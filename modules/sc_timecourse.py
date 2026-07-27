"""Sample-aware temporal analysis for single-cell RNA-seq experiments.

Measured sampling time and pseudotime answer different questions.  This module
uses the former: it aggregates cells at the biological-sample level before
testing cell-type composition and gene-expression dynamics.  It deliberately
does not treat cells from one sample as independent replicates.
"""

import os
import re

import numpy as np

from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping


def _natural_key(value):
    """Natural-sort a label without an extra dependency."""
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", str(value))]


def _time_value(label):
    """Return a sortable value for common labels such as D0, day7 and 12h."""
    text = str(label).strip().lower()
    match = re.match(
        r"^(?:(day|d|hour|hr|h|week|wk|w|minute|min|m)\s*)?"
        r"(-?\d+(?:\.\d+)?)\s*"
        r"(?:(day|d|hour|hr|h|week|wk|w|minute|min|m))?$",
        text,
    )
    if not match:
        return None
    prefix, number, suffix = match.groups()
    unit = prefix or suffix or ""
    multipliers = {
        "": 1.0, "m": 1.0 / 60.0, "min": 1.0 / 60.0,
        "h": 1.0, "hr": 1.0, "hour": 1.0,
        "d": 24.0, "day": 24.0,
        "w": 24.0 * 7, "wk": 24.0 * 7, "week": 24.0 * 7,
    }
    return float(number) * multipliers[unit]


def order_timepoints(values, requested_order=""):
    """Return ordered string labels and numeric plotting positions.

    Explicit order takes precedence.  Otherwise numeric columns and familiar
    labels (D0/D3/D7, 0h/12h) are sorted chronologically; arbitrary labels use
    natural sorting and are treated as ordered experimental categories.
    """
    import pandas as pd

    series = pd.Series(values).dropna()
    if series.empty:
        raise ValueError("时间列没有有效值")
    observed = list(dict.fromkeys(series.astype(str).tolist()))
    if len(observed) < 3:
        raise ValueError(f"至少需要 3 个不同时间点，当前只有 {len(observed)} 个")

    explicit = [item.strip() for item in str(requested_order or "").split(",") if item.strip()]
    if explicit:
        unknown = [item for item in explicit if item not in observed]
        missing = [item for item in observed if item not in explicit]
        if unknown or missing:
            detail = []
            if unknown:
                detail.append("不存在: " + ", ".join(unknown))
            if missing:
                detail.append("未包含: " + ", ".join(missing))
            raise ValueError("时间顺序必须覆盖所有时间点（" + "；".join(detail) + "）")
        return explicit, {label: float(index) for index, label in enumerate(explicit)}

    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().all():
        numeric_by_label = {
            str(label): float(pd.to_numeric(pd.Series([label]), errors="coerce").iloc[0])
            for label in observed
        }
        labels = sorted(observed, key=lambda item: numeric_by_label[item])
        return labels, {label: numeric_by_label[label] for label in labels}

    parsed = {label: _time_value(label) for label in observed}
    if all(value is not None for value in parsed.values()):
        labels = sorted(observed, key=lambda item: parsed[item])
        return labels, {label: float(parsed[label]) for label in labels}

    labels = sorted(observed, key=_natural_key)
    return labels, {label: float(index) for index, label in enumerate(labels)}


def benjamini_hochberg(pvalues):
    """BH-adjust p values while preserving NaN for unavailable inference."""
    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted
    indices = np.where(valid)[0]
    ordered = indices[np.argsort(values[indices])]
    n_valid = len(ordered)
    running = 1.0
    for rank in range(n_valid, 0, -1):
        index = ordered[rank - 1]
        running = min(running, values[index] * n_valid / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def _temporal_statistics(time_values, values, time_labels, min_replicates, allow_inference):
    """Return monotonic direction plus a sample-level, non-parametric screen."""
    from scipy.stats import kruskal, spearmanr

    time_values = np.asarray(time_values, dtype=float)
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(time_values) & np.isfinite(values)
    time_values = time_values[finite]
    values = values[finite]
    labels = np.asarray(time_labels, dtype=str)[finite]
    if len(values) < 3:
        return np.nan, np.nan, "insufficient_units"

    try:
        rho = float(spearmanr(time_values, values).statistic)
    except Exception:
        rho = np.nan

    unique_labels = list(dict.fromkeys(labels.tolist()))
    group_values = [values[labels == label] for label in unique_labels]
    has_replicates = all(len(group) >= int(min_replicates) for group in group_values)
    if not allow_inference or len(unique_labels) < 3 or not has_replicates:
        return rho, np.nan, "descriptive_only"
    try:
        p_value = float(kruskal(*group_values).pvalue)
    except Exception:
        return rho, np.nan, "test_unavailable"
    if not np.isfinite(p_value):
        return rho, np.nan, "test_unavailable"
    return rho, p_value, "sample_level_kruskal"


def _matrix_variance(matrix):
    """Column variance without turning a sparse count matrix dense."""
    from scipy import sparse

    if sparse.issparse(matrix):
        mean = np.asarray(matrix.mean(axis=0)).ravel()
        second = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
    else:
        dense = np.asarray(matrix, dtype=float)
        mean = np.nanmean(dense, axis=0)
        second = np.nanmean(np.square(dense), axis=0)
    return np.maximum(second - np.square(mean), 0)


def _is_count_like_matrix(matrix):
    """Return whether a matrix is safe to describe as raw non-negative counts.

    Import compatibility may populate ``layers['counts']`` by copying ``X``.
    The layer name alone therefore cannot justify count-based pseudobulk
    inference.  Check a bounded set of non-zero values instead of densifying a
    large sparse matrix.
    """
    from scipy import sparse

    if sparse.issparse(matrix):
        values = np.asarray(matrix.data, dtype=float)
    else:
        values = np.asarray(matrix, dtype=float).reshape(-1)
    if values.size == 0:
        return False
    if values.size > 200_000:
        step = max(1, values.size // 200_000)
        values = values[::step][:200_000]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return False
    nonnegative = float(np.mean(values >= 0))
    integer_like = float(np.mean(np.isclose(values, np.round(values))))
    return nonnegative >= 0.999 and integer_like >= 0.995


def select_temporal_genes(adata, matrix, max_genes):
    """Prefer existing HVGs, then rank them by observed variance."""
    max_genes = max(10, int(max_genes))
    candidates = np.arange(adata.n_vars)
    if "highly_variable" in adata.var.columns:
        hvg = np.flatnonzero(np.asarray(adata.var["highly_variable"], dtype=bool))
        if len(hvg):
            candidates = hvg
    variances = _matrix_variance(matrix)
    ranked = candidates[np.argsort(variances[candidates])[::-1]]
    return ranked[:min(max_genes, len(ranked))]


def _matrix_row_expression(matrix, indices, gene_indices, use_counts):
    """Aggregate a sample × celltype unit onto a stable expression scale."""
    selected = matrix[indices, :]
    if use_counts:
        total = float(np.asarray(selected.sum()).ravel()[0])
        gene_sums = np.asarray(selected[:, gene_indices].sum(axis=0)).ravel()
        return np.log1p(gene_sums / max(total, 1.0) * 1_000_000.0)
    return np.asarray(selected[:, gene_indices].mean(axis=0)).ravel()


class SCTimecourseAnalysis(BaseAnalysis):
    """Analyze measured temporal dynamics using biological samples as units."""

    MODULE_NAME = "sc_timecourse"
    DISPLAY_NAME = "单细胞时序动态"
    DESCRIPTION = "多时间点细胞组成、样本级伪 bulk 基因动态与时间点 UMAP"
    # A user-supplied celltype/annotation column is preferred; leiden is only
    # a fallback, so it must not be a hard AnnData input contract.
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        timepoint_key = str(self.params.get("timepoint_key", "timepoint") or "").strip()
        if not timepoint_key:
            return "需要指定真实采样时间列（如 timepoint、day、hour）"
        if timepoint_key not in adata.obs.columns:
            return f"时间列 '{timepoint_key}' 不在 adata.obs 中"
        return None

    def _resolve_celltype_key(self, adata):
        requested = str(self.params.get("celltype_key", "celltype") or "").strip()
        key, info = resolve_obs_grouping(
            adata, requested,
            fallbacks=("celltype", "annotation", "leiden"),
            max_categories=60, max_numeric_categories=30, require_multiple=True,
        )
        if key is None:
            raise ValueError(
                f"细胞类型列 '{requested}' 不可用于时序分析：{info.get('reason', '')}"
            )
        if key != requested:
            self.progress(-1, f"细胞类型列已改用 '{key}'：{info.get('requested_reason', '')}")
        return key

    def _resolve_optional_group(self, adata, requested, label):
        requested = str(requested or "").strip()
        if not requested:
            return None
        if requested not in adata.obs.columns:
            raise ValueError(f"{label}列 '{requested}' 不在 adata.obs 中")
        raw = adata.obs[requested]
        missing = raw.isna() | raw.astype(str).str.strip().eq("")
        if bool(missing.any()):
            raise ValueError(
                f"{label}列 '{requested}' 有 {int(missing.sum())} 个缺失值；"
                "请先补齐或不进行该分层，不能把缺失值作为一个条件。"
            )
        key, info = resolve_obs_grouping(
            adata, requested, max_categories=60,
            max_numeric_categories=30, require_multiple=True,
        )
        if key is None:
            raise ValueError(f"{label}列 '{requested}' 不可用：{info.get('reason', '')}")
        return key

    def _sample_column(self, adata, requested, timepoint_key, condition_key=None):
        import pandas as pd

        requested = str(requested or "").strip()
        if not requested:
            return None, "未提供 sample_id；所有时间动态均为描述性，不报告推断性 p 值。"
        if requested not in adata.obs.columns:
            return None, f"未找到样本列 '{requested}'；仅输出描述性时间趋势，不报告 p 值。"

        technical_batch_names = {
            "batch", "technical_batch", "sequencing_batch", "library_batch",
        }
        if (requested.lower() in technical_batch_names
                and not bool(self.params.get("confirm_batch_is_biological_sample", False))):
            return None, (
                f"样本列 '{requested}' 看起来是技术 batch；未确认其代表独立生物学样本，"
                "已关闭推断性 p 值/FDR。"
            )

        raw = adata.obs[requested]
        missing = raw.isna() | raw.astype(str).str.strip().eq("")
        if bool(missing.any()):
            return None, (
                f"样本列 '{requested}' 有 {int(missing.sum())} 个缺失值；"
                "不能将缺失值合并为同一伪样本，已关闭推断性 p 值/FDR。"
            )

        values = raw.astype(str).str.strip()
        per_sample = values.value_counts()
        if values.nunique() < 2:
            return None, f"样本列 '{requested}' 的独立样本数不足；仅输出描述性时间趋势。"
        if int(per_sample.median()) < 2 or values.nunique() > max(20, int(0.5 * adata.n_obs)):
            return None, (
                f"样本列 '{requested}' 接近一细胞一个 ID，不能作为生物学重复；"
                "仅输出描述性时间趋势。"
            )

        design = pd.DataFrame({
            "sample": values,
            "timepoint": adata.obs[timepoint_key].astype(str),
        })
        timepoints_per_sample = design.groupby("sample", observed=True)["timepoint"].nunique()
        if bool((timepoints_per_sample > 1).any()):
            return None, (
                f"样本列 '{requested}' 有同一样本跨多个时间点；当前模型不支持配对/混合效应推断，"
                "已降级为描述性结果。"
            )
        if condition_key:
            design["condition"] = adata.obs[condition_key].astype(str)
            conditions_per_sample = design.groupby("sample", observed=True)["condition"].nunique()
            if bool((conditions_per_sample > 1).any()):
                return None, (
                    f"样本列 '{requested}' 有同一样本跨多个条件；当前模型不支持配对/交互推断，"
                "已降级为描述性结果。"
            )
        return requested, None

    def _longitudinal_subject_warning(self, adata, timepoint_key, sample_key):
        """Detect common subject IDs that make library-level units paired.

        A distinct library/sample ID at each visit does not by itself make
        repeated biopsies from the same donor independent.  This module has no
        mixed-effects or paired test, so retain the sample-level tables but
        disable inferential p values when a conventional subject column shows
        repeated measurements across timepoints.
        """
        import pandas as pd

        if not sample_key:
            return None
        subject_names = {
            "donor_id", "donor", "subject_id", "subject", "patient_id",
            "patient", "participant_id", "participant", "individual_id",
        }
        candidates = [
            column for column in adata.obs.columns
            if str(column).lower() in subject_names and str(column) != sample_key
        ]
        for column in candidates:
            values = adata.obs[column]
            missing = values.isna() | values.astype(str).str.strip().eq("")
            if bool(missing.any()):
                continue
            design = pd.DataFrame({
                "subject": values.astype(str).str.strip(),
                "timepoint": adata.obs[timepoint_key].astype(str),
            })
            visits = design.groupby("subject", observed=True)["timepoint"].nunique()
            if bool((visits > 1).any()):
                return (
                    f"检测到 '{column}' 显示同一受试者跨多个时间点；"
                    "当前模型不支持配对/混合效应推断，已保留样本级表格但关闭 p 值/FDR。"
                )
        return None

    def _save_table(self, frame, results_dir, filename, label):
        path = os.path.join(results_dir, filename)
        frame.to_csv(path, index=False)  # result_files: caller registers returned artifact
        return {"file_path": path, "file_type": "csv", "category": "table", "label": label}

    def _composition_rows(self, obs, time_labels, time_values, celltype_key, sample_key, condition_key):
        import pandas as pd

        work = pd.DataFrame({
            "timepoint": obs["_sc_timepoint"].astype(str),
            "celltype": obs[celltype_key].astype(str),
            "condition": obs[condition_key].astype(str) if condition_key else "All",
        }, index=obs.index)
        work["unit"] = obs[sample_key].astype(str) if sample_key else work["timepoint"]
        unit_columns = ["condition", "timepoint", "unit"]
        units = work.groupby(unit_columns, observed=True).size().rename("n_cells").reset_index()
        counts = work.groupby(unit_columns + ["celltype"], observed=True).size().rename("cell_count").reset_index()

        celltypes = sorted(work["celltype"].unique().tolist(), key=_natural_key)
        full = units.merge(pd.DataFrame({"celltype": celltypes}), how="cross")
        composition = full.merge(counts, on=unit_columns + ["celltype"], how="left")
        composition["cell_count"] = composition["cell_count"].fillna(0).astype(int)
        composition["proportion"] = composition["cell_count"] / composition["n_cells"].clip(lower=1)
        composition["time_position"] = composition["timepoint"].map(time_values).astype(float)
        return composition, units, celltypes

    def _composition_trends(self, composition, time_labels, sample_key, min_replicates,
                            allow_inference=True):
        import pandas as pd

        rows = []
        for (condition, celltype), frame in composition.groupby(["condition", "celltype"], observed=True):
            counts_by_time = frame.groupby("timepoint", observed=True)["unit"].nunique()
            inference_ready = bool(sample_key) and bool(allow_inference) and all(
                int(counts_by_time.get(label, 0)) >= int(min_replicates) for label in time_labels
            )
            rho, pvalue, status = _temporal_statistics(
                frame["time_position"], frame["proportion"], frame["timepoint"],
                min_replicates, inference_ready,
            )
            rows.append({
                "condition": str(condition), "celltype": str(celltype),
                "spearman_rho": rho, "kruskal_pvalue": pvalue,
                "inference_status": status,
                "n_units": int(frame["unit"].nunique()),
                "min_replicates_per_timepoint": int(counts_by_time.min()) if len(counts_by_time) else 0,
            })
        trends = pd.DataFrame(rows)
        if not trends.empty:
            trends["fdr_bh"] = benjamini_hochberg(trends["kruskal_pvalue"].values)
            trends["direction"] = np.select(
                [trends["spearman_rho"] >= 0.25, trends["spearman_rho"] <= -0.25],
                ["increasing", "decreasing"], default="non_monotonic_or_flat",
            )
        return trends

    def _expression_trends(self, adata, obs, celltype_key, sample_key, condition_key,
                           time_labels, time_values, min_cells, min_replicates,
                           max_genes, use_counts, allow_inference=True):
        import pandas as pd

        matrix = adata.layers["counts"] if use_counts else adata.X
        gene_indices = select_temporal_genes(adata, matrix, max_genes)
        gene_names = adata.var_names[gene_indices].astype(str).tolist()
        work = pd.DataFrame({
            "timepoint": obs["_sc_timepoint"].astype(str),
            "celltype": obs[celltype_key].astype(str),
            "condition": obs[condition_key].astype(str) if condition_key else "All",
        }, index=obs.index)
        work["unit"] = obs[sample_key].astype(str) if sample_key else work["timepoint"]
        group_columns = ["condition", "timepoint", "unit", "celltype"]
        records = []
        for keys, cell_indices in work.groupby(group_columns, observed=True).indices.items():
            if len(cell_indices) < int(min_cells):
                continue
            values = _matrix_row_expression(matrix, cell_indices, gene_indices, use_counts)
            record = dict(zip(group_columns, [str(item) for item in keys]))
            record["n_cells"] = int(len(cell_indices))
            record.update(dict(zip(gene_names, values)))
            records.append(record)
        if not records:
            return pd.DataFrame(), pd.DataFrame(), {}, use_counts

        pseudobulk = pd.DataFrame(records)
        pseudobulk["time_position"] = pseudobulk["timepoint"].map(time_values).astype(float)
        trend_rows = []
        profiles = {}
        for (condition, celltype), frame in pseudobulk.groupby(["condition", "celltype"], observed=True):
            replicate_counts = frame.groupby("timepoint", observed=True)["unit"].nunique()
            inference_ready = bool(sample_key) and bool(allow_inference) and use_counts and all(
                int(replicate_counts.get(label, 0)) >= int(min_replicates) for label in time_labels
            )
            mean_by_time = frame.groupby("timepoint", observed=True)[gene_names].mean().reindex(time_labels)
            for gene in gene_names:
                rho, pvalue, status = _temporal_statistics(
                    frame["time_position"], frame[gene], frame["timepoint"],
                    min_replicates, inference_ready,
                )
                trend_rows.append({
                    "condition": str(condition), "celltype": str(celltype), "gene": gene,
                    "spearman_rho": rho, "kruskal_pvalue": pvalue,
                    "inference_status": status,
                    "n_units": int(frame["unit"].nunique()),
                    "min_replicates_per_timepoint": int(replicate_counts.min()) if len(replicate_counts) else 0,
                })
                profiles[(str(condition), str(celltype), gene)] = mean_by_time[gene].to_numpy(dtype=float)
        trends = pd.DataFrame(trend_rows)
        if not trends.empty:
            trends["fdr_bh"] = np.nan
            for _, indexes in trends.groupby(["condition", "celltype"], observed=True).groups.items():
                trends.loc[indexes, "fdr_bh"] = benjamini_hochberg(
                    trends.loc[indexes, "kruskal_pvalue"].values
                )
            trends["direction"] = np.select(
                [trends["spearman_rho"] >= 0.25, trends["spearman_rho"] <= -0.25],
                ["increasing", "decreasing"], default="non_monotonic_or_flat",
            )
        return pseudobulk, trends, profiles, use_counts

    def run(self, input_path):
        import json
        import pandas as pd
        from modules.native_figures import heatmap_figure, line_figure, umap_figure

        self.progress(5, "读取单细胞数据与时间元数据...")
        adata = self.load_adata(input_path)
        timepoint_key = str(self.params.get("timepoint_key", "timepoint") or "").strip()
        if timepoint_key not in adata.obs.columns:
            raise ValueError(f"时间列 '{timepoint_key}' 不在 adata.obs 中；请先在 h5ad.obs 中提供真实采样时间。")
        missing_timepoints = adata.obs[timepoint_key].isna() | adata.obs[timepoint_key].astype(str).str.strip().eq("")
        if bool(missing_timepoints.any()):
            raise ValueError(
                f"时间列 '{timepoint_key}' 有 {int(missing_timepoints.sum())} 个缺失值；"
                "请先补齐真实采样时间，不能把缺失值当作一个时间点。"
            )
        time_labels, time_values = order_timepoints(
            adata.obs[timepoint_key], self.params.get("time_order", "")
        )
        adata.obs["_sc_timepoint"] = pd.Categorical(
            adata.obs[timepoint_key].astype(str), categories=time_labels, ordered=True,
        )
        celltype_key = self._resolve_celltype_key(adata)
        condition_key = self._resolve_optional_group(adata, self.params.get("condition_key", ""), "条件")
        sample_key, sample_warning = self._sample_column(
            adata, self.params.get("sample_key", "sample_id"), timepoint_key, condition_key,
        )
        min_cells = max(1, int(self.params.get("min_cells_per_celltype_sample", 20)))
        # A single sample per timepoint is not a biological replicate and
        # cannot support a between-timepoint inferential screen.  Clamp here
        # as well as in the UI because tasks can be submitted through the API
        # or AI tools without the browser form.
        min_replicates = max(2, int(self.params.get("min_replicates_per_timepoint", 2)))
        max_genes = max(20, int(self.params.get("max_genes", 300)))
        warnings = []
        if sample_warning:
            warnings.append(sample_warning)
        technical_batch_names = {
            "batch", "technical_batch", "sequencing_batch", "library_batch",
        }
        technical_timepoint = timepoint_key.lower() in technical_batch_names
        if (technical_timepoint
                and not bool(self.params.get("confirm_batch_is_biological_timepoint", False))):
            warnings.append(
                f"时间列 '{timepoint_key}' 看起来是技术 batch；未确认其代表真实采样时间，"
                "已关闭推断性 p 值/FDR。"
            )
        elif technical_timepoint:
            warnings.append(
                f"已人工确认 '{timepoint_key}' 代表真实采样时间；"
                "请勿再将同一列用于批次校正。"
            )
        longitudinal_warning = self._longitudinal_subject_warning(
            adata, timepoint_key, sample_key,
        )
        if longitudinal_warning:
            warnings.append(longitudinal_warning)
        allow_sample_inference = (
            bool(sample_key)
            and not bool(longitudinal_warning)
            and (not technical_timepoint or bool(self.params.get("confirm_batch_is_biological_timepoint", False)))
        )
        if "batch" in adata.obs.columns and timepoint_key != "batch":
            overlap = pd.crosstab(adata.obs[timepoint_key].astype(str), adata.obs["batch"].astype(str))
            if len(overlap) and (overlap.gt(0).sum(axis=1) <= 1).all():
                warnings.append("每个时间点仅对应一个 technical batch，时间与技术批次完全混杂；需谨慎解释时间差异。")
        if condition_key:
            warnings.append("condition 当前仅作分层展示；本模块不进行 time × condition 交互检验。")

        self.progress(20, "汇总样本层面的细胞组成...")
        composition, units, celltypes = self._composition_rows(
            adata.obs, time_labels, time_values, celltype_key, sample_key, condition_key
        )
        composition_trends = self._composition_trends(
            composition, time_labels, sample_key, min_replicates, allow_sample_inference,
        )
        if (not composition_trends.empty
                and composition_trends["inference_status"].eq("test_unavailable").any()):
            warnings.append(
                "至少一个细胞类型的组成比例恒定或无法计算 Kruskal 检验；"
                "对应 p 值/FDR 未报告。"
            )
        if sample_key:
            design_index = pd.MultiIndex.from_product(
                [sorted(units["condition"].astype(str).unique()), time_labels],
                names=["condition", "timepoint"],
            )
            replicate_counts = units.groupby(["condition", "timepoint"], observed=True)["unit"].nunique()
            replicate_counts = replicate_counts.reindex(design_index, fill_value=0)
            if bool((replicate_counts < min_replicates).any()):
                warnings.append(
                    f"至少一个 condition × timepoint 少于 {min_replicates} 个独立样本；"
                    "对应趋势仅为描述性，不报告 p 值/FDR。"
                )

        results_dir = os.path.join(self.project_dir, "results")
        os.makedirs(results_dir, exist_ok=True)
        plots_dir = self.ensure_plots_dir()
        result_files = [
            self._save_table(composition, results_dir, "sc_timecourse_composition.csv", "Sample-level cell-type proportions"),
            self._save_table(composition_trends, results_dir, "sc_timecourse_composition_trends.csv", "Cell-type temporal composition trends"),
        ]
        units_export = units.copy()
        if "total_counts" in adata.obs.columns:
            unit_meta = pd.DataFrame({
                "condition": adata.obs[condition_key].astype(str) if condition_key else "All",
                "timepoint": adata.obs["_sc_timepoint"].astype(str),
                "unit": adata.obs[sample_key].astype(str) if sample_key else adata.obs["_sc_timepoint"].astype(str),
                "total_counts": pd.to_numeric(adata.obs["total_counts"], errors="coerce"),
            })
            medians = unit_meta.groupby(["condition", "timepoint", "unit"], observed=True)["total_counts"].median().rename("median_total_counts").reset_index()
            units_export = units_export.merge(medians, on=["condition", "timepoint", "unit"], how="left")
        result_files.append(self._save_table(units_export, results_dir, "sc_timecourse_sample_summary.csv", "Timepoint sample summary"))

        if self.params.get("show_timepoint_umap", True) and "X_umap" in adata.obsm:
            fig = umap_figure(adata, "_sc_timepoint", title="UMAP by measured timepoint")
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, "sc_timecourse_umap.png", "umap", "Timepoint UMAP", formats=("png", "svg"), dpi=300,
            ))

        mean_composition = composition.groupby(["condition", "timepoint", "celltype"], observed=True)["proportion"].mean().reset_index()
        top_celltypes = (mean_composition.groupby("celltype", observed=True)["proportion"].mean()
                         .sort_values(ascending=False).head(max(2, int(self.params.get("top_celltypes", 8)))).index.tolist())
        if self.params.get("show_composition_trajectory", True) and top_celltypes:
            series, labels = [], []
            time_positions = [time_values[label] for label in time_labels]
            for condition in mean_composition["condition"].unique():
                for celltype in top_celltypes:
                    profile = mean_composition[(mean_composition["condition"] == condition) & (mean_composition["celltype"] == celltype)]
                    profile = profile.set_index("timepoint").reindex(time_labels)["proportion"].fillna(0)
                    series.append(profile.values)
                    labels.append(str(celltype) if condition == "All" else f"{condition} · {celltype}")
            fig = line_figure(time_positions, series, title="Cell-type composition across time", x_label="Timepoint", y_label="Mean proportion", labels=labels, markers=["o"] * len(series))
            fig.axes[0].set_xticks(time_positions, time_labels)
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, "sc_timecourse_composition_trajectory.png", "bar", "Cell-type composition trajectories", formats=("png", "svg"), dpi=300,
            ))
        if self.params.get("show_composition_heatmap", True) and top_celltypes:
            rows, row_labels = [], []
            for condition in mean_composition["condition"].unique():
                for timepoint in time_labels:
                    row = mean_composition[(mean_composition["condition"] == condition) & (mean_composition["timepoint"] == timepoint)].set_index("celltype")["proportion"]
                    rows.append(row.reindex(top_celltypes).fillna(0).values)
                    row_labels.append(timepoint if condition == "All" else f"{condition} · {timepoint}")
            fig = heatmap_figure(np.asarray(rows), x_labels=top_celltypes, y_labels=row_labels,
                                 title="Cell-type composition heatmap", x_label="Cell type", y_label="Timepoint / condition",
                                 colorbar_label="Proportion", vmin=0, vmax=max(1e-6, float(np.max(rows))))
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, "sc_timecourse_composition_heatmap.png", "heatmap", "Cell-type composition heatmap", formats=("png", "svg"), dpi=300,
            ))

        pseudobulk = pd.DataFrame()
        gene_trends = pd.DataFrame()
        use_counts = bool("counts" in adata.layers and _is_count_like_matrix(adata.layers["counts"]))
        if self.params.get("enable_gene_trends", True):
            self.progress(55, "构建 sample × celltype 伪 bulk 并检测基因动态...")
            pseudobulk, gene_trends, profiles, use_counts = self._expression_trends(
                adata, adata.obs, celltype_key, sample_key, condition_key,
                time_labels, time_values, min_cells, min_replicates, max_genes, use_counts,
                allow_sample_inference,
            )
            if not use_counts:
                warnings.append("未找到非负近似整数的原始 counts 层；基因曲线仅基于标准化细胞均值，已禁用基因推断性 p 值。")
            elif not gene_trends.empty:
                warnings.append(
                    f"基因 FDR 仅在每个 condition × celltype 的最多 {max_genes} 个高变/高方差候选基因内校正，"
                    "属于探索性筛选，不代表全转录组 FDR。"
                )
            if (not gene_trends.empty
                    and gene_trends["inference_status"].eq("test_unavailable").any()):
                warnings.append(
                    "至少一个候选基因的样本级 Kruskal 检验不可计算；"
                    "对应 p 值/FDR 未报告。"
                )
            if not pseudobulk.empty:
                result_files.append(self._save_table(pseudobulk, results_dir, "sc_timecourse_pseudobulk_expression.csv", "Sample × cell-type pseudobulk expression"))
            if not gene_trends.empty:
                result_files.append(self._save_table(gene_trends, results_dir, "sc_timecourse_gene_trends.csv", "Cell-type temporal gene trends"))
                top_n = max(5, int(self.params.get("top_dynamic_genes", 30)))
                ranking = gene_trends.copy()
                ranking["rank_fdr"] = ranking["fdr_bh"].fillna(2.0)
                ranking["rank_effect"] = ranking["spearman_rho"].abs().fillna(0.0)
                top = ranking.sort_values(["rank_fdr", "rank_effect"], ascending=[True, False]).head(top_n)
                matrix_rows, gene_labels = [], []
                for row in top.itertuples(index=False):
                    values = np.asarray(profiles.get((row.condition, row.celltype, row.gene), []), dtype=float)
                    if len(values) != len(time_labels):
                        continue
                    std = float(np.nanstd(values))
                    scaled = (values - np.nanmean(values)) / std if std > 1e-12 else np.zeros_like(values)
                    matrix_rows.append(scaled)
                    gene_labels.append(f"{row.condition} · {row.celltype} · {row.gene}")
                if matrix_rows and self.params.get("show_gene_heatmap", True):
                    fig = heatmap_figure(np.asarray(matrix_rows), x_labels=time_labels, y_labels=gene_labels,
                                         title="Top cell-type temporal gene dynamics", x_label="Timepoint", y_label="Cell type · gene",
                                         colorbar_label="Row z-score", vmin=-2, vmax=2)
                    result_files.extend(self.save_matplotlib_figure(
                        fig, plots_dir, "sc_timecourse_gene_dynamics_heatmap.png", "heatmap", "Temporal pseudobulk gene dynamics", formats=("png", "svg"), dpi=300,
                    ))

        composition_inference_available = bool(
            not composition_trends.empty
            and composition_trends["inference_status"].eq("sample_level_kruskal").any()
        )
        gene_inference_available = bool(
            not gene_trends.empty
            and gene_trends["inference_status"].eq("sample_level_kruskal").any()
        )
        metadata = {
            "timepoint_key": timepoint_key, "timepoint_order": time_labels,
            "time_values": time_values, "sample_key": sample_key,
            "condition_key": condition_key, "celltype_key": celltype_key,
            "inference_unit": "sample" if sample_key else "none (descriptive only)",
            "sample_inference_allowed": bool(allow_sample_inference),
            "gene_expression_source": "counts_pseudobulk_cpm" if use_counts else "normalized_cell_mean",
            "counts_layer_is_raw": bool(use_counts),
            "composition_inference_available": composition_inference_available,
            "gene_inference_available": gene_inference_available,
            "warnings": warnings,
        }
        metadata_path = os.path.join(results_dir, "sc_timecourse_metadata.json")
        with open(metadata_path, "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, ensure_ascii=False, indent=2)
        result_files.append({"file_path": metadata_path, "file_type": "json", "category": "info", "label": "Time-course analysis metadata and cautions"})

        self.progress(90, "保存带时间顺序的 AnnData 输出...")
        output_path = self.save_output(adata, "sc_timecourse")
        n_composition_sig = int((composition_trends.get("fdr_bh", pd.Series(dtype=float)) < 0.05).sum()) if not composition_trends.empty else 0
        n_gene_sig = int((gene_trends.get("fdr_bh", pd.Series(dtype=float)) < 0.05).sum()) if not gene_trends.empty else 0
        self.progress(100, "单细胞时序动态分析完成")
        return {
            "output_adata": output_path,
            "result_files": result_files,
            'summary': {
                "n_cells": int(adata.n_obs), "n_timepoints": int(len(time_labels)),
                "n_celltypes": int(len(celltypes)), "n_samples": int(units["unit"].nunique()) if sample_key else 0,
                "n_composition_significant": int(n_composition_sig),
                "n_gene_significant": int(n_gene_sig),
                "timepoint_key": timepoint_key, "sample_key": sample_key or "",
                "celltype_key": celltype_key, "condition_key": condition_key or "",
                "inference_available": bool(composition_inference_available or gene_inference_available),
                "composition_inference_available": composition_inference_available,
                "gene_inference_available": gene_inference_available,
                "counts_layer_is_raw": bool(use_counts),
                "sample_inference_allowed": bool(allow_sample_inference),
                "gene_expression_source": "counts_pseudobulk_cpm" if use_counts else "normalized_cell_mean",
                "warnings": warnings,
            },
        }
