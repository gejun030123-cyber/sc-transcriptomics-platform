"""Standard CSV export and sample-level pseudobulk DEG for batched scRNA-seq."""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from modules.base import BaseAnalysis
from modules.sc_batch import (
    _available_path,
    aggregate_pseudobulk_counts,
    build_comparison_plan,
    sample_design_from_obs,
)
from modules.sc_de_utils import resolve_cell_grouping


CSV_PACKAGE_FOLDERS = {
    "proportions": "01_group_proportions",
    "expression": "02_gene_expression",
    "deg": "03_differential_expression",
    "go": "04_go_enrichment",
}


def _as_bool(value, default=True):
    if value is None:
        return bool(default)
    if isinstance(value, str):
        return value.strip().lower() not in {"", "0", "false", "no", "off"}
    return bool(value)


def _pseudobulk_design_is_identifiable(meta, batch_key="batch"):
    """Check whether ``~ batch + condition`` has full column rank.

    A batch perfectly aligned with condition cannot be adjusted by DESeq2;
    reporting a coefficient in that case would be a statistical fiction.
    """
    import numpy as np

    if batch_key not in meta.columns or meta[batch_key].nunique(dropna=True) <= 1:
        return True, ""
    if meta["condition"].nunique(dropna=True) < 2:
        return False, "当前 pseudobulk 单元少于两个条件"
    batch = pd.get_dummies(meta[batch_key].astype(str), prefix="batch", drop_first=True)
    condition = pd.get_dummies(
        meta["condition"].astype(str), prefix="condition", drop_first=True,
    )
    design = pd.concat(
        [pd.Series(1.0, index=meta.index, name="intercept"), batch, condition], axis=1,
    ).astype(float)
    rank = int(np.linalg.matrix_rank(design.to_numpy(dtype=float)))
    if rank < design.shape[1]:
        return False, "batch 与 condition 完全混杂，~ batch + condition 不可识别"
    return True, ""


def _safe_folder_name(value, fallback="sc_batch_results"):
    text = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in str(value or ""))
    return text.strip("._-") or fallback


def batch_csv_package_dirs(project_dir, export_folder="sc_batch_results", included_keys=None):
    """Create only the requested DOCX-defined directories for one result package."""
    root = os.path.join(project_dir, "results", _safe_folder_name(export_folder))
    selected = tuple(included_keys or CSV_PACKAGE_FOLDERS.keys())
    unknown = [key for key in selected if key not in CSV_PACKAGE_FOLDERS]
    if unknown:
        raise ValueError("未知 CSV 结果目录: " + ", ".join(unknown))
    directories = {
        key: os.path.join(root, CSV_PACKAGE_FOLDERS[key])
        for key in selected
    }
    for path in directories.values():
        os.makedirs(path, exist_ok=True)
    return root, directories


def _gene_metadata(adata):
    gene_id = (
        adata.var["gene_id"].astype(str).to_numpy()
        if "gene_id" in adata.var.columns else adata.var_names.astype(str).to_numpy()
    )
    gene_name = (
        adata.var["gene_name"].astype(str).to_numpy()
        if "gene_name" in adata.var.columns else adata.var_names.astype(str).to_numpy()
    )
    return gene_id, gene_name


def _logcpm(counts):
    counts = np.asarray(counts, dtype=float)
    library_sizes = counts.sum(axis=1, keepdims=True)
    return np.log2(counts / np.maximum(library_sizes, 1.0) * 1_000_000.0 + 1.0)


def _welch_logcpm_deg(counts, group_1_mask, group_2_mask, gene_ids, gene_names,
                      comparison, group_1, group_2, cluster, min_samples):
    """A transparent fallback when DESeq2 cannot be fit on pseudobulk counts."""
    from scipy.stats import ttest_ind
    from statsmodels.stats.multitest import multipletests

    counts = np.asarray(counts, dtype=float)
    group_1_counts = counts[group_1_mask]
    group_2_counts = counts[group_2_mask]
    if len(group_1_counts) < int(min_samples) or len(group_2_counts) < int(min_samples):
        return pd.DataFrame()
    transformed = _logcpm(counts)
    test = ttest_ind(
        transformed[group_1_mask], transformed[group_2_mask], axis=0,
        equal_var=False, nan_policy="omit",
    )
    mean_cpm_1 = np.expm1(transformed[group_1_mask] * np.log(2)).mean(axis=0)
    mean_cpm_2 = np.expm1(transformed[group_2_mask] * np.log(2)).mean(axis=0)
    log2fc = np.log2((mean_cpm_1 + 0.5) / (mean_cpm_2 + 0.5))
    pvalues = np.asarray(test.pvalue, dtype=float)
    padj = np.full(pvalues.shape, np.nan, dtype=float)
    valid = np.isfinite(pvalues)
    if bool(valid.any()):
        padj[valid] = multipletests(pvalues[valid], method="fdr_bh")[1]
    return pd.DataFrame({
        "comparison": comparison,
        "group_1": group_1,
        "group_2": group_2,
        "cluster": cluster,
        "gene_id": gene_ids,
        "gene": gene_names,
        "base_mean_count": counts.mean(axis=0),
        "log2FC": log2fc,
        "statistic": np.asarray(test.statistic, dtype=float),
        "pvalue": pvalues,
        "padj": padj,
        "n_samples_group_1": int(group_1_mask.sum()),
        "n_samples_group_2": int(group_2_mask.sum()),
        "method": "welch_t_log2cpm",
    })


def _deseq2_deg(counts, conditions, group_1, group_2, gene_ids, gene_names,
                comparison, cluster, batch_values=None, batch_key="batch"):
    """Fit PyDESeq2 on pseudobulk counts; callers may fall back safely."""
    from pydeseq2.ds import DeseqStats

    dds = _make_deseq2_dataset(
        counts, conditions, batch_values=batch_values, batch_key=batch_key,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", str(group_1), str(group_2)])
    stats.summary()
    result = stats.results_df.copy()
    result.index = result.index.astype(str)
    columns = [f"g{index}" for index in range(np.asarray(counts).shape[1])]
    if result.index.is_unique and set(columns).issubset(set(result.index)):
        result = result.reindex(columns)
    else:
        raise RuntimeError(
            "PyDESeq2 结果缺少可与输入基因一一对应的索引，已停止以避免错配基因。"
        )
    if len(result) != len(gene_ids) or len(result) != len(gene_names):
        raise RuntimeError("PyDESeq2 结果行数与输入基因元数据不一致。")
    result = result.reset_index(drop=True)
    return pd.DataFrame({
        "comparison": comparison,
        "group_1": group_1,
        "group_2": group_2,
        "cluster": cluster,
        "gene_id": gene_ids,
        "gene": gene_names,
        "base_mean_count": result.get("baseMean", np.nan),
        "log2FC": result.get("log2FoldChange", np.nan),
        "statistic": result.get("stat", np.nan),
        "pvalue": result.get("pvalue", np.nan),
        "padj": result.get("padj", np.nan),
        "n_samples_group_1": int((np.asarray(conditions) == group_1).sum()),
        "n_samples_group_2": int((np.asarray(conditions) == group_2).sum()),
        "method": "pydeseq2",
        "design": dds.uns.get("design_formula", "~ condition"),
    })


def _make_deseq2_dataset(counts, conditions, batch_values=None, batch_key="batch"):
    """Construct one PyDESeq2 dataset with the same design used for DEG."""
    from pydeseq2.dds import DeseqDataSet

    integer_counts = np.asarray(np.rint(counts), dtype=np.int64)
    columns = [f"g{index}" for index in range(integer_counts.shape[1])]
    count_frame = pd.DataFrame(integer_counts, columns=columns)
    metadata = pd.DataFrame({"condition": pd.Categorical(conditions)})
    has_batch = batch_values is not None and len(np.unique(batch_values)) > 1
    if has_batch:
        metadata[batch_key] = pd.Categorical(np.asarray(batch_values).astype(str))
    design_factors = [batch_key, "condition"] if has_batch else ["condition"]
    design_formula = "~ " + " + ".join(design_factors)
    try:
        dds = DeseqDataSet(
            counts=count_frame, metadata=metadata, design_factors=design_factors,
            refit_cooks=True, quiet=True,
        )
    except TypeError:
        # PyDESeq2 >= 0.5 changed ``design_factors`` to a formula API.
        dds = DeseqDataSet(
            counts=count_frame, metadata=metadata, design=design_formula,
            refit_cooks=True, quiet=True,
        )
    dds.uns["design_formula"] = design_formula
    return dds


def _median_ratio_normalize(counts):
    """Return DESeq2 median-ratio normalized counts with a robust local fallback."""
    values = np.asarray(counts, dtype=float)
    try:
        from pydeseq2.preprocessing import deseq2_norm

        normalized, size_factors = deseq2_norm(values)
        normalized = np.asarray(normalized, dtype=float)
        size_factors = np.asarray(size_factors, dtype=float).reshape(-1)
        if (normalized.shape == values.shape and len(size_factors) == len(values)
                and np.isfinite(size_factors).all() and np.all(size_factors > 0)):
            return normalized, size_factors
    except Exception:
        pass

    # Sparse pseudobulk units may contain a zero in every gene.  The ordinary
    # median-ratio estimate is then undefined, so use the positive-count
    # version and retain an explicit transform label in the audit instead of
    # inventing missing sample-level QC values.
    with np.errstate(divide="ignore", invalid="ignore"):
        log_values = np.where(values > 0, np.log(values), np.nan)
        log_geomean = np.nanmean(log_values, axis=0)
        log_ratios = log_values - log_geomean[None, :]
        log_size_factors = np.nanmedian(log_ratios, axis=1)
    library_sizes = values.sum(axis=1)
    fallback = library_sizes / max(float(np.nanmedian(library_sizes[library_sizes > 0])), 1.0)
    log_size_factors = np.where(
        np.isfinite(log_size_factors), log_size_factors,
        np.log(np.maximum(fallback, 1e-8)),
    )
    size_factors = np.exp(log_size_factors)
    size_factors[~np.isfinite(size_factors) | (size_factors <= 0)] = 1.0
    positive = size_factors[size_factors > 0]
    if len(positive):
        size_factors /= np.exp(np.mean(np.log(positive)))
    return values / size_factors[:, None], size_factors


def _pseudobulk_vst(counts, conditions, batch_values=None, batch_key="batch"):
    """Fit a blind PyDESeq2 VST for sample QC, with an honest safe fallback.

    The QC transform intentionally uses ``use_design=False``: PCA and sample
    correlations must show the full donor-to-donor structure rather than a
    condition-adjusted representation.  The DESeq2 model used for inference is
    still fitted separately with its requested design.
    """
    try:
        dds = _make_deseq2_dataset(
            counts, conditions, batch_values=batch_values, batch_key=batch_key,
        )
        dds.vst(use_design=False)
        transformed = np.asarray(dds.layers["vst_counts"], dtype=float)
        size_factors = np.asarray(dds.obs["size_factors"], dtype=float).reshape(-1)
        if (transformed.shape != np.asarray(counts).shape or len(size_factors) != len(counts)
                or not np.isfinite(transformed).all() or not np.isfinite(size_factors).all()
                or np.any(size_factors <= 0)):
            raise RuntimeError("PyDESeq2 VST 返回了无效的变换矩阵或 size factor")
        return transformed, size_factors, "pydeseq2_blind_vst", ""
    except Exception as exc:
        normalized, size_factors = _median_ratio_normalize(counts)
        transformed = np.log2(np.maximum(normalized, 0.0) + 1.0)
        return (
            transformed,
            size_factors,
            "deseq2_median_ratio_log2_fallback",
            f"PyDESeq2 VST 未完成，QC 使用 DESeq2 median-ratio normalized log2 fallback: {exc}",
        )


def _sample_pca(values):
    """Small, dependency-light PCA for an n_samples × n_genes QC matrix."""
    matrix = np.nan_to_num(np.asarray(values, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    coordinates = np.zeros((matrix.shape[0], 2), dtype=float)
    explained = np.zeros(2, dtype=float)
    try:
        left, singular, _ = np.linalg.svd(centered, full_matrices=False)
        n_components = min(2, len(singular))
        coordinates[:, :n_components] = left[:, :n_components] * singular[:n_components]
        total = float(np.sum(singular ** 2))
        if total > 0:
            explained[:n_components] = singular[:n_components] ** 2 / total
    except np.linalg.LinAlgError:
        pass
    return coordinates, explained


def _sample_correlation(values):
    """Pearson sample × sample correlation from transformed pseudobulk counts."""
    matrix = np.nan_to_num(np.asarray(values, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = np.corrcoef(matrix)
    correlation = np.asarray(correlation, dtype=float)
    correlation[~np.isfinite(correlation)] = 0.0
    np.fill_diagonal(correlation, 1.0)
    return correlation


def _top_deg_positions(table, top_n):
    """Select visible genes without requiring an arbitrary significance cutoff."""
    ranked = table.copy()
    ranked["_padj_rank"] = pd.to_numeric(ranked.get("padj"), errors="coerce").fillna(1.0)
    ranked["_abs_log2fc_rank"] = pd.to_numeric(
        ranked.get("log2FC"), errors="coerce",
    ).abs().fillna(0.0)
    ranked = ranked.sort_values(
        ["_padj_rank", "_abs_log2fc_rank"], ascending=[True, False], kind="mergesort",
    )
    return ranked.head(max(1, int(top_n))).index.to_numpy(dtype=int)


def _requested_gene_labels(value):
    """Parse the one shared Volcano/MA gene-label control deterministically."""
    if value is None:
        return ()
    if isinstance(value, (list, tuple, set)):
        raw_values = value
    else:
        raw_values = str(value).replace("\n", ",").replace(";", ",").replace("；", ",").split(",")
    return tuple(dict.fromkeys(
        str(gene).strip() for gene in raw_values if str(gene).strip()
    ))


def _pseudobulk_library_size_figure(qc):
    """Show the three quantities that determine pseudobulk sample support."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    table = qc["sample_table"]
    labels = table["sample_id"].astype(str).tolist()
    conditions = table["condition"].astype(str).tolist()
    unique_conditions = list(dict.fromkeys(conditions))
    palette = ["#4C78A8", "#E45756", "#54A24B", "#B279A2", "#F2CF5B"]
    color_map = {
        condition: palette[index % len(palette)]
        for index, condition in enumerate(unique_conditions)
    }
    colors = [color_map[condition] for condition in conditions]
    x = np.arange(len(table))
    panels = (
        ("total_raw_counts", "Total raw counts"),
        ("n_cells", "Number of cells"),
        ("deseq2_size_factor", "DESeq2 size factor"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.7), dpi=180, sharex=True)
    for axis, (column, ylabel) in zip(axes, panels):
        values = pd.to_numeric(table[column], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        bars = axis.bar(x, values, color=colors, edgecolor="white", linewidth=0.45)
        axis.set_title(ylabel, loc="left", pad=4)
        axis.set_ylabel(ylabel)
        axis.set_xticks(x, labels, rotation=42, ha="right")
        axis.grid(axis="y", color="#D9DEE7", linewidth=0.55, alpha=0.8)
        axis.set_axisbelow(True)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
        if len(values) <= 8:
            offset = max(float(np.nanmax(values)) * 0.018, 0.02) if len(values) else 0.02
            for bar, value in zip(bars, values):
                if abs(value) >= 1_000_000:
                    value_label = f"{value / 1_000_000:.1f}M"
                elif abs(value) >= 1_000:
                    value_label = f"{value / 1_000:.1f}k"
                else:
                    value_label = f"{value:.3g}"
                axis.text(
                    bar.get_x() + bar.get_width() / 2, bar.get_height() + offset,
                    value_label, ha="center", va="bottom", fontsize=7.2,
                )
    fig.suptitle("Pseudobulk sample support", x=0.08, ha="left", fontsize=12, fontweight="semibold")
    fig.legend(
        handles=[Patch(facecolor=color_map[value], edgecolor="none", label=value)
                 for value in unique_conditions],
        loc="upper right", bbox_to_anchor=(0.985, 0.995), frameon=False,
        title="Condition", fontsize=8, title_fontsize=8,
    )
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.28, top=0.78, wspace=0.36)
    return fig


def _pseudobulk_qc_data(counts, meta):
    """Prepare one cluster/all-cells pseudobulk QC unit and its audit table."""
    batch_values = meta["batch"].astype(str).to_numpy() if "batch" in meta.columns else None
    transformed, size_factors, transform_name, warning = _pseudobulk_vst(
        counts, meta["condition"].astype(str).to_numpy(), batch_values=batch_values,
    )
    coordinates, explained = _sample_pca(transformed)
    sample_table = meta.copy()
    sample_table = sample_table.rename(columns={"library_size": "total_raw_counts"})
    sample_table["deseq2_size_factor"] = size_factors
    sample_table["expression_transform"] = transform_name
    desired_columns = [
        column for column in (
            "sample_id", "condition", "batch", "cluster", "total_raw_counts", "n_cells",
            "deseq2_size_factor", "expression_transform",
        ) if column in sample_table.columns
    ]
    return {
        "transformed": transformed,
        "coordinates": coordinates,
        "explained_variance": explained,
        "correlation": _sample_correlation(transformed),
        "sample_table": sample_table.loc[:, desired_columns],
        "transform": transform_name,
        "warning": warning,
    }


class SCBatchCSVExport(BaseAnalysis):
    """Export a stable, documented CSV package without modifying the h5ad."""

    MODULE_NAME = "sc_csv_export"
    DISPLAY_NAME = "批量单细胞 CSV 导出"
    DESCRIPTION = "导出细胞元数据、样本组成、pseudobulk 表达和比较注册表"

    def run(self, input_path):
        adata = self.load_adata(input_path)
        sample_key = str(self.params.get("sample_key", "sample_id") or "").strip()
        condition_key = str(self.params.get("condition_key", "condition") or "").strip()
        cluster_key = str(
            self.params.get("cluster_key", self.params.get("celltype_key", "leiden")) or ""
        ).strip()
        min_samples = max(2, int(self.params.get("min_samples_per_group", 2)))
        comparisons = str(self.params.get("comparisons", "") or "")
        comparison_mode = str(self.params.get("comparison_mode", "all_pairwise") or "all_pairwise")
        if comparison_mode not in {"all_pairwise", "vs_reference"}:
            raise ValueError("comparison_mode 必须为 all_pairwise 或 vs_reference")
        reference_group = str(self.params.get("reference_group", "") or "")
        prefix = str(self.params.get("export_prefix", "sc_batch") or "sc_batch").strip()
        prefix = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in prefix).strip("._-") or "sc_batch"
        export_folder = str(self.params.get("export_folder", "sc_batch_results") or "sc_batch_results")
        export_scope = str(self.params.get("export_scope", "full") or "full")
        if export_scope not in {"full", "proportions_expression"}:
            raise ValueError("export_scope 必须为 full 或 proportions_expression")
        included_keys = (
            ("proportions", "expression")
            if export_scope == "proportions_expression" else tuple(CSV_PACKAGE_FOLDERS)
        )

        self.progress(8, "正在验证样本级设计...")
        design = sample_design_from_obs(adata, sample_key, condition_key)
        package_root, package_dirs = batch_csv_package_dirs(
            self.project_dir, export_folder, included_keys=included_keys,
        )
        result_files = []

        # One cell per row preserves a future platform's option to re-aggregate
        # by user-selected metadata without requiring the original h5ad.
        self.progress(18, "正在导出细胞和样本元数据...")
        cell_metadata = adata.obs.copy()
        cell_metadata.insert(0, "cell_id", adata.obs_names.astype(str))
        # UMAP coordinates are cell-level data, so keep them beside the
        # corresponding cell_id rather than attempting to aggregate them into
        # the sample-level tables.  Preserve every available dimension for
        # nonstandard embeddings while maintaining the familiar UMAP_1/2 names.
        if "X_umap" in adata.obsm:
            coordinates = np.asarray(adata.obsm["X_umap"])
            if coordinates.ndim == 2 and coordinates.shape[0] == adata.n_obs:
                for index in range(coordinates.shape[1]):
                    cell_metadata[f"UMAP_{index + 1}"] = coordinates[:, index]
        cell_path = _available_path(package_dirs["proportions"], f"{prefix}_cell_metadata", ".csv")
        cell_metadata.to_csv(cell_path, index=False)
        result_files.append({"file_path": cell_path, "file_type": "csv", "category": "table",
                             "label": "Single-cell metadata with sample provenance"})

        sample_path = _available_path(package_dirs["proportions"], f"{prefix}_sample_metadata", ".csv")
        design.to_csv(sample_path, index=False)
        result_files.append({"file_path": sample_path, "file_type": "csv", "category": "table",
                             "label": "Biological sample design"})

        self.progress(35, "正在聚合样本级原始 counts...")
        pseudo_meta, pseudo_counts = aggregate_pseudobulk_counts(
            adata, sample_key=sample_key, condition_key=condition_key, min_cells=1,
        )
        gene_id, gene_name = _gene_metadata(adata)
        # Match the DOCX gene-expression convention (one gene per row and one
        # sample per column) while retaining the stable source identifier.
        wide_counts = pd.DataFrame({"GeneID": gene_id, "GeneName": gene_name})
        wide_logcpm = pd.DataFrame({"GeneID": gene_id, "GeneName": gene_name})
        for row_index, (_, row) in enumerate(pseudo_meta.iterrows()):
            sample_id = str(row["sample_id"])
            wide_counts[sample_id] = np.rint(pseudo_counts[row_index]).astype(np.int64)
            wide_logcpm[sample_id] = _logcpm(pseudo_counts[row_index:row_index + 1])[0]
        count_path = _available_path(package_dirs["expression"], f"{prefix}_pseudobulk_counts", ".csv")
        logcpm_path = _available_path(package_dirs["expression"], f"{prefix}_pseudobulk_log2cpm", ".csv")
        wide_counts.to_csv(count_path, index=False)
        wide_logcpm.to_csv(logcpm_path, index=False)
        result_files.extend([
            {"file_path": count_path, "file_type": "csv", "category": "table",
             "label": "Sample-level pseudobulk raw counts (gene × sample)"},
            {"file_path": logcpm_path, "file_type": "csv", "category": "table",
             "label": "Sample-level pseudobulk log2(CPM+1) (gene × sample)"},
        ])

        if cluster_key and cluster_key in adata.obs.columns:
            self.progress(52, "正在导出每样本聚类组成...")
            composition = pd.DataFrame({
                "sample_id": adata.obs[sample_key].astype(str),
                "condition": adata.obs[condition_key].astype(str),
                "cluster": adata.obs[cluster_key].astype(str),
            })
            composition = composition.groupby(
                ["sample_id", "condition", "cluster"], observed=True
            ).size().rename("cell_count").reset_index()
            totals = composition.groupby("sample_id", observed=True)["cell_count"].sum().rename("total_cells")
            composition = composition.join(totals, on="sample_id")
            composition["proportion"] = composition["cell_count"] / composition["total_cells"].clip(lower=1)
            composition_path = _available_path(package_dirs["proportions"], f"{prefix}_sample_cluster_proportions", ".csv")
            composition.to_csv(composition_path, index=False)
            result_files.append({"file_path": composition_path, "file_type": "csv", "category": "table",
                                 "label": "Sample-level cluster proportions"})
            # The downstream DOCX example is wide: one row per sample, one
            # component/cluster per column.  Keep the long table above for
            # traceability and add this directly plottable representation.
            proportion_wide = composition.pivot_table(
                index=["sample_id", "condition"], columns="cluster", values="proportion",
                aggfunc="first", fill_value=0,
            ).reset_index()
            proportion_wide.columns = [
                "SampleID" if column == "sample_id" else
                "SampleGroup" if column == "condition" else
                f"Cluster_{column}"
                for column in proportion_wide.columns
            ]
            if "experiment_id" in adata.obs.columns:
                experiment_by_sample = pd.DataFrame({
                    "SampleID": adata.obs[sample_key].astype(str),
                    "ExperimentID": adata.obs["experiment_id"].astype(str),
                }).drop_duplicates("SampleID")
                proportion_wide = proportion_wide.merge(experiment_by_sample, on="SampleID", how="left")
                first_columns = ["SampleID", "SampleGroup", "ExperimentID"]
                proportion_wide = proportion_wide[
                    [column for column in first_columns if column in proportion_wide.columns]
                    + [column for column in proportion_wide.columns if column not in first_columns]
                ]
            wide_composition_path = _available_path(package_dirs["proportions"], f"{prefix}_group_proportions", ".csv")
            proportion_wide.to_csv(wide_composition_path, index=False)
            result_files.append({"file_path": wide_composition_path, "file_type": "csv", "category": "table",
                                 "label": "DOCX-compatible group proportions (wide)"})

        comparison_plan = pd.DataFrame()
        if export_scope == "full":
            self.progress(70, "正在生成可运行比较注册表...")
            comparison_plan = build_comparison_plan(
                design, requested_comparisons=comparisons, mode=comparison_mode,
                reference_group=reference_group, min_samples_per_group=min_samples,
            )
            comparison_path = _available_path(package_dirs["deg"], f"{prefix}_comparison_manifest", ".csv")
            comparison_plan.to_csv(comparison_path, index=False)
            result_files.append({"file_path": comparison_path, "file_type": "csv", "category": "table",
                                 "label": "Pseudobulk comparison registry"})

            # GO is intentionally never fabricated from a DEG table: record
            # whether it is pending or blocked rather than mistaking an empty
            # directory for a completed enrichment analysis.
            go_status = comparison_plan.copy()
            if go_status.empty:
                go_status = pd.DataFrame(columns=["comparison_id", "comparison", "go_status", "reason"])
            else:
                go_status["go_status"] = np.where(
                    go_status["status"].eq("runnable"), "not_run", "not_runnable"
                )
                go_status["reason"] = np.where(
                    go_status["status"].eq("runnable"),
                    "尚未运行 GO 富集；需先使用相同 comparison_id 的有效 DEG CSV。",
                    go_status["reason"],
                )
            go_status_path = _available_path(package_dirs["go"], f"{prefix}_go_enrichment_status", ".csv")
            go_status.to_csv(go_status_path, index=False)
            result_files.append({"file_path": go_status_path, "file_type": "csv", "category": "table",
                                 "label": "GO enrichment status by comparison"})
        else:
            self.progress(70, "已按要求跳过 DEG、GO 与比较注册表导出...")

        metadata = {
            "schema_version": "sc-batch-csv-v1",
            "input_h5ad": os.path.basename(input_path),
            "sample_key": sample_key,
            "condition_key": condition_key,
            "cluster_key": cluster_key if cluster_key in adata.obs.columns else "",
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_samples": int(len(design)),
            "csv_package_dir": package_root,
            "csv_package_folders": {
                key: CSV_PACKAGE_FOLDERS[key] for key in included_keys
            },
            "export_scope": export_scope,
            "pseudobulk_count_source": "layers[counts]",
            "notes": [
                "cell_metadata.csv 可供下游按任意注释字段重新汇总。",
                "pseudobulk_counts.csv 为所有细胞的样本级原始 counts。",
            ],
        }
        if export_scope == "full":
            metadata["comparison_inference_unit"] = "biological_sample_pseudobulk"
            metadata["notes"].extend([
                "簇特异 DEG 由 sc_pseudobulk_deg 单独导出。",
                "comparison_manifest 中 not_runnable 的比较不应报告样本级统计推断。",
            ])
        else:
            metadata["notes"].append("本结果包只包含 01_group_proportions 与 02_gene_expression。")
        metadata_path = _available_path(package_dirs["proportions"], f"{prefix}_export_metadata", ".json")
        with open(metadata_path, "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, ensure_ascii=False, indent=2)
        result_files.append({"file_path": metadata_path, "file_type": "json", "category": "info",
                             "label": "Batch CSV export contract metadata"})
        self.progress(100, "批量单细胞 CSV 导出完成")
        summary = {
            "schema_version": metadata["schema_version"],
            "n_cells": int(adata.n_obs), "n_samples": int(len(design)),
            "n_comparisons": int(len(comparison_plan)),
            "csv_package_dir": package_root,
            "n_runnable_comparisons": int((comparison_plan.get("status", pd.Series(dtype=str)) == "runnable").sum()),
            "cluster_proportions_exported": bool(cluster_key and cluster_key in adata.obs.columns),
            "export_scope": export_scope,
        }
        return {
            "output_adata": input_path,
            "result_files": result_files,
            "summary": summary,
        }


class SCPseudobulkDEG(BaseAnalysis):
    """Run all requested condition contrasts on sample-level pseudobulk counts."""

    MODULE_NAME = "sc_pseudobulk_deg"
    DISPLAY_NAME = "样本级 pseudobulk 差异表达"
    DESCRIPTION = "以生物学样本为统计单位的条件 DEG；可导出每个比较的完整基因 CSV"
    # pseudobulk 以 celltype/cluster × sample 为单元，输入必须已经带分组列。
    INPUT_REQUIRES = ['leiden']

    def run(self, input_path):
        adata = self.load_adata(input_path)
        adata = self.apply_scope(adata)
        sample_key = str(self.params.get("sample_key", "sample_id") or "").strip()
        condition_key = str(self.params.get("condition_key", "condition") or "").strip()
        batch_key = str(self.params.get("batch_key", "") or "").strip()
        scope = str(self.params.get("analysis_scope", "per_cluster") or "per_cluster")
        if scope not in {"all_cells", "per_cluster"}:
            raise ValueError("analysis_scope 必须为 all_cells 或 per_cluster")
        if scope == "per_cluster":
            grouping = resolve_cell_grouping(self.params, adata.obs.columns)
            cluster_key = grouping["key"]
        else:
            grouping = {
                "mode": "all_cells", "key": "", "label": "全部细胞",
                "value_label": "全部细胞",
            }
            cluster_key = ""
        min_cells = max(1, int(self.params.get("min_cells_per_sample_celltype", 20)))
        min_samples = max(2, int(self.params.get("min_samples_per_group", 2)))
        requested_method = str(self.params.get("method", "deseq2") or "deseq2")
        if requested_method not in {"deseq2", "welch_logcpm"}:
            raise ValueError(
                f"不支持的统计方法 '{requested_method}'；仅支持 deseq2 与 welch_logcpm。"
            )
        comparison_mode = str(self.params.get("comparison_mode", "all_pairwise") or "all_pairwise")
        if comparison_mode not in {"all_pairwise", "vs_reference"}:
            raise ValueError("comparison_mode 必须为 all_pairwise 或 vs_reference")
        reference_group = str(self.params.get("reference_group", "") or "")
        requested_comparisons = str(self.params.get("comparisons", "") or "")
        prefix = str(self.params.get("export_prefix", "sc_pseudobulk") or "sc_pseudobulk")
        prefix = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in prefix).strip("._-") or "sc_pseudobulk"
        export_folder = str(self.params.get("export_folder", "sc_batch_results") or "sc_batch_results")
        # Gene-level tables are registered as task outputs by default so the
        # result page can offer a direct CSV download for every contrast.  The
        # opt-out remains useful for compact, internal-only pipeline runs.
        export_full_tables = _as_bool(self.params.get("export_full_tables"), default=True)
        analysis_id = _safe_folder_name(self.params.get("_analysis_id"), fallback="pseudobulk_deg")
        show_deg_figures = self.params.get("show_deg_figures", True)
        if isinstance(show_deg_figures, str):
            show_deg_figures = show_deg_figures.strip().lower() not in {"", "0", "false", "no", "off"}
        plot_cluster_limit = int(self.params.get("plot_max_clusters", 0))
        if plot_cluster_limit < 0:
            raise ValueError("plot_max_clusters 必须为 0 或正整数")
        volcano_top_n = min(20, max(0, int(self.params.get("volcano_top_n", 8))))
        requested_label_genes = _requested_gene_labels(
            self.params.get("label_genes") or self.params.get("volcano_label_genes", "")
        )
        pseudobulk_heatmap_top_n = min(
            50, max(4, int(self.params.get("pseudobulk_heatmap_top_n", 20))),
        )
        padj_cutoff = float(self.params.get("padj_cutoff", 0.05))
        log2fc_cutoff = float(self.params.get("log2fc_cutoff", 1.0))
        if not np.isfinite(padj_cutoff) or not 0 <= padj_cutoff <= 1:
            raise ValueError("padj_cutoff 必须位于 [0, 1]")
        if not np.isfinite(log2fc_cutoff) or log2fc_cutoff < 0:
            raise ValueError("log2fc_cutoff 必须是非负有限数值")

        self.progress(8, "正在验证样本与条件设计...")
        design = sample_design_from_obs(adata, sample_key, condition_key)
        plan = build_comparison_plan(
            design, requested_comparisons=requested_comparisons,
            mode=comparison_mode, reference_group=reference_group,
            min_samples_per_group=min_samples,
        )
        self.progress(20, "正在按样本聚合原始 counts...")
        pseudo_meta, pseudo_counts = aggregate_pseudobulk_counts(
            adata, sample_key=sample_key, condition_key=condition_key,
            celltype_key=cluster_key, min_cells=min_cells, batch_key=batch_key,
        )
        pseudo_meta = pseudo_meta.rename(columns={"celltype": "cluster"})
        gene_ids, gene_names = _gene_metadata(adata)
        package_root, package_dirs = batch_csv_package_dirs(
            self.project_dir, export_folder, included_keys=("deg",),
        )
        internal_dir = os.path.join(package_dirs["deg"], ".internal")
        os.makedirs(internal_dir, exist_ok=True)
        warnings = []
        unit_statuses = []
        all_results = []
        internal_sources = []
        result_files = []
        # These sample-level diagnostics are mandatory evidence for a
        # pseudobulk inference result.  ``show_deg_figures`` only controls the
        # optional Volcano/MA pair; it must not hide donor outliers, library
        # support or a DEG driven by a single sample.
        plots_dir = self.ensure_plots_dir()
        from figure_engine import NatureFigureDirector

        director = NatureFigureDirector()
        pseudobulk_qc_by_cluster = {}
        pseudobulk_qc_units = []

        runnable_indices = plan.index[plan["status"] == "runnable"].tolist() if not plan.empty else []
        for iteration, plan_index in enumerate(runnable_indices, start=1):
            item = plan.loc[plan_index]
            group_1, group_2 = str(item["group_1"]), str(item["group_2"])
            comparison = str(item["comparison"])
            self.progress(25 + int((iteration - 1) * 60 / max(len(runnable_indices), 1)),
                          f"正在运行 {comparison} 的样本级 pseudobulk DEG...")
            comparison_rows = []
            for cluster in sorted(pseudo_meta["cluster"].unique().tolist()):
                unit_mask = pseudo_meta["cluster"].eq(cluster)
                meta = pseudo_meta.loc[unit_mask].reset_index(drop=True)
                counts = pseudo_counts[np.asarray(unit_mask)]
                group_1_mask = meta["condition"].eq(group_1).to_numpy()
                group_2_mask = meta["condition"].eq(group_2).to_numpy()
                if int(group_1_mask.sum()) < min_samples or int(group_2_mask.sum()) < min_samples:
                    unit_statuses.append({
                        "comparison_id": item["comparison_id"], "comparison": comparison,
                        "cluster": str(cluster), "status": "not_runnable",
                        "reason": f"每组至少需要 {min_samples} 个样本×{grouping['label']} 单元",
                    })
                    continue
                qc_key = str(cluster)
                if qc_key not in pseudobulk_qc_by_cluster:
                    self.progress(-1, f"正在生成 cluster {cluster} 的 pseudobulk 样本 QC...")
                    try:
                        qc = _pseudobulk_qc_data(counts, meta)
                        pseudobulk_qc_by_cluster[qc_key] = qc
                        safe_cluster = _safe_folder_name(cluster, fallback="All")
                        plot_context = (
                            "all cells" if not cluster_key
                            else f"{grouping['label']} {cluster}"
                        )
                        qc_table_path = _available_path(
                            package_dirs["deg"],
                            f"{prefix}_task_{analysis_id}_sample_qc_{safe_cluster}", ".csv",
                        )
                        qc["sample_table"].to_csv(qc_table_path, index=False)
                        result_files.append({
                            "file_path": qc_table_path, "file_type": "csv", "category": "qc",
                            "label": f"Pseudobulk sample QC: {plot_context}",
                        })
                        pseudobulk_qc_units.append({
                            "cluster": str(cluster),
                            "sample_qc_file": qc_table_path,
                            "n_samples": int(len(qc["sample_table"])),
                            "transform": qc["transform"],
                        })
                        if qc["warning"]:
                            warnings.append(f"cluster {cluster}: {qc['warning']}")
                        try:
                            pca_spec = director.spec_from_params(
                                "pca", self.params, width="single",
                                title=f"Pseudobulk sample PCA · {plot_context}",
                                evidence_role="quality_control", show_sample_labels=True,
                                show_legend=True, confidence_ellipse=False,
                            )
                            fig_pca = director.render(pca_spec, {
                                "coordinates": qc["coordinates"],
                                "groups": qc["sample_table"]["condition"].astype(str).tolist(),
                                "samples": qc["sample_table"]["sample_id"].astype(str).tolist(),
                                "explained_variance": qc["explained_variance"],
                            })
                            result_files.extend(self.save_matplotlib_figure(
                                fig_pca, plots_dir,
                                f"sc_pseudobulk_sample_pca_{prefix}_{safe_cluster}.png",
                                "pca", f"Pseudobulk sample PCA: {plot_context}",
                            ))

                            correlation_spec = director.spec_from_params(
                                "correlation_heatmap", self.params, width="single",
                                title=f"Pseudobulk sample correlation · {plot_context}",
                                evidence_role="quality_control", correlation_method="pearson",
                                col_cluster=False, mask_diagonal=False, annotate_cells=True,
                                max_col_labels=12,
                            )
                            fig_correlation = director.render(correlation_spec, {
                                "correlation_matrix": qc["correlation"],
                                "sample_labels": qc["sample_table"]["sample_id"].astype(str).tolist(),
                                "groups": qc["sample_table"]["condition"].astype(str).tolist(),
                            })
                            result_files.extend(self.save_matplotlib_figure(
                                fig_correlation, plots_dir,
                                f"sc_pseudobulk_sample_correlation_{prefix}_{safe_cluster}.png",
                                "correlation", f"Pseudobulk sample correlation: {plot_context}",
                            ))

                            fig_library = _pseudobulk_library_size_figure(qc)
                            result_files.extend(self.save_matplotlib_figure(
                                fig_library, plots_dir,
                                f"sc_pseudobulk_library_size_{prefix}_{safe_cluster}.png",
                                "qc", f"Pseudobulk library size: {plot_context}",
                            ))
                        except Exception as exc:
                            warnings.append(
                                f"cluster {cluster}: pseudobulk 样本 QC 绘图失败（{exc}）"
                            )
                    except Exception as exc:
                        warnings.append(
                            f"cluster {cluster}: 无法生成 pseudobulk 样本 QC（{exc}）"
                        )
                if requested_method == "deseq2":
                    identifiable, design_reason = _pseudobulk_design_is_identifiable(meta)
                    if not identifiable:
                        unit_statuses.append({
                            "comparison_id": item["comparison_id"], "comparison": comparison,
                            "cluster": str(cluster), "status": "not_runnable",
                            "reason": design_reason,
                        })
                        warnings.append(f"{comparison} / cluster {cluster}: {design_reason}")
                        continue
                    try:
                        batch_values = (
                            meta["batch"].astype(str).to_numpy()
                            if "batch" in meta.columns else None
                        )
                        table = _deseq2_deg(
                            counts, meta["condition"].astype(str).to_numpy(), group_1, group_2,
                            gene_ids, gene_names, comparison, str(cluster),
                            batch_values=batch_values, batch_key="batch",
                        )
                    except Exception as exc:
                        # Do not silently swap a negative-binomial model for a
                        # different test.  A failed formal model is a failed
                        # unit, not evidence from an unannounced fallback.
                        reason = f"PyDESeq2 拟合失败，未回退到 Welch: {exc}"
                        warnings.append(
                            f"{comparison} / cluster {cluster}: {reason}"
                        )
                        unit_statuses.append({
                            "comparison_id": item["comparison_id"], "comparison": comparison,
                            "cluster": str(cluster), "status": "model_failed", "reason": reason,
                        })
                        continue
                else:
                    table = _welch_logcpm_deg(
                        counts, group_1_mask, group_2_mask, gene_ids, gene_names,
                        comparison, group_1, group_2, str(cluster), min_samples,
                    )
                    if batch_key:
                        warnings.append(
                            f"{comparison} / cluster {cluster}: 当前选择 Welch，未建模 batch_key={batch_key}。"
                        )
                if not table.empty:
                    if "design" not in table.columns:
                        table["design"] = (
                            "~ condition (batch not modeled)" if batch_key
                            else "~ condition"
                        )
                    table.insert(0, "comparison_id", item["comparison_id"])
                    # Canonical direction aliases shared with the cell-level
                    # contract and downstream enrichment.  Positive log2FC is
                    # always experimental/group_1 over control/group_2.
                    table["experimental_group"] = group_1
                    table["control_group"] = group_2
                    table["comparison_type"] = "condition"
                    table["deg_scope"] = "per_cluster" if cluster_key else "all_cells"
                    table["analysis_grouping"] = grouping["mode"]
                    table["analysis_group_key"] = cluster_key
                    table["analysis_group_label"] = grouping["label"]
                    table["inference_unit"] = "biological_sample"
                    table["statistical_status"] = (
                        "biological_replicate_model"
                        if requested_method == "deseq2" else "sample_level_approximate"
                    )
                    table["significant"] = (table["padj"] < padj_cutoff) & (
                        np.abs(table["log2FC"]) >= log2fc_cutoff
                    )
                    # The expression heatmap shares the exact sample ordering
                    # with PCA/correlation and uses the blind VST (or an
                    # explicitly labelled median-ratio log2 fallback).  This
                    # makes a one-donor DEG immediately visible rather than
                    # allowing the group statistic to conceal it.
                    qc = pseudobulk_qc_by_cluster.get(str(cluster))
                    if qc is not None:
                        try:
                            top_positions = _top_deg_positions(table, pseudobulk_heatmap_top_n)
                            top_positions = top_positions[
                                (top_positions >= 0) & (top_positions < qc["transformed"].shape[1])
                            ]
                            if len(top_positions):
                                heatmap_table = table.loc[top_positions]
                                gene_labels = heatmap_table["gene"].astype(str).tolist()
                                if len(set(gene_labels)) != len(gene_labels):
                                    gene_labels = [
                                        f"{gene} ({gene_id})"
                                        for gene, gene_id in zip(
                                            heatmap_table["gene"].astype(str),
                                            heatmap_table["gene_id"].astype(str),
                                        )
                                    ]
                                safe_cluster = _safe_folder_name(cluster, fallback="All")
                                plot_context = (
                                    "all cells" if not cluster_key
                                    else f"{grouping['label']} {cluster}"
                                )
                                heatmap_spec = director.spec_from_params(
                                    "heatmap", self.params, width="single",
                                    title=(f"Top DEG expression · {comparison}"
                                           + (f" · {cluster}" if cluster_key else "")),
                                    evidence_role="quality_control", zscore="row",
                                    row_cluster=True, col_cluster=False,
                                    max_row_labels=pseudobulk_heatmap_top_n,
                                    max_col_labels=12,
                                )
                                fig_heatmap = director.render(heatmap_spec, {
                                    # ``transformed`` is samples × genes; the
                                    # heatmap contract is genes × samples.
                                    "matrix": qc["transformed"][:, top_positions].T,
                                    "gene_labels": gene_labels,
                                    "sample_labels": qc["sample_table"]["sample_id"].astype(str).tolist(),
                                    "annotations": {
                                        "Condition": qc["sample_table"]["condition"].astype(str).tolist(),
                                    },
                                    "colorbar_label": "Gene-wise z-score of VST expression",
                                })
                                result_files.extend(self.save_matplotlib_figure(
                                    fig_heatmap, plots_dir,
                                    (f"sc_pseudobulk_top_deg_expression_{prefix}_"
                                     f"{item['comparison_id']}_{safe_cluster}.png"),
                                    "heatmap",
                                    f"Pseudobulk Top DEG expression: {comparison} · {plot_context}",
                                ))
                        except Exception as exc:
                            warnings.append(
                                f"{comparison} / cluster {cluster}: Top DEG 样本表达热图失败（{exc}）"
                            )
                    comparison_rows.append(table)
                    unit_statuses.append({
                        "comparison_id": item["comparison_id"], "comparison": comparison,
                        "cluster": str(cluster), "status": "completed",
                        "reason": "" if requested_method == "deseq2" else "Welch log2CPM 为显式选择的近似样本级分析。",
                    })
            if comparison_rows:
                frame = pd.concat(comparison_rows, ignore_index=True)
                internal_path = _available_path(
                    internal_dir, f"{prefix}_task_{analysis_id}_deg_{item['comparison_id']}", ".csv",
                )
                frame.to_csv(internal_path, index=False)
                internal_sources.append(internal_path)
                if export_full_tables:
                    output_path = _available_path(
                        package_dirs["deg"], f"{prefix}_deg_{item['comparison_id']}", ".csv"
                    )
                    frame.to_csv(output_path, index=False)
                    plan.loc[plan_index, "output_file"] = os.path.basename(output_path)
                    result_files.append({"file_path": output_path, "file_type": "csv", "category": "table",
                                         "label": f"Pseudobulk DEG: {comparison}"})
                plan.loc[plan_index, "n_rows"] = int(len(frame))
                plan.loc[plan_index, "n_significant"] = int(frame["significant"].sum())
                if show_deg_figures:
                    cluster_priority = (
                        frame.groupby("cluster", observed=True)["significant"]
                        .sum().sort_values(ascending=False).index.astype(str).tolist()
                    )
                    clusters_to_plot = (
                        cluster_priority if plot_cluster_limit == 0
                        else cluster_priority[:plot_cluster_limit]
                    )
                    if len(cluster_priority) > len(clusters_to_plot):
                        warnings.append(
                            f"{comparison}: pseudobulk 图仅展示显著 DEG 最多的前 "
                            f"{plot_cluster_limit} 个 cluster；完整结果保留在受控内部表。"
                        )
                    for cluster in clusters_to_plot:
                        view = frame.loc[frame["cluster"].astype(str) == cluster].copy()
                        if view.empty:
                            continue
                        safe_cluster = _safe_folder_name(cluster, fallback="All")
                        figure_stem = (
                            f"{prefix}_{item['comparison_id']}_{safe_cluster}"
                        )
                        title_suffix = (
                            comparison if not cluster_key
                            else f"{comparison} · {grouping['label']} {cluster}"
                        )
                        plot_context = (
                            "all cells" if not cluster_key
                            else f"{grouping['label']} {cluster}"
                        )
                        try:
                            volcano_spec = director.spec_from_params(
                                "volcano", self.params, width="single",
                                title=f"{comparison}\nPseudobulk Volcano · {plot_context}",
                                evidence_role="comparison", fc_threshold=log2fc_cutoff,
                                fdr_threshold=padj_cutoff, label_n=volcano_top_n,
                                label_genes=requested_label_genes, show_legend=True,
                            )
                            fig_volcano = director.render(volcano_spec, view)
                            result_files.extend(self.save_matplotlib_figure(
                                fig_volcano, plots_dir,
                                f"sc_pseudobulk_volcano_{figure_stem}.png",
                                "volcano", f"Pseudobulk Volcano: {title_suffix}",
                            ))

                            ma_view = view.rename(columns={"base_mean_count": "mean_expression"})
                            ma_spec = director.spec_from_params(
                                "ma", self.params, width="single",
                                title=f"{comparison}\nPseudobulk MA · {plot_context}",
                                evidence_role="comparison", fc_threshold=log2fc_cutoff,
                                fdr_threshold=padj_cutoff, label_n=min(volcano_top_n, 6),
                                label_genes=requested_label_genes, show_legend=False,
                            )
                            fig_ma = director.render(ma_spec, ma_view)
                            result_files.extend(self.save_matplotlib_figure(
                                fig_ma, plots_dir,
                                f"sc_pseudobulk_ma_{figure_stem}.png",
                                "ma", f"Pseudobulk MA: {title_suffix}",
                            ))
                        except Exception as exc:
                            warnings.append(
                                f"{comparison} / {grouping['label']} {cluster}: pseudobulk DEG 绘图失败（{exc}）"
                            )
                all_results.append(frame)
            else:
                plan.loc[plan_index, "status"] = "not_runnable"
                plan.loc[plan_index, "reason"] = (
                    "过滤后没有满足每组最少样本数的 sample × "
                    f"{grouping['label']} 单元"
                )

        self.progress(90, "正在保存比较审计与结果摘要...")
        internal_table = internal_sources[0] if internal_sources else ""
        if all_results:
            combined = pd.concat(all_results, ignore_index=True)
            if export_full_tables:
                plan_path = _available_path(package_dirs["deg"], f"{prefix}_comparison_manifest", ".csv")
                plan.to_csv(plan_path, index=False)
                result_files.append({"file_path": plan_path, "file_type": "csv", "category": "table",
                                     "label": "Pseudobulk DEG comparison registry"})
                combined_path = _available_path(package_dirs["deg"], f"{prefix}_deg_all_comparisons", ".csv")
                combined.to_csv(combined_path, index=False)
                result_files.append({"file_path": combined_path, "file_type": "csv", "category": "table",
                                     "label": "Pseudobulk DEG combined across comparisons"})
        if all_results:
            statistical_status = (
                "biological_replicate_model" if requested_method == "deseq2"
                else "sample_level_approximate"
            )
        elif unit_statuses:
            statuses = {str(row.get("status", "")) for row in unit_statuses}
            statistical_status = (
                "model_failed" if "model_failed" in statuses
                else "no_runnable_units"
            )
        else:
            statistical_status = "no_comparison_plan"

        audit = {
            "analysis_level": "sample_level_pseudobulk",
            "inference_unit": "biological_sample",
            "method_requested": requested_method,
            "statistical_status": statistical_status,
            "design": {"sample_key": sample_key, "condition_key": condition_key,
                       "cluster_key": cluster_key, "batch_key": batch_key,
                       "analysis_grouping": grouping["mode"],
                       "analysis_group_label": grouping["label"]},
            "comparison_plan": plan.to_dict(orient="records"),
            "unit_statuses": unit_statuses,
            "sample_qc": {
                "mandatory": True,
                "units": pseudobulk_qc_units,
                "pca_coloring": "condition",
                "correlation": "Pearson on blind VST / explicit median-ratio log2 fallback",
                "top_deg_heatmap": {
                    "top_n": pseudobulk_heatmap_top_n,
                    "sample_order": "pseudobulk aggregation order; not clustered",
                    "display": "gene-wise z-score of the transformed pseudobulk expression",
                },
            },
            "warnings": warnings,
            "full_table_exported": export_full_tables,
        }
        audit_path = _available_path(internal_dir, f"{prefix}_task_{analysis_id}_audit", ".json")
        Path(audit_path).write_text(json.dumps(audit, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
        result_files.append({"file_path": audit_path, "file_type": "json", "category": "qc",
                             "label": "Pseudobulk DEG statistical audit"})
        self.progress(100, "样本级 pseudobulk DEG 完成")
        summary = {
            "n_comparisons": int(len(plan)),
            "n_completed_comparisons": int(len({
                str(frame["comparison_id"].iloc[0]) for frame in all_results if not frame.empty
            })),
            "analysis_scope": "per_cluster" if cluster_key else "all_cells",
            "sample_key": sample_key, "condition_key": condition_key,
            "scope_key": str(self.params.get("scope_key", "") or "").strip() or None,
            "scope_values": ([v.strip() for v in str(self.params.get("scope_values", "") or "").split(",") if v.strip()] or None),
            "cluster_key": cluster_key,
            "analysis_grouping": grouping["mode"],
            "analysis_group_key": cluster_key,
            "analysis_group_label": grouping["label"],
            "batch_key": batch_key,
            "csv_package_dir": package_root,
            "method_requested": requested_method,
            "statistical_status": audit["statistical_status"],
            "analysis_status": (
                "completed_with_warning"
                if audit["statistical_status"] in {"no_runnable_units", "model_failed", "no_comparison_plan"}
                else "completed"
            ),
            "n_model_failed_units": int(sum(row["status"] == "model_failed" for row in unit_statuses)),
            "n_pseudobulk_qc_units": int(len(pseudobulk_qc_units)),
            "pseudobulk_qc_files": [row["sample_qc_file"] for row in pseudobulk_qc_units],
            "pseudobulk_qc_mandatory": True,
            "deg_source_file": internal_table,
            "deg_source_files": internal_sources,
            "deg_source_level": "pseudobulk",
            "deg_source_task_contract": "condition",
            "warnings": warnings[:20],
        }
        if unit_statuses and not all_results and any(
            row.get("status") == "model_failed" for row in unit_statuses
        ):
            summary["error"] = "所有可运行 pseudobulk 单元的 DESeq2 模型均拟合失败；未自动回退。"
        return {
            "output_adata": input_path,
            "result_files": result_files,
            "summary": summary,
        }
