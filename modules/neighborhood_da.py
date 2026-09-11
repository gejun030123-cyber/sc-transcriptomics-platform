"""Sample-level differential abundance for local scRNA-seq neighbourhoods.

This is intentionally a small, transparent neighbourhood analysis for the
platform's usual 1--2 concurrent laboratory runs.  It follows the important
Milo-style contract -- local, overlapping KNN neighbourhoods and *samples*
as the independent observations -- without presenting itself as a replacement
for the full R/Milo negative-binomial implementation.

For every selected seed cell, the module forms a KNN neighbourhood in an
integrated/PCA representation.  Its abundance in a biological sample is the
fraction of that sample's scoped cells falling in the neighbourhood.  A
condition comparison is then made across samples, with BH correction across
neighbourhoods for each requested comparison.  Cell counts are never used as
replicates.
"""

from __future__ import annotations

import itertools
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from modules.base import BaseAnalysis


_TECHNICAL_BATCH_NAMES = {
    "batch", "technical_batch", "sequencing_batch", "library_batch",
    "lane", "run", "sequencing_run",
}


def _as_bool(value, default=False):
    if value is None:
        return bool(default)
    if isinstance(value, str):
        return value.strip().lower() not in {"", "0", "false", "no", "off"}
    return bool(value)


def _safe_name(value, fallback="comparison"):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._-")
    return text or fallback


def _bh_adjust(pvalues):
    values = np.asarray(pvalues, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not bool(valid.any()):
        return result
    indices = np.where(valid)[0]
    ordered = indices[np.argsort(values[indices])]
    running = 1.0
    for rank in range(len(ordered), 0, -1):
        index = ordered[rank - 1]
        running = min(running, float(values[index]) * len(ordered) / rank)
        result[index] = min(1.0, running)
    return result


def _parse_pairs(value):
    pairs = []
    for raw_item in str(value or "").replace("\n", ";").split(";"):
        item = raw_item.strip()
        if not item:
            continue
        if "-vs-" not in item:
            raise ValueError("指定比较格式应为 IBD-vs-Control；多个比较用分号分隔。")
        experimental, reference = (part.strip() for part in item.split("-vs-", 1))
        if not experimental or not reference or experimental == reference:
            raise ValueError("指定比较中的两个条件必须不同且非空。")
        pair = (experimental, reference)
        if pair not in pairs:
            pairs.append(pair)
    return pairs


def _validate_design(adata, sample_key, condition_key, *, allow_technical_sample,
                     allow_technical_condition, min_cells_per_sample):
    """Validate one biological sample -> one condition and derive denominators."""
    sample_key = str(sample_key or "").strip()
    condition_key = str(condition_key or "").strip()
    if not sample_key or sample_key not in adata.obs.columns:
        raise ValueError("邻域差异丰度需要有效的生物学样本列（如 sample_id）。")
    if not condition_key or condition_key not in adata.obs.columns:
        raise ValueError("邻域差异丰度需要有效的条件列（如 condition）。")
    if sample_key == condition_key:
        raise ValueError("生物学样本列与条件列不能相同。")
    if sample_key.lower() in _TECHNICAL_BATCH_NAMES and not allow_technical_sample:
        raise ValueError(
            f"样本列 '{sample_key}' 看起来是技术 batch；请填写真实 sample_id，"
            "或明确确认它确为独立生物学样本。"
        )
    if condition_key.lower() in _TECHNICAL_BATCH_NAMES and not allow_technical_condition:
        raise ValueError(
            f"条件列 '{condition_key}' 看起来是技术 batch；请填写真实 condition，"
            "或明确确认它确为生物学条件。"
        )

    obs = adata.obs
    sample_values = obs[sample_key]
    condition_values = obs[condition_key]
    invalid = (
        sample_values.isna() | condition_values.isna()
        | sample_values.astype(str).str.strip().eq("")
        | condition_values.astype(str).str.strip().eq("")
    )
    if bool(invalid.any()):
        raise ValueError("样本列或条件列包含缺失/空值，无法进行样本级邻域比较。")
    design = pd.DataFrame({
        "sample_id": sample_values.astype(str),
        "condition": condition_values.astype(str),
    })
    n_conditions = design.groupby("sample_id", observed=True)["condition"].nunique()
    if bool((n_conditions > 1).any()):
        raise ValueError("同一个样本属于多个条件；当前邻域模块不支持配对/交互模型。")
    design = design.groupby(["sample_id", "condition"], observed=True).size().rename(
        "n_scoped_cells"
    ).reset_index()
    design["eligible_for_inference"] = (
        design["n_scoped_cells"] >= max(1, int(min_cells_per_sample))
    )
    return design


def _select_representation(adata, requested):
    requested = str(requested or "auto").strip()
    if requested != "auto":
        if requested not in adata.obsm:
            raise ValueError(f"指定的邻域表示 '{requested}' 不在 adata.obsm 中。")
        return requested
    for key in (
        "X_pca_harmony", "X_scanorama", "X_scvi", "X_scVI", "X_sysvi",
        "X_pca_combat", "X_mnn", "X_pca",
    ):
        if key in adata.obsm:
            return key
    raise ValueError(
        "未找到可用于构建邻域的高维表示（优先 X_harmony/X_scvi/X_pca）。"
        "请先运行降维或批次校正。"
    )


def _neighbourhood_memberships(representation, n_neighbourhoods, neighbourhood_size,
                               random_state):
    """Return deterministic, unique overlapping KNN neighbourhoods.

    Sampling a bounded number of seed cells keeps the module tractable for the
    internal platform while making the exact seed set reproducible.
    """
    from sklearn.neighbors import NearestNeighbors

    values = np.asarray(representation, dtype=float)
    if values.ndim != 2 or values.shape[0] < 3:
        raise ValueError("邻域分析至少需要 3 个细胞及一个二维以上的表示矩阵。")
    if not np.isfinite(values).all():
        raise ValueError("邻域表示含非有限数值，无法构建 KNN 邻域。")
    n_cells = values.shape[0]
    n_seed = max(1, min(int(n_neighbourhoods), n_cells))
    n_neighbors = max(2, min(int(neighbourhood_size), n_cells))
    rng = np.random.default_rng(int(random_state))
    seed_indices = np.sort(rng.choice(n_cells, size=n_seed, replace=False))
    model = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    model.fit(values)
    indices = model.kneighbors(values[seed_indices], return_distance=False)

    unique = {}
    for seed_index, members in zip(seed_indices, indices):
        membership = np.sort(np.asarray(members, dtype=int))
        unique.setdefault(tuple(membership.tolist()), int(seed_index))
    memberships = [np.asarray(key, dtype=int) for key in unique]
    seeds = [unique[tuple(item.tolist())] for item in memberships]
    if not memberships:
        raise ValueError("未构建出有效 KNN 邻域。")
    return np.asarray(seeds, dtype=int), memberships


def _requested_comparisons(design, explicit_pairs):
    conditions = sorted(design.loc[design["eligible_for_inference"], "condition"].unique())
    if len(conditions) < 2:
        raise ValueError("满足最低细胞数的样本不足以覆盖两个条件。")
    pairs = explicit_pairs or list(itertools.combinations(conditions, 2))
    missing = [
        f"{experimental}-vs-{reference}" for experimental, reference in pairs
        if experimental not in conditions or reference not in conditions
    ]
    if missing:
        raise ValueError("指定比较中含没有合格样本的条件: " + "、".join(missing))
    return pairs


def _neighbourhood_sample_table(memberships, adata, design):
    """Materialize zero counts for every (neighbourhood, sample) pair."""
    sample_by_cell = adata.obs["_neighborhood_da_sample"].astype(str).to_numpy()
    lookup = design.set_index("sample_id")
    sample_ids = lookup.index.astype(str).tolist()
    rows = []
    for neighbourhood_id, members in enumerate(memberships):
        observed = pd.Series(sample_by_cell[members]).value_counts()
        for sample_id in sample_ids:
            n_cells = int(observed.get(sample_id, 0))
            total = int(lookup.loc[sample_id, "n_scoped_cells"])
            rows.append({
                "neighbourhood_id": int(neighbourhood_id),
                "sample_id": str(sample_id),
                "condition": str(lookup.loc[sample_id, "condition"]),
                "n_neighbourhood_cells": n_cells,
                "n_scoped_cells": total,
                "eligible_for_inference": bool(lookup.loc[sample_id, "eligible_for_inference"]),
                "proportion": n_cells / max(total, 1),
            })
    return pd.DataFrame(rows)


def _compare_neighbourhoods(sample_table, pairs, min_samples_per_condition,
                            proportion_pseudocount):
    from scipy.stats import mannwhitneyu

    rows = []
    for experimental, reference in pairs:
        comparison = f"{experimental} vs {reference}"
        pair_table = sample_table[
            sample_table["condition"].isin([experimental, reference])
            & sample_table["eligible_for_inference"]
        ]
        for neighbourhood_id, values in pair_table.groupby("neighbourhood_id", observed=True):
            exp_values = values.loc[
                values["condition"].eq(experimental), "proportion"
            ].to_numpy(dtype=float)
            ref_values = values.loc[
                values["condition"].eq(reference), "proportion"
            ].to_numpy(dtype=float)
            eligible = (
                len(exp_values) >= int(min_samples_per_condition)
                and len(ref_values) >= int(min_samples_per_condition)
            )
            if eligible:
                try:
                    statistic, pvalue = mannwhitneyu(
                        exp_values, ref_values, alternative="two-sided",
                    )
                    test_status = "tested_sample_level"
                except ValueError:
                    statistic, pvalue, test_status = np.nan, np.nan, "test_unavailable"
            else:
                statistic, pvalue, test_status = np.nan, np.nan, "insufficient_biological_samples"
            mean_exp = float(np.mean(exp_values)) if len(exp_values) else np.nan
            mean_ref = float(np.mean(ref_values)) if len(ref_values) else np.nan
            rows.append({
                "neighbourhood_id": int(neighbourhood_id),
                "comparison": comparison,
                "experimental_condition": experimental,
                "reference_condition": reference,
                "n_samples_experimental": int(len(exp_values)),
                "n_samples_reference": int(len(ref_values)),
                "mean_proportion_experimental": mean_exp,
                "mean_proportion_reference": mean_ref,
                "mean_proportion_difference": mean_exp - mean_ref,
                "log2_proportion_fold_change": float(np.log2(
                    (mean_exp + proportion_pseudocount)
                    / (mean_ref + proportion_pseudocount)
                )),
                "statistic": float(statistic) if np.isfinite(statistic) else np.nan,
                "p_value": float(pvalue) if np.isfinite(pvalue) else np.nan,
                "test": "mann_whitney_u",
                "test_status": test_status,
                "inference_unit": "biological_sample",
            })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr_bh"] = np.nan
        for comparison, frame in result.groupby("comparison", observed=True):
            result.loc[frame.index, "fdr_bh"] = _bh_adjust(frame["p_value"].to_numpy(dtype=float))
    return result


class NeighborhoodDAAnalysis(BaseAnalysis):
    MODULE_NAME = "neighborhood_da"
    DISPLAY_NAME = "邻域差异丰度"
    DESCRIPTION = "以样本为重复的局部细胞状态差异丰度分析"
    INPUT_REQUIRES = ["X_umap"]

    def validate_input(self, adata):
        if "X_umap" not in adata.obsm:
            return "缺少 UMAP 嵌入（obsm['X_umap']）；请先运行降维分析（dimred）。"
        try:
            _select_representation(adata, self.params.get("representation", "auto"))
        except ValueError as exc:
            return str(exc)
        return None

    def run(self, input_path):
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm

        self.progress(5, "加载并验证样本级邻域设计...")
        adata = self.load_adata(input_path)
        adata = self.apply_scope(adata)
        if adata.n_obs < 3:
            raise ValueError("范围内细胞少于 3 个，无法进行邻域差异丰度分析。")

        sample_key = str(self.params.get("sample_key", "sample_id") or "").strip()
        condition_key = str(self.params.get("condition_key", "condition") or "").strip()
        min_cells_per_sample = int(self.params.get("min_cells_per_sample", 20) or 20)
        min_samples_per_condition = int(self.params.get("min_samples_per_condition", 2) or 2)
        design = _validate_design(
            adata, sample_key, condition_key,
            allow_technical_sample=_as_bool(
                self.params.get("confirm_batch_is_biological_sample"), False,
            ),
            allow_technical_condition=_as_bool(
                self.params.get("confirm_batch_is_biological_condition"), False,
            ),
            min_cells_per_sample=min_cells_per_sample,
        )
        representation_key = _select_representation(adata, self.params.get("representation", "auto"))
        representation = np.asarray(adata.obsm[representation_key], dtype=float)
        explicit_pairs = _parse_pairs(self.params.get("comparisons", ""))
        comparisons = _requested_comparisons(design, explicit_pairs)
        adata.obs["_neighborhood_da_sample"] = adata.obs[sample_key].astype(str).to_numpy()

        self.progress(22, "构建重叠 KNN 邻域...")
        seeds, memberships = _neighbourhood_memberships(
            representation,
            n_neighbourhoods=int(self.params.get("n_neighbourhoods", 100) or 100),
            neighbourhood_size=int(self.params.get("neighbourhood_size", 50) or 50),
            random_state=int(self.params.get("random_state", 0) or 0),
        )
        sample_table = _neighbourhood_sample_table(memberships, adata, design)

        self.progress(48, "按生物学样本比较邻域丰度...")
        pseudocount = float(self.params.get("proportion_pseudocount", 1e-4) or 1e-4)
        if pseudocount <= 0:
            raise ValueError("比例 fold-change 的伪计数必须大于 0。")
        tests = _compare_neighbourhoods(
            sample_table, comparisons, min_samples_per_condition, pseudocount,
        )
        if tests.empty:
            raise ValueError("未生成可比较的邻域结果。")

        celltype_key = str(self.params.get("celltype_key", "celltype") or "").strip()
        umap = np.asarray(adata.obsm["X_umap"], dtype=float)
        neighbourhood_rows = []
        for neighbourhood_id, (seed, members) in enumerate(zip(seeds, memberships)):
            row = {
                "neighbourhood_id": int(neighbourhood_id),
                "seed_cell_id": str(adata.obs_names[int(seed)]),
                "seed_index": int(seed),
                "n_cells": int(len(members)),
                "umap_1": float(umap[int(seed), 0]),
                "umap_2": float(umap[int(seed), 1]),
            }
            if celltype_key and celltype_key in adata.obs.columns:
                labels = adata.obs[celltype_key].astype(str).iloc[members]
                leading = labels.value_counts()
                row["dominant_celltype"] = str(leading.index[0]) if len(leading) else ""
                row["dominant_celltype_fraction"] = float(leading.iloc[0] / len(members)) if len(leading) else np.nan
            neighbourhood_rows.append(row)
        neighbourhoods = pd.DataFrame(neighbourhood_rows)
        tests = tests.merge(neighbourhoods, on="neighbourhood_id", how="left", validate="many_to_one")

        # The data are useful for downstream plotting, but the temporary
        # sample helper must never leak into a persisted H5AD.
        adata.obs.drop(columns=["_neighborhood_da_sample"], inplace=True)
        adata.uns["neighborhood_da"] = {
            "method": "sample_level_overlapping_knn_neighbourhoods",
            "statistical_contract": (
                "Mann-Whitney comparison of per-sample neighbourhood proportions; "
                "BH correction within each requested condition comparison."
            ),
            "representation": representation_key,
            "sample_key": sample_key,
            "condition_key": condition_key,
            "n_neighbourhoods": int(len(memberships)),
            "neighbourhood_size": int(len(memberships[0])),
            "comparisons": [f"{experimental} vs {reference}" for experimental, reference in comparisons],
            "min_cells_per_sample": min_cells_per_sample,
            "min_samples_per_condition": min_samples_per_condition,
            "proportion_pseudocount": pseudocount,
        }

        self.progress(68, "导出邻域统计和空间图...")
        run_key = _safe_name(self.params.get("_analysis_id"), "manual_run")
        results_dir = Path(self.project_dir) / "results" / "neighborhood_da" / run_key
        results_dir.mkdir(parents=True, exist_ok=True)
        plots_dir = self.ensure_plots_dir()
        result_files = []

        def save_table(frame, filename, label):
            path = results_dir / filename
            frame.to_csv(path, index=False)
            result_files.append({"file_path": str(path), "file_type": "csv", "category": "table", "label": label})

        save_table(design, "sample_design.csv", "Neighborhood DA sample design")
        save_table(neighbourhoods, "neighbourhood_metadata.csv", "Neighborhood metadata")
        save_table(sample_table, "neighbourhood_sample_proportions.csv", "Neighborhood sample-level proportions")
        save_table(tests, "neighbourhood_differential_abundance.csv", "Neighborhood differential-abundance statistics")
        manifest_path = results_dir / "analysis_manifest.json"
        manifest_path.write_text(json.dumps(adata.uns["neighborhood_da"], ensure_ascii=False, indent=2), encoding="utf-8")
        result_files.append({"file_path": str(manifest_path), "file_type": "json", "category": "metadata", "label": "Neighborhood DA analysis manifest"})

        for comparison, table in tests.groupby("comparison", observed=True):
            values = table["log2_proportion_fold_change"].to_numpy(dtype=float)
            limit = max(float(np.nanmax(np.abs(values))) if len(values) else 0.0, 0.25)
            figure, axis = plt.subplots(figsize=(7.0, 5.8), dpi=180)
            axis.scatter(umap[:, 0], umap[:, 1], s=2.2, c="#D0D5DD", alpha=0.28, linewidths=0)
            points = axis.scatter(
                table["umap_1"], table["umap_2"], s=32,
                c=values, cmap="coolwarm", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                edgecolors="white", linewidths=0.35, zorder=2,
            )
            significant = table["fdr_bh"].le(float(self.params.get("fdr_cutoff", 0.05) or 0.05)).fillna(False)
            if bool(significant.any()):
                axis.scatter(
                    table.loc[significant, "umap_1"], table.loc[significant, "umap_2"],
                    s=72, facecolors="none", edgecolors="#111827", linewidths=0.9,
                    zorder=3, label=f"FDR ≤ {float(self.params.get('fdr_cutoff', 0.05) or 0.05):.2g}",
                )
                axis.legend(frameon=False, loc="best", fontsize=8)
            colorbar = figure.colorbar(points, ax=axis, pad=0.02)
            colorbar.set_label("log2 sample-proportion fold change", fontsize=8)
            axis.set_title(f"Neighbourhood differential abundance: {comparison}", loc="left", fontsize=10, fontweight="semibold")
            axis.set_xlabel("UMAP1")
            axis.set_ylabel("UMAP2")
            axis.grid(False)
            filename = f"neighbourhood_da_{_safe_name(comparison)}.png"
            result_files.extend(self.save_matplotlib_figure(
                figure, plots_dir, filename, "umap",
                f"Neighborhood DA: {comparison}", formats=("png", "svg"), dpi=300,
            ))

        self.progress(90, "保存邻域差异丰度结果...")
        output_path = self.save_output(adata, "neighborhood_da")
        self.progress(100, "完成")
        n_tested = int(tests["test_status"].eq("tested_sample_level").sum())
        return {
            "output_adata": output_path,
            "result_files": result_files,
            "summary": {
                "n_cells": int(adata.n_obs),
                "n_neighbourhoods": int(len(memberships)),
                "neighbourhood_size": int(len(memberships[0])),
                "representation": representation_key,
                "comparisons": [f"{experimental} vs {reference}" for experimental, reference in comparisons],
                "n_tested_neighbourhoods": n_tested,
                "n_significant_neighbourhoods": int(tests["fdr_bh"].le(float(self.params.get("fdr_cutoff", 0.05) or 0.05)).sum()),
                "inference_unit": "biological_sample",
                "scope_key": str(self.params.get("scope_key", "") or "").strip() or None,
            },
        }
