"""Exploratory cell-level condition DEG exports for a downstream plot platform.

This module is intentionally distinct from ``sc_pseudobulk_deg``.  It supplies
complete gene-level log2FC/FDR tables for a selected control-versus-treatment
pair even when an experiment has only one library per condition.  Such a
comparison is useful for exploration and visualization, but is marked in every
row as cell-level rather than biological-sample-level inference.
"""

import os
import re
from itertools import combinations

import numpy as np
import pandas as pd

from modules.base import BaseAnalysis
from modules.sc_batch import _available_path
from modules.sc_batch_export import batch_csv_package_dirs


def _safe_name(value, fallback="comparison"):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._-")
    return text or fallback


def _comparison_pairs(groups, requested_comparisons="", mode="all_pairwise", reference_group=""):
    """Return ordered (experimental, control) pairs without assuming sample IDs."""
    groups = sorted({str(value) for value in groups})
    requested = str(requested_comparisons or "").strip()
    if requested:
        pairs = []
        for raw in re.split(r"[;\n]+", requested):
            item = raw.strip()
            if not item:
                continue
            parts = re.split(r"\s*(?:-vs-|\s+vs\s+)\s*", item, maxsplit=1,
                             flags=re.IGNORECASE)
            if len(parts) != 2 or not all(part.strip() for part in parts):
                raise ValueError(f"比较格式无效: '{item}'；请使用 Experimental-vs-Control")
            experimental, control = (part.strip() for part in parts)
            if experimental not in groups or control not in groups:
                raise ValueError(f"比较组不存在: '{item}'")
            if experimental != control and (experimental, control) not in pairs:
                pairs.append((experimental, control))
        return pairs
    if mode == "vs_reference":
        reference = str(reference_group or "").strip()
        if reference not in groups:
            raise ValueError("vs_reference 模式需要填写存在的 reference_group")
        return [(group, reference) for group in groups if group != reference]
    return list(combinations(groups, 2))


class SCCellLevelDEG(BaseAnalysis):
    """Write exploratory DEG tables at whole-dataset and cluster resolution."""

    MODULE_NAME = "sc_cell_deg"
    DISPLAY_NAME = "单细胞级条件 DEG 导出"
    DESCRIPTION = "按对照/实验组导出完整 log2FC、p.adjust 的探索性单细胞 DEG CSV"

    def run(self, input_path):
        import scanpy as sc

        adata = self.load_adata(input_path)
        condition_key = str(self.params.get("condition_key", "condition") or "").strip()
        cluster_key = str(self.params.get("cluster_key", "leiden") or "").strip()
        scope = str(self.params.get("analysis_scope", "both") or "both")
        min_cells = max(2, int(self.params.get("min_cells_per_group", 10)))
        method = str(self.params.get("method", "wilcoxon") or "wilcoxon")
        comparison_mode = str(self.params.get("comparison_mode", "all_pairwise") or "all_pairwise")
        reference_group = str(self.params.get("reference_group", "") or "")
        requested_comparisons = str(self.params.get("comparisons", "") or "")
        export_folder = str(self.params.get("export_folder", "sc_batch_results") or "sc_batch_results")
        prefix = _safe_name(self.params.get("export_prefix"), "sc_cell_level")

        if condition_key not in adata.obs.columns:
            raise ValueError(f"adata.obs 中缺少条件列 '{condition_key}'")
        conditions = adata.obs[condition_key].astype(str).str.strip()
        if bool(conditions.eq("").any()):
            raise ValueError(f"条件列 '{condition_key}' 存在空值")
        groups = sorted(conditions.unique().tolist())
        if len(groups) < 2:
            raise ValueError("条件列至少需要两个组")
        pairs = _comparison_pairs(
            groups, requested_comparisons=requested_comparisons,
            mode=comparison_mode, reference_group=reference_group,
        )
        if scope not in {"all_cells", "per_cluster", "both"}:
            raise ValueError("analysis_scope 必须为 all_cells、per_cluster 或 both")
        scopes = (["all_cells", "per_cluster"] if scope == "both" else [scope])
        if "per_cluster" in scopes:
            if cluster_key not in adata.obs.columns:
                raise ValueError(f"按簇 DEG 需要聚类列 '{cluster_key}'")
            cluster_values = sorted(adata.obs[cluster_key].astype(str).unique().tolist())

        package_root, package_dirs = batch_csv_package_dirs(self.project_dir, export_folder)
        result_files = []
        manifest_rows = []
        units_per_comparison = sum(
            1 if item == "all_cells" else len(cluster_values) for item in scopes
        )
        total_units = max(1, len(pairs) * units_per_comparison)
        current_unit = 0
        for experimental, control in pairs:
            comparison = f"{experimental} vs {control}"
            comparison_id = f"{_safe_name(experimental)}_vs_{_safe_name(control)}"
            for scope_name in scopes:
                clusters = ["All"] if scope_name == "all_cells" else cluster_values
                result_rows = []
                comparison_status = "exploratory_runnable"
                comparison_reason = (
                    "细胞级 Wilcoxon；细胞不是独立生物学重复，结果仅供探索性可视化。"
                )
                for cluster in clusters:
                    current_unit += 1
                    unit_name = "全部细胞" if cluster == "All" else f"cluster {cluster}"
                    self.progress(
                        5 + int((current_unit - 1) * 80 / total_units),
                        f"正在计算 {comparison} / {unit_name} 的单细胞 DEG...",
                    )
                    mask = conditions.isin([experimental, control]).to_numpy()
                    if cluster != "All":
                        mask &= adata.obs[cluster_key].astype(str).eq(cluster).to_numpy()
                    subset = adata[mask].copy()
                    labels = subset.obs[condition_key].astype(str)
                    n_experimental = int(labels.eq(experimental).sum())
                    n_control = int(labels.eq(control).sum())
                    status_base = {
                        "comparison_id": comparison_id, "comparison": comparison,
                        "experimental_group": experimental, "control_group": control,
                        "deg_scope": scope_name, "cluster": cluster,
                        "n_cells_experimental": n_experimental,
                        "n_cells_control": n_control,
                    }
                    if n_experimental < min_cells or n_control < min_cells:
                        manifest_rows.append({
                            **status_base, "status": "not_runnable",
                            "reason": f"每组该簇至少需要 {min_cells} 个细胞",
                        })
                        continue
                    subset.obs["_cell_deg_condition"] = labels.astype("category")
                    sc.tl.rank_genes_groups(
                        subset, groupby="_cell_deg_condition", groups=[experimental],
                        reference=control, method=method, n_genes=subset.n_vars, pts=True,
                    )
                    table = sc.get.rank_genes_groups_df(subset, group=experimental)
                    table = table.rename(columns={
                        "names": "gene", "logfoldchanges": "log2FC", "pvals": "p_val",
                        "pvals_adj": "p.adjust", "scores": "statistic",
                        "pct_nz_group": "pct.1", "pct_nz_reference": "pct.2",
                    })
                    table["comparison_id"] = comparison_id
                    table["comparison"] = comparison
                    table["experimental_group"] = experimental
                    table["control_group"] = control
                    table["deg_scope"] = scope_name
                    table["cluster"] = cluster
                    table["avg_log2FC"] = table["log2FC"]
                    table["p_val_adj"] = table["p.adjust"]
                    table["n_cells_experimental"] = n_experimental
                    table["n_cells_control"] = n_control
                    table["inference_unit"] = "cell"
                    table["statistical_status"] = "exploratory_no_biological_replicates"
                    table["method"] = method
                    result_rows.append(table)
                    manifest_rows.append({
                        **status_base, "status": "exploratory_runnable",
                        "reason": comparison_reason,
                    })

                if result_rows:
                    frame = pd.concat(result_rows, ignore_index=True)
                    output_path = _available_path(
                        package_dirs["deg"],
                        f"{prefix}_deg_{scope_name}_{comparison_id}", ".csv",
                    )
                    frame.to_csv(output_path, index=False)
                    result_files.append({
                        "file_path": output_path, "file_type": "csv", "category": "table",
                        "label": (
                            f"Cell-level exploratory DEG ({scope_name}): {comparison}"
                        ),
                    })
                else:
                    comparison_status = "not_runnable"
                    comparison_reason = (
                        f"{scope_name} 没有任何比较单元同时满足每组 {min_cells} 个细胞"
                    )
                # Comparison-wide status complements the per-cluster status rows.
                manifest_rows.append({
                    "comparison_id": comparison_id, "comparison": comparison,
                    "experimental_group": experimental, "control_group": control,
                    "deg_scope": scope_name, "cluster": "__comparison__",
                    "n_cells_experimental": int(conditions.eq(experimental).sum()),
                    "n_cells_control": int(conditions.eq(control).sum()),
                    "status": comparison_status, "reason": comparison_reason,
                })

        self.progress(90, "正在保存单细胞 DEG 比较注册表...")
        manifest = pd.DataFrame(manifest_rows)
        manifest_path = _available_path(
            package_dirs["deg"], f"{prefix}_comparison_manifest", ".csv"
        )
        manifest.to_csv(manifest_path, index=False)
        result_files.append({
            "file_path": manifest_path, "file_type": "csv", "category": "table",
            "label": "Cell-level DEG comparison registry",
        })
        self.progress(100, "单细胞级探索性 DEG CSV 导出完成")
        summary = {
            "analysis_level": "cell_level_exploratory",
            "inference_unit": "cell",
            "statistical_status": "exploratory_no_biological_replicates",
            "condition_key": condition_key,
            "cluster_key": cluster_key if "per_cluster" in scopes else "",
            "analysis_scopes": scopes,
            "n_comparisons": int(len(pairs)), "csv_package_dir": package_root,
            "warning": "细胞级 p.adjust 不能替代有生物学重复的样本级推断。",
        }
        return {
            "output_adata": input_path,
            "result_files": result_files,
            "summary": summary,
        }
