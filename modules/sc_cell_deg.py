"""Exploratory cell-level differential expression with explicit comparison contracts.

This module deliberately does *not* call cells biological replicates.  It is
for three useful descriptive questions: condition contrasts, clusters within
one selected sample, and samples compared within a selected cluster (or across
all cells).  Confirmatory condition inference belongs to ``sc_pseudobulk_deg``.
"""

from __future__ import annotations

import json
import os
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from modules.base import BaseAnalysis
from modules.sc_batch import (
    _available_path,
    _looks_like_control,
    comparison_id_for_groups,
)
from modules.sc_batch_export import batch_csv_package_dirs
from modules.sc_de_utils import log1p_adata_for_cell_level_de, normalization_semantics
from modules.native_figures import diverging_bar_figure


def _safe_name(value, fallback="comparison"):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._-")
    return text or fallback


def _as_bool(value, default=True):
    if value is None:
        return bool(default)
    if isinstance(value, str):
        return value.strip().lower() not in {"", "0", "false", "no", "off"}
    return bool(value)


def _comparison_pairs(groups, requested_comparisons="", mode="all_pairwise", reference_group=""):
    """Return ordered (experimental, control) pairs from explicit categories."""
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
        if not pairs:
            raise ValueError("未解析到有效比较")
        return pairs
    if mode == "vs_reference":
        reference = str(reference_group or "").strip()
        if reference not in groups:
            raise ValueError("vs_reference 模式需要填写存在的 reference_group")
        return [(group, reference) for group in groups if group != reference]
    # all_pairwise 无固有方向：control 风格标签固定放在对照组位置，
    # 保证 log2FC 正值 = 实验组上调；其余按字母序（与 pseudobulk 一致）。
    ordered = []
    for left, right in combinations(groups, 2):
        if _looks_like_control(left) and not _looks_like_control(right):
            ordered.append((right, left))
        else:
            ordered.append((left, right))
    return ordered


class SCCellLevelDEG(BaseAnalysis):
    """Run descriptive cell-level contrasts on a safe log1p expression scale."""

    MODULE_NAME = "sc_cell_deg"
    DISPLAY_NAME = "探索性细胞级比较"
    DESCRIPTION = "细胞级探索性比较；不将细胞当作生物学重复"

    def _comparison_contract(self, adata):
        """Build comparison units without conflating their statistical meaning."""
        comparison_type = str(
            self.params.get("comparison_type", "condition") or "condition"
        ).strip().lower()
        aliases = {
            "condition": "condition",
            "within_sample_clusters": "within_sample_clusters",
            "cluster_within_sample": "within_sample_clusters",
            "between_samples_within_cluster": "between_samples_within_cluster",
            "sample_within_cluster": "between_samples_within_cluster",
            "between_samples_all_cells": "between_samples_all_cells",
            "sample_all_cells": "between_samples_all_cells",
        }
        if comparison_type not in aliases:
            raise ValueError("comparison_type 必须为 condition、within_sample_clusters、between_samples_within_cluster 或 between_samples_all_cells")
        comparison_type = aliases[comparison_type]
        condition_key = str(self.params.get("condition_key", "condition") or "").strip()
        cluster_key = str(self.params.get("cluster_key", "leiden") or "").strip()
        sample_key = str(self.params.get("sample_key", "sample_id") or "").strip()
        selected_sample = str(self.params.get("selected_sample", "") or "").strip()
        selected_cluster = str(self.params.get("selected_cluster", "") or "").strip()
        requested_scope = str(self.params.get("analysis_scope", "both") or "both")

        if comparison_type == "condition":
            if condition_key not in adata.obs.columns:
                available = ", ".join(str(item) for item in adata.obs.columns[:10])
                raise ValueError(
                    f"adata.obs 中缺少条件列 '{condition_key}'；当前 obs 列：{available}"
                    + ("…" if len(adata.obs.columns) > 10 else "")
                    + "。条件间比较需要真实条件元数据（如 condition/treatment/group），"
                    "请先通过批量 10x manifest 导入或在 h5ad.obs 中补充后重试。"
                )
            labels = adata.obs[condition_key].astype(str).str.strip()
            if bool(labels.eq("").any()):
                raise ValueError(f"条件列 '{condition_key}' 存在空值")
            groups = sorted(labels.unique().tolist())
            if len(groups) < 2:
                raise ValueError("条件列至少需要两个组")
            if requested_scope not in {"all_cells", "per_cluster", "both"}:
                raise ValueError("analysis_scope 必须为 all_cells、per_cluster 或 both")
            scopes = ["all_cells", "per_cluster"] if requested_scope == "both" else [requested_scope]
            units = []
            for scope in scopes:
                if scope == "all_cells":
                    units.append((scope, "All", pd.Series(True, index=adata.obs_names), condition_key, groups))
                    continue
                if cluster_key not in adata.obs.columns:
                    raise ValueError(f"按簇 DEG 需要聚类列 '{cluster_key}'")
                for cluster in sorted(adata.obs[cluster_key].astype(str).unique().tolist()):
                    units.append((scope, cluster, adata.obs[cluster_key].astype(str).eq(cluster), condition_key, groups))
            return {
                "comparison_type": comparison_type, "condition_key": condition_key,
                "cluster_key": cluster_key, "sample_key": sample_key,
                "selected_sample": "", "selected_cluster": "", "units": units,
                "comparison_unit": "cell", "description": "条件间细胞级探索性比较；不能替代样本级推断。",
            }

        if sample_key not in adata.obs.columns:
            available = ", ".join(str(item) for item in adata.obs.columns[:10])
            raise ValueError(
                f"该比较需要样本列 '{sample_key}'；当前 obs 列：{available}"
                + ("…" if len(adata.obs.columns) > 10 else "")
                + "。请先通过批量 10x manifest 导入写入 sample_id，"
                "或在 h5ad.obs 中补充样本元数据后重试。"
            )
        samples = adata.obs[sample_key].astype(str).str.strip()
        if bool(samples.eq("").any()):
            raise ValueError(f"样本列 '{sample_key}' 存在空值")

        if comparison_type == "within_sample_clusters":
            if cluster_key not in adata.obs.columns:
                raise ValueError(f"同一样本不同簇比较需要聚类列 '{cluster_key}'")
            if not selected_sample:
                raise ValueError("同一样本不同簇比较必须填写 selected_sample")
            sample_mask = samples.eq(selected_sample)
            if not bool(sample_mask.any()):
                raise ValueError(f"selected_sample '{selected_sample}' 不在 '{sample_key}' 中")
            clusters = adata.obs.loc[sample_mask, cluster_key].astype(str).str.strip()
            if len(clusters.unique()) < 2:
                raise ValueError("所选样本至少需要两个非空 cluster")
            return {
                "comparison_type": comparison_type, "condition_key": condition_key,
                "cluster_key": cluster_key, "sample_key": sample_key,
                "selected_sample": selected_sample, "selected_cluster": "", "units": [
                    ("within_sample_clusters", "SelectedSample", sample_mask, cluster_key,
                     sorted(clusters.unique().tolist()))
                ], "comparison_unit": "cell",
                "description": "同一生物样本内的 cluster marker 对比；仅作细胞级描述。",
            }

        if comparison_type == "between_samples_within_cluster":
            if cluster_key not in adata.obs.columns:
                raise ValueError(f"不同样本同一簇比较需要聚类列 '{cluster_key}'")
            if not selected_cluster:
                raise ValueError("不同样本同一簇比较必须填写 selected_cluster")
            cluster_mask = adata.obs[cluster_key].astype(str).eq(selected_cluster)
            available = sorted(samples.loc[cluster_mask].unique().tolist())
            if len(available) < 2:
                raise ValueError("所选 cluster 至少需要来自两个样本的细胞")
            return {
                "comparison_type": comparison_type, "condition_key": condition_key,
                "cluster_key": cluster_key, "sample_key": sample_key,
                "selected_sample": "", "selected_cluster": selected_cluster, "units": [
                    ("within_cluster", selected_cluster, cluster_mask, sample_key, available)
                ], "comparison_unit": "cell",
                "description": "不同样本在同一 cluster 内的细胞级探索性比较；样本对只有描述性意义。",
            }

        available = sorted(samples.unique().tolist())
        if len(available) < 2:
            raise ValueError("不同样本比较至少需要两个样本")
        return {
            "comparison_type": comparison_type, "condition_key": condition_key,
            "cluster_key": cluster_key, "sample_key": sample_key,
            "selected_sample": "", "selected_cluster": "", "units": [
                ("all_cells", "All", pd.Series(True, index=adata.obs_names), sample_key, available)
            ], "comparison_unit": "cell",
            "description": "不同样本整体细胞级探索性比较；正式条件结论请使用样本级 pseudobulk。",
        }

    def run(self, input_path):
        import scanpy as sc

        adata = self.load_adata(input_path)
        adata = self.apply_filters(adata, self.MODULE_NAME)
        if adata.n_obs == 0:
            raise ValueError("筛选后没有细胞可用于差异分析")
        contract = self._comparison_contract(adata)
        comparison_mode = str(self.params.get("comparison_mode", "all_pairwise") or "all_pairwise")
        if comparison_mode not in {"all_pairwise", "vs_reference"}:
            raise ValueError("comparison_mode 必须为 all_pairwise 或 vs_reference")
        reference_group = str(self.params.get("reference_group", "") or "")
        requested_comparisons = str(self.params.get("comparisons", "") or "")
        min_cells = max(2, int(self.params.get("min_cells_per_group", 10)))
        method = str(self.params.get("method", "wilcoxon") or "wilcoxon")
        if method not in {"wilcoxon", "t-test", "t-test_overestim_var"}:
            raise ValueError(
                "细胞级 method 必须为 wilcoxon、t-test 或 t-test_overestim_var"
            )
        export_folder = str(self.params.get("export_folder", "sc_batch_results") or "sc_batch_results")
        prefix = _safe_name(self.params.get("export_prefix"), "sc_cell_level")
        export_full_tables = _as_bool(self.params.get("export_full_tables"), default=False)
        analysis_id = _safe_name(self.params.get("_analysis_id"), "cell_deg")
        show_deg_figures = _as_bool(self.params.get("show_deg_figures"), default=True)
        max_plot_units = min(24, max(1, int(self.params.get("plot_max_units", 12))))
        volcano_top_n = min(20, max(0, int(self.params.get("volcano_top_n", 8))))
        padj_cutoff = float(self.params.get("padj_cutoff", 0.05))
        log2fc_cutoff = float(self.params.get("log2fc_cutoff", 1.0))
        if not np.isfinite(padj_cutoff) or not 0 <= padj_cutoff <= 1:
            raise ValueError("padj_cutoff 必须位于 [0, 1]")
        if not np.isfinite(log2fc_cutoff) or log2fc_cutoff < 0:
            raise ValueError("log2fc_cutoff 必须是非负有限数值")

        self.progress(5, "正在从原始 counts 构建 log1p 表达层...")
        expression = log1p_adata_for_cell_level_de(adata)
        package_root, package_dirs = batch_csv_package_dirs(
            self.project_dir, export_folder, included_keys=("deg",),
        )
        internal_dir = os.path.join(package_dirs["deg"], ".internal")
        os.makedirs(internal_dir, exist_ok=True)

        result_files = []
        manifest_rows = []
        frames_by_output = {}
        all_frames = []
        total_units = sum(len(_comparison_pairs(
            groups, requested_comparisons, comparison_mode, reference_group
        )) for _, _, _, _, groups in contract["units"])
        completed_units = 0

        for scope_name, cluster, base_mask, grouping_key, groups in contract["units"]:
            pairs = _comparison_pairs(groups, requested_comparisons, comparison_mode, reference_group)
            labels = expression.obs[grouping_key].astype(str).str.strip()
            for experimental, control in pairs:
                completed_units += 1
                comparison = f"{experimental} vs {control}"
                comparison_id = comparison_id_for_groups(experimental, control)
                self.progress(
                    8 + int((completed_units - 1) * 75 / max(total_units, 1)),
                    f"正在计算 {comparison} / {scope_name} 的探索性细胞级 DEG...",
                )
                mask = base_mask.to_numpy() & labels.isin([experimental, control]).to_numpy()
                subset = expression[mask].copy()
                subset_labels = subset.obs[grouping_key].astype(str).str.strip()
                n_experimental = int(subset_labels.eq(experimental).sum())
                n_control = int(subset_labels.eq(control).sum())
                status_base = {
                    "comparison_id": comparison_id, "comparison": comparison,
                    "experimental_group": experimental, "control_group": control,
                    "deg_scope": scope_name, "cluster": cluster,
                    "comparison_type": contract["comparison_type"],
                    "grouping_key": grouping_key,
                    "selected_sample": contract["selected_sample"],
                    "selected_cluster": contract["selected_cluster"],
                    "n_cells_experimental": n_experimental,
                    "n_cells_control": n_control,
                }
                if n_experimental < min_cells or n_control < min_cells:
                    manifest_rows.append({
                        **status_base, "status": "not_runnable",
                        "reason": f"每组该比较单元至少需要 {min_cells} 个细胞",
                    })
                    continue
                subset.obs["_cell_deg_group"] = subset_labels.astype("category")
                sc.tl.rank_genes_groups(
                    subset, groupby="_cell_deg_group", groups=[experimental], reference=control,
                    method=method, n_genes=subset.n_vars, pts=True, use_raw=False,
                )
                table = sc.get.rank_genes_groups_df(subset, group=experimental).rename(columns={
                    "names": "gene", "logfoldchanges": "log2FC", "pvals": "p_val",
                    "pvals_adj": "p.adjust", "scores": "statistic",
                    "pct_nz_group": "pct.1", "pct_nz_reference": "pct.2",
                })
                # scanpy >=1.11 在显式 groups+reference 时只返回 pct_nz_group，
                # pct.2（对照组检出率）需要从 uns['rank_genes_groups']['pts'] 按
                # 基因名补回，否则 Seurat 风格 pct.1/pct.2 契约残缺。
                if "pct.2" not in table.columns:
                    pts = subset.uns.get("rank_genes_groups", {}).get("pts")
                    ref_pct = {}
                    if isinstance(pts, dict) and control in pts:
                        ref_values = np.asarray(pts[control], dtype=float)
                        ref_pct = {
                            str(gene): float(value)
                            for gene, value in zip(subset.var_names, ref_values)
                            if np.isfinite(value)
                        }
                    elif hasattr(pts, "loc") and control in getattr(pts, "columns", []):
                        ref_pct = {
                            str(gene): float(value)
                            for gene, value in pts[control].items()
                            if np.isfinite(value)
                        }
                    table["pct.2"] = table["gene"].astype(str).map(ref_pct).fillna(0.0)
                # MA 图的横轴需要每个基因在本比较单元中的平均表达量。这里
                # 使用从受保护 counts 重建的 log1p 表达层的均值；不把原始
                # UMI 总数误当作已标准化表达，也不引入额外的 DE 统计检验。
                mean_expression = np.asarray(subset.X.mean(axis=0)).reshape(-1)
                mean_by_gene = pd.Series(
                    mean_expression,
                    index=pd.Index(subset.var_names).astype(str),
                )
                table["mean_expression"] = table["gene"].astype(str).map(mean_by_gene)
                for key, value in status_base.items():
                    table[key] = value
                table["avg_log2FC"] = table["log2FC"]
                table["p_val_adj"] = table["p.adjust"]
                table["inference_unit"] = "cell"
                table["statistical_status"] = "exploratory_no_biological_replicates"
                table["method"] = method
                table["expression_scale"] = "log1p_from_layers[counts]"
                all_frames.append(table)
                frames_by_output.setdefault((comparison_id, scope_name), []).append(table)
                manifest_rows.append({
                    **status_base, "status": "exploratory_runnable",
                    "reason": contract["description"],
                })

        internal_table = ""
        internal_sources = []
        if all_frames:
            # Keep one internal source per contrast/scope.  This is not a
            # user-facing CSV export; it lets downstream enrichment bind to a
            # precise DEG unit without guessing from file modification time.
            for (comparison_id, scope_name), frames in sorted(frames_by_output.items()):
                frame = pd.concat(frames, ignore_index=True)
                internal_path = _available_path(
                    internal_dir,
                    f"{prefix}_task_{analysis_id}_deg_{scope_name}_{comparison_id}", ".csv",
                )
                frame.to_csv(internal_path, index=False)
                internal_sources.append(internal_path)
            internal_table = internal_sources[0]
            if export_full_tables:
                for (comparison_id, scope_name), frames in sorted(frames_by_output.items()):
                    frame = pd.concat(frames, ignore_index=True)
                    output_path = _available_path(
                        package_dirs["deg"], f"{prefix}_deg_{scope_name}_{comparison_id}", ".csv",
                    )
                    frame.to_csv(output_path, index=False)
                    result_files.append({
                        "file_path": output_path, "file_type": "csv", "category": "table",
                        "label": f"Cell-level exploratory DEG ({scope_name}): {frame['comparison'].iloc[0]}",
                    })

        plot_warnings = []
        plot_units_generated = 0
        if all_frames and show_deg_figures:
            self.progress(88, "正在生成细胞级 Volcano、MA 与 DEG 计数图...")
            from figure_engine import NatureFigureDirector

            plots_dir = self.ensure_plots_dir()
            director = NatureFigureDirector()
            plotted = pd.concat(all_frames, ignore_index=True)
            plot_groups = []
            for keys, frame in plotted.groupby(
                ["comparison_id", "comparison", "deg_scope", "cluster"],
                observed=True,
                sort=False,
            ):
                cells = int(frame["n_cells_experimental"].iloc[0]) + int(
                    frame["n_cells_control"].iloc[0]
                )
                # 全细胞比较优先；其余按可用细胞数排序，避免每次运行导出
                # 过多的 cluster 图并淹没结果浏览页。
                is_all_cells = str(keys[2]) == "all_cells"
                plot_groups.append((not is_all_cells, -cells, keys, frame.copy()))
            plot_groups.sort(key=lambda item: (item[0], item[1], tuple(map(str, item[2]))))
            if len(plot_groups) > max_plot_units:
                plot_warnings.append(
                    f"细胞级 DEG 图仅展示优先级最高的 {max_plot_units} 个比较单元；"
                    "完整结果保留在受控内部表。"
                )

            bar_labels, up_counts, down_counts = [], [], []
            for _, _, keys, frame in plot_groups[:max_plot_units]:
                comparison_id, comparison, scope_name, cluster = keys
                view = frame.rename(columns={"p.adjust": "padj"}).copy()
                view["padj"] = pd.to_numeric(view["padj"], errors="coerce")
                view["log2FC"] = pd.to_numeric(view["log2FC"], errors="coerce")
                significant = view["padj"].lt(padj_cutoff)
                up_counts.append(int((significant & view["log2FC"].ge(log2fc_cutoff)).sum()))
                down_counts.append(int((significant & view["log2FC"].le(-log2fc_cutoff)).sum()))
                suffix = "all cells" if str(scope_name) == "all_cells" else f"cluster {cluster}"
                bar_labels.append(f"{comparison} · {suffix}")

                safe_stem = _safe_name(
                    f"{prefix}_{comparison_id}_{scope_name}_{cluster}", "comparison"
                )
                title_suffix = f"{comparison} · {suffix}"
                try:
                    volcano_spec = director.spec_from_params(
                        "volcano", self.params, width="single",
                        title=f"Cell-level exploratory DEG · {title_suffix}",
                        evidence_role="discovery", fc_threshold=log2fc_cutoff,
                        fdr_threshold=padj_cutoff, label_n=volcano_top_n,
                        show_legend=True,
                    )
                    fig_volcano = director.render(volcano_spec, view)
                    fig_volcano._nature_semantic_warnings.append(
                        "细胞级探索性比较；细胞不是独立生物学重复，不能替代样本级 pseudobulk 推断。"
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_volcano, plots_dir, f"sc_cell_deg_volcano_{safe_stem}.png",
                        "volcano", f"Cell-level exploratory Volcano: {title_suffix}",
                    ))

                    ma_spec = director.spec_from_params(
                        "ma", self.params, width="single",
                        title=f"Cell-level exploratory MA · {title_suffix}",
                        evidence_role="discovery", fc_threshold=log2fc_cutoff,
                        fdr_threshold=padj_cutoff, label_n=min(volcano_top_n, 6),
                        show_legend=False,
                    )
                    fig_ma = director.render(ma_spec, view)
                    fig_ma._nature_semantic_warnings.append(
                        "MA 图横轴为本比较单元中 log1p 标准化表达的平均值；仅作细胞级探索性展示。"
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_ma, plots_dir, f"sc_cell_deg_ma_{safe_stem}.png",
                        "ma", f"Cell-level exploratory MA: {title_suffix}",
                    ))
                    plot_units_generated += 1
                except Exception as exc:
                    plot_warnings.append(f"{title_suffix}: 细胞级 DEG 绘图失败（{exc}）")

            if bar_labels:
                try:
                    fig_bar = diverging_bar_figure(
                        bar_labels, up_counts, down_counts,
                        title=(
                            "Cell-level exploratory significant DEG counts "
                            f"(FDR<{padj_cutoff}, |log2FC|≥{log2fc_cutoff})"
                        ),
                        x_label="comparison unit", y_label="number of genes", rotation=35,
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_bar, plots_dir, "sc_cell_deg_significant_counts_bar.png",
                        "bar", "Cell-level exploratory significant DEG counts",
                        formats=("png", "svg"), dpi=300,
                    ))
                except Exception as exc:
                    plot_warnings.append(f"细胞级 DEG 计数柱状图绘制失败（{exc}）")

        audit = {
            "analysis_level": "cell_level_exploratory",
            "inference_unit": "cell",
            "statistical_status": "exploratory_no_biological_replicates",
            "comparison_contract": {key: value for key, value in contract.items() if key != "units"},
            "expression_scale": normalization_semantics(adata),
            "min_cells_per_group": min_cells,
            "method": method,
            "n_runnable_units": int(sum(row["status"] == "exploratory_runnable" for row in manifest_rows)),
            "units": manifest_rows,
            "full_table_exported": export_full_tables,
            "plot_units_generated": plot_units_generated,
            "plot_warnings": plot_warnings,
        }
        audit_path = _available_path(internal_dir, f"{prefix}_task_{analysis_id}_audit", ".json")
        Path(audit_path).write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
        # Keep the statistical audit first: callers that need the task contract
        # must not accidentally select a figure-engine readiness sidecar.
        result_files.insert(0, {
            "file_path": audit_path, "file_type": "json", "category": "qc",
            "label": "Cell-level DEG statistical audit",
        })
        if export_full_tables:
            manifest_path = _available_path(package_dirs["deg"], f"{prefix}_comparison_manifest", ".csv")
            pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)
            result_files.append({
                "file_path": manifest_path, "file_type": "csv", "category": "table",
                "label": "Cell-level DEG comparison registry",
            })

        self.progress(100, "单细胞级探索性 DEG 完成")
        return {
            "output_adata": input_path,
            "result_files": result_files,
            "summary": {
                "analysis_level": "cell_level_exploratory",
                "inference_unit": "cell",
                "statistical_status": "exploratory_no_biological_replicates",
                "comparison_type": contract["comparison_type"],
                "condition_key": contract["condition_key"],
                "sample_key": contract["sample_key"],
                "cluster_key": contract["cluster_key"],
                "selected_sample": contract["selected_sample"],
                "selected_cluster": contract["selected_cluster"],
                "n_runnable_units": audit["n_runnable_units"],
                "deg_source_file": internal_table,
                "deg_source_files": internal_sources,
                "deg_source_level": "cell_level",
                "deg_source_task_contract": contract["comparison_type"],
                "csv_package_dir": package_root,
                "plot_units_generated": plot_units_generated,
                "plot_warnings": plot_warnings,
                "warning": "细胞级 p.adjust 不能替代有生物学重复的样本级推断。",
            },
        }
