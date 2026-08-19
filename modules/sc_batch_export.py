"""Standard CSV export and sample-level pseudobulk DEG for batched scRNA-seq."""

import json
import os

import numpy as np
import pandas as pd

from modules.base import BaseAnalysis
from modules.sc_batch import (
    _available_path,
    aggregate_pseudobulk_counts,
    build_comparison_plan,
    sample_design_from_obs,
)


CSV_PACKAGE_FOLDERS = {
    "proportions": "01_group_proportions",
    "expression": "02_gene_expression",
    "deg": "03_differential_expression",
    "go": "04_go_enrichment",
}


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
                comparison, cluster):
    """Fit PyDESeq2 on pseudobulk counts; callers may fall back safely."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    integer_counts = np.asarray(np.rint(counts), dtype=np.int64)
    columns = [f"g{index}" for index in range(integer_counts.shape[1])]
    count_frame = pd.DataFrame(integer_counts, columns=columns)
    metadata = pd.DataFrame({"condition": pd.Categorical(conditions)})
    try:
        dds = DeseqDataSet(
            counts=count_frame, metadata=metadata, design_factors="condition",
            refit_cooks=True,
        )
    except TypeError:
        # PyDESeq2 >= 0.5 changed ``design_factors`` to a formula API.
        dds = DeseqDataSet(
            counts=count_frame, metadata=metadata, design="~condition", refit_cooks=True,
        )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", str(group_1), str(group_2)])
    stats.summary()
    result = stats.results_df.copy().reset_index(drop=True)
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
    })


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
    DESCRIPTION = "以生物学样本而非单细胞为统计单位，批量导出每个条件比较的 DEG CSV"

    def run(self, input_path):
        adata = self.load_adata(input_path)
        sample_key = str(self.params.get("sample_key", "sample_id") or "").strip()
        condition_key = str(self.params.get("condition_key", "condition") or "").strip()
        cluster_key = str(
            self.params.get("cluster_key", self.params.get("celltype_key", "leiden")) or ""
        ).strip()
        scope = str(self.params.get("analysis_scope", "per_cluster") or "per_cluster")
        if scope == "all_cells":
            cluster_key = ""
        min_cells = max(1, int(self.params.get("min_cells_per_sample_celltype", 20)))
        min_samples = max(2, int(self.params.get("min_samples_per_group", 2)))
        requested_method = str(self.params.get("method", "deseq2") or "deseq2")
        comparison_mode = str(self.params.get("comparison_mode", "all_pairwise") or "all_pairwise")
        reference_group = str(self.params.get("reference_group", "") or "")
        requested_comparisons = str(self.params.get("comparisons", "") or "")
        prefix = str(self.params.get("export_prefix", "sc_pseudobulk") or "sc_pseudobulk")
        prefix = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in prefix).strip("._-") or "sc_pseudobulk"
        export_folder = str(self.params.get("export_folder", "sc_batch_results") or "sc_batch_results")

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
            celltype_key=cluster_key, min_cells=min_cells,
        )
        pseudo_meta = pseudo_meta.rename(columns={"celltype": "cluster"})
        gene_ids, gene_names = _gene_metadata(adata)
        package_root, package_dirs = batch_csv_package_dirs(self.project_dir, export_folder)
        warnings = []
        all_results = []
        result_files = []

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
                    continue
                if requested_method == "deseq2":
                    try:
                        table = _deseq2_deg(
                            counts, meta["condition"].astype(str).to_numpy(), group_1, group_2,
                            gene_ids, gene_names, comparison, str(cluster),
                        )
                    except Exception as exc:
                        warnings.append(
                            f"{comparison} / cluster {cluster}: PyDESeq2 不可用，已回退 Welch log2CPM（{exc}）"
                        )
                        table = _welch_logcpm_deg(
                            counts, group_1_mask, group_2_mask, gene_ids, gene_names,
                            comparison, group_1, group_2, str(cluster), min_samples,
                        )
                else:
                    table = _welch_logcpm_deg(
                        counts, group_1_mask, group_2_mask, gene_ids, gene_names,
                        comparison, group_1, group_2, str(cluster), min_samples,
                    )
                if not table.empty:
                    table.insert(0, "comparison_id", item["comparison_id"])
                    table["significant"] = (table["padj"] < float(self.params.get("padj_cutoff", 0.05))) & (
                        np.abs(table["log2FC"]) >= float(self.params.get("log2fc_cutoff", 1.0))
                    )
                    comparison_rows.append(table)
            if comparison_rows:
                frame = pd.concat(comparison_rows, ignore_index=True)
                output_path = _available_path(
                    package_dirs["deg"], f"{prefix}_deg_{item['comparison_id']}", ".csv"
                )
                frame.to_csv(output_path, index=False)
                plan.loc[plan_index, "output_file"] = os.path.basename(output_path)
                plan.loc[plan_index, "n_rows"] = int(len(frame))
                plan.loc[plan_index, "n_significant"] = int(frame["significant"].sum())
                result_files.append({"file_path": output_path, "file_type": "csv", "category": "table",
                                     "label": f"Pseudobulk DEG: {comparison}"})
                all_results.append(frame)
            else:
                plan.loc[plan_index, "status"] = "not_runnable"
                plan.loc[plan_index, "reason"] = "过滤后没有满足每组最少样本数的 sample × cluster 单元"

        self.progress(90, "正在保存比较注册表和整合结果...")
        plan_path = _available_path(package_dirs["deg"], f"{prefix}_comparison_manifest", ".csv")
        plan.to_csv(plan_path, index=False)
        result_files.append({"file_path": plan_path, "file_type": "csv", "category": "table",
                             "label": "Pseudobulk DEG comparison registry"})
        if all_results:
            combined = pd.concat(all_results, ignore_index=True)
            combined_path = _available_path(package_dirs["deg"], f"{prefix}_deg_all_comparisons", ".csv")
            combined.to_csv(combined_path, index=False)
            result_files.append({"file_path": combined_path, "file_type": "csv", "category": "table",
                                 "label": "Pseudobulk DEG combined across comparisons"})
        self.progress(100, "样本级 pseudobulk DEG 完成")
        summary = {
            "n_comparisons": int(len(plan)),
            "n_completed_comparisons": int(plan.get("output_file", pd.Series(dtype=str)).notna().sum()) if "output_file" in plan else 0,
            "analysis_scope": "per_cluster" if cluster_key else "all_cells",
            "sample_key": sample_key, "condition_key": condition_key,
            "cluster_key": cluster_key,
            "csv_package_dir": package_root,
            "method_requested": requested_method,
            "warnings": warnings[:20],
        }
        return {
            "output_adata": input_path,
            "result_files": result_files,
            "summary": summary,
        }
