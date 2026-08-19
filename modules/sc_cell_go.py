"""GO ORA tables derived from the platform's cell-level DEG CSV contract."""

import os
import re
from pathlib import Path

import pandas as pd

from config import Config
from modules.base import BaseAnalysis
from modules.sc_batch import _available_path
from modules.sc_batch_export import batch_csv_package_dirs


GO_RESULT_COLUMNS = [
    "comparison_id", "comparison", "experimental_group", "control_group",
    "analysis_level", "statistical_status", "deg_scope", "cluster", "gene_set", "Term", "Overlap",
    "P-value", "Adjusted P-value", "Odds Ratio", "Combined Score", "Genes",
    "n_input_genes", "status", "reason",
]


LOCAL_GO_ASPECTS = {
    "BP": "GO_Biological_Process_2023",
    "CC": "GO_Cellular_Component_2023",
    "MF": "GO_Molecular_Function_2023",
}


def _safe_name(value, fallback="gene_set"):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._-")
    return text or fallback


def _requested_gene_sets(params):
    """Resolve one or more GO libraries without silently widening an analysis."""
    requested = str(params.get("gene_sets", "") or "").strip()
    if requested:
        libraries = [item.strip() for item in re.split(r"[,;\n]+", requested) if item.strip()]
    else:
        aspect_mode = str(params.get("go_aspects", "single") or "single")
        if aspect_mode == "BP_CC_MF":
            libraries = list(LOCAL_GO_ASPECTS.values())
        else:
            libraries = [str(params.get("gene_set", LOCAL_GO_ASPECTS["BP"]) or LOCAL_GO_ASPECTS["BP"])]
    deduplicated = []
    for library in libraries:
        if library not in deduplicated:
            deduplicated.append(library)
    return deduplicated


def _local_gene_set_path(library, directory):
    """Return a platform-local GMT/TXT file for a configured GO library."""
    base = Path(directory)
    for suffix in (".gmt", ".txt"):
        candidate = base / f"{library}{suffix}"
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"本地 GO 基因集缺失: {library}；请在 {base} 放置 {library}.gmt 或 .txt"
    )


def _resolve_local_gene_set_dir(value, project_dir):
    """Permit only the managed global cache or files inside this project."""
    managed_dir = (Path(Config.DATA_DIR) / "go_gene_sets").resolve()
    requested = Path(value).resolve() if str(value or "").strip() else managed_dir
    project_root = Path(project_dir).resolve()
    if not any(requested == root or root in requested.parents for root in (managed_dir, project_root)):
        raise ValueError("本地 GO 基因集目录必须位于平台缓存或当前项目目录内")
    return str(requested)


def _read_local_gene_sets(path):
    """Read Enrichr-style GMT/TXT without any network request."""
    gene_sets = {}
    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            term = fields[0].strip()
            genes = [gene.strip() for gene in fields[2:] if gene.strip()]
            if term and genes:
                gene_sets[term] = genes
    if not gene_sets:
        raise ValueError(f"本地 GO 基因集文件为空或格式无效: {path}")
    return gene_sets


def _run_ora(gp, selected, library, organism, execution_mode, local_gene_set_dir):
    """Run ORA locally when requested; Enrichr remains an explicit opt-in."""
    if execution_mode == "local":
        gene_sets = _read_local_gene_sets(_local_gene_set_path(library, local_gene_set_dir))
        return gp.enrich(
            gene_list=selected, gene_sets=gene_sets, outdir=None, no_plot=True,
        ).results
    return gp.enrichr(
        gene_list=selected, gene_sets=[library], organism=organism,
        outdir=None, no_plot=True,
    ).results


def _latest_deg_files(deg_dir, prefix):
    """Return the latest CSV for each (comparison, DEG scope) pair.

    Older exports did not carry ``deg_scope``.  They remain readable by
    inferring whole-dataset versus per-cluster scope from the cluster column.
    """
    candidates = sorted(Path(deg_dir).glob(f"{prefix}_deg_*.csv"))
    latest = {}
    for path in candidates:
        try:
            frame = pd.read_csv(path, nrows=2)
        except (OSError, pd.errors.EmptyDataError):
            continue
        if frame.empty or "comparison_id" not in frame.columns:
            continue
        comparison_id = str(frame["comparison_id"].iloc[0])
        if "deg_scope" in frame.columns:
            deg_scope = str(frame["deg_scope"].iloc[0])
        elif "cluster" in frame.columns and str(frame["cluster"].iloc[0]) == "All":
            deg_scope = "all_cells"
        else:
            deg_scope = "per_cluster"
        key = (comparison_id, deg_scope)
        existing = latest.get(key)
        if existing is None or path.stat().st_mtime > existing.stat().st_mtime:
            latest[key] = path
    return latest


class SCCellGOEnrichment(BaseAnalysis):
    """Run GO ORA separately for whole-dataset and per-cluster DEG exports."""

    MODULE_NAME = "sc_cell_go"
    DISPLAY_NAME = "单细胞 DEG GO 富集导出"
    DESCRIPTION = "基于单细胞级 DEG CSV，为每个比较导出一份 GO 富集表"

    def run(self, input_path):
        # Keep this module chained like other analyses, even though it consumes
        # the persisted DEG CSV package rather than expression from the h5ad.
        export_folder = str(self.params.get("export_folder", "sc_batch_results") or "sc_batch_results")
        deg_prefix = str(self.params.get("deg_prefix", "sc_cell_level") or "sc_cell_level")
        output_prefix = str(self.params.get("export_prefix", "sc_cell_go") or "sc_cell_go")
        libraries = _requested_gene_sets(self.params)
        organism = str(self.params.get("organism", "Human") or "Human")
        execution_mode = str(self.params.get("execution_mode", "local") or "local").lower()
        if execution_mode not in {"local", "enrichr"}:
            raise ValueError("execution_mode 必须为 local 或 enrichr")
        local_gene_set_dir = _resolve_local_gene_set_dir(
            self.params.get("local_gene_set_dir", ""), self.project_dir,
        )
        padj_cutoff = float(self.params.get("padj_cutoff", 0.05))
        log2fc_cutoff = float(self.params.get("log2fc_cutoff", 0.25))
        min_genes = max(3, int(self.params.get("min_genes", 5)))
        top_n = max(1, int(self.params.get("top_n", 30)))

        package_root, package_dirs = batch_csv_package_dirs(self.project_dir, export_folder)
        source_files = _latest_deg_files(package_dirs["deg"], deg_prefix)
        if not source_files:
            raise ValueError(
                f"未在 03_differential_expression 找到 '{deg_prefix}_deg_*.csv'；"
                "请先运行 sc_cell_deg，并保持相同的结果包文件夹。"
            )
        try:
            import gseapy as gp
        except ImportError as exc:
            raise ImportError("缺少 gseapy，无法运行 GO 富集") from exc

        result_files = []
        statuses = []
        total_outputs = max(1, len(source_files) * len(libraries))
        output_index = 0
        for (comparison_id, deg_scope), source_path in sorted(source_files.items()):
            deg = pd.read_csv(source_path)
            comparison = str(deg.get("comparison", pd.Series([comparison_id])).iloc[0])
            experimental = str(deg.get("experimental_group", pd.Series([""])).iloc[0])
            control = str(deg.get("control_group", pd.Series([""])).iloc[0])
            gene_column = "gene" if "gene" in deg.columns else "names"
            padj_column = "p.adjust" if "p.adjust" in deg.columns else "p_val_adj"
            if gene_column not in deg.columns or padj_column not in deg.columns or "log2FC" not in deg.columns:
                raise ValueError(f"DEG 文件缺少 gene/log2FC/{padj_column} 列: {source_path.name}")
            if deg_scope == "per_cluster" and "cluster" in deg.columns:
                units = [(str(cluster), table.copy()) for cluster, table in deg.groupby("cluster", sort=True)]
            else:
                units = [("All", deg)]
            for library in libraries:
                output_index += 1
                self.progress(
                    5 + int((output_index - 1) * 80 / total_outputs),
                    f"正在运行 {comparison_id} / {deg_scope} / {library} 的 GO 富集...",
                )
                scope_results = []
                for cluster, unit_deg in units:
                    selected = unit_deg.loc[
                        (pd.to_numeric(unit_deg[padj_column], errors="coerce") <= padj_cutoff)
                        & (pd.to_numeric(unit_deg["log2FC"], errors="coerce") >= log2fc_cutoff),
                        gene_column,
                    ].dropna().astype(str).drop_duplicates().tolist()
                    base = {
                        "comparison_id": comparison_id, "comparison": comparison,
                        "experimental_group": experimental, "control_group": control,
                        "analysis_level": "cell_level_exploratory",
                        "statistical_status": "exploratory_no_biological_replicates",
                        "deg_scope": deg_scope, "cluster": cluster,
                        "gene_set": library, "n_input_genes": len(selected),
                    }
                    if len(selected) < min_genes:
                        result = pd.DataFrame([{
                            **base, "status": "not_runnable",
                            "reason": f"显著上调基因少于 {min_genes} 个，无法进行稳定 GO ORA",
                        }])
                    else:
                        try:
                            enrichment = _run_ora(
                                gp, selected, library, organism, execution_mode, local_gene_set_dir,
                            )
                            if enrichment is None or enrichment.empty:
                                result = pd.DataFrame([{
                                    **base, "status": "empty", "reason": "未返回显著/匹配 GO term",
                                }])
                            else:
                                result = enrichment.head(top_n).copy()
                                for key, value in base.items():
                                    result[key] = value
                                result["status"] = "completed"
                                result["reason"] = ""
                        except Exception as exc:
                            source = "本地 GO 基因集" if execution_mode == "local" else "GO 服务"
                            result = pd.DataFrame([{
                                **base, "status": "failed", "reason": f"{source}不可用: {exc}",
                            }])
                    scope_results.append(result)
                result = pd.concat(scope_results, ignore_index=True)
                for column in GO_RESULT_COLUMNS:
                    if column not in result.columns:
                        result[column] = ""
                result = result[GO_RESULT_COLUMNS + [
                    column for column in result.columns if column not in GO_RESULT_COLUMNS
                ]]
                output_path = _available_path(
                    package_dirs["go"],
                    f"{output_prefix}_go_{_safe_name(library)}_{deg_scope}_{comparison_id}", ".csv",
                )
                result.to_csv(output_path, index=False)
                result_files.append({
                    "file_path": output_path, "file_type": "csv", "category": "table",
                    "label": f"Cell-level GO enrichment ({library}, {deg_scope}): {comparison}",
                })
                statuses.append({
                    "comparison_id": comparison_id, "comparison": comparison,
                    "gene_set": library, "deg_scope": deg_scope,
                    "execution_mode": execution_mode,
                    "status": "completed" if result["status"].eq("completed").any() else str(result["status"].iloc[0]),
                    "output_file": os.path.basename(output_path),
                    "n_clusters": int(result["cluster"].nunique()),
                })

        manifest_path = _available_path(
            package_dirs["go"], f"{output_prefix}_go_comparison_manifest", ".csv"
        )
        pd.DataFrame(statuses).to_csv(manifest_path, index=False)
        result_files.append({
            "file_path": manifest_path, "file_type": "csv", "category": "table",
            "label": "Cell-level GO enrichment comparison registry",
        })
        self.progress(100, "单细胞 DEG GO 富集导出完成")
        summary = {
            "n_comparisons": int(len(statuses)), "csv_package_dir": package_root,
            "analysis_level": "cell_level_exploratory", "gene_sets": libraries,
            "organism": organism, "execution_mode": execution_mode,
        }
        return {
            "output_adata": input_path,
            "result_files": result_files,
            "summary": summary,
        }
