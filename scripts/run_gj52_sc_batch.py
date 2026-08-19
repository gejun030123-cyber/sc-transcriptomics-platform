#!/usr/bin/env python3
"""Run GJ 52 organoid/tissue scRNA workflow with batch correction."""

import json
import os
import sys
import traceback
from datetime import datetime


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

PROJECT_ID = "gj52_organoid_tissue_batch"
PROJECT_NAME = "GJ52 organoid vs tissue scRNA"
REPORT_PATH = "GJ52_SC_BATCH_RESULTS_20260707.md"
GALLERY_NAME = "gj52_plot_gallery.html"

SOURCE_DIRS = [
    {
        "path": "/home/oelab/data/GJ/000000/52alg_singlet_matrix",
        "batch": "organoid",
        "sample": "52alg_singlet_matrix",
        "description": "类器官",
    },
    {
        "path": "/home/oelab/data/GJ/000000/52tissue_singlet_matrix",
        "batch": "tissue",
        "sample": "52tissue_singlet_matrix",
        "description": "对应组织",
    },
]

PIPELINE = [
    ("qc", {
        "mito_perc": 0.20,
        "nUMIs": 500,
        "detected_genes": 200,
        "max_detected_genes": 8000,
        "ribo_perc": 0,
        "hb_perc": 0,
        "batch_key": "batch",
        "batch_adaptive_qc": True,
        "mad_multiplier": 4.0,
        "save_counts_layer": True,
        "show_qc_filter_summary": True,
        "show_doublet_histogram": True,
    }),
    ("normalize", {
        "method": "log1p",
        "target_sum": 10000,
        "show_expression_distribution": True,
    }),
    ("hvg", {
        "n_top_genes": 3000,
        "hvg_flavor": "seurat_v3",
        "batch_key": "batch",
        "batch_hvg_strategy": "union",
        "exclude_mt_genes": True,
        "show_hvg_rank_plot": True,
    }),
    ("dimred", {
        "n_comps": 50,
        "umap_n_neighbors": 20,
        "umap_min_dist": 0.4,
        "umap_metric": "euclidean",
        "enable_tsne": False,
        "show_pca_scatter": True,
    }),
    ("batch_correct", {
        "method": "combat",
        "batch_key": "batch",
        "n_pcs": 50,
        "evaluate_correction": True,
    }),
    ("clustering", {
        "resolutions": "0.4,0.6,0.8,1.0,1.2",
        "primary_resolution": 0.8,
        "n_neighbors": 20,
        "clustering_method": "leiden",
        "n_iterations": 2,
        "distance_metric": "euclidean",
        "use_corrected": True,
        "batch_key": "batch",
        "show_cluster_size_bar": True,
        "show_cluster_batch_composition": True,
        "show_labeled_umap": True,
        "show_resolution_sankey": True,
    }),
    ("qc_reassess", {
        "cluster_key": "leiden",
        "doublet_threshold": 0.3,
        "mt_threshold": 15.0,
        "ribosomal_threshold": 0,
        "min_cells_per_cluster": 10,
        "auto_remove": False,
        "show_cluster_qc_bar": True,
        "show_qc_umap_panel": True,
    }),
    ("annotation", {
        "method": "auto_marker",
        "cluster_key": "leiden",
        "resolution": "0.8",
        "marker_set": "TME",
        "confidence_method": "score_margin",
        "mark_unknown": True,
        "show_celltype_composition": True,
        "show_marker_score_heatmap": True,
        "show_marker_expression_violin": True,
        "show_annotation_score_umap": True,
    }),
    ("deg", {
        "groupby": "leiden",
        "reference": "rest",
        "method": "wilcoxon",
        "n_genes": 20,
        "show_dotplot": True,
        "show_deg_counts_bar": True,
        "show_top_marker_umap_panel": True,
        "top_marker_umap_genes": 9,
        "show_marker_heatmap": True,
        "marker_heatmap_top_n": 3,
        "plot_genes_umap": "EPCAM,KRT8,KRT18,CD3D,NKG7,MS4A1,LYZ,COL1A1,PECAM1,MKI67",
        "pval_cutoff": 0.05,
        "logfc_cutoff": 0.25,
        "min_pct": 0.1,
    }),
    ("proportion", {
        "groupby": "celltype",
        "batch_key": "batch",
        "compare_groups": "organoid-vs-tissue",
        "stat_test": "chi_square",
        "n_permutations": 1000,
        "min_cells_per_group": 1,
        "show_proportion_heatmap": True,
    }),
]


def ensure_project():
    from config import Config
    from database import init_db, get_conn
    from models import Project

    init_db()
    for path in [
        Config.project_dir(PROJECT_ID),
        Config.uploads_dir(PROJECT_ID),
        Config.intermediate_dir(PROJECT_ID),
        Config.results_dir(PROJECT_ID),
        Config.plots_dir(PROJECT_ID),
        Config.branches_dir(PROJECT_ID),
    ]:
        os.makedirs(path, exist_ok=True)

    conn = get_conn()
    try:
        conn.execute("DELETE FROM result_files WHERE project_id=?", (PROJECT_ID,))
        conn.execute("DELETE FROM analysis_tasks WHERE project_id=?", (PROJECT_ID,))
        conn.execute("DELETE FROM pipeline_runs WHERE project_id=?", (PROJECT_ID,))
        conn.commit()
    finally:
        conn.close()

    project = Project(
        id=PROJECT_ID,
        name=PROJECT_NAME,
        description="GJ 52 organoid and matched tissue 10x singlet matrices, integrated by batch.",
        status="processing",
        metadata_json=json.dumps({"source_dirs": SOURCE_DIRS}, ensure_ascii=False),
    )
    project.save()
    return project


def build_merged_input():
    from config import Config
    import anndata as ad
    import scanpy as sc

    out_path = os.path.join(Config.uploads_dir(PROJECT_ID), "gj52_organoid_tissue_raw.h5ad")
    sc.settings.cachedir = Config.CACHE_DIR
    os.makedirs(sc.settings.cachedir, exist_ok=True)
    adatas = []
    raw_stats = []
    for item in SOURCE_DIRS:
        matrix_dir = item["path"]
        if not os.path.isdir(matrix_dir):
            raise FileNotFoundError(f"10x directory not found: {matrix_dir}")
        adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=True)
        adata.var_names_make_unique()
        adata.obs["barcode"] = adata.obs_names.astype(str)
        adata.obs["batch"] = item["batch"]
        adata.obs["sample"] = item["sample"]
        adata.obs["source_type"] = item["batch"]
        adata.obs["source_description"] = item["description"]
        adata.obs_names = [f"{item['batch']}:{bc}" for bc in adata.obs["barcode"].astype(str)]
        adata.layers["counts"] = adata.X.copy()
        raw_stats.append({
            "batch": item["batch"],
            "sample": item["sample"],
            "path": matrix_dir,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        })
        adatas.append(adata)

    merged = ad.concat(adatas, join="outer", merge="same", fill_value=0)
    merged.obs["batch"] = merged.obs["batch"].astype("category")
    merged.obs["sample"] = merged.obs["sample"].astype("category")
    merged.obs["source_type"] = merged.obs["source_type"].astype("category")
    merged.uns["source_dirs_json"] = json.dumps(SOURCE_DIRS, ensure_ascii=False)
    merged.uns["raw_batch_stats_json"] = json.dumps(raw_stats, ensure_ascii=False)
    merged.write_h5ad(out_path)
    return out_path, raw_stats


def run_module(module_name, params, input_path):
    from config import Config
    from models import AnalysisTask, ResultFile
    from modules import MODULE_REGISTRY

    task = AnalysisTask(
        project_id=PROJECT_ID,
        module_name=module_name,
        params_json=json.dumps(params, ensure_ascii=False),
    )
    task.save()
    task.mark_running()

    progress_log = []

    def progress_cb(pct, message):
        progress_log.append({
            "time": datetime.now().strftime("%H:%M:%S"),
            "pct": pct,
            "message": message,
        })
        try:
            task.update_progress(pct, message, json.dumps(progress_log, ensure_ascii=False))
        except Exception:
            pass

    cls = MODULE_REGISTRY[module_name]
    module = cls(
        project_dir=Config.project_dir(PROJECT_ID),
        params=params,
        progress_callback=progress_cb,
    )

    try:
        result = module.run(input_path)
        for rf in result.get("result_files") or []:
            ResultFile.create(
                task_id=task.id,
                project_id=PROJECT_ID,
                file_type=rf.get("file_type", ""),
                category=rf.get("category", ""),
                label=rf.get("label", ""),
                file_path=rf.get("file_path", ""),
            )
        output_adata = result.get("output_adata")
        result_json = json.dumps(result.get("summary", {}), ensure_ascii=False)
        task.mark_completed(output_adata, result_json)
        task.status = "completed"
        task.progress = 100
        task.progress_message = "已完成"
        task.output_adata_path = output_adata
        task.result_json = result_json
        return task, result
    except Exception:
        task.mark_failed(traceback.format_exc())
        raise


def make_gallery(result_files):
    from config import Config
    import plotly.io as pio

    out_path = os.path.join(Config.results_dir(PROJECT_ID), GALLERY_NAME)
    sections = []
    first = True
    for rf in result_files:
        fpath = rf.get("file_path", "")
        if not fpath.endswith(".json") or not os.path.isfile(fpath):
            continue
        try:
            with open(fpath, encoding="utf-8") as handle:
                fig = pio.from_json(handle.read())
            html = pio.to_html(fig, include_plotlyjs=("include" if first else False), full_html=False)
            first = False
            sections.append(f"<section><h2>{rf.get('label', os.path.basename(fpath))}</h2>{html}</section>")
        except Exception as exc:
            sections.append(
                f"<section><h2>{rf.get('label', os.path.basename(fpath))}</h2>"
                f"<pre>Could not render {fpath}: {exc}</pre></section>"
            )

    doc = """<!doctype html>
<html><head><meta charset="utf-8"><title>GJ52 scRNA batch-corrected gallery</title>
<style>
body{font-family:Arial,sans-serif;margin:24px;background:#f8f9fb;color:#1f2933}
section{background:#fff;border:1px solid #dde3ea;border-radius:8px;margin:0 0 20px;padding:16px}
h1{margin:0 0 8px} h2{font-size:16px;margin:0 0 12px}
.meta{color:#667085;margin-bottom:20px}
</style></head><body>
<h1>GJ52 scRNA Batch-Corrected Plot Gallery</h1>
<div class="meta">Organoid and matched tissue integrated by batch. Generated by scripts/run_gj52_sc_batch.py</div>
""" + "\n".join(sections or ["<p>No plot JSON files were produced.</p>"]) + "\n</body></html>\n"
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(doc)
    return out_path


def summarize(final_path, module_results, raw_path, raw_stats, gallery_path):
    import scanpy as sc
    from config import Config
    from models import Project

    adata = sc.read_h5ad(final_path)
    batch_counts = adata.obs["batch"].astype(str).value_counts().to_dict() if "batch" in adata.obs.columns else {}
    cluster_counts = adata.obs["leiden"].astype(str).value_counts().sort_index().to_dict() if "leiden" in adata.obs.columns else {}
    celltypes = adata.obs["celltype"].astype(str).value_counts().to_dict() if "celltype" in adata.obs.columns else {}

    rows = []
    batch_summary = None
    for module_name, task, result in module_results:
        summary = result.get("summary", {})
        if module_name == "batch_correct":
            batch_summary = summary
        rows.append(
            f"| `{module_name}` | `{task.status}` | `{task.output_adata_path}` | "
            f"`{json.dumps(summary, ensure_ascii=False)}` |"
        )

    report = f"""# GJ52 类器官与对应组织单细胞分析结果

生成时间：{datetime.now().isoformat(timespec='seconds')}

## 数据

- 项目 ID：`{PROJECT_ID}`
- 项目目录：`{Config.project_dir(PROJECT_ID)}`
- 原始合并 h5ad：`{raw_path}`
- 最终 h5ad：`{final_path}`
- 图表 HTML：`{gallery_path}`

## 输入批次

```json
{json.dumps(raw_stats, ensure_ascii=False, indent=2)}
```

## 流程

| 模块 | 状态 | 输出 h5ad | 摘要 |
|---|---|---|---|
{os.linesep.join(rows)}

## 最终数据概览

- 细胞数：{adata.n_obs}
- 基因数：{adata.n_vars}
- obs 列：`{', '.join(map(str, adata.obs.columns[:40]))}`
- obsm：`{', '.join(map(str, adata.obsm.keys()))}`

## 最终 batch 细胞数

```json
{json.dumps(batch_counts, ensure_ascii=False, indent=2)}
```

## Leiden cluster 细胞数

```json
{json.dumps(cluster_counts, ensure_ascii=False, indent=2)}
```

## 自动注释结果

```json
{json.dumps(celltypes, ensure_ascii=False, indent=2)}
```

## 批次效应整合说明

- 本次将 `52alg_singlet_matrix` 标记为 `batch=organoid`，将 `52tissue_singlet_matrix` 标记为 `batch=tissue`。
- 本次运行参数为 `method=combat`。当前环境中 OmicVerse ComBat 入口会触发 `unhashable type: 'list'`，因此平台批次校正模块自动回退到 Scanpy ComBat PCA，并在 `batch_correct` summary 的 `fallback` 字段记录。
- 下游聚类使用校正后的 embedding：`{(batch_summary or {}).get('embedding_key', 'unknown')}`。
- `organoid` 与 `tissue` 同时也是生物来源差异，不只是技术批次。批次校正后的 UMAP/cluster 用于共同嵌入和分群，不应单独作为消除全部生物差异的证据。
- 自动注释使用内置 `TME` marker set，只能作为初步参考；关键细胞类型应结合 marker heatmap、dotplot、DEG 和原始文献/实验背景复核。

## 使用说明

网页项目页：`/projects/{PROJECT_ID}`

离线总画廊：

`{gallery_path}`
"""
    with open(REPORT_PATH, "w", encoding="utf-8") as handle:
        handle.write(report)

    project = Project.get_by_id(PROJECT_ID)
    if project:
        project.status = "completed"
        project.metadata_json = json.dumps({
            "source_dirs": SOURCE_DIRS,
            "raw_path": raw_path,
            "final_adata_path": final_path,
            "gallery_path": gallery_path,
            "batch_counts": batch_counts,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        }, ensure_ascii=False)
        project.save()

    return {
        "report": os.path.abspath(REPORT_PATH),
        "gallery": gallery_path,
        "final_path": final_path,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "batch_counts": batch_counts,
        "clusters": cluster_counts,
        "celltypes": celltypes,
    }


def main():
    ensure_project()
    raw_path, raw_stats = build_merged_input()

    current_input = raw_path
    module_results = []
    all_result_files = []
    for module_name, params in PIPELINE:
        task, result = run_module(module_name, params, current_input)
        module_results.append((module_name, task, result))
        all_result_files.extend(result.get("result_files") or [])
        current_input = result["output_adata"]

    gallery_path = make_gallery(all_result_files)
    summary = summarize(current_input, module_results, raw_path, raw_stats, gallery_path)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Project URL: /projects/{PROJECT_ID}")
    print(f"Report: {os.path.abspath(REPORT_PATH)}")
    print(f"Gallery: {gallery_path}")


if __name__ == "__main__":
    main()
