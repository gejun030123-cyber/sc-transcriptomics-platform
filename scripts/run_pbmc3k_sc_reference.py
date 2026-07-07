#!/usr/bin/env python3
"""Download classic PBMC3k and run the platform scRNA core workflow."""

import json
import os
import sys
import traceback
from datetime import datetime


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

PROJECT_ID = "pbmc3k_classic"
PROJECT_NAME = "Classic PBMC3k scRNA"
REPORT_PATH = "PBMC3K_SC_RESULTS_20260707.md"
GALLERY_NAME = "pbmc3k_plot_gallery.html"


PIPELINE = [
    ("qc", {
        "mito_perc": 0.20,
        "nUMIs": 500,
        "detected_genes": 200,
        "max_detected_genes": 6000,
        "ribo_perc": 0,
        "hb_perc": 0,
        "batch_key": "",
        "save_counts_layer": True,
    }),
    ("normalize", {
        "method": "log1p",
        "target_sum": 10000,
    }),
    ("hvg", {
        "n_top_genes": 2000,
        "hvg_flavor": "seurat_v3",
        "batch_key": "",
        "exclude_mt_genes": True,
    }),
    ("dimred", {
        "n_comps": 50,
        "umap_n_neighbors": 15,
        "umap_min_dist": 0.5,
        "umap_metric": "euclidean",
        "enable_tsne": False,
    }),
    ("clustering", {
        "resolutions": "0.8,0.4,1.2",
        "primary_resolution": 0.8,
        "n_neighbors": 15,
        "clustering_method": "leiden",
        "n_iterations": 2,
        "distance_metric": "euclidean",
        "use_corrected": False,
    }),
    ("annotation", {
        "method": "auto_marker",
        "cluster_key": "leiden",
        "resolution": "0.8",
        "marker_set": "PBMC",
        "confidence_method": "score_margin",
        "mark_unknown": True,
    }),
    ("deg", {
        "groupby": "leiden",
        "reference": "rest",
        "method": "wilcoxon",
        "n_genes": 20,
        "show_dotplot": True,
        "pval_cutoff": 0.05,
        "logfc_cutoff": 0.25,
        "min_pct": 0.1,
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
        description="Classic 10x PBMC3k dataset downloaded via scanpy.datasets.pbmc3k()",
        status="processing",
        metadata_json=json.dumps({"source": "scanpy.datasets.pbmc3k"}, ensure_ascii=False),
    )
    project.save()
    return project


def download_pbmc3k():
    from config import Config
    import scanpy as sc

    raw_path = os.path.join(Config.uploads_dir(PROJECT_ID), "pbmc3k_raw.h5ad")
    if os.path.isfile(raw_path):
        return raw_path, "cached"

    sc.settings.datasetdir = os.path.join(Config.project_dir(PROJECT_ID), "download_cache")
    os.makedirs(sc.settings.datasetdir, exist_ok=True)
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    adata.write_h5ad(raw_path)
    return raw_path, "downloaded"


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
        tb = traceback.format_exc()
        task.mark_failed(tb)
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
            sections.append(
                f"<section><h2>{rf.get('label', os.path.basename(fpath))}</h2>{html}</section>"
            )
        except Exception as exc:
            sections.append(
                f"<section><h2>{rf.get('label', os.path.basename(fpath))}</h2>"
                f"<pre>Could not render {fpath}: {exc}</pre></section>"
            )

    doc = """<!doctype html>
<html><head><meta charset="utf-8"><title>PBMC3k scRNA plot gallery</title>
<style>
body{font-family:Arial,sans-serif;margin:24px;background:#f8f9fb;color:#1f2933}
section{background:#fff;border:1px solid #dde3ea;border-radius:8px;margin:0 0 20px;padding:16px}
h1{margin:0 0 8px} h2{font-size:16px;margin:0 0 12px}
.meta{color:#667085;margin-bottom:20px}
</style></head><body>
<h1>PBMC3k scRNA Plot Gallery</h1>
<div class="meta">Generated by scripts/run_pbmc3k_sc_reference.py</div>
""" + "\n".join(sections or ["<p>No plot JSON files were produced.</p>"]) + "\n</body></html>\n"
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(doc)
    return out_path


def summarize(final_path, module_results, raw_path, download_status, gallery_path):
    import scanpy as sc
    from config import Config
    from models import Project

    adata = sc.read_h5ad(final_path)
    cluster_counts = {}
    if "leiden" in adata.obs.columns:
        cluster_counts = adata.obs["leiden"].astype(str).value_counts().sort_index().to_dict()
    celltypes = {}
    if "celltype" in adata.obs.columns:
        celltypes = adata.obs["celltype"].astype(str).value_counts().head(20).to_dict()

    rows = []
    for module_name, task, result in module_results:
        rows.append(
            f"| `{module_name}` | `{task.status}` | `{task.output_adata_path}` | "
            f"`{json.dumps(result.get('summary', {}), ensure_ascii=False)}` |"
        )

    report = f"""# PBMC3k 单细胞分析结果

生成时间：{datetime.now().isoformat(timespec='seconds')}

## 数据

- 项目 ID：`{PROJECT_ID}`
- 项目目录：`{Config.project_dir(PROJECT_ID)}`
- 原始数据：`{raw_path}`
- 下载状态：`{download_status}`
- 最终 h5ad：`{final_path}`
- 图表 HTML：`{gallery_path}`

## 流程

| 模块 | 状态 | 输出 h5ad | 摘要 |
|---|---|---|---|
{os.linesep.join(rows)}

## 最终数据概览

- 细胞数：{adata.n_obs}
- 基因数：{adata.n_vars}
- obs 列：`{', '.join(map(str, adata.obs.columns[:30]))}`
- obsm：`{', '.join(map(str, adata.obsm.keys()))}`

## Leiden cluster 细胞数

```json
{json.dumps(cluster_counts, ensure_ascii=False, indent=2)}
```

## 注释结果（如有）

```json
{json.dumps(celltypes, ensure_ascii=False, indent=2)}
```

## Scanpy 教程对照与语义审计

- 官方 Scanpy PBMC3k 教程使用同一经典 10x PBMC3k 数据，原始规模为 `2700 x 32738`，基础过滤后为 `2700 x 13714`，按 `n_genes_by_counts < 2500`、`n_genes_by_counts > 200`、`pct_counts_mt < 5` 后为 `2638 x 13714`。
- 本平台流程不是逐行复刻官方教程：平台 QC 默认包含 Scrublet 双细胞处理、细胞周期评分和更宽松的 MT 阈值，因此最终为 `{adata.n_obs} x {adata.n_vars}`。报告不能声称与官方教程完全一致，只能声称使用同一数据并获得相近 PBMC 生物学结构。
- 官方教程 Leiden `resolution=0.7` 得到 8 个群；本平台主聚类使用 `primary_resolution=0.8`，得到 `{len(cluster_counts)}` 个群，粒度接近官方教程。
- 当前注释使用 `PBMC` marker set，覆盖官方教程中的 CD4 T、CD14+ Monocytes、B、CD8 T、NK、FCGR3A+ Monocytes、Dendritic 和 Megakaryocytes marker 结构。
- QC summary 中 `cells_removed_by_qc_and_doublet` 表示 QC 与 doublet 过滤后的综合减少量，不再声称是纯 doublet 数。
- annotation summary 中 `score_margin` 是 marker 分数差距，不是概率置信度；若使用该方法，报告字段应为 `mean_score_margin`。
- `Unknown` 代表 marker 分数不足或冲突的细胞，不应解读为新的细胞类型。

## 使用说明

打开项目页面后进入 `Classic PBMC3k scRNA` 项目即可查看任务与结果。图表也可以直接打开：

`{gallery_path}`
"""
    with open(REPORT_PATH, "w", encoding="utf-8") as handle:
        handle.write(report)

    project = Project.get_by_id(PROJECT_ID)
    if project:
        project.status = "completed"
        project.metadata_json = json.dumps({
            "source": "scanpy.datasets.pbmc3k",
            "raw_path": raw_path,
            "final_adata_path": final_path,
            "gallery_path": gallery_path,
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
        "clusters": cluster_counts,
        "celltypes": celltypes,
    }


def main():
    from config import Config

    ensure_project()
    raw_path, download_status = download_pbmc3k()

    current_input = raw_path
    module_results = []
    all_result_files = []
    for module_name, params in PIPELINE:
        task, result = run_module(module_name, params, current_input)
        module_results.append((module_name, task, result))
        all_result_files.extend(result.get("result_files") or [])
        current_input = result["output_adata"]

    gallery_path = make_gallery(all_result_files)
    summary = summarize(current_input, module_results, raw_path, download_status, gallery_path)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Project URL: /projects/{PROJECT_ID}")
    print(f"Report: {os.path.abspath(REPORT_PATH)}")
    print(f"Gallery: {gallery_path}")


if __name__ == "__main__":
    main()
