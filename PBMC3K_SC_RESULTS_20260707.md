# PBMC3k 单细胞分析结果

生成时间：2026-07-07T13:40:14

## 数据

- 项目 ID：`pbmc3k_classic`
- 项目目录：`/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic`
- 原始数据：`/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/uploads/pbmc3k_raw.h5ad`
- 下载状态：`cached`
- 最终 h5ad：`/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/intermediate/deg_output.h5ad`
- 图表 HTML：`/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/results/pbmc3k_plot_gallery.html`

## 流程

| 模块 | 状态 | 输出 h5ad | 摘要 |
|---|---|---|---|
| `qc` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/intermediate/qc_output.h5ad` | `{"cells_before": 2700, "cells_after": 2661, "cells_removed": 39, "pct_removed": 1.4, "n_genes": 13714, "cells_removed_by_qc_and_doublet": 39, "novelty_median": 0.3669, "cell_cycle_available": true, "phase_counts": {"G1": 2487, "S": 135, "G2M": 39}, "s_genes_found": 40, "g2m_genes_found": 52}` |
| `normalize` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/intermediate/normalize_output.h5ad` | `{"n_cells": 2661, "n_genes": 13714, "method": "log1p", "target_sum": 10000.0}` |
| `hvg` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/intermediate/hvg_output.h5ad` | `{"n_cells": 2661, "n_genes_total": 13714, "n_hvgs": 1998, "hvg_flavor": "seurat_v3", "force_include_count": 0}` |
| `dimred` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/intermediate/dimred_output.h5ad` | `{"n_cells": 2661, "n_pcs": 50, "pca_variance_ratio_top5": 0.063, "embedding_method": "umap", "tsne_enabled": false}` |
| `clustering` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/intermediate/clustering_output.h5ad` | `{"n_clusters_0.8": 9, "n_clusters_0.4": 6, "n_clusters_1.2": 10, "resolutions": [0.8, 0.4, 1.2], "best_resolution": 0.8, "primary_resolution": 0.8}` |
| `annotation` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/intermediate/annotation_output.h5ad` | `{"n_celltypes": 9, "celltype_counts": {"CD4 T cells": 1149, "CD14+ Monocytes": 428, "NK cells": 356, "B cells": 328, "Unknown": 213, "FCGR3A+ Monocytes": 133, "Dendritic cells": 37, "Megakaryocytes": 12, "CD8 T cells": 5}, "cluster_column": "leiden_0.8", "method_used": "auto_marker", "mean_score_margin": 1.007}` |
| `deg` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/intermediate/deg_output.h5ad` | `{"n_groups": 9, "groups": ["0", "1", "2", "3", "4", "5", "6", "7", "8"], "method": "wilcoxon", "reference": "rest", "total_deg_genes": 180, "custom_dotplot_genes": null}` |

## 最终数据概览

- 细胞数：2661
- 基因数：13714
- obs 列：`n_genes_by_counts, log1p_n_genes_by_counts, total_counts, log1p_total_counts, pct_counts_in_top_20_genes, total_counts_mt, log1p_total_counts_mt, pct_counts_mt, total_counts_ribo, log1p_total_counts_ribo, pct_counts_ribo, total_counts_hb, log1p_total_counts_hb, pct_counts_hb, novelty_score, S_score, G2M_score, phase, nUMIs, mito_perc, ribo_perc, hb_perc, detected_genes, cell_complexity, n_counts, n_genes, passing_mt, passing_nUMIs, passing_ngenes, doublet_score`
- obsm：`X_pca, X_umap`

## Leiden cluster 细胞数

```json
{
  "0": 645,
  "1": 570,
  "2": 338,
  "3": 418,
  "4": 153,
  "5": 6,
  "6": 482,
  "7": 37,
  "8": 12
}
```

## 注释结果（如有）

```json
{
  "CD4 T cells": 1149,
  "CD14+ Monocytes": 428,
  "NK cells": 356,
  "B cells": 328,
  "Unknown": 213,
  "FCGR3A+ Monocytes": 133,
  "Dendritic cells": 37,
  "Megakaryocytes": 12,
  "CD8 T cells": 5
}
```

## Scanpy 教程对照与语义审计

- 官方 Scanpy PBMC3k 教程使用同一经典 10x PBMC3k 数据，原始规模为 `2700 x 32738`，基础过滤后为 `2700 x 13714`，按 `n_genes_by_counts < 2500`、`n_genes_by_counts > 200`、`pct_counts_mt < 5` 后为 `2638 x 13714`。
- 本平台流程不是逐行复刻官方教程：平台 QC 默认包含 Scrublet 双细胞处理、细胞周期评分和更宽松的 MT 阈值，因此最终为 `2661 x 13714`。报告不能声称与官方教程完全一致，只能声称使用同一数据并获得相近 PBMC 生物学结构。
- 官方教程 Leiden `resolution=0.7` 得到 8 个群；本平台主聚类使用 `primary_resolution=0.8`，得到 `9` 个群，粒度接近官方教程。
- 当前注释使用 `PBMC` marker set，覆盖官方教程中的 CD4 T、CD14+ Monocytes、B、CD8 T、NK、FCGR3A+ Monocytes、Dendritic 和 Megakaryocytes marker 结构。
- QC summary 中 `cells_removed_by_qc_and_doublet` 表示 QC 与 doublet 过滤后的综合减少量，不再声称是纯 doublet 数。
- annotation summary 中 `score_margin` 是 marker 分数差距，不是概率置信度；若使用该方法，报告字段应为 `mean_score_margin`。
- `Unknown` 代表 marker 分数不足或冲突的细胞，不应解读为新的细胞类型。

## 使用说明

打开项目页面后进入 `Classic PBMC3k scRNA` 项目即可查看任务与结果。图表也可以直接打开：

`/home/oelab/data/sc-transcriptomics-platform/data/projects/pbmc3k_classic/results/pbmc3k_plot_gallery.html`
