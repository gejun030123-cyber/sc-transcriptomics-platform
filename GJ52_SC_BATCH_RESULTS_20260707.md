# GJ52 类器官与对应组织单细胞分析结果

生成时间：2026-07-07T18:02:34

## 数据

- 项目 ID：`gj52_organoid_tissue_batch`
- 项目目录：`/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch`
- 原始合并 h5ad：`/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/uploads/gj52_organoid_tissue_raw.h5ad`
- 最终 h5ad：`/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/proportion_output.h5ad`
- 图表 HTML：`/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/results/gj52_plot_gallery.html`

## 输入批次

```json
[
  {
    "batch": "organoid",
    "sample": "52alg_singlet_matrix",
    "path": "/home/oelab/data/GJ/000000/52alg_singlet_matrix",
    "n_cells": 17435,
    "n_genes": 16572
  },
  {
    "batch": "tissue",
    "sample": "52tissue_singlet_matrix",
    "path": "/home/oelab/data/GJ/000000/52tissue_singlet_matrix",
    "n_cells": 13877,
    "n_genes": 16921
  }
]
```

## 流程

| 模块 | 状态 | 输出 h5ad | 摘要 |
|---|---|---|---|
| `qc` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/qc_output.h5ad` | `{"cells_before": 31312, "cells_after": 26294, "cells_removed": 5018, "pct_removed": 16.0, "n_genes": 16937, "cells_removed_by_qc_and_doublet": 386, "novelty_median": 0.6519, "cell_cycle_available": true, "phase_counts": {"G1": 14069, "G2M": 7199, "S": 5026}, "s_genes_found": 42, "g2m_genes_found": 53}` |
| `normalize` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/normalize_output.h5ad` | `{"n_cells": 26294, "n_genes": 15147, "method": "log1p", "target_sum": 10000.0}` |
| `hvg` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/hvg_output.h5ad` | `{"n_cells": 26294, "n_genes_total": 15147, "n_hvgs": 4342, "hvg_flavor": "seurat_v3", "force_include_count": 0}` |
| `dimred` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/dimred_output.h5ad` | `{"n_cells": 26294, "n_pcs": 50, "pca_variance_ratio_top5": 0.143, "embedding_method": "umap", "tsne_enabled": false}` |
| `batch_correct` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/batch_correct_output.h5ad` | `{"method": "combat", "n_cells": 26294, "embedding_key": "X_pca_combat", "fallback": {"requested_method": "combat", "method_used": "scanpy_combat_pca", "reason": "unhashable type: 'list'", "n_hvg_for_combat": 4342}, "eval_metrics": {"asw_batch": 0.0139, "graph_components": 1}}` |
| `clustering` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/clustering_output.h5ad` | `{"n_clusters_0.4": 20, "n_clusters_0.6": 26, "n_clusters_0.8": 33, "n_clusters_1.0": 35, "n_clusters_1.2": 41, "resolutions": [0.4, 0.6, 0.8, 1.0, 1.2], "best_resolution": 0.8, "primary_resolution": 0.8}` |
| `qc_reassess` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/qc_reassess_output.h5ad` | `{"n_clusters": 33, "n_low_quality": 0, "low_quality_clusters": [], "doublet_threshold": 0.3, "mt_threshold": 15.0}` |
| `annotation` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/annotation_output.h5ad` | `{"n_celltypes": 9, "celltype_counts": {"Tumor Epithelial": 14905, "Unknown": 4613, "T cells": 2624, "CAF": 1731, "Plasma cells": 967, "Monocyte/Macrophage": 801, "Endothelial": 350, "NK cells": 157, "B cells": 146}, "cluster_column": "leiden_0.8", "method_used": "auto_marker", "mean_score_margin": 0.861}` |
| `deg` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/deg_output.h5ad` | `{"n_groups": 33, "groups": ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32"], "method": "wilcoxon", "reference": "rest", "total_deg_genes": 660, "custom_dotplot_genes": null}` |
| `proportion` | `completed` | `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/proportion_output.h5ad` | `{"chi2": 6275.35, "p_value": 0.0, "n_batches": 2, "n_groups": 9, "stat_test": "chi_square"}` |

## 最终数据概览

- 细胞数：26294
- 基因数：15147
- obs 列：`barcode, batch, sample, source_type, source_description, n_genes_by_counts, log1p_n_genes_by_counts, total_counts, log1p_total_counts, pct_counts_in_top_20_genes, total_counts_mt, log1p_total_counts_mt, pct_counts_mt, total_counts_ribo, log1p_total_counts_ribo, pct_counts_ribo, total_counts_hb, log1p_total_counts_hb, pct_counts_hb, novelty_score, S_score, G2M_score, phase, nUMIs, mito_perc, ribo_perc, hb_perc, detected_genes, cell_complexity, n_counts, n_genes, passing_mt, passing_nUMIs, passing_ngenes, doublet_score, predicted_doublet, leiden_0.4, leiden_0.6, leiden_0.8, leiden_1.0`
- obsm：`X_pca, X_pca_combat, X_umap`

## 最终 batch 细胞数

```json
{
  "tissue": 13501,
  "organoid": 12793
}
```

## Leiden cluster 细胞数

```json
{
  "0": 149,
  "1": 2112,
  "10": 703,
  "11": 623,
  "12": 21,
  "13": 914,
  "14": 356,
  "15": 196,
  "16": 1693,
  "17": 1172,
  "18": 40,
  "19": 2938,
  "2": 19,
  "20": 183,
  "21": 3385,
  "22": 222,
  "23": 344,
  "24": 175,
  "25": 835,
  "26": 987,
  "27": 506,
  "28": 58,
  "29": 1407,
  "3": 968,
  "30": 360,
  "31": 650,
  "32": 1202,
  "4": 2189,
  "5": 25,
  "6": 565,
  "7": 207,
  "8": 504,
  "9": 586
}
```

## 自动注释结果

```json
{
  "Tumor Epithelial": 14905,
  "Unknown": 4613,
  "T cells": 2624,
  "CAF": 1731,
  "Plasma cells": 967,
  "Monocyte/Macrophage": 801,
  "Endothelial": 350,
  "NK cells": 157,
  "B cells": 146
}
```

## 批次效应整合说明

- 本次将 `52alg_singlet_matrix` 标记为 `batch=organoid`，将 `52tissue_singlet_matrix` 标记为 `batch=tissue`。
- 本次运行参数为 `method=combat`。当前环境中 OmicVerse ComBat 入口会触发 `unhashable type: 'list'`，因此平台批次校正模块自动回退到 Scanpy ComBat PCA，并在 `batch_correct` summary 的 `fallback` 字段记录。
- 下游聚类使用校正后的 embedding：`X_pca_combat`。
- `organoid` 与 `tissue` 同时也是生物来源差异，不只是技术批次。批次校正后的 UMAP/cluster 用于共同嵌入和分群，不应单独作为消除全部生物差异的证据。
- 自动注释使用内置 `TME` marker set，只能作为初步参考；关键细胞类型应结合 marker heatmap、dotplot、DEG 和原始文献/实验背景复核。

## 使用说明

网页项目页：`/projects/gj52_organoid_tissue_batch`

离线总画廊：

`/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/results/gj52_plot_gallery.html`
