# GJ52 batch method comparison

- Input: `/home/oelab/data/sc-transcriptomics-platform/data/projects/gj52_organoid_tissue_batch/intermediate/hvg_output.h5ad`
- Batch key: `batch`
- Resolution: `0.8`
- Gallery: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/gj52_batch_method_comparison.html`
- Metrics CSV: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/batch_method_metrics.csv`

## Completed methods

| n_cells | n_clusters | asw_batch | weighted_batch_entropy | weighted_max_batch_fraction | min_cluster_size | median_cluster_size | small_clusters_lt_50 | asw_reference_celltype | method | embedding_key | runtime_sec |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 26294 | 28 | 0.122 | 0.0556 | 0.9872 | 26 | 501.0 | 1 | 0.1276 | baseline | X_pca | 53.1 |
| 26294 | 27 | 0.0156 | 0.0802 | 0.9833 | 26 | 648.0 | 1 | 0.093 | combat | X_pca_combat | 38.8 |
| 26294 | 27 | 0.109 | 0.2259 | 0.9296 | 25 | 501.0 | 1 | 0.1298 | harmony | X_pca_harmony | 25.8 |
| 26294 | 31 | 0.122 | 0.3084 | 0.915 | 6 | 626.0 | 4 | 0.1276 | bbknn | X_pca | 17.3 |
| 26294 | 18 | 0.0804 | 0.3634 | 0.8934 | 87 | 1387.0 | 0 | 0.0978 | scanorama | X_scanorama | 31.9 |
| 26294 | 25 | 0.0417 | 0.6455 | 0.7745 | 111 | 1178.0 | 0 | 0.2134 | sysvi | X_sysvi | 190.9 |

## Skipped or failed methods

None.

## Interpretation notes

- Lower absolute batch ASW suggests stronger batch mixing, but for organoid vs tissue this can also mean biological differences were compressed.
- Higher weighted batch entropy means clusters contain both organoid and tissue more evenly.
- Lower weighted max batch fraction means fewer clusters are dominated by a single source.
- Reference celltype labels come from the previous pipeline annotation and are used only as a visual guide.

## Output files

- `baseline` `batch`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/baseline_umap_batch.json`
- `baseline` `cluster`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/baseline_umap_cluster.json`
- `baseline` `reference_celltype`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/baseline_umap_reference_celltype.json`
- `baseline` `h5ad`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/baseline_comparison.h5ad`
- `combat` `batch`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/combat_umap_batch.json`
- `combat` `cluster`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/combat_umap_cluster.json`
- `combat` `reference_celltype`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/combat_umap_reference_celltype.json`
- `combat` `h5ad`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/combat_comparison.h5ad`
- `harmony` `batch`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/harmony_umap_batch.json`
- `harmony` `cluster`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/harmony_umap_cluster.json`
- `harmony` `reference_celltype`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/harmony_umap_reference_celltype.json`
- `harmony` `h5ad`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/harmony_comparison.h5ad`
- `bbknn` `batch`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/bbknn_umap_batch.json`
- `bbknn` `cluster`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/bbknn_umap_cluster.json`
- `bbknn` `reference_celltype`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/bbknn_umap_reference_celltype.json`
- `bbknn` `h5ad`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/bbknn_comparison.h5ad`
- `scanorama` `batch`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/scanorama_umap_batch.json`
- `scanorama` `cluster`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/scanorama_umap_cluster.json`
- `scanorama` `reference_celltype`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/scanorama_umap_reference_celltype.json`
- `scanorama` `h5ad`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/scanorama_comparison.h5ad`
- `sysvi` `batch`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/sysvi_umap_batch.json`
- `sysvi` `cluster`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/sysvi_umap_cluster.json`
- `sysvi` `reference_celltype`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/sysvi_umap_reference_celltype.json`
- `sysvi` `h5ad`: `data/projects/gj52_organoid_tissue_batch/results/batch_method_comparison_20260707_v3_sysvi/sysvi_comparison.h5ad`
