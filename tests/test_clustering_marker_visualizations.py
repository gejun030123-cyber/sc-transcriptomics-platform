"""Regression coverage for cluster marker heatmap and violin exports."""

import os

import numpy as np
import pandas as pd
import pytest


def test_marker_preview_exports_heatmap_and_one_violin_per_cluster(tmp_path):
    pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    import anndata
    import scanpy as sc

    from modules.clustering import ClusteringAnalysis

    rng = np.random.default_rng(11)
    genes = ['CL0_MARK', 'CL1_MARK', 'CL2_MARK'] + [f'G{index}' for index in range(9)]
    counts = rng.poisson(0.15, size=(72, len(genes))).astype(float)
    for cluster_index in range(3):
        start = cluster_index * 24
        counts[start:start + 24, cluster_index] += 8
    adata = anndata.AnnData(
        counts.copy(),
        obs=pd.DataFrame({
            'leiden_0.5': pd.Categorical(['0'] * 24 + ['1'] * 24 + ['2'] * 24),
        }),
        var=pd.DataFrame(index=genes),
    )
    adata.layers['counts'] = counts.copy()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    analysis = ClusteringAnalysis(
        str(tmp_path),
        {
            'compute_marker_preview': True,
            'marker_min_per_cluster': 1,
            'marker_max_per_cluster': 2,
            'marker_heatmap_top_n': 1,
            'marker_violin_top_n': 1,
            'marker_violin_max_cells_per_group': 50,
        },
        lambda *_: None,
    )
    result_files = []
    selection = analysis._run_marker_preview(
        adata, 'leiden_0.5', analysis.ensure_plots_dir(), result_files,
    )

    png_names = {os.path.basename(item['file_path']) for item in result_files
                 if item['file_type'] == 'png'}
    assert 'cluster_marker_dotplot.png' in png_names
    assert 'cluster_marker_heatmap.png' in png_names
    violin_pngs = [name for name in png_names if name.startswith('cluster_marker_violin_')]
    assert len(violin_pngs) == 3
    assert selection['visualizations']['heatmap'] is True
    assert selection['visualizations']['heatmap_n_genes'] == 3
    assert selection['visualizations']['violin_clusters'] == ['0', '1', '2']
    assert selection['visualizations']['violin_expression_source'].startswith('layers["counts"]')

