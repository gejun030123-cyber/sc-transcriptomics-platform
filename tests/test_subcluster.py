"""Focused tests for the targeted single-cell subcluster module."""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def _analysis(tmp_path, **params):
    from modules.subcluster import SubclusterAnalysis
    defaults = {'source_cluster_key': 'leiden', 'target_cluster': '0', 'min_cells': 5,
                'n_neighbors': 4, 'resolution': 0.4, 'run_enrichment': False,
                'n_genes': 8, 'marker_heatmap_top_n': 2}
    defaults.update(params)
    return SubclusterAnalysis(str(tmp_path), defaults, lambda *_: None)


def _adata():
    ad = pytest.importorskip('anndata')
    rng = np.random.default_rng(3)
    x = rng.poisson(1.2, size=(24, 16)).astype(float)
    # Give the first parent cluster two detectable sub-populations.
    x[:6, :3] += 5
    x[6:12, 3:6] += 5
    adata = ad.AnnData(x, obs=pd.DataFrame({'leiden': ['0'] * 12 + ['1'] * 12}),
                       var=pd.DataFrame(index=[f'G{i}' for i in range(16)]))
    import scanpy as sc
    sc.pp.pca(adata)
    return adata


def test_validate_input_requires_existing_target(tmp_path):
    adata = _adata()
    assert _analysis(tmp_path).validate_input(adata) is None
    assert 'does not exist' in _analysis(tmp_path, target_cluster='99').validate_input(adata)
    assert 'not found' in _analysis(tmp_path, source_cluster_key='missing').validate_input(adata)


def test_subcluster_run_writes_reusable_outputs(tmp_path):
    pytest.importorskip('scanpy')
    adata = _adata()
    input_path = tmp_path / 'input.h5ad'
    adata.write_h5ad(input_path)
    result = _analysis(tmp_path).run(str(input_path))

    assert os.path.exists(result['output_adata'])
    assert result['summary']['n_cells'] == 12
    assert result['summary']['n_subclusters'] >= 1
    result_names = {os.path.basename(item['file_path']) for item in result['result_files']}
    assert {'subcluster_umap.json', 'subcluster_size_bar.json', 'subcluster_deg_results.csv'} <= result_names
    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert {'parent_cluster', 'subcluster'} <= set(output.obs.columns)
