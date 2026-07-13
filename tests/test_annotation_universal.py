"""End-to-end smoke coverage for universal first-pass annotation."""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def test_universal_annotation_emits_coverage_evidence(tmp_path):
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['Cd3d', 'Trac', 'Cd3e', 'Lck', 'Nkg7', 'Klrd1', 'Pecam1', 'Vwf',
             'Col1a1', 'Dcn'] + [f'Gene{i}' for i in range(70)]
    rng = np.random.default_rng(12)
    x = rng.poisson(1.0, size=(20, len(genes))).astype(float)
    x[:10, :4] += 5
    x[10:, 4:6] += 5
    adata = ad.AnnData(
        x,
        obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 10 + ['1'] * 10)}),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = rng.normal(size=(20, 2))
    input_path = tmp_path / 'input.h5ad'
    adata.write_h5ad(input_path)

    analysis = AnnotationAnalysis(
        str(tmp_path),
        {'method': 'auto_marker', 'marker_set': 'Universal', 'cluster_key': 'leiden',
         'min_markers_per_type': 2, 'show_celltype_composition': False,
         'show_marker_score_heatmap': False, 'show_marker_expression_violin': False,
         'show_annotation_score_umap': False},
        lambda *_: None,
    )
    result = analysis.run(str(input_path))

    assert result['summary']['marker_set'] == 'Universal'
    assert result['summary']['marker_coverage']['T cells']['matched'] >= 2
    assert any(os.path.basename(item['file_path']) == 'annotation_marker_coverage.json'
               for item in result['result_files'])


def test_multi_evidence_adds_hierarchical_and_audit_fields(tmp_path):
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK', 'NKG7', 'KLRD1'] + [f'G{i}' for i in range(60)]
    adata = ad.AnnData(np.random.default_rng(7).poisson(1, (12, len(genes))).astype(float),
                       obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 6 + ['1'] * 6)}),
                       var=pd.DataFrame(index=genes))
    adata.X[:6, :4] += 4
    adata.obsm['X_umap'] = np.random.default_rng(8).normal(size=(12, 2))
    path = tmp_path / 'multi.h5ad'
    adata.write_h5ad(path)
    result = AnnotationAnalysis(str(tmp_path), {
        'method': 'multi_evidence', 'marker_set': 'Universal', 'cluster_key': 'leiden',
        'min_markers_per_type': 2, 'cluster_agreement_threshold': 0.5,
        'show_celltype_composition': False, 'show_marker_score_heatmap': False,
        'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
    }, lambda *_: None).run(str(path))
    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert result['summary']['method_used'] == 'multi_evidence'
    assert {'marker_label', 'final_annotation', 'annotation_status', 'cell_type_l1',
            'cell_type_l2', 'cell_type_l3', 'cell_ontology_id'} <= set(output.obs.columns)
