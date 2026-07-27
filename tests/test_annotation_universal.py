"""End-to-end smoke coverage for universal first-pass annotation."""
import os
import sys
import types

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
    assert 'marker_selection' in result['summary']
    assert result['summary']['marker_selection']['parameters']['method'] == 'wilcoxon'
    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert 'marker_selection' in output.uns
    assert any(os.path.basename(item['file_path']) == 'annotation_marker_coverage.png'
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
        'confidence_method': 'entropy', 'mark_unknown': True,
        'show_celltype_composition': False, 'show_marker_score_heatmap': False,
        'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
    }, lambda *_: None).run(str(path))
    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert result['summary']['method_used'] == 'multi_evidence'
    assert {'marker_label', 'final_annotation', 'annotation_status', 'cell_type_l1',
            'cell_type_l2', 'cell_type_l3', 'cell_ontology_id'} <= set(output.obs.columns)
    assert (output.obs['final_annotation'].astype(str) == output.obs['celltype'].astype(str)).all()
    # Low entropy alone must not erase a cluster whose marker label is coherent.
    assert (output.obs['celltype'].astype(str) != 'Unknown').all()
    assert any(os.path.basename(item['file_path']) == 'annotation_cluster_review.csv'
               for item in result['result_files'])


def test_cluster_review_marks_unknown_and_inconsistent_cluster_for_review():
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_annotation_cluster_review

    adata = ad.AnnData(
        np.ones((6, 2)),
        obs=pd.DataFrame({
            'leiden': ['0', '0', '0', '1', '1', '1'],
            'celltype': ['T', 'T', 'B', 'Unknown', 'Unknown', 'Unknown'],
            'final_annotation': ['T', 'T', 'B', 'Unknown', 'Unknown', 'Unknown'],
            'annotation_status': ['review', 'review', 'review', 'Unknown', 'Unknown', 'Unknown'],
            'annotation_score_margin': [0.2, 0.2, 0.2, 0.01, 0.01, 0.01],
        }),
        var=pd.DataFrame(index=['g1', 'g2']),
    )
    review = build_annotation_cluster_review(
        adata, 'leiden', {'cluster_markers': {'0': ['TRAC'], '1': ['NKG7']}},
    ).set_index('cluster')

    assert review.loc['0', 'needs_review']
    assert review.loc['1', 'needs_review']
    assert review.loc['1', 'top_markers'] == 'NKG7'


def test_marker_validation_prefers_counts_over_scaled_x():
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_marker_validation_matrix

    counts = np.asarray([
        [10, 0], [8, 0], [0, 12], [0, 9],
    ], dtype=float)
    adata = ad.AnnData(
        np.asarray([[-2, 0], [-1, 0], [0, -3], [0, -2]], dtype=float),
        obs=pd.DataFrame({'celltype': ['T', 'T', 'B', 'B']}),
        var=pd.DataFrame(index=['CD3D', 'MS4A1']),
    )
    adata.layers['counts'] = counts
    validation = build_marker_validation_matrix(adata, ['CD3D', 'MS4A1'])

    assert 'counts' in validation['expression_source']
    assert (validation['mean_expression'] >= 0).all()
    assert validation['categories'] == ['T', 'B']
    assert validation['detection_fraction'][0, 0] == 1.0


def test_celltypist_reference_is_optional_and_keeps_marker_labels(tmp_path, monkeypatch):
    """The reference adapter writes evidence without replacing Marker labels."""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules import annotation as annotation_module

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK', 'MS4A1', 'CD79A'] + [f'G{i}' for i in range(20)]
    counts = np.random.default_rng(101).poisson(1, (8, len(genes))).astype(float)
    counts[:4, :4] += 5
    adata = ad.AnnData(
        counts.copy(),
        obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 4 + ['1'] * 4)}),
        var=pd.DataFrame(index=genes),
    )
    adata.layers['counts'] = counts
    adata.obsm['X_umap'] = np.random.default_rng(102).normal(size=(8, 2))
    input_path = tmp_path / 'celltypist_reference.h5ad'
    adata.write_h5ad(input_path)

    model_path = tmp_path / 'fake_model.pkl'
    model_path.write_bytes(b'local test model')
    monkeypatch.setattr(annotation_module, 'CELLTYPIST_MODEL_DIR', str(tmp_path))

    fake_celltypist = types.ModuleType('celltypist')
    fake_celltypist.__version__ = 'test'

    class FakeResult:
        predicted_labels = pd.DataFrame({
            'predicted_labels': ['T cells'] * 4 + ['Unassigned'] * 4,
            'conf_score': [0.9] * 4 + [0.2] * 4,
        })
        probability_matrix = pd.DataFrame({
            'T cells': [0.9] * 4 + [0.2] * 4,
            'B cells': [0.1] * 4 + [0.1] * 4,
        })

    def fake_annotate(query, **kwargs):
        assert kwargs['model'] == str(model_path)
        assert kwargs['mode'] == 'prob match'
        assert query.n_obs == 8
        assert query.X.min() >= 0
        return FakeResult()

    fake_celltypist.annotate = fake_annotate
    monkeypatch.setitem(sys.modules, 'celltypist', fake_celltypist)

    result = annotation_module.AnnotationAnalysis(
        str(tmp_path),
        {
            'method': 'multi_evidence', 'marker_set': 'Universal', 'cluster_key': 'leiden',
            'min_markers_per_type': 2, 'use_celltypist_reference': True,
            'celltypist_model': 'fake_model.pkl', 'celltypist_mode': 'prob match',
            'celltypist_p_threshold': 0.5,
            'show_celltype_composition': False, 'show_marker_score_heatmap': False,
            'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
        },
        lambda *_: None,
    ).run(str(input_path))

    assert result['summary']['celltypist_reference']['enabled'] is True
    assert result['summary']['celltypist_reference']['model_name'] == 'fake_model.pkl'
    assert result['summary']['celltypist_reference']['status_counts']['unassigned'] == 4
    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert {'celltypist_label', 'celltypist_confidence', 'celltypist_comparison'} <= set(output.obs.columns)
    assert 'celltype' in output.obs.columns
    assert (output.obs['final_annotation'].astype(str) == output.obs['celltype'].astype(str)).all()


def test_organoid_annotation_uses_selected_tissue_panel(tmp_path):
    """Organoid + organoid_type must reach scoring and coverage evidence."""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['NPHS1', 'NPHS2', 'PODXL', 'WT1', 'SIX2', 'CITED1'] + [f'Gene{i}' for i in range(30)]
    rng = np.random.default_rng(32)
    x = rng.poisson(0.2, size=(12, len(genes))).astype(float)
    x[:6, :4] += 7
    x[6:, 4:6] += 5
    adata = ad.AnnData(
        x,
        obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 6 + ['1'] * 6)}),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = rng.normal(size=(12, 2))
    input_path = tmp_path / 'kidney_organoid.h5ad'
    adata.write_h5ad(input_path)

    result = AnnotationAnalysis(
        str(tmp_path),
        {'method': 'auto_marker', 'marker_set': 'Organoid', 'organoid_type': 'kidney',
         'cluster_key': 'leiden', 'min_markers_per_type': 2,
         'show_celltype_composition': False, 'show_marker_score_heatmap': False,
         'show_marker_expression_violin': False, 'show_annotation_score_umap': False},
        lambda *_: None,
    ).run(str(input_path))

    assert result['summary']['marker_set'] == 'Organoid'
    assert result['summary']['organoid_type'] == 'kidney'
    assert result['summary']['marker_coverage']['Podocytes']['matched'] == 4
    assert result['summary']['maturity_evidence'] is not None
    assert result['summary']['annotation_version'] == 'v1'
    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert {'cell_lineage', 'cell_type_l1', 'cell_type_l2', 'cell_type_l3',
            'cell_state', 'cell_state_flags', 'developmental_state'} <= set(output.obs.columns)


def test_organoid_maturity_evidence_and_quality_evidence_are_auditable():
    ad = pytest.importorskip('anndata')
    from modules.annotation import (
        compute_annotation_quality_evidence,
        compute_organoid_maturity_evidence,
    )

    maturity_adata = ad.AnnData(
        np.asarray([[8, 0, 0, 0], [7, 0, 0, 0], [0, 8, 0, 0], [0, 7, 0, 0]], dtype=float),
        obs=pd.DataFrame({'culture_day': [3, 3, 21, 21]}),
        var=pd.DataFrame(index=['LGR5', 'ALPI', 'MKI67', 'VIL1']),
    )
    maturity = compute_organoid_maturity_evidence(
        maturity_adata, 'intestinal', time_key='culture_day'
    )
    assert maturity['time_key'] == 'culture_day'
    assert set(maturity['maturity_state']) >= {'progenitor_like', 'mature_like'}
    assert maturity['marker_coverage']['progenitor']['matched'] == 1
    assert maturity['marker_coverage']['mature']['matched'] == 2

    adata = ad.AnnData(
        np.asarray([
            [4, 0, 0, 3], [3, 0, 0, 2],  # balanced epithelial / neural evidence
            [4, 0, 0, 0], [3, 0, 0, 0],
        ], dtype=float),
        obs=pd.DataFrame({
            'score_Epithelial': [1.0, 0.9, 1.2, 1.1],
            'score_Neuron': [0.95, 0.85, -0.2, -0.1],
        }),
        var=pd.DataFrame(index=['EPCAM', 'TUBB3', 'KRT8', 'LYZ']),
    )
    evidence = compute_annotation_quality_evidence(
        adata,
        ['score_Epithelial', 'score_Neuron'],
        {'Epithelial': ['EPCAM', 'KRT8'], 'Neuron': ['TUBB3']},
        doublet_threshold=0.30,
        ambient_threshold=0.35,
        ambient_prevalence=0.5,
    )
    assert evidence['doublet_status'][0] == 'suspect_doublet'
    assert evidence['top1'][0] == 'Epithelial'
    assert evidence['top2'][0] == 'Neuron'
    assert 'EPCAM' in evidence['ambient_genes']

    single = compute_annotation_quality_evidence(
        adata, ['score_Epithelial'], {'Epithelial': ['EPCAM']}
    )
    assert set(single['doublet_status']) == {'not_evaluated'}
