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


def test_auto_marker_panel_uses_colorectal_only_with_expression_evidence():
    """不能仅凭 10x 全基因清单切换到特异组织面板。"""
    ad = pytest.importorskip('anndata')
    from modules.annotation import resolve_auto_marker_set

    intestinal = ad.AnnData(
        np.asarray([
            [5, 4, 3, 2, 4, 3],
            [4, 3, 2, 3, 5, 4],
            [3, 4, 4, 3, 3, 2],
            [4, 5, 3, 4, 2, 3],
        ], dtype=float),
        var=pd.DataFrame(index=['ALPI', 'VIL1', 'MUC2', 'TFF3', 'MKI67', 'TOP2A']),
    )
    panel, evidence = resolve_auto_marker_set(intestinal, 'Auto')

    assert panel == 'Colorectal'
    assert evidence['n_supported_programmes'] >= 3
    assert evidence['programme_evidence']['Goblet-like epithelial']['detected_genes'] == 2

    broad_epithelial = ad.AnnData(
        np.ones((4, 3), dtype=float),
        var=pd.DataFrame(index=['EPCAM', 'KRT8', 'KRT18']),
    )
    panel, evidence = resolve_auto_marker_set(broad_epithelial, 'Auto')

    assert panel == 'Universal'
    assert evidence['n_supported_programmes'] == 0


def test_refined_colorectal_panel_preserves_confirmed_label_contract():
    """已确认的本地精细标签必须有稳定 Marker 与层级字段。"""
    from modules.annotation import (
        COLORECTAL_REFINED_MARKERS,
        REFINED_PANEL_VERSION,
        get_marker_set,
        hierarchy_for_label,
    )

    panel = get_marker_set('Colorectal_refined')
    assert panel == COLORECTAL_REFINED_MARKERS
    # v2 面板：17 个可独立判别的 programme，而不是 11 个宽标签。
    assert REFINED_PANEL_VERSION == 'colorectal_refined_v2'
    assert set(panel) == {
        'CYP3A5+ enterocytes',
        'KRT20/MALRD1+ absorptive enterocytes',
        'HSD17B2+ absorptive enterocytes',
        'MTTP/RBP2+ absorptive enterocytes',
        'CKB/FABP1+ metabolic enterocytes',
        'NDRG1/ANKRD37+ hypoxia-response enterocytes',
        'MAML3/CHRM3+ epithelial cells (review)',
        'TFF1/REG4+ secretory cells',
        'MUC2/FCGBP+ Goblet cells',
        'SPINK4/CA4+ secretory cells',
        'DDIT3/ATF3+ ER-stress enterocytes',
        'TA/stem-like cells',
        'S-phase enterocytes',
        'MIR924HG/DIAPH3+ cycling cells',
        'Histone-rich cycling enterocytes',
        'BIRC5/CDKN3+ cycling enterocytes',
        'CDC20/CCNB1+ cycling enterocytes',
    }
    # 别名 Colorectal_refined_v2 指向同一个稳定面板。
    assert get_marker_set('Colorectal_refined_v2') == COLORECTAL_REFINED_MARKERS
    assert hierarchy_for_label(
        'CYP3A5+ enterocytes', 'Colorectal_refined'
    )[:3] == ('Epithelial lineage', 'Absorptive epithelium', 'CYP3A5+ enterocytes')
    assert hierarchy_for_label(
        'HSD17B2+ absorptive enterocytes', 'Colorectal_refined'
    )[:3] == (
        'Epithelial lineage', 'Absorptive epithelium',
        'HSD17B2+ absorptive enterocytes',
    )
    assert hierarchy_for_label(
        'Unresolved epithelial programme: GENE1/GENE2 (cluster 3; review)',
        'Colorectal_refined',
    )[:3] == (
        'Epithelial lineage', 'Unresolved epithelial programme',
        'Unresolved epithelial programme: GENE1/GENE2 (cluster 3; review)',
    )


def test_refined_display_aliases_preserve_marker_programme_provenance():
    ad = pytest.importorskip('anndata')
    from modules.annotation import (
        apply_curated_annotation_display_labels,
        hierarchy_for_label,
    )

    programmes = [
        'SPINK4/CA4+ secretory cells',
        'CKB/FABP1+ metabolic enterocytes',
        'Unresolved epithelial programme: CLDN2/NME1 (cluster 3; review)',
        'KRT20/MALRD1+ absorptive enterocytes',
    ]
    adata = ad.AnnData(
        np.ones((len(programmes), 1)),
        obs=pd.DataFrame({'celltype': programmes}),
        var=pd.DataFrame(index=['GENE']),
    )

    changed = apply_curated_annotation_display_labels(adata, 'Colorectal_refined')

    assert changed == 3
    assert adata.obs['celltype'].astype(str).tolist() == [
        'Goblet cells', 'High metabolic cells', 'Inflammatory stress cells',
        'KRT20/MALRD1+ absorptive enterocytes',
    ]
    assert adata.obs['annotation_marker_programme'].astype(str).tolist() == programmes
    assert adata.obs['annotation_display_label_source'].astype(str).tolist() == [
        'curated_current_run_map', 'curated_current_run_map',
        'curated_current_run_map', 'marker_programme',
    ]
    assert hierarchy_for_label(
        'High metabolic cells', 'Colorectal_refined',
    )[:3] == ('Epithelial lineage', 'Absorptive epithelium', 'High metabolic cells')


def test_refined_colorectal_prefers_cluster_specific_signature_over_shared_score():
    """共享高表达模块不能把有明确 Top marker 的精细簇合并掉。"""
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_marker_decisions

    genes = ['CKB', 'FABP1', 'CYP3A5', 'MAML3', 'HSD17B2', 'LAMA1']
    x = np.asarray([
        [8, 7, 0, 0, 0, 0], [7, 8, 0, 0, 0, 0],
        [0, 0, 8, 7, 0, 0], [0, 0, 7, 8, 0, 0],
        [0, 0, 0, 0, 8, 7], [0, 0, 0, 0, 7, 8],
    ], dtype=float)
    labels = ['0', '0', '1', '1', '2', '2']
    obs = pd.DataFrame({'leiden': labels})
    # Simulate a shared high-metabolic score that would otherwise win every
    # cluster.  The independently ranked cluster signatures must take priority
    # for the explicitly local refined panel.
    for label, values in {
        'CKB/FABP1+ metabolic enterocytes': [10.0] * 6,
        'CYP3A5+ enterocytes': [0.2, 0.2, 1.0, 1.0, 0.2, 0.2],
        'HSD17B2+ absorptive enterocytes': [0.2, 0.2, 0.2, 0.2, 1.0, 1.0],
    }.items():
        obs[f'score_{label}'] = values
        obs[f'raw_marker_{label}'] = values
    adata = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=genes))
    markers = {
        'CKB/FABP1+ metabolic enterocytes': ['CKB', 'FABP1'],
        # 面板内 CYP3A5 是唯一锚点；MAML3 属于另一个 programme，因此
        # 单个锚点重叠走 anchor-supported，而不是被误算成两个 signature。
        'CYP3A5+ enterocytes': ['CYP3A5'],
        'HSD17B2+ absorptive enterocytes': ['HSD17B2'],
    }
    decisions = build_marker_decisions(
        adata,
        'leiden',
        markers,
        marker_set_name='Colorectal_refined',
        multi_evidence=True,
        agreement_threshold=0.6,
        cluster_marker_selection={
            'cluster_markers': {
                '0': ['CKB', 'FABP1'],
                '1': ['CYP3A5', 'MAML3'],
                '2': ['HSD17B2', 'LAMA1'],
            },
        },
    )

    assert decisions['cluster_annotations'] == {
        '0': 'CKB/FABP1+ metabolic enterocytes',
        '1': 'CYP3A5+ enterocytes',
        '2': 'HSD17B2+ absorptive enterocytes',
    }
    evidence = decisions['cluster_evidence']['1']
    assert evidence['refined_signature_mode_applied'] is True
    # CYP3A5 是人工锚点：单个强特异 marker 足以覆盖共享高分模块，但等级为
    # anchor-supported，而不是把未验证的 MAML3 当作第二签名。
    assert evidence['refined_signature_overlap'] == ['CYP3A5']
    assert evidence['refined_signature_anchor_overlap'] == ['CYP3A5']
    assert evidence['evidence_tier'] == 'anchor-supported'
    assert decisions['cluster_evidence']['0']['evidence_tier'] == 'confirmed'
    assert decisions['cluster_evidence']['2']['evidence_tier'] == 'anchor-supported'


def test_annotation_exports_condition_stacked_composition_without_donor(tmp_path):
    """The Control/IBD composition plot needs only a condition column."""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK', 'NKG7', 'KLRD1'] + [f'Gene{i}' for i in range(36)]
    rng = np.random.default_rng(71)
    x = rng.poisson(1.0, size=(12, len(genes))).astype(float)
    x[:6, :4] += 5
    x[6:, 4:6] += 5
    adata = ad.AnnData(
        x,
        obs=pd.DataFrame({
            'leiden': pd.Categorical(['0'] * 6 + ['1'] * 6),
            'condition': ['IBD'] * 6 + ['Control'] * 6,
        }),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = rng.normal(size=(12, 2))
    input_path = tmp_path / 'condition_only.h5ad'
    adata.write_h5ad(input_path)

    result = AnnotationAnalysis(
        str(tmp_path),
        {
            'method': 'auto_marker', 'marker_set': 'Universal', 'cluster_key': 'leiden',
            'condition_key': 'condition', 'min_markers_per_type': 2,
            'show_celltype_composition': True, 'show_marker_score_heatmap': False,
            'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
        },
        lambda *_: None,
    ).run(str(input_path))

    assert any(
        os.path.basename(item['file_path']) == 'annotation_condition_celltype_composition.png'
        for item in result['result_files']
    )


def test_annotation_requires_explicit_map_for_cluster_merge(tmp_path):
    """A legacy threshold is ignored; an explicit map preserves source IDs."""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK'] + [f'G{i}' for i in range(30)]
    x = np.random.default_rng(83).poisson(1.0, size=(12, len(genes))).astype(float)
    x[:, :4] += 5
    adata = ad.AnnData(
        x,
        obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 6 + ['1'] * 6)}),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = np.random.default_rng(84).normal(size=(12, 2))
    input_path = tmp_path / 'unmerged_source_clusters.h5ad'
    adata.write_h5ad(input_path)

    result = AnnotationAnalysis(
        str(tmp_path),
        {
            'method': 'auto_marker', 'marker_set': 'Universal', 'cluster_key': 'leiden',
            'merge_similar_threshold': 0.01,
            'cluster_merge_map': 'Merged T cells=0,1', 'min_markers_per_type': 2,
            'show_celltype_composition': False, 'show_marker_score_heatmap': False,
            'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
        },
        lambda *_: None,
    ).run(str(input_path))

    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert set(output.obs['leiden'].astype(str)) == {'0', '1'}
    assert set(output.obs['annotation_source_cluster'].astype(str)) == {
        'leiden=0', 'leiden=1',
    }
    assert set(output.obs['celltype'].astype(str)) == {'Merged T cells'}
    assert output.uns['annotation_metadata']['cluster_merge_policy'] == 'manual_confirmed_only'
    assert result['summary']['n_source_clusters'] == 2
    assert result['summary']['cluster_merge_policy'] == 'manual_confirmed_only'
    assert result['summary']['confirmed_cluster_merges'] == {'Merged T cells': ['0', '1']}
    assert any('已忽略旧版“相似簇合并阈值”' in warning
               for warning in result['summary']['runtime_warnings'])


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


def test_cluster_review_flags_zero_evidence_cell_vote_recovery():
    """仅靠逐细胞投票恢复、无正向 marker 证据的簇必须进入人工复核。"""
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_annotation_cluster_review

    adata = ad.AnnData(
        np.ones((6, 2)),
        obs=pd.DataFrame({
            'leiden': ['0', '0', '0', '1', '1', '1'],
            # cluster 0：投票一致且 margin 充足，但完全没有 marker 证据
            'celltype': ['T cells'] * 3 + ['B cells'] * 3,
            'final_annotation': ['T cells'] * 3 + ['B cells'] * 3,
            'annotation_status': ['research'] * 6,
            'annotation_evidence_tier': ['provisional'] * 6,
            'annotation_score_margin': [0.4] * 6,
            'annotation_decision_reason': (
                ['unknown_recovered_by_cell_vote'] * 3
                + ['cluster_mean_marker_decision'] * 3
            ),
        }),
        var=pd.DataFrame(index=['g1', 'g2']),
    )
    review = build_annotation_cluster_review(
        adata, 'leiden', {'cluster_markers': {'0': ['TRAC'], '1': ['MS4A1']}},
    ).set_index('cluster')

    assert review.loc['0', 'final_decision_reason'] == 'unknown_recovered_by_cell_vote'
    assert review.loc['0', 'needs_review']
    assert not review.loc['1', 'needs_review']


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
        transformed = query.X.toarray() if hasattr(query.X, 'toarray') else np.asarray(query.X)
        nonempty = counts.sum(axis=1) > 0
        np.testing.assert_allclose(
            np.expm1(transformed[nonempty]).sum(axis=1), 10_000.0, rtol=1e-6,
        )
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


def test_celltypist_model_must_be_regular_file_inside_managed_dir(tmp_path):
    from modules.annotation import resolve_celltypist_model

    managed = tmp_path / 'managed_models'
    managed.mkdir()
    valid = managed / 'custom.pkl'
    valid.write_bytes(b'test')
    outside = tmp_path / 'outside.pkl'
    outside.write_bytes(b'outside')

    assert resolve_celltypist_model('custom.pkl', str(managed)) == str(valid.resolve())
    with pytest.raises(ValueError, match='文件名'):
        resolve_celltypist_model(str(outside), str(managed))
    with pytest.raises(ValueError, match='文件名'):
        resolve_celltypist_model('../outside.pkl', str(managed))

    link = managed / 'linked.pkl'
    try:
        link.symlink_to(outside)
    except OSError:
        pytest.skip('当前文件系统不支持符号链接测试')
    with pytest.raises(ValueError, match='符号链接'):
        resolve_celltypist_model('linked.pkl', str(managed))


def test_annotation_uses_counts_for_signed_x_and_preserves_input_representation(tmp_path):
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK', 'NKG7', 'KLRD1'] + [f'G{i}' for i in range(40)]
    counts = np.random.default_rng(301).poisson(1, (12, len(genes))).astype(float)
    counts[:6, :4] += 6
    residuals = np.linspace(-2.5, 2.5, counts.size).reshape(counts.shape)
    adata = ad.AnnData(
        residuals.copy(),
        obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 6 + ['1'] * 6)}),
        var=pd.DataFrame(index=genes),
    )
    adata.layers['counts'] = counts
    adata.uns['normalization'] = {'x_contains': 'pearson_residuals'}
    adata.obsm['X_umap'] = np.random.default_rng(302).normal(size=(12, 2))
    input_path = tmp_path / 'signed_annotation.h5ad'
    adata.write_h5ad(input_path)

    result = AnnotationAnalysis(str(tmp_path), {
        'method': 'auto_marker', 'marker_set': 'Universal', 'cluster_key': 'leiden',
        'min_markers_per_type': 2, 'show_celltype_composition': False,
        'show_marker_score_heatmap': False, 'show_marker_expression_violin': False,
        'show_annotation_score_umap': False,
    }, lambda *_: None).run(str(input_path))

    output = __import__('scanpy').read_h5ad(result['output_adata'])
    np.testing.assert_allclose(np.asarray(output.X), residuals)
    assert 'layers["counts"]' in result['summary']['annotation_expression_source']
    assert output.uns['annotation_metadata']['annotation_expression_source'] == (
        result['summary']['annotation_expression_source']
    )


def test_llm_assisted_annotation_keeps_marker_and_llm_audit_fields(tmp_path, monkeypatch):
    """LLM mode labels clusters but never gives the model raw AnnData data."""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules import annotation as annotation_module

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK', 'MS4A1', 'CD79A'] + [f'G{i}' for i in range(30)]
    counts = np.random.default_rng(202).poisson(1, (10, len(genes))).astype(float)
    counts[:5, :4] += 5
    counts[5:, 4:6] += 5
    adata = ad.AnnData(
        counts,
        obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 5 + ['1'] * 5)}),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = np.random.default_rng(203).normal(size=(10, 2))
    input_path = tmp_path / 'llm_input.h5ad'
    adata.write_h5ad(input_path)
    captured = {}

    def fake_llm_annotation(cluster_summaries, **kwargs):
        captured['summaries'] = cluster_summaries
        captured['kwargs'] = kwargs
        assert all('input_path' not in row and 'donor_id' not in row for row in cluster_summaries)
        return {
            'enabled': True,
            'provider': 'openai',
            'model': 'fake-annotation-model',
            'prompt_version': 'test-v1',
            'request_hash': 'a' * 64,
            'cluster_input': cluster_summaries,
            'tissue_context_provided': True,
            'sent_data_policy': 'cluster-only',
            'global_note': 'Review all clusters.',
            'annotations': [
                {'cluster': '0', 'cell_type': 'T cells', 'confidence': 'high',
                 'rationale': 'CD3D and TRAC.', 'review_note': ''},
                {'cluster': '1', 'cell_type': 'B cells', 'confidence': 'medium',
                 'rationale': 'MS4A1 and CD79A.', 'review_note': ''},
            ],
        }

    monkeypatch.setattr(annotation_module, 'run_llm_cluster_annotation', fake_llm_annotation)
    result = annotation_module.AnnotationAnalysis(
        str(tmp_path),
        {
            'method': 'llm_assisted', 'marker_set': 'Universal', 'cluster_key': 'leiden',
            'min_markers_per_type': 2, 'llm_species': 'human',
            'llm_tissue_context': 'human PBMC', 'llm_max_clusters': 2,
            'show_celltype_composition': False, 'show_marker_score_heatmap': False,
            'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
        },
        lambda *_: None,
    ).run(str(input_path))

    assert captured['kwargs']['species'] == 'human'
    assert result['summary']['method_used'] == 'llm_assisted'
    assert result['summary']['llm_annotation']['model'] == 'fake-annotation-model'
    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert {'marker_label', 'llm_annotation', 'llm_confidence', 'llm_rationale'} <= set(output.obs.columns)
    assert set(output.obs['celltype'].astype(str)) == {'T cells', 'B cells'}
    assert set(output.obs['annotation_source'].astype(str)) == {'llm_cluster_evidence'}


def test_organoid_annotation_uses_selected_tissue_panel(tmp_path):
    """Organoid + organoid_type must reach scoring and coverage evidence."""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['NPHS1', 'NPHS2', 'PODXL', 'WT1', 'SIX2', 'CITED1', 'MKI67'] + [f'Gene{i}' for i in range(30)]
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
            'cell_state', 'dominant_cell_state', 'cell_state_flags',
            'state_high_cycling', 'developmental_state'} <= set(output.obs.columns)
    assert 'library-size normalized log1p' in result['summary']['cell_state_evidence']['expression_source']
    assert 'dominant_state_counts' in result['summary']['cell_state_evidence']
    assert 'multi_label_state_counts' in result['summary']['cell_state_evidence']


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
    assert evidence['doublet_status'][0] == 'not_evaluated'
    assert evidence['lineage_mixture_status'][0] == 'mixed_lineage_review'
    assert evidence['top1'][0] == 'Epithelial'
    assert evidence['top2'][0] == 'Neuron'
    assert 'EPCAM' in evidence['ambient_genes']

    single = compute_annotation_quality_evidence(
        adata, ['score_Epithelial'], {'Epithelial': ['EPCAM']}
    )
    assert set(single['doublet_status']) == {'not_evaluated'}


def test_colorectal_hierarchy_recovers_goblet_and_supported_myeloid_cluster():
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_marker_decisions

    genes = ['EPCAM', 'KRT8', 'FCGBP', 'CLCA1', 'MUC2', 'ITLN1',
             'LST1', 'CSF1R', 'STAB1', 'MS4A6A']
    x = np.zeros((10, len(genes)), dtype=float)
    x[:5, :6] = 5
    x[5:, 6:] = 6
    obs = pd.DataFrame({'leiden': ['12'] * 5 + ['15'] * 5})
    # Cluster 15 has only 3/5 per-cell Myeloid top votes.  Strong, specific
    # LST1/CSF1R/STAB1/MS4A6A support must keep it out of Unknown.
    obs['score_Epithelial'] = [2.0] * 5 + [0.1, 0.1, 0.1, 1.8, 1.8]
    obs['score_Myeloid'] = [0.1] * 5 + [2.0, 2.0, 2.0, 0.2, 0.2]
    obs['score_Goblet-like epithelial'] = [3.0] * 5 + [0.0] * 5
    obs['score_Absorptive/enterocyte-like epithelial'] = [0.2] * 10
    for label in ('Epithelial', 'Myeloid', 'Goblet-like epithelial',
                  'Absorptive/enterocyte-like epithelial'):
        obs[f'raw_marker_{label}'] = obs[f'score_{label}']
    adata = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=genes))
    markers = {
        'Epithelial': ['EPCAM', 'KRT8'],
        'Myeloid': ['LST1', 'CSF1R', 'STAB1', 'MS4A6A'],
        'Goblet-like epithelial': ['FCGBP', 'CLCA1', 'MUC2', 'ITLN1'],
        'Absorptive/enterocyte-like epithelial': ['EPCAM', 'KRT8'],
    }

    decisions = build_marker_decisions(
        adata, 'leiden', markers, marker_set_name='Colorectal',
        multi_evidence=True, agreement_threshold=0.8,
        marker_min_pct=0.1, marker_min_delta_pct=0.05,
    )

    assert decisions['cluster_annotations']['12'] == 'Goblet-like epithelial'
    assert decisions['cluster_annotations']['15'] == 'Myeloid'
    myeloid = decisions['cluster_evidence']['15']
    assert myeloid['cell_prediction_agreement'] == pytest.approx(0.6)
    assert myeloid['strong_marker_support'] is True
    assert myeloid['decision_reason'] == 'low_cell_agreement_but_cluster_markers_support_candidate'
    assert {'LST1', 'CSF1R'} <= set(myeloid['detected_positive_markers'])


def test_doublet_status_is_inherited_from_qc_not_marker_mixture():
    ad = pytest.importorskip('anndata')
    from modules.annotation import compute_annotation_quality_evidence

    adata = ad.AnnData(
        np.ones((4, 2)),
        obs=pd.DataFrame({
            'predicted_doublet': [False, True, False, False],
            'doublet_score': [0.1, 0.7, 0.2, 0.1],
            'score_Epithelial': [1.0, 1.0, 1.0, 1.0],
            'score_Myeloid': [0.95, 0.95, 0.95, 0.95],
        }),
        var=pd.DataFrame(index=['EPCAM', 'LST1']),
    )
    evidence = compute_annotation_quality_evidence(
        adata, ['score_Epithelial', 'score_Myeloid'],
        {'Epithelial': ['EPCAM'], 'Myeloid': ['LST1']},
        assigned_labels=['Epithelial'] * 4,
    )

    assert list(evidence['doublet_status']).count('qc_predicted_doublet') == 1
    assert list(evidence['lineage_mixture_status']).count('mixed_lineage_review') == 4
    assert evidence['doublet_source'] == 'qc_predicted_doublet'


def test_ambient_signal_requires_foreign_lineage_markers():
    ad = pytest.importorskip('anndata')
    from modules.annotation import compute_annotation_quality_evidence

    counts = np.asarray([
        [5, 5, 0, 0], [4, 4, 0, 0],  # expected epithelial expression
        [3, 3, 5, 5], [3, 3, 4, 4],  # epithelial signal inside T cells
    ], dtype=float)
    adata = ad.AnnData(
        counts.copy(),
        obs=pd.DataFrame({
            'score_Epithelial': [2, 2, 0, 0],
            'score_T cells': [0, 0, 2, 2],
        }),
        var=pd.DataFrame(index=['EPCAM', 'KRT8', 'CD3D', 'TRAC']),
    )
    adata.layers['counts'] = counts
    evidence = compute_annotation_quality_evidence(
        adata, ['score_Epithelial', 'score_T cells'],
        {'Epithelial': ['EPCAM', 'KRT8'], 'T cells': ['CD3D', 'TRAC']},
        assigned_labels=['Epithelial', 'Epithelial', 'T cells', 'T cells'],
        ambient_prevalence=0.25, ambient_threshold=0.5,
    )

    assert set(evidence['ambient_status'][:2]) == {'no_strong_ambient_signal'}
    assert set(evidence['ambient_status'][2:]) == {'ambient_signal'}


def test_cell_state_scores_are_log_normalized_and_multi_label():
    ad = pytest.importorskip('anndata')
    from modules.annotation import compute_cell_state_evidence

    genes = ['ISG15', 'IFIT1', 'FOS', 'JUN', 'OTHER']
    counts = np.asarray([
        [10, 10, 10, 10, 60],
        [100, 100, 100, 100, 600],
    ], dtype=float)
    adata = ad.AnnData(counts.copy(), var=pd.DataFrame(index=genes))
    adata.layers['counts'] = counts
    evidence = compute_cell_state_evidence(adata, threshold=0.5)

    assert 'library-size normalized log1p' in evidence['expression_source']
    assert evidence['scores']['Interferon response'][0] == pytest.approx(
        evidence['scores']['Interferon response'][1]
    )
    assert evidence['high_flags']['Interferon response'].all()
    assert evidence['high_flags']['Stress response'].all()
    assert all('Interferon response' in flags and 'Stress response' in flags
               for flags in evidence['state_flags'])


def test_fine_annotation_splits_t_cell_clusters_and_falls_back_to_broad(tmp_path):
    """细分注释：T 谱系内输出 CD4+/CD8+ 亚型，无亚型证据时回退大谱系。"""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK', 'CD4', 'IL7R', 'LEF1', 'CCR7', 'TCF7',
             'CD8A', 'CD8B', 'GZMK', 'GZMA', 'CCL5'] + [f'G{i}' for i in range(60)]
    rng = np.random.default_rng(41)
    x = rng.poisson(1.0, size=(30, len(genes))).astype(float)
    x[:10, 0:4] += 6
    x[10:20, 0:4] += 6
    x[10:20, 4:9] += 6
    x[20:30, 0:4] += 6
    x[20:30, 9:14] += 6
    adata = ad.AnnData(
        x,
        obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 10 + ['1'] * 10 + ['2'] * 10)}),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = rng.normal(size=(30, 2))
    input_path = tmp_path / 'fine.h5ad'
    adata.write_h5ad(input_path)

    result = AnnotationAnalysis(str(tmp_path), {
        'method': 'multi_evidence', 'marker_set': 'Universal', 'cluster_key': 'leiden',
        'min_markers_per_type': 2, 'fine_annotation': True, 'annotate_all': True,
        'show_celltype_composition': False, 'show_marker_score_heatmap': False,
        'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
    }, lambda *_: None).run(str(input_path))

    output = __import__('scanpy').read_h5ad(result['output_adata'])
    by_cluster = output.obs.groupby('leiden', observed=True)['celltype'].apply(
        lambda s: s.astype(str).mode().iloc[0]
    )
    assert by_cluster['0'] == 'T cells'          # 无亚型 marker -> 大谱系回退
    assert by_cluster['1'] == 'CD4+ T cells'
    assert by_cluster['2'] == 'CD8+ T cells'
    assert (output.obs['celltype'].astype(str) != 'Unknown').all()
    assert bool(output.uns['annotation_metadata']['fine_annotation']) is True
    assert bool(output.uns['annotation_metadata']['fine_mode_applied']) is True
    assert bool(output.uns['annotation_metadata']['annotate_all']) is True
    assert output.obs['cell_type_l3'].isin(['CD4+ T cell', 'CD8+ T cell']).any()
    assert result['summary']['agreement_source'] == 'hierarchical_prediction_1'


def test_annotate_all_keeps_low_agreement_label_for_review():
    """尽量注释所有细胞：低一致率且无 Marker 支持的簇保留标签并标记复核。"""
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_marker_decisions

    genes = ['CD3D', 'TRAC', 'MS4A1', 'CD79A']
    x = np.zeros((5, len(genes)), dtype=float)
    x[0, 0] = 5
    x[1, 0] = 5
    x[2, 2] = 5
    x[3, 2] = 5
    x[4, 2] = 5
    obs = pd.DataFrame({
        'leiden': ['0'] * 5,
        'score_T cells': [3.0, 3.0, 0.5, 0.5, 0.5],
        'score_B cells': [0.3, 0.3, 2.0, 2.0, 2.0],
        'raw_marker_T cells': [3.0, 3.0, 0.5, 0.5, 0.5],
        'raw_marker_B cells': [0.3, 0.3, 2.0, 2.0, 2.0],
    })
    adata = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=genes))
    markers = {'T cells': ['CD3D', 'TRAC'], 'B cells': ['MS4A1', 'CD79A']}

    decisions = build_marker_decisions(
        adata, 'leiden', markers, multi_evidence=True, agreement_threshold=0.6,
    )
    assert decisions['cluster_annotations']['0'] == 'T cells'
    assert decisions['cluster_evidence']['0']['decision_reason'] == (
        'low_cell_agreement_label_kept_for_review'
    )
    assert decisions['cluster_evidence']['0']['strong_marker_support'] is False

    strict = build_marker_decisions(
        adata, 'leiden', markers, multi_evidence=True,
        agreement_threshold=0.6, annotate_all=False,
    )
    assert strict['cluster_annotations']['0'] == 'Unknown'
    assert strict['cluster_evidence']['0']['decision_reason'] == (
        'low_cell_agreement_without_specific_marker_support'
    )


def test_annotate_all_keeps_unknown_when_no_positive_evidence_at_all():
    """完全没有正向 Marker 证据时仍标 Unknown，而不是硬猜。"""
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_marker_decisions

    adata = ad.AnnData(
        np.zeros((4, 2), dtype=float),
        obs=pd.DataFrame({
            'leiden': ['0'] * 4,
            'score_T cells': [-0.2, -0.1, -0.3, -0.2],
            'score_B cells': [-0.1, -0.2, -0.1, -0.3],
            'raw_marker_T cells': [-0.2, -0.1, -0.3, -0.2],
            'raw_marker_B cells': [-0.1, -0.2, -0.1, -0.3],
        }),
        var=pd.DataFrame(index=['CD3D', 'MS4A1']),
    )
    decisions = build_marker_decisions(
        adata, 'leiden', {'T cells': ['CD3D'], 'B cells': ['MS4A1']},
        multi_evidence=True,
    )
    assert decisions['cluster_annotations']['0'] == 'Unknown'
    assert decisions['cluster_evidence']['0']['decision_reason'] == 'no_positive_marker_evidence'


def test_annotation_accepts_high_resolution_many_clusters(tmp_path):
    """高分辨率聚类（60 个簇）也能进入初步注释，而不是被分组上限拦住。"""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK'] + [f'G{i}' for i in range(26)]
    n_clusters = 60
    cells_per_cluster = 5
    n_cells = n_clusters * cells_per_cluster
    rng = np.random.default_rng(5)
    x = rng.poisson(1.0, size=(n_cells, len(genes))).astype(float)
    x[:, :4] += 3
    adata = ad.AnnData(
        x,
        obs=pd.DataFrame({'leiden': pd.Categorical(
            np.repeat(np.arange(n_clusters).astype(str), cells_per_cluster),
        )}),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = rng.normal(size=(n_cells, 2))
    input_path = tmp_path / 'many_clusters.h5ad'
    adata.write_h5ad(input_path)

    result = AnnotationAnalysis(str(tmp_path), {
        'method': 'auto_marker', 'marker_set': 'Universal', 'cluster_key': 'leiden',
        'min_markers_per_type': 2,
        'show_celltype_composition': False, 'show_marker_score_heatmap': False,
        'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
    }, lambda *_: None).run(str(input_path))

    output = __import__('scanpy').read_h5ad(result['output_adata'])
    assert output.obs['leiden'].nunique() == n_clusters
    assert result['summary']['n_celltypes'] >= 1
    assert output.obs['celltype'].notna().all()


def test_annotation_honors_explicit_cluster_key_over_resolution_field(tmp_path):
    """A selected cluster column must not be overwritten by the UI default resolution."""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules.annotation import AnnotationAnalysis

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK'] + [f'G{i}' for i in range(24)]
    x = np.random.default_rng(25).poisson(1, (12, len(genes))).astype(float)
    x[:, :4] += 4
    adata = ad.AnnData(
        x,
        obs=pd.DataFrame({
            'leiden_0.6': pd.Categorical(['0'] * 6 + ['1'] * 6),
            'leiden_0.8': pd.Categorical(['0', '1', '2'] * 4),
            'leiden': pd.Categorical(['0', '1', '2'] * 4),
        }),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = np.random.default_rng(26).normal(size=(12, 2))
    input_path = tmp_path / 'explicit_cluster_key.h5ad'
    adata.write_h5ad(input_path)

    result = AnnotationAnalysis(str(tmp_path), {
        'method': 'auto_marker', 'marker_set': 'Universal',
        'cluster_key': 'leiden_0.6', 'resolution': '0.8',
        'min_markers_per_type': 2,
        'show_celltype_composition': False, 'show_marker_score_heatmap': False,
        'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
    }, lambda *_: None).run(str(input_path))

    assert result['summary']['cluster_column'] == 'leiden_0.6'


def test_colorectal_second_pass_does_not_let_broad_epithelial_mask_subtype():
    """Generic EPCAM/KRT evidence cannot out-rank positive subtype evidence."""
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_marker_decisions

    epithelial_subtypes = [
        'Stem/crypt-like epithelial', 'TA/S-phase epithelial',
        'TA/G2M epithelial', 'Goblet-like epithelial',
        'Absorptive/enterocyte-like epithelial', 'BEST4+ absorptive epithelial',
        'Enteroendocrine-like epithelial', 'Paneth/LYZ+ secretory epithelial',
        'Inflammatory epithelial', 'Regenerative/stress epithelial',
    ]
    obs = pd.DataFrame({'leiden': ['0'] * 4})
    obs['score_Epithelial'] = [3.0] * 4
    obs['raw_marker_Epithelial'] = [3.0] * 4
    for label in epithelial_subtypes:
        score = 1.2 if label == 'Goblet-like epithelial' else -0.1
        obs[f'score_{label}'] = [score] * 4
        obs[f'raw_marker_{label}'] = [score] * 4
    adata = ad.AnnData(
        np.ones((4, 2)), obs=obs, var=pd.DataFrame(index=['MUC2', 'EPCAM']),
    )
    decisions = build_marker_decisions(
        adata, 'leiden',
        {
            'Epithelial': ['EPCAM'],
            'Goblet-like epithelial': ['MUC2'],
            **{label: ['EPCAM'] for label in epithelial_subtypes if label != 'Goblet-like epithelial'},
        },
        marker_set_name='Colorectal',
    )
    assert decisions['cluster_annotations']['0'] == 'Goblet-like epithelial'


def test_unknown_recovered_by_cell_vote_when_cells_agree():
    """零正向证据但逐细胞投票一致的簇也给出标签（review），而不是 Unknown。"""
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_marker_decisions

    adata = ad.AnnData(
        np.zeros((4, 2), dtype=float),
        obs=pd.DataFrame({
            'leiden': ['0'] * 4,
            'score_T cells': [-0.1, -0.2, -0.1, -0.3],
            'score_B cells': [-0.5, -0.6, -0.4, -0.7],
            'raw_marker_T cells': [-0.1, -0.2, -0.1, -0.3],
            'raw_marker_B cells': [-0.5, -0.6, -0.4, -0.7],
        }),
        var=pd.DataFrame(index=['CD3D', 'MS4A1']),
    )
    decisions = build_marker_decisions(
        adata, 'leiden', {'T cells': ['CD3D'], 'B cells': ['MS4A1']},
        multi_evidence=True, agreement_threshold=0.6,
    )
    evidence = decisions['cluster_evidence']['0']
    assert decisions['cluster_annotations']['0'] == 'T cells'
    assert evidence['decision_reason'] == 'unknown_recovered_by_cell_vote'
    assert evidence['cell_prediction_agreement'] == 1.0
    assert evidence['score_margin'] > 0


def test_llm_unknown_falls_back_to_marker_candidate(tmp_path, monkeypatch):
    """LLM 判 Unknown 的簇回退到 Marker 规则候选标签，llm_annotation 保留原判。"""
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    from modules import annotation as annotation_module

    genes = ['CD3D', 'TRAC', 'CD3E', 'LCK', 'MS4A1', 'CD79A'] + [f'G{i}' for i in range(30)]
    counts = np.random.default_rng(202).poisson(1, (10, len(genes))).astype(float)
    counts[:5, :4] += 5
    counts[5:, 4:6] += 5
    adata = ad.AnnData(
        counts,
        obs=pd.DataFrame({'leiden': pd.Categorical(['0'] * 5 + ['1'] * 5)}),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm['X_umap'] = np.random.default_rng(203).normal(size=(10, 2))
    input_path = tmp_path / 'llm_fallback.h5ad'
    adata.write_h5ad(input_path)

    def fake_llm_annotation(cluster_summaries, **kwargs):
        return {
            'enabled': True,
            'provider': 'openai',
            'model': 'fake-annotation-model',
            'prompt_version': 'test-v1',
            'request_hash': 'b' * 64,
            'cluster_input': cluster_summaries,
            'tissue_context_provided': True,
            'sent_data_policy': 'cluster-only',
            'global_note': 'Review all clusters.',
            'annotations': [
                {'cluster': '0', 'cell_type': 'T cells', 'confidence': 'high',
                 'rationale': 'CD3D and TRAC.', 'review_note': ''},
                {'cluster': '1', 'cell_type': 'Unknown', 'confidence': 'low',
                 'rationale': 'mixed markers.', 'review_note': 'Check B cell markers.'},
            ],
        }

    monkeypatch.setattr(annotation_module, 'run_llm_cluster_annotation', fake_llm_annotation)
    result = annotation_module.AnnotationAnalysis(
        str(tmp_path),
        {
            'method': 'llm_assisted', 'marker_set': 'Universal', 'cluster_key': 'leiden',
            'min_markers_per_type': 2, 'llm_species': 'human',
            'llm_tissue_context': 'human PBMC', 'llm_max_clusters': 2,
            'show_celltype_composition': False, 'show_marker_score_heatmap': False,
            'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
        },
        lambda *_: None,
    ).run(str(input_path))

    output = __import__('scanpy').read_h5ad(result['output_adata'])
    cluster1 = output.obs[output.obs['leiden'].astype(str) == '1']
    assert set(cluster1['celltype'].astype(str)) == {'B cells'}   # Marker 回退
    assert set(cluster1['llm_annotation'].astype(str)) == {'Unknown'}  # 原判保留
    assert cluster1['llm_unknown_fallback'].all()
    assert result['summary']['n_llm_unknown_fallback'] == 5
    assert output.uns['annotation_metadata']['unknown_recovery'][
        'llm_unknown_marker_fallback_cells'] == 5
