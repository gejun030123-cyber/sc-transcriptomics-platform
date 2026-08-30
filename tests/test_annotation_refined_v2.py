"""Cluster-marker-first fine annotation (colorectal_refined_v2) tests.

覆盖规格中的验收点：
- Auto 三路路由（Universal / Colorectal / Colorectal refined v2）
- 17 个细分 programme 面板与面板版本 manifest
- evidence tier（confirmed / anchor-supported / provisional / unresolved-review）
- 来源 cluster 永不自动合并
- 可下载复核 CSV、人工审阅模板、marker panel manifest
- Organoid 面板自动派生锚点后的细标签优先判定
"""
import os

import numpy as np
import pandas as pd
import pytest

sys_path_inserted = False


def _syspath():
    global sys_path_inserted
    if not sys_path_inserted:
        import sys
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        sys_path_inserted = True


def _build_refined_adata(n_cells_per_cluster=4, seed=7):
    _syspath()
    ad = pytest.importorskip('anndata')
    programmes = {
        '0': ['CYP3A5', 'A1'],        # 单锚点 CYP3A5 -> anchor-supported
        '1': ['HSD17B2', 'C11orf86'], # confirmed
        '2': ['MUC2', 'FCGBP'],       # confirmed
        '3': ['ASCL2', 'LGR5'],       # confirmed
        '4': ['DDIT3', 'ATF3'],       # confirmed
        '5': ['FEN1', 'MCM10'],       # confirmed
        '6': ['UNRES1', 'UNRES2'],    # 证据不足 -> unresolved-review
    }
    panel_genes = [
        'CYP3A5', 'PTPRR', 'HSD17B2', 'C11orf86', 'MUC2', 'FCGBP',
        'ASCL2', 'LGR5', 'DDIT3', 'ATF3', 'FEN1', 'MCM10',
        'A1', 'UNRES1', 'UNRES2',
    ]
    rng = np.random.default_rng(seed)
    rows, labels = [], []
    for cluster, genes in programmes.items():
        for _ in range(n_cells_per_cluster):
            rows.append([6.0 if gene in genes else 0.0 for gene in panel_genes])
            labels.append(cluster)
    adata = ad.AnnData(
        np.asarray(rows, dtype=float),
        obs=pd.DataFrame({'leiden': pd.Categorical(labels)}),
        var=pd.DataFrame(index=panel_genes),
    )
    adata.obsm['X_umap'] = rng.normal(size=(len(labels), 2))
    return adata


def _run_refined(tmp_path, adata, extra_params=None, name='refined'):
    from modules.annotation import AnnotationAnalysis
    input_path = tmp_path / f'{name}_input.h5ad'
    adata.write_h5ad(input_path)
    params = {
        'method': 'auto_marker', 'marker_set': 'Colorectal_refined',
        'cluster_key': 'leiden', 'min_markers_per_type': 2,
        'show_celltype_composition': False, 'show_marker_score_heatmap': False,
        'show_marker_expression_violin': False, 'show_annotation_score_umap': False,
    }
    params.update(extra_params or {})
    return AnnotationAnalysis(
        str(tmp_path), params, lambda *_: None,
    ).run(str(input_path))


def test_auto_routing_reaches_refined_panel():
    _syspath()
    ad = pytest.importorskip('anndata')
    from modules.annotation import resolve_auto_marker_set

    genes = [
        'CYP3A5', 'PTPRR', 'HSD17B2', 'C11orf86', 'MUC2', 'FCGBP',
        'ASCL2', 'LGR5', 'ALPI', 'VIL1', 'TYMS', 'PCNA',
    ]
    refined_adata = ad.AnnData(
        np.ones((8, len(genes)), dtype=float) * 5,
        var=pd.DataFrame(index=genes),
    )
    panel, evidence = resolve_auto_marker_set(refined_adata, 'Auto')
    assert panel == 'Colorectal_refined'
    assert evidence['panel_version'] == 'colorectal_refined_v2'
    assert evidence['n_supported_refined_programmes'] >= 3
    assert evidence['n_supported_programmes'] >= 3

    broad_only = ad.AnnData(
        np.ones((6, 6), dtype=float) * 5,
        var=pd.DataFrame(index=['MUC2', 'TFF3', 'MKI67', 'TOP2A', 'ALPI', 'VIL1']),
    )
    panel, evidence = resolve_auto_marker_set(broad_only, 'Auto')
    assert panel == 'Colorectal'
    # MUC2+TFF3 只支撑 1 个细分 programme，不足以进入 refined v2。
    assert evidence['n_supported_refined_programmes'] <= 2

    universal = ad.AnnData(
        np.ones((4, 3), dtype=float),
        var=pd.DataFrame(index=['EPCAM', 'KRT8', 'KRT18']),
    )
    panel, evidence = resolve_auto_marker_set(universal, 'Auto')
    assert panel == 'Universal'
    assert evidence['n_supported_programmes'] == 0


def test_refined_end_to_end_emits_tiers_and_review_files(tmp_path):
    _syspath()
    adata = _build_refined_adata()
    result = _run_refined(tmp_path, adata)

    summary = result['summary']
    assert summary['n_source_clusters'] == 7
    assert summary['n_final_labels'] >= 6
    assert summary['label_reduction_warning'] is None
    assert summary['panel_version'] == 'colorectal_refined_v2'

    basenames = {os.path.basename(item['file_path']) for item in result['result_files']}
    assert 'annotation_cluster_review.csv' in basenames
    assert 'annotation_manual_map_template.csv' in basenames
    assert 'annotation_marker_panel_manifest.json' in basenames

    output = __import__('scanpy').read_h5ad(result['output_adata'])
    required_columns = {
        'annotation_source_cluster', 'celltype_l1', 'celltype_l2', 'celltype_l3',
        'cell_type_l1', 'cell_type_l2', 'cell_type_l3', 'cell_state_flags',
        'annotation_evidence_tier', 'annotation_decision_reason',
    }
    assert required_columns <= set(output.obs.columns)

    by_cluster = output.obs.groupby('leiden', observed=True)['celltype'].agg(
        lambda series: series.astype(str).mode().iloc[0]
    )
    assert by_cluster['0'] == 'CYP3A5+ Enterocytes'
    assert (output.obs.loc[output.obs['leiden'].astype(str) == '0',
                            'annotation_marker_programme'].astype(str)
            == 'CYP3A5+ enterocytes').all()
    assert by_cluster['6'].startswith('Unresolved epithelial programme:')
    assert 'UNRES1/UNRES2' in by_cluster['6']

    tiers = output.obs.groupby('leiden', observed=True)['annotation_evidence_tier'].agg(
        lambda series: series.astype(str).mode().iloc[0]
    )
    assert tiers['0'] == 'anchor-supported'
    assert tiers['1'] == 'confirmed'
    assert tiers['6'] == 'unresolved-review'

    assert output.uns['marker_panel_manifest']['panel_version'] == 'colorectal_refined_v2'
    manifest_path = os.path.join(
        tmp_path, 'results', 'annotation_marker_panel_manifest.json',
    )
    assert os.path.isfile(manifest_path)
    manifest = __import__('json').load(open(manifest_path, encoding='utf-8'))
    assert manifest['panel_version'] == 'colorectal_refined_v2'
    assert manifest['thresholds']['signature_overlap_min'] == 2

    review = pd.read_csv(
        os.path.join(tmp_path, 'results', 'annotation_cluster_review.csv'), dtype=str,
    )
    assert 'evidence_tier' in review.columns
    assert 'umap_adjacent_clusters' in review.columns
    assert 'knn_adjacent_clusters' in review.columns
    row6 = review.set_index('cluster').loc['6']
    assert row6['evidence_tier'] == 'unresolved-review'
    assert row6['needs_review'] == 'True'

    template = pd.read_csv(
        os.path.join(tmp_path, 'results', 'annotation_manual_map_template.csv'), dtype=str,
    )
    assert {'cluster', 'recommended_label', 'manual_label'} <= set(template.columns)
    assert template.set_index('cluster').loc['6', 'manual_label'] == 'KEEP'
    assert template.set_index('cluster').loc['6', 'recommended_label'].startswith(
        'Unresolved epithelial programme:'
    )


def test_manual_map_csv_is_parsed_and_applied(tmp_path):
    _syspath()
    from modules.annotation import parse_annotation_manual_map_csv

    adata = _build_refined_adata()
    first = _run_refined(tmp_path, adata, name='first')
    template_path = os.path.join(tmp_path, 'results', 'annotation_manual_map_template.csv')
    template = pd.read_csv(template_path, dtype=str)
    template.loc[template['cluster'] == '6', 'manual_label'] = 'Human-confirmed TA cells'
    edited_path = tmp_path / 'edited_manual_map.csv'
    template.to_csv(edited_path, index=False)

    parsed = parse_annotation_manual_map_csv(str(edited_path))
    assert parsed['6'] == 'Human-confirmed TA cells'

    with pytest.raises(ValueError, match='cluster 与 manual_label'):
        parse_annotation_manual_map_csv('cluster,n_cells\n0,4\n')

    second = _run_refined(tmp_path, adata, extra_params={
        'manual_map_csv': str(edited_path),
    }, name='second')
    output = __import__('scanpy').read_h5ad(second['output_adata'])
    by_cluster = output.obs.groupby('leiden', observed=True)['celltype'].agg(
        lambda series: series.astype(str).mode().iloc[0]
    )
    assert by_cluster['6'] == 'Human-confirmed TA cells'
    assert output.obs.loc[
        output.obs['leiden'].astype(str) == '6', 'annotation_evidence_tier'
    ].eq('confirmed').all()
    # 来源 cluster 永远保留，即使人工改了显示标签。
    assert set(output.obs['annotation_source_cluster'].astype(str)) == {
        'leiden=0', 'leiden=1', 'leiden=2', 'leiden=3',
        'leiden=4', 'leiden=5', 'leiden=6',
    }
    assert second['summary']['n_source_clusters'] == 7
    assert second['summary']['manual_map_csv_applied']['n_clusters'] == 1


def test_organoid_uses_derived_anchor_signature():
    _syspath()
    ad = pytest.importorskip('anndata')
    from modules.annotation import build_marker_decisions, derive_panel_anchor_markers

    markers = {
        'Podocytes': ['NPHS1', 'NPHS2', 'PODXL', 'WT1'],
        'Proximal tubule': ['SLC34A1', 'CUBN'],
        'Stromal cells': ['COL1A1', 'DCN'],
    }
    anchors = derive_panel_anchor_markers(markers)
    assert 'NPHS1' in anchors['Podocytes']
    assert 'SLC34A1' in anchors['Proximal tubule']

    genes = ['NPHS1', 'UNIQ', 'SLC34A1', 'COL1A1', 'DCN']
    x = np.asarray([
        [8, 7, 0, 0, 0], [7, 8, 0, 0, 0],
        [0, 0, 8, 0, 0], [0, 0, 7, 0, 0],
        [0, 0, 0, 8, 8], [0, 0, 0, 7, 7],
    ], dtype=float)
    obs = pd.DataFrame({'leiden': ['0', '0', '1', '1', '2', '2']})
    # 共享 Stromal 高分不能覆盖有自身 signature 的簇。
    for label, values in {
        'Stromal cells': [2.0] * 6,
        'Podocytes': [1.0, 1.0, 0.2, 0.2, 0.2, 0.2],
        'Proximal tubule': [0.2, 0.2, 1.0, 1.0, 0.2, 0.2],
    }.items():
        obs[f'score_{label}'] = values
        obs[f'raw_marker_{label}'] = values
    adata = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=genes))
    decisions = build_marker_decisions(
        adata, 'leiden', markers, marker_set_name='Organoid',
        cluster_marker_selection={'cluster_markers': {
            '0': ['NPHS1', 'UNIQ'],
            '1': ['SLC34A1'],
            '2': ['COL1A1', 'DCN'],
        }},
    )
    assert decisions['cluster_annotations']['0'] == 'Podocytes'
    assert decisions['cluster_annotations']['1'] == 'Proximal tubule'
    assert decisions['cluster_annotations']['2'] == 'Stromal cells'
    assert decisions['cluster_evidence']['0']['evidence_tier'] == 'anchor-supported'
    assert decisions['cluster_evidence']['1']['evidence_tier'] == 'anchor-supported'
    assert decisions['cluster_evidence']['2']['evidence_tier'] == 'confirmed'


def test_neighborhood_evidence_is_review_only():
    _syspath()
    ad = pytest.importorskip('anndata')
    from modules.annotation import compute_cluster_neighborhood_evidence

    rng = np.random.default_rng(11)
    coords = np.asarray([
        [0, 0], [0.1, 0.1], [9, 9], [9.1, 9.1], [-9, -9], [-9.1, -9.1],
    ], dtype=float)
    adata = ad.AnnData(
        np.ones((6, 2)),
        obs=pd.DataFrame({'leiden': pd.Categorical(['0', '0', '1', '1', '2', '2'])}),
        var=pd.DataFrame(index=['g1', 'g2']),
    )
    adata.obsm['X_umap'] = coords
    evidence = compute_cluster_neighborhood_evidence(adata, 'leiden', k=2)
    assert evidence['0']['umap_adjacent_clusters'] == ['1', '2']
    assert evidence['1']['umap_adjacent_clusters'] == ['0', '2']
    assert evidence['2']['umap_adjacent_clusters'] == ['0', '1']
    assert len(evidence['0']['knn_adjacent_clusters']) == 0
