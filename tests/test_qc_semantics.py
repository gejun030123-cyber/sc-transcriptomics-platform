import anndata as ad
import numpy as np
import pandas as pd


def test_ribosomal_gene_mask_excludes_s6_kinases_and_pseudogenes():
    from modules.qc import _ribosomal_gene_mask

    genes = [
        'RPS18', 'RPSA', 'RPS4Y1', 'RPS27L', 'RPL13A', 'RPLP0', 'RPL36AL',
        'RPS6KA1', 'RPS6KB1', 'RPS18P3', 'RPL13AP5', 'GAPDH',
    ]

    assert _ribosomal_gene_mask(genes) == [
        True, True, True, True, True, True, True,
        False, False, False, False, False,
    ]


def test_gene_name_mapping_keeps_first_duplicate_deterministically():
    from modules.qc import _first_gene_name_to_var_name

    mapping = _first_gene_name_to_var_name(
        ['ENSG_FIRST', 'ENSG_SECOND', 'ENSG_OTHER'],
        ['MCM5', 'MCM5', 'PCNA'],
    )

    assert mapping == {'MCM5': 'ENSG_FIRST', 'PCNA': 'ENSG_OTHER'}


def test_qc_overview_grid_has_room_for_all_six_metrics():
    from modules.qc import _qc_overview_grid

    assert _qc_overview_grid(6) == (3, 2)
    assert np.prod(_qc_overview_grid(6)) == 6


def test_scdblfinder_dependency_requires_supported_installed_version():
    from importlib.metadata import PackageNotFoundError

    from modules.qc import _scdblfinder_dependency_status

    def missing(_package):
        raise PackageNotFoundError

    assert _scdblfinder_dependency_status(missing) == {
        'available': False,
        'package': 'pyscdblfinder',
        'minimum_version': '0.2.0',
        'installed_version': None,
        'reason': 'pyscdblfinder_missing_or_version_below_0.2.0',
    }
    old = _scdblfinder_dependency_status(lambda _package: '0.1.0')
    assert old['available'] is False
    assert old['installed_version'] == '0.1.0'
    assert old['reason'] == 'pyscdblfinder_missing_or_version_below_0.2.0'
    assert _scdblfinder_dependency_status(lambda _package: '0.2.1')['available'] is True


def test_qc_stops_before_analysis_when_scdblfinder_dependency_is_unavailable(tmp_path, monkeypatch):
    import pytest

    from modules.qc import QCAnalysis

    input_path = tmp_path / 'counts.h5ad'
    ad.AnnData(X=np.ones((2, 2))).write_h5ad(input_path)
    monkeypatch.setattr(
        'modules.qc._scdblfinder_dependency_status',
        lambda: {
            'available': False,
            'package': 'pyscdblfinder',
            'minimum_version': '0.2.0',
            'installed_version': None,
            'reason': 'pyscdblfinder_missing_or_version_below_0.2.0',
        },
    )

    with pytest.raises(ValueError, match='无法运行 scDblFinder'):
        QCAnalysis(
            project_dir=str(tmp_path),
            params={'doublets_method': 'scdblfinder'},
            progress_callback=lambda *_: None,
        ).run(str(input_path))


def test_doublet_method_provenance_records_the_effective_omicverse_method():
    from modules.qc import (
        _doublet_method_provenance,
        _effective_omicverse_doublets_method,
    )

    adata = ad.AnnData(X=np.ones((1, 1)))
    adata.uns['status_args'] = {'qc': {'doublets_method': 'scrublet'}}
    effective = _effective_omicverse_doublets_method(adata)
    assert effective == 'scrublet'
    provenance = _doublet_method_provenance(
        'scdblfinder', effective,
        verification_source='omicverse.status_args.qc.doublets_method',
    )
    assert provenance['requested_doublets_method'] == 'scdblfinder'
    assert provenance['effective_doublets_method'] == 'scrublet'
    assert provenance['fallback'] is True


def test_qc_removal_counts_include_every_filter_stage():
    from modules.qc import _qc_removal_counts

    counts = _qc_removal_counts(
        100, 90, 80,
        {'max_genes': 2, 'ribo': 3, 'hb': 1, 'batch_adaptive_qc': 4},
    )

    assert counts['cells_removed_by_qc_and_doublet'] == 20
    assert counts['cells_removed_by_initial_qc_and_doublet'] == 10
    assert counts['cells_removed_by_additional_qc'] == 10
    assert sum(counts[key] for key in [
        'cells_removed_by_max_genes', 'cells_removed_by_ribo',
        'cells_removed_by_hb', 'cells_removed_by_batch_adaptive_qc',
    ]) == 10


def test_qc_removal_reason_summary_keeps_unique_and_overlapping_counts():
    from modules.qc import _qc_removal_reason_summary

    qc_before = pd.DataFrame({
        'total_counts': [100, 400, 1000, 1000],
        'n_genes_by_counts': [100, 200, 900, 1000],
        'pct_counts_mt': [30, 1, 1, 1],
        'pct_counts_ribo': [1, 1, 40, 1],
        'pct_counts_hb': [1, 1, 1, 20],
    }, index=['a', 'b', 'c', 'd'])
    summary = _qc_removal_reason_summary(
        qc_before,
        ['c', 'd'],
        ad.AnnData(
            X=np.ones((1, 1)),
            obs=pd.DataFrame(index=['d']),
            var=pd.DataFrame(index=['gene']),
        ),
        mito_perc=0.2,
        n_umis=500,
        n_genes_min=250,
        n_genes_max=800,
        ribo_perc_max=30,
        hb_perc_max=10,
        doublet_index=['b'],
        stage_removed={'batch_adaptive_qc': 0},
    )

    assert summary['unique_removed_total'] == 3
    assert summary['removed_low_umi'] == 2
    assert summary['removed_high_mt'] == 1
    assert summary['removed_doublets'] == 1
    assert summary['removed_high_ribo'] == 1
    assert summary['removed_high_hb'] == 0
    assert summary['reason_overlap_cells'] >= 0


def test_scrublet_summary_reports_global_and_batch_rates():
    from modules.qc import _scrublet_summary

    adata = ad.AnnData(
        X=np.ones((4, 1)),
        obs=pd.DataFrame({
            'batch': pd.Categorical(['A', 'A', 'B', 'B']),
            'predicted_doublet': [True, False, False, True],
        }, index=['a', 'b', 'c', 'd']),
        var=pd.DataFrame(index=['gene']),
    )
    summary = _scrublet_summary(adata, 'batch', 0.25)

    assert summary['predicted_doublets'] == 2
    assert summary['doublet_rate'] == 0.5
    assert summary['by_batch']['A']['predicted_doublets'] == 1
    assert summary['by_batch']['B']['threshold'] == 0.25


def test_scdblfinder_summary_uses_classifier_evidence_not_scrublet_simulation():
    from modules.qc import _scdblfinder_summary

    adata = ad.AnnData(
        X=np.ones((4, 1)),
        obs=pd.DataFrame({
            'batch': pd.Categorical(['A', 'A', 'B', 'B']),
            'predicted_doublet': [True, False, False, True],
            'scdblfinder_doublet': [True, False, False, True],
            'scdblfinder_score': [0.92, 0.04, 0.11, 0.88],
        }, index=['a', 'b', 'c', 'd']),
        var=pd.DataFrame(index=['gene']),
    )

    summary = _scdblfinder_summary(adata, 'batch')

    assert summary['method'] == 'scdblfinder'
    assert summary['score_column'] == 'scdblfinder_score'
    assert summary['class_column'] == 'scdblfinder_doublet'
    assert summary['threshold_method'] == 'scdblfinder_model_classification'
    assert summary['threshold'] is None
    assert 'n_simulated_doublets' not in summary
    assert summary['by_batch']['A']['predicted_doublets'] == 1


def test_scrublet_summary_preserves_independent_batch_thresholds():
    from modules.qc import _scrublet_summary

    adata = ad.AnnData(
        X=np.ones((4, 1)),
        obs=pd.DataFrame({
            'batch': pd.Categorical(['A', 'A', 'B', 'B']),
            'predicted_doublet': [True, False, False, True],
        }, index=['a', 'b', 'c', 'd']),
        var=pd.DataFrame(index=['gene']),
    )
    summary = _scrublet_summary(
        adata,
        'batch',
        scrublet_payload={
            'batches': {
                'A': {'threshold': 0.12},
                'B': {'threshold': 0.31},
            },
        },
    )

    assert summary['threshold_by_batch'] == {'A': 0.12, 'B': 0.31}
    assert summary['by_batch']['A']['threshold'] == 0.12
    assert summary['by_batch']['B']['threshold'] == 0.31

def test_scrublet_summary_reports_detected_rate_and_sim_evidence():
    from modules.qc import _scrublet_summary

    adata = ad.AnnData(
        X=np.ones((6, 1)),
        obs=pd.DataFrame({
            'batch': pd.Categorical(['A', 'A', 'A', 'B', 'B', 'B']),
            'predicted_doublet': [True, False, False, False, True, False],
        }, index=['a1', 'a2', 'a3', 'b1', 'b2', 'b3']),
        var=pd.DataFrame(index=['gene']),
    )
    summary = _scrublet_summary(
        adata,
        'batch',
        scrublet_payload={
            'batches': {
                'A': {
                    'threshold': 0.45,
                    'parameters': {'expected_doublet_rate': 0.05, 'n_neighbors': 12},
                    'doublet_scores_sim': [0.1, 0.3, 0.5, 0.7, 0.9],
                },
                'B': {
                    'threshold': 0.60,
                    'parameters': {'expected_doublet_rate': 0.05, 'n_neighbors': 10},
                    'doublet_scores_sim': [0.2, 0.4, 0.6, 0.8, 1.0],
                },
            },
        },
        run_params={
            'expected_doublet_rate': 0.05,
            'sim_doublet_ratio': 2.0,
            'n_prin_comps': 30,
            'random_state': 1234,
            'threshold_method': 'automatic_simulated_distribution',
        },
    )

    import pytest

    assert summary['detected_doublet_rate'] == summary['doublet_rate']
    assert summary['detected_doublet_rate'] == pytest.approx(1 / 3)
    assert summary['threshold_method'] == 'automatic_simulated_distribution'
    a = summary['by_batch']['A']
    assert a['detected_doublet_rate'] == pytest.approx(1 / 3)
    assert a['n_simulated_doublets'] == 5
    assert a['sim_score_quantiles']['max'] == 0.9
    # 0.5, 0.7 and 0.9 are strictly above threshold 0.45 -> 3/5
    assert a['detectable_doublet_fraction'] == 0.6
    # estimate uses the already-rounded detected rate (0.333333 / 0.6)
    assert a['estimated_overall_doublet_rate'] == 0.555555
    assert a['n_neighbors'] == 12
    b = summary['by_batch']['B']
    # only 0.8 and 1.0 above 0.60
    assert b['detectable_doublet_fraction'] == 0.4
    # 1/3 detected rate is far above the 0.005 review cut-off
    assert b['status'] == 'pass'
    assert summary['status'] == 'pass'


def test_scrublet_summary_low_detected_rate_is_review_not_failure():
    from modules.qc import _scrublet_summary

    adata = ad.AnnData(
        X=np.ones((4, 1)),
        obs=pd.DataFrame({
            'batch': pd.Categorical(['A', 'A', 'B', 'B']),
            'predicted_doublet': [True, False, False, False],
        }, index=['a1', 'a2', 'b1', 'b2']),
        var=pd.DataFrame(index=['gene']),
    )
    summary = _scrublet_summary(
        adata,
        'batch',
        scrublet_payload={
            'batches': {
                'A': {'threshold': 0.5, 'doublet_scores_sim': [0.1, 0.9, 0.9, 0.9]},
                'B': {'threshold': 0.5, 'doublet_scores_sim': [0.1, 0.1, 0.1, 0.1]},
            },
        },
        run_params={'threshold_method': 'automatic_simulated_distribution'},
    )

    # A: 1/2 = 0.5 detected -> pass; B: 0/2 -> review with estimated rate None
    assert summary['by_batch']['A']['status'] == 'pass'
    assert summary['by_batch']['B']['status'] == 'review'
    assert summary['by_batch']['B']['detected_doublet_rate'] == 0.0
    # simulated scores never reach the threshold -> no estimate
    assert summary['by_batch']['B']['estimated_overall_doublet_rate'] is None
    # Global detected rate 1/4 = 0.25 is far above the 0.005 review cut-off
    assert summary['status'] == 'pass'


def test_scrublet_summary_unbatched_keeps_top_level_sim_scores():
    from modules.qc import _scrublet_summary

    adata = ad.AnnData(
        X=np.ones((3, 1)),
        obs=pd.DataFrame({
            'predicted_doublet': [True, False, False],
        }, index=['c1', 'c2', 'c3']),
        var=pd.DataFrame(index=['gene']),
    )
    summary = _scrublet_summary(
        adata,
        None,
        doublet_threshold=0.5,
        scrublet_payload={
            'threshold': 0.5,
            'doublet_scores_sim': [0.1, 0.6, 0.8, 0.9],
            'parameters': {'expected_doublet_rate': 0.05, 'n_neighbors': 8},
        },
        run_params={'threshold_method': 'automatic_simulated_distribution'},
    )

    import pytest

    assert summary['n_simulated_doublets'] == 4
    assert summary['detectable_doublet_fraction'] == 0.75  # 0.6/0.8/0.9 above 0.5
    assert summary['expected_doublet_rate'] == 0.05
    assert summary['n_neighbors'] == 8
    assert summary['estimated_overall_doublet_rate'] == pytest.approx((1 / 3) / 0.75)


def test_canonicalize_qc_obs_fields_removes_duplicate_aliases():
    from modules.qc import canonicalize_qc_obs_fields

    adata = ad.AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({
            'total_counts': [10, 20],
            'nUMIs': [10, 20],
            'n_counts': [10, 20],
            'n_genes_by_counts': [3, 4],
            'detected_genes': [3, 4],
            'pct_counts_mt': [1.0, 2.0],
            'mito_perc': [0.01, 0.02],
        }, index=['a', 'b']),
        var=pd.DataFrame(index=['gene']),
    )

    canonicalize_qc_obs_fields(adata)

    assert set(['total_counts', 'n_genes_by_counts', 'pct_counts_mt']).issubset(adata.obs.columns)
    assert not set(['nUMIs', 'n_counts', 'detected_genes', 'mito_perc']) & set(adata.obs.columns)


def test_dimred_embedding_uses_configured_batch_key_not_literal_batch():
    from modules.dimred import _embedding_color_keys

    keys = _embedding_color_keys('sample_id', 'leiden', 'n_genes_by_counts')

    assert keys == ['sample_id', 'leiden', 'n_genes_by_counts']
    assert 'batch' not in keys


def test_batch_cluster_overlap_reports_batch_dominated_clusters():
    from modules.io_utils import batch_cluster_overlap

    adata = ad.AnnData(
        X=np.ones((8, 1)),
        obs=pd.DataFrame({
            'leiden': ['0', '0', '0', '0', '1', '1', '1', '1'],
            'batch': ['A', 'A', 'A', 'A', 'A', 'A', 'B', 'B'],
        }, index=[f'cell_{i}' for i in range(8)]),
        var=pd.DataFrame(index=['GENE1']),
    )

    result = batch_cluster_overlap(adata, 'leiden', 'batch')

    assert result['valid'] is True
    assert result['batch_dominated_cluster_count'] == 1
    assert result['batch_dominated_cell_fraction'] == 0.5
    assert result['warnings']


def test_clustering_resolution_parser_rejects_empty_input_and_deduplicates():
    import pytest
    from modules.clustering import cluster_key_for_method, parse_resolutions

    assert parse_resolutions('0.8,0.8,1.0') == [0.8, 1.0]
    assert cluster_key_for_method('louvain', 0.8) == 'louvain_0.8'
    assert cluster_key_for_method('leiden', 0.8) == 'leiden_0.8'
    with pytest.raises(ValueError):
        parse_resolutions('')


def test_metric_subset_is_bounded_and_keeps_each_cluster():
    from modules.clustering import metric_subset

    rep = np.arange(60, dtype=float).reshape(20, 3)
    labels = np.asarray(['a'] * 10 + ['b'] * 6 + ['c'] * 4)
    subset_rep, subset_labels = metric_subset(rep, labels, max_cells=8)

    assert subset_rep.shape == (8, 3)
    assert set(subset_labels) == {'a', 'b', 'c'}


def test_dimred_schema_exposes_batch_key():
    from modules.schemas import PARAM_SCHEMAS

    dimred_keys = [item['key'] for item in PARAM_SCHEMAS['dimred']]

    assert 'batch_key' in dimred_keys


def test_dimred_defaults_use_hvg_elbow_and_empty_batch_key():
    from modules.schemas import PARAM_SCHEMAS

    dimred = {item['key']: item for item in PARAM_SCHEMAS['dimred']}

    assert dimred['n_comps']['default'] == 25
    assert dimred['auto_n_comps']['default'] == 'elbow'
    assert dimred['batch_key']['default'] == ''
    assert dimred['pca_hvg_only']['default'] is True


def test_dimred_elbow_selection_never_exceeds_small_candidate_space():
    from modules.dimred import _select_elbow_n_comps

    assert _select_elbow_n_comps(np.ones(2), 2) == 2
    assert _select_elbow_n_comps(np.ones(3), 3) == 3
    assert _select_elbow_n_comps(np.ones(10), 10) == 5


def test_dimred_pca_fingerprint_changes_with_coordinates():
    from modules.dimred import _pca_fingerprint

    first = np.arange(30, dtype=float).reshape(6, 5)
    second = first.copy()
    second[0, 0] += 1

    assert _pca_fingerprint(first) != _pca_fingerprint(second)


def test_restore_scanpy_qc_percentages_from_omicverse_fractions():
    from modules.io_utils import restore_scanpy_qc_percentages

    adata = ad.AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({
            'mito_perc': [0.0203, 0.1045],
            'ribo_perc': [0.368, 0.250],
            'hb_perc': [0.0, 0.0136],
            'pct_counts_mt': [0.0203, 0.1045],
            'pct_counts_ribo': [0.368, 0.250],
            'pct_counts_hb': [0.0, 0.0136],
        }, index=['cell_a', 'cell_b']),
        var=pd.DataFrame(index=['GENE1']),
    )

    returned = restore_scanpy_qc_percentages(adata)

    assert returned is adata
    np.testing.assert_allclose(adata.obs['pct_counts_mt'], [2.03, 10.45])
    np.testing.assert_allclose(adata.obs['pct_counts_ribo'], [36.8, 25.0])
    np.testing.assert_allclose(adata.obs['pct_counts_hb'], [0.0, 1.36])
    np.testing.assert_allclose(adata.obs['mito_perc'], [0.0203, 0.1045])

    # UI thresholds are percentages and compare directly to pct_counts_*.
    assert (adata.obs['pct_counts_ribo'] <= 30.0).tolist() == [False, True]
    assert (adata.obs['pct_counts_hb'] <= 1.0).tolist() == [True, False]


def test_restore_scanpy_qc_percentages_accepts_existing_percent_aliases():
    from modules.io_utils import restore_scanpy_qc_percentages

    adata = ad.AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({
            'mito_perc': [2.0, 12.0],
            'pct_counts_mt': [2.0, 12.0],
        }, index=['cell_a', 'cell_b']),
        var=pd.DataFrame(index=['GENE1']),
    )

    restore_scanpy_qc_percentages(adata)

    np.testing.assert_allclose(adata.obs['pct_counts_mt'], [2.0, 12.0])


def test_cell_cycle_copy_is_normalized_even_when_qc_log_columns_exist():
    import scanpy as sc
    from modules.qc import _normalized_cell_cycle_copy

    counts = np.asarray([
        [10.0, 0.0, 5.0],
        [1.0, 7.0, 2.0],
    ])
    adata = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(
            {'log1p_total_counts': np.log1p(counts.sum(axis=1))},
            index=['cell_a', 'cell_b'],
        ),
        var=pd.DataFrame(index=['GENE1', 'GENE2', 'GENE3']),
    )

    normalized = _normalized_cell_cycle_copy(adata, sc)

    np.testing.assert_allclose(adata.X, counts)
    np.testing.assert_allclose(
        np.expm1(normalized.X).sum(axis=1), np.full(adata.n_obs, 1e4),
        rtol=1e-6,
    )
    assert normalized.uns['log1p']['base'] is None


def test_qc_reassess_uses_percentage_scale_for_legacy_qc_output(tmp_path):
    from modules.qc_reassess import QCReassessAnalysis

    adata = ad.AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({
            'leiden': pd.Categorical(['0', '0', '1', '1']),
            'mito_perc': [0.20, 0.20, 0.02, 0.02],
            'pct_counts_mt': [0.20, 0.20, 0.02, 0.02],
            'total_counts': [1000.0] * 4,
            'n_genes_by_counts': [500] * 4,
            'predicted_doublet': [False] * 4,
        }, index=[f'cell_{i}' for i in range(4)]),
        var=pd.DataFrame(index=['GENE1', 'GENE2']),
    )
    input_path = tmp_path / 'legacy_qc.h5ad'
    adata.write_h5ad(input_path)
    module = QCReassessAnalysis(
        project_dir=str(tmp_path),
        params={'mt_threshold': 15.0, 'min_cells_per_cluster': 1},
        progress_callback=lambda *_: None,
    )
    module.save_matplotlib_figure = lambda *_args, **_kwargs: []

    result = module.run(str(input_path))
    table = pd.read_csv(tmp_path / 'results' / 'low_quality_clusters.csv')
    cluster_zero = table.loc[table['cluster'].astype(str) == '0'].iloc[0]

    assert result['summary']['low_quality_clusters'] == ['0']
    assert cluster_zero['mean_pct_mt'] == 20.0
    assert bool(cluster_zero['low_quality']) is True


def test_qc_reassess_does_not_treat_mean_score_as_doublet_fraction(tmp_path):
    from modules.qc_reassess import QCReassessAnalysis

    adata = ad.AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({
            'leiden': pd.Categorical(['0', '0', '1', '1']),
            'predicted_doublet': [False, False, False, False],
            'doublet_score': [0.9, 0.8, 0.1, 0.2],
            'total_counts': [1000.0] * 4,
            'n_genes_by_counts': [500] * 4,
        }, index=[f'cell_{i}' for i in range(4)]),
        var=pd.DataFrame(index=['GENE1', 'GENE2']),
    )
    input_path = tmp_path / 'scores.h5ad'
    adata.write_h5ad(input_path)
    module = QCReassessAnalysis(
        project_dir=str(tmp_path),
        params={'doublet_threshold': 0.3, 'min_cells_per_cluster': 1},
        progress_callback=lambda *_: None,
    )
    module.save_matplotlib_figure = lambda *_args, **_kwargs: []

    result = module.run(str(input_path))
    table = pd.read_csv(tmp_path / 'results' / 'low_quality_clusters.csv')
    table['cluster'] = table['cluster'].astype(str)
    table = table.set_index('cluster')

    assert result['summary']['doublet_fraction_source'] == 'predicted_doublet'
    assert table.loc['0', 'doublet_fraction'] == 0.0
    assert table.loc['0', 'mean_doublet_score'] == 0.85
    assert 'high_doublet' not in str(table.loc['0', 'low_reasons'])
