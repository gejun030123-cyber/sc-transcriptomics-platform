# tests/test_bulk_qc_helpers.py
import os
import numpy as np
import pandas as pd
import pytest


def test_gini_uniform():
    """完全均匀分布 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([100.0, 100.0, 100.0, 100.0])
    assert abs(_gini(vals)) < 1e-6


def test_gini_concentrated():
    """一个值占主导 → Gini 接近 1"""
    from modules.bulk_qc import _gini
    vals = np.array([0.0] * 20 + [1000.0])
    assert _gini(vals) > 0.9


def test_gini_all_zeros():
    """全零输入 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([0.0, 0.0, 0.0])
    assert _gini(vals) == 0.0


def test_gini_empty():
    """空数组 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([])
    assert _gini(vals) == 0.0


def test_bulk_mt_detection_recognizes_versioned_human_ensembl_ids():
    """Ensembl-only human counts must not bypass the MT% QC threshold."""
    import anndata as ad
    from modules.bulk_qc import _detect_mitochondrial_genes

    adata = ad.AnnData(
        X=np.ones((2, 3)),
        var=pd.DataFrame(index=[
            'ENSG00000198888.2',  # MT-ND1
            'ENSG00000141510',    # TP53
            'ENSG00000111640',    # GAPDH
        ]),
    )

    result = _detect_mitochondrial_genes(adata)

    assert result['mask'].tolist() == [True, False, False]
    assert result['n_mitochondrial_genes'] == 1
    assert result['sources'] == {'human_ensembl_id': 1}


def test_bulk_mt_detection_recognizes_chromosome_annotation_without_symbol():
    """chrM/MT annotation is a species-agnostic fallback when symbols are absent."""
    import anndata as ad
    from modules.bulk_qc import _detect_mitochondrial_genes

    adata = ad.AnnData(
        X=np.ones((2, 3)),
        var=pd.DataFrame(
            {'Chr': ['chrM', '1', 'MT']},
            index=['feature_1', 'feature_2', 'feature_3'],
        ),
    )

    result = _detect_mitochondrial_genes(adata)

    assert result['mask'].tolist() == [True, False, True]
    assert result['n_mitochondrial_genes'] == 2
    assert result['sources'] == {'chromosome': 2}


def test_bulk_qc_applies_mt_threshold_for_ensembl_only_input(tmp_path):
    """A high-MT Ensembl-only sample is actually filtered, not silently kept."""
    import anndata as ad
    from modules.bulk_qc import BulkQCAnalysis

    adata = ad.AnnData(
        X=np.array([
            [30, 70, 50],   # MT-ND1 = 20%; threshold is inclusive
            [25, 75, 50],
            [5, 100, 30],
            [50, 50, 10],   # MT-ND1 > 20%; must fail
        ], dtype=int),
        obs=pd.DataFrame(index=['S1', 'S2', 'S3', 'S4']),
        var=pd.DataFrame(index=[
            'ENSG00000198888', 'ENSG00000141510', 'ENSG00000111640',
        ]),
    )
    input_path = tmp_path / 'ensembl_only_counts.h5ad'
    adata.write_h5ad(input_path)

    result = BulkQCAnalysis(
        str(tmp_path), {
            'input_measurement': 'raw_counts',
            'min_counts': 0, 'min_genes': 0,
            'max_mt_pct': 20, 'max_ribo_pct': 100,
            'min_sample_expr': 0, 'detect_outliers': False,
        }, lambda *_: None,
    ).run(str(input_path))

    mitochondrial_qc = result['summary']['mitochondrial_qc']
    assert result['summary']['samples_after'] == 3
    assert mitochondrial_qc['status'] == 'identified'
    assert mitochondrial_qc['threshold_applied'] is True
    assert mitochondrial_qc['sources']['human_ensembl_id'] == 1
    metrics = pd.read_csv(tmp_path / 'results' / 'bulk_qc_sample_metrics.csv')
    assert metrics.loc[metrics['sample'] == 'S4', 'fail_reasons'].item() == 'mt_pct>20.0'


def test_bulk_mt_detection_covers_release_drifted_mt_tp_id():
    """tRNA stable IDs drift between GTF releases; the union set must cover both."""
    import anndata as ad
    from modules.bulk_qc import _detect_mitochondrial_genes

    adata = ad.AnnData(
        X=np.ones((2, 2)),
        var=pd.DataFrame(index=['ENSG00000210196', 'ENSG00000141510']),
    )

    result = _detect_mitochondrial_genes(adata)

    assert result['mask'].tolist() == [True, False]
    assert result['sources'] == {'human_ensembl_id': 1}


def test_bulk_qc_flags_identified_mt_genes_with_zero_reads(tmp_path):
    """Correctly identified MT genes without reads must not be reported as clean 0%."""
    import anndata as ad
    from modules.bulk_qc import BulkQCAnalysis

    adata = ad.AnnData(
        X=np.array([
            [0, 100, 50],
            [0, 120, 60],
            [0, 90, 45],
            [0, 110, 55],
        ], dtype=int),
        obs=pd.DataFrame(index=['S1_1', 'S1_2', 'S2_1', 'S2_2']),
        var=pd.DataFrame(
            {'gene_name': ['MT-ND1', 'GAPDH', 'ACTB']},
            index=['ENSG00000198888', 'ENSG00000111640', 'ENSG00000149257'],
        ),
    )
    input_path = tmp_path / 'zero_mt_counts.h5ad'
    adata.write_h5ad(input_path)

    result = BulkQCAnalysis(
        str(tmp_path), {
            'input_measurement': 'raw_counts',
            'min_counts': 0, 'min_genes': 0,
            'max_mt_pct': 20, 'max_ribo_pct': 100,
            'min_sample_expr': 0, 'detect_outliers': False,
        }, lambda *_: None,
    ).run(str(input_path))

    summary = result['summary']
    mitochondrial_qc = summary['mitochondrial_qc']
    # The MT feature is identified, but an all-zero block cannot support MT QC.
    assert mitochondrial_qc['status'] == 'identified_no_reads'
    assert mitochondrial_qc['informative'] is False
    assert mitochondrial_qc['threshold_applied'] is False
    assert mitochondrial_qc['total_mt_counts'] == 0
    assert mitochondrial_qc['mt_genes_with_counts'] == 0
    assert mitochondrial_qc['n_mitochondrial_genes'] == 1
    # 0% stays the numeric truth; it is flagged instead of presented as evidence.
    assert summary['median_mt_pct'] == 0.0
    assert any('MT counts 均为 0' in warning for warning in summary['warnings'])
    # QC intentionally keeps every gene, so genes_removed=0 must be explained.
    assert any('genes_removed=0' in warning for warning in summary['warnings'])
    assert summary['gene_expression_filter']['n_genes_removed_in_qc'] == 0
    assert summary['gene_expression_filter']['applied_in_qc'] is False


def test_outlier_diagnostic_explains_invariant_qc_metric_threshold():
    """An invariant metric yields a null threshold; the reason must be recorded."""
    from modules.bulk_qc import _detect_within_replicate_outliers

    samples = ['Ctr_1', 'Ctr_2', 'Ctr_3', 'Treat_1', 'Treat_2', 'Treat_3']
    groups = ['Ctr'] * 3 + ['Treat'] * 3
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(6, 4))
    corr = np.full((6, 6), 0.99)
    np.fill_diagonal(corr, 1.0)
    qc_metrics = pd.DataFrame(
        {
            'mt_pct': [0.0] * 6,
            'ribo_pct': [4.0, 4.2, 3.9, 5.0, 5.1, 4.9],
        },
        index=samples,
    )

    _, details = _detect_within_replicate_outliers(
        coords, corr, samples, groups, qc_metrics=qc_metrics,
    )

    assert details['qc_metric_deviation_thresholds']['mt_pct'] is None
    assert details['qc_metric_deviation_threshold_notes']['mt_pct'] == 'invariant_metric_no_threshold'
    assert 'mt_pct' in details['invariant_qc_metrics']
    assert details['qc_metric_deviation_thresholds']['ribo_pct'] is not None
    assert 'ribo_pct' not in details['qc_metric_deviation_threshold_notes']


def test_infer_groups_dash():
    """样本名含 '-' → 取前缀"""
    from modules.bulk_qc import _infer_groups
    names = ['ctrl-1', 'ctrl-2', 'hmc3-1', 'hmc3-2']
    assert _infer_groups(names) == ['ctrl', 'ctrl', 'hmc3', 'hmc3']


def test_infer_groups_underscore():
    """样本名含 '_' → 取前缀"""
    from modules.bulk_qc import _infer_groups
    names = ['ctrl_1', 'rapa_3']
    assert _infer_groups(names) == ['ctrl', 'rapa']


def test_infer_groups_no_sep():
    """无分隔符 → 样本名本身"""
    from modules.bulk_qc import _infer_groups
    names = ['sampleA', 'sampleB']
    assert _infer_groups(names) == ['sampleA', 'sampleB']


def test_infer_multifactor_groups_preserves_combined_group():
    """Ctr_B_1 应默认分到 Ctr_B，不应与 Ctr_En 合并。"""
    from modules.bulk_qc import _infer_groups

    names = ['Ctr_B_1', 'Ctr_B_2', 'Ctr_En_1', 'Ctr_En_2']
    assert _infer_groups(names) == ['Ctr_B', 'Ctr_B', 'Ctr_En', 'Ctr_En']


def test_multifactor_group_candidates_include_combined_and_factors():
    from modules.io_utils import infer_sample_group_candidates

    names = [
        'Ctr_B_1', 'Ctr_B_2', 'Ctr_En_1', 'Ctr_En_2',
        'PEA_B_1', 'PEA_B_2', 'PEA_En_1', 'PEA_En_2',
    ]
    candidates = infer_sample_group_candidates(names)

    assert [c['key'] for c in candidates] == ['combined', 'factor_1', 'factor_2']
    assert candidates[0]['values'] == ['Ctr_B', 'Ctr_En', 'PEA_B', 'PEA_En']
    assert candidates[0]['group_sizes'] == {
        'Ctr_B': 2, 'Ctr_En': 2, 'PEA_B': 2, 'PEA_En': 2,
    }
    assert candidates[1]['values'] == ['Ctr', 'PEA']
    assert candidates[2]['values'] == ['B', 'En']


def test_simple_group_names_do_not_duplicate_candidate():
    from modules.io_utils import infer_sample_group_candidates

    candidates = infer_sample_group_candidates(['Ctrl_1', 'Ctrl_2', 'Drug_1', 'Drug_2'])
    assert len(candidates) == 1
    assert candidates[0]['values'] == ['Ctrl', 'Drug']


def test_unstructured_sample_names_have_no_auto_group():
    from modules.io_utils import infer_sample_group_candidates

    assert infer_sample_group_candidates(['sampleA', 'sampleB']) == []


def test_obs_columns_api_returns_multifactor_candidates(test_project):
    from app import create_app
    from config import Config

    path = os.path.join(Config.uploads_dir(test_project), 'bulk.tsv')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(
            'gene\tCtr_B_1\tCtr_B_2\tCtr_En_1\tCtr_En_2\tPEA_B_1\tPEA_B_2\tPEA_En_1\tPEA_En_2\n'
            'G1\t1\t2\t3\t4\t5\t6\t7\t8\n'
            'G2\t8\t7\t6\t5\t4\t3\t2\t1\n'
        )

    app = create_app()
    app.config['TESTING'] = True
    response = app.test_client().get('/api/obs-columns', query_string={'file_path': path})
    payload = response.get_json()

    assert response.status_code == 200
    assert payload['sample_groups']['auto_group']['values'] == ['Ctr_B', 'Ctr_En', 'PEA_B', 'PEA_En']
    assert len(payload['sample_groups']['auto_group_candidates']) == 3


def test_obs_columns_api_separates_group_batch_and_cluster_suggestions(test_project):
    import anndata
    import pandas as pd
    from app import create_app
    from config import Config

    n_obs = 24
    obs = pd.DataFrame({
        # Deliberately put technical columns first: this reproduces the old
        # physical-order bug in the web form.
        'barcode': [f'cell_{i}' for i in range(n_obs)],
        'batch': pd.Categorical(['B1'] * 12 + ['B2'] * 12),
        'sample': pd.Categorical([f'S{i // 6 + 1}' for i in range(n_obs)]),
        'condition': pd.Categorical(['Control'] * 12 + ['Treatment'] * 12),
        'leiden': pd.Categorical(['0'] * 12 + ['1'] * 12),
    }, index=[f'cell_{i}' for i in range(n_obs)])
    adata = anndata.AnnData(
        np.ones((n_obs, 4)),
        obs=obs,
        var=pd.DataFrame(index=['G1', 'G2', 'G3', 'G4']),
    )
    path = os.path.join(Config.uploads_dir(test_project), 'grouping.h5ad')
    adata.write_h5ad(path)

    app = create_app()
    app.config['TESTING'] = True
    response = app.test_client().get(
        '/api/obs-columns',
        query_string={'file_path': path, 'module_name': 'deg'},
    )
    payload = response.get_json()

    assert response.status_code == 200
    assert payload['suggestions']['groupby'] == 'condition'
    assert payload['suggestions']['color_by'] == 'condition'
    assert payload['suggestions']['batch_key'] == 'batch'
    assert payload['suggestions']['cluster_key'] == 'leiden'
    assert payload['suggestions']['sample_key'] == 'sample'
    assert payload['grouping_candidates'][0] == 'condition'


def test_detect_outliers_no_outlier():
    """正常聚类数据 → 无离群点"""
    from modules.bulk_qc import _detect_outliers_mahal
    np.random.seed(42)
    coords = np.random.randn(20, 3)
    outliers = _detect_outliers_mahal(coords, [f's{i}' for i in range(20)])
    assert len(outliers) == 0


def test_detect_outliers_with_outlier():
    """含明显离群点 → 检测到"""
    from modules.bulk_qc import _detect_outliers_mahal
    np.random.seed(42)
    coords = np.random.randn(20, 3)
    coords = np.vstack([coords, [[50, 50, 50]]])  # 离群点
    samples = [f's{i}' for i in range(21)]
    outliers = _detect_outliers_mahal(coords, samples)
    assert 's20' in outliers


def test_detect_outliers_uses_later_principal_components():
    """A PC4-only anomaly must not disappear behind a PC1/PC2-only check."""
    from modules.bulk_qc import _detect_outliers_mahal

    rng = np.random.default_rng(7)
    coords = rng.normal(size=(24, 5))
    coords[-1, 3] = 30.0
    samples = [f's{i}' for i in range(len(coords))]

    outliers, details = _detect_outliers_mahal(coords, samples, return_details=True)

    assert details['n_components'] >= 4
    assert 's23' in outliers


def test_detect_outliers_too_few():
    """样本数 < 4 → 返回空"""
    from modules.bulk_qc import _detect_outliers_mahal
    coords = np.array([[0, 0], [1, 1], [2, 2]])
    outliers = _detect_outliers_mahal(coords, ['a', 'b', 'c'])
    assert outliers == []


def test_within_replicate_outlier_diagnostic_does_not_flag_coherent_factor_stratum():
    """A B/En-wide shift must not turn every B replicate into a global-PCA outlier."""
    from modules.bulk_qc import _detect_within_replicate_outliers

    rng = np.random.default_rng(19)
    groups = ['Ctr_B'] * 3 + ['Ctr_En'] * 3 + ['PEA_B'] * 3 + ['PEA_En'] * 3
    samples = [f'S{i}' for i in range(len(groups))]
    centers = {'Ctr_B': 10.0, 'Ctr_En': -10.0, 'PEA_B': 11.0, 'PEA_En': -9.0}
    coords = np.vstack([rng.normal(loc=centers[group], scale=0.08, size=4) for group in groups])
    corr = np.full((len(groups), len(groups)), 0.35)
    np.fill_diagonal(corr, 1.0)
    for i, group_i in enumerate(groups):
        for j, group_j in enumerate(groups):
            if group_i == group_j and i != j:
                corr[i, j] = 0.98

    table, details = _detect_within_replicate_outliers(coords, corr, samples, groups)

    assert not table['outlier_flag'].any()
    assert details['method'] == 'within_replicate_group_agreement'


def test_within_replicate_outlier_diagnostic_flags_one_aberrant_sample_not_its_group():
    from modules.bulk_qc import _detect_within_replicate_outliers

    rng = np.random.default_rng(23)
    groups = ['Ctr_B'] * 3 + ['Ctr_En'] * 3 + ['PEA_B'] * 3 + ['PEA_En'] * 3
    samples = [f'S{i}' for i in range(len(groups))]
    coords = rng.normal(scale=0.12, size=(len(groups), 5))
    # One Ctr_B sample is far from its two own replicates; the other strata
    # retain their own coherent local centers.
    for group_index, group in enumerate(('Ctr_B', 'Ctr_En', 'PEA_B', 'PEA_En')):
        indices = [i for i, value in enumerate(groups) if value == group]
        coords[indices] += group_index * 4.0
    coords[2] += 12.0
    corr = np.full((len(groups), len(groups)), 0.4)
    np.fill_diagonal(corr, 1.0)
    for i, group_i in enumerate(groups):
        for j, group_j in enumerate(groups):
            if group_i == group_j and i != j:
                corr[i, j] = 0.985
    corr[2, 0] = corr[0, 2] = 0.35
    corr[2, 1] = corr[1, 2] = 0.35

    table, _ = _detect_within_replicate_outliers(coords, corr, samples, groups)

    flagged = table.loc[table['outlier_flag'], 'sample_id'].tolist()
    assert flagged == ['S2']


def test_log2_normalization_omits_library_size_comparison(tmp_path):
    """A pure log2 transform must not export an unchanged library-size panel."""
    from modules.bulk_normalize import BulkNormalizeAnalysis

    expression = pd.DataFrame(
        {
            'Ctr_1': [118.5, 202.2, 3.1, 0.0, 11.4],
            'Ctr_2': [121.8, 190.4, 2.9, 0.0, 10.8],
            'Treat_1': [225.7, 111.3, 8.8, 1.0, 19.2],
            'Treat_2': [219.2, 108.7, 7.9, 1.2, 18.5],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5'],
    )
    input_path = tmp_path / 'continuous_expression.tsv'
    expression.to_csv(input_path, sep='\t', index_label='gene')
    module = BulkNormalizeAnalysis(
        str(tmp_path), {
            'method': 'log2', 'input_measurement': 'continuous_expression',
            'min_expr_samples': 0,
        }, lambda *_: None,
    )

    result = module.run(str(input_path))

    stems = {os.path.basename(item['file_path']) for item in result['result_files']}
    assert not any(stem.startswith('bulk_norm_libsize') for stem in stems)
    assert any(stem.startswith('bulk_norm_boxplot_compare') for stem in stems)
    assert result['summary']['n_genes_after_filter'] == 5
    assert result['summary']['max_zero_pct_filter_enabled'] is False
    assert result['summary']['max_zero_pct_semantics'].startswith('0 disables')


def test_bulk_normalize_requires_explicit_scale_before_log2_for_ambiguous_values(tmp_path):
    """Non-integer values alone cannot prove that another log is safe."""
    from modules.bulk_qc import BulkQCAnalysis
    from modules.bulk_normalize import BulkNormalizeAnalysis

    expression = pd.DataFrame(
        [[2.2, 2.4, 3.1, 3.3], [4.8, 5.0, 4.9, 5.1]],
        index=['G1', 'G2'], columns=['Ctr_1', 'Ctr_2', 'Treat_1', 'Treat_2'],
    )
    input_path = tmp_path / 'unknown_expression.tsv'
    expression.to_csv(input_path, sep='\t', index_label='gene')

    with pytest.raises(ValueError, match='无法仅靠数值区分'):
        BulkNormalizeAnalysis(
            str(tmp_path), {'method': 'log2', 'min_expr_samples': 0}, lambda *_: None,
        ).run(str(input_path))
    with pytest.raises(ValueError, match='无法仅靠数值区分'):
        BulkQCAnalysis(
            str(tmp_path), {'min_sample_expr': 0}, lambda *_: None,
        ).run(str(input_path))


def test_bulk_normalize_applies_nonzero_zero_filter_even_when_expression_filter_is_off(tmp_path):
    """A positive zero threshold is active independently of min_expr_samples."""
    from modules.bulk_normalize import BulkNormalizeAnalysis

    expression = pd.DataFrame(
        {
            'Ctr_1': [5.1, 3.2, 0.0], 'Ctr_2': [5.2, 3.4, 0.0],
            'Treat_1': [7.1, 2.8, 1.0], 'Treat_2': [7.2, 2.9, 1.1],
        }, index=['G1', 'G2', 'Zero_in_half'],
    )
    input_path = tmp_path / 'linear_values.tsv'
    expression.to_csv(input_path, sep='\t', index_label='gene')

    result = BulkNormalizeAnalysis(
        str(tmp_path), {
            'method': 'log2', 'input_measurement': 'continuous_expression',
            'min_expr_samples': 0, 'max_zero_pct': 25,
        }, lambda *_: None,
    ).run(str(input_path))

    assert result['summary']['n_genes_after_filter'] == 2
    assert result['summary']['max_zero_pct_filter_enabled'] is True


def test_bulk_normalize_preserves_explicit_log_expression_without_double_log(tmp_path):
    """An explicitly declared log matrix must round-trip through ``none``."""
    import anndata
    from modules.bulk_normalize import BulkNormalizeAnalysis

    original = np.array([
        [2.0, 3.0, 4.0, 5.0],
        [2.2, 3.2, 4.2, 5.2],
        [3.0, 2.0, 5.0, 4.0],
        [3.2, 2.2, 5.2, 4.2],
    ])
    input_path = tmp_path / 'already_log.h5ad'
    anndata.AnnData(
        original,
        obs=pd.DataFrame(index=['Ctr_1', 'Ctr_2', 'Treat_1', 'Treat_2']),
        var=pd.DataFrame(index=['G1', 'G2', 'G3', 'G4']),
    ).write_h5ad(input_path)

    result = BulkNormalizeAnalysis(
        str(tmp_path), {
            'method': 'none', 'input_measurement': 'log_transformed',
            'min_expr_samples': 0,
        }, lambda *_: None,
    ).run(str(input_path))

    output = anndata.read_h5ad(result['output_adata'])
    np.testing.assert_allclose(output.X, original)
    assert result['summary']['input_measurement'] == 'log_transformed'
    assert result['summary']['method'] == 'none'


def test_bulk_qc_normalize_and_pca_share_the_same_log_pca_convention(tmp_path):
    """The QC and normalized-output PCA spectra must match for one matrix."""
    from modules.bulk_qc import BulkQCAnalysis
    from modules.bulk_normalize import BulkNormalizeAnalysis
    from modules.bulk_pca import BulkPCAAnalysis

    rng = np.random.default_rng(71)
    base = rng.gamma(shape=2.5, scale=5.0, size=(40, 6))
    base[:, 3:] *= np.linspace(1.4, 3.0, 40)[:, None]
    expression = pd.DataFrame(
        base,
        index=[f'G{i}' for i in range(base.shape[0])],
        columns=['Ctr_1', 'Ctr_2', 'Ctr_3', 'Treat_1', 'Treat_2', 'Treat_3'],
    )
    input_path = tmp_path / 'linear_expression.tsv'
    expression.to_csv(input_path, sep='\t', index_label='gene')

    qc = BulkQCAnalysis(
        str(tmp_path), {
            'input_measurement': 'continuous_expression', 'min_sample_expr': 0,
            'detect_outliers': False,
        }, lambda *_: None,
    ).run(str(input_path))
    normalized = BulkNormalizeAnalysis(
        str(tmp_path), {
            'method': 'log2', 'min_expr_samples': 0,
        }, lambda *_: None,
    ).run(qc['output_adata'])
    pca = BulkPCAAnalysis(str(tmp_path), {'n_comps': 4}, lambda *_: None).run(
        normalized['output_adata']
    )

    norm_pc1 = normalized['summary']['pca_compare_variance_ratio']['output'][0] * 100
    assert qc['summary']['pca_preprocessing']['feature_scaling'] == 'none'
    assert normalized['summary']['pca_preprocessing']['feature_scaling'] == 'none'
    assert pca['summary']['pca_preprocessing']['feature_scaling'] == 'none'
    assert qc['summary']['pc1_variance_pct'] == pytest.approx(norm_pc1, abs=0.01)
    assert pca['summary']['pc1_variance_pct'] == pytest.approx(norm_pc1, abs=0.01)


def test_bulk_pca_uses_auto_factor_marker_and_exports_scores(tmp_path):
    from modules.bulk_pca import BulkPCAAnalysis

    sample_names = [
        'Ctr_B_1', 'Ctr_B_2', 'Ctr_En_1', 'Ctr_En_2',
        'PEA_B_1', 'PEA_B_2', 'PEA_En_1', 'PEA_En_2',
    ]
    rng = np.random.default_rng(13)
    expression = pd.DataFrame(
        rng.poisson(lam=40, size=(12, len(sample_names))),
        index=[f'G{i}' for i in range(12)], columns=sample_names,
    )
    input_path = tmp_path / 'multifactor_counts.tsv'
    expression.to_csv(input_path, sep='\t', index_label='gene')
    module = BulkPCAAnalysis(
        str(tmp_path),
        {'n_comps': 5, 'color_by': '_auto_group_', 'batch_by': '第2因素：B, En'},
        lambda *_: None,
    )

    result = module.run(str(input_path))

    summary = result['summary']
    assert summary['color_by_used'] == '_auto_group'
    assert summary['batch_by_used'] == '_auto_factor2'
    scores_path = next(
        item['file_path'] for item in result['result_files']
        if item['label'] == 'PCA sample scores and metadata'
    )
    scores = pd.read_csv(scores_path)
    assert {'sample_id', 'PC1', 'PC3', 'color_group', 'marker_group', 'factor1', 'factor2'} <= set(scores.columns)
    assert set(scores['marker_group']) == {'B', 'En'}


def test_continuous_bulk_qc_uses_total_expression_terms(tmp_path):
    from modules.bulk_qc import BulkQCAnalysis

    expression = pd.DataFrame(
        {
            'Ctr_1': [118.5, 202.2, 3.1, 0.0, 11.4],
            'Ctr_2': [121.8, 190.4, 2.9, 0.0, 10.8],
            'Treat_1': [225.7, 111.3, 8.8, 1.0, 19.2],
            'Treat_2': [219.2, 108.7, 7.9, 1.2, 18.5],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5'],
    )
    input_path = tmp_path / 'continuous_expression.tsv'
    expression.to_csv(input_path, sep='\t', index_label='gene')
    module = BulkQCAnalysis(
        str(tmp_path), {
            'input_measurement': 'continuous_expression',
            'min_sample_expr': 0, 'detect_outliers': False,
        }, lambda *_: None,
    )

    result = module.run(str(input_path))

    assert result['summary']['input_measurement'] == 'continuous_expression'
    assert 'median_total_expression' in result['summary']
    assert 'median_library_size' not in result['summary']
    metrics_path = next(item['file_path'] for item in result['result_files'] if item['label'] == '样本 QC 指标')
    assert 'total_expression' in pd.read_csv(metrics_path).columns
    assert result['summary']['qc_transform_for_correlation_and_pca'].startswith('log2(')


def test_log_transformed_bulk_qc_does_not_apply_a_second_log(tmp_path):
    import anndata
    from modules.bulk_qc import BulkQCAnalysis

    adata = anndata.AnnData(
        np.log2(np.array([
            [101.0, 99.0, 3.0, 0.0, 9.0],
            [100.0, 98.0, 4.0, 0.0, 8.0],
            [210.0, 105.0, 8.0, 1.0, 18.0],
            [205.0, 102.0, 9.0, 1.0, 17.0],
        ]) + 1.0),
        obs=pd.DataFrame(index=['Ctr_1', 'Ctr_2', 'Treat_1', 'Treat_2']),
        var=pd.DataFrame(index=['G1', 'G2', 'G3', 'G4', 'G5']),
    )
    adata.uns['normalization'] = {'is_log_transformed': True}
    input_path = tmp_path / 'log_expression.h5ad'
    adata.write_h5ad(input_path)

    result = BulkQCAnalysis(
        str(tmp_path), {'min_sample_expr': 0, 'detect_outliers': False}, lambda *_: None,
    ).run(str(input_path))

    assert result['summary']['input_measurement'] == 'log_transformed'
    assert result['summary']['qc_transform_for_correlation_and_pca'] == (
        'input log-transformed expression (no second log)'
    )
    assert 'median_total_transformed_expression' in result['summary']


def test_bulk_qc_uses_combined_replicate_groups_and_factor_aware_correlation_summary(tmp_path):
    import anndata
    from modules.bulk_qc import BulkQCAnalysis

    treatments = ('Ctr', 'NH4Cl', 'PEA', 'TMAO')
    sample_names = [
        f'{treatment}_{factor}_{replicate}'
        for treatment in treatments for factor in ('B', 'En') for replicate in range(1, 4)
    ]
    rng = np.random.default_rng(41)
    expression = np.empty((30, len(sample_names)), dtype=float)
    # Make factor 2 dominate globally while retaining tight repeats inside
    # each treatment x factor group: this must not become 12 false outliers.
    for treatment in treatments:
        for factor in ('B', 'En'):
            base = rng.uniform(10.0, 30.0, size=30)
            if factor == 'B':
                base[:12] += 90.0
            indices = [index for index, name in enumerate(sample_names)
                       if name.startswith(f'{treatment}_{factor}_')]
            expression[:, indices] = base[:, None] + rng.normal(scale=0.08, size=(30, len(indices)))
    frame = pd.DataFrame(expression, index=[f'G{i}' for i in range(30)], columns=sample_names)
    input_path = tmp_path / 'multifactor_continuous.tsv'
    frame.to_csv(input_path, sep='\t', index_label='gene')

    result = BulkQCAnalysis(
        str(tmp_path), {
            'input_measurement': 'continuous_expression',
            'min_sample_expr': 0, 'detect_outliers': True,
        }, lambda *_: None,
    ).run(str(input_path))

    assert result['summary']['outlier_samples'] == []
    assert result['summary']['outlier_detection']['replicate_group_column'] == '_auto_group'
    assert {'factor2_within:B', 'factor2_within:En', 'factor2_between:B_vs_En'} <= set(
        result['summary']['correlation_medians']
    )
    outlier_path = next(
        item['file_path'] for item in result['result_files']
        if item['label'] == '重复组内离群诊断（相关性、PCA 与 QC 指标）'
    )
    outlier_table = pd.read_csv(outlier_path)
    assert set(outlier_table['replicate_group']) == {
        f'{treatment}_{factor}' for treatment in treatments for factor in ('B', 'En')
    }
    output = anndata.read_h5ad(result['output_adata'])
    assert {'_auto_group', '_auto_factor1', '_auto_factor2'} <= set(output.obs.columns)
