"""Regression tests for safe obs grouping semantics across SC/Bulk modules."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def _adata(n_obs=24, n_vars=12):
    anndata = pytest.importorskip('anndata')
    rng = np.random.default_rng(11)
    obs = pd.DataFrame({
        'leiden': pd.Categorical(['0'] * (n_obs // 2) + ['1'] * (n_obs - n_obs // 2)),
        'batch': pd.Categorical(['A'] * (n_obs // 2) + ['B'] * (n_obs - n_obs // 2)),
        'continuous_qc': np.arange(n_obs, dtype=float),
    })
    return anndata.AnnData(
        rng.poisson(2.0, size=(n_obs, n_vars)).astype(float),
        obs=obs,
        var=pd.DataFrame(index=[f'G{i}' for i in range(n_vars)]),
    )


def test_obs_grouping_rejects_continuous_and_accepts_categories():
    from modules.io_utils import obs_grouping_info

    adata = _adata()
    bad = obs_grouping_info(adata, 'continuous_qc', require_multiple=True)
    assert bad['valid'] is False
    assert '连续' in bad['reason'] or '高基数' in bad['reason']

    adata.obs['cell_like_id'] = [f'cell_{i}' for i in range(adata.n_obs)]
    id_info = obs_grouping_info(adata, 'cell_like_id', require_multiple=True)
    assert id_info['valid'] is False

    good = obs_grouping_info(adata, 'leiden', require_multiple=True)
    assert good['valid'] is True
    assert good['n_unique'] == 2


def test_obs_grouping_accepts_many_sample_ids_in_large_study():
    """A 131-sample study must keep sample_id usable for batch/sample roles.

    The small categorical cap exists to keep plotting groupings readable; it
    must not disqualify a real library/donor identifier just because the study
    contains more than 50 samples.
    """
    from modules.io_utils import obs_grouping_info, rank_obs_grouping_candidates

    n_obs, n_samples = 2600, 131
    adata = _adata(n_obs=n_obs, n_vars=6)
    adata.obs['sample_id'] = pd.Categorical(
        [f'S{index % n_samples:03d}' for index in range(n_obs)]
    )
    adata.obs['Sample'] = pd.Categorical(
        [f'CAP{index % n_samples:03d}' for index in range(n_obs)]
    )
    adata.obs['donor'] = pd.Categorical(
        [f'D{index % 40:02d}' for index in range(n_obs)]
    )

    info = obs_grouping_info(adata, 'sample_id', max_categories=50, require_multiple=True)
    assert info['valid'] is True
    assert info['n_unique'] == n_samples
    assert info['identifier_column'] is True
    assert obs_grouping_info(
        adata, 'Sample', max_categories=50, require_multiple=True,
    )['valid'] is True

    # Sample identifiers stay selectable for batch/sample parameters, but must
    # not be auto-suggested as the biological groupby.
    assert 'sample_id' in rank_obs_grouping_candidates(adata, 'dimred', purpose='batch_key')
    assert 'sample_id' in rank_obs_grouping_candidates(adata, 'dimred', purpose='sample_key')
    assert 'sample_id' not in rank_obs_grouping_candidates(adata, 'deg', purpose='groupby')
    assert 'donor' not in rank_obs_grouping_candidates(adata, 'deg', purpose='groupby')


def test_obs_grouping_still_rejects_one_id_per_cell_and_huge_identifiers():
    from modules.io_utils import obs_grouping_info

    n_obs = 400
    adata = _adata(n_obs=n_obs, n_vars=4)
    # One "sample" per cell is a barcode, not a replicate design.
    adata.obs['sample_id'] = [f'S{index}' for index in range(n_obs)]
    one_per_cell = obs_grouping_info(adata, 'sample_id', require_multiple=True)
    assert one_per_cell['valid'] is False
    assert one_per_cell['identifier_column'] is True

    # Identifier fields are documented but not unbounded.
    large = _adata(n_obs=3000, n_vars=4)
    large.obs['donor_id'] = [f'D{index}' for index in range(large.n_obs)]
    assert obs_grouping_info(large, 'donor_id', require_multiple=True)['valid'] is False


def test_obs_grouping_accepts_numeric_sample_ids():
    """Numeric 1..131 IDs are sample identifiers, not a continuous QC gradient."""
    from modules.io_utils import obs_grouping_info

    n_obs, n_samples = 2620, 131
    adata = _adata(n_obs=n_obs, n_vars=4)
    adata.obs['sample_id'] = [index % n_samples + 1 for index in range(n_obs)]

    info = obs_grouping_info(adata, 'sample_id', max_categories=50, require_multiple=True)
    assert info['valid'] is True
    assert info['identifier_column'] is True
    assert info['n_unique'] == n_samples

    # A continuous numeric field is still rejected.
    assert obs_grouping_info(adata, 'continuous_qc', require_multiple=True)['valid'] is False


def test_per_sample_figure_status_suppresses_only_large_sample_columns():
    from modules.io_utils import (
        MAX_PER_SAMPLE_FIGURE_LEVELS,
        per_sample_figure_status,
    )

    n_obs = 2620
    adata = _adata(n_obs=n_obs, n_vars=4)
    adata.obs['sample_id'] = [f'S{index % 131:03d}' for index in range(n_obs)]
    adata.obs['cluster'] = [f'K{index % 131:03d}' for index in range(n_obs)]
    adata.obs['condition'] = [f'C{index % 40:02d}' for index in range(n_obs)]

    assert per_sample_figure_status(adata, 'sample_id') == (True, 131)
    # An ordinary grouping with the same cardinality is not a sample figure.
    assert per_sample_figure_status(adata, 'cluster') == (False, 131)
    assert per_sample_figure_status(adata, 'condition') == (False, 40)
    assert per_sample_figure_status(adata, 'not_a_column') == (False, 0)
    assert MAX_PER_SAMPLE_FIGURE_LEVELS == 20


def test_per_sample_figure_status_keeps_small_studies():
    from modules.io_utils import per_sample_figure_status

    n_obs = 40
    adata = _adata(n_obs=n_obs, n_vars=4)
    adata.obs['sample_id'] = [f'S{index % 5}' for index in range(n_obs)]

    assert per_sample_figure_status(adata, 'sample_id') == (False, 5)
    # The caller may raise the limit for an explicitly requested wide figure.
    assert per_sample_figure_status(adata, 'sample_id', max_levels=3) == (True, 5)


def test_obs_grouping_candidates_are_parameter_aware():
    from modules.io_utils import rank_obs_grouping_candidates

    adata = _adata()
    adata.obs['barcode'] = [f'cell_{i}' for i in range(adata.n_obs)]
    adata.obs['sample'] = pd.Categorical(
        [f'S{i // 6 + 1}' for i in range(adata.n_obs)]
    )
    adata.obs['condition'] = pd.Categorical(
        ['Control'] * (adata.n_obs // 2) + ['Treatment'] * (adata.n_obs - adata.n_obs // 2)
    )

    assert rank_obs_grouping_candidates(adata, 'deg', purpose='groupby')[0] == 'condition'
    assert rank_obs_grouping_candidates(adata, 'deg', purpose='batch_key')[0] == 'batch'
    assert rank_obs_grouping_candidates(adata, 'deg', purpose='cluster_key')[0] == 'leiden'
    assert rank_obs_grouping_candidates(adata, 'deg', purpose='sample_key')[0] == 'sample'


def test_batch_correct_rejects_continuous_batch_before_method_run(tmp_path):
    from modules.batch_correct import BatchCorrectAnalysis

    adata = _adata()
    adata.obsm['X_pca'] = np.zeros((adata.n_obs, 3), dtype=float)
    analysis = BatchCorrectAnalysis(
        str(tmp_path), {'batch_key': 'continuous_qc'}, lambda *_: None,
    )
    message = analysis.validate_input(adata)
    assert message is not None
    assert 'batch_key' in message
    assert '分类' in message or '连续' in message


def test_proportion_does_not_build_crosstab_from_continuous_batch(tmp_path):
    from modules.proportion import ProportionAnalysis

    adata = _adata()
    input_path = tmp_path / 'input.h5ad'
    adata.write_h5ad(input_path)
    analysis = ProportionAnalysis(
        str(tmp_path),
        {'groupby': 'leiden', 'batch_key': 'continuous_qc'},
        lambda *_: None,
    )
    with pytest.raises(ValueError, match='batch_key'):
        analysis.run(str(input_path))


def test_bulk_deg_missing_group_never_splits_samples_by_file_order(tmp_path):
    from modules.bulk_deg import BulkDEGAnalysis

    # A tiny expression matrix with no ``condition`` metadata.  The module must
    # ask for a real group column instead of inventing Group1/Group2 by position.
    matrix = pd.DataFrame(
        np.arange(30, dtype=float).reshape(5, 6),
        index=[f'G{i}' for i in range(5)],
        columns=[f'S{i}' for i in range(6)],
    )
    input_path = tmp_path / 'bulk.tsv'
    matrix.to_csv(input_path, sep='\t')
    analysis = BulkDEGAnalysis(
        str(tmp_path), {'method': 't-test', 'groupby': 'condition'}, lambda *_: None,
    )
    with pytest.raises(ValueError, match='分组列'):
        analysis.run(str(input_path))


def test_bulk_deg_rejects_comparisons_from_a_different_grouping_column(tmp_path):
    from modules.bulk_deg import BulkDEGAnalysis

    adata = _adata(n_obs=6, n_vars=5)
    adata.obs['condition'] = ['Ctrl_B'] * 3 + ['Drug_B'] * 3
    adata.obs['treatment'] = ['Ctrl'] * 3 + ['Drug'] * 3
    input_path = tmp_path / 'input.h5ad'
    adata.write_h5ad(input_path)

    analysis = BulkDEGAnalysis(
        str(tmp_path),
        {
            'method': 't-test', 'groupby': 'treatment',
            'comparisons': 'Drug_B-vs-Ctrl_B',
        },
        lambda *_: None,
    )
    with pytest.raises(ValueError, match='指定比较与分组列'):
        analysis.run(str(input_path))


def test_plotly_umap_keeps_continuous_qc_as_one_colorscale_trace():
    from modules.visualization import umap_scatter

    adata = _adata()
    adata.obsm['X_umap'] = np.column_stack([
        np.arange(adata.n_obs, dtype=float),
        np.zeros(adata.n_obs, dtype=float),
    ])
    figure = umap_scatter(adata, 'continuous_qc')
    assert len(figure['data']) == 1
    assert figure['data'][0]['marker'].get('colorscale') is not None


def test_bulk_pca_uses_sample_name_groups_when_obs_group_is_missing_or_constant():
    from modules.bulk_pca import _resolve_pca_color_groups

    adata = _adata(n_obs=6, n_vars=5)
    adata.obs.index = ['Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'Drug_1', 'Drug_2', 'Drug_3']
    adata.obs['condition'] = 'all_samples'

    groups, info = _resolve_pca_color_groups(adata, 'condition')

    assert groups == ['Ctrl', 'Ctrl', 'Ctrl', 'Drug', 'Drug', 'Drug']
    assert info['used'] == '_auto_group'
    assert info['source'] == 'sample_name_inference'
    assert info['group_counts'] == {'Ctrl': 3, 'Drug': 3}
    assert '_auto_group' in adata.obs


def test_bulk_pca_materializes_and_uses_second_auto_factor_from_ui_text():
    from modules.bulk_pca import _resolve_pca_groups

    adata = _adata(n_obs=8, n_vars=5)
    adata.obs.index = [
        'Ctr_B_1', 'Ctr_B_2', 'Ctr_En_1', 'Ctr_En_2',
        'PEA_B_1', 'PEA_B_2', 'PEA_En_1', 'PEA_En_2',
    ]

    colors, color_info = _resolve_pca_groups(adata, '_auto_group_', role='color')
    markers, marker_info = _resolve_pca_groups(
        adata, '第2因素：B, En', role='marker',
    )

    assert color_info['used'] == '_auto_group'
    assert marker_info['used'] == '_auto_factor2'
    assert marker_info['source'] == 'sample_name_factor_inference'
    assert set(colors) == {'Ctr_B', 'Ctr_En', 'PEA_B', 'PEA_En'}
    assert set(markers) == {'B', 'En'}
    assert {'_auto_group', '_auto_factor1', '_auto_factor2'} <= set(adata.obs.columns)


def test_sample_name_groups_support_attached_replicate_suffixes():
    from modules.io_utils import infer_sample_group_candidates

    sample_names = [
        'EC_CT_BSA1', 'EC_CT_BSA2', 'EC_CT_BSA3',
        'Hep_KO_FFA1', 'Hep_KO_FFA2', 'Hep_KO_FFA3',
    ]
    candidates = infer_sample_group_candidates(sample_names)

    assert candidates[0]['key'] == 'combined'
    assert candidates[0]['mapping']['EC_CT_BSA1'] == 'EC_CT_BSA'
    assert candidates[0]['mapping']['Hep_KO_FFA3'] == 'Hep_KO_FFA'
