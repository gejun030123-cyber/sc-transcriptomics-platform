"""Focused tests for bounded before/after batch-integration evaluation."""

from types import SimpleNamespace
import sys

import numpy as np
import pandas as pd


def _analysis(params=None):
    from modules.batch_correct import BatchCorrectAnalysis

    return BatchCorrectAnalysis(
        project_dir='/tmp/test-batch-evaluation',
        params=params or {},
        progress_callback=lambda *_: None,
    )


def test_evaluation_sample_is_bounded_and_keeps_rare_batches():
    labels = pd.Series(['a'] * 8 + ['b'] + ['c'])

    indices, metadata = _analysis()._evaluation_sample_indices(labels, requested_size=5)

    assert len(indices) == 5
    assert set(labels.iloc[indices]) == {'a', 'b', 'c'}
    assert metadata['strategy'] == 'batch_stratified'
    assert metadata['n_batches_in_sample'] == 3
    assert metadata['used_size'] == len(indices)


def test_shared_pre_post_evaluation_emits_comparable_deltas():
    obs = pd.DataFrame({
        'batch': ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B'],
        'celltype': ['T', 'B', 'T', 'B', 'T', 'B', 'T', 'B'],
    })
    adata = SimpleNamespace(obs=obs, obsm={})
    analysis = _analysis({'evaluation_sample_size': 8, 'bio_label_key': 'celltype'})
    indices, sampling = analysis._evaluation_sample_indices(obs['batch'], requested_size=8)

    before_representation = np.array([
        [-10, 0], [-9, 0], [-8, 0], [-7, 0],
        [7, 0], [8, 0], [9, 0], [10, 0],
    ], dtype=float)
    # After correction, cells are organized by biological label while batches
    # are interleaved within each label.
    after_representation = np.array([
        [0, 0], [10, 0], [0, .1], [10, .1],
        [.1, 0], [10.1, 0], [.1, .1], [10.1, .1],
    ], dtype=float)

    before = analysis._compute_evaluation(
        adata, before_representation, 'batch', indices, sampling,
    )
    after = analysis._compute_evaluation(
        adata, after_representation, 'batch', indices, sampling,
    )
    comparison = analysis._build_evaluation_comparison(before, after, sampling, 'harmony')

    assert before['n_evaluation_cells'] == after['n_evaluation_cells'] == 8
    assert comparison['sampling']['used_size'] == 8
    assert comparison['delta']['abs_asw_batch'] < 0
    assert comparison['delta']['asw_bio'] > 0
    assert {row['metric'] for row in analysis._comparison_rows(comparison)} >= {
        'abs_asw_batch', 'asw_bio',
    }


def test_tiny_evaluation_records_asw_warning_instead_of_raising():
    obs = pd.DataFrame({'batch': ['A', 'B']})
    adata = SimpleNamespace(obs=obs, obsm={})
    analysis = _analysis({'evaluation_sample_size': 100})
    indices, sampling = analysis._evaluation_sample_indices(obs['batch'], requested_size=100)

    metrics = analysis._compute_evaluation(
        adata, np.array([[0.0, 0.0], [1.0, 1.0]]), 'batch', indices, sampling,
    )

    assert metrics['n_evaluation_cells'] == 2
    assert 'asw_batch_warning' in metrics
    assert 'neighbor_mixing_warning' not in metrics


def test_pre_post_comparison_writes_downloadable_delta_table(tmp_path):
    analysis = _analysis()
    analysis.project_dir = str(tmp_path)
    comparison = {
        'before': {
            'abs_asw_batch': 0.7,
            'mean_neighbor_batch_entropy': 0.2,
            'mean_neighbor_same_batch_fraction': 0.8,
        },
        'after': {
            'abs_asw_batch': 0.3,
            'mean_neighbor_batch_entropy': 0.6,
            'mean_neighbor_same_batch_fraction': 0.4,
        },
        'delta': {
            'abs_asw_batch': -0.4,
            'mean_neighbor_batch_entropy': 0.4,
            'mean_neighbor_same_batch_fraction': -0.4,
        },
        'metric_directions': {
            'abs_asw_batch': 'lower',
            'mean_neighbor_batch_entropy': 'higher',
            'mean_neighbor_same_batch_fraction': 'lower',
        },
        'sampling': {'used_size': 40, 'strategy': 'batch_stratified'},
        'comparison_scope': 'embedding_and_graph',
    }
    metrics = dict(comparison['after'])
    adata = SimpleNamespace(obs=pd.DataFrame({'batch': ['A', 'B']}))

    result_files = analysis._make_eval_plots(
        adata, metrics, 'batch', 'harmony', analysis.ensure_plots_dir(), comparison,
    )

    csv_files = [item for item in result_files if item['file_type'] == 'csv']
    assert len(csv_files) == 1
    table = pd.read_csv(csv_files[0]['file_path'])
    assert set(table['metric']) == {
        'abs_asw_batch', 'mean_neighbor_batch_entropy',
        'mean_neighbor_same_batch_fraction',
    }


def test_scvi_records_the_actual_training_layer(monkeypatch):
    """The scVI completion path must not reference a missing ``layer`` name."""
    setup_calls = []

    class FakeSCVI:
        @classmethod
        def setup_anndata(cls, adata, layer, batch_key):
            setup_calls.append((adata, layer, batch_key))

        def __init__(self, _adata, **_kwargs):
            pass

        def train(self, **_kwargs):
            return None

        def get_latent_representation(self):
            return np.ones((3, 2), dtype=np.float32)

    fake_scvi = SimpleNamespace(model=SimpleNamespace(SCVI=FakeSCVI))
    monkeypatch.setitem(sys.modules, 'scvi', fake_scvi)
    work = SimpleNamespace(
        layers={'counts': np.array([[1, 0], [0, 2], [3, 1]])},
        n_obs=3, n_vars=2, obs_names=['a', 'b', 'c'], obsm={},
    )
    analysis = _analysis({'max_epochs': 1})
    monkeypatch.setattr(analysis, '_hvg_adata', lambda _adata: work)

    representation, info = analysis._run_scvi(work, 'batch')

    assert representation == 'X_scVI'
    assert setup_calls[0][1:] == ('counts', 'batch')
    assert info['layer'] == 'counts'
    assert work.obsm['X_scVI'].shape == (3, 2)
