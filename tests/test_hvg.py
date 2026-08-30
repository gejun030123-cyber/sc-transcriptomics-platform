"""Regression tests for HVG input scale, selection semantics, and diagnostics."""

import numpy as np
import pandas as pd
import pytest


def _expression_adata(n_obs=24, n_vars=40):
    anndata = pytest.importorskip('anndata')
    rng = np.random.default_rng(17)
    counts = rng.poisson(3.0, size=(n_obs, n_vars)).astype(np.float32)
    counts[:, 0] += np.arange(n_obs) % 5
    totals = counts.sum(axis=1, keepdims=True)
    normalized = np.log1p(counts / np.maximum(totals, 1) * 1e4).astype(np.float32)
    obs = pd.DataFrame({
        'batch': pd.Categorical(['A'] * (n_obs // 2) + ['B'] * (n_obs - n_obs // 2)),
    })
    var_names = ['MT-ND1', *[f'Gene{i}' for i in range(1, n_vars)]]
    adata = anndata.AnnData(
        X=normalized,
        obs=obs,
        var=pd.DataFrame(index=var_names),
    )
    adata.layers['counts'] = counts
    adata.uns['log1p'] = {'base': None}
    return adata


def _install_fake_hvg(monkeypatch, calls):
    import scanpy as sc

    def fake_hvg(adata, **kwargs):
        calls.append(dict(kwargs))
        n_top = int(kwargs['n_top_genes'])
        n_vars = adata.n_vars
        selected = np.zeros(n_vars, dtype=bool)
        selected[:n_top] = True
        adata.var['highly_variable'] = selected
        adata.var['means'] = np.linspace(0.01, 2.0, n_vars)
        if kwargs['flavor'] == 'seurat_v3':
            adata.var['variances'] = np.linspace(2.0, 0.1, n_vars)
            adata.var['variances_norm'] = np.linspace(3.0, 0.5, n_vars)
            rank = np.full(n_vars, np.nan, dtype=float)
            rank[:n_top] = np.arange(n_top, dtype=float)
            adata.var['highly_variable_rank'] = rank
        else:
            adata.var['dispersions'] = np.linspace(2.0, 0.1, n_vars)
            adata.var['dispersions_norm'] = np.linspace(3.0, 0.5, n_vars)
        if kwargs.get('batch_key'):
            adata.var['highly_variable_nbatches'] = np.where(selected, 2, 1)
            adata.var['highly_variable_intersection'] = selected
        adata.uns['hvg'] = {'flavor': kwargs['flavor']}

    monkeypatch.setattr(sc.pp, 'highly_variable_genes', fake_hvg)


@pytest.mark.parametrize(
    ('flavor', 'expected_layer', 'expected_input'),
    [
        ('seurat_v3', 'counts', 'counts'),
        ('seurat', None, 'X'),
        ('cell_ranger', None, 'X'),
    ],
)
def test_hvg_flavors_use_the_correct_expression_source(
    flavor, expected_layer, expected_input,
):
    from modules.hvg import _resolve_hvg_input

    layer, source = _resolve_hvg_input(_expression_adata(), flavor)
    assert layer == expected_layer
    assert source == expected_input


def test_dispersion_flavor_rejects_raw_counts_x():
    from modules.hvg import _resolve_hvg_input

    adata = _expression_adata()
    adata.X = adata.layers['counts'].copy()
    with pytest.raises(ValueError, match='log-normalized'):
        _resolve_hvg_input(adata, 'seurat')


def test_seurat_v3_rejects_noninteger_counts_layer():
    from modules.hvg import _resolve_hvg_input

    adata = _expression_adata()
    adata.layers['counts'] = adata.X.copy()
    with pytest.raises(ValueError, match='原始非负整数 counts'):
        _resolve_hvg_input(adata, 'seurat_v3')


def test_previous_hvg_columns_are_cleared_before_rerun():
    from modules.hvg import _clear_previous_hvg_results

    var = pd.DataFrame({
        'gene_ids': ['id1', 'id2'],
        'highly_variable': [True, False],
        'variances_norm': [2.0, 1.0],
        'dispersions_norm': [9.0, 8.0],
        'highly_variable_rank_scanpy': [0.0, np.nan],
    })
    _clear_previous_hvg_results(var)
    assert list(var.columns) == ['gene_ids']


def test_force_include_and_exclusion_preserve_requested_hvg_count():
    from modules.hvg import _finalize_hvg_selection

    n_vars = 8
    var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_vars)])
    var['highly_variable'] = [True, True, True, False, False, False, False, False]
    var['variances_norm'] = np.arange(n_vars, 0, -1, dtype=float)
    var['highly_variable_rank'] = [0.0, 1.0, 2.0, np.nan, np.nan, np.nan, np.nan, np.nan]
    excluded = np.array([True, False, False, False, False, False, False, False])
    forced = np.array([False, False, False, False, False, True, False, False])

    _finalize_hvg_selection(var, 3, 'seurat_v3', excluded, forced)

    assert int(var['highly_variable'].sum()) == 3
    assert not bool(var.loc['Gene0', 'highly_variable'])
    assert bool(var.loc['Gene5', 'highly_variable'])
    assert var['highly_variable_rank'].notna().equals(var['highly_variable'])
    assert 'highly_variable_rank_scanpy' in var.columns


def test_too_many_forced_genes_fail_instead_of_exceeding_target():
    from modules.hvg import _finalize_hvg_selection

    var = pd.DataFrame({
        'highly_variable': [True, True, False, False],
        'variances_norm': [4.0, 3.0, 2.0, 1.0],
        'highly_variable_rank': [0.0, 1.0, np.nan, np.nan],
    })
    with pytest.raises(ValueError, match='超过目标 HVG 数'):
        _finalize_hvg_selection(
            var, 2, 'seurat_v3', np.zeros(4, dtype=bool),
            np.array([True, True, True, False]),
        )


def test_regress_cc_fails_loudly_without_scores(tmp_path):
    from modules.hvg import HVGAnalysis

    analysis = HVGAnalysis(
        str(tmp_path),
        {'hvg_flavor': 'seurat_v3', 'regress_cc': True, 'cc_scoring': False},
        lambda *_: None,
    )
    message = analysis.validate_input(_expression_adata())
    assert message is not None
    assert 'S_score/G2M_score' in message


def test_regress_cc_rejects_raw_counts_even_when_scores_exist(tmp_path):
    from modules.hvg import HVGAnalysis

    adata = _expression_adata()
    adata.X = adata.layers['counts'].copy()
    adata.obs['S_score'] = 0.0
    adata.obs['G2M_score'] = 0.0
    analysis = HVGAnalysis(
        str(tmp_path),
        {'hvg_flavor': 'seurat_v3', 'regress_cc': True},
        lambda *_: None,
    )
    message = analysis.validate_input(adata)
    assert message is not None
    assert '不能作用于原始 counts' in message


def test_invalid_requested_batch_is_not_silently_ignored(tmp_path):
    from modules.hvg import HVGAnalysis

    analysis = HVGAnalysis(
        str(tmp_path),
        {'hvg_flavor': 'seurat_v3', 'batch_key': 'missing_batch'},
        lambda *_: None,
    )
    message = analysis.validate_input(_expression_adata())
    assert message is not None
    assert 'batch_key' in message
    assert 'missing_batch' in message


@pytest.mark.parametrize(
    ('flavor', 'expected_layer'),
    [('seurat_v3', 'counts'), ('seurat', None), ('cell_ranger', None)],
)
def test_run_calls_scanpy_once_and_plots_all_genes(
    tmp_path, monkeypatch, flavor, expected_layer,
):
    from modules.hvg import HVGAnalysis
    import modules.native_figures as native_figures

    adata = _expression_adata()
    calls = []
    plotted_lengths = []
    saved = {}
    _install_fake_hvg(monkeypatch, calls)

    def fake_scatter(x, y, **kwargs):
        plotted_lengths.append((len(x), len(y), len(kwargs.get('labels') or [])))
        return object()

    monkeypatch.setattr(native_figures, 'scatter_figure', fake_scatter)
    analysis = HVGAnalysis(
        str(tmp_path),
        {
            'n_top_genes': 10,
            'hvg_flavor': flavor,
            'batch_key': 'batch',
            'show_hvg_rank_plot': True,
        },
        lambda *_: None,
    )
    analysis.load_adata = lambda _: adata.copy()
    analysis.ensure_plots_dir = lambda: str(tmp_path)
    analysis.save_matplotlib_figure = lambda *args, **kwargs: []

    def fake_save_output(output, _module_name):
        saved['adata'] = output
        return str(tmp_path / 'hvg_output.h5ad')

    analysis.save_output = fake_save_output
    result = analysis.run('unused.h5ad')

    assert len(calls) == 1
    assert calls[0].get('layer') == expected_layer
    assert calls[0]['batch_key'] == 'batch'
    assert result['summary']['n_hvgs'] == 10
    assert int(saved['adata'].var['highly_variable'].sum()) == 10
    assert [entry[:2] for entry in plotted_lengths] == [
        (adata.n_vars, adata.n_vars),
        (adata.n_vars, adata.n_vars),
    ]
    assert plotted_lengths[1][2] == adata.n_vars


def test_obsolete_manual_batch_strategy_is_not_exposed():
    from modules.schemas import PARAM_SCHEMAS

    keys = {field['key'] for field in PARAM_SCHEMAS['hvg']}
    assert 'batch_hvg_strategy' not in keys
