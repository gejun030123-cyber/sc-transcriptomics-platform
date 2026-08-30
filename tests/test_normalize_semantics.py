from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def _counts_adata(*, transformed_x=False, with_counts=True, n_obs=12, n_vars=18):
    rng = np.random.default_rng(42)
    counts = rng.poisson(4, size=(n_obs, n_vars)).astype(np.float32)
    counts[:, 0] += 1
    if transformed_x:
        row_sums = counts.sum(axis=1, keepdims=True)
        x = np.log1p(counts * (10000.0 / row_sums))
    else:
        x = counts.copy()
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f'cell_{i}' for i in range(n_obs)]),
        var=pd.DataFrame(index=[f'gene_{i}' for i in range(n_vars)]),
    )
    if with_counts:
        adata.layers['counts'] = counts.copy()
    return adata, counts


def _run_normalize(tmp_path, adata, params):
    from modules.normalize import NormalizeAnalysis

    input_path = tmp_path / 'input.h5ad'
    adata.write_h5ad(input_path)
    captured_filenames = []
    module = NormalizeAnalysis(
        project_dir=str(tmp_path),
        params=params,
        progress_callback=lambda *_: None,
    )

    def _capture_figure(_fig, _plots_dir, filename, *_args, **_kwargs):
        captured_filenames.append(filename)
        return []

    module.save_matplotlib_figure = _capture_figure
    result = module.run(str(input_path))
    return ad.read_h5ad(result['output_adata']), result, captured_filenames


def test_log1p_preserves_existing_counts_and_normalizes_from_them(tmp_path):
    adata, counts = _counts_adata(transformed_x=True, with_counts=True)

    out, result, filenames = _run_normalize(
        tmp_path, adata, {'method': 'log1p', 'target_sum': 10000},
    )

    expected = np.log1p(counts * (10000.0 / counts.sum(axis=1, keepdims=True)))
    np.testing.assert_array_equal(out.layers['counts'], counts)
    np.testing.assert_allclose(out.X, expected, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(
        out.layers['normalized'],
        counts * (10000.0 / counts.sum(axis=1, keepdims=True)),
        rtol=1e-6, atol=1e-6,
    )
    assert result['summary']['counts_source'] == 'layers[counts]'
    assert result['summary']['normalized_layer_matches_x'] is False
    assert result['summary']['layer_semantics']['X'].startswith('log1p')
    assert 'normalize_libsize.png' in filenames
    assert 'highly_variable' not in out.var.columns
    assert out.shape == adata.shape


def test_missing_counts_rejects_already_normalized_x(tmp_path):
    from modules.normalize import NormalizeAnalysis

    adata, _ = _counts_adata(transformed_x=True, with_counts=False)
    input_path = tmp_path / 'normalized_without_counts.h5ad'
    adata.write_h5ad(input_path)
    module = NormalizeAnalysis(
        project_dir=str(tmp_path), params={'method': 'log1p'},
        progress_callback=lambda *_: None,
    )

    with pytest.raises(ValueError, match=r"未找到 layers\['counts'\]"):
        module.run(str(input_path))


def test_raw_x_is_promoted_to_counts_only_when_count_like(tmp_path):
    adata, counts = _counts_adata(transformed_x=False, with_counts=False)

    out, result, _ = _run_normalize(tmp_path, adata, {'method': 'log1p'})

    np.testing.assert_array_equal(out.layers['counts'], counts)
    assert result['summary']['counts_source'] == 'adata.X'


def test_pearson_keeps_x_and_normalized_in_sync_without_hvg_side_effects(tmp_path):
    adata, counts = _counts_adata(transformed_x=True, with_counts=True)
    counts[:, -1] = 0
    adata.layers['counts'][:, -1] = 0

    out, result, filenames = _run_normalize(
        tmp_path, adata,
        {'method': 'pearson_residuals', 'clip_values': True},
    )

    np.testing.assert_array_equal(out.layers['counts'], counts)
    np.testing.assert_allclose(out.layers['normalized'], out.X, rtol=0, atol=0)
    assert np.isfinite(out.X).all()
    assert np.any(out.X < 0)
    assert not np.allclose(out.X, counts)
    assert np.max(np.abs(out.X)) <= np.sqrt(out.n_obs) + 1e-8
    np.testing.assert_array_equal(out.X[:, -1], np.zeros(out.n_obs))
    assert out.shape == adata.shape
    assert list(out.obs_names) == list(adata.obs_names)
    assert list(out.var_names) == list(adata.var_names)
    assert 'highly_variable' not in out.var.columns
    assert 'pearson_residuals_normalization' in out.uns
    assert out.uns['normalization']['x_contains'] == 'pearson_residuals'
    assert result['summary']['target_sum'] is None
    assert result['summary']['library_size_plot'] is False
    assert result['summary']['zero_count_genes_set_to_zero'] == 1
    assert 'normalize_libsize.png' not in filenames
    assert filenames == ['normalize_expression_distribution.png']


def test_missing_pearson_api_raises_instead_of_copying_counts():
    from modules.normalize import _require_pearson_normalizer

    fake_scanpy = SimpleNamespace(
        experimental=SimpleNamespace(pp=SimpleNamespace()),
    )

    with pytest.raises(RuntimeError, match='不能把原始 counts 静默写成 normalized'):
        _require_pearson_normalizer(fake_scanpy)


def test_pearson_clip_checkbox_uses_scanpy_semantics():
    from modules.normalize import _pearson_clip_setting

    assert _pearson_clip_setting(True) == (True, None)
    assert _pearson_clip_setting('true') == (True, None)
    enabled, clip = _pearson_clip_setting(False)
    assert enabled is False and np.isinf(clip)
    enabled, clip = _pearson_clip_setting('false')
    assert enabled is False and np.isinf(clip)


def test_expression_distribution_uses_only_the_selected_output_scale():
    from modules.normalize import _expression_distribution_payload

    log_values, log_label, _ = _expression_distribution_payload(
        'log1p', np.asarray([[0.0, 1.2], [2.3, 0.0]]),
    )
    pearson_values, pearson_label, _ = _expression_distribution_payload(
        'pearson_residuals', np.asarray([[-1.0, 0.0], [0.5, 1.2]]),
    )

    np.testing.assert_allclose(log_values, [1.2, 2.3])
    np.testing.assert_allclose(pearson_values, [-1.0, 0.0, 0.5, 1.2])
    assert 'raw' not in log_label.lower()
    assert pearson_label == 'Pearson residuals'


def test_sparse_total_normalization_preserves_sparse_shape_and_zero_cells():
    from modules.normalize import _normalize_total_matrix

    counts = sparse.csr_matrix(
        np.asarray([[1, 2, 0], [0, 0, 0], [4, 0, 1]], dtype=np.float32)
    )

    normalized = _normalize_total_matrix(counts, 10000)

    assert sparse.isspmatrix_csr(normalized)
    assert normalized.shape == counts.shape
    np.testing.assert_allclose(
        np.asarray(normalized.sum(axis=1)).ravel(), [10000, 0, 10000],
        rtol=1e-6,
    )
