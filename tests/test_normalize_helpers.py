import numpy as np
import pytest


def test_tmm_equal_samples():
    """完全相同的样本 → 因子全为 1"""
    from modules.bulk_normalize import _tmm_normalize
    counts = np.array([
        [10, 20, 30, 40, 50],
        [10, 20, 30, 40, 50],
        [10, 20, 30, 40, 50],
    ], dtype=float)
    factors = _tmm_normalize(counts)
    assert len(factors) == 3
    np.testing.assert_allclose(factors, 1.0, atol=0.1)


def test_tmm_scaled_sample():
    """组成差异大的样本 → TMM 因子有明显差异"""
    from modules.bulk_normalize import _tmm_normalize
    np.random.seed(42)
    # Sample 0: uniform expression, Sample 1: strong composition shift (half genes up, half down)
    base = np.random.randint(100, 500, size=200).astype(float)
    shifted = base.copy()
    shifted[:100] *= 5   # half genes strongly up
    shifted[100:] /= 5   # half genes strongly down (keep > 0)
    shifted = np.maximum(shifted, 1)
    counts = np.vstack([base, shifted])
    factors = _tmm_normalize(counts)
    # Composition shift should produce factors that differ
    assert abs(factors[1] - factors[0]) > 0.05


def test_tmm_no_crash_zeros():
    """含零值不崩溃"""
    from modules.bulk_normalize import _tmm_normalize
    counts = np.array([
        [0, 10, 20, 0, 50],
        [5, 0, 30, 40, 0],
        [10, 20, 0, 0, 100],
    ], dtype=float)
    factors = _tmm_normalize(counts)
    assert len(factors) == 3
    assert all(f > 0 for f in factors)


def test_vst_shape():
    """VST 输出形状与输入一致"""
    from modules.bulk_normalize import _vst_transform
    counts = np.array([[100, 200, 300], [150, 250, 350]], dtype=float)
    sf = np.array([1.0, 1.0])
    result = _vst_transform(counts, sf)
    assert result.shape == counts.shape


def test_vst_variance_stabilized():
    """VST 后方差应比原始 log2 更稳定"""
    from modules.bulk_normalize import _vst_transform
    np.random.seed(42)
    counts = np.array([
        [10, 100, 1000, 10000],
        [12, 120, 1200, 12000],
        [8, 80, 800, 8000],
    ], dtype=float)
    sf = np.array([1.0, 1.0, 1.0])
    vst = _vst_transform(counts, sf)
    variances = vst.var(axis=0)
    assert variances.max() / (variances.min() + 1e-10) < 10


def test_rlog_shape():
    """rlog 输出形状与输入一致"""
    from modules.bulk_normalize import _rlog_transform
    counts = np.array([[100, 200, 300], [150, 250, 350]], dtype=float)
    sf = np.array([1.0, 1.0])
    result = _rlog_transform(counts, sf)
    assert result.shape == counts.shape


def test_rlog_small_sample_regularity():
    """小样本时 rlog 收缩效应"""
    from modules.bulk_normalize import _rlog_transform
    counts = np.array([
        [10, 100, 1000],
        [50, 50, 500],
        [5, 200, 2000],
    ], dtype=float)
    sf = np.array([1.0, 1.0, 1.0])
    rlog = _rlog_transform(counts, sf)
    direct_log = np.log2(counts + 1)
    assert rlog.var(axis=0).mean() <= direct_log.var(axis=0).mean() * 1.5


def test_tmm_large_sample():
    """较大样本集不崩溃"""
    from modules.bulk_normalize import _tmm_normalize
    np.random.seed(42)
    counts = np.random.randint(5, 5000, size=(12, 100)).astype(float)
    factors = _tmm_normalize(counts)
    assert len(factors) == 12
    assert all(f > 0 for f in factors)


def test_bulk_normalize_schema_supports_log2_for_continuous_expression():
    from modules.schemas import PARAM_SCHEMAS

    method = next(p for p in PARAM_SCHEMAS['bulk_normalize'] if p['key'] == 'method')
    assert 'log2' in method['options']


def test_measurement_type_uses_filename_hint_for_fpkm():
    from modules.bulk_normalize import _infer_measurement_type
    assert _infer_measurement_type(np.array([[1., 2.], [3., 4.]]), 'sample.fpkm.tsv') == 'continuous_expression'


def test_measurement_type_accepts_integer_counts():
    from modules.bulk_normalize import _infer_measurement_type
    assert _infer_measurement_type(np.array([[1, 2], [3, 4]]), 'counts.tsv') == 'raw_counts'
