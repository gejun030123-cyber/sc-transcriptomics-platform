# tests/test_bulk_qc_helpers.py
import numpy as np
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


def test_detect_outliers_too_few():
    """样本数 < 4 → 返回空"""
    from modules.bulk_qc import _detect_outliers_mahal
    coords = np.array([[0, 0], [1, 1], [2, 2]])
    outliers = _detect_outliers_mahal(coords, ['a', 'b', 'c'])
    assert outliers == []
