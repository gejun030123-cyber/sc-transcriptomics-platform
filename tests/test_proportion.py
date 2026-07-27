# tests/test_proportion.py
"""Tests for modules/proportion.py — _run_stat_test."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
import pytest
import anndata as ad


class TestRunStatTest:
    """测试 _run_stat_test 统计检验函数。"""

    def _make_contingency(self, data, rows=None, cols=None):
        """创建列联表 DataFrame。"""
        rows = rows or [f'r{i}' for i in range(len(data))]
        cols = cols or [f'c{i}' for i in range(len(data[0]))]
        return pd.DataFrame(data, index=rows, columns=cols)

    def test_chi_square_independent(self):
        """独立分组 → p 值应较大（不显著）。"""
        from modules.proportion import _run_stat_test
        ct = self._make_contingency([[50, 50], [50, 50]])
        stat, pval = _run_stat_test(ct, 'chi_square')
        assert stat == pytest.approx(0.0)
        assert pval > 0.05

    def test_chi_square_dependent(self):
        """强烈关联 → p 值应极小。"""
        from modules.proportion import _run_stat_test
        ct = self._make_contingency([[100, 10], [10, 100]])
        stat, pval = _run_stat_test(ct, 'chi_square')
        assert stat > 50
        assert pval < 0.001

    def test_fisher_exact_2x2(self):
        """2x2 列联表 → 使用 Fisher 精确检验。"""
        from modules.proportion import _run_stat_test
        ct = self._make_contingency([[10, 0], [0, 10]])
        stat, pval = _run_stat_test(ct, 'fisher_exact')
        assert pval < 0.05  # 完全分离，应显著

    def test_fisher_exact_fallback_to_chi2(self):
        """非 2x2 列联表 + fisher_exact → 降级为卡方。"""
        from modules.proportion import _run_stat_test
        ct = self._make_contingency([[10, 20, 30], [40, 50, 60]])
        stat, pval = _run_stat_test(ct, 'fisher_exact')
        assert stat > 0
        assert 0 < pval <= 1

    def test_permutation_test(self):
        """置换检验 → p 值合理。"""
        from modules.proportion import _run_stat_test
        np.random.seed(42)
        ct = self._make_contingency([[100, 10], [10, 100]])
        stat, pval = _run_stat_test(ct, 'permutation', n_permutations=500)
        assert pval < 0.05  # 强关联应显著

    def test_permutation_independent(self):
        """独立分组 + 置换检验 → p 值应较大。"""
        from modules.proportion import _run_stat_test
        np.random.seed(42)
        ct = self._make_contingency([[50, 50], [50, 50]])
        stat, pval = _run_stat_test(ct, 'permutation', n_permutations=500)
        assert pval > 0.01  # 不显著

    def test_unknown_test_defaults_to_chi2(self):
        """未知检验类型 → 降级为卡方。"""
        from modules.proportion import _run_stat_test
        ct = self._make_contingency([[10, 20], [30, 40]])
        stat, pval = _run_stat_test(ct, 'unknown_test')
        assert stat > 0
        assert 0 < pval <= 1

    def test_single_row_contingency(self):
        """单行列联表 → 卡方值为 0。"""
        from modules.proportion import _run_stat_test
        ct = self._make_contingency([[10, 20, 30]])
        stat, pval = _run_stat_test(ct, 'chi_square')
        assert stat == pytest.approx(0.0)

    def test_return_types(self):
        """返回值类型正确。"""
        from modules.proportion import _run_stat_test
        ct = self._make_contingency([[10, 20], [30, 40]])
        stat, pval = _run_stat_test(ct, 'chi_square')
        assert isinstance(stat, (int, float))
        assert isinstance(pval, (int, float))


def _sample_level_adata(condition_key='condition'):
    rows = [
        ('s1', 'Ctrl', 'T'), ('s1', 'Ctrl', 'T'), ('s1', 'Ctrl', 'B'),
        ('s2', 'Ctrl', 'T'), ('s2', 'Ctrl', 'B'), ('s2', 'Ctrl', 'B'),
        ('s3', 'Treat', 'T'), ('s3', 'Treat', 'T'), ('s3', 'Treat', 'T'),
        ('s4', 'Treat', 'B'), ('s4', 'Treat', 'B'), ('s4', 'Treat', 'B'),
    ]
    obs = pd.DataFrame(rows, columns=['sample_id', condition_key, 'celltype'])
    obs.index = [f'cell_{index}' for index in range(len(obs))]
    return ad.AnnData(X=np.ones((len(obs), 2)), obs=obs, var=pd.DataFrame(index=['g1', 'g2']))


def test_sample_level_composition_uses_samples_and_keeps_zero_proportions():
    from modules.proportion import _sample_level_composition

    result = _sample_level_composition(
        _sample_level_adata(), 'celltype', 'sample_id', 'condition',
        min_samples_per_condition=2, min_cells_per_sample=2,
    )

    assert result['available'] is True
    assert result['inference_ready'] is True
    assert result['condition_counts'] == {'Ctrl': 2, 'Treat': 2}
    # s3 has no B cells, but it must remain an explicit zero rather than be
    # silently omitted from the treatment distribution.
    s3_b = result['proportions'].query("sample == 's3' and celltype == 'B'")
    assert len(s3_b) == 1
    assert s3_b.iloc[0]['proportion'] == 0
    assert set(result['tests']['test']) == {'mann_whitney_u'}
    assert result['tests']['fdr_bh'].notna().all()


def test_sample_level_composition_refuses_technical_batch_as_condition_without_confirmation():
    from modules.proportion import _sample_level_composition

    result = _sample_level_composition(
        _sample_level_adata(condition_key='batch'), 'celltype', 'sample_id', 'batch',
        min_samples_per_condition=2, min_cells_per_sample=2,
    )

    assert result['available'] is False
    assert result['inference_ready'] is False
    assert any('技术 batch' in warning for warning in result['warnings'])
