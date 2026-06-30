# tests/test_proportion.py
"""Tests for modules/proportion.py — _run_stat_test."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
import pytest


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
