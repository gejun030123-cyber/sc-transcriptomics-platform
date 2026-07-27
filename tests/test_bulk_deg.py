# tests/test_bulk_deg.py
"""Tests for modules/bulk_deg.py — comparison parsing helpers."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from modules.bulk_deg import (
    _parse_comparisons,
    _parse_custom_groups,
    _welch_ttest_log_expression,
)


class TestParseComparisons:
    """测试 _parse_comparisons 函数。"""

    def test_single_comparison(self):
        """单个比较。"""
        assert _parse_comparisons('A-vs-B') == [('A', 'B')]

    def test_multiple_comparisons(self):
        """分号分隔的多个比较。"""
        result = _parse_comparisons('A-vs-B;C-vs-D')
        assert result == [('A', 'B'), ('C', 'D')]

    def test_newline_separator(self):
        """换行符分隔。"""
        result = _parse_comparisons('A-vs-B\nC-vs-D')
        assert result == [('A', 'B'), ('C', 'D')]

    def test_mixed_separators(self):
        """混合分隔符。"""
        result = _parse_comparisons('A-vs-B;C-vs-D\nE-vs-F')
        assert len(result) == 3

    def test_empty_string(self):
        """空字符串 → 空列表。"""
        assert _parse_comparisons('') == []
        assert _parse_comparisons(None) == []

    def test_whitespace_handling(self):
        """带空格的输入。"""
        result = _parse_comparisons('  A -vs- B  ;  C -vs- D  ')
        assert result == [('A', 'B'), ('C', 'D')]

    def test_no_vs_separator(self):
        """不含 -vs- 的项被跳过。"""
        result = _parse_comparisons('A-vs-B;invalid;C-vs-D')
        assert result == [('A', 'B'), ('C', 'D')]

    def test_empty_part_skipped(self):
        """不完整的比较被跳过。"""
        result = _parse_comparisons('A-vs-;-vs-B;C-vs-D')
        assert result == [('C', 'D')]

    def test_group_names_with_special_chars(self):
        """组名含下划线、数字。"""
        result = _parse_comparisons('Treated_3h-vs-Control_0h')
        assert result == [('Treated_3h', 'Control_0h')]

    def test_multiple_vs_in_name(self):
        """组名中含多个 -vs- → 只按第一个分割。"""
        result = _parse_comparisons('A-vs-B-vs-C')
        assert result == [('A', 'B-vs-C')]


class TestParseCustomGroups:
    """测试 _parse_custom_groups 函数。"""

    def test_single_group(self):
        """单个自定义组。"""
        result = _parse_custom_groups('High=Treated_1h+Treated_3h')
        assert result == {'High': ['Treated_1h', 'Treated_3h']}

    def test_multiple_groups(self):
        """多个组。"""
        result = _parse_custom_groups('High=Treated_1h+Treated_3h\nLow=Ctrl')
        assert result == {'High': ['Treated_1h', 'Treated_3h'], 'Low': ['Ctrl']}

    def test_single_member(self):
        """单个成员。"""
        result = _parse_custom_groups('Ctrl=Control')
        assert result == {'Ctrl': ['Control']}

    def test_empty_string(self):
        """空字符串 → 空 dict。"""
        assert _parse_custom_groups('') == {}
        assert _parse_custom_groups(None) == {}

    def test_no_equals_skipped(self):
        """不含 = 的行被跳过。"""
        result = _parse_custom_groups('High=A+B\ninvalid\nLow=C')
        assert 'High' in result
        assert 'Low' in result
        assert len(result) == 2

    def test_whitespace_handling(self):
        """带空格的输入。"""
        result = _parse_custom_groups('  High = A + B  ')
        assert result == {'High': ['A', 'B']}

    def test_empty_name_skipped(self):
        """空组名被跳过。"""
        result = _parse_custom_groups('=A+B')
        assert result == {}

    def test_empty_members_skipped(self):
        """空成员被跳过。"""
        result = _parse_custom_groups('High=')
        assert result == {}


def test_welch_ttest_log_expression_uses_mean_difference_as_log2fc():
    """Log-scale input must use a difference, not a ratio of log values."""
    import pandas as pd

    data = pd.DataFrame({
        't1': [4.0, 2.0], 't2': [4.1, 2.1], 't3': [3.9, 1.9],
        'c1': [2.0, 2.0], 'c2': [2.1, 2.1], 'c3': [1.9, 1.9],
    }, index=['changed', 'stable'])

    result = _welch_ttest_log_expression(
        data, ['t1', 't2', 't3'], ['c1', 'c2', 'c3'])

    assert result.loc['changed', 'log2FC'] == pytest.approx(2.0)
    assert result.loc['stable', 'log2FC'] == pytest.approx(0.0)
    assert result.loc['changed', 'qvalue'] < 0.05


def test_volcano_display_ceiling_handles_padj_underflow():
    """Extreme/zero padj values must not stretch the panel to y=300."""
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from modules.bulk_deg import _draw_bulk_volcano, _volcano_y_limit

    padj = np.array([0.0, 1e-80, 0.01, 0.5])
    assert _volcano_y_limit(padj, 0.05) <= 30

    fig, ax = plt.subplots()
    info = _draw_bulk_volcano(
        ax,
        pd.DataFrame({
            'gene': ['A', 'B', 'C', 'D'],
            'log2FC': [2.0, -1.5, 0.2, 0.0],
            'padj': padj,
            'regulation': ['Up', 'Down', 'NS', 'NS'],
        }),
        pval_threshold=0.05,
        fc_threshold=2.0,
        top_n=0,
    )
    assert info['n_clipped'] == 2
    assert ax.get_ylim()[1] <= 30
    plt.close(fig)


def test_volcano_gene_labels_use_spaced_boxes_and_arrows():
    """Dense top genes stay readable: cap labels and separate each side."""
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from modules.bulk_deg import _draw_bulk_volcano

    n_each = 8
    deg_df = pd.DataFrame({
        'gene': [f'UP{i}' for i in range(n_each)] + [f'DOWN{i}' for i in range(n_each)],
        'log2FC': np.r_[np.linspace(1.2, 2.0, n_each), np.linspace(-1.2, -2.0, n_each)],
        'padj': np.r_[np.linspace(1e-9, 8e-9, n_each), np.linspace(1e-9, 8e-9, n_each)],
        'regulation': ['Up'] * n_each + ['Down'] * n_each,
    })

    fig, ax = plt.subplots()
    info = _draw_bulk_volcano(ax, deg_df, top_n=50, show_legend=False)
    fig.canvas.draw()
    labels = [text for text in ax.texts if text.get_text().startswith(('UP', 'DOWN'))]

    assert len(labels) == 10  # visual labels are capped; the Top-DEG CSV is unchanged
    assert all(label.get_bbox_patch() is not None for label in labels)
    assert all(label.arrow_patch is not None for label in labels)
    for alignment in ('left', 'right'):
        y_positions = sorted(
            label.get_position()[1] for label in labels
            if label.get_horizontalalignment() == alignment
        )
        assert all(
            second - first >= max(0.55, info['y_limit'] * 0.04) - 1e-7
            for first, second in zip(y_positions, y_positions[1:])
        )
    plt.close(fig)
