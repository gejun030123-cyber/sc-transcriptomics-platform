# tests/test_bulk_deg.py
"""Tests for modules/bulk_deg.py — comparison parsing helpers."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from modules.bulk_deg import (
    _build_comparison_matrix,
    _format_volcano_contrast_title,
    _parse_comparisons,
    _parse_custom_groups,
    _welch_ttest_log_expression,
)


def test_volcano_contrast_title_formats_common_condition_strata():
    assert _format_volcano_contrast_title('NH4Cl_B', 'Ctr_B') == 'NH₄Cl vs Control (B)'
    assert _format_volcano_contrast_title('PEA_En', 'Ctr_En') == 'PEA vs Control (En)'
    assert _format_volcano_contrast_title('treated', 'vehicle') == 'treated vs vehicle'


def test_comparison_matrix_joins_repeated_gene_symbols_by_feature_id():
    """Repeated display symbols must not create a Cartesian-product merge."""
    import pandas as pd

    first = pd.DataFrame({
        'gene_id': ['id_1', 'id_2'], 'gene': ['DUP', 'DUP'],
        'gene_name': ['DUP', 'DUP'], 'log2FC': [1.0, 2.0],
        'padj': [0.01, 0.02], 'regulation': ['Up', 'Up'],
    })
    second = pd.DataFrame({
        'gene_id': ['id_1', 'id_2'], 'gene': ['DUP', 'DUP'],
        'gene_name': ['DUP', 'DUP'], 'log2FC': [-1.0, -2.0],
        'padj': [0.03, 0.04], 'regulation': ['Down', 'Down'],
    })

    result = _build_comparison_matrix([first, second], ['A-vs-B', 'C-vs-D'])

    assert result.shape[0] == 2
    assert result['gene_id'].tolist() == ['id_1', 'id_2']
    assert 'A-vs-B_log2FC' in result
    assert 'C-vs-D_log2FC' in result


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
    assert _volcano_y_limit(padj, 0.05) == 18

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
    assert ax.get_ylim()[1] == 18
    assert 'log' in ax.get_ylabel()
    plt.close(fig)


def test_volcano_gene_labels_use_short_leaders_without_boxes():
    """Dense top genes stay readable without turning into an annotation map."""
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

    assert len(labels) == 6  # visual labels are capped; the Top-DEG CSV is unchanged
    assert all(label.get_bbox_patch() is None for label in labels)
    assert all(label.arrow_patch is not None for label in labels)
    renderer = fig.canvas.get_renderer()
    from matplotlib.text import Text
    boxes = [Text.get_window_extent(label, renderer) for label in labels]
    assert not any(
        left.overlaps(right)
        for index, left in enumerate(boxes) for right in boxes[index + 1:]
    )
    plt.close(fig)


def test_volcano_keeps_extreme_fold_changes_at_compact_x_axis_edges():
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from modules.bulk_deg import _draw_bulk_volcano

    fig, ax = plt.subplots()
    info = _draw_bulk_volcano(
        ax,
        pd.DataFrame({
            'gene': ['far_down', 'center', 'far_up'],
            'log2FC': [-6.0, 0.0, 4.0],
            'padj': [0.002, 0.8, 0.003],
            'regulation': ['Down', 'NS', 'Up'],
        }),
        top_n=0,
    )

    points = np.vstack([
        collection.get_offsets() for collection in ax.collections
        if len(collection.get_offsets())
    ])
    assert info['x_limit'] == 2.5
    assert info['n_x_clipped'] == 2
    assert np.sum(np.isclose(np.abs(points[:, 0]), 2.5 * 0.985)) == 2
    plt.close(fig)


def test_volcano_custom_view_range_retains_clipped_edge_markers():
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from modules.bulk_deg import _draw_bulk_volcano

    fig, ax = plt.subplots()
    info = _draw_bulk_volcano(
        ax,
        pd.DataFrame({
            'gene': ['far_down', 'center', 'far_up', 'high_fdr'],
            'log2FC': [-6.0, 0.0, 4.0, 0.5],
            'padj': [0.002, 0.8, 0.003, 1e-12],
            'regulation': ['Down', 'NS', 'Up', 'Up'],
        }),
        top_n=0, x_min=-1.5, x_max=2.0, y_min=0.0, y_limit=5.0,
    )

    points = np.vstack([
        collection.get_offsets() for collection in ax.collections
        if len(collection.get_offsets())
    ])
    assert ax.get_xlim() == (-1.5, 2.0)
    assert ax.get_ylim() == (0.0, 5.0)
    assert info['n_x_clipped'] == 2
    assert info['n_clipped'] == 1
    assert np.any(np.isclose(points[:, 0], -1.5 + 3.5 * 0.0075))
    assert np.any(np.isclose(points[:, 0], 2.0 - 3.5 * 0.0075))
    assert np.any(np.isclose(points[:, 1], 5.0 * 0.985))
    plt.close(fig)
