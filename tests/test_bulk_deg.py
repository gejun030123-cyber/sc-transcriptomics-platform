# tests/test_bulk_deg.py
"""Tests for modules/bulk_deg.py — comparison parsing helpers."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from modules.bulk_deg import _parse_comparisons, _parse_custom_groups


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
