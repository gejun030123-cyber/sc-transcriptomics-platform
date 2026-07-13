# tests/test_base.py
"""Tests for modules/base.py — BaseAnalysis helper methods."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
import pytest
import anndata
from modules.base import BaseAnalysis


# --- 测试用具体子类 ---

class _StubAnalysis(BaseAnalysis):
    """用于测试的最小化具体子类。"""
    MODULE_NAME = '_stub'

    def run(self, input_path: str) -> dict:
        return {}


def _make_analysis(params=None):
    """创建 _StubAnalysis 实例，带可选参数。"""
    return _StubAnalysis(
        project_dir='/tmp/test_project',
        params=params or {},
        progress_callback=lambda pct, msg: None,
    )


def test_progress_allows_missing_callback():
    """Branch/agent runners may omit progress callbacks; progress should be a no-op."""
    analysis = _StubAnalysis(
        project_dir='/tmp/test_project',
        params={},
        progress_callback=None,
    )
    analysis.progress(10, 'running')


def test_save_matplotlib_figure_writes_publication_formats(tmp_path):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    files = _make_analysis().save_matplotlib_figure(
        fig, str(tmp_path), 'quality_plot.png', 'qc', 'Quality plot'
    )
    plt.close(fig)

    assert {item['file_type'] for item in files} == {'png', 'svg'}
    assert (tmp_path / 'quality_plot.png').is_file()
    assert (tmp_path / 'quality_plot.svg').is_file()
    assert all(item['label'] == 'Quality plot' for item in files)


def _make_adata(n_obs=10, obs_dict=None):
    """创建测试用 AnnData。"""
    X = np.random.rand(n_obs, 5)
    obs = pd.DataFrame(obs_dict or {}, index=[f'cell_{i}' for i in range(n_obs)])
    adata = anndata.AnnData(X=X, obs=obs)
    adata.var_names = [f'Gene{i}' for i in range(5)]
    return adata


# --- apply_filters tests ---

class TestApplyFilters:
    """测试 apply_filters 方法。"""

    def test_no_filters(self):
        """无过滤规则 → adata 不变。"""
        a = _make_analysis()
        adata = _make_adata(10)
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 10

    def test_ge_filter(self):
        """>= 过滤：移除低于阈值的行。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'n_genes', 'op': '>=', 'value': 5}
        ]}})
        adata = _make_adata(10, {'n_genes': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 6  # 5,6,7,8,9,10

    def test_le_filter(self):
        """<= 过滤。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'mt_pct', 'op': '<=', 'value': 10}
        ]}})
        adata = _make_adata(5, {'mt_pct': [5.0, 10.0, 10.1, 15.0, 20.0]})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 2

    def test_eq_filter(self):
        """== 过滤。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'batch', 'op': '==', 'value': 'A'}
        ]}})
        adata = _make_adata(4, {'batch': ['A', 'B', 'A', 'C']})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 2

    def test_neq_filter(self):
        """!= 过滤。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'batch', 'op': '!=', 'value': 'B'}
        ]}})
        adata = _make_adata(4, {'batch': ['A', 'B', 'A', 'C']})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 3

    def test_in_filter(self):
        """in 过滤：保留列表中的值。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'celltype', 'op': 'in', 'value': ['T', 'B']}
        ]}})
        adata = _make_adata(5, {'celltype': ['T', 'B', 'NK', 'T', 'Macro']})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 3

    def test_not_in_filter(self):
        """not_in 过滤。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'celltype', 'op': 'not_in', 'value': ['NK', 'Macro']}
        ]}})
        adata = _make_adata(5, {'celltype': ['T', 'B', 'NK', 'T', 'Macro']})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 3

    def test_between_filter(self):
        """between 过滤：保留 [min, max] 范围内的值。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'score', 'op': 'between', 'value': [0.3, 0.7]}
        ]}})
        adata = _make_adata(5, {'score': [0.1, 0.3, 0.5, 0.7, 0.9]})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 3  # 0.3, 0.5, 0.7

    def test_missing_column_skipped(self):
        """列不存在 → 跳过该规则，不崩溃。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'nonexistent', 'op': '>=', 'value': 5}
        ]}})
        adata = _make_adata(5)
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 5

    def test_unsupported_operator(self):
        """不支持的操作符 → 跳过，不崩溃。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'score', 'op': 'LIKE', 'value': 'test'}
        ]}})
        adata = _make_adata(5, {'score': [1, 2, 3, 4, 5]})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 5

    def test_chained_filters(self):
        """多条规则串联：先过滤 n_genes，再过滤 mt_pct。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'n_genes', 'op': '>=', 'value': 3},
            {'column': 'mt_pct', 'op': '<=', 'value': 15},
        ]}})
        adata = _make_adata(5, {
            'n_genes': [1, 2, 3, 4, 5],
            'mt_pct': [5.0, 10.0, 20.0, 8.0, 12.0],
        })
        result = a.apply_filters(adata, 'qc')
        # n_genes >= 3: cell_2, cell_3, cell_4
        # mt_pct <= 15: cell_2(20→排除), cell_3(8), cell_4(12)
        assert result.n_obs == 2

    def test_in_with_string_value(self):
        """in 操作符收到字符串而非列表 → 跳过，不崩溃。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'batch', 'op': 'in', 'value': 'A'}
        ]}})
        adata = _make_adata(3, {'batch': ['A', 'B', 'C']})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 3  # 跳过，不变

    def test_between_with_invalid_value(self):
        """between 操作符收到非列表 → 跳过，不崩溃。"""
        a = _make_analysis({'_filters': {'qc': [
            {'column': 'score', 'op': 'between', 'value': 0.5}
        ]}})
        adata = _make_adata(3, {'score': [0.1, 0.5, 0.9]})
        result = a.apply_filters(adata, 'qc')
        assert result.n_obs == 3


# --- get_plotly_layout tests ---

class TestGetPlotlyLayout:
    """测试 get_plotly_layout 方法。"""

    def test_default_layout(self):
        """无 _visualization 参数 → 使用默认主题。"""
        a = _make_analysis()
        layout = a.get_plotly_layout('Test Title')
        assert layout['title'] == 'Test Title'
        assert layout['width'] == 800
        assert layout['height'] == 500
        assert layout['font']['family'] == 'Arial'

    def test_custom_theme(self):
        """指定 nature 主题。"""
        a = _make_analysis({'_visualization': {'theme': 'nature'}})
        layout = a.get_plotly_layout()
        assert layout['font']['family'] == 'Helvetica'

    def test_overrides(self):
        """overrides 覆盖默认值。"""
        a = _make_analysis()
        layout = a.get_plotly_layout(width=1200)
        assert layout['width'] == 1200


# --- get_viz_params tests ---

class TestGetVizParams:
    """测试 get_viz_params 方法。"""

    def test_default_params(self):
        """无 _visualization → 使用默认值。"""
        a = _make_analysis()
        viz = a.get_viz_params()
        assert viz['umap_point_size'] == 5
        assert viz['umap_opacity'] == 0.7

    def test_custom_params(self):
        """自定义可视化参数。"""
        a = _make_analysis({'_visualization': {'umap_point_size': 10}})
        viz = a.get_viz_params()
        assert viz['umap_point_size'] == 10


# --- load_adata / save_output 测试 ---

class TestLoadAdata:
    """load_adata 方法测试。"""

    def test_load_valid_h5ad(self, tmp_path):
        """加载有效 h5ad 文件应返回 AnnData。"""
        import scanpy as sc
        adata = _make_adata()
        h5ad_path = str(tmp_path / "test.h5ad")
        adata.write_h5ad(h5ad_path)

        a = _make_analysis()
        result = a.load_adata(h5ad_path)
        assert isinstance(result, anndata.AnnData)
        assert result.n_obs == 10
        assert result.n_vars == 5

    def test_load_nonexistent_file(self):
        """加载不存在的文件应抛出异常。"""
        a = _make_analysis()
        with pytest.raises((FileNotFoundError, OSError)):
            a.load_adata("/nonexistent/path/file.h5ad")

    def test_load_non_h5ad_raises(self, tmp_path):
        """加载非 h5ad 文件应抛出 OSError。"""
        csv_path = str(tmp_path / "data.csv")
        with open(csv_path, 'w') as f:
            f.write("gene,s1,s2\nG1,1,2\n")

        a = _make_analysis()
        with pytest.raises((OSError, ValueError)):
            a.load_adata(csv_path)


class TestSaveOutput:
    """save_output 方法测试。"""

    def test_save_creates_file(self, tmp_path):
        """save_output 应创建 h5ad 文件并返回路径。"""
        adata = _make_adata()
        a = _make_analysis(params=None)
        a.project_dir = str(tmp_path)

        output_path = a.save_output(adata, 'test_module')
        assert os.path.exists(output_path)
        assert output_path.endswith('test_module_output.h5ad')
        assert 'intermediate' in output_path

    def test_save_creates_intermediate_dir(self, tmp_path):
        """save_output 应自动创建 intermediate 目录。"""
        adata = _make_adata()
        a = _make_analysis(params=None)
        a.project_dir = str(tmp_path)

        intermediate = os.path.join(str(tmp_path), 'intermediate')
        assert not os.path.exists(intermediate)

        a.save_output(adata, 'mod')
        assert os.path.exists(intermediate)

    def test_save_roundtrip(self, tmp_path):
        """保存后重新加载应得到相同数据。"""
        import scanpy as sc
        adata = _make_adata()
        a = _make_analysis(params=None)
        a.project_dir = str(tmp_path)

        output_path = a.save_output(adata, 'roundtrip')
        loaded = sc.read_h5ad(output_path)
        assert loaded.n_obs == adata.n_obs
        assert loaded.n_vars == adata.n_vars
