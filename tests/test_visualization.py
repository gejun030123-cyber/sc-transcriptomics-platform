# tests/test_visualization.py
"""Tests for modules/visualization.py — data transformation and clustering helpers."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest


class TestTransformHeatmapData:
    """测试 transform_heatmap_data 函数。"""

    def test_zscore_basic(self):
        """z-score 标准化后均值接近 0。"""
        from modules.visualization import transform_heatmap_data
        data = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        result = transform_heatmap_data(data, row_scaling='zscore')
        assert result.shape == data.shape
        np.testing.assert_allclose(result.mean(axis=0), 0, atol=1e-6)

    def test_center(self):
        """center 后均值为 0。"""
        from modules.visualization import transform_heatmap_data
        data = np.array([[10.0, 20.0], [30.0, 40.0]])
        result = transform_heatmap_data(data, row_scaling='center')
        np.testing.assert_allclose(result.mean(axis=0), 0, atol=1e-6)

    def test_none_scaling(self):
        """无标准化 → 数据不变（除 clip 外）。"""
        from modules.visualization import transform_heatmap_data
        data = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = transform_heatmap_data(data, row_scaling='none', clip_range=None)
        np.testing.assert_allclose(result, data)

    def test_clip_range(self):
        """截断到 [-2, 2]。"""
        from modules.visualization import transform_heatmap_data
        data = np.array([[-10.0, -1.0], [0.0, 1.0], [10.0, 2.0]])
        result = transform_heatmap_data(data, row_scaling='none', clip_range=(-2, 2))
        assert result.max() <= 2.0
        assert result.min() >= -2.0

    def test_missing_value_mean_fill(self):
        """NaN 值用列均值填充。"""
        from modules.visualization import transform_heatmap_data
        data = np.array([[1.0, np.nan], [3.0, 4.0], [5.0, 6.0]])
        result = transform_heatmap_data(data, row_scaling='none', missing_value='mean_fill',
                                        clip_range=None)
        assert not np.isnan(result).any()
        # NaN 应被填充为 (4+6)/2 = 5
        assert result[0, 1] == pytest.approx(5.0)

    def test_missing_value_zero_fill(self):
        """NaN 值用 0 填充。"""
        from modules.visualization import transform_heatmap_data
        data = np.array([[1.0, np.nan], [3.0, 4.0]])
        result = transform_heatmap_data(data, row_scaling='none', missing_value='zero_fill',
                                        clip_range=None)
        assert result[0, 1] == 0.0

    def test_winsorize(self):
        """Winsorize 截断极端值。"""
        from modules.visualization import transform_heatmap_data
        np.random.seed(42)
        data = np.random.randn(100, 3)
        data[0, 0] = 100  # 极端值
        result = transform_heatmap_data(data, row_scaling='none', winsorize='1pct',
                                        clip_range=None)
        assert result[0, 0] < 100  # 极端值应被截断


def test_plotly_umap_uses_shared_style_and_existing_category_colours():
    import anndata
    import pandas as pd

    from modules.visualization import umap_scatter

    adata = anndata.AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({'cell_type': ['T', 'B', 'T', 'B']}),
    )
    adata.obsm['X_umap'] = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    adata.uns['cell_type_colors'] = ['#123456', '#ABCDEF']

    figure = umap_scatter(
        adata, 'cell_type',
        viz_params={'figure_width': 1000, 'figure_height': 700},
    )

    assert figure['layout']['width'] == 1000
    assert figure['layout']['height'] == 700
    assert figure['layout']['xaxis']['showticklabels'] is False
    assert figure['data'][0]['marker']['color'] == '#123456'

class TestClusterHeatmap:
    """测试 cluster_heatmap 函数。"""

    def test_basic_clustering(self):
        """基本聚类返回正确数量的索引。"""
        from modules.visualization import cluster_heatmap
        data = np.random.rand(10, 5)
        order = cluster_heatmap(data)
        assert len(order) == 10
        assert set(order) == set(range(10))

    def test_single_sample(self):
        """单样本 → 返回 [0]。"""
        from modules.visualization import cluster_heatmap
        data = np.array([[1.0, 2.0, 3.0]])
        order = cluster_heatmap(data)
        assert order == [0]

    def test_two_samples(self):
        """两个样本 → 返回 2 个索引。"""
        from modules.visualization import cluster_heatmap
        data = np.array([[1.0, 2.0], [10.0, 20.0]])
        order = cluster_heatmap(data)
        assert len(order) == 2

    def test_pearson_metric(self):
        """Pearson 距离度量。"""
        from modules.visualization import cluster_heatmap
        data = np.random.rand(8, 5)
        order = cluster_heatmap(data, metric='pearson')
        assert len(order) == 8

    def test_cosine_metric(self):
        """Cosine 距离度量。"""
        from modules.visualization import cluster_heatmap
        data = np.random.rand(8, 5)
        order = cluster_heatmap(data, metric='cosine')
        assert len(order) == 8


class TestComputeVariability:
    """测试 compute_gene_variability 函数（已在 test_heatmap_helpers.py 中有部分测试）。"""

    def test_zero_variance_gene(self):
        """方差为零的基因 → 返回 0。"""
        from modules.visualization import compute_gene_variability
        data = np.array([[5.0, 1.0], [5.0, 2.0], [5.0, 3.0]])
        result = compute_gene_variability(data, metric='var')
        assert result[0] == pytest.approx(0.0)
        assert result[1] > 0

    def test_unknown_metric_raises(self):
        """未知度量 → 抛出 ValueError。"""
        from modules.visualization import compute_gene_variability
        with pytest.raises(ValueError, match='未知变异度量'):
            compute_gene_variability(np.array([[1.0]]), metric='invalid')
