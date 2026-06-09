import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
import pandas as pd
import tempfile
import json
from modules.visualization import (
    compute_gene_variability, transform_heatmap_data, cluster_heatmap,
    build_annotation_bar, save_plotly_json,
)


class TestComputeGeneVariability:
    @pytest.fixture
    def sample_data(self):
        return np.array([
            [1.0, 10.0, 0.1],
            [2.0, 20.0, 0.2],
            [3.0, 30.0, 0.3],
            [4.0, 40.0, 0.4],
            [5.0, 50.0, 0.5],
        ])

    def test_var(self, sample_data):
        result = compute_gene_variability(sample_data, metric='var')
        assert len(result) == 3
        np.testing.assert_almost_equal(result[0], np.var([1, 2, 3, 4, 5]))
        np.testing.assert_almost_equal(result[1], np.var([10, 20, 30, 40, 50]))

    def test_mad(self, sample_data):
        result = compute_gene_variability(sample_data, metric='mad')
        assert len(result) == 3
        assert all(v >= 0 for v in result)
        uniform = np.array([[1, 5], [1, 5], [1, 5], [1, 5]])
        result_u = compute_gene_variability(uniform, metric='mad')
        assert result_u[0] == 0.0

    def test_cv(self, sample_data):
        result = compute_gene_variability(sample_data, metric='cv')
        assert len(result) == 3
        assert all(v >= 0 for v in result)

    def test_range(self, sample_data):
        result = compute_gene_variability(sample_data, metric='range')
        np.testing.assert_almost_equal(result[0], 4.0)
        np.testing.assert_almost_equal(result[1], 40.0)

    def test_unknown_metric_raises(self, sample_data):
        with pytest.raises(ValueError):
            compute_gene_variability(sample_data, metric='unknown')


class TestTransformHeatmapData:
    @pytest.fixture
    def sample_data(self):
        return np.array([
            [10.0, 100.0, 5.0],
            [20.0, 200.0, np.nan],
            [30.0, 300.0, 15.0],
            [40.0, 400.0, 20.0],
        ])

    def test_zscore(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='zscore')
        np.testing.assert_almost_equal(result[:, 0].mean(), 0.0, decimal=10)
        assert abs(result[:, 0].std() - 1.0) < 0.01

    def test_center(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='center')
        np.testing.assert_almost_equal(result[:, 0].mean(), 0.0, decimal=10)
        assert result[:, 0].std() != 1.0

    def test_none_scaling(self, sample_data):
        # row_scaling='none' + default clip_range=(-3,3) clips values; NaN stays NaN
        result = transform_heatmap_data(sample_data.copy(), row_scaling='none', clip_range=None)
        assert np.any(np.isnan(result))  # NaN preserved when missing_value='ignore'

    def test_winsorize(self):
        data = np.array([[1, 100], [2, 200], [3, 300], [1000, 400], [5, 500]])
        result = transform_heatmap_data(data.copy(), row_scaling='none', winsorize='5pct')
        assert result[3, 0] < 1000

    def test_missing_mean_fill(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='none',
                                        missing_value='mean_fill', clip_range=None)
        assert not np.any(np.isnan(result))
        col2_mean = np.nanmean(sample_data[:, 2])
        assert abs(result[1, 2] - col2_mean) < 0.01

    def test_missing_zero_fill(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='none', missing_value='zero_fill')
        assert result[1, 2] == 0.0

    def test_clip_range(self, sample_data):
        # Use zero_fill to avoid NaN issues in comparison
        result = transform_heatmap_data(sample_data.copy(), row_scaling='zscore',
                                        clip_range=(-1, 1), missing_value='zero_fill')
        assert np.all(result >= -1.0)
        assert np.all(result <= 1.0)

    def test_no_clip(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='zscore', clip_range=None)
        assert np.any(np.abs(result) > 0)


class TestClusterHeatmap:
    @pytest.fixture
    def cluster_data(self):
        return np.array([
            [1, 1, -1, -1],
            [1, 1, -1, -1],
            [1, 1, -1, -1],
            [-1, -1, 1, 1],
            [-1, -1, 1, 1],
            [-1, -1, 1, 1],
        ], dtype=float)

    def test_euclidean_ward(self, cluster_data):
        order = cluster_heatmap(cluster_data, method='ward', metric='euclidean')
        assert len(order) == 6
        assert set(order) == {0, 1, 2, 3, 4, 5}
        a_positions = [order.index(i) for i in [0, 1, 2]]
        b_positions = [order.index(i) for i in [3, 4, 5]]
        assert max(a_positions) - min(a_positions) <= 2
        assert max(b_positions) - min(b_positions) <= 2

    def test_pearson_metric(self, cluster_data):
        order = cluster_heatmap(cluster_data, method='complete', metric='pearson')
        assert len(order) == 6
        assert set(order) == {0, 1, 2, 3, 4, 5}

    def test_cosine_metric(self, cluster_data):
        order = cluster_heatmap(cluster_data, method='average', metric='cosine')
        assert len(order) == 6

    def test_single_sample(self):
        data = np.array([[1, 2, 3]])
        order = cluster_heatmap(data)
        assert order == [0]

    def test_two_samples(self):
        data = np.array([[1, 2], [3, 4]])
        order = cluster_heatmap(data)
        assert len(order) == 2


class TestBuildAnnotationBar:
    def test_single_column(self):
        obs = pd.DataFrame({'group': ['A', 'A', 'B', 'B']}, index=['s1', 's2', 's3', 's4'])
        fig_data, unique_groups = build_annotation_bar(obs, ['group'], sample_order=[0, 1, 2, 3])
        assert 'group' in fig_data
        assert len(unique_groups['group']) == 2

    def test_multiple_columns(self):
        obs = pd.DataFrame({
            'group': ['A', 'A', 'B', 'B'],
            'batch': ['b1', 'b2', 'b1', 'b2']
        }, index=['s1', 's2', 's3', 's4'])
        fig_data, unique_groups = build_annotation_bar(obs, ['group', 'batch'])
        assert 'group' in fig_data
        assert 'batch' in fig_data

    def test_empty_columns(self):
        obs = pd.DataFrame({'group': ['A', 'B']}, index=['s1', 's2'])
        fig_data, unique_groups = build_annotation_bar(obs, [])
        assert len(fig_data) == 0

    def test_missing_column_skipped(self):
        obs = pd.DataFrame({'group': ['A', 'B']}, index=['s1', 's2'])
        fig_data, _ = build_annotation_bar(obs, ['group', 'nonexistent'])
        assert 'group' in fig_data
        assert 'nonexistent' not in fig_data


class TestSavePlotlyJson:
    def test_save_creates_file(self):
        import plotly.graph_objects as go
        fig = go.Figure(go.Heatmap(z=[[1, 2], [3, 4]]))
        with tempfile.TemporaryDirectory() as tmpdir:
            result_files = []
            save_plotly_json(fig, tmpdir, 'test_heatmap.json', result_files,
                            category='heatmap', label='测试热图')
            assert os.path.exists(os.path.join(tmpdir, 'test_heatmap.json'))
            assert len(result_files) == 1
            assert result_files[0]['category'] == 'heatmap'
            assert result_files[0]['label'] == '测试热图'

    def test_save_valid_json(self):
        import plotly.graph_objects as go
        fig = go.Figure(go.Heatmap(z=[[1, 2], [3, 4]]))
        with tempfile.TemporaryDirectory() as tmpdir:
            result_files = []
            save_plotly_json(fig, tmpdir, 'test.json', result_files)
            with open(os.path.join(tmpdir, 'test.json')) as f:
                data = json.load(f)
            assert 'data' in data
