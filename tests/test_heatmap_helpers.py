import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from modules.visualization import compute_gene_variability, transform_heatmap_data


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
