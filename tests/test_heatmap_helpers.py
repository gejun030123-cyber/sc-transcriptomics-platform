import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from modules.visualization import compute_gene_variability


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
