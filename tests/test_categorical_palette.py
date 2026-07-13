"""Color palettes should be categorical, stable, and sufficiently distinct."""
import os
import sys

import pandas as pd
import pytest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def test_categorical_color_map_uses_one_color_per_category():
    ad = pytest.importorskip('anndata')
    from modules.visualization import categorical_color_map

    adata = ad.AnnData(X=np.array([[1], [2], [3], [4]]),
                       obs=pd.DataFrame({'leiden': pd.Categorical(['0', '1', '2', '1'])}))
    palette = categorical_color_map(adata, 'leiden')
    assert set(palette) == {'0', '1', '2'}
    assert len(set(palette.values())) == 3
    assert 'leiden_colors' in adata.uns
