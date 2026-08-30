"""Regression coverage for the static multi-resolution cluster flow figure."""

from types import SimpleNamespace

import pandas as pd


def test_resolution_flow_uses_proportional_alluvial_bands(tmp_path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import PathPatch, Rectangle
    from modules.clustering import ClusteringAnalysis

    # The middle resolution splits both coarse clusters, and the last one
    # splits only one of those branches.  This makes it possible to verify
    # that adjacent-resolution links, rather than cluster IDs, define flows.
    adata = SimpleNamespace(
        obs=pd.DataFrame({
            'leiden_0.2': ['0'] * 6 + ['1'] * 6,
            'leiden_0.6': ['0'] * 4 + ['1'] * 2 + ['1'] * 3 + ['2'] * 3,
            'leiden_1.0': ['0'] * 3 + ['1'] + ['2'] * 2 + ['1'] * 3 + ['3'] * 3,
        }),
        n_obs=12,
    )
    analysis = ClusteringAnalysis(
        str(tmp_path), {'show_resolution_sankey': True}, lambda *_: None,
    )
    captured = {}

    def capture_figure(figure, *args, **kwargs):
        captured['figure'] = figure
        captured['filename'] = args[1]
        return []

    analysis.save_matplotlib_figure = capture_figure
    result_files = []
    analysis._plot_resolution_sankey(
        adata, [0.2, 0.6, 1.0], 'leiden', str(tmp_path), result_files,
    )

    figure = captured['figure']
    try:
        axis = figure.axes[0]
        node_patches = [patch for patch in axis.patches if isinstance(patch, Rectangle)]
        band_patches = [patch for patch in axis.patches if isinstance(patch, PathPatch)]

        # One node per observed cluster and one proportional band per nonzero
        # crosstab entry across each adjacent resolution pair.
        assert len(node_patches) == 2 + 3 + 4
        assert len(band_patches) == 4 + 5
        assert all(patch.get_path().codes.tolist().count(4) == 6 for patch in band_patches)
        assert max(patch.get_height() for patch in node_patches) > min(
            patch.get_height() for patch in node_patches
        )

        assert captured['filename'] == 'cluster_resolution_sankey.png'
        assert figure._native_layout_rect == (0, 0.08, 1, 0.93)
        assert any('Band width = shared cells' in text.get_text() for text in figure.texts)
    finally:
        plt.close(figure)
