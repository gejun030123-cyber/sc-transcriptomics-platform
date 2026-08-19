"""Regression tests for the shared Nature-style figure contract."""

import json
import os


def test_plotly_save_applies_nature_style(tmp_path):
    import plotly.graph_objects as go

    from modules.base import BaseAnalysis
    from modules.figure_style import NATURE_PALETTE

    class Stub(BaseAnalysis):
        MODULE_NAME = 'style_stub'

        def run(self, input_path):
            return {}

    figure = go.Figure(go.Scatter(x=[1, 2], y=[2, 3], mode='markers'))
    result = Stub(str(tmp_path), {}, lambda *_: None).save_plotly_json(
        figure, str(tmp_path), 'styled.json', 'scatter', 'Styled plot'
    )
    payload = json.loads((tmp_path / 'styled.json').read_text(encoding='utf-8'))

    assert result['file_type'] == 'plotly_json'
    assert payload['layout']['paper_bgcolor'] == 'white'
    assert payload['layout']['colorway'][0] == NATURE_PALETTE[0]
    assert payload['layout']['xaxis']['showline'] is True
    assert payload['layout']['yaxis']['showgrid'] is True


def test_visualization_params_normalize_publication_controls():
    from modules.figure_style import normalize_visualization_params

    params = normalize_visualization_params({
        'theme': 'dark', 'font_size': 100, 'figure_width': 10,
        'umap_opacity': 2, 'umap_label_categories': 'true',
        'umap_hide_axes': 'false',
    })

    assert params['bg_color'] == '#1A1A2E'
    assert params['font_size'] == 28
    assert params['figure_width'] == 320
    assert params['umap_opacity'] == 1.0
    assert params['umap_label_categories'] is True
    assert params['umap_hide_axes'] is False


def test_native_umap_labels_replace_redundant_legend():
    import anndata
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    from modules.native_figures import umap_figure

    adata = anndata.AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame({'cell_type': ['T', 'B', 'T', 'B']}),
    )
    adata.obsm['X_umap'] = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    fig = umap_figure(adata, 'cell_type', label_categories=True)

    assert fig.axes[0].get_legend() is None
    assert {text.get_text() for text in fig.axes[0].texts} == {'T', 'B'}
    plt.close(fig)


def test_static_renderer_separates_violin_and_groups_bars(tmp_path):
    from modules.reporting.static_rendering import render_plotly_payload

    payload = {
        'data': [
            {'type': 'violin', 'y': [1, 2, 3, 4], 'name': 'Metric A'},
            {'type': 'violin', 'y': [2, 3, 4, 5], 'name': 'Metric B'},
            {'type': 'bar', 'x': ['Before', 'After'], 'y': [5, 4], 'name': 'Cells'},
            {'type': 'bar', 'x': ['Before', 'After'], 'y': [3, 2], 'name': 'Genes'},
        ],
        'layout': {'title': {'text': 'Style smoke'}, 'width': 800, 'height': 500},
    }
    rendered = render_plotly_payload(payload, str(tmp_path / 'style_smoke'), 'Style smoke')

    assert {item['file_type'] for item in rendered} == {'png', 'svg'}
    assert all(os.path.isfile(item['file_path']) for item in rendered)
    assert all(os.path.getsize(item['file_path']) > 0 for item in rendered)


def test_native_grouped_bars_pool_high_cardinality_series():
    import numpy as np

    from modules.native_figures import grouped_bar_figure

    series = [(f'batch_{i}', np.array([i + 1, 2 * (i + 1)])) for i in range(30)]
    fig = grouped_bar_figure(['cluster_0', 'cluster_1'], series, stacked=True)

    # 11 dominant batches + explicit Other, rather than a legend thousands of
    # pixels tall when the input annotation column is high-cardinality.
    assert len(fig.axes[0].get_legend_handles_labels()[0]) == 12
    assert tuple(round(value, 1) for value in fig.get_size_inches()) == (7.0, 5.0)


def test_native_save_preserves_layout_dimensions(tmp_path):
    import matplotlib.pyplot as plt
    from PIL import Image

    from modules.base import BaseAnalysis

    class Stub(BaseAnalysis):
        MODULE_NAME = 'style_stub'

        def run(self, input_path):
            return {}

    figure, ax = plt.subplots(figsize=(14, 8), dpi=150)
    ax.plot([0, 1], [0, 1])
    result = Stub(str(tmp_path), {}, lambda *_: None).save_matplotlib_figure(
        figure, str(tmp_path), 'wide.png', 'umap', 'Wide figure',
        formats=('png',), dpi=100,
    )
    assert result[0]['file_type'] == 'png'
    assert Image.open(tmp_path / 'wide.png').size[0] >= 1200


def test_native_save_can_preserve_panel_aspect_with_visualization_defaults(tmp_path):
    import matplotlib.pyplot as plt
    from PIL import Image

    from modules.base import BaseAnalysis

    class Stub(BaseAnalysis):
        MODULE_NAME = 'style_stub'

        def run(self, input_path):
            return {}

    figure, _ = plt.subplots(figsize=(13.2, 4.0), dpi=150)
    result = Stub(
        str(tmp_path),
        {'_visualization': {'figure_width': 800, 'figure_height': 500}},
        lambda *_: None,
    ).save_matplotlib_figure(
        figure, str(tmp_path), 'panel.png', 'umap', 'Panel figure',
        formats=('png',), dpi=100, preserve_aspect=True,
    )
    assert result[0]['file_type'] == 'png'
    width, height = Image.open(tmp_path / 'panel.png').size
    assert width / height > 2.5


def test_correlation_heatmap_masks_diagonal_and_uses_readable_scale():
    import matplotlib.pyplot as plt
    import numpy as np

    from modules.native_figures import correlation_heatmap_figure

    corr = np.array([
        [1.0, 0.992, 0.955, 0.951],
        [0.992, 1.0, 0.958, 0.954],
        [0.955, 0.958, 1.0, 0.990],
        [0.951, 0.954, 0.990, 1.0],
    ])
    fig, metadata = correlation_heatmap_figure(
        corr, ['ctrl-1', 'ctrl-2', 'drug-1', 'drug-2'],
        group_labels=['ctrl', 'ctrl', 'drug', 'drug'],
    )

    assert metadata['display_vmin'] == 0.95
    assert metadata['display_vmax'] == 1.0
    assert metadata['off_diagonal']['min'] == 0.951
    assert set(metadata['group_color_map']) == {'ctrl', 'drug'}
    heatmap_ax = next(axis for axis in fig.axes if axis.images and axis.get_xlabel() == 'Sample')
    assert np.all(np.ma.getmaskarray(heatmap_ax.images[0].get_array()).diagonal())
    plt.close(fig)


def test_correlation_pairwise_tables_distinguish_within_and_between_groups():
    import numpy as np

    from modules.native_figures import correlation_pairwise_table, summarize_correlation_pairs

    corr = np.array([
        [1.0, 0.99, 0.90, 0.89],
        [0.99, 1.0, 0.91, 0.90],
        [0.90, 0.91, 1.0, 0.98],
        [0.89, 0.90, 0.98, 1.0],
    ])
    pairs = correlation_pairwise_table(
        corr, ['a1', 'a2', 'b1', 'b2'], ['A', 'A', 'B', 'B'])
    summary = summarize_correlation_pairs(pairs)

    assert len(pairs) == 6
    assert set(pairs['pair_type']) == {'within_group', 'between_group'}
    within = summary.loc[summary['pair_set'] == 'within_group', 'median_r'].iloc[0]
    between = summary.loc[summary['pair_set'] == 'between_group', 'median_r'].iloc[0]
    assert within > between


def test_correlation_heatmap_centres_negative_values_at_zero():
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import TwoSlopeNorm

    from modules.native_figures import correlation_heatmap_figure

    corr = np.array([
        [1.0, -0.30, 0.15],
        [-0.30, 1.0, 0.25],
        [0.15, 0.25, 1.0],
    ])
    fig, metadata = correlation_heatmap_figure(corr, ['a', 'b', 'c'])
    heatmap_ax = next(axis for axis in fig.axes if axis.images and axis.get_xlabel() == 'Sample')

    assert metadata['display_vmin'] < 0
    assert isinstance(heatmap_ax.images[0].norm, TwoSlopeNorm)
    assert heatmap_ax.images[0].norm.vcenter == 0
    plt.close(fig)


def test_correlation_heatmap_thins_sample_ticks_for_36_samples():
    import matplotlib.pyplot as plt
    import numpy as np

    from modules.native_figures import correlation_heatmap_figure

    n_samples = 36
    corr = np.full((n_samples, n_samples), 0.96)
    np.fill_diagonal(corr, 1.0)
    fig, metadata = correlation_heatmap_figure(
        corr, [f'long_sample_name_{index + 1:02d}' for index in range(n_samples)],
    )
    heatmap_ax = next(axis for axis in fig.axes if axis.images and axis.get_xlabel().startswith('Sample'))

    assert metadata['sample_label_stride'] == 2
    assert metadata['sample_labels_shown'] == 19
    assert len(heatmap_ax.get_xticklabels()) == 19
    assert len(heatmap_ax.get_yticklabels()) == 19
    assert 'every 2 samples' in heatmap_ax.get_xlabel()
    plt.close(fig)


def test_bulk_pca_embedding_uses_group_legend_and_boxed_sample_labels():
    import matplotlib.pyplot as plt
    import numpy as np

    from modules.bulk_pca import _pca_embedding_figure

    coords = np.array([
        [-2.0, -1.0], [-1.5, -0.8], [-1.7, -1.3],
        [1.3, 1.0], [1.7, 1.3], [1.5, 0.7],
    ])
    fig = _pca_embedding_figure(
        coords, ['Ctrl'] * 3 + ['Drug'] * 3,
        [f'S{i}' for i in range(6)], 'PCA', 'PC1 (20%)', 'PC2 (18%)',
        group_label='_auto_group', subtitle='Color: _auto_group · PC1 + PC2: 38.0%',
    )
    ax = fig.axes[0]

    assert [text.get_text() for text in ax.get_legend().get_texts()] == ['Ctrl (n=3)', 'Drug (n=3)']
    labels = [text for text in ax.texts if text.get_text().startswith('S')]
    assert len(labels) == 6
    assert all(label.get_bbox_patch() is not None for label in labels)
    plt.close(fig)


def test_bulk_pca_embedding_hides_replicate_labels_in_dense_panels():
    import matplotlib.pyplot as plt
    import numpy as np

    from modules.bulk_pca import _pca_embedding_figure

    coords = np.column_stack([np.arange(9, dtype=float), np.arange(9, dtype=float)])
    fig = _pca_embedding_figure(
        coords, ['Ctrl'] * 4 + ['Drug'] * 5,
        [f'S{i}' for i in range(9)], 'PCA', 'PC1', 'PC2',
        subtitle='Replicate names are intentionally omitted.',
    )
    ax = fig.axes[0]

    assert all(not text.get_text().startswith('S') for text in ax.texts)
    plt.close(fig)
