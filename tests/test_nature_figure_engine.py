"""Phase 1 regression tests for the deterministic Nature Figure Engine."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def test_figure_spec_uses_physical_submission_widths_and_scientific_controls():
    from figure_engine import FigureSpec

    single = FigureSpec(plot_type="volcano", width="single", fc_threshold=1.25)
    double = FigureSpec(plot_type="heatmap", width="double")

    assert single.width_mm == 89.0
    assert double.width_mm == 183.0
    assert single.fc_threshold == 1.25
    assert FigureSpec.from_mapping(single.to_dict()) == single
    with pytest.raises(ValueError, match="fdr_threshold"):
        FigureSpec(plot_type="volcano", fdr_threshold=0)


def test_director_rejects_arbitrary_visual_options_and_maps_legacy_params():
    from figure_engine import NatureFigureDirector

    director = NatureFigureDirector()
    spec = director.spec_from_params(
        "volcano",
        {
            "_visualization": {"theme": "nature"},
            "fc_threshold": 1.0,
            "fdr_threshold": 0.01,
            "label_genes": "TP53, MYC",
        },
    )

    assert spec.mode == "nature_portfolio"
    assert spec.label_genes == ("TP53", "MYC")
    with pytest.raises(TypeError, match="Unknown FigureSpec"):
        director.create_spec("pca", arbitrary_font_size=18)


def test_nature_pca_encodes_group_and_batch_without_default_labels():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    director = NatureFigureDirector()
    spec = director.create_spec("pca", show_legend=True)
    figure = director.render(
        spec,
        {
            "coordinates": np.array([
                [-2.1, 0.1], [-1.8, -0.2], [-2.0, 0.3],
                [1.8, 0.0], [2.2, 0.2], [2.0, -0.3],
            ]),
            "groups": ["Control"] * 3 + ["Treatment"] * 3,
            "batches": ["B1", "B2", "B1", "B1", "B2", "B1"],
            "samples": [f"sample_{index}" for index in range(6)],
            "explained_variance": [0.52, 0.21],
        },
    )
    axis = figure.axes[0]

    assert axis.get_xlabel() == "PC1 (52.0%)"
    assert axis.get_ylabel() == "PC2 (21.0%)"
    assert not axis.xaxis.get_gridlines()[0].get_visible()
    assert not any(text.get_text().startswith("sample_") for text in axis.texts)
    assert figure._nature_encodings == {"color": "group", "marker": "batch"}
    plt.close(figure)


def test_nature_volcano_applies_thresholds_and_controlled_gene_labels():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame({
        "gene": ["TP53", "MYC", "STAT1", "ACTB", "GAPDH", "FOXA2"],
        "log2FC": [2.4, -2.1, 1.6, 0.1, -0.2, -1.8],
        "padj": [1e-7, 2e-6, 0.002, 0.8, 0.7, 0.004],
    })
    director = NatureFigureDirector()
    spec = director.create_spec(
        "volcano", fc_threshold=1.0, fdr_threshold=0.05,
        label_n=3, label_genes=("TP53",),
    )
    figure = director.render(spec, frame)
    axis = figure.axes[0]
    labels = {text.get_text() for text in axis.texts}

    assert "TP53" in labels
    assert len([line for line in axis.lines if line.get_linestyle() != "None"]) == 3
    assert axis.spines["top"].get_visible() is False
    assert axis.spines["right"].get_visible() is False
    assert not any(line.get_visible() for line in axis.xaxis.get_gridlines())
    plt.close(figure)


def test_nature_volcano_auto_labels_balance_directions_and_skip_bare_feature_ids():
    from figure_engine import NatureFigureDirector
    from figure_engine.templates.volcano import _select_labels

    frame = pd.DataFrame({
        'gene': ['NFKBIA', 'BHLHE40', 'ENSG00000261600', 'IFI44L', 'MX1', 'ENSMUSG000001'],
        'log2FC': [1.0, 0.8, 2.2, -0.9, -0.7, -2.0],
        'padj': [1e-12, 1e-8, 1e-15, 1e-14, 1e-9, 1e-16],
    })
    frame['_gene'] = frame['gene']
    frame['_log2fc'] = frame['log2FC']
    frame['_fdr'] = frame['padj']
    regulation = pd.Series(['Up', 'Up', 'Up', 'Down', 'Down', 'Down'])

    selected = _select_labels(
        frame, NatureFigureDirector().create_spec('volcano', label_n=4), regulation,
    )

    assert frame.loc[selected, 'gene'].tolist() == ['NFKBIA', 'BHLHE40', 'IFI44L', 'MX1']


def test_nature_volcano_caps_extreme_fdr_tail_with_top_triangles():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame({
        'gene': [f'GENE{index}' for index in range(6)],
        'log2FC': [-2.4, -1.3, -0.4, 0.4, 1.3, 2.4],
        'padj': [0.2, 1e-4, 1e-12, 1e-30, 1e-90, 1e-250],
    })
    figure = NatureFigureDirector().render(
        NatureFigureDirector().create_spec('volcano', label_n=0), frame,
    )
    axis = figure.axes[0]
    points = np.vstack([
        collection.get_offsets()
        for collection in axis.collections
        if len(collection.get_offsets())
    ])
    points = points[np.argsort(points[:, 0])]

    cap = axis.get_ylim()[1] * 0.985
    # The display remains in raw −log10(FDR) units; extreme FDR values are
    # retained as top-edge triangles instead of a nonlinear tail transform.
    assert axis.get_ylim()[1] == 18
    assert np.sum(np.isclose(points[:, 1], cap)) == 3
    assert 'compressed tail' not in axis.get_ylabel()
    assert figure._nature_encodings['y'].startswith('-log10(FDR); top triangles')
    plt.close(figure)


def test_nature_volcano_keeps_extreme_points_and_labels_inside_axes():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    rng = np.random.default_rng(11)
    log2fc = np.concatenate([
        rng.normal(-9.5, 1.2, 30),
        rng.normal(0.0, 0.8, 200),
        rng.normal(8.2, 1.0, 30),
    ])
    padj = np.concatenate([
        rng.uniform(1e-8, 1e-3, 30),
        rng.uniform(0.1, 0.9, 200),
        rng.uniform(1e-8, 1e-3, 30),
    ])
    frame = pd.DataFrame({
        'gene': [f'GENE{index}' for index in range(len(log2fc))],
        'log2FC': log2fc,
        'padj': padj,
    })
    figure = NatureFigureDirector().render(
        NatureFigureDirector().create_spec(
            'volcano', fc_threshold=1.0, fdr_threshold=0.05, label_n=8,
        ),
        frame,
    )
    axis = figure.axes[0]
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    axes_box = axis.get_window_extent(renderer)

    # The x-axis stays compact, and extreme effects remain visibly encoded at
    # both edges rather than stretching the sparse central cloud.
    points = np.vstack([
        collection.get_offsets()
        for collection in axis.collections
        if len(collection.get_offsets())
    ])
    assert axis.get_xlim() == (-2.5, 2.5)
    assert np.sum(np.isclose(np.abs(points[:, 0]), 2.5 * 0.985)) >= 60

    # Gene labels must stay inside the plotting area on both sides (1 px
    # tolerance for renderer rounding).
    from matplotlib.text import Text
    for text in [item for item in axis.texts if item.get_visible()]:
        box = Text.get_window_extent(text, renderer)
        assert box.x0 >= axes_box.x0 - 1
        assert box.x1 <= axes_box.x1 + 1
        assert box.y0 >= axes_box.y0 - 1
        assert box.y1 <= axes_box.y1 + 1
    plt.close(figure)


def test_nature_volcano_separates_dense_labels_with_leader_lines_and_legend():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame({
        'gene': ['UP_A', 'UP_B', 'UP_C', 'DOWN_A', 'DOWN_B', 'DOWN_C'],
        'log2FC': [1.05, 1.08, 1.12, -1.05, -1.08, -1.12],
        'padj': [1e-12, 2e-12, 3e-12, 1e-12, 2e-12, 3e-12],
    })
    figure = NatureFigureDirector().render(
        NatureFigureDirector().create_spec(
            'volcano', fc_threshold=1.0, label_n=6, show_legend=True,
        ), frame,
    )
    axis = figure.axes[0]
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    axes_box = axis.get_window_extent(renderer)
    from matplotlib.text import Text
    labels = [
        text for text in axis.texts
        if text.get_visible() and text.get_text().startswith(('UP_', 'DOWN_'))
    ]

    # A very dense cluster can safely retain four to six of the requested
    # labels after the final collision guard.  The full list remains in CSV.
    assert 4 <= len(labels) <= 6
    assert all(label.arrow_patch is not None for label in labels)
    for index, label in enumerate(labels):
        assert label.get_bbox_patch() is None
        box = Text.get_window_extent(label, renderer)
        assert box.x0 >= axes_box.x0 - 1
        assert box.x1 <= axes_box.x1 + 1
        assert box.y0 >= axes_box.y0 - 1
        assert box.y1 <= axes_box.y1 + 1
        assert not any(
                box.overlaps(Text.get_window_extent(other, renderer))
            for other in labels[index + 1:]
        )
    assert axis.get_legend() is not None
    assert axis.get_legend().get_title().get_text() == ''
    assert {text.get_text() for text in axis.get_legend().get_texts()} == {
        'Up  3', 'Down  3', 'NS  0',
    }
    plt.close(figure)


def test_multi_volcano_promotes_proxy_legends_to_one_shared_key():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureComposer, NatureFigureDirector
    from figure_engine.composer import FigurePanel

    frame = pd.DataFrame({
        'gene': ['UP', 'DOWN', 'NS'],
        'log2FC': [1.5, -1.5, 0.1],
        'padj': [0.001, 0.002, 0.8],
    })
    director = NatureFigureDirector()
    outer_spec = director.create_spec('diagnostic', width='double').with_updates(
        plot_type='composite', height_mm=82.0,
    )
    panel_spec = director.create_spec(
        'volcano', width='single', label_n=0, show_legend=True,
        extra={'legend_counts': False},
    )
    panels = [
        FigurePanel(plot_type='volcano', data=frame, spec=panel_spec,
                    label='a', row=0, column=0),
        FigurePanel(plot_type='volcano', data=frame, spec=panel_spec,
                    label='b', row=0, column=1),
    ]
    figure = NatureFigureComposer(director).compose(
        panels, outer_spec, nrows=1, ncols=2, shared_legend=True,
    )

    assert len(figure.legends) == 1
    assert {text.get_text() for text in figure.legends[0].get_texts()} == {
        'Up', 'Down', 'NS',
    }
    assert all(axis.get_legend() is None for axis in figure.axes)
    plt.close(figure)


def test_nature_ma_keeps_extreme_fold_changes_uncompressed_and_labels_requested_genes():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame({
        'gene': ['EDGE_DOWN', 'TARGET_NS', 'EDGE_UP', 'MID_UP', 'MID_DOWN'],
        'mean_expression': [8, 20, 64, 180, 350],
        'log2FC': [-25.0, 0.15, 24.0, 2.2, -2.1],
        'padj': [1e-8, 0.9, 1e-9, 0.002, 0.003],
    })
    figure = NatureFigureDirector().render(
        NatureFigureDirector().create_spec(
            'ma', fc_threshold=1.0, fdr_threshold=0.05,
            label_n=1, label_genes=('TARGET_NS',),
        ), frame,
    )
    axis = figure.axes[0]
    points = np.vstack([
        collection.get_offsets()
        for collection in axis.collections
        if len(collection.get_offsets())
    ])
    y_min, y_max = axis.get_ylim()
    assert np.isclose(points[:, 1].min(), -25.0)
    assert np.isclose(points[:, 1].max(), 24.0)
    assert points[:, 1].min() > y_min
    assert points[:, 1].max() < y_max
    assert {text.get_text() for text in axis.texts if text.get_visible()} == {'TARGET_NS'}
    plt.close(figure)


def test_nature_heatmap_integrates_annotations_and_controls_label_density():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    rng = np.random.default_rng(7)
    matrix = rng.normal(size=(18, 12))
    matrix[:6, 6:] += 1.8
    director = NatureFigureDirector()
    spec = director.create_spec(
        "heatmap", max_row_labels=7, max_col_labels=6,
    )
    figure = director.render(
        spec,
        {
            "matrix": matrix,
            "gene_labels": [f"Gene{index:02d}" for index in range(18)],
            "sample_labels": [f"Sample{index:02d}" for index in range(12)],
            "annotations": {
                "Group": ["Control"] * 6 + ["Treatment"] * 6,
                "Batch": ["B1", "B2"] * 6,
            },
        },
    )
    heat_axis = max((axis for axis in figure.axes if axis.images),
                    key=lambda axis: np.asarray(axis.images[0].get_array()).size)
    annotation_axis = min((axis for axis in figure.axes if axis.images),
                          key=lambda axis: np.asarray(axis.images[0].get_array()).size)

    values = np.asarray(heat_axis.images[0].get_array())
    assert abs(float(np.nanmean(values))) < 1e-8
    assert len(heat_axis.get_yticklabels()) <= 7
    assert len(heat_axis.get_xticklabels()) <= 6
    assert {text.get_text() for text in annotation_axis.get_yticklabels()} == {"Group", "Batch"}
    assert "gene-wise z-score" in figure._nature_encodings["color"]
    plt.close(figure)


def test_enrichment_and_gsea_use_fixed_dotplot_encodings():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    director = NatureFigureDirector()
    ora = pd.DataFrame({
        "Term": ["Cellular response to lipid", "Mitochondrial organization", "Inflammation"],
        "FDR": [0.002, 0.01, 0.04],
        "GeneRatio": ["8/120", "5/120", "3/120"],
        "Count": [8, 5, 3],
    })
    ora_spec = director.create_spec("enrichment", top_n=3)
    ora_figure = director.render(ora_spec, ora)
    assert ora_figure.axes[0].get_xlabel() == "Gene ratio"
    assert ora_figure._nature_encodings == {
        "x": "Gene ratio", "size": "Count", "color": "FDR",
    }

    gsea = pd.DataFrame({
        "Term": ["Fatty acid metabolism", "Interferon response", "Cell cycle"],
        "NES": [2.1, -1.8, 1.4],
        "FDR": [0.004, 0.009, 0.03],
        "setSize": [86, 72, 49],
    })
    gsea_spec = director.create_spec("gsea", top_n=3)
    gsea_figure = director.render(gsea_spec, gsea)
    assert "NES" in gsea_figure.axes[0].get_xlabel()
    assert any(np.allclose(line.get_xdata(), [0, 0]) for line in gsea_figure.axes[0].lines)
    plt.close(ora_figure)
    plt.close(gsea_figure)


def test_bulk_enrichment_adapter_prefers_dotplot_for_complete_ora_tables():
    import matplotlib.pyplot as plt

    from modules.bulk_enrichment import _enrichment_figure

    frame = pd.DataFrame({
        "Term": ["Lipid transport", "Fatty acid metabolism"],
        "Adjusted P-value": [0.002, 0.01],
        "Overlap": ["5/100", "3/100"],
        "Count": [5, 3],
    })
    figure = _enrichment_figure(frame, "GO enrichment", database="GO_BP")

    assert getattr(figure, "_nature_spec_object").plot_type == "enrichment_dotplot"
    assert figure.axes[0].collections
    assert not figure.axes[0].patches
    plt.close(figure)


def test_export_and_validator_preserve_editable_89mm_outputs(tmp_path):
    import matplotlib.pyplot as plt
    from PIL import Image

    from figure_engine import NatureFigureDirector, export_figure

    rng = np.random.default_rng(13)
    coordinates = np.vstack([
        rng.normal((-1.8, 0), 0.25, size=(6, 2)),
        rng.normal((1.8, 0), 0.25, size=(6, 2)),
    ])
    director = NatureFigureDirector()
    spec = director.create_spec(
        "pca", width="single", show_legend=True,
        formats=("svg", "pdf", "png"),
    )
    figure = director.render(
        spec,
        {
            "coordinates": coordinates,
            "groups": ["Control"] * 6 + ["Treatment"] * 6,
            "explained_variance": [0.48, 0.24],
        },
    )
    paths, report = export_figure(figure, tmp_path / "pca", spec)

    assert set(paths) == {"svg", "pdf", "png"}
    assert all(Path(path).stat().st_size > 0 for path in paths.values())
    assert "<text" in Path(paths["svg"]).read_text(encoding="utf-8")
    assert Path(paths["pdf"]).read_bytes().startswith(b"%PDF-")
    with Image.open(paths["png"]) as image:
        assert image.width == pytest.approx(89 / 25.4 * 600, abs=2)
    assert report.ready is True
    assert report.score >= 90
    assert report.metrics["canvas_mm"][0] == pytest.approx(89.0, abs=0.01)
    plt.close(figure)


def test_nature_validator_fails_on_one_real_text_overlap(tmp_path):
    import matplotlib.pyplot as plt

    from figure_engine import FigureSpec, export_figure

    fig, ax = plt.subplots(figsize=(89 / 25.4, 55 / 25.4), dpi=600)
    ax.plot([0, 1], [0, 1], color='#355F8A')
    ax.text(.5, .5, 'first label', transform=ax.transAxes, fontsize=7)
    ax.text(.5, .5, 'second label', transform=ax.transAxes, fontsize=7)
    spec = FigureSpec(plot_type='pca', width='single', height_mm=55,
                      formats=('svg', 'pdf', 'png'))
    _, report = export_figure(fig, tmp_path / 'overlap', spec)

    assert report.metrics['label_overlap_count'] >= 1
    assert report.ready is False
    assert any(issue.code == 'label_overlap' and issue.severity == 'error'
               for issue in report.issues)
    plt.close(fig)
