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
