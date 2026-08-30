import json
from pathlib import Path

import numpy as np


def _embedding_data():
    rng = np.random.default_rng(42)
    coordinates = np.vstack([
        rng.normal(loc=(-2.0, 0.0), scale=0.55, size=(120, 2)),
        rng.normal(loc=(2.0, 0.0), scale=0.55, size=(120, 2)),
    ])
    categories = np.asarray(["Cluster A"] * 120 + ["Cluster B"] * 120)
    return coordinates, categories


def test_single_cell_embedding_has_dedicated_contract_and_rasterized_cell_layer():
    import matplotlib.pyplot as plt
    from figure_engine import NatureFigureDirector

    coordinates, categories = _embedding_data()
    director = NatureFigureDirector()
    spec = director.spec_from_params(
        "umap", {"_figure_mode": "nature_portfolio"},
        title="UMAP by cluster",
    )
    figure = director.render(spec, {
        "coordinates": coordinates,
        "values": categories,
        "value_type": "categorical",
        "value_label": "cluster",
    })

    assert spec.plot_type == "embedding"
    assert figure._nature_spec["plot_type"] == "embedding"
    assert figure._nature_encodings["point"] == "cell"
    assert any(collection.get_rasterized() for collection in figure.axes[0].collections)
    plt.close(figure)


def test_base_embedding_reads_gene_expression_and_uses_exact_export_contract(tmp_path):
    import anndata as ad
    import pandas as pd
    from modules.base import BaseAnalysis

    class DummyAnalysis(BaseAnalysis):
        def run(self, input_path):
            return {}

    coordinates, categories = _embedding_data()
    expression = np.column_stack([
        np.linspace(0.0, 5.0, len(coordinates)),
        np.linspace(3.0, 0.0, len(coordinates)),
    ])
    adata = ad.AnnData(
        X=expression,
        obs=pd.DataFrame({"cluster": pd.Categorical(categories)}),
        var=pd.DataFrame(index=["GENE1", "GENE2"]),
    )
    adata.obsm["X_umap"] = coordinates
    module = DummyAnalysis(
        project_dir=str(tmp_path),
        params={"_figure_mode": "nature_portfolio"},
        progress_callback=lambda *_: None,
    )

    figure = module.build_publication_umap(adata, "GENE1", title="GENE1 expression")
    assert figure._nature_encodings["color"] == "continuous GENE1"
    outputs = module.save_matplotlib_figure(
        figure, module.ensure_plots_dir(), "gene1_umap.png", "umap", "GENE1 UMAP",
    )

    by_type = {item["file_type"]: item for item in outputs}
    assert {"svg", "pdf", "png", "json"}.issubset(by_type)
    assert all(Path(item["file_path"]).is_file() for item in outputs)
    report = json.loads(Path(by_type["json"]["file_path"]).read_text(encoding="utf-8"))
    assert report["nature_readiness_score"] >= 90
    assert not any(issue["severity"] == "error" for issue in report["issues"])


def test_embedding_continuous_scale_records_percentile_clipping_warning():
    import matplotlib.pyplot as plt
    from figure_engine import NatureFigureDirector

    coordinates, _ = _embedding_data()
    values = np.linspace(0.0, 1.0, len(coordinates))
    values[-1] = 1000.0
    director = NatureFigureDirector()
    spec = director.spec_from_params("embedding", {}, title="QC burden")
    figure = director.render(spec, {
        "coordinates": coordinates,
        "values": values,
        "value_type": "continuous",
        "value_label": "MT%",
    })

    assert figure._nature_encodings["color"] == "continuous MT%"
    assert any("1–99" in warning for warning in figure._nature_semantic_warnings)
    plt.close(figure)


def test_dense_celltype_labels_do_not_overlap_at_final_canvas_size():
    import matplotlib.pyplot as plt
    from figure_engine import FigureValidator, NatureFigureDirector

    rng = np.random.default_rng(2026)
    coordinates = []
    labels = []
    for index in range(20):
        center = rng.normal(0.0, 0.03, size=2)
        coordinates.append(rng.normal(center, 0.005, size=(20, 2)))
        labels.extend([f'Cell type {index} with a very long biological label'] * 20)

    director = NatureFigureDirector()
    spec = director.spec_from_params('embedding', {})
    figure = director.render(spec, {
        'coordinates': np.vstack(coordinates),
        'values': np.asarray(labels),
        'value_type': 'categorical',
        'value_label': 'celltype',
    })
    try:
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        visible = [text for text in figure.axes[0].texts if text.get_visible()]
        overlaps = sum(
            left.get_window_extent(renderer).overlaps(right.get_window_extent(renderer))
            for index, left in enumerate(visible)
            for right in visible[index + 1:]
        )

        assert overlaps == 0
        report = FigureValidator().validate(figure, spec)
        assert report.metrics['label_overlap_count'] == 0
    finally:
        plt.close(figure)


def test_single_column_volcano_wraps_long_pseudobulk_title_inside_canvas():
    import matplotlib.pyplot as plt
    import pandas as pd
    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame({
        "gene": [f"GENE{i}" for i in range(40)],
        "log2FC": np.linspace(-2.0, 2.0, 40),
        "padj": np.geomspace(1e-5, 0.8, 40),
    })
    director = NatureFigureDirector()
    spec = director.spec_from_params(
        "volcano", {"_figure_mode": "nature_portfolio"}, width="single",
        title=(
            "Pseudobulk DEG · Control vs Long_treatment_condition · "
            "cluster 0"
        ),
    )
    figure = director.render(spec, frame)
    figure.canvas.draw()
    title = figure.axes[0]._left_title
    title_box = title.get_window_extent(figure.canvas.get_renderer())

    assert "\n" in title.get_text()
    assert title_box.x0 >= figure.bbox.x0
    assert title_box.x1 <= figure.bbox.x1
    assert title_box.y1 <= figure.bbox.y1
    plt.close(figure)
