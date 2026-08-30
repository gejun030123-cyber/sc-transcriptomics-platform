import os
import sys
import types

import numpy as np
import pandas as pd
import pytest


def _make_10x_dir(path):
    path.mkdir(parents=True, exist_ok=True)
    for name in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
        (path / name).write_text("", encoding="utf-8")


def _tiny_batch_adata(tmp_path):
    import anndata as ad

    samples = ["Ctrl_1", "Ctrl_2", "Trt_1", "Trt_2"]
    obs = []
    rows = []
    for sample in samples:
        condition = "Control" if sample.startswith("Ctrl") else "Treatment"
        for cell in range(2):
            obs.append({
                "sample_id": sample, "condition": condition,
                "celltype": "TypeA" if cell == 0 else "TypeB",
                "leiden": str(cell),
                "replicate": sample.rsplit("_", 1)[1],
            })
            rows.append([100 if condition == "Treatment" else 5, 12, 3, 0])
    adata = ad.AnnData(
        X=np.asarray(rows, dtype=np.float32),
        obs=pd.DataFrame(obs, index=[f"cell_{index}" for index in range(len(obs))]),
        var=pd.DataFrame({"gene_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]},
                         index=["G1", "G2", "G3", "G4"]),
    )
    adata.layers["counts"] = np.asarray(rows, dtype=np.int64)
    input_path = tmp_path / "input.h5ad"
    adata.write_h5ad(input_path)
    return input_path


def test_pseudobulk_count_validation_checks_entire_matrix():
    import anndata as ad
    from modules.sc_batch import aggregate_pseudobulk_counts

    counts = np.ones((2, 100_001), dtype=np.float32)
    counts[-1, -1] = 0.5
    adata = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame({
            'sample_id': ['S1', 'S2'], 'condition': ['A', 'B'],
        }),
        var=pd.DataFrame(index=[f'G{i}' for i in range(counts.shape[1])]),
    )
    adata.layers['counts'] = counts

    with pytest.raises(ValueError, match='非负整数'):
        aggregate_pseudobulk_counts(adata)


def test_pseudobulk_enrichment_rejects_online_mode_before_sending_genes(tmp_path):
    from modules.sc_cell_go import SCCellGOEnrichment

    module = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={'source_level': 'pseudobulk', 'execution_mode': 'enrichr'},
        progress_callback=lambda *_: None,
    )
    with pytest.raises(ValueError, match='tested-gene background'):
        module.run(str(tmp_path / 'unused.h5ad'))


def test_local_gene_set_paths_reject_traversal_and_resolve_project_relative(tmp_path):
    from modules.sc_cell_go import (
        _local_gene_set_path, _requested_gene_sets, _resolve_local_gene_set_dir,
    )

    local = tmp_path / 'resources' / 'go'
    local.mkdir(parents=True)
    (local / 'Safe.gmt').write_text('Term\tID\tGENE\n', encoding='utf-8')

    assert _resolve_local_gene_set_dir('resources/go', str(tmp_path)) == str(local.resolve())
    assert _local_gene_set_path('Safe', local) == (local / 'Safe.gmt').resolve()
    with pytest.raises(ValueError, match='名称'):
        _requested_gene_sets({'gene_sets': '../secret'})

    outside = tmp_path / 'outside.gmt'
    outside.write_text('Term\tID\tGENE\n', encoding='utf-8')
    link = local / 'Linked.gmt'
    try:
        link.symlink_to(outside)
    except OSError:
        pytest.skip('当前文件系统不支持符号链接测试')
    with pytest.raises(ValueError, match='符号链接'):
        _local_gene_set_path('Linked', local)


def test_discovery_manifest_does_not_guess_biology(tmp_path, monkeypatch):
    from config import Config
    from modules.sc_batch import discovered_manifest_frame, read_sample_manifest

    root = tmp_path / "server_data"
    _make_10x_dir(root / "nested" / "Control_rep1")
    _make_10x_dir(root / "nested" / "Treatment_rep1")
    monkeypatch.setenv("SC_BATCH_SOURCE_ROOTS", str(root))

    manifest = discovered_manifest_frame(str(root))
    assert len(manifest) == 2
    assert manifest["condition"].eq("").all()
    assert manifest["replicate"].eq("").all()

    manifest["condition"] = ["Control", "Treatment"]
    manifest["replicate"] = ["1", "1"]
    path = tmp_path / "manifest.csv"
    manifest.to_csv(path, index=False)
    parsed = read_sample_manifest(str(path), str(root))
    assert all(os.path.isabs(item) for item in parsed["matrix_path"])
    assert tuple(Config.sc_batch_source_roots()) == (os.path.realpath(root),)


def test_batch_csv_export_writes_standard_contract(tmp_path):
    from modules.sc_batch_export import CSV_PACKAGE_FOLDERS, SCBatchCSVExport

    input_path = _tiny_batch_adata(tmp_path)
    module = SCBatchCSVExport(
        project_dir=str(tmp_path),
        params={"sample_key": "sample_id", "condition_key": "condition"},
        progress_callback=lambda *_: None,
    )
    result = module.run(str(input_path))
    labels = {item["label"] for item in result["result_files"]}
    assert "Biological sample design" in labels
    assert "Sample-level pseudobulk raw counts (gene × sample)" in labels
    assert "Sample-level cluster proportions" in labels
    assert "Pseudobulk comparison registry" in labels
    assert "GO enrichment status by comparison" in labels
    assert "DOCX-compatible group proportions (wide)" in labels
    comparison_file = next(item["file_path"] for item in result["result_files"]
                           if item["label"] == "Pseudobulk comparison registry")
    plan = pd.read_csv(comparison_file)
    # all_pairwise 方向约定：control 风格标签固定放在对照（右侧），
    # 保证 log2FC 正值 = 实验组（Treatment）上调。
    assert plan.loc[0, "comparison"] == "Treatment vs Control"
    assert plan.loc[0, "group_1"] == "Treatment"
    assert plan.loc[0, "group_2"] == "Control"
    assert plan.loc[0, "direction_basis"] == "control_name_heuristic"
    assert plan.loc[0, "status"] == "runnable"
    proportions_file = next(item["file_path"] for item in result["result_files"]
                            if item["label"] == "DOCX-compatible group proportions (wide)")
    proportions = pd.read_csv(proportions_file)
    assert {"SampleID", "SampleGroup", "Cluster_0", "Cluster_1"}.issubset(proportions.columns)
    package_root = tmp_path / "results" / "sc_batch_results"
    assert {path.name for path in package_root.iterdir() if path.is_dir()} == set(CSV_PACKAGE_FOLDERS.values())


def test_batch_csv_export_can_write_only_first_two_folders(tmp_path):
    from modules.sc_batch_export import CSV_PACKAGE_FOLDERS, SCBatchCSVExport

    input_path = _tiny_batch_adata(tmp_path)
    result = SCBatchCSVExport(
        project_dir=str(tmp_path),
        params={
            "sample_key": "sample_id", "condition_key": "condition",
            "export_folder": "first_two_only", "export_scope": "proportions_expression",
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))
    labels = {item["label"] for item in result["result_files"]}
    assert "Pseudobulk comparison registry" not in labels
    assert "GO enrichment status by comparison" not in labels
    package_root = tmp_path / "results" / "first_two_only"
    assert {path.name for path in package_root.iterdir() if path.is_dir()} == {
        CSV_PACKAGE_FOLDERS["proportions"], CSV_PACKAGE_FOLDERS["expression"],
    }
    assert result["summary"]["export_scope"] == "proportions_expression"


def test_batch_csv_export_includes_umap_coordinates_in_cell_metadata(tmp_path):
    import anndata as ad
    from modules.sc_batch_export import SCBatchCSVExport

    input_path = _tiny_batch_adata(tmp_path)
    adata = ad.read_h5ad(input_path)
    adata.obsm["X_umap"] = np.column_stack((
        np.arange(adata.n_obs, dtype=float),
        np.arange(adata.n_obs, dtype=float) + 0.5,
    ))
    adata.write_h5ad(input_path)
    result = SCBatchCSVExport(
        project_dir=str(tmp_path),
        params={"sample_key": "sample_id", "condition_key": "condition"},
        progress_callback=lambda *_: None,
    ).run(str(input_path))
    cell_file = next(item["file_path"] for item in result["result_files"]
                     if item["label"] == "Single-cell metadata with sample provenance")
    exported = pd.read_csv(cell_file)
    assert {"UMAP_1", "UMAP_2"}.issubset(exported.columns)
    assert exported.loc[0, "UMAP_1"] == 0
    assert exported.loc[0, "UMAP_2"] == 0.5


def test_pseudobulk_welch_exports_one_file_per_comparison(tmp_path):
    from modules.sc_batch_export import SCPseudobulkDEG

    input_path = _tiny_batch_adata(tmp_path)
    module = SCPseudobulkDEG(
        project_dir=str(tmp_path),
        params={
            "sample_key": "sample_id", "condition_key": "condition",
            "analysis_scope": "all_cells", "method": "welch_logcpm",
            "min_samples_per_group": 2, "min_cells_per_sample_celltype": 1,
            "export_prefix": "demo", "export_full_tables": True,
        },
        progress_callback=lambda *_: None,
    )
    result = module.run(str(input_path))
    per_comparison = [item for item in result["result_files"] if item["label"].startswith("Pseudobulk DEG:")]
    assert len(per_comparison) == 1
    assert "03_differential_expression" in per_comparison[0]["file_path"]
    table = pd.read_csv(per_comparison[0]["file_path"])
    assert {"comparison_id", "gene", "log2FC", "padj", "method"}.issubset(table.columns)
    assert table["method"].eq("welch_t_log2cpm").all()


def test_pseudobulk_per_cluster_figures_cover_all_completed_clusters_by_default(tmp_path, monkeypatch):
    """Per-cluster pseudobulk plots must not silently stop at the old top-eight cap."""
    import anndata as ad
    import figure_engine
    from modules.sc_batch_export import SCPseudobulkDEG

    samples = ["Ctrl_1", "Ctrl_2", "Trt_1", "Trt_2"]
    obs, rows = [], []
    for sample in samples:
        condition = "Control" if sample.startswith("Ctrl") else "Treatment"
        for cluster in range(9):
            obs.append({
                "sample_id": sample,
                "condition": condition,
                "leiden": str(cluster),
            })
            rows.append([
                10 + cluster + (20 if condition == "Treatment" else 0),
                5 + cluster,
            ])
    adata = ad.AnnData(
        X=np.asarray(rows, dtype=np.float32),
        obs=pd.DataFrame(obs, index=[f"cell_{index}" for index in range(len(obs))]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    adata.layers["counts"] = np.asarray(rows, dtype=np.int64)
    input_path = tmp_path / "nine_clusters.h5ad"
    adata.write_h5ad(input_path)

    heatmap_shapes = []
    figure_spec_calls = []

    class FakeNatureFigureDirector:
        def spec_from_params(self, plot_type, *_args, **_kwargs):
            figure_spec_calls.append((plot_type, _kwargs))
            return plot_type

        def render(self, spec, frame):
            if spec == "heatmap":
                heatmap_shapes.append(np.asarray(frame["matrix"]).shape)
            return object()

    saved = []

    def fake_save(self, _fig, _plots_dir, filename, category, label, **_kwargs):
        saved.append((filename, category, label))
        return []

    monkeypatch.setattr(figure_engine, "NatureFigureDirector", FakeNatureFigureDirector)
    monkeypatch.setattr(SCPseudobulkDEG, "save_matplotlib_figure", fake_save)

    result = SCPseudobulkDEG(
        project_dir=str(tmp_path),
        params={
            "sample_key": "sample_id", "condition_key": "condition",
            "analysis_scope": "per_cluster", "cluster_key": "leiden",
            "method": "welch_logcpm", "min_samples_per_group": 2,
            "min_cells_per_sample_celltype": 1,
            "label_genes": "G1; G0",
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))

    # Every completed pseudobulk unit now carries donor-level evidence in
    # addition to the optional DEG pair: PCA, correlation, library support
    # and Top DEG sample expression make a single-donor signal visible.
    assert len(saved) == 54  # Nine clusters × (PCA, correlation, library, heatmap, Volcano, MA).
    assert {category for _, category, _ in saved} == {
        "pca", "correlation", "qc", "heatmap", "volcano", "ma",
    }
    assert sum(category == "pca" for _, category, _ in saved) == 9
    assert sum(category == "correlation" for _, category, _ in saved) == 9
    assert sum(category == "qc" for _, category, _ in saved) == 9
    assert sum(category == "heatmap" for _, category, _ in saved) == 9
    assert heatmap_shapes == [(2, 4)] * 9  # genes × samples, never samples × genes.
    for plot_type in ("volcano", "ma"):
        calls = [kwargs for kind, kwargs in figure_spec_calls if kind == plot_type]
        assert len(calls) == 9
        assert all(call["label_genes"] == ("G1", "G0") for call in calls)
    assert result["summary"]["n_pseudobulk_qc_units"] == 9
    assert result["summary"]["pseudobulk_qc_mandatory"] is True
    assert not any("仅展示显著 DEG 最多的前" in warning for warning in result["summary"]["warnings"])


def test_pseudobulk_exports_gene_csv_by_default(tmp_path):
    """The analysis form defaults to a task-visible, full tested-gene table."""
    from modules.sc_batch_export import SCPseudobulkDEG

    input_path = _tiny_batch_adata(tmp_path)
    result = SCPseudobulkDEG(
        project_dir=str(tmp_path),
        params={
            "sample_key": "sample_id", "condition_key": "condition",
            "analysis_scope": "all_cells", "method": "welch_logcpm",
            "min_samples_per_group": 2, "min_cells_per_sample_celltype": 1,
            "show_deg_figures": False,
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))

    output = next(item for item in result["result_files"]
                  if item["label"].startswith("Pseudobulk DEG:"))
    table = pd.read_csv(output["file_path"])
    assert len(table) == 4
    assert {"gene", "log2FC", "pvalue", "padj", "significant"}.issubset(table.columns)
    assert result["summary"]["deg_source_file"] != output["file_path"]
    # Disabling Volcano/MA does not suppress sample-level pseudobulk QC.
    qc_output = next(item for item in result["result_files"]
                     if item["label"] == "Pseudobulk sample QC: all cells")
    qc_table = pd.read_csv(qc_output["file_path"])
    assert {
        "sample_id", "condition", "total_raw_counts", "n_cells", "deseq2_size_factor",
        "expression_transform",
    }.issubset(qc_table.columns)
    assert not any(item["label"].startswith("Pseudobulk Volcano:")
                   for item in result["result_files"])
    assert any(item["label"].startswith("Pseudobulk sample PCA:")
               for item in result["result_files"])
    assert any(item["label"].startswith("Pseudobulk sample correlation:")
               for item in result["result_files"])
    assert any(item["label"].startswith("Pseudobulk library size:")
               for item in result["result_files"])
    assert any(item["label"].startswith("Pseudobulk Top DEG expression:")
               for item in result["result_files"])


def test_pseudobulk_enrichment_uses_full_background_and_directional_ora(tmp_path):
    from modules.sc_batch import _available_path
    from modules.sc_batch_export import batch_csv_package_dirs
    from modules.sc_cell_go import SCCellGOEnrichment

    input_path = _tiny_batch_adata(tmp_path)
    _, package_dirs = batch_csv_package_dirs(str(tmp_path), "pseudobulk_go_package")
    deg = pd.DataFrame({
        "comparison_id": ["Treatment_vs_Control"] * 6,
        "comparison": ["Treatment vs Control"] * 6,
        "experimental_group": ["Treatment"] * 6,
        "control_group": ["Control"] * 6,
        "deg_scope": ["all_cells"] * 6,
        "cluster": ["All"] * 6,
        "gene": [f"GENE{index}" for index in range(6)],
        "log2FC": [2.0, 2.0, 2.0, -2.0, -2.0, -2.0],
        "statistic": [4.0, 3.0, 2.5, -4.0, -3.0, -2.5],
        "padj": [0.01] * 6,
        "method": ["pydeseq2"] * 6,
    })
    deg.to_csv(_available_path(
        package_dirs["deg"], "sc_pseudobulk_deg_Treatment_vs_Control", ".csv"
    ), index=False)
    local_dir = tmp_path / "go_gene_sets"
    local_dir.mkdir()
    (local_dir / "GO_Biological_Process_2023.gmt").write_text(
        "Up term\tGO:1\tGENE0\tGENE1\tGENE2\n"
        "Down term\tGO:2\tGENE3\tGENE4\tGENE5\n", encoding="utf-8"
    )

    result = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={
            "source_level": "pseudobulk", "deg_prefix": "sc_pseudobulk",
            "export_folder": "pseudobulk_go_package", "execution_mode": "local",
            "direction_mode": "up_down", "min_genes": 3,
            "ora_min_size": 1,
            "local_gene_set_dir": str(local_dir), "show_enrichment_plots": False,
            "export_full_tables": True,
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))
    output = next(item["file_path"] for item in result["result_files"]
                   if item["label"].startswith("Pseudobulk ORA enrichment ("))
    table = pd.read_csv(output)
    assert set(table["direction"]) == {"Up", "Down"}
    assert table["analysis_level"].eq("sample_level_pseudobulk").all()
    assert table["n_background_genes"].eq(6).all()
    assert set(table["Term"]) == {"Up term", "Down term"}


def test_pseudobulk_per_cluster_enrichment_plots_cover_all_significant_clusters_by_default(
        tmp_path, monkeypatch):
    """Every completed, significant pseudobulk cluster gets its enrichment plot."""
    import modules.bulk_enrichment as bulk_enrichment
    import modules.sc_cell_go as sc_cell_go
    from modules.sc_batch import _available_path
    from modules.sc_batch_export import batch_csv_package_dirs
    from modules.sc_cell_go import SCCellGOEnrichment

    input_path = _tiny_batch_adata(tmp_path)
    _, package_dirs = batch_csv_package_dirs(str(tmp_path), "all_cluster_go_package")
    deg_rows = []
    for cluster in range(9):
        for gene in ("GENE1", "GENE2", "GENE3"):
            deg_rows.append({
                "comparison_id": "Treatment_vs_Control",
                "comparison": "Treatment vs Control",
                "experimental_group": "Treatment",
                "control_group": "Control",
                "deg_scope": "per_cluster",
                "cluster": str(cluster),
                "gene": gene,
                "log2FC": 2.0,
                "padj": 0.01,
                "inference_unit": "biological_sample",
            })
    deg_path = _available_path(
        package_dirs["deg"], "sc_pseudobulk_deg_Treatment_vs_Control", ".csv",
    )
    pd.DataFrame(deg_rows).to_csv(deg_path, index=False)
    local_dir = tmp_path / "all_cluster_gene_sets"
    local_dir.mkdir()
    (local_dir / "GO_Biological_Process_2023.gmt").write_text(
        "Shared term\tGO:1\tGENE1\tGENE2\tGENE3\n", encoding="utf-8",
    )

    def fake_run_ora_full(*_args, **_kwargs):
        result = pd.DataFrame({
            "Term": ["Shared term"],
            "P-value": [0.001],
            "Adjusted P-value": [0.001],
            "Genes": ["GENE1;GENE2;GENE3"],
            "Significant": [True],
        })
        result.attrs["testing_family"] = {"n_terms": 1}
        return result

    rendered_clusters = []

    def fake_render(_analysis, frame, **_kwargs):
        rendered_clusters.append(str(frame["cluster"].iloc[0]))
        return [], []

    monkeypatch.setattr(sc_cell_go, "run_ora_full", fake_run_ora_full)
    monkeypatch.setattr(bulk_enrichment, "_render_enrichment_variants", fake_render)

    result = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={
            "source_level": "pseudobulk", "execution_mode": "local",
            "deg_source_files": [deg_path], "method": "ORA", "direction_mode": "up",
            "min_genes": 3, "ora_min_size": 1,
            "local_gene_set_dir": str(local_dir), "export_full_tables": True,
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))

    output = next(item["file_path"] for item in result["result_files"]
                  if item["label"].startswith("Pseudobulk ORA enrichment ("))
    table = pd.read_csv(output)
    assert set(table["cluster"].astype(str)) == {str(cluster) for cluster in range(9)}
    assert set(rendered_clusters) == {str(cluster) for cluster in range(9)}
    assert not any("GO 图仅展示最小 FDR 最优的前" in warning
                   for warning in result["summary"]["plot_warnings"])


def test_pseudobulk_enrichment_can_run_a_single_selected_cluster(tmp_path, monkeypatch):
    """An explicit cluster selection limits both the statistical table and plots."""
    import modules.bulk_enrichment as bulk_enrichment
    import modules.sc_cell_go as sc_cell_go
    from modules.sc_batch import _available_path
    from modules.sc_batch_export import batch_csv_package_dirs
    from modules.sc_cell_go import SCCellGOEnrichment

    input_path = _tiny_batch_adata(tmp_path)
    _, package_dirs = batch_csv_package_dirs(str(tmp_path), "selected_cluster_go_package")
    rows = []
    for cluster in ("0", "1"):
        for gene in ("GENE1", "GENE2", "GENE3"):
            rows.append({
                "comparison_id": "Treatment_vs_Control", "comparison": "Treatment vs Control",
                "experimental_group": "Treatment", "control_group": "Control",
                "deg_scope": "per_cluster", "cluster": cluster, "gene": gene,
                "log2FC": 2.0, "padj": 0.01, "inference_unit": "biological_sample",
            })
    deg_path = _available_path(
        package_dirs["deg"], "sc_pseudobulk_deg_Treatment_vs_Control", ".csv",
    )
    pd.DataFrame(rows).to_csv(deg_path, index=False)
    local_dir = tmp_path / "selected_cluster_gene_sets"
    local_dir.mkdir()
    (local_dir / "GO_Biological_Process_2023.gmt").write_text(
        "Shared term\tGO:1\tGENE1\tGENE2\tGENE3\n", encoding="utf-8",
    )

    def fake_run_ora_full(*_args, **_kwargs):
        result = pd.DataFrame({
            "Term": ["Shared term"], "P-value": [0.001],
            "Adjusted P-value": [0.001], "Genes": ["GENE1;GENE2;GENE3"],
            "Significant": [True],
        })
        result.attrs["testing_family"] = {"n_terms": 1}
        return result

    rendered_clusters = []

    def fake_render(_analysis, frame, **_kwargs):
        rendered_clusters.append(str(frame["cluster"].iloc[0]))
        return [], []

    monkeypatch.setattr(sc_cell_go, "run_ora_full", fake_run_ora_full)
    monkeypatch.setattr(bulk_enrichment, "_render_enrichment_variants", fake_render)
    result = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={
            "source_level": "pseudobulk", "execution_mode": "local",
            "deg_source_files": [deg_path], "method": "ORA", "direction_mode": "up",
            "target_clusters": "1", "min_genes": 3, "ora_min_size": 1,
            "local_gene_set_dir": str(local_dir), "export_full_tables": True,
        }, progress_callback=lambda *_: None,
    ).run(str(input_path))

    output = next(item["file_path"] for item in result["result_files"]
                  if item["label"].startswith("Pseudobulk ORA enrichment ("))
    table = pd.read_csv(output)
    assert set(table["cluster"].astype(str)) == {"1"}
    assert rendered_clusters == ["1"]
    assert result["summary"]["selected_clusters"] == ["1"]


def test_pseudobulk_ora_maps_selected_ids_and_excludes_untested_background(tmp_path):
    from modules.sc_batch import _available_path
    from modules.sc_batch_export import batch_csv_package_dirs
    from modules.sc_cell_go import SCCellGOEnrichment

    input_path = _tiny_batch_adata(tmp_path)
    _, package_dirs = batch_csv_package_dirs(str(tmp_path), 'mapped_go_package')
    deg = pd.DataFrame({
        'comparison_id': ['Treatment_vs_Control'] * 5,
        'comparison': ['Treatment vs Control'] * 5,
        # Older pseudobulk exports used group_1/group_2 only.
        'group_1': ['Treatment'] * 5,
        'group_2': ['Control'] * 5,
        'deg_scope': ['all_cells'] * 5,
        'cluster': ['All'] * 5,
        'gene': [f'ENSG{i}' for i in range(1, 6)],
        'gene_symbol': ['A', 'B', 'C', 'D', 'E'],
        'log2FC': [2.0, 2.0, 2.0, 0.1, 3.0],
        'statistic': [4.0, 3.0, 2.0, 0.1, 5.0],
        # ENSG5 did not form a usable adjusted test and is not in ORA N.
        'padj': [0.01, 0.01, 0.01, 0.8, np.nan],
        'inference_unit': ['biological_sample'] * 5,
    })
    deg.to_csv(_available_path(
        package_dirs['deg'], 'sc_pseudobulk_deg_Treatment_vs_Control', '.csv',
    ), index=False)
    local_dir = tmp_path / 'mapped_gene_sets'
    local_dir.mkdir()
    (local_dir / 'GO_Biological_Process_2023.gmt').write_text(
        'Mapped term\tGO:1\tA\tB\tC\nOther term\tGO:2\tD\n', encoding='utf-8',
    )

    result = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={
            'source_level': 'pseudobulk', 'execution_mode': 'local',
            'export_folder': 'mapped_go_package', 'direction_mode': 'up',
            'min_genes': 3, 'ora_min_size': 1,
            'local_gene_set_dir': str(local_dir), 'show_enrichment_plots': False,
            'export_full_tables': True,
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))

    output = next(item['file_path'] for item in result['result_files']
                  if item['label'].startswith('Pseudobulk ORA enrichment ('))
    table = pd.read_csv(output)
    assert table['status'].eq('completed').all()
    assert table['experimental_group'].eq('Treatment').all()
    assert table['control_group'].eq('Control').all()
    assert table['n_background_genes'].eq(4).all()
    assert table.loc[table['Term'].eq('Mapped term'), 'Genes'].iloc[0] == 'A;B;C'


def test_pseudobulk_deg_records_requested_batch_design(tmp_path):
    import anndata as ad
    from modules.sc_batch_export import SCPseudobulkDEG

    input_path = _tiny_batch_adata(tmp_path)
    adata = ad.read_h5ad(input_path)
    batch_by_sample = {"Ctrl_1": "B1", "Ctrl_2": "B2", "Trt_1": "B1", "Trt_2": "B2"}
    adata.obs["batch"] = adata.obs["sample_id"].map(batch_by_sample)
    adata.write_h5ad(input_path)
    result = SCPseudobulkDEG(
        project_dir=str(tmp_path),
        params={
            "sample_key": "sample_id", "condition_key": "condition",
            "batch_key": "batch", "analysis_scope": "all_cells",
            "method": "welch_logcpm", "min_samples_per_group": 2,
            "min_cells_per_sample_celltype": 1, "show_deg_figures": False,
            "export_full_tables": True,
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))
    output = next(item["file_path"] for item in result["result_files"]
                   if item["label"].startswith("Pseudobulk DEG:"))
    table = pd.read_csv(output)
    assert table["design"].eq("~ condition (batch not modeled)").all()
    assert any("未建模 batch_key=batch" in warning for warning in result["summary"]["warnings"])


def test_pseudobulk_enrichment_runs_local_preranked_gsea(tmp_path):
    from modules.sc_batch import _available_path
    from modules.sc_batch_export import batch_csv_package_dirs
    from modules.sc_cell_go import SCCellGOEnrichment

    input_path = _tiny_batch_adata(tmp_path)
    _, package_dirs = batch_csv_package_dirs(str(tmp_path), "pseudobulk_gsea_package")
    deg = pd.DataFrame({
        "comparison_id": ["Treatment_vs_Control"] * 6,
        "comparison": ["Treatment vs Control"] * 6,
        "experimental_group": ["Treatment"] * 6,
        "control_group": ["Control"] * 6,
        "deg_scope": ["all_cells"] * 6,
        "cluster": ["All"] * 6,
        "gene": [f"GENE{index}" for index in range(6)],
        "log2FC": [3.0, 2.0, 1.0, -1.0, -2.0, -3.0],
        "statistic": [6.0, 4.0, 2.0, -2.0, -4.0, -6.0],
        # GSEA uses the complete finite ranking even when independent
        # filtering leaves an adjusted P value missing for one gene.
        "padj": [0.01] * 5 + [np.nan],
        "method": ["pydeseq2"] * 6,
    })
    deg.to_csv(_available_path(
        package_dirs["deg"], "sc_pseudobulk_deg_Treatment_vs_Control", ".csv"
    ), index=False)
    local_dir = tmp_path / "go_gene_sets"
    local_dir.mkdir()
    (local_dir / "GO_Biological_Process_2023.gmt").write_text(
        "Positive term\tGO:1\tGENE0\tGENE1\tGENE2\n"
        "Negative term\tGO:2\tGENE3\tGENE4\tGENE5\n", encoding="utf-8"
    )

    result = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={
            "source_level": "pseudobulk", "deg_prefix": "sc_pseudobulk",
            "export_folder": "pseudobulk_gsea_package", "execution_mode": "local",
            "method": "GSEA", "min_genes": 3,
            "gsea_min_size": 2, "gsea_max_size": 6, "permutation_num": 100,
            "local_gene_set_dir": str(local_dir), "show_enrichment_plots": False,
            "export_full_tables": True,
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))
    output = next(item["file_path"] for item in result["result_files"]
                   if item["label"].startswith("Pseudobulk GSEA enrichment ("))
    table = pd.read_csv(output)
    assert table["method"].eq("GSEA").all()
    assert table["n_input_genes"].eq(6).all()
    assert table["gene_universe"].eq("complete_valid_ranking").all()
    assert {"Adjusted P-value", "nes"}.issubset(table.columns)
    assert set(table["Term"]) == {"Positive term", "Negative term"}


def test_cell_level_deg_exports_log2fc_and_adjusted_pvalues(tmp_path):
    from modules.sc_cell_deg import SCCellLevelDEG

    input_path = _tiny_batch_adata(tmp_path)
    module = SCCellLevelDEG(
        project_dir=str(tmp_path),
        params={
            "condition_key": "condition", "cluster_key": "leiden",
            "analysis_scope": "both", "comparisons": "Treatment-vs-Control",
            "min_cells_per_group": 2, "export_folder": "cell_deg_package",
            "export_full_tables": True,
        },
        progress_callback=lambda *_: None,
    )
    result = module.run(str(input_path))
    result_files = [item for item in result["result_files"]
                    if item["label"].startswith("Cell-level exploratory DEG (")]
    assert len(result_files) == 2
    assert {"all_cells", "per_cluster"} == {
        item["label"].split("(", 1)[1].split(")", 1)[0] for item in result_files
    }
    result_file = next(item for item in result_files if "(all_cells)" in item["label"])
    table = pd.read_csv(result_file["file_path"])
    assert "03_differential_expression" in result_file["file_path"]
    assert {"comparison", "experimental_group", "control_group", "cluster", "log2FC", "p.adjust", "avg_log2FC", "p_val_adj"}.issubset(table.columns)
    assert table["deg_scope"].eq("all_cells").all()
    assert table["cluster"].eq("All").all()
    assert table["statistical_status"].eq("exploratory_no_biological_replicates").all()


def test_cell_level_go_writes_one_table_per_deg_comparison(tmp_path, monkeypatch):
    from modules.sc_batch import _available_path
    from modules.sc_batch_export import batch_csv_package_dirs
    from modules.sc_cell_go import SCCellGOEnrichment

    input_path = _tiny_batch_adata(tmp_path)
    _, package_dirs = batch_csv_package_dirs(str(tmp_path), "go_package")
    deg = pd.DataFrame({
        "comparison_id": ["G1_Treatment_vs_G1_Control"] * 6,
        "comparison": ["G1_Treatment vs G1_Control"] * 6,
        "experimental_group": ["G1_Treatment"] * 6,
        "control_group": ["G1_Control"] * 6,
        "gene": [f"GENE{index}" for index in range(6)],
        "log2FC": [1.0] * 6,
        "p.adjust": [0.01] * 6,
    })
    deg["deg_scope"] = "all_cells"
    deg["cluster"] = "All"
    deg.to_csv(_available_path(package_dirs["deg"], "sc_cell_level_deg_all_cells_G1_Treatment_vs_G1_Control", ".csv"), index=False)
    per_cluster = deg.copy()
    per_cluster["deg_scope"] = "per_cluster"
    per_cluster["cluster"] = "0"
    per_cluster.to_csv(_available_path(package_dirs["deg"], "sc_cell_level_deg_per_cluster_G1_Treatment_vs_G1_Control", ".csv"), index=False)

    fake_results = pd.DataFrame({
        "Gene_set": ["GO_Biological_Process_2023"], "Term": ["Example pathway"],
        "Overlap": ["3/100"], "P-value": [0.001], "Adjusted P-value": [0.01],
        "Odds Ratio": [2.0], "Combined Score": [6.0], "Genes": ["GENE1;GENE2;GENE3"],
    })
    fake_gp = types.SimpleNamespace(
        enrichr=lambda **_kwargs: types.SimpleNamespace(results=fake_results)
    )
    monkeypatch.setitem(sys.modules, "gseapy", fake_gp)

    result = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={"export_folder": "go_package", "min_genes": 5, "execution_mode": "enrichr",
                "export_full_tables": True},
        progress_callback=lambda *_: None,
    ).run(str(input_path))
    go_files = [item for item in result["result_files"]
                if item["label"].startswith("Cell-level GO enrichment (")]
    assert len(go_files) == 2
    go_file = next(item for item in go_files if ", per_cluster)" in item["label"])
    table = pd.read_csv(go_file["file_path"])
    assert "04_go_enrichment" in go_file["file_path"]
    assert table.loc[0, "comparison_id"] == "G1_Treatment_vs_G1_Control"
    assert table.loc[0, "deg_scope"] == "per_cluster"
    assert table.loc[0, "cluster"] == 0
    assert table.loc[0, "status"] == "completed"


def test_cell_level_go_exports_pathway_csv_by_default(tmp_path, monkeypatch):
    """Enrichment tables are task-visible without a separate export flag."""
    from modules.sc_batch import _available_path
    from modules.sc_batch_export import batch_csv_package_dirs
    from modules.sc_cell_go import SCCellGOEnrichment

    input_path = _tiny_batch_adata(tmp_path)
    _, package_dirs = batch_csv_package_dirs(str(tmp_path), "default_go_export")
    pd.DataFrame({
        "comparison_id": ["Treatment_vs_Control"] * 5,
        "comparison": ["Treatment vs Control"] * 5,
        "experimental_group": ["Treatment"] * 5,
        "control_group": ["Control"] * 5,
        "deg_scope": ["all_cells"] * 5,
        "cluster": ["All"] * 5,
        "gene": [f"GENE{index}" for index in range(5)],
        "log2FC": [1.0] * 5,
        "p.adjust": [0.01] * 5,
    }).to_csv(_available_path(
        package_dirs["deg"], "sc_cell_level_deg_all_cells_Treatment_vs_Control", ".csv",
    ), index=False)
    fake_results = pd.DataFrame({
        "Term": ["Example pathway"], "Overlap": ["3/100"],
        "P-value": [0.001], "Adjusted P-value": [0.01],
        "Odds Ratio": [2.0], "Combined Score": [6.0],
        "Genes": ["GENE0;GENE1;GENE2"],
    })
    monkeypatch.setitem(sys.modules, "gseapy", types.SimpleNamespace(
        enrichr=lambda **_kwargs: types.SimpleNamespace(results=fake_results),
    ))

    result = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={"export_folder": "default_go_export", "min_genes": 5,
                "execution_mode": "enrichr"},
        progress_callback=lambda *_: None,
    ).run(str(input_path))

    output = next(item for item in result["result_files"]
                  if item["label"].startswith("Cell-level GO enrichment ("))
    assert "04_go_enrichment" in output["file_path"]
    assert {"Term", "Adjusted P-value", "Genes", "status"}.issubset(
        pd.read_csv(output["file_path"]).columns
    )
    assert result["summary"]["full_table_exported"] is True


def test_cell_level_go_runs_bp_cc_mf_from_local_gene_sets(tmp_path, monkeypatch):
    from modules.sc_batch import _available_path
    from modules.sc_batch_export import batch_csv_package_dirs
    from modules.sc_cell_go import LOCAL_GO_ASPECTS, SCCellGOEnrichment

    input_path = _tiny_batch_adata(tmp_path)
    _, package_dirs = batch_csv_package_dirs(str(tmp_path), "local_go_package")
    deg = pd.DataFrame({
        "comparison_id": ["Treatment_vs_Control"] * 6,
        "comparison": ["Treatment vs Control"] * 6,
        "experimental_group": ["Treatment"] * 6,
        "control_group": ["Control"] * 6,
        "deg_scope": ["all_cells"] * 6,
        "cluster": ["All"] * 6,
        "gene": [f"GENE{index}" for index in range(6)],
        "log2FC": [1.0] * 6,
        "p.adjust": [0.01] * 6,
    })
    deg.to_csv(_available_path(
        package_dirs["deg"], "sc_cell_level_deg_all_cells_Treatment_vs_Control", ".csv"
    ), index=False)
    local_dir = tmp_path / "local_go_gene_sets"
    local_dir.mkdir()
    for library in LOCAL_GO_ASPECTS.values():
        (local_dir / f"{library}.gmt").write_text(
            "Example local GO term\tGO:0000000\tGENE0\tGENE1\tGENE2\n", encoding="utf-8"
        )

    fake_results = pd.DataFrame({
        "Term": ["Example local GO term"], "Overlap": ["3/100"],
        "P-value": [0.001], "Adjusted P-value": [0.01],
        "Odds Ratio": [2.0], "Combined Score": [6.0], "Genes": ["GENE0;GENE1;GENE2"],
    })
    calls = []
    fake_gp = types.SimpleNamespace(
        enrich=lambda **kwargs: calls.append(kwargs) or types.SimpleNamespace(results=fake_results)
    )
    monkeypatch.setitem(sys.modules, "gseapy", fake_gp)

    result = SCCellGOEnrichment(
        project_dir=str(tmp_path),
        params={
            "export_folder": "local_go_package", "execution_mode": "local",
            "go_aspects": "BP_CC_MF", "local_gene_set_dir": str(local_dir),
            "min_genes": 5, "export_full_tables": True,
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))
    go_files = [item for item in result["result_files"]
                if item["label"].startswith("Cell-level GO enrichment (")]
    assert len(go_files) == 3
    assert len(calls) == 3
    assert all(isinstance(call["gene_sets"], dict) for call in calls)
    exported_sets = {pd.read_csv(item["file_path"]).loc[0, "gene_set"] for item in go_files}
    assert exported_sets == set(LOCAL_GO_ASPECTS.values())


def test_server_manifest_endpoints_scan_then_submit(test_project, tmp_path, monkeypatch):
    import io
    import worker
    from app import create_app

    root = tmp_path / "server_data"
    _make_10x_dir(root / "run" / "ctrl_1")
    _make_10x_dir(root / "run" / "trt_1")
    monkeypatch.setenv("SC_BATCH_SOURCE_ROOTS", str(root))
    submitted = []
    monkeypatch.setattr(worker, "submit_task", lambda *args: submitted.append(args) or True)
    app = create_app()
    app.config["TESTING"] = True

    with app.test_client() as client:
        scan = client.post(
            f"/projects/{test_project}/upload/discover-10x-directory",
            data={"source_root": str(root)},
        )
        assert scan.status_code == 200
        payload = scan.get_json()
        assert payload["n_samples"] == 2
        assert client.get(payload["download_url"]).status_code == 200

        manifest = pd.DataFrame({
            "sample_id": ["ctrl_1", "trt_1"],
            "matrix_dir": ["run/ctrl_1", "run/trt_1"],
            "condition": ["Control", "Treatment"],
            "replicate": ["1", "1"],
        }).to_csv(index=False).encode()
        submit = client.post(
            f"/projects/{test_project}/upload/import-10x-manifest",
            data={
                "source_root": str(root), "dataset_name": "demo",
                "manifest_file": (io.BytesIO(manifest), "design.csv"),
            },
            content_type="multipart/form-data",
        )
    assert submit.status_code == 200
    assert submit.get_json()["n_samples"] == 2
    assert submitted and submitted[0][2] == "sc_batch_import"
