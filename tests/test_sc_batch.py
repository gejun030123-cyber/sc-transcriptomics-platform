import os
import sys
import types

import numpy as np
import pandas as pd


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
    assert plan.loc[0, "comparison"] == "Control vs Treatment"
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
            "export_prefix": "demo",
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


def test_cell_level_deg_exports_log2fc_and_adjusted_pvalues(tmp_path):
    from modules.sc_cell_deg import SCCellLevelDEG

    input_path = _tiny_batch_adata(tmp_path)
    module = SCCellLevelDEG(
        project_dir=str(tmp_path),
        params={
            "condition_key": "condition", "cluster_key": "leiden",
            "analysis_scope": "both", "comparisons": "Treatment-vs-Control",
            "min_cells_per_group": 2, "export_folder": "cell_deg_package",
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
        params={"export_folder": "go_package", "min_genes": 5, "execution_mode": "enrichr"},
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
            "min_genes": 5,
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
