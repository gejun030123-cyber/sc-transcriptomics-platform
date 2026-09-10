"""Regression tests for single-cell DEG statistical contracts and compact output."""

import json
import os

import numpy as np
import pandas as pd


def _comparison_adata(tmp_path):
    import anndata as ad

    rows = np.asarray([
        [30, 2, 1], [25, 3, 1], [1, 28, 2], [2, 24, 2],
        [18, 3, 1], [16, 4, 1], [2, 17, 2], [1, 15, 2],
    ], dtype=np.int64)
    obs = pd.DataFrame({
        "sample_id": ["S1"] * 4 + ["S2"] * 4,
        "condition": ["Control"] * 4 + ["Treatment"] * 4,
        "leiden": ["A", "A", "B", "B"] * 2,
    }, index=[f"cell_{i}" for i in range(len(rows))])
    adata = ad.AnnData(
        X=rows.astype(np.float32), obs=obs,
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    adata.layers["counts"] = rows
    path = tmp_path / "comparison.h5ad"
    adata.write_h5ad(path)
    return path


def test_within_sample_cluster_deg_is_compact_and_uses_count_log1p(tmp_path):
    from modules.sc_cell_deg import SCCellLevelDEG

    input_path = _comparison_adata(tmp_path)
    result = SCCellLevelDEG(
        project_dir=str(tmp_path),
        params={
            "comparison_type": "within_sample_clusters", "sample_key": "sample_id",
            "selected_sample": "S1", "cluster_key": "leiden",
            "comparisons": "A-vs-B", "min_cells_per_group": 2,
            "export_full_tables": False,
        }, progress_callback=lambda *_: None,
    ).run(str(input_path))

    assert not [item for item in result["result_files"] if item["file_type"] == "csv"]
    assert result["summary"]["comparison_type"] == "within_sample_clusters"
    source = result["summary"]["deg_source_file"]
    table = pd.read_csv(source)
    assert table["comparison_type"].eq("within_sample_clusters").all()
    assert table["selected_sample"].eq("S1").all()
    assert table["expression_scale"].eq("log1p_from_layers[counts]").all()
    assert np.isfinite(pd.to_numeric(table["log2FC"])).all()


def test_pearson_residual_input_reconstructs_log1p_before_cell_deg(tmp_path):
    import anndata as ad
    from modules.sc_cell_deg import SCCellLevelDEG

    input_path = _comparison_adata(tmp_path)
    adata = ad.read_h5ad(input_path)
    # Simulate a completed Pearson-residual normalization: X has signed values
    # while protected counts are still available for a valid DE scale.
    adata.X = np.linspace(-3, 3, adata.n_obs * adata.n_vars).reshape(adata.shape)
    adata.uns["normalization"] = {"x_contains": "pearson_residuals"}
    adata.write_h5ad(input_path)

    result = SCCellLevelDEG(
        project_dir=str(tmp_path),
        params={
            "condition_key": "condition", "analysis_scope": "all_cells",
            "comparisons": "Treatment-vs-Control", "min_cells_per_group": 2,
            "export_full_tables": False,
        }, progress_callback=lambda *_: None,
    ).run(str(input_path))

    table = pd.read_csv(result["summary"]["deg_source_file"])
    assert np.isfinite(pd.to_numeric(table["log2FC"])).all()
    audit = next(item["file_path"] for item in result["result_files"] if item["file_type"] == "json")
    assert json.loads(open(audit, encoding="utf-8").read())["expression_scale"]["input_x_contains"] == "pearson_residuals"


def test_pseudobulk_deseq2_failure_never_silently_falls_back_to_welch(tmp_path, monkeypatch):
    import anndata as ad
    import modules.sc_batch_export as export

    input_path = _comparison_adata(tmp_path)
    original = ad.read_h5ad(input_path)
    replicate = original.copy()
    replicate.obs_names = [f"rep_{name}" for name in replicate.obs_names]
    replicate.obs["sample_id"] = replicate.obs["sample_id"].map({"S1": "S3", "S2": "S4"})
    merged = ad.concat([original, replicate], merge="same")
    merged.layers["counts"] = np.vstack([original.layers["counts"], replicate.layers["counts"]])
    merged.write_h5ad(input_path)

    def fail_deseq(*args, **kwargs):
        raise RuntimeError("intentional model failure")

    def unexpected_welch(*args, **kwargs):
        raise AssertionError("Welch fallback must not run")

    monkeypatch.setattr(export, "_deseq2_deg", fail_deseq)
    monkeypatch.setattr(export, "_welch_logcpm_deg", unexpected_welch)
    result = export.SCPseudobulkDEG(
        project_dir=str(tmp_path),
        params={
            "sample_key": "sample_id", "condition_key": "condition",
            "analysis_scope": "all_cells", "method": "deseq2",
            "min_samples_per_group": 2, "min_cells_per_sample_celltype": 1,
            "show_deg_figures": False, "export_full_tables": False,
        }, progress_callback=lambda *_: None,
    ).run(str(input_path))

    assert result["summary"]["n_model_failed_units"] == 1
    assert result["summary"]["deg_source_file"] == ""
    # Formal DEG tables remain absent after a model failure, while mandatory
    # donor-level QC remains available to diagnose the failed unit.
    assert not [item for item in result["result_files"]
                if item["label"].startswith("Pseudobulk DEG:")]
    assert any(item["label"].startswith("Pseudobulk sample QC:")
               for item in result["result_files"])


def test_pseudobulk_rejects_batch_condition_perfect_confounding():
    from modules.sc_batch_export import _pseudobulk_design_is_identifiable

    meta = pd.DataFrame({
        "condition": ["Control", "Control", "Treatment", "Treatment"],
        "batch": ["B1", "B1", "B2", "B2"],
    })
    identifiable, reason = _pseudobulk_design_is_identifiable(meta)

    assert not identifiable
    assert "完全混杂" in reason


def test_pseudobulk_uses_annotated_celltypes_as_distinct_deg_units(tmp_path):
    """Post-annotation DEG must aggregate sample × celltype, not Leiden IDs."""
    import anndata as ad
    from modules.sc_batch_export import SCPseudobulkDEG

    rows = []
    metadata = []
    for sample, condition in (
        ("C1", "Control"), ("C2", "Control"),
        ("T1", "Treatment"), ("T2", "Treatment"),
    ):
        for celltype, counts in (("T cells", [20, 2, 1]), ("B cells", [2, 20, 1])):
            for index in range(2):
                rows.append(counts)
                metadata.append({
                    "sample_id": sample, "condition": condition,
                    # Deliberately collapse Leiden labels: only the reviewed
                    # annotation can recover the two biological units here.
                    "leiden": "0", "celltype": celltype,
                })
    values = np.asarray(rows, dtype=np.int64)
    adata = ad.AnnData(
        X=values.astype(np.float32), obs=pd.DataFrame(metadata),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    adata.layers["counts"] = values
    input_path = tmp_path / "annotated.h5ad"
    adata.write_h5ad(input_path)

    result = SCPseudobulkDEG(
        project_dir=str(tmp_path),
        params={
            "sample_key": "sample_id", "condition_key": "condition",
            "analysis_scope": "per_cluster",
            "grouping_mode": "annotated_celltype", "celltype_key": "celltype",
            "method": "welch_logcpm", "min_samples_per_group": 2,
            "min_cells_per_sample_celltype": 2,
            "show_deg_figures": False, "export_full_tables": False,
        }, progress_callback=lambda *_: None,
    ).run(str(input_path))

    assert result["summary"]["analysis_grouping"] == "annotated_celltype"
    assert result["summary"]["analysis_group_key"] == "celltype"
    table = pd.read_csv(result["summary"]["deg_source_file"])
    assert set(table["cluster"]) == {"T cells", "B cells"}
    assert table["analysis_group_label"].eq("细胞类型").all()
    assert table["inference_unit"].eq("biological_sample").all()


def test_celltype_marker_deg_creates_task_bound_enrichment_sources(tmp_path):
    import anndata as ad
    from modules.deg import DEGAnalysis

    values = np.asarray(
        [[18, 2, 1]] * 6 + [[2, 18, 1]] * 6,
        dtype=np.int64,
    )
    adata = ad.AnnData(
        X=values.astype(np.float32),
        obs=pd.DataFrame({"celltype": ["T cells"] * 6 + ["B cells"] * 6}),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    adata.layers["counts"] = values
    input_path = tmp_path / "celltype_markers.h5ad"
    adata.write_h5ad(input_path)

    result = DEGAnalysis(
        project_dir=str(tmp_path),
        params={
            "groupby": "celltype", "method": "wilcoxon", "n_genes": 2,
            "show_dotplot": False, "show_marker_heatmap": False,
            "show_top_marker_umap_panel": False, "show_deg_counts_bar": False,
        }, progress_callback=lambda *_: None,
    ).run(str(input_path))

    summary = result["summary"]
    assert summary["analysis_grouping"] == "annotated_celltype"
    assert summary["deg_source_task_contract"] == "celltype_marker_vs_rest"
    assert len(summary["deg_source_files"]) == 2
    source = pd.read_csv(summary["deg_source_files"][0])
    assert source["analysis_group_label"].eq("细胞类型").all()
    assert source["inference_unit"].eq("cell").all()
    assert {"comparison_id", "log2FC", "padj", "cluster"}.issubset(source.columns)


def test_go_route_binds_the_selected_deg_task_instead_of_scanning_csv(tmp_path, monkeypatch):
    from app import create_app
    from config import Config
    from models import AnalysisTask, Project
    import worker

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "DB_PATH", str(tmp_path / "test.db"))
    monkeypatch.setattr(Config, "PLATFORM_ACCESS_PASSWORD", "")
    app = create_app()
    app.config["TESTING"] = True
    project_id = "sc_deg_source"
    Project(id=project_id, name="source").save()
    uploads = Config.uploads_dir(project_id)
    os.makedirs(uploads, exist_ok=True)
    input_path = os.path.join(uploads, "input.h5ad")
    _comparison_adata(tmp_path).replace(input_path)
    source_dir = os.path.join(Config.results_dir(project_id), "internal")
    os.makedirs(source_dir, exist_ok=True)
    source_path = os.path.join(source_dir, "selected_deg.csv")
    pd.DataFrame({
        "comparison_id": ["Treatment_vs_Control"], "deg_scope": ["all_cells"],
        "comparison": ["Treatment vs Control"], "gene": ["G1"],
        "log2FC": [1.0], "padj": [0.01], "inference_unit": ["cell"],
    }).to_csv(source_path, index=False)
    source_task = AnalysisTask(
        project_id=project_id, module_name="deg", status="completed",
        output_adata_path=input_path,
        result_json=json.dumps({
            "deg_source_files": [source_path], "deg_source_level": "cell_level",
            "deg_source_task_contract": "celltype_marker_vs_rest",
            "analysis_grouping": "annotated_celltype",
            "analysis_group_label": "细胞类型",
        }),
    )
    source_task.save()
    submitted = []
    monkeypatch.setattr(
        worker, "submit_task",
        lambda *args: submitted.append(args) or True,
    )
    with app.test_client() as client:
        page = client.get(f"/projects/{project_id}/analyze/sc_cell_go")
        response = client.post(
            f"/projects/{project_id}/analyze/sc_cell_go",
            data={"input_path": input_path, "deg_source_task_id": source_task.id},
        )

    assert page.status_code == 200
    assert 'id="sc-go-cluster-picker"' in page.get_data(as_text=True)
    assert 'data-clusters=' in page.get_data(as_text=True)
    assert 'data-group-label="细胞类型"' in page.get_data(as_text=True)
    assert response.status_code == 302
    assert submitted
    params = submitted[0][3]
    assert params["source_task_id"] == source_task.id
    assert params["deg_source_files"] == [source_path]
    assert params["source_level"] == "cell_level"


def test_go_route_rejects_a_cluster_not_in_the_selected_deg_task(tmp_path, monkeypatch):
    from app import create_app
    from config import Config
    from models import AnalysisTask, Project
    import worker

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "DB_PATH", str(tmp_path / "test.db"))
    monkeypatch.setattr(Config, "PLATFORM_ACCESS_PASSWORD", "")
    app = create_app()
    app.config["TESTING"] = True
    project_id = "sc_go_cluster_validation"
    Project(id=project_id, name="source").save()
    uploads = Config.uploads_dir(project_id)
    os.makedirs(uploads, exist_ok=True)
    input_path = os.path.join(uploads, "input.h5ad")
    _comparison_adata(tmp_path).replace(input_path)
    source_dir = os.path.join(Config.results_dir(project_id), "internal")
    os.makedirs(source_dir, exist_ok=True)
    source_path = os.path.join(source_dir, "selected_deg.csv")
    pd.DataFrame({
        "comparison_id": ["Treatment_vs_Control"], "deg_scope": ["per_cluster"],
        "comparison": ["Treatment vs Control"], "cluster": ["0"], "gene": ["G1"],
        "log2FC": [1.0], "padj": [0.01], "inference_unit": ["cell"],
    }).to_csv(source_path, index=False)
    source_task = AnalysisTask(
        project_id=project_id, module_name="sc_cell_deg", status="completed",
        output_adata_path=input_path,
        result_json=json.dumps({
            "deg_source_files": [source_path], "deg_source_level": "cell_level",
            "deg_source_task_contract": "condition",
        }),
    )
    source_task.save()
    submitted = []
    monkeypatch.setattr(worker, "submit_task", lambda *args: submitted.append(args) or True)

    with app.test_client() as client:
        response = client.post(
            f"/projects/{project_id}/analyze/sc_cell_go",
            data={
                "input_path": input_path, "deg_source_task_id": source_task.id,
                "target_clusters": "does-not-exist",
            }, follow_redirects=True,
        )

    assert response.status_code == 200
    assert "不存在以下 cluster" in response.get_data(as_text=True)
    assert not submitted
