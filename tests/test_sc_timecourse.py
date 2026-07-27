"""Tests for sample-aware measured-time single-cell analysis."""

import os

import numpy as np
import pandas as pd
import pytest


def test_order_timepoints_handles_common_day_labels():
    from modules.sc_timecourse import order_timepoints

    labels, values = order_timepoints(["D10", "D2", "D0", "D2"])

    assert labels == ["D0", "D2", "D10"]
    assert values["D0"] < values["D2"] < values["D10"]


def test_benjamini_hochberg_preserves_missing_values():
    from modules.sc_timecourse import benjamini_hochberg

    adjusted = benjamini_hochberg([0.01, 0.04, np.nan, 0.03])

    assert adjusted[0] == pytest.approx(0.03)
    assert np.isnan(adjusted[2])
    assert np.all(adjusted[[0, 1, 3]] <= 1)


def _make_temporal_adata(
    tmp_path,
    *,
    timepoints=("D0", "D3", "D7"),
    replicates=2,
    counts_are_raw=True,
    vary_composition=True,
    path=None,
):
    """Create a compact multi-timepoint dataset with real sample replicates."""
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(42)
    genes = ["EARLY", "LATE", "B_STATE", "G3", "G4", "G5"]
    counts, obs_rows = [], []
    for time_index, timepoint in enumerate(timepoints):
        for replicate in range(replicates):
            sample_id = f"{timepoint}_R{replicate + 1}"
            for celltype in ["A", "B"]:
                n_cells = (
                    12 + time_index * 2 if celltype == "A" and vary_composition
                    else 12 - time_index * 2 if celltype == "B" and vary_composition
                    else 12
                )
                for _ in range(n_cells):
                    row = rng.poisson(2, len(genes)).astype(float)
                    if celltype == "A":
                        row[0] += 10 - time_index * 3
                        row[1] += time_index * 4
                    else:
                        row[2] += time_index * 3
                    counts.append(row)
                    obs_rows.append({
                        "timepoint": timepoint,
                        "sample_id": sample_id,
                        "celltype": celltype,
                        "batch": "technical_1" if replicate == 0 else "technical_2",
                    })
    counts = np.asarray(counts)
    normalized = np.log1p(counts / np.maximum(counts.sum(axis=1, keepdims=True), 1) * 10_000)
    adata = anndata.AnnData(
        normalized,
        obs=pd.DataFrame(obs_rows),
        var=pd.DataFrame(index=genes),
    )
    adata.layers["counts"] = counts if counts_are_raw else normalized.copy()
    adata.obsm["X_umap"] = rng.normal(size=(adata.n_obs, 2))
    path = path or tmp_path / "temporal_input.h5ad"
    adata.write_h5ad(path)
    return path


def _run_timecourse(tmp_path, path, **params):
    from modules.sc_timecourse import SCTimecourseAnalysis

    module_params = {
        "timepoint_key": "timepoint",
        "sample_key": "sample_id",
        "celltype_key": "celltype",
        "min_replicates_per_timepoint": 2,
        "min_cells_per_celltype_sample": 5,
        "max_genes": 6,
        "top_dynamic_genes": 5,
    }
    module_params.update(params)
    module = SCTimecourseAnalysis(
        project_dir=str(tmp_path),
        params=module_params,
        progress_callback=lambda *_: None,
    )
    return module.run(str(path))


def test_sc_timecourse_uses_sample_level_units_and_exports_outputs(tmp_path):
    pytest.importorskip("scanpy")

    path = _make_temporal_adata(tmp_path)
    result = _run_timecourse(tmp_path, path)

    assert result["summary"]["inference_available"] is True
    assert result["summary"]["n_timepoints"] == 3
    assert os.path.exists(result["output_adata"])
    labels = {item["label"] for item in result["result_files"]}
    assert "Sample-level cell-type proportions" in labels
    assert "Cell-type temporal gene trends" in labels
    trends = pd.read_csv(tmp_path / "results" / "sc_timecourse_gene_trends.csv")
    assert set(trends["inference_status"]) == {"sample_level_kruskal"}
    assert {"A", "B"}.issubset(set(trends["celltype"]))


def test_sc_timecourse_requires_real_timepoint_column(tmp_path):
    pytest.importorskip("scanpy")
    from modules.sc_timecourse import SCTimecourseAnalysis

    path = _make_temporal_adata(tmp_path)
    module = SCTimecourseAnalysis(
        project_dir=str(tmp_path),
        params={"timepoint_key": "not_a_column"},
        progress_callback=lambda *_: None,
    )

    with pytest.raises(ValueError, match="时间列"):
        module.run(str(path))


def test_sc_timecourse_requires_nonmissing_timepoints(tmp_path):
    pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")
    from modules.sc_timecourse import SCTimecourseAnalysis

    path = _make_temporal_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.obs["timepoint"] = adata.obs["timepoint"].astype(object)
    adata.obs.loc[adata.obs.index[0], "timepoint"] = np.nan
    adata.write_h5ad(path)
    module = SCTimecourseAnalysis(
        project_dir=str(tmp_path),
        params={"timepoint_key": "timepoint"},
        progress_callback=lambda *_: None,
    )

    with pytest.raises(ValueError, match="缺失值"):
        module.run(str(path))


def test_sc_timecourse_rejects_missing_condition_values(tmp_path):
    pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")

    path = _make_temporal_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.obs["condition"] = np.where(
        adata.obs["sample_id"].astype(str).str.endswith("_R1"), "Control", "Treatment"
    )
    adata.obs["condition"] = adata.obs["condition"].astype(object)
    adata.obs.loc[adata.obs.index[0], "condition"] = ""
    adata.write_h5ad(path)

    with pytest.raises(ValueError, match="条件列.*缺失值"):
        _run_timecourse(tmp_path, path, condition_key="condition")


def test_sc_timecourse_clamps_single_replicate_to_descriptive(tmp_path):
    pytest.importorskip("scanpy")

    path = _make_temporal_adata(tmp_path, replicates=1)
    result = _run_timecourse(
        tmp_path, path,
        # Direct/API callers may send 1.  It must never activate p values.
        min_replicates_per_timepoint=1,
    )

    assert result["summary"]["inference_available"] is False
    assert result["summary"]["composition_inference_available"] is False
    assert result["summary"]["gene_inference_available"] is False
    assert any("少于 2 个独立样本" in warning for warning in result["summary"]["warnings"])
    composition = pd.read_csv(tmp_path / "results" / "sc_timecourse_composition_trends.csv")
    assert set(composition["inference_status"]) == {"descriptive_only"}


def test_sc_timecourse_rejects_technical_batch_as_sample_unit(tmp_path):
    pytest.importorskip("scanpy")

    result = _run_timecourse(tmp_path, _make_temporal_adata(tmp_path), sample_key="batch")

    assert result["summary"]["sample_key"] == ""
    assert result["summary"]["n_samples"] == 0
    assert result["summary"]["inference_available"] is False
    assert any("技术 batch" in warning for warning in result["summary"]["warnings"])


def test_sc_timecourse_rejects_cell_barcode_like_sample_column(tmp_path):
    pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")

    path = _make_temporal_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.obs["barcode_like"] = [f"cell_{index}" for index in range(adata.n_obs)]
    adata.write_h5ad(path)
    result = _run_timecourse(tmp_path, path, sample_key="barcode_like")

    assert result["summary"]["sample_key"] == ""
    assert result["summary"]["inference_available"] is False
    assert any("接近一细胞一个 ID" in warning for warning in result["summary"]["warnings"])


def test_sc_timecourse_rejects_repeated_sample_across_timepoints(tmp_path):
    pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")

    path = _make_temporal_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.obs["longitudinal_sample"] = [
        "donor_1" if sample_id.endswith("_R1") else f"unique_{timepoint}"
        for sample_id, timepoint in zip(adata.obs["sample_id"], adata.obs["timepoint"])
    ]
    adata.write_h5ad(path)
    result = _run_timecourse(tmp_path, path, sample_key="longitudinal_sample")

    assert result["summary"]["sample_key"] == ""
    assert result["summary"]["inference_available"] is False
    assert any("同一样本跨多个时间点" in warning for warning in result["summary"]["warnings"])


def test_sc_timecourse_disables_inference_for_repeated_donor_measurements(tmp_path):
    pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")

    path = _make_temporal_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.obs["donor_id"] = [
        "donor_1" if sample_id.endswith("_R1") else f"donor_{timepoint}_r2"
        for sample_id, timepoint in zip(adata.obs["sample_id"], adata.obs["timepoint"])
    ]
    adata.write_h5ad(path)
    result = _run_timecourse(tmp_path, path)

    assert result["summary"]["sample_key"] == "sample_id"
    assert result["summary"]["sample_inference_allowed"] is False
    assert result["summary"]["inference_available"] is False
    assert any("同一受试者跨多个时间点" in warning for warning in result["summary"]["warnings"])


def test_sc_timecourse_disables_gene_inference_for_non_count_layer(tmp_path):
    pytest.importorskip("scanpy")

    result = _run_timecourse(
        tmp_path,
        _make_temporal_adata(tmp_path, counts_are_raw=False),
    )

    assert result["summary"]["counts_layer_is_raw"] is False
    assert result["summary"]["composition_inference_available"] is True
    assert result["summary"]["gene_inference_available"] is False
    assert result["summary"]["inference_available"] is True
    assert any("已禁用基因推断性 p 值" in warning for warning in result["summary"]["warnings"])
    gene_trends = pd.read_csv(tmp_path / "results" / "sc_timecourse_gene_trends.csv")
    assert set(gene_trends["inference_status"]) == {"descriptive_only"}


def test_sc_timecourse_preserves_real_time_spacing_in_exports(tmp_path):
    pytest.importorskip("scanpy")

    result = _run_timecourse(
        tmp_path,
        _make_temporal_adata(tmp_path, timepoints=("D0", "D1", "D10")),
    )

    assert result["summary"]["n_timepoints"] == 3
    composition = pd.read_csv(tmp_path / "results" / "sc_timecourse_composition.csv")
    positions = composition.groupby("timepoint")["time_position"].first().to_dict()
    assert positions["D0"] == pytest.approx(0.0)
    assert positions["D1"] == pytest.approx(24.0)
    assert positions["D10"] == pytest.approx(240.0)


def test_sc_timecourse_does_not_mark_constant_composition_as_inferential(tmp_path):
    pytest.importorskip("scanpy")

    result = _run_timecourse(
        tmp_path,
        _make_temporal_adata(tmp_path, vary_composition=False),
        enable_gene_trends=False,
    )

    assert result["summary"]["composition_inference_available"] is False
    assert result["summary"]["inference_available"] is False
    assert any("无法计算 Kruskal" in warning for warning in result["summary"]["warnings"])
    trends = pd.read_csv(tmp_path / "results" / "sc_timecourse_composition_trends.csv")
    assert set(trends["inference_status"]) == {"test_unavailable"}


def test_sc_timecourse_rejects_batch_timepoint_without_explicit_confirmation(tmp_path):
    pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")

    path = _make_temporal_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.obs["batch"] = adata.obs["timepoint"].astype(str)
    adata.write_h5ad(path)
    result = _run_timecourse(tmp_path, path, timepoint_key="batch")

    assert result["summary"]["sample_inference_allowed"] is False
    assert result["summary"]["inference_available"] is False
    assert any("看起来是技术 batch" in warning for warning in result["summary"]["warnings"])


def test_sc_timecourse_result_interpretation_separates_descriptive_outputs():
    from modules.reporting.review_evidence import build_result_interpretation

    interpretation = build_result_interpretation("sc_timecourse", {
        "n_timepoints": 3,
        "sample_key": "",
        "inference_available": False,
        "composition_inference_available": False,
        "gene_inference_available": False,
    })

    assert any("描述性" in item for item in interpretation["evidence"])
    assert "没有满足设计与重复要求" in interpretation["cautions"][0]


def test_obs_columns_api_does_not_auto_select_technical_batch(test_project):
    pytest.importorskip("anndata")
    from app import create_app
    from config import Config

    path = os.path.join(Config.uploads_dir(test_project), "temporal_metadata.h5ad")
    _make_temporal_adata(tmp_path=None, path=path)
    app = create_app()
    app.config["TESTING"] = True
    response = app.test_client().get("/api/obs-columns", query_string={"file_path": path})
    payload = response.get_json()

    assert response.status_code == 200
    assert "timepoint" in payload["time_candidates"]
    assert "sample_id" in payload["sample_candidates"]
    assert "batch" not in payload["sample_candidates"]
    assert "batch" not in payload["time_candidates"]
