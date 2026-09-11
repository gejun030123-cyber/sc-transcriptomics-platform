"""Focused contracts for sample-level KNN neighbourhood abundance analysis."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _ibd_neighbourhood_adata(path):
    ad = pytest.importorskip("anndata")
    rng = np.random.default_rng(12)
    rows, coordinates = [], []
    # The two conditions occupy different regions.  Each sample contributes
    # multiple cells, but inference must still count only eight samples.
    for condition, center in (("Control", -3.0), ("IBD", 3.0)):
        for sample_index in range(4):
            sample = f"{condition}_{sample_index + 1}"
            for cell_index in range(24):
                coordinates.append([center + rng.normal(0, 0.35), rng.normal(0, 0.35), rng.normal(0, 0.2)])
                rows.append({
                    "sample_id": sample,
                    "condition": condition,
                    "celltype": "Stress epithelial" if condition == "IBD" else "Absorptive enterocyte",
                    "cell": f"{sample}_cell_{cell_index}",
                })
    obs = pd.DataFrame(rows).set_index("cell")
    values = rng.poisson(2, size=(len(obs), 5)).astype(float)
    adata = ad.AnnData(values, obs=obs, var=pd.DataFrame(index=[f"G{index}" for index in range(5)]))
    representation = np.asarray(coordinates)
    adata.obsm["X_pca"] = representation
    adata.obsm["X_umap"] = representation[:, :2]
    adata.write_h5ad(path)
    return path


def test_neighborhood_da_uses_sample_level_proportions_and_exports_audit(tmp_path):
    from modules.neighborhood_da import NeighborhoodDAAnalysis

    input_path = _ibd_neighbourhood_adata(tmp_path / "ibd.h5ad")
    result = NeighborhoodDAAnalysis(
        project_dir=str(tmp_path),
        params={
            "sample_key": "sample_id", "condition_key": "condition",
            "representation": "X_pca", "comparisons": "IBD-vs-Control",
            "n_neighbourhoods": 20, "neighbourhood_size": 16,
            "min_cells_per_sample": 10, "min_samples_per_condition": 3,
            "random_state": 1,
        },
        progress_callback=lambda *_: None,
    ).run(str(input_path))

    assert result["summary"]["inference_unit"] == "biological_sample"
    assert result["summary"]["n_neighbourhoods"] <= 20
    root = tmp_path / "results" / "neighborhood_da" / "manual_run"
    sample_proportions = pd.read_csv(root / "neighbourhood_sample_proportions.csv")
    statistics = pd.read_csv(root / "neighbourhood_differential_abundance.csv")
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))

    # Four donors in each group, not 96 cells in each group, form the tests.
    first = statistics.iloc[0]
    assert first["n_samples_experimental"] == 4
    assert first["n_samples_reference"] == 4
    assert first["inference_unit"] == "biological_sample"
    assert set(sample_proportions["sample_id"]) == {
        "Control_1", "Control_2", "Control_3", "Control_4",
        "IBD_1", "IBD_2", "IBD_3", "IBD_4",
    }
    assert statistics["fdr_bh"].notna().all()
    assert manifest["representation"] == "X_pca"
    assert manifest["statistical_contract"].startswith("Mann-Whitney")
    assert any(item["label"].startswith("Neighborhood DA:") for item in result["result_files"])


def test_neighborhood_da_refuses_technical_batch_without_explicit_confirmation(tmp_path):
    from modules.neighborhood_da import NeighborhoodDAAnalysis
    anndata = pytest.importorskip("anndata")

    input_path = _ibd_neighbourhood_adata(tmp_path / "ibd.h5ad")
    adata = anndata.read_h5ad(input_path)
    adata.obs["batch"] = adata.obs["sample_id"].astype(str)
    adata.write_h5ad(input_path)
    with pytest.raises(ValueError, match="技术 batch"):
        NeighborhoodDAAnalysis(
            project_dir=str(tmp_path),
            params={"sample_key": "batch", "condition_key": "condition"},
            progress_callback=lambda *_: None,
        ).run(str(input_path))


def test_neighborhood_da_preflight_requires_sample_design_and_high_dimensional_representation(tmp_path):
    from anndata import read_h5ad
    from modules.design_preflight import build_design_preflight

    input_path = _ibd_neighbourhood_adata(tmp_path / "ibd.h5ad")
    adata = read_h5ad(input_path)
    ready = build_design_preflight(
        adata, "neighborhood_da",
        {"sample_key": "sample_id", "condition_key": "condition", "representation": "X_pca"},
    )
    assert ready["status"] == "ready"
    assert {item["name"] for item in ready["checks"]} >= {
        "邻域丰度样本设计", "UMAP 嵌入", "邻域高维表示",
    }

    del adata.obsm["X_pca"]
    blocked = build_design_preflight(
        adata, "neighborhood_da",
        {"sample_key": "sample_id", "condition_key": "condition", "representation": "auto"},
    )
    assert blocked["status"] == "blocked"
