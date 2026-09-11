"""Regression coverage for local precomputed-SCENIC result panels."""

import json
import os

import numpy as np
import pandas as pd
import pytest


def _scenic_adata(tmp_path, include_targets=True):
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(31)
    groups = np.repeat(["Epithelial", "Myeloid", "T cell"], 18)
    adata = anndata.AnnData(
        X=rng.poisson(2, size=(len(groups), 8)).astype(float),
        obs=pd.DataFrame({"celltype": pd.Categorical(groups)}),
        var=pd.DataFrame(index=[f"GENE{i}" for i in range(8)]),
    )
    regulons = ["HNF4A (+)", "SPI1 (+)", "TCF7 (+)", "RELA (+)", "JUN_extended"]
    auc = rng.uniform(0.01, 0.08, size=(len(groups), len(regulons)))
    auc[groups == "Epithelial", 0] += 0.55
    auc[groups == "Myeloid", 1] += 0.50
    auc[groups == "T cell", 2] += 0.45
    adata.obsm["X_aucell"] = pd.DataFrame(auc, index=adata.obs_names, columns=regulons)
    if include_targets:
        adata.uns["scenic_regulon_targets"] = {
            "HNF4A (+)": ["KRT8", "KRT18", "EPCAM"],
            "SPI1 (+)": ["LYZ", "TYROBP"],
            "TCF7 (+)": ["CCR7", "LTB"],
        }
    path = tmp_path / "scenic_input.h5ad"
    adata.write_h5ad(path)
    return path


def test_scenic_module_is_registered_and_ordered():
    from modules import MODULE_REGISTRY, PIPELINE_DEPS, PIPELINE_ORDER
    from modules.schemas import PARAM_SCHEMAS, SC_MODULE_LIST

    assert MODULE_REGISTRY["scenic"].MODULE_NAME == "scenic"
    assert any(item["name"] == "scenic" for item in SC_MODULE_LIST)
    assert "scenic" in PARAM_SCHEMAS
    assert PIPELINE_ORDER.index("annotation") < PIPELINE_ORDER.index("scenic")
    assert PIPELINE_DEPS["scenic"] == ["annotation"]


def test_scenic_allows_an_explicit_cluster_grouping_without_celltype(tmp_path):
    """RSS accepts a reviewed cluster column when no celltype column exists."""
    import anndata
    from modules.scenic import ScenicAnalysis

    path = _scenic_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.obs["leiden"] = pd.Categorical(
        ["0"] * 18 + ["1"] * 18 + ["2"] * 18,
    )
    del adata.obs["celltype"]
    adata.write_h5ad(path)

    result = ScenicAnalysis(str(tmp_path), {
        "groupby": "leiden", "max_heatmap_cells": 40,
        "show_network": False, "show_correlation": False,
    }, None).run(str(path))

    assert result["summary"]["groupby"] == "leiden"
    assert result["summary"]["n_groups"] == 3


def test_scenic_safe_name_keeps_unicode_analysis_ids_distinct():
    from modules.scenic import _safe_name

    assert _safe_name("分析 结果") == "分析_结果"
    assert _safe_name("分析/结果") == "分析_结果"


def test_scenic_rss_is_chunk_invariant_and_matches_scipy_definition():
    from scipy.spatial.distance import jensenshannon
    from modules.scenic import ScenicAnalysis

    rng = np.random.default_rng(4)
    matrix = rng.random((17, 9))
    matrix[:, -1] = 0.0
    labels = np.array(["A"] * 7 + ["B"] * 6 + ["C"] * 4)
    groups = ["A", "B", "C"]
    chunked = ScenicAnalysis._rss(matrix, labels, groups, chunk_size=2)
    larger_chunk = ScenicAnalysis._rss(matrix, labels, groups, chunk_size=16)
    assert np.allclose(chunked, larger_chunk, atol=1e-14)
    for group_index, group in enumerate(groups):
        q = (labels == group).astype(float)
        q /= q.sum()
        for regulon_index in range(matrix.shape[1] - 1):
            p = matrix[:, regulon_index] / matrix[:, regulon_index].sum()
            expected = 1.0 - jensenshannon(p, q, base=2.0)
            assert chunked[group_index, regulon_index] == pytest.approx(expected)
    assert np.all(chunked[:, -1] == 0.0)


def test_scenic_rejects_explicit_zero_display_limits(tmp_path):
    from modules.scenic import ScenicAnalysis

    with pytest.raises(ValueError, match="AUCell 热图展示细胞数"):
        ScenicAnalysis(str(tmp_path), {
            "groupby": "celltype", "max_heatmap_cells": 0,
        }, None).run(str(_scenic_adata(tmp_path)))


def test_scenic_generates_core_aucell_rss_network_and_correlation_outputs(tmp_path):
    from modules.scenic import ScenicAnalysis

    result = ScenicAnalysis(str(tmp_path), {
        "groupby": "celltype", "auc_obsm_key": "X_aucell", "top_n_regulons": 4,
        "max_heatmap_cells": 40, "top_n_rss": 3, "show_network": True,
        "show_correlation": True, "top_n_correlation_regulons": 4,
    }, lambda *_: None).run(str(_scenic_adata(tmp_path)))

    assert os.path.exists(result["output_adata"])
    assert result["summary"]["n_regulons"] == 5
    assert result["summary"]["n_groups"] == 3
    assert result["summary"]["n_network_edges"] > 0
    output_names = {os.path.basename(item["file_path"]) for item in result["result_files"]}
    assert {"scenic_aucell_heatmap.png", "scenic_rss_bubble.png", "scenic_rss_bar.png",
            "scenic_regulon_network.png", "scenic_regulon_correlation.png",
            "regulon_specificity_score.csv", "group_regulon_aucell_summary.csv",
            "regulon_network_edges.csv", "regulon_coactivity_correlation.csv",
            "analysis_manifest.json"}.issubset(output_names)

    root = tmp_path / "results" / "scenic" / "manual_run"
    rss = pd.read_csv(root / "regulon_specificity_score.csv")
    assert rss.shape[0] == 15
    assert rss["rss"].between(0, 1).all()
    epithelial_top = rss[rss["cell_group"].eq("Epithelial")].sort_values("rss", ascending=False).iloc[0]
    assert epithelial_top["regulon"] == "HNF4A (+)"
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["activity_method"] == "precomputed_AUCell"
    assert manifest["network"]["status"] == "available"


def test_scenic_requires_named_precomputed_aucell_matrix(tmp_path):
    from modules.scenic import ScenicAnalysis

    path = _scenic_adata(tmp_path)
    anndata = pytest.importorskip("anndata")
    adata = anndata.read_h5ad(path)
    del adata.obsm["X_aucell"]
    adata.write_h5ad(path)
    with pytest.raises(ValueError, match="AUCell"):
        ScenicAnalysis(str(tmp_path), {"groupby": "celltype"}, None).run(str(path))


def test_scenic_keeps_core_outputs_when_targets_are_not_available(tmp_path):
    from modules.scenic import ScenicAnalysis

    result = ScenicAnalysis(str(tmp_path), {
        "groupby": "celltype", "max_heatmap_cells": 40, "show_network": True,
    }, None).run(str(_scenic_adata(tmp_path, include_targets=False)))
    assert result["summary"]["n_network_edges"] == 0
    assert any("跳过网络图" in warning for warning in result["summary"]["warnings"])
