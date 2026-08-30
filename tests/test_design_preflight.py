"""Focused checks for the non-mutating experiment-design preflight."""

import numpy as np
import pandas as pd
import anndata as ad


def _bulk_adata(groups=("Ctrl", "Ctrl", "Treat", "Treat"), values=None):
    if values is None:
        values = [[index + 1, index + 2] for index in range(len(groups))]
    values = np.asarray(values, dtype=float)
    return ad.AnnData(
        X=values,
        obs=pd.DataFrame({"condition": list(groups)}, index=[f"s{i}" for i in range(len(groups))]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )


def test_bulk_preflight_builds_group_sizes_and_default_contrast():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _bulk_adata(), "bulk_deg", {"groupby": "condition", "method": "t-test"}
    )

    assert result["status"] == "ready"
    assert result["contrast"]["selected"] == {"group1": "Ctrl", "group2": "Treat"}
    assert result["contrast"]["groups"] == [{"name": "Ctrl", "n": 2}, {"name": "Treat", "n": 2}]
    assert result["recommended_params"]["groupby"] == "condition"


def test_bulk_preflight_blocks_count_method_for_continuous_expression():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _bulk_adata(values=[[0.2, 1.1], [0.3, 1.2], [0.4, 1.3], [0.5, 1.4]]),
        "bulk_deg", {"groupby": "condition", "method": "deseq2"},
    )

    assert result["status"] == "blocked"
    assert any("原始计数" in message for message in result["blockers"])


def test_bulk_preflight_blocks_singleton_group():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _bulk_adata(groups=("Ctrl", "Ctrl", "Treat")),
        "bulk_deg", {"groupby": "condition", "method": "t-test"},
    )

    assert result["status"] == "blocked"
    assert any(check["name"] == "生物学重复" for check in result["checks"])


def test_proportion_preflight_marks_missing_sample_replicates_exploratory():
    from modules.design_preflight import build_design_preflight

    adata = ad.AnnData(
        X=np.ones((8, 2)),
        obs=pd.DataFrame({
            "celltype": ["T", "B"] * 4,
            "condition": ["Ctrl"] * 4 + ["Treat"] * 4,
            "sample_id": [f"cell_{i}" for i in range(8)],
        }, index=[f"cell_{i}" for i in range(8)]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )

    result = build_design_preflight(
        adata,
        "proportion",
        {"groupby": "celltype", "condition_key": "condition", "sample_key": "sample_id"},
    )

    assert result["status"] == "exploratory"
    sample_check = next(check for check in result["checks"] if check["name"] == "样本级组成设计")
    assert sample_check["status"] == "warning"
def _sc_deg_adata():
    """Single-cell AnnData without biological sample/condition metadata."""
    return ad.AnnData(
        X=np.ones((12, 3)),
        obs=pd.DataFrame({
            "batch": ["b1"] * 6 + ["b2"] * 6,
            "leiden": ["0", "1"] * 6,
            "celltype": ["T", "B"] * 6,
        }, index=[f"cell_{i}" for i in range(12)]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )


def test_sc_cell_deg_preflight_blocks_missing_condition_column():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _sc_deg_adata(), "sc_cell_deg",
        {"comparison_type": "condition", "condition_key": "condition"},
    )

    assert result["status"] == "blocked"
    assert any(check["name"] == "条件列" and check["status"] == "blocked"
               for check in result["checks"])
    assert any("缺少条件列" in message for message in result["blockers"])
    assert "当前可用分组列" in result["blockers"][0] or "可用分组列" in result["blockers"][0]


def test_sc_cell_deg_preflight_blocks_missing_sample_column_for_sample_modes():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _sc_deg_adata(), "sc_cell_deg",
        {"comparison_type": "between_samples_all_cells", "sample_key": "sample_id"},
    )

    assert result["status"] == "blocked"
    assert any(check["name"] == "样本列" and check["status"] == "blocked"
               for check in result["checks"])


def test_sc_pseudobulk_preflight_blocks_missing_sample_and_condition():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _sc_deg_adata(), "sc_pseudobulk_deg",
        {"sample_key": "sample_id", "condition_key": "condition",
         "cluster_key": "leiden"},
    )

    assert result["status"] == "blocked"
    names = {check["name"] for check in result["checks"] if check["status"] == "blocked"}
    assert {"样本列", "条件列"} <= names
    assert any("缺少样本列" in message for message in result["blockers"])
    assert any("缺少条件列" in message for message in result["blockers"])


def test_sc_pseudobulk_preflight_warns_when_sample_key_is_technical_batch():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _sc_deg_adata(), "sc_pseudobulk_deg",
        {"sample_key": "batch", "condition_key": "condition",
         "cluster_key": "leiden"},
    )

    sample_checks = [check for check in result["checks"] if check["name"] == "样本列"]
    assert any(check["status"] == "warning" and "技术 batch" in check["message"]
               for check in sample_checks)
    assert result["status"] == "blocked"  # condition still missing


def test_sc_pseudobulk_preflight_passes_with_real_sample_design():
    from modules.design_preflight import build_design_preflight

    adata = ad.AnnData(
        X=np.ones((20, 3)),
        obs=pd.DataFrame({
            "sample_id": ["S1"] * 5 + ["S2"] * 5 + ["S3"] * 5 + ["S4"] * 5,
            "condition": ["Ctrl"] * 10 + ["Treat"] * 10,
            "leiden": ["0", "1"] * 10,
        }, index=[f"cell_{i}" for i in range(20)]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )

    result = build_design_preflight(
        adata, "sc_pseudobulk_deg",
        {"sample_key": "sample_id", "condition_key": "condition",
         "cluster_key": "leiden", "analysis_scope": "per_cluster"},
    )

    assert result["status"] == "ready"
    assert result["recommended_params"]["sample_key"] == "sample_id"

