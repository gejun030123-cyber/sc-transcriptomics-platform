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
