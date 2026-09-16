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


def test_bulk_deg_preflight_warns_when_low_expression_filter_is_absent():
    """Count-based DEG must not silently accept an unfiltered gene set."""
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _bulk_adata(), "bulk_deg", {"groupby": "condition", "method": "deseq2"}
    )

    check = next(item for item in result["checks"] if item["name"] == "低表达基因过滤")
    assert check["status"] == "warning"
    assert any("低表达基因过滤" in message for message in result["warnings"])


def test_bulk_deg_preflight_passes_after_normalization_gene_filter():
    from modules.design_preflight import build_design_preflight

    adata = _bulk_adata()
    adata.uns["normalization"] = {"gene_expression_filter": {
        "applied": True, "n_genes_before": 100, "n_genes_after": 40,
    }}

    result = build_design_preflight(
        adata, "bulk_deg", {"groupby": "condition", "method": "deseq2"}
    )

    check = next(item for item in result["checks"] if item["name"] == "低表达基因过滤")
    assert check["status"] == "pass"
    assert "100 → 40" in check["message"]


def test_bulk_deg_preflight_accepts_normalized_x_with_verified_raw_layer():
    """A normalizer's log-scale X must not hide its DESeq2 raw-count layer."""
    from modules.design_preflight import build_design_preflight

    raw = np.asarray([[10, 2], [12, 3], [20, 5], [24, 6]], dtype=int)
    adata = _bulk_adata(values=np.log2(raw + 1))
    adata.layers['raw'] = raw
    adata.uns['normalization'] = {'is_log_transformed': True, 'gene_expression_filter': {
        'applied': True, 'n_genes_before': 10, 'n_genes_after': 2,
    }}

    result = build_design_preflight(
        adata, "bulk_deg", {"groupby": "condition", "method": "deseq2"}
    )

    assert result["status"] == "ready"
    scale_check = next(item for item in result["checks"] if item["name"] == "表达量尺度")
    assert scale_check["status"] == "pass"
    assert scale_check["value"]["count_layer"] == "raw"


def test_bulk_preflight_blocks_singleton_group():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _bulk_adata(groups=("Ctrl", "Ctrl", "Treat")),
        "bulk_deg", {"groupby": "condition", "method": "t-test"},
    )

    assert result["status"] == "blocked"
    assert any(check["name"] == "生物学重复" for check in result["checks"])


def test_bulk_preflight_blocks_comparisons_from_another_groupby_column():
    from modules.design_preflight import build_design_preflight

    adata = _bulk_adata(groups=("Ctrl_B", "Ctrl_B", "Drug_B", "Drug_B"))
    adata.obs["treatment"] = ["Ctrl", "Ctrl", "Drug", "Drug"]
    result = build_design_preflight(
        adata,
        "bulk_deg",
        {
            "groupby": "treatment", "method": "t-test",
            "comparisons": "Drug_B-vs-Ctrl_B",
        },
    )

    assert result["status"] == "blocked"
    assert any("不属于当前分组列" in message for message in result["blockers"])


def test_bulk_preflight_allows_comparison_using_a_custom_group():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _bulk_adata(groups=("Ctrl", "Ctrl", "Low", "Low", "High", "High")),
        "bulk_deg",
        {
            "groupby": "condition", "method": "t-test",
            "custom_groups": "Combined=Low+High",
            "comparisons": "Combined-vs-Ctrl",
        },
    )

    assert result["status"] == "ready"


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
    adata.layers["counts"] = adata.X.copy()

    result = build_design_preflight(
        adata, "sc_pseudobulk_deg",
        {"sample_key": "sample_id", "condition_key": "condition",
         "cluster_key": "leiden", "analysis_scope": "per_cluster"},
    )

    assert result["status"] == "ready"
    assert result["recommended_params"]["sample_key"] == "sample_id"


def test_sc_pseudobulk_preflight_requires_annotation_for_celltype_mode():
    from modules.design_preflight import build_design_preflight

    result = build_design_preflight(
        _sc_deg_adata(), "sc_pseudobulk_deg",
        {
            "sample_key": "batch", "condition_key": "batch",
            "analysis_scope": "per_cluster",
            "grouping_mode": "annotated_celltype", "celltype_key": "missing_celltype",
        },
    )

    annotation_check = next(
        check for check in result["checks"] if check["name"] == "细胞类型注释列"
    )
    assert annotation_check["status"] == "blocked"
    assert "先完成细胞注释" in annotation_check["message"]


def test_functional_state_preflight_previews_scoped_sample_celltype_units():
    from modules.design_preflight import build_design_preflight

    adata = ad.AnnData(
        X=np.ones((24, 2)),
        obs=pd.DataFrame({
            "sample_id": ["H1"] * 6 + ["H2"] * 6 + ["D1"] * 6 + ["D2"] * 6,
            "condition": ["Healthy"] * 12 + ["IBD"] * 12,
            "celltype": ["Epithelial"] * 4 + ["Immune"] * 2 + ["Epithelial"] * 4 + ["Immune"] * 2
                        + ["Epithelial"] * 4 + ["Immune"] * 2 + ["Epithelial"] * 4 + ["Immune"] * 2,
        }, index=[f"cell_{i}" for i in range(24)]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )

    result = build_design_preflight(
        adata, "functional_state",
        {
            "sample_key": "sample_id", "condition_key": "condition", "celltype_key": "celltype",
            "scope_key": "celltype", "scope_values": "Epithelial",
            "min_cells_per_sample_celltype": 3,
        },
    )

    check = next(item for item in result["checks"] if item["name"] == "功能状态统计预览")
    assert check["status"] == "pass"
    assert check["value"] == {
        "selected_cells": 16, "total_cells": 24,
        "samples_per_condition": {"Healthy": 2, "IBD": 2},
        "celltype_cells": {"Epithelial": 16},
        "valid_sample_x_celltype_units": 4, "total_sample_x_celltype_units": 4,
        "units_below_min_cells": 0, "min_cells_per_sample_celltype": 3,
    }
