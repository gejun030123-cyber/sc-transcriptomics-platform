"""Presentation guards for compact task-result summaries."""

from routes.results import _result_summary_items, _review_evidence_for_display


def test_result_summary_omits_long_nested_and_technical_fields():
    items, omitted = _result_summary_items({
        "n_cells": 1250,
        "comparison": "Treatment vs Control",
        "input_path": "/protected/project/uploads/source.h5ad",
        "per_sample_qc": {"sample_a": {"cells": 600}},
        "warning": "x" * 97,
        "comparisons": ["Treatment vs Control", "Dose vs Control"],
    })

    displayed = {item["key"]: item["value"] for item in items}
    assert displayed == {
        "n_cells": "1250",
        "comparison": "Treatment vs Control",
        "comparisons": "Treatment vs Control；Dose vs Control",
    }
    assert omitted == 3


def test_result_summary_caps_the_number_of_cards():
    items, omitted = _result_summary_items({f"metric_{i}": i for i in range(17)})

    assert len(items) == 16
    assert items[-1] == {"key": "metric_15", "value": "15"}
    assert omitted == 1


def test_review_evidence_only_exposes_compact_values():
    evidence = {
        "title": "QC 审批证据",
        "status": "pass",
        "checks": [
            {
                "name": "删除原因拆分",
                "status": "pass",
                "message": "已记录。",
                "value": {"unique_removed_total": 12, "reason_overlap_cells": 2},
            },
            {
                "name": "逐样本统计",
                "status": "pass",
                "message": "已记录。",
                "value": {"sample_a": {"detected_rate": 0.02}},
            },
        ],
    }

    displayed = _review_evidence_for_display(evidence)
    compact, nested = displayed["checks"]

    assert compact["display_value"] == "unique removed total: 12；reason overlap cells: 2"
    assert "display_value" not in nested
    assert "value" not in nested
