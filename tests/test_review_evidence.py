from types import SimpleNamespace


def test_qc_review_evidence_flags_high_filtering():
    from modules.reporting.review_evidence import build_review_evidence

    evidence = build_review_evidence(
        "qc",
        {"cells_before": 100, "cells_after": 20, "pct_removed": 80},
        [
            SimpleNamespace(file_path="/tmp/qc_filter_summary.json", label="QC filter summary", category="qc"),
            SimpleNamespace(file_path="/tmp/qc_doublet_score_histogram.json", label="Doublet score", category="histogram"),
        ],
    )

    assert evidence["module_name"] == "qc"
    assert evidence["status"] == "warning"
    assert any(c["name"] == "细胞保留率" and c["status"] == "warning" for c in evidence["checks"])
    assert any(c["name"] == "双细胞证据" and c["status"] == "pass" for c in evidence["checks"])


def test_annotation_review_evidence_uses_unknown_ratio_and_confidence():
    from modules.reporting.review_evidence import build_review_evidence

    evidence = build_review_evidence(
        "annotation",
        {
            "n_celltypes": 4,
            "celltype_counts": {"T cell": 40, "B cell": 30, "Unknown": 30},
            "mean_confidence": 0.62,
        },
        [
            SimpleNamespace(file_path="/tmp/annotation_marker_expression_box.json", label="Marker expression", category="boxplot"),
            SimpleNamespace(file_path="/tmp/annotation_score_umap.json", label="Annotation score UMAP", category="umap"),
        ],
    )

    assert evidence["status"] == "review"
    assert any(c["name"] == "Unknown 比例" and c["status"] == "review" for c in evidence["checks"])
    assert any(c["name"] == "平均置信度" and c["status"] == "pass" for c in evidence["checks"])


def test_non_single_cell_review_module_returns_none():
    from modules.reporting.review_evidence import build_review_evidence

    assert build_review_evidence("bulk_qc", {"n_samples": 6}, []) is None
