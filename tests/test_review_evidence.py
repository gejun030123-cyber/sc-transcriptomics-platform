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


def test_annotation_optional_confidence_is_not_reported_as_warning():
    from modules.reporting.review_evidence import build_review_evidence

    evidence = build_review_evidence(
        "annotation",
        {
            "n_celltypes": 3,
            "celltype_counts": {"T cell": 40, "B cell": 30},
            "marker_selection": {"warnings": ["cluster 2 最终仅保留 1 个 marker"]},
        },
        [
            SimpleNamespace(file_path="/tmp/annotation_marker_expression_dotplot.png", label="Marker expression", category="dotplot"),
        ],
    )

    confidence = next(c for c in evidence["checks"] if c["name"] == "注释置信度")
    confidence_umap = next(c for c in evidence["checks"] if c["name"] == "置信度 UMAP")
    assert confidence["status"] == "pass"
    assert confidence_umap["status"] == "pass"
    assert any(c["name"] == "Marker 选择提示" for c in evidence["checks"])


def test_proportion_review_evidence_requires_sample_level_statistics():
    from modules.reporting.review_evidence import build_review_evidence

    evidence = build_review_evidence(
        "proportion",
        {
            "sample_level_inference_ready": True,
            "sample_key": "sample_id",
            "condition_key": "condition",
            "sample_level_n_tests": 4,
        },
        [
            SimpleNamespace(file_path="/tmp/sample_level_cell_proportions.csv", label="Sample-level Cell Proportions", category="table"),
            SimpleNamespace(file_path="/tmp/sample_level_proportion_tests.csv", label="Sample-level Proportion Tests", category="table"),
        ],
    )

    assert evidence["module_name"] == "proportion"
    assert any(c["name"] == "统计单位" and c["status"] == "pass" for c in evidence["checks"])


def test_batch_correct_review_evidence_uses_pre_post_deltas():
    from modules.reporting.review_evidence import build_review_evidence

    evidence = build_review_evidence(
        "batch_correct",
        {
            "evaluation_comparison": {
                "before": {
                    "abs_asw_batch": 0.62,
                    "asw_bio": 0.48,
                    "mean_neighbor_batch_entropy": 0.31,
                    "mean_neighbor_same_batch_fraction": 0.82,
                    "bio_label_key": "celltype",
                },
                "after": {
                    "abs_asw_batch": 0.18,
                    "asw_bio": 0.51,
                    "mean_neighbor_batch_entropy": 0.66,
                    "mean_neighbor_same_batch_fraction": 0.35,
                    "bio_label_key": "celltype",
                },
                "delta": {
                    "abs_asw_batch": -0.44,
                    "asw_bio": 0.03,
                    "mean_neighbor_batch_entropy": 0.35,
                    "mean_neighbor_same_batch_fraction": -0.47,
                },
                "sampling": {"used_size": 400, "strategy": "batch_stratified"},
                "comparison_scope": "embedding_and_graph",
                "warnings": [],
            }
        },
        [
            SimpleNamespace(
                file_path="/tmp/batch_evaluation_pre_post_harmony.csv",
                label="Batch integration pre/post metrics",
                category="table",
            )
        ],
    )

    assert evidence["module_name"] == "batch_correct"
    assert evidence["status"] == "pass"
    assert any(c["name"] == "Batch ASW" and c["status"] == "pass" for c in evidence["checks"])
    assert any(c["name"] == "生物结构保留" and c["status"] == "pass" for c in evidence["checks"])


def test_non_single_cell_review_module_returns_none():
    from modules.reporting.review_evidence import build_review_evidence

    assert build_review_evidence("bulk_qc", {"n_samples": 6}, []) is None


def test_result_interpretation_has_safe_generic_fallback():
    from modules.reporting.review_evidence import build_result_interpretation
    text = build_result_interpretation('hvg', {'n_cells': 100, 'n_hvgs': 2000}, [])
    assert text['title'] == '结果解读'
    assert '100' in text['conclusion']


def test_result_interpretation_uses_review_evidence_when_available():
    from modules.reporting.review_evidence import build_result_interpretation
    text = build_result_interpretation('annotation', {'n_celltypes': 1, 'celltype_counts': {'Unknown': 10}}, [])
    assert text['cautions']
