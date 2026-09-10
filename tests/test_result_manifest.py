import json


def test_build_task_manifest_records_files(tmp_path):
    from models import AnalysisTask
    from modules.reporting.result_manifest import build_task_manifest

    output_adata = tmp_path / "output.h5ad"
    output_adata.write_text("h5ad placeholder", encoding="utf-8")
    csv_path = tmp_path / "result.csv"
    csv_path.write_text("gene,score\nA,1\n", encoding="utf-8")

    task = AnalysisTask(
        id="task-123",
        project_id="proj-1",
        module_name="qc",
        params_json=json.dumps({"mito_perc": 0.2}),
    )
    result = {
        "output_adata": str(output_adata),
        "summary": {"cells_after": 42},
        "result_files": [
            {
                "file_type": "csv",
                "category": "qc",
                "label": "QC table",
                "file_path": str(csv_path),
            }
        ],
    }

    manifest = build_task_manifest(task, result)

    assert manifest["manifest_version"] == 1
    assert manifest["project_id"] == "proj-1"
    assert manifest["task_id"] == "task-123"
    assert manifest["module_name"] == "qc"
    assert manifest["params"] == {"mito_perc": 0.2}
    assert manifest["summary"] == {"cells_after": 42}
    assert manifest["output_adata"]["exists"] is True
    assert manifest["result_files"][0]["file"]["exists"] is True
    assert manifest["result_files"][0]["file"]["size_bytes"] > 0
    assert manifest["review_evidence"]["module_name"] == "qc"


def test_write_task_manifest_creates_manifest_file(tmp_path):
    from models import AnalysisTask
    from modules.reporting.result_manifest import write_task_manifest

    project_dir = tmp_path / "project"
    project_dir.mkdir()
    output_adata = project_dir / "intermediate" / "qc_output.h5ad"
    output_adata.parent.mkdir()
    output_adata.write_text("h5ad placeholder", encoding="utf-8")

    task = AnalysisTask(
        id="task-456",
        project_id="proj-1",
        module_name="qc",
        params_json="{}",
    )
    result = {
        "output_adata": str(output_adata),
        "summary": {"cells_before": 50, "cells_after": 45},
        "result_files": [],
    }

    rf = write_task_manifest(str(project_dir), task, result)

    assert rf["file_type"] == "json"
    assert rf["category"] == "manifest"
    manifest_path = project_dir / "results" / "manifests" / "qc_task-456_analysis_manifest.json"
    assert rf["file_path"] == str(manifest_path)
    assert manifest_path.exists()

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["summary"]["cells_after"] == 45
    assert manifest["output_adata"]["exists"] is True
    assert manifest["review_evidence"]["status"] in {"pass", "review", "warning"}


def test_manifest_records_module_declared_failure(tmp_path):
    from models import AnalysisTask
    from modules.reporting.result_manifest import build_task_manifest

    task = AnalysisTask(
        id="task-failed", project_id="proj-1", module_name="qc", params_json="{}",
    )
    manifest = build_task_manifest(task, {
        "output_adata": str(tmp_path / "input.h5ad"),
        "result_files": [],
        "summary": {"error": "required dependency is unavailable"},
    })

    assert manifest["status"] == "failed"
    assert manifest["error"] == "required dependency is unavailable"


def test_failed_result_does_not_snapshot_its_input_h5ad(tmp_path):
    import worker

    input_h5ad = tmp_path / "input.h5ad"
    input_h5ad.write_bytes(b"input")
    task = type("Task", (), {"id": "task-failed", "module_name": "qc"})()
    result = worker._snapshot_task_artifacts(task, str(tmp_path), {
        "output_adata": str(input_h5ad), "result_files": [], "error": "failed",
    })

    assert result["output_adata"] is None
    assert not list((tmp_path / "results" / "task_artifacts").rglob("*.h5ad"))
