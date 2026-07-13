import json
import os


def test_pipeline_report_manifest_and_resume_plan(tmp_path):
    from models import AnalysisTask, PipelineRun, ResultFile
    from modules.reporting.pipeline_report import write_pipeline_report

    project_dir = tmp_path / "project"
    project_dir.mkdir()
    input_path = project_dir / "input.h5ad"
    input_path.write_text("input", encoding="utf-8")
    qc_output = project_dir / "intermediate" / "qc_output.h5ad"
    qc_output.parent.mkdir()
    qc_output.write_text("qc", encoding="utf-8")

    run = PipelineRun(
        id="run-1",
        project_id="proj-1",
        name="SC pipeline",
        analysis_type="sc",
        status="failed",
        current_module="normalize",
        progress=50,
        input_path=str(input_path),
        modules_json=json.dumps(["qc", "normalize", "hvg"]),
        params_json=json.dumps({"normalize": {"method": "log1p"}}),
        task_ids_json=json.dumps(["task-qc", "task-norm"]),
        error_traceback="normalize failed",
    )
    qc_task = AnalysisTask(
        id="task-qc",
        project_id="proj-1",
        module_name="qc",
        status="completed",
        output_adata_path=str(qc_output),
        result_json=json.dumps({"cells_after": 10}),
    )
    norm_task = AnalysisTask(
        id="task-norm",
        project_id="proj-1",
        module_name="normalize",
        status="failed",
        error_traceback="bad input",
    )
    result_file = ResultFile(
        task_id="task-qc",
        project_id="proj-1",
        file_type="csv",
        category="table",
        label="QC table",
        file_path=str(project_dir / "results" / "qc.csv"),
    )

    artifacts = write_pipeline_report(
        str(project_dir),
        run,
        [qc_task, norm_task],
        {"task-qc": [result_file], "task-norm": []},
        {"qc": "质控", "normalize": "标准化"},
    )

    assert os.path.isfile(artifacts["manifest_path"])
    assert os.path.isfile(artifacts["report_path"])
    manifest = artifacts["manifest"]
    assert manifest["resume_plan"]["can_resume"] is True
    assert manifest["resume_plan"]["failed_module"] == "normalize"
    assert manifest["resume_plan"]["remaining_modules"] == ["normalize", "hvg"]
    assert manifest["resume_plan"]["resume_input_path"] == str(qc_output)


def test_resume_pipeline_run_api_creates_new_run(tmp_path, monkeypatch):
    from app import create_app
    from config import Config
    from database import get_conn
    from models import AnalysisTask, PipelineRun

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    db_path = tmp_path / "test.db"
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "DB_PATH", str(db_path))

    submitted = {}

    def fake_submit_pipeline_run(**kwargs):
        submitted.update(kwargs)
        return True

    monkeypatch.setattr("worker.submit_pipeline_run", fake_submit_pipeline_run)

    app = create_app()
    app.config["TESTING"] = True
    pid = "resume_project"
    project_dir = Config.project_dir(pid)
    os.makedirs(project_dir, exist_ok=True)
    input_path = os.path.join(project_dir, "input.h5ad")
    qc_output = os.path.join(project_dir, "qc_output.h5ad")
    with open(input_path, "w", encoding="utf-8") as fh:
        fh.write("input")
    with open(qc_output, "w", encoding="utf-8") as fh:
        fh.write("qc")

    conn = get_conn()
    try:
        conn.execute("INSERT INTO projects (id, name, status) VALUES (?, ?, ?)", (pid, "Resume Project", "active"))
        conn.commit()
    finally:
        conn.close()

    qc_task = AnalysisTask(
        id="task-qc",
        project_id=pid,
        module_name="qc",
        status="completed",
        output_adata_path=qc_output,
    )
    qc_task.save()
    failed_task = AnalysisTask(
        id="task-norm",
        project_id=pid,
        module_name="normalize",
        status="failed",
    )
    failed_task.save()
    run = PipelineRun(
        id="run-failed",
        project_id=pid,
        name="Failed run",
        analysis_type="sc",
        status="failed",
        current_module="normalize",
        input_path=input_path,
        modules_json=json.dumps(["qc", "normalize"]),
        params_json=json.dumps({"normalize": {"method": "log1p"}}),
        task_ids_json=json.dumps(["task-qc", "task-norm"]),
    )
    run.save()

    with app.test_client() as client:
        resp = client.post(f"/api/projects/{pid}/pipeline-runs/{run.id}/resume")

    assert resp.status_code == 201
    data = resp.get_json()
    assert data["resume_plan"]["remaining_modules"] == ["normalize"]
    assert submitted["modules"] == ["normalize"]
    assert submitted["input_path"] == qc_output
