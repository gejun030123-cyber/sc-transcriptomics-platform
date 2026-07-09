"""Result route and analysis input path security tests."""
import json
import os

import pytest

from app import create_app


@pytest.fixture
def client(tmp_path, monkeypatch):
    from config import Config

    data_dir = tmp_path / "data"
    data_dir.mkdir(exist_ok=True)
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "DB_PATH", str(tmp_path / "test.db"))

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as c:
        yield c


def _make_project(pid):
    from database import get_conn
    from config import Config

    os.makedirs(Config.project_dir(pid), exist_ok=True)
    conn = get_conn()
    try:
        conn.execute(
            "INSERT OR REPLACE INTO projects (id, name, status) VALUES (?,?,?)",
            (pid, pid, "active"),
        )
        conn.commit()
    finally:
        conn.close()


def _make_task_with_result(pid, filename="plot.json"):
    from config import Config
    from models import AnalysisTask, ResultFile

    project_dir = Config.project_dir(pid)
    plots_dir = os.path.join(project_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    path = os.path.join(plots_dir, filename)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"data": [], "layout": {}}, fh)

    task = AnalysisTask(project_id=pid, module_name="qc", status="completed")
    task.save()
    result_file = ResultFile.create(
        task_id=task.id,
        project_id=pid,
        file_type="plotly_json",
        category="qc",
        label="QC plot",
        file_path=path,
    )
    return task, result_file


def test_result_page_rejects_cross_project_file(client):
    _make_project("project_a")
    _make_project("project_b")
    _, result_file = _make_task_with_result("project_b")

    resp = client.get(f"/projects/project_a/results/file/{result_file.id}")

    assert resp.status_code == 302
    assert b"data" not in resp.data


def test_project_result_file_api_requires_matching_project(client):
    _make_project("project_a")
    _make_project("project_b")
    _, result_file = _make_task_with_result("project_b")

    wrong = client.get(f"/api/projects/project_a/result-file/{result_file.id}")
    assert wrong.status_code == 403

    legacy_without_project = client.get(f"/api/result-file/{result_file.id}")
    assert legacy_without_project.status_code == 400

    ok = client.get(f"/api/projects/project_b/result-file/{result_file.id}")
    assert ok.status_code == 200
    assert ok.get_json() == {"data": [], "layout": {}}


def test_analysis_submit_rejects_input_outside_project(client, tmp_path):
    from models import AnalysisTask

    _make_project("project_a")
    external_input = tmp_path / "outside.h5ad"
    external_input.write_text("not real h5ad", encoding="utf-8")

    resp = client.post(
        "/projects/project_a/analyze/qc",
        data={"input_path": str(external_input)},
    )

    assert resp.status_code == 302
    assert AnalysisTask.get_by_project("project_a") == []
