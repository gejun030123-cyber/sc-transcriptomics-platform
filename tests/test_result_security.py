"""Result route and analysis input path security tests."""
import json
import io
import os
import zipfile

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


def _make_csv_result(pid, filename="table.csv"):
    from config import Config
    from models import AnalysisTask, ResultFile

    results_dir = Config.results_dir(pid)
    os.makedirs(results_dir, exist_ok=True)
    path = os.path.join(results_dir, filename)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("gene,log2FC,padj\nTP53,2.1,0.01\nGAPDH,0.2,0.9\nEGFR,-1.4,0.04\n")
    task = AnalysisTask(project_id=pid, module_name="bulk_deg", status="completed")
    task.save()
    return ResultFile.create(
        task_id=task.id,
        project_id=pid,
        file_type="csv",
        category="table",
        label="DEG table",
        file_path=path,
    )


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
    assert ok.status_code == 410
    assert "retired" in ok.get_json()["error"]


def test_csv_result_preview_is_project_scoped_and_filterable(client):
    _make_project("project_a")
    _make_project("project_b")
    result_file = _make_csv_result("project_b")

    wrong = client.get(f"/api/projects/project_a/result-file/{result_file.id}/table-preview")
    assert wrong.status_code == 403

    response = client.get(
        f"/api/projects/project_b/result-file/{result_file.id}/table-preview",
        query_string={"search": "tp", "column": "padj", "max": "0.05"},
    )
    assert response.status_code == 200
    data = response.get_json()
    assert data["matched_rows"] == 1
    assert data["rows"][0]["gene"] == "TP53"


def test_design_preflight_api_returns_safe_contrast_preview(client):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from config import Config

    _make_project("project_a")
    uploads = Config.uploads_dir("project_a")
    os.makedirs(uploads, exist_ok=True)
    input_path = os.path.join(uploads, "bulk_counts.h5ad")
    ad.AnnData(
        X=np.asarray([[1, 3], [2, 4], [7, 9], [8, 10]], dtype=float),
        obs=pd.DataFrame({"condition": ["Ctrl", "Ctrl", "Treat", "Treat"]}, index=["s1", "s2", "s3", "s4"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    ).write_h5ad(input_path)

    response = client.post(
        "/api/projects/project_a/design-preflight",
        json={
            "file_path": input_path,
            "module_name": "bulk_deg",
            "params": {"groupby": "condition", "method": "t-test"},
        },
    )

    assert response.status_code == 200
    data = response.get_json()
    assert data["status"] == "ready"
    assert data["contrast"]["selected"] == {"group1": "Ctrl", "group2": "Treat"}


def test_bulk_deg_form_blocks_insufficient_replicates_before_task_submission(client):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from config import Config
    from models import AnalysisTask

    _make_project("project_a")
    uploads = Config.uploads_dir("project_a")
    os.makedirs(uploads, exist_ok=True)
    input_path = os.path.join(uploads, "singleton_group.h5ad")
    ad.AnnData(
        X=np.asarray([[1, 3], [2, 4], [7, 9]], dtype=float),
        obs=pd.DataFrame({"condition": ["Ctrl", "Ctrl", "Treat"]}, index=["s1", "s2", "s3"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    ).write_h5ad(input_path)

    response = client.post(
        "/projects/project_a/analyze/bulk_deg",
        data={"input_path": input_path, "groupby": "condition", "method": "t-test"},
        follow_redirects=True,
    )

    assert response.status_code == 200
    assert '分析前检查未通过'.encode() in response.data
    assert AnalysisTask.get_by_project("project_a") == []


def test_project_plot_archive_contains_gallery_and_sources(client):
    _make_project("project_a")
    _, result_file = _make_task_with_result("project_a")

    resp = client.get('/projects/project_a/results/plots-archive')

    assert resp.status_code == 200
    assert resp.mimetype == 'application/zip'
    with zipfile.ZipFile(io.BytesIO(resp.data)) as archive:
        names = archive.namelist()
        assert 'plot_gallery.html' in names
        assert 'manifest.json' in names
        plot_sources = [name for name in names if name.startswith('static_images/')]
        assert len(plot_sources) == 2
        manifest = json.loads(archive.read('manifest.json'))
        assert manifest[0]['label'] == 'QC plot'
        assert manifest[0]['archive_path'] == plot_sources[0]


def test_task_result_page_has_individual_and_batch_plot_exports(client):
    _make_project("project_a")
    task, _ = _make_task_with_result("project_a")

    resp = client.get(f'/projects/project_a/task/{task.id}')

    assert resp.status_code == 200
    assert '科研图（默认展示）'.encode() in resp.data
    assert '下载 PNG'.encode() in resp.data
    assert b'result-figure-header' in resp.data
    assert b'result-figure-preview' in resp.data
    assert b'Plotly.newPlot' not in resp.data


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
