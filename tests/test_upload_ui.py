"""Regression tests for upload file-selection feedback and AJAX responses."""

import io
from pathlib import Path


def test_ajax_upload_returns_machine_readable_success(test_project):
    from app import create_app
    from config import Config

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        response = client.post(
            f"/projects/{test_project}/upload",
            data={"file": (io.BytesIO(b"gene,sample\nA,1\n"), "counts.csv")},
            content_type="multipart/form-data",
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

    assert response.status_code == 201
    assert response.get_json() == {
        "ok": True,
        "filename": "counts.csv",
        "size_bytes": 16,
    }
    assert Path(Config.uploads_dir(test_project), "counts.csv").read_bytes() == b"gene,sample\nA,1\n"


def test_ajax_upload_reports_invalid_extension_instead_of_redirecting(test_project):
    from app import create_app

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        response = client.post(
            f"/projects/{test_project}/upload",
            data={"file": (io.BytesIO(b"not supported"), "payload.exe")},
            content_type="multipart/form-data",
            headers={"X-Requested-With": "XMLHttpRequest"},
        )

    assert response.status_code == 400
    assert "不支持的文件格式" in response.get_json()["error"]


def test_upload_page_has_feedback_for_each_file_selector(test_project):
    from app import create_app

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        response = client.get(f"/projects/{test_project}/upload")

    html = response.get_data(as_text=True)
    assert "batch-zip-input').addEventListener('change'" in html
    assert "sc-batch-manifest').addEventListener('change'" in html
    assert "X-Requested-With" in html
    assert "已选择 manifest" in html
    assert "id=\"btn-choose-files\"" in html
    assert "bulk-count-import-card" in html
    assert "import-bulk-counts" in html
    assert "sample_metadata" in html
    assert "bulk-gene-annotation" in html
    assert "gene_annotation" in html
    assert "GTF/GFF" in html
    assert "scrollIntoView" in html
    assert html.index('id="upload-progress"') < html.index("多个 10x ZIP 合并导入")
    assert html.index('id="upload-result"') < html.index("多个 10x ZIP 合并导入")
