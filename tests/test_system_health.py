def test_dependency_status_has_groups_and_module_availability():
    from modules.platform.system_health import dependency_status
    from modules import MODULE_REGISTRY

    status = dependency_status()

    assert "groups" in status
    assert "core" in status["groups"]
    assert "modules" in status
    assert "qc" in status["modules"]
    assert "available" in status["modules"]["qc"]
    assert "summary" in status
    assert set(MODULE_REGISTRY).issubset(status["modules"])
    assert "optional_missing" in status["modules"]["subcluster"]
    assert "gseapy" not in status["modules"]["bulk_enrichment"]["missing"]
    assert {"hap.py", "som.py"} <= set(status["external_tools"])


def test_dependency_status_honors_configured_nextflow(tmp_path, monkeypatch):
    from config import Config
    from modules.platform.system_health import dependency_status

    nextflow = tmp_path / "nextflow"
    nextflow.write_text("#!/bin/sh\n", encoding="utf-8")
    nextflow.chmod(0o755)
    monkeypatch.setattr(Config, "WES_NEXTFLOW_BIN", str(nextflow))

    status = dependency_status()
    assert status["external_tools"]["nextflow"] == {
        "installed": True,
        "path": str(nextflow),
    }


def test_system_dependencies_api(tmp_path, monkeypatch):
    from app import create_app
    from config import Config

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    db_path = tmp_path / "test.db"
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "DB_PATH", str(db_path))

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        resp = client.get("/api/system/dependencies")

    assert resp.status_code == 200
    data = resp.get_json()
    assert "groups" in data
    assert "modules" in data
