from config import Config
from modules.ai_config import get_effective_ai_config
import pytest


@pytest.fixture
def client(test_project):
    from app import create_app

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as test_client:
        yield test_client


def test_ai_settings_page_and_masked_status(client, test_project, monkeypatch):
    monkeypatch.setattr(Config, "AI_API_KEY", "env-secret-key")
    monkeypatch.setattr(Config, "AI_API_URL", "https://env.example/v1")
    monkeypatch.setattr(Config, "AI_MODEL", "env-model")
    monkeypatch.setattr(Config, "AI_API_TOKEN", "")

    page = client.get("/settings/ai")
    assert page.status_code == 200
    assert "AI API 设置" in page.get_data(as_text=True)

    response = client.get("/api/settings/ai")
    assert response.status_code == 200
    data = response.get_json()
    assert data["configured"] is True
    assert data["api_key_masked"] != "env-secret-key"
    assert "env-secret-key" not in response.get_data(as_text=True)


def test_ai_settings_save_changes_runtime_without_restart(client, test_project, monkeypatch):
    monkeypatch.setattr(Config, "AI_API_KEY", "")
    monkeypatch.setattr(Config, "AI_API_URL", "https://env.example/v1")
    monkeypatch.setattr(Config, "AI_MODEL", "env-model")
    monkeypatch.setattr(Config, "AI_API_TOKEN", "")

    response = client.post("/api/settings/ai", json={
        "provider": "openai",
        "api_url": "https://api.example.test/v1",
        "model": "new-model",
        "api_key": "new-secret-key",
    })
    assert response.status_code == 200
    assert response.get_json()["settings"]["api_key_masked"] != "new-secret-key"

    effective = get_effective_ai_config()
    assert effective == {
        "provider": "openai",
        "api_url": "https://api.example.test/v1",
        "model": "new-model",
        "api_key": "new-secret-key",
    }

    # Blank API key preserves the current secret while changing the model.
    response = client.post("/api/settings/ai", json={
        "provider": "auto",
        "api_url": "https://api.example.test/v1",
        "model": "second-model",
        "api_key": "",
    })
    assert response.status_code == 200
    assert get_effective_ai_config()["api_key"] == "new-secret-key"
    assert get_effective_ai_config()["model"] == "second-model"


def test_ai_settings_rejects_invalid_url_and_can_reset(client, test_project, monkeypatch):
    monkeypatch.setattr(Config, "AI_API_KEY", "env-key")
    monkeypatch.setattr(Config, "AI_API_URL", "https://env.example/v1")
    monkeypatch.setattr(Config, "AI_MODEL", "env-model")
    monkeypatch.setattr(Config, "AI_API_TOKEN", "")

    response = client.post("/api/settings/ai", json={
        "api_url": "not-a-url", "model": "model", "api_key": "key"
    })
    assert response.status_code == 400

    client.post("/api/settings/ai", json={
        "api_url": "https://api.example.test/v1", "model": "new", "api_key": "key"
    })
    response = client.post("/api/settings/ai/reset")
    assert response.status_code == 200
    assert get_effective_ai_config()["api_url"] == "https://env.example/v1"
    assert get_effective_ai_config()["api_key"] == "env-key"


def test_ai_settings_connection_test_does_not_save(client, test_project, monkeypatch):
    monkeypatch.setattr(Config, "AI_API_KEY", "saved-key")
    monkeypatch.setattr(Config, "AI_API_URL", "https://env.example/v1")
    monkeypatch.setattr(Config, "AI_MODEL", "env-model")
    monkeypatch.setattr(Config, "AI_API_TOKEN", "")

    class Response:
        status_code = 200
        text = "{}"

    calls = {}

    def fake_get(url, **kwargs):
        calls["url"] = url
        calls["headers"] = kwargs["headers"]
        return Response()

    monkeypatch.setattr("requests.get", fake_get)
    response = client.post("/api/settings/ai/test", json={
        "provider": "openai",
        "api_url": "https://api.example.test/v1",
        "model": "test-model",
        "api_key": "temporary-key",
    })
    assert response.status_code == 200
    assert calls["url"] == "https://api.example.test/v1/models"
    assert calls["headers"]["Authorization"] == "Bearer temporary-key"
    assert get_effective_ai_config()["api_key"] == "saved-key"
