import json


def test_chat_config_route_hides_api_key(tmp_path, monkeypatch):
    from app import create_app
    from config import Config

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    db_path = tmp_path / "test.db"
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "DB_PATH", str(db_path))
    monkeypatch.setattr(Config, "AI_API_KEY", "secret-value")
    monkeypatch.setattr(Config, "AI_API_URL", "https://api.deepseek.com/anthropic")
    monkeypatch.setattr(Config, "AI_MODEL", "deepseek-chat")
    monkeypatch.setattr(Config, "AI_API_TOKEN", "")

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        resp = client.get("/api/chat/config")

    assert resp.status_code == 200
    data = resp.get_json()
    assert data["configured"] is True
    assert data["provider"] == "deepseek-anthropic"
    assert data["model"] == "deepseek-chat"
    assert "secret-value" not in json.dumps(data)


def test_chat_endpoint_rejects_invalid_project_id(tmp_path, monkeypatch):
    from app import create_app
    from config import Config

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    db_path = tmp_path / "test.db"
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "DB_PATH", str(db_path))
    monkeypatch.setattr(Config, "AI_API_KEY", "secret-value")
    monkeypatch.setattr(Config, "AI_API_TOKEN", "")

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        resp = client.post(
            "/api/chat",
            json={"message": "hello", "project_id": "../bad"},
        )

    assert resp.status_code == 400
    assert "Invalid project ID" in resp.get_json()["error"]


def test_chat_endpoint_uses_form_context_without_storing_it_as_visible_history(tmp_path, monkeypatch):
    from app import create_app
    from config import Config
    from modules import ai_adapter
    from routes.chat import _chat_histories

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "DB_PATH", str(tmp_path / "test.db"))
    monkeypatch.setattr(Config, "AI_API_KEY", "secret-value")
    monkeypatch.setattr(Config, "AI_API_TOKEN", "")
    _chat_histories.clear("context-project")

    captured = {}

    def fake_chat(messages, project_id=None):
        captured["message"] = messages[-1]["content"]
        return {
            "reply": "建议先检查分组。",
            "tool_calls": [],
            "proposed_tools": [],
            "messages": messages + [{"role": "assistant", "content": "建议先检查分组。"}],
        }

    monkeypatch.setattr(ai_adapter, "chat", fake_chat)
    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        response = client.post("/api/chat", json={
            "project_id": "context-project",
            "message": "推荐当前参数",
            "context": {"module": "bulk_deg", "params": {"method": "welch_ttest"}},
        })
        history = client.get("/api/chat/history/context-project").get_json()

    assert response.status_code == 200
    assert "bulk_deg" in captured["message"]
    assert history["messages"][0]["content"] == "推荐当前参数"


def test_deepseek_anthropic_endpoint_is_messages_url(monkeypatch):
    from config import Config
    from modules.ai_adapter import _anthropic_messages_endpoint, get_ai_config_status

    monkeypatch.setattr(Config, "AI_API_URL", "https://api.deepseek.com/anthropic")
    monkeypatch.setattr(Config, "AI_API_KEY", "secret-value")
    monkeypatch.setattr(Config, "AI_MODEL", "deepseek-chat")

    assert _anthropic_messages_endpoint() == "https://api.deepseek.com/anthropic/v1/messages"
    status = get_ai_config_status()
    assert status["provider"] == "deepseek-anthropic"
    assert status["mode"] == "anthropic_messages"
