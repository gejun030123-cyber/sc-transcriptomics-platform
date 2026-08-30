import os

import pytest

from config import Config
from database import get_conn
from models import AnalysisTask, ResultFile


@pytest.fixture
def client(test_project):
    from app import create_app

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as test_client:
        yield test_client


def test_task_results_expose_only_project_owned_image_attachments(test_project):
    from modules.ai_tools import _get_task_results

    plot_dir = os.path.join(Config.project_dir(test_project), "plots")
    os.makedirs(plot_dir, exist_ok=True)
    image_path = os.path.join(plot_dir, "pca.png")
    with open(image_path, "wb") as handle:
        handle.write(b"not-a-real-png-but-a-project-owned-file")

    task = AnalysisTask(project_id=test_project, module_name="bulk_pca", status="completed")
    task.save()
    image = ResultFile.create(
        task_id=task.id, project_id=test_project, file_type="png",
        category="pca", label="PCA 图", file_path=image_path,
    )
    table_path = os.path.join(plot_dir, "summary.csv")
    with open(table_path, "w", encoding="utf-8") as handle:
        handle.write("gene\nTP53\n")
    ResultFile.create(
        task_id=task.id, project_id=test_project, file_type="csv",
        category="table", label="摘要表", file_path=table_path,
    )

    result = _get_task_results(task.id)
    assert len(result["attachments"]) == 1
    attachment = result["attachments"][0]
    assert attachment["id"] == image.id
    assert attachment["url"] == f"/api/projects/{test_project}/result-file/{image.id}"
    assert result["result_files"][0]["id"] == image.id


def test_chat_response_and_history_keep_image_attachments(client, test_project, monkeypatch):
    from modules import ai_adapter
    from routes.chat import _chat_histories

    monkeypatch.setattr(Config, "AI_API_KEY", "test-key")
    monkeypatch.setattr(Config, "AI_API_TOKEN", "")
    _chat_histories.clear(test_project)
    attachments = [{
        "id": "figure-1", "task_id": "task-1", "label": "火山图",
        "category": "volcano", "type": "png",
        "url": f"/api/projects/{test_project}/result-file/figure-1",
    }]

    def fake_chat(messages, project_id=None):
        return {
            "reply": "这是差异表达火山图。",
            "tool_calls": [{"name": "get_task_results", "args": {"task_id": "task-1"}}],
            "proposed_tools": [],
            "attachments": attachments,
            "messages": messages + [{"role": "assistant", "content": "这是差异表达火山图。"}],
        }

    monkeypatch.setattr(ai_adapter, "chat", fake_chat)
    response = client.post("/api/chat", json={"project_id": test_project, "message": "展示结果图"})
    assert response.status_code == 200
    assert response.get_json()["attachments"] == attachments

    history = client.get(f"/api/chat/history/{test_project}").get_json()
    assert history["messages"][-1]["attachments"] == attachments
