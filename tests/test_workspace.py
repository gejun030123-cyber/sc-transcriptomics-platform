"""AI 主工作台的页面契约测试。"""


def test_workspace_page_is_project_scoped_and_reserves_wes_atac(test_project):
    from app import create_app

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        response = client.get(f"/projects/{test_project}/workspace")

    assert response.status_code == 200
    html = response.get_data(as_text=True)
    assert "AI 分析执行台" in html
    assert "WES 外显子组" in html
    assert "Bulk ATAC" in html
    assert "规划中" in html
    assert f"/projects/{test_project}/sc-analysis" in html
    assert f"/projects/{test_project}/bulk-analysis" in html
    assert f'window.WORKSPACE_CONFIG = {{ projectId: "{test_project}" }}' in html


def test_workspace_assets_are_served(test_project):
    from app import create_app

    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as client:
        css = client.get("/static/css/workspace.css")
        js = client.get("/static/js/workspace.js")

    assert css.status_code == 200
    assert ".workspace-grid" in css.get_data(as_text=True)
    assert js.status_code == 200
    assert "workspaceSend" in js.get_data(as_text=True)
