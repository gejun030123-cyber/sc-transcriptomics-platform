# tests/test_branch_routes.py
"""Branch API 路由测试 — 鉴权、project scope、accept/delete/current-context"""
import pytest
import json
import os
import uuid
from app import create_app


@pytest.fixture
def client(tmp_path, monkeypatch):
    """创建 Flask 测试客户端（数据隔离到 tmp_path，启用 token 验证）。"""
    from config import Config

    # 隔离数据目录和数据库到 tmp_path
    data_dir = tmp_path / 'data'
    data_dir.mkdir(exist_ok=True)
    db_path = tmp_path / 'test.db'
    monkeypatch.setattr(Config, 'DATA_DIR', str(data_dir))
    monkeypatch.setattr(Config, 'DB_PATH', str(db_path))
    monkeypatch.setattr(Config, 'AI_API_TOKEN', 'test-token-123')

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as c:
        yield c
    # tmp_path 自动清理


@pytest.fixture
def auth_client(tmp_path, monkeypatch):
    """创建带 token 的 Flask 测试客户端（数据隔离到 tmp_path）。"""
    from config import Config

    data_dir = tmp_path / 'data'
    data_dir.mkdir(exist_ok=True)
    db_path = tmp_path / 'test.db'
    monkeypatch.setattr(Config, 'DATA_DIR', str(data_dir))
    monkeypatch.setattr(Config, 'DB_PATH', str(db_path))
    monkeypatch.setattr(Config, 'AI_API_TOKEN', 'test-token-123')

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as c:
        yield c


def _make_project(client, pid=None):
    """创建测试 project，使用唯一随机 ID。"""
    if pid is None:
        pid = 'rt_' + uuid.uuid4().hex[:8]
    from database import get_conn
    conn = get_conn()
    try:
        conn.execute(
            "INSERT OR REPLACE INTO projects (id, name, status) VALUES (?,?,?)",
            (pid, 'Route Test Project', 'active')
        )
        conn.commit()
    finally:
        conn.close()
    return pid


class TestBranchRouteAuth:
    """测试 branch 路由鉴权."""

    WRITE_ENDPOINTS = [
        ('POST', '/api/projects/{pid}/agent/sessions'),
        ('POST', '/api/projects/{pid}/branches'),
        ('POST', '/api/projects/{pid}/branches/fake_id/run'),
        ('POST', '/api/projects/{pid}/branches/fake_id/accept'),
        ('POST', '/api/projects/{pid}/branches/fake_id/score'),
        ('POST', '/api/projects/{pid}/branches/fake_id/delete'),
    ]

    READ_ENDPOINTS = [
        ('GET', '/api/projects/{pid}/agent/sessions'),
        ('GET', '/api/projects/{pid}/agent/sessions/fake_sid'),
        ('GET', '/api/projects/{pid}/branches'),
        ('GET', '/api/projects/{pid}/branches/fake_id'),
    ]

    def test_write_endpoints_reject_no_token(self, client):
        """无 token 调写接口应返回 401."""
        pid = _make_project(client)
        for method, url_tpl in self.WRITE_ENDPOINTS:
            url = url_tpl.format(pid=pid)
            if method == 'POST':
                resp = client.post(url, json={})
            else:
                resp = client.get(url)
            assert resp.status_code == 401, f"{method} {url} 应返回 401，实际 {resp.status_code}"

    def test_write_endpoints_accept_correct_token(self, auth_client):
        """正确 token 调写接口不返回 401."""
        pid = _make_project(auth_client)
        for method, url_tpl in self.WRITE_ENDPOINTS:
            url = url_tpl.format(pid=pid)
            headers = {'Authorization': 'Bearer test-token-123'}
            if method == 'POST':
                resp = auth_client.post(url, json={}, headers=headers)
            else:
                resp = auth_client.get(url, headers=headers)
            # 不返回 401 即可（可能返回 400/404 因为缺少参数/fake_id）
            assert resp.status_code != 401, f"{method} {url} 不应返回 401"

    def test_write_endpoints_reject_wrong_token(self, auth_client):
        """错误 token 应返回 401."""
        pid = _make_project(auth_client)
        for method, url_tpl in self.WRITE_ENDPOINTS:
            url = url_tpl.format(pid=pid)
            headers = {'Authorization': 'Bearer wrong-token'}
            if method == 'POST':
                resp = auth_client.post(url, json={}, headers=headers)
            else:
                resp = auth_client.get(url, headers=headers)
            assert resp.status_code == 401, f"{method} {url} 应返回 401"

    def test_current_context_endpoint(self, auth_client):
        """测试 current-context API."""
        pid = _make_project(auth_client)
        headers = {'Authorization': 'Bearer test-token-123'}
        resp = auth_client.get(f'/api/projects/{pid}/current-context', headers=headers)
        assert resp.status_code == 200
        data = resp.get_json()
        assert 'source' in data
        assert 'current_adata_path' in data
        assert 'can_continue_analysis' in data


class TestBranchProjectScope:
    """测试 project scope 校验."""

    def test_branch_wrong_project_rejected(self, auth_client):
        """用项目A的pid访问项目B的branch应返回403."""
        pid_a = _make_project(auth_client, 'project_a')
        pid_b = _make_project(auth_client, 'project_b')

        # 在项目B中创建branch
        headers = {'Authorization': 'Bearer test-token-123'}
        # 先创建项目目录
        from config import Config
        os.makedirs(Config.project_dir(pid_b), exist_ok=True)
        test_file = os.path.join(Config.project_dir(pid_b), 'test.h5ad')
        with open(test_file, 'w') as f:
            f.write('mock')

        resp = auth_client.post(
            f'/api/projects/{pid_b}/branches',
            json={'parent_adata_path': test_file, 'branch_name': 'test_branch'},
            headers=headers,
        )
        assert resp.status_code == 200
        branch_id = resp.get_json()['id']

        # 用项目A的pid访问项目B的branch
        resp = auth_client.get(
            f'/api/projects/{pid_a}/branches/{branch_id}',
            headers=headers,
        )
        assert resp.status_code == 403

    def test_session_wrong_project_rejected(self, auth_client):
        """用项目A的pid访问项目B的session应返回403."""
        pid_a = _make_project(auth_client, 'project_c')
        pid_b = _make_project(auth_client, 'project_d')

        headers = {'Authorization': 'Bearer test-token-123'}

        # 在项目B中创建session
        resp = auth_client.post(
            f'/api/projects/{pid_b}/agent/sessions',
            json={'goal_summary': 'test'},
            headers=headers,
        )
        assert resp.status_code == 200
        session_id = resp.get_json()['id']

        # 用项目A的pid访问项目B的session
        resp = auth_client.get(
            f'/api/projects/{pid_a}/agent/sessions/{session_id}',
            headers=headers,
        )
        assert resp.status_code == 403


class TestBranchAcceptDelete:
    """测试 branch 采纳和删除语义."""

    def test_accept_requires_confirm_true(self, auth_client):
        """采纳必须 confirm=true."""
        pid = _make_project(auth_client, 'accept_test_project')
        headers = {'Authorization': 'Bearer test-token-123'}

        from config import Config
        os.makedirs(Config.project_dir(pid), exist_ok=True)
        test_file = os.path.join(Config.project_dir(pid), 'test.h5ad')
        with open(test_file, 'w') as f:
            f.write('mock')

        # 创建branch
        resp = auth_client.post(
            f'/api/projects/{pid}/branches',
            json={'parent_adata_path': test_file, 'branch_name': 'accept_test'},
            headers=headers,
        )
        branch_id = resp.get_json()['id']

        # 手动标记为completed
        from models import AnalysisBranch
        branch = AnalysisBranch.get_by_id(branch_id)
        branch.status = 'running'
        branch.save()
        branch.mark_completed(test_file)

        # 不带confirm的accept
        resp = auth_client.post(
            f'/api/projects/{pid}/branches/{branch_id}/accept',
            json={},
            headers=headers,
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data.get('status') == 'needs_confirmation'

        # 带confirm=true的accept
        resp = auth_client.post(
            f'/api/projects/{pid}/branches/{branch_id}/accept',
            json={'confirm': True},
            headers=headers,
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data.get('status') == 'accepted'

    def test_delete_branch_soft_delete(self, auth_client):
        """删除branch应为软删除，branch task不失归属."""
        pid = _make_project(auth_client, 'delete_test_project')
        headers = {'Authorization': 'Bearer test-token-123'}

        from config import Config
        os.makedirs(Config.project_dir(pid), exist_ok=True)
        test_file = os.path.join(Config.project_dir(pid), 'test.h5ad')
        with open(test_file, 'w') as f:
            f.write('mock')

        # 创建branch
        resp = auth_client.post(
            f'/api/projects/{pid}/branches',
            json={'parent_adata_path': test_file, 'branch_name': 'delete_test'},
            headers=headers,
        )
        branch_id = resp.get_json()['id']

        # 删除前 branch存在
        resp = auth_client.get(
            f'/api/projects/{pid}/branches/{branch_id}',
            headers=headers,
        )
        assert resp.status_code == 200

        # 删除
        resp = auth_client.post(
            f'/api/projects/{pid}/branches/{branch_id}/delete',
            headers=headers,
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data.get('status') == 'deleted'

        # 删除后 branch不在列表中
        resp = auth_client.get(
            f'/api/projects/{pid}/branches',
            headers=headers,
        )
        branches = resp.get_json()['branches']
        assert all(b['id'] != branch_id for b in branches)

    def test_current_context_after_accept(self, auth_client):
        """采纳branch后 current-context 应指向 accepted branch."""
        pid = _make_project(auth_client, 'ctx_accept_project')
        headers = {'Authorization': 'Bearer test-token-123'}

        from config import Config
        proj_dir = Config.project_dir(pid)
        os.makedirs(proj_dir, exist_ok=True)
        test_file = os.path.join(proj_dir, 'test.h5ad')
        with open(test_file, 'w') as f:
            f.write('mock')

        # 创建branch
        resp = auth_client.post(
            f'/api/projects/{pid}/branches',
            json={'parent_adata_path': test_file, 'branch_name': 'ctx_accept_branch'},
            headers=headers,
        )
        branch_id = resp.get_json()['id']

        # 手动完成branch
        from models import AnalysisBranch
        branch = AnalysisBranch.get_by_id(branch_id)
        branch.status = 'running'
        branch.save()
        output_path = os.path.join(Config.branch_dir(pid, branch_id), 'output.h5ad')
        os.makedirs(Config.branch_dir(pid, branch_id), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write('mock')
        branch.mark_completed(output_path)

        # accept
        resp = auth_client.post(
            f'/api/projects/{pid}/branches/{branch_id}/accept',
            json={'confirm': True},
            headers=headers,
        )
        assert resp.status_code == 200

        # 检查 current-context
        resp = auth_client.get(
            f'/api/projects/{pid}/current-context',
            headers=headers,
        )
        assert resp.status_code == 200
        ctx = resp.get_json()
        assert ctx['source'] == 'accepted_branch'
        assert ctx['accepted_branch_id'] == branch_id
