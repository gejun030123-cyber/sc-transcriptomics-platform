# tests/test_agent_jobs.py
"""Agent Job 异步任务测试 — 成功/失败/路由认证"""
import pytest
import json
import os
import time
from app import create_app


@pytest.fixture
def agent_client(tmp_path, monkeypatch):
    """创建 Flask 测试客户端（数据隔离 + token）。"""
    from config import Config
    data_dir = tmp_path / 'data'
    data_dir.mkdir(exist_ok=True)
    monkeypatch.setattr(Config, 'DATA_DIR', str(data_dir))
    monkeypatch.setattr(Config, 'DB_PATH', str(tmp_path / 'test.db'))
    monkeypatch.setattr(Config, 'AI_API_TOKEN', 'test-token-123')

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as c:
        yield c


def _make_project_and_goal(pid='job_test_proj'):
    """创建测试项目 + session + goal + 有效 checkpoint 文件。返回 (pid, goal, valid_path)."""
    from database import get_conn
    from config import Config

    conn = get_conn()
    try:
        conn.execute(
            "INSERT OR REPLACE INTO projects (id, name, status) VALUES (?,?,?)",
            (pid, 'Job Test', 'active')
        )
        conn.commit()
    finally:
        conn.close()

    from models import AgentSession, AgentGoal

    session = AgentSession(project_id=pid, status='active', goal_summary='test sweep')
    session.save()

    goal = AgentGoal(
        session_id=session.id,
        project_id=pid,
        goal_type='target_cluster_refinement',
        target_entity='microglia',
        goal_json='{"target_cell_type":"microglia"}',
    )
    goal.save()

    proj_dir = Config.project_dir(pid)
    os.makedirs(proj_dir, exist_ok=True)
    valid_path = os.path.join(proj_dir, 'test.h5ad')
    with open(valid_path, 'w') as f:
        f.write('mock')

    return pid, goal, valid_path


class TestAgentJobSubmission:
    """测试 job 提交."""

    def test_submit_returns_job_id(self, agent_client):
        """提交 sweep 立即返回 job_id."""
        pid, goal, valid_path = _make_project_and_goal()

        headers = {'Authorization': 'Bearer test-token-123'}
        resp = agent_client.post(
            f'/api/projects/{pid}/agent/sessions',
            json={'goal_summary': 'test'},
            headers=headers,
        )

        # 直接用函数测试
        from modules.agent_orchestrator import run_parameter_sweep
        candidates = [{'name': 'c1', 'modules': ['clustering'], 'params': {}}]
        result = run_parameter_sweep(goal.id, candidates, valid_path, pid)

        assert result['status'] == 'submitted'
        assert 'job_id' in result
        assert 'poll_url' in result
        assert result['n_candidates'] == 1

    def test_session_id_set_on_job(self, agent_client):
        """Job 创建时包含 session_id."""
        pid, goal, valid_path = _make_project_and_goal()

        from modules.agent_orchestrator import run_parameter_sweep
        from models import AgentJob
        candidates = [{'name': 'c1', 'modules': ['clustering'], 'params': {}}]
        result = run_parameter_sweep(goal.id, candidates, valid_path, pid)

        job = AgentJob.get_by_id(result['job_id'])
        assert job.session_id is not None
        assert job.session_id == goal.session_id


class TestAgentJobRoute:
    """测试 job 路由."""

    def test_job_status_route(self, agent_client):
        """GET job status 返回正确."""
        pid, goal, valid_path = _make_project_and_goal()

        from modules.agent_orchestrator import run_parameter_sweep
        candidates = [{'name': 'c1', 'modules': ['clustering'], 'params': {}}]
        result = run_parameter_sweep(goal.id, candidates, valid_path, pid)

        headers = {'Authorization': 'Bearer test-token-123'}
        resp = agent_client.get(
            f'/api/projects/{pid}/agent/jobs/{result["job_id"]}',
            headers=headers,
        )
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['id'] == result['job_id']
        assert data['job_type'] == 'parameter_sweep'
        assert data['status'] in ('pending', 'running', 'completed', 'failed')

    def test_job_status_rejects_no_token(self, agent_client):
        """无 token 调 job status 返回 401."""
        pid, goal, valid_path = _make_project_and_goal('job_auth_test')
        from modules.agent_orchestrator import run_parameter_sweep
        result = run_parameter_sweep(goal.id, [{'name': 'c1', 'modules': ['clustering'], 'params': {}}], valid_path, pid)

        resp = agent_client.get(f'/api/projects/{pid}/agent/jobs/{result["job_id"]}')
        assert resp.status_code == 401

    def test_job_status_cross_project_rejected(self, agent_client):
        """跨项目访问 job 返回 403."""
        pid1, goal1, path1 = _make_project_and_goal('job_proj_a')
        pid2, goal2, path2 = _make_project_and_goal('job_proj_b')

        from modules.agent_orchestrator import run_parameter_sweep
        result = run_parameter_sweep(goal1.id, [{'name': 'c1', 'modules': ['clustering'], 'params': {}}], path1, pid1)

        headers = {'Authorization': 'Bearer test-token-123'}
        resp = agent_client.get(
            f'/api/projects/{pid2}/agent/jobs/{result["job_id"]}',
            headers=headers,
        )
        assert resp.status_code == 403


class TestAgentJobBackground:
    """测试后台 job 执行."""

    def test_job_session_id_is_set_by_submit(self, test_project):
        """submit_sweep_job 从 goal 读取 session_id 并写入 job."""
        from models import AgentSession, AgentGoal, AgentJob
        from config import Config
        import os

        pid = test_project
        session = AgentSession(project_id=pid, status='active', goal_summary='test')
        session.save()

        goal = AgentGoal(
            session_id=session.id, project_id=pid,
            goal_type='target_cluster_refinement', target_entity='microglia',
            goal_json='{"target_cell_type":"microglia"}',
        )
        goal.save()

        valid_path = os.path.join(Config.project_dir(pid), 'test.h5ad')
        with open(valid_path, 'w') as f:
            f.write('mock')

        from modules.agent_orchestrator import run_parameter_sweep
        candidates = [{'name': 'c1', 'modules': ['clustering'], 'params': {}}]
        result = run_parameter_sweep(goal.id, candidates, valid_path, pid)

        job = AgentJob.get_by_id(result['job_id'])
        assert job.session_id is not None
        assert job.session_id == session.id

    def test_job_marks_failed_on_error(self, test_project, monkeypatch):
        """后台线程异常时 job 进入 failed 而非永久 running."""
        from models import AgentSession, AgentGoal, AgentJob
        from config import Config
        import os

        pid = test_project
        session = AgentSession(project_id=pid, status='active', goal_summary='test')
        session.save()

        goal = AgentGoal(
            session_id=session.id, project_id=pid,
            goal_type='target_cluster_refinement', target_entity='microglia',
            goal_json='{"target_cell_type":"microglia"}',
        )
        goal.save()

        valid_path = os.path.join(Config.project_dir(pid), 'test.h5ad')
        with open(valid_path, 'w') as f:
            f.write('mock')

        from modules import MODULE_REGISTRY

        class BadModule:
            DISPLAY_NAME = 'Bad Module'
            DESCRIPTION = 'Always fails'
            def __init__(self, project_dir=None, params=None, progress_callback=None):
                pass
            def run(self, input_path):
                raise RuntimeError('simulated module crash')

        monkeypatch.setitem(MODULE_REGISTRY, 'clustering', BadModule)

        from modules.agent_orchestrator import run_parameter_sweep
        candidates = [{'name': 'fail_c1', 'modules': ['clustering'], 'params': {}}]
        result = run_parameter_sweep(goal.id, candidates, valid_path, pid)

        job_id = result['job_id']

        for _ in range(50):
            job = AgentJob.get_by_id(job_id)
            if job.status in ('completed', 'failed'):
                break
            time.sleep(0.2)

        job = AgentJob.get_by_id(job_id)
        assert job.status != 'running', f"job should not be stuck running, got {job.status}"
        assert job.finished_at is not None
