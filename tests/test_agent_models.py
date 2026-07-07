# tests/test_agent_models.py
"""Agent 数据模型测试"""
import pytest
import json
from database import init_db


class TestAgentModels:
    """测试 Agent 模型的 CRUD 操作."""

    def test_init_db_creates_agent_tables(self):
        """验证数据库初始化创建所有 agent 表."""
        init_db()
        from database import get_conn
        conn = get_conn()
        try:
            tables = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
            table_names = {t['name'] for t in tables}
        finally:
            conn.close()

        expected = {'agent_sessions', 'agent_goals', 'agent_steps',
                     'analysis_branches', 'candidate_scores'}
        assert expected.issubset(table_names), f"缺少表: {expected - table_names}"

    def test_agent_session_create_and_get(self, test_project):
        """测试 AgentSession 的创建和查询."""
        from models import AgentSession
        session = AgentSession(
            project_id=test_project,
            status='active',
            goal_summary='测试目标',
        )
        session.save()

        fetched = AgentSession.get_by_id(session.id)
        assert fetched is not None
        assert fetched.project_id == test_project
        assert fetched.status == 'active'
        assert fetched.goal_summary == '测试目标'

    def test_agent_session_list_by_project(self, test_project):
        """测试按项目列出 sessions."""
        from models import AgentSession
        session = AgentSession(
            project_id=test_project,
            status='active',
        )
        session.save()
        sessions = AgentSession.get_by_project(test_project)
        assert len(sessions) > 0

    def test_agent_goal_create_and_get(self, test_project):
        """测试 AgentGoal 的创建和查询."""
        from models import AgentSession, AgentGoal
        session = AgentSession(project_id=test_project, status='active')
        session.save()

        goal = AgentGoal(
            session_id=session.id,
            project_id=test_project,
            goal_type='target_cluster_refinement',
            target_entity='microglia',
            goal_json=json.dumps({
                'target_cell_type': 'microglia',
                'positive_markers': ['P2RY12', 'TMEM119'],
            }),
        )
        goal.save()

        goals = AgentGoal.get_by_session(session.id)
        assert len(goals) == 1
        assert goals[0].goal_type == 'target_cluster_refinement'

    def test_agent_step_create_and_index(self, test_project):
        """测试 AgentStep 创建和 step_index 自动递增."""
        from models import AgentSession, AgentStep
        session = AgentSession(project_id=test_project, status='active')
        session.save()

        step1 = AgentStep(
            session_id=session.id,
            step_index=0,
            step_type='test',
            thought_summary='第一步',
            action_name='test_action',
        )
        step1.save()

        next_idx = AgentStep.next_index(session.id)
        assert next_idx == 1

        step2 = AgentStep(
            session_id=session.id,
            step_index=next_idx,
            step_type='test',
            thought_summary='第二步',
            action_name='test_action',
        )
        step2.save()

        steps = AgentStep.get_by_session(session.id)
        assert len(steps) == 2

    def test_analysis_branch_lifecycle(self, test_project):
        """测试 AnalysisBranch 生命周期：create → run → complete → accept."""
        from models import AnalysisBranch
        branch = AnalysisBranch(
            project_id=test_project,
            parent_adata_path='/tmp/test/test.h5ad',
            branch_name='test_candidate',
            purpose='测试候选',
            params_json=json.dumps({'clustering': {'resolutions': '1.0'}}),
            status='pending',
        )
        branch.save()

        assert branch.status == 'pending'
        assert branch.accepted == 0

        branch.mark_running()
        fetched = AnalysisBranch.get_by_id(branch.id)
        assert fetched.status == 'running'

        branch.mark_completed('/tmp/test/output.h5ad')
        fetched = AnalysisBranch.get_by_id(branch.id)
        assert fetched.status == 'completed'
        assert fetched.output_adata_path == '/tmp/test/output.h5ad'

        branch.accept()
        fetched = AnalysisBranch.get_by_id(branch.id)
        assert fetched.accepted == 1

    def test_analysis_branch_mark_failed(self, test_project):
        """测试失败分支保存错误信息."""
        from models import AnalysisBranch
        branch = AnalysisBranch(
            project_id=test_project,
            parent_adata_path='/tmp/test/test2.h5ad',
            branch_name='test_fail',
            status='pending',
        )
        branch.save()
        branch.mark_failed('测试错误信息')
        fetched = AnalysisBranch.get_by_id(branch.id)
        assert fetched.status == 'failed'
        assert '测试错误信息' in fetched.error_traceback

    def test_candidate_score_create_and_get(self, test_project):
        """测试 CandidateScore 的创建和查询."""
        from models import AnalysisBranch, CandidateScore
        branch = AnalysisBranch(
            project_id=test_project,
            parent_adata_path='/tmp/test/test3.h5ad',
            branch_name='test_score',
            status='completed',
        )
        branch.save()

        score = CandidateScore(
            branch_id=branch.id,
            project_id=test_project,
            evaluator_name='sc_cluster_signature',
            target_label='microglia:7',
            score_json=json.dumps({'total_score': 0.85}),
            total_score=0.85,
            recommendation='最佳 cluster: 7, 置信度: high',
        )
        score.save()

        scores = CandidateScore.get_by_branch(branch.id)
        assert len(scores) == 1
        assert scores[0].total_score == 0.85


class TestBranchIsolation:
    """测试未采纳 branch 不会污染主线."""

    def test_branch_task_excluded_from_latest_adata(self, test_project):
        """验证分支任务的 h5ad 不会成为项目 latest h5ad."""
        from models import AnalysisTask, AnalysisBranch
        from modules.ai_tools import resolve_current_adata_path
        from config import Config
        import os

        project_dir = Config.project_dir(test_project)
        intermediate_dir = os.path.join(project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)

        # 创建主线任务
        main_path = os.path.join(intermediate_dir, 'main_output.h5ad')
        with open(main_path, 'w') as f:
            f.write('mock')
        main_task = AnalysisTask(
            project_id=test_project,
            module_name='clustering',
            status='completed',
            output_adata_path=main_path,
            branch_id=None,
        )
        main_task.save()

        # 创建分支任务
        branch = AnalysisBranch(
            project_id=test_project,
            parent_adata_path=main_path,
            branch_name='test_branch',
            status='completed',
        )
        branch.save()

        branch_path = os.path.join(intermediate_dir, 'branch_output.h5ad')
        with open(branch_path, 'w') as f:
            f.write('mock')
        branch_task = AnalysisTask(
            project_id=test_project,
            module_name='clustering',
            status='completed',
            output_adata_path=branch_path,
            branch_id=branch.id,
        )
        branch_task.save()

        # 解析 current context 应该返回主线
        ctx = resolve_current_adata_path(test_project)
        assert ctx['path'] is not None
        assert 'main_output' in ctx['path']

    def test_project_latest_adata_excludes_branches(self, test_project):
        """验证 Project.get_latest_adata_path() 排除分支任务."""
        from models import Project, AnalysisTask, AnalysisBranch
        from config import Config
        import os

        project_dir = Config.project_dir(test_project)
        intermediate_dir = os.path.join(project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)

        # 主线任务
        main_path = os.path.join(intermediate_dir, 'main_v2.h5ad')
        with open(main_path, 'w') as f:
            f.write('mock')
        main_task = AnalysisTask(
            project_id=test_project,
            module_name='normalize',
            status='completed',
            output_adata_path=main_path,
            branch_id=None,
        )
        main_task.save()

        # 分支任务
        branch = AnalysisBranch(
            project_id=test_project,
            parent_adata_path=main_path,
            branch_name='test_iso',
            status='completed',
        )
        branch.save()

        branch_dir = Config.branch_dir(test_project, branch.id)
        os.makedirs(branch_dir, exist_ok=True)
        branch_path = os.path.join(branch_dir, 'output.h5ad')
        with open(branch_path, 'w') as f:
            f.write('mock')
        branch_task = AnalysisTask(
            project_id=test_project,
            module_name='clustering',
            status='completed',
            output_adata_path=branch_path,
            branch_id=branch.id,
        )
        branch_task.save()

        project = Project.get_by_id(test_project)
        latest = project.get_latest_adata_path()
        assert latest is not None
        assert 'branches' not in latest

    def test_run_analysis_uses_accepted_branch_by_default(self, test_project):
        """验证采纳 branch 后 resolve_current_adata_path + _run_analysis 使用 accepted branch."""
        from models import AnalysisBranch
        from modules.ai_tools import resolve_current_adata_path, _run_analysis
        from config import Config
        import os

        project_dir = Config.project_dir(test_project)
        os.makedirs(project_dir, exist_ok=True)

        # 先确认无数据时返回 none
        ctx_before = resolve_current_adata_path(test_project)
        assert ctx_before['source'] == 'none'

        # 创建 accepted branch
        branch = AnalysisBranch(
            project_id=test_project,
            parent_adata_path=os.path.join(project_dir, 'input.h5ad'),
            branch_name='accepted_candidate',
            status='pending',
        )
        branch.save()
        branch.mark_running()

        branch_dir = Config.branch_dir(test_project, branch.id)
        os.makedirs(branch_dir, exist_ok=True)
        accepted_output = os.path.join(branch_dir, 'output.h5ad')
        with open(accepted_output, 'w') as f:
            f.write('mock')
        branch.mark_completed(accepted_output)
        branch.accept()

        # 1. 验证 context resolver 返回 accepted branch
        ctx = resolve_current_adata_path(test_project)
        assert ctx['source'] == 'accepted_branch'
        assert ctx['path'] == accepted_output

        # 2. 验证 _run_analysis 空 input 时使用 accepted branch 作为默认输入
        # 用 try/except 包裹，因为 submit_task 会被触发但我们只关心 input_source
        import worker as worker_mod
        captured = {}
        original_submit = worker_mod.submit_task

        def mock_submit(task_id, p_id, mod_name, params, proj_dir, input_path):
            captured['input_path'] = input_path
            # 不真正提交任务

        try:
            worker_mod.submit_task = mock_submit
            result = _run_analysis(
                {"module_name": "clustering", "params": {"resolutions": "1.0"}},
                test_project,
            )
            assert result.get('input_source') == 'accepted_branch'
            assert captured.get('input_path') == accepted_output
        finally:
            worker_mod.submit_task = original_submit

    def test_run_analysis_input_source_returned(self, test_project):
        """验证 _run_analysis 返回 input_path 和 input_source."""
        from modules.ai_tools import _run_analysis
        from config import Config
        import os

        project_dir = Config.project_dir(test_project)
        intermediate_dir = os.path.join(project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        test_input = os.path.join(intermediate_dir, 'test_input.h5ad')
        with open(test_input, 'w') as f:
            f.write('mock')

        # 手动指定 input_path
        result = _run_analysis(
            {"module_name": "clustering", "input_path": test_input},
            test_project,
        )

        assert 'input_path' in result
        assert 'input_source' in result
        assert result['input_source'] == 'manual'


class TestConfigValidation:
    """测试路径验证安全规则."""

    def test_validate_pid_rejects_invalid(self):
        """测试非法 project ID 被拒绝."""
        from config import Config
        with pytest.raises(ValueError):
            Config._validate_pid('../../../etc')

    def test_validate_pid_accepts_valid(self):
        """测试合法 project ID 通过验证."""
        from config import Config
        assert Config._validate_pid('abc123') == 'abc123'

    def test_branch_dir_path(self, tmp_path, monkeypatch):
        """测试分支目录路径在项目目录内."""
        from config import Config
        monkeypatch.setattr(Config, 'DATA_DIR', str(tmp_path))

        proj_dir = Config.project_dir('test_validation')
        branch_dir = Config.branch_dir('test_validation', 'branch_001')

        import os
        os.makedirs(branch_dir, exist_ok=True)

        real_proj = os.path.realpath(proj_dir)
        real_branch = os.path.realpath(branch_dir)
        assert real_branch.startswith(real_proj + os.sep) or real_branch == real_proj
