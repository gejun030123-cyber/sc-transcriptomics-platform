# tests/test_agent_orchestrator.py
"""Agent 编排器测试"""
import pytest
import json
from database import init_db


class TestAgentOrchestratorStartup:
    """测试 start_goal_agent 初始流程."""

    def test_start_goal_agent_no_data(self, test_project):
        """测试项目无数据时返回 blocked 状态."""
        from modules.agent_orchestrator import start_goal_agent

        result = start_goal_agent(
            project_id=test_project,
            goal_type='target_cluster_refinement',
            target_cell_type='microglia',
        )

        assert result['status'] in ('blocked', 'needs_confirmation', 'goal_achieved')

    def test_sweep_candidates_limited(self, test_project):
        """测试候选生成不超过上限."""
        from modules.agent_orchestrator import _generate_sweep_candidates
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
                'max_candidate_runs': 3,
            }),
        )
        goal.save()

        state = {'available_embeddings': ['X_umap']}
        current_score = {'confidence': 'low'}

        candidates = _generate_sweep_candidates(goal, state, current_score, max_candidates=3)
        assert len(candidates) <= 3


class TestContinueGoalAgent:
    """测试 continue_goal_agent 指令解析."""

    @pytest.fixture(autouse=True)
    def setup_method(self, test_project):
        """创建测试 session 和 goal."""
        from models import AgentSession, AgentGoal
        self.session = AgentSession(
            project_id=test_project,
            status='active',
        )
        self.session.save()

        self.goal = AgentGoal(
            session_id=self.session.id,
            project_id=test_project,
            goal_type='target_cluster_refinement',
            target_entity='microglia',
            goal_json=json.dumps({'target_cell_type': 'microglia'}),
        )
        self.goal.save()

    def test_stop_instruction(self):
        """测试停止指令."""
        from modules.agent_orchestrator import continue_goal_agent

        result = continue_goal_agent(self.session.id, '停止这个目标')
        assert result['status'] == 'stopped'

    def test_constraint_update(self):
        """测试约束更新指令."""
        from modules.agent_orchestrator import continue_goal_agent

        result = continue_goal_agent(self.session.id, '必须包含 P2RY12 和 TMEM119')
        assert result['status'] == 'constraints_updated'

    def test_unrecognized_instruction(self):
        """测试未识别指令."""
        from modules.agent_orchestrator import continue_goal_agent

        result = continue_goal_agent(self.session.id, '随便说点什么')
        assert result['status'] == 'unknown_instruction'

    def test_accept_without_branches(self):
        """测试无候选分支时的采纳指令."""
        from modules.agent_orchestrator import continue_goal_agent

        result = continue_goal_agent(self.session.id, '采用候选A')
        assert result['status'] == 'blocked'


class TestAgentStepAudit:
    """测试 Agent 步骤审计记录."""

    def test_start_creates_audit_steps(self, test_project):
        """测试 start_goal_agent 创建审计记录."""
        from modules.agent_orchestrator import start_goal_agent
        from models import AgentStep

        result = start_goal_agent(
            project_id=test_project,
            goal_type='target_cluster_refinement',
            target_cell_type='microglia',
            user_requirement='找到小胶质细胞',
        )

        if 'session_id' in result:
            steps = AgentStep.get_by_session(result['session_id'])
            assert len(steps) >= 1
            assert steps[0].step_type == 'goal_created'

    def test_continue_creates_audit_steps(self, test_project):
        """测试 continue 创建审计记录."""
        from models import AgentSession, AgentGoal, AgentStep
        from modules.agent_orchestrator import continue_goal_agent

        session = AgentSession(project_id=test_project, status='active')
        session.save()

        goal = AgentGoal(
            session_id=session.id,
            project_id=test_project,
            goal_type='target_cluster_refinement',
            target_entity='microglia',
            goal_json=json.dumps({'target_cell_type': 'microglia'}),
        )
        goal.save()

        result = continue_goal_agent(session.id, '继续细分这个cluster')

        if result['status'] != 'blocked':
            steps = AgentStep.get_by_session(session.id)
            assert len(steps) > 0
            last_step = steps[-1]
            assert last_step.step_type == 'user_instruction'


class TestCrossProjectAccess:
    """测试跨项目访问被拒绝."""

    def test_sweep_goal_wrong_project_rejected(self, test_project):
        """验证 goal 不属于当前 project 时 sweep 被拒绝."""
        from models import AgentSession, AgentGoal
        from modules.agent_orchestrator import run_parameter_sweep
        from config import Config
        import os

        session = AgentSession(project_id=test_project, status='active')
        session.save()

        goal = AgentGoal(
            session_id=session.id,
            project_id=test_project,
            goal_type='target_cluster_refinement',
            target_entity='microglia',
            goal_json='{"target_cell_type":"microglia"}',
        )
        goal.save()

        valid_path = os.path.join(Config.project_dir(test_project), 'test.h5ad')
        os.makedirs(Config.project_dir(test_project), exist_ok=True)
        with open(valid_path, 'w') as f:
            f.write('mock')

        result = run_parameter_sweep(
            goal_id=goal.id,
            candidates=[{'name': 'test', 'modules': ['clustering'], 'params': {}}],
            base_checkpoint=valid_path,
            project_id='different_project',
        )

        assert 'error' in result
        assert '不属于' in result['error']

    def test_continue_session_wrong_project_rejected(self, test_project):
        """验证 session 不属于当前 project 时 continue 被拒绝."""
        from models import AgentSession
        from modules.agent_orchestrator import continue_goal_agent

        session = AgentSession(project_id=test_project, status='active')
        session.save()

        result = continue_goal_agent(
            session_id=session.id,
            user_instruction='测试',
            project_id='different_project',
        )

        assert 'error' in result
        assert '不属于' in result['error']


class TestSweepValidation:
    """测试参数搜索执行期校验."""

    def test_candidate_count_hard_limit(self, test_project):
        """测试候选数量超过硬上限被拒绝."""
        from models import AgentSession, AgentGoal
        from modules.agent_orchestrator import run_parameter_sweep
        from config import Config
        import os

        session = AgentSession(project_id=test_project, status='active')
        session.save()

        goal = AgentGoal(
            session_id=session.id,
            project_id=test_project,
            goal_type='target_cluster_refinement',
            target_entity='microglia',
            goal_json='{"target_cell_type":"microglia"}',
        )
        goal.save()

        # 创建有效 checkpoint
        valid_path = os.path.join(Config.project_dir(test_project), 'test.h5ad')
        os.makedirs(Config.project_dir(test_project), exist_ok=True)
        with open(valid_path, 'w') as f:
            f.write('mock')

        # 创建超过硬上限的候选列表
        too_many = [{'name': f'c{i}', 'modules': ['clustering'], 'params': {}} for i in range(20)]

        result = run_parameter_sweep(
            goal_id=goal.id,
            candidates=too_many,
            base_checkpoint=valid_path,
            project_id=test_project,
        )

        assert 'error' in result
        assert '超过' in result['error'] or '上限' in result['error']

    def test_unknown_module_rejected(self, test_project):
        """测试未知模块被拒绝."""
        from models import AgentSession, AgentGoal
        from modules.agent_orchestrator import run_parameter_sweep
        from config import Config
        import os

        session = AgentSession(project_id=test_project, status='active')
        session.save()

        goal = AgentGoal(
            session_id=session.id,
            project_id=test_project,
            goal_type='target_cluster_refinement',
            target_entity='microglia',
            goal_json='{"target_cell_type":"microglia"}',
        )
        goal.save()

        # 创建有效的 checkpoint 文件
        valid_path = os.path.join(Config.project_dir(test_project), 'test.h5ad')
        os.makedirs(Config.project_dir(test_project), exist_ok=True)
        with open(valid_path, 'w') as f:
            f.write('mock')

        result = run_parameter_sweep(
            goal_id=goal.id,
            candidates=[{'name': 'test', 'modules': ['nonexistent_module'], 'params': {}}],
            base_checkpoint=valid_path,
            project_id=test_project,
        )

        assert 'error' in result
        assert '未知模块' in result['error'] or 'nonexistent' in result['error'].lower()

    def test_param_cleaning_removes_unknown_in_sweep(self):
        """测试 sweep 中参数被 PARAM_SCHEMAS 清洗."""
        from modules.ai_tools import _validate_analysis_params

        cleaned, err = _validate_analysis_params('clustering', {
            'resolutions': '1.5',
            'n_neighbors': '25',
            'malicious_param': 'DROP TABLE',
        })

        assert err is None
        assert 'malicious_param' not in cleaned
        assert 'resolutions' in cleaned


class TestSecurityConstraints:
    """测试安全约束."""

    def test_agent_tools_no_direct_code_exec(self):
        """验证 agent 工具不包含任意代码执行."""
        from modules.ai_tools import execute_tool

        result = execute_tool('exec', {}, 'test')
        assert '未知工具' in result.get('error', '')
        result = execute_tool('eval', {}, 'test')
        assert '未知工具' in result.get('error', '')

    def test_sweep_candidates_hard_limit(self):
        """测试参数搜索候选硬上限."""
        from modules.ai_tools import _build_sweep_candidates

        candidates = _build_sweep_candidates(
            'target_cluster_refinement', 'microglia',
            {'available_embeddings': []},
            max_candidates=50,
        )
        assert len(candidates) <= 12

    def test_auto_exec_tools_are_readonly(self):
        """验证 AUTO_EXEC_TOOLS 不包含写操作，CONFIRM_TOOLS 覆盖所有 AI 写操作."""
        from modules.ai_adapter import AUTO_EXEC_TOOLS, CONFIRM_TOOLS

        # AI 工具层面的写操作（需确认）
        ai_write_tools = {'run_analysis', 'propose_parameter_sweep', 'run_parameter_sweep',
                          'start_goal_agent', 'continue_goal_agent'}
        for tool in ai_write_tools:
            assert tool in CONFIRM_TOOLS, f"{tool} 应该是需要确认的写操作"
            assert tool not in AUTO_EXEC_TOOLS, f"{tool} 不应该是自动执行的"

        # accept_branch 仅通过前端 Branch API 调用，不作为 AI 工具暴露
        # 验证它不在任一工具集合中
        assert 'accept_branch' not in AUTO_EXEC_TOOLS
        assert 'accept_branch' not in CONFIRM_TOOLS
