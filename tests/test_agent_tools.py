# tests/test_agent_tools.py
"""Agent 工具函数测试 — 只读工具和路径安全"""
import pytest
import os
import json


class TestExecuteToolRegistry:
    """测试工具注册表."""

    def test_all_readonly_tools_registered(self):
        """验证所有只读工具在 execute_tool 中注册."""
        from modules.ai_tools import execute_tool

        readonly_tools = [
            'get_project_status', 'get_task_results', 'list_modules',
            'inspect_analysis_state', 'inspect_adata', 'get_cluster_summary',
            'score_cell_type_signature', 'list_builtin_markers',
        ]

        for tool in readonly_tools:
            # 测试调用一个不存在的任务也应该返回合理的错误（不会抛 KeyError）
            result = execute_tool(tool, {}, 'nonexistent_project')
            # 对于项目不存在的错误是可以接受的，但不应该是"未知工具"
            assert '未知工具' not in str(result.get('error', '')), f"{tool} 未注册"

    def test_confirm_tools_registered(self):
        """验证所有需要确认的工具在 execute_tool 中注册."""
        from modules.ai_tools import execute_tool

        confirm_tools = [
            'run_analysis', 'propose_parameter_sweep', 'run_parameter_sweep',
            'start_goal_agent', 'continue_goal_agent',
        ]

        for tool in confirm_tools:
            result = execute_tool(tool, {}, 'nonexistent_project')
            assert '未知工具' not in str(result.get('error', '')), f"{tool} 未注册"

    def test_unknown_tool_returns_error(self):
        """测试未知工具名返回 error."""
        from modules.ai_tools import execute_tool
        result = execute_tool('nonexistent_tool', {}, 'test')
        assert '未知工具' in result['error']


class TestPathValidation:
    """测试路径安全校验."""

    def test_validate_project_path_rejects_symlink(self, tmp_path, monkeypatch):
        """测试符号链接被拒绝."""
        from modules.ai_tools import _validate_project_path
        from config import Config

        monkeypatch.setattr(Config, 'DATA_DIR', str(tmp_path))

        real_file = tmp_path / 'real.h5ad'
        real_file.write_text('test')

        symlink_path = str(tmp_path / 'link.h5ad')
        os.symlink(str(real_file), symlink_path)

        err = _validate_project_path(symlink_path, 'test_pid')
        assert err is not None  # 应该返回错误

    def test_validate_project_path_rejects_outside_path(self, tmp_path, monkeypatch):
        """测试项目目录外的路径被拒绝."""
        from modules.ai_tools import _validate_project_path
        from config import Config

        monkeypatch.setattr(Config, 'DATA_DIR', str(tmp_path))
        outside_path = '/tmp/outside_file.h5ad'

        err = _validate_project_path(outside_path, 'test_pid')
        # 可能因为文件不存在而返回错误，但至少不应该接受
        assert err is not None

    def test_find_latest_adata_no_project(self):
        """测试项目不存在时返回 None."""
        from modules.ai_tools import _find_latest_adata
        result = _find_latest_adata('nonexistent_project')
        assert result is None


class TestAIToolsCoverage:
    """测试 AI tools 核心功能覆盖."""

    def test_list_modules_bulk(self):
        """测试列出 Bulk 模块."""
        from modules.ai_tools import _list_modules
        result = _list_modules('bulk')
        assert 'modules' in result
        for name in result['modules']:
            assert name.startswith('bulk_')

    def test_list_modules_sc(self):
        """测试列出单细胞模块."""
        from modules.ai_tools import _list_modules
        result = _list_modules('sc')
        assert 'modules' in result

    def test_list_builtin_markers(self):
        """测试列出内置 marker."""
        from modules.ai_tools import _list_builtin_markers
        result = _list_builtin_markers()
        assert 'cell_types' in result
        assert 'microglia' in result['cell_types']

    def test_validate_analysis_params_removes_unknown(self):
        """测试参数校验移除未知键."""
        from modules.ai_tools import _validate_analysis_params
        cleaned, err = _validate_analysis_params('clustering', {
            'resolutions': '1.0',
            'unknown_param': 'value',
        })
        assert err is None
        assert 'unknown_param' not in cleaned

    def test_validate_analysis_params_type_check(self):
        """测试参数类型校验."""
        from modules.ai_tools import _validate_analysis_params
        from modules.schemas import PARAM_SCHEMAS

        # 找到一个有 number 参数的模块
        for mod, schemas in PARAM_SCHEMAS.items():
            for s in schemas:
                if s.get('type') == 'number':
                    cleaned, err = _validate_analysis_params(mod, {s['key']: 'not_a_number'})
                    assert err is not None
                    return

        pytest.skip("没有找到带 number 类型参数的模块")


class TestProposeParameterSweep:
    """测试参数搜索候选生成."""

    def test_respects_max_candidates(self):
        """测试候选数不超过 max_candidates."""
        from modules.ai_tools import _build_sweep_candidates

        candidates = _build_sweep_candidates(
            'target_cluster_refinement', 'microglia',
            {'available_embeddings': ['X_umap']},
            max_candidates=3,
        )
        assert len(candidates) <= 3

    def test_candidate_structure(self):
        """测试候选结构完整."""
        from modules.ai_tools import _build_sweep_candidates

        candidates = _build_sweep_candidates(
            'target_cluster_refinement', 'microglia',
            {'available_embeddings': []},
            max_candidates=2,
        )

        for c in candidates:
            assert 'name' in c
            assert 'modules' in c
            assert 'params' in c
            assert 'rationale' in c
            assert 'clustering' in c['modules']
