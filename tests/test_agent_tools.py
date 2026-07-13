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
            'recommend_analysis_config',
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


class TestBulkUploadDiscovery:
    """Bulk 原始表格应能被 AI 发现和检查。"""

    @staticmethod
    def _write_bulk_table(path):
        path.write_text(
            "Geneid\tGeneName\tStart\tCtrl_1\tTreat_1\n"
            "ENSG1\tGENE1\t1\t10\t20\n"
            "ENSG2\tGENE2\t2\t30\t40\n",
            encoding='utf-8',
        )

    def test_project_status_lists_xls_and_tsv(self, test_project):
        from config import Config
        from modules.ai_tools import _get_project_status

        uploads = Config.uploads_dir(test_project)
        self._write_bulk_table(__import__('pathlib').Path(uploads) / 'counts.xls')
        self._write_bulk_table(__import__('pathlib').Path(uploads) / 'counts.tsv')

        result = _get_project_status(test_project)
        assert 'counts.xls' in result['uploaded_files']
        assert 'counts.tsv' in result['uploaded_files']
        xls_detail = next(item for item in result['uploaded_file_details']
                          if item['name'] == 'counts.xls')
        assert xls_detail['path'].endswith('/uploads/counts.xls')
        assert xls_detail['extension'] == '.xls'

    def test_bulk_analysis_resolver_falls_back_to_xls(self, test_project):
        from config import Config
        from modules.ai_tools import resolve_current_analysis_input

        path = __import__('pathlib').Path(Config.uploads_dir(test_project)) / 'counts.xls'
        self._write_bulk_table(path)

        result = resolve_current_analysis_input(test_project, 'bulk_qc')
        assert result['source'] == 'bulk_upload'
        assert result['path'] == str(path)

    def test_inspect_adata_reads_tab_delimited_xls(self, test_project):
        from config import Config
        from modules.ai_tools import _inspect_adata

        path = __import__('pathlib').Path(Config.uploads_dir(test_project)) / 'counts.xls'
        self._write_bulk_table(path)

        result = _inspect_adata({}, test_project)
        assert 'error' not in result
        assert result['data_type'] == 'bulk_expression_matrix'
        assert result['input_file'] == 'counts.xls'
        assert result['n_obs'] == 2
        assert result['n_vars'] == 2
        assert result['sample_names'] == ['Ctrl_1', 'Treat_1']


class TestAnalysisConfigRecommendation:
    """对话 AI 应基于数据证据决策方法和参数。"""

    @staticmethod
    def _write_multifactor_table(path, continuous=True):
        samples = [
            'Ctr_B_1', 'Ctr_B_2', 'Ctr_En_1', 'Ctr_En_2',
            'PEA_B_1', 'PEA_B_2', 'PEA_En_1', 'PEA_En_2',
        ]
        rows = ['gene\t' + '\t'.join(samples)]
        values = [
            [1.1, 1.2, 2.1, 2.2, 3.1, 3.2, 4.1, 4.2],
            [2.2, 2.3, 1.2, 1.3, 5.2, 5.3, 3.2, 3.3],
            [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
        ]
        if not continuous:
            values = [[int(round(value * 10)) for value in row] for row in values]
        for idx, row in enumerate(values, 1):
            rows.append(f'G{idx}\t' + '\t'.join(map(str, row)))
        path.write_text('\n'.join(rows) + '\n', encoding='utf-8')

    def test_continuous_bulk_recommends_log2(self, test_project):
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _recommend_analysis_config

        path = Path(Config.uploads_dir(test_project)) / 'continuous.tsv'
        self._write_multifactor_table(path, continuous=True)
        result = _recommend_analysis_config({
            'module_name': 'bulk_normalize', 'input_path': str(path),
        }, test_project)

        assert result['data_profile']['measurement_type'] == 'continuous_expression'
        assert result['recommended_params']['method'] == 'log2'
        assert result['should_run'] is True

    def test_raw_counts_recommend_deseq2_normalization(self, test_project):
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _recommend_analysis_config

        path = Path(Config.uploads_dir(test_project)) / 'counts.tsv'
        self._write_multifactor_table(path, continuous=False)
        result = _recommend_analysis_config({
            'module_name': 'bulk_normalize', 'input_path': str(path),
        }, test_project)

        assert result['data_profile']['measurement_type'] == 'raw_counts'
        assert result['recommended_params']['method'] == 'deseq2'

    def test_bulk_deg_recommends_stratified_comparisons_and_prerequisite(self, test_project):
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _recommend_analysis_config

        path = Path(Config.uploads_dir(test_project)) / 'continuous.tsv'
        self._write_multifactor_table(path, continuous=True)
        result = _recommend_analysis_config({
            'module_name': 'bulk_deg',
            'input_path': str(path),
            'objective': '保留 B/En 分层差异',
        }, test_project)

        params = result['recommended_params']
        assert params['method'] == 't-test'
        assert params['groupby'] == '_auto_group_'
        assert params['comparisons'] == 'PEA_B-vs-Ctr_B;PEA_En-vs-Ctr_En'
        assert result['should_run'] is False
        assert result['prerequisites'][0]['params'] == {'method': 'log2'}

    def test_strict_objective_tightens_deg_thresholds(self, test_project):
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _recommend_analysis_config

        path = Path(Config.uploads_dir(test_project)) / 'counts.tsv'
        self._write_multifactor_table(path, continuous=False)
        result = _recommend_analysis_config({
            'module_name': 'bulk_deg', 'input_path': str(path), 'objective': '严格验证',
        }, test_project)

        assert result['recommended_params']['fc_threshold'] == 2
        assert result['recommended_params']['pval_threshold'] == 0.01

    def test_compatibility_gate_rejects_deseq2_for_continuous_values(self, test_project):
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _validate_method_compatibility

        path = Path(Config.uploads_dir(test_project)) / 'continuous.tsv'
        self._write_multifactor_table(path, continuous=True)

        error = _validate_method_compatibility(
            'bulk_deg', str(path), {'method': 'deseq2'})
        assert '原始整数 counts' in error

    def test_run_analysis_does_not_create_incompatible_task(self, test_project):
        from pathlib import Path
        from config import Config
        from models import AnalysisTask
        from modules.ai_tools import _run_analysis

        path = Path(Config.uploads_dir(test_project)) / 'continuous.tsv'
        self._write_multifactor_table(path, continuous=True)
        result = _run_analysis({
            'module_name': 'bulk_normalize',
            'input_path': str(path),
            'params': {'method': 'deseq2'},
        }, test_project)

        assert '方法与数据不兼容' in result['error']
        assert AnalysisTask.get_by_project(test_project) == []

    def test_auto_group_mapping_for_ai_uses_combined_groups(self, test_project):
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _infer_auto_group_mapping

        path = Path(Config.uploads_dir(test_project)) / 'continuous.tsv'
        self._write_multifactor_table(path, continuous=True)
        mapping = _infer_auto_group_mapping(str(path))

        assert mapping['Ctr_B_1'] == 'Ctr_B'
        assert mapping['Ctr_En_1'] == 'Ctr_En'
        assert mapping['PEA_B_1'] == 'PEA_B'

    def test_recommendation_tool_is_read_only_and_sweep_proposal_requires_confirmation(self):
        from modules.ai_adapter import AUTO_EXEC_TOOLS, CONFIRM_TOOLS, TOOLS_ANTHROPIC

        tool_names = {tool['name'] for tool in TOOLS_ANTHROPIC}
        assert 'recommend_analysis_config' in tool_names
        assert 'recommend_analysis_config' in AUTO_EXEC_TOOLS
        assert 'recommend_analysis_config' not in CONFIRM_TOOLS
        assert 'propose_parameter_sweep' not in AUTO_EXEC_TOOLS
        assert 'propose_parameter_sweep' in CONFIRM_TOOLS


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
