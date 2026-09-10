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
            'run_analysis', 'run_pipeline', 'propose_parameter_sweep', 'run_parameter_sweep',
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
        assert 'functional_state' in result['modules']
        assert 'sc_pseudobulk_deg' in result['modules']
        assert 'neighborhood_da' in result['modules']

    def test_ibd_downstream_parameter_contract_matches_web_schema(self):
        """AI must expose the same form fields as the web entry point."""
        from modules.ai_tools import _get_module_parameters

        functional = _get_module_parameters({'module_name': 'functional_state'})
        functional_keys = {field['key'] for field in functional['parameters']}
        assert {'analysis_focus', 'sample_key', 'condition_key', 'celltype_key'} <= functional_keys

        neighbourhood = _get_module_parameters({'module_name': 'neighborhood_da'})
        neighbourhood_keys = {field['key'] for field in neighbourhood['parameters']}
        assert {'representation', 'n_neighbourhoods', 'sample_key', 'condition_key'} <= neighbourhood_keys

    def test_list_builtin_markers(self):
        """测试列出内置 marker."""
        from modules.ai_tools import _list_builtin_markers
        result = _list_builtin_markers()
        assert 'cell_types' in result
        assert 'microglia' in result['cell_types']
        assert 'organoid_panels' in result
        assert 'kidney' in result['organoid_panels']
        assert 'NPHS2' in result['organoid_panels']['kidney']['cell_types']['Podocytes']

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

    def test_run_pipeline_submits_single_sequential_job(self, test_project, monkeypatch):
        from pathlib import Path
        from config import Config
        from models import PipelineRun
        from modules.ai_tools import execute_tool

        input_path = Path(Config.uploads_dir(test_project)) / 'input.h5ad'
        input_path.write_bytes(b'placeholder')
        submitted = {}

        def fake_submit_pipeline_run(**kwargs):
            submitted.update(kwargs)
            return True

        monkeypatch.setattr('worker.submit_pipeline_run', fake_submit_pipeline_run)
        result = execute_tool('run_pipeline', {
            'analysis_type': 'sc',
            'modules': ['qc', 'normalize', 'hvg'],
            'input_path': str(input_path),
            'params': {'hvg': {'n_top_genes': '3000'}},
        }, test_project)

        assert result['status'] == 'submitted'
        assert result['pipeline_url'].endswith('/pipeline-runs/' + result['pipeline_run_id'])
        assert submitted['modules'] == ['qc', 'normalize', 'hvg']
        assert submitted['input_path'] == str(input_path)
        run = PipelineRun.get_by_id(result['pipeline_run_id'])
        assert run is not None
        assert json.loads(run.modules_json) == ['qc', 'normalize', 'hvg']


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

    @staticmethod
    def _write_ibd_organoid_h5ad(path):
        """Write a minimal, replicated IBD/control epithelial input for AI recommendations."""
        anndata = pytest.importorskip('anndata')
        import numpy as np
        import pandas as pd

        rows = []
        for condition, sample_id in (
            ('Control', 'CTRL_1'), ('Control', 'CTRL_2'),
            ('IBD', 'IBD_1'), ('IBD', 'IBD_2'),
        ):
            for celltype in ('TA cells', 'Metabolic enterocytes'):
                for cell_index in range(3):
                    rows.append({
                        'sample_id': sample_id,
                        'condition': condition,
                        'celltype': celltype,
                        'leiden': '0' if celltype == 'TA cells' else '1',
                        'barcode_note': f'{sample_id}_{cell_index}',
                    })
        counts = np.tile(np.array([4, 8, 12, 16], dtype=float), (len(rows), 1))
        adata = anndata.AnnData(
            counts,
            obs=pd.DataFrame(rows),
            var=pd.DataFrame(index=['HNF4A', 'PPARA', 'MKI67', 'NFKB1']),
        )
        adata.layers['counts'] = counts.copy()
        adata.obsm['X_pca'] = np.arange(len(rows) * 3, dtype=float).reshape(len(rows), 3)
        adata.obsm['X_umap'] = np.arange(len(rows) * 2, dtype=float).reshape(len(rows), 2)
        # The recommender only needs the persisted graph marker, not graph data.
        adata.uns['neighbors'] = {'params': {'n_neighbors': 5}}
        adata.write_h5ad(path)

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
        assert result['recommended_params']['input_measurement'] == 'auto'
        assert result['should_run'] is False
        assert any('不能证明' in warning for warning in result['warnings'])

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

    def test_run_analysis_blocks_singleton_bulk_group_before_task_creation(self, test_project):
        import anndata as ad
        import numpy as np
        import pandas as pd
        from config import Config
        from models import AnalysisTask
        from modules.ai_tools import _run_analysis

        path = os.path.join(Config.uploads_dir(test_project), 'singleton.h5ad')
        ad.AnnData(
            X=np.asarray([[1, 2], [2, 3], [5, 6]], dtype=float),
            obs=pd.DataFrame({'condition': ['Ctrl', 'Ctrl', 'Treat']}, index=['s1', 's2', 's3']),
            var=pd.DataFrame(index=['G1', 'G2']),
        ).write_h5ad(path)
        result = _run_analysis({
            'module_name': 'bulk_deg',
            'input_path': path,
            'params': {'groupby': 'condition', 'method': 't-test'},
        }, test_project)

        assert '分析前检查未通过' in result['error']
        assert AnalysisTask.get_by_project(test_project) == []

    def test_ai_sc_go_binds_completed_deg_task_to_fixed_internal_sources(self, test_project):
        from pathlib import Path
        from models import AnalysisTask
        from config import Config
        from modules.ai_tools import _bind_sc_cell_go_source

        source_path = Path(Config.results_dir(test_project)) / 'deg_source.csv'
        source_path.parent.mkdir(parents=True, exist_ok=True)
        source_path.write_text('comparison_id,cluster\ndemo,9\n', encoding='utf-8')
        source_task = AnalysisTask(
            project_id=test_project, module_name='sc_pseudobulk_deg', status='completed',
            result_json=json.dumps({
                'deg_source_files': [str(source_path)],
                'deg_source_level': 'pseudobulk',
            }),
        )
        source_task.save()

        params = {'deg_source_task_id': source_task.id, 'source_level': 'cell_level'}
        assert _bind_sc_cell_go_source(params, test_project) is None
        assert params['deg_source_files'] == [str(source_path)]
        assert params['source_level'] == 'pseudobulk'
        assert params['source_task_id'] == source_task.id

    def test_ai_can_query_every_web_adjustable_parameter(self):
        from modules import MODULE_REGISTRY
        from modules.ai_tools import _get_module_parameters
        from modules.schemas import PARAM_SCHEMAS

        for module_name in MODULE_REGISTRY:
            result = _get_module_parameters({'module_name': module_name})
            assert 'error' not in result
            returned_keys = {field['key'] for field in result['parameters']}
            expected_keys = {field['key'] for field in PARAM_SCHEMAS.get(module_name, [])}
            assert returned_keys == expected_keys

    def test_ai_preserves_active_custom_go_display_parameters(self):
        from modules.ai_tools import _validate_analysis_params

        requested = {
            'go_priority_allocation_mode': 'custom',
            'go_inflammation_slots': 3,
            'go_lipid_slots': 2,
            'go_confirmed_theme_slots': 1,
            'go_top_pathway_slots': 2,
            'plot_top_n': 8,
        }
        cleaned, error = _validate_analysis_params('sc_cell_go', requested)

        assert error is None
        assert cleaned == requested

    def test_ai_activates_an_unambiguous_dependent_web_parameter(self):
        from modules.ai_tools import _validate_analysis_params

        cleaned, error = _validate_analysis_params(
            'sc_cell_go', {'go_inflammation_slots': 3},
        )

        assert error is None
        assert cleaned == {
            'go_priority_allocation_mode': 'custom',
            'go_inflammation_slots': 3,
        }

    def test_ai_rejects_a_number_outside_the_web_parameter_range(self):
        from modules.ai_tools import _validate_analysis_params

        cleaned, error = _validate_analysis_params(
            'sc_cell_go', {'plot_top_n': 31},
        )

        assert cleaned is None
        assert '不能大于 30' in error

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

    def test_sc_timecourse_recommendation_uses_safe_metadata_candidates(self, test_project):
        """AI profiling must exercise the NumPy/Pandas path and avoid batch IDs."""
        anndata = pytest.importorskip('anndata')
        import numpy as np
        import pandas as pd
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _recommend_analysis_config

        rows = []
        for timepoint in ('D0', 'D3', 'D7'):
            for replicate in (1, 2):
                for celltype in ('A', 'B'):
                    rows.append({
                        'timepoint': timepoint,
                        'sample_id': f'{timepoint}_R{replicate}',
                        'celltype': celltype,
                        'batch': f'technical_{replicate}',
                    })
        counts = np.arange(1, len(rows) * 4 + 1, dtype=float).reshape(len(rows), 4)
        adata = anndata.AnnData(
            counts,
            obs=pd.DataFrame(rows),
            var=pd.DataFrame(index=['G1', 'G2', 'G3', 'G4']),
        )
        adata.layers['counts'] = counts.copy()
        path = Path(Config.uploads_dir(test_project)) / 'temporal.h5ad'
        adata.write_h5ad(path)

        result = _recommend_analysis_config({
            'module_name': 'sc_timecourse', 'input_path': str(path),
        }, test_project)

        assert 'error' not in result
        assert result['should_run'] is True
        assert result['recommended_params']['timepoint_key'] == 'timepoint'
        assert result['recommended_params']['sample_key'] == 'sample_id'
        assert result['data_profile']['counts_layer_is_raw'] is True
        assert 'batch' not in result['data_profile']['sample_candidates']
        assert 'batch' not in result['data_profile']['time_candidates']

    def test_ibd_downstream_recommendations_choose_sample_level_modules(self, test_project):
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _recommend_analysis_config

        path = Path(Config.uploads_dir(test_project)) / 'ibd_organoid.h5ad'
        self._write_ibd_organoid_h5ad(path)

        functional = _recommend_analysis_config({
            'module_name': 'functional_state', 'input_path': str(path),
            'objective': 'IBD 类器官中炎症应激、TA 增殖和代谢 enterocyte 分化状态',
        }, test_project)
        assert functional['should_run'] is True
        assert functional['recommended_params']['analysis_focus'] == 'ibd_organoid_epithelial'
        assert functional['recommended_params']['sample_key'] == 'sample_id'
        assert functional['recommended_params']['condition_key'] == 'condition'

        pseudobulk = _recommend_analysis_config({
            'module_name': 'sc_pseudobulk_deg', 'input_path': str(path),
            'objective': '比较 IBD 与对照中各 cell type 的表达变化',
        }, test_project)
        assert pseudobulk['should_run'] is True
        assert pseudobulk['recommended_params']['grouping_mode'] == 'annotated_celltype'
        assert pseudobulk['recommended_params']['celltype_key'] == 'celltype'
        assert pseudobulk['recommended_params']['comparisons'] == 'IBD-vs-Control'

        neighbourhood = _recommend_analysis_config({
            'module_name': 'neighborhood_da', 'input_path': str(path),
            'objective': '查找 IBD 富集的连续 epithelial state',
        }, test_project)
        assert neighbourhood['should_run'] is True
        assert neighbourhood['recommended_params']['representation'] == 'X_pca'
        assert neighbourhood['recommended_params']['comparisons'] == 'IBD-vs-Control'

        proportion = _recommend_analysis_config({
            'module_name': 'proportion', 'input_path': str(path),
            'objective': '检验 TA 和 metabolic enterocyte 的比例变化',
        }, test_project)
        assert proportion['should_run'] is True
        assert proportion['recommended_params']['analysis_unit'] == 'sample'
        assert proportion['recommended_params']['groupby'] == 'celltype'

    def test_organoid_annotation_recommends_matching_marker_panel(self, test_project):
        """自然语言中的类器官类型应映射到对应的 annotation 参数。"""
        anndata = pytest.importorskip('anndata')
        import numpy as np
        import pandas as pd
        from pathlib import Path
        from config import Config
        from modules.ai_tools import _recommend_analysis_config

        adata = anndata.AnnData(
            np.ones((6, 4), dtype=float),
            obs=pd.DataFrame({'leiden': ['0', '0', '1', '1', '1', '0']}),
            var=pd.DataFrame(index=['NPHS1', 'NPHS2', 'PODXL', 'WT1']),
        )
        path = Path(Config.uploads_dir(test_project)) / 'kidney_organoid.h5ad'
        adata.write_h5ad(path)

        result = _recommend_analysis_config({
            'module_name': 'annotation', 'input_path': str(path),
            'objective': '肾类器官细胞注释',
        }, test_project)

        assert result['recommended_params']['marker_set'] == 'Organoid'
        assert result['recommended_params']['organoid_type'] == 'kidney'
        assert result['recommended_params']['method'] == 'multi_evidence'

    def test_recommendation_tool_is_read_only_and_sweep_proposal_requires_confirmation(self):
        from modules.ai_adapter import (
            AUTO_EXEC_TOOLS, CONFIRM_TOOLS, TOOLS_ANTHROPIC, SYSTEM_PROMPT,
        )

        tool_names = {tool['name'] for tool in TOOLS_ANTHROPIC}
        assert 'recommend_analysis_config' in tool_names
        assert 'get_module_parameters' in tool_names
        assert 'recommend_analysis_config' in AUTO_EXEC_TOOLS
        assert 'get_module_parameters' in AUTO_EXEC_TOOLS
        assert 'recommend_analysis_config' not in CONFIRM_TOOLS
        assert 'propose_parameter_sweep' not in AUTO_EXEC_TOOLS
        assert 'propose_parameter_sweep' in CONFIRM_TOOLS
        assert 'IBD / 对照类器官单细胞下游流程' in SYSTEM_PROMPT
        assert 'neighborhood_da' in SYSTEM_PROMPT


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
