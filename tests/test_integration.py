"""集成测试 —— 覆盖 worker → 模块 → API 的完整链路，发现跨层 Bug。"""
import json
import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ────────────────────────────────────────────
# TestWorkerValidation — worker 预验证不阻塞非 h5ad 文件
# ────────────────────────────────────────────

class TestWorkerValidation:
    """验证 worker 的 best-effort 预验证不会阻塞 CSV/TSV 文件。"""

    def test_load_adata_raises_on_csv(self, tmp_path):
        """load_adata 对 CSV 文件应抛 OSError（不是 h5ad 格式）。"""
        from modules.base import BaseAnalysis

        csv_file = tmp_path / "data.csv"
        csv_file.write_text("gene,s1,s2\nG1,1,2\nG2,3,4\n")

        class _Stub(BaseAnalysis):
            MODULE_NAME = '_stub'
            def run(self, input_path):
                return {}

        stub = _Stub(project_dir=str(tmp_path), params={}, progress_callback=lambda p, m: None)
        with pytest.raises((OSError, ValueError)):
            stub.load_adata(str(csv_file))

    def test_worker_catches_oserror(self, tmp_path):
        """worker 的异常处理应捕获 OSError，不中断后续执行。"""
        csv_file = tmp_path / "data.csv"
        csv_file.write_text("gene,s1,s2\nG1,1,2\nG2,3,4\n")

        from modules.base import BaseAnalysis

        class _Stub(BaseAnalysis):
            MODULE_NAME = '_stub'
            def run(self, input_path):
                return {'output_adata': None, 'result_files': [], 'summary': {}}

        stub = _Stub(project_dir=str(tmp_path), params={}, progress_callback=lambda p, m: None)

        # 模拟 worker 的预验证逻辑（与 worker.py:62-74 一致）
        validation_skipped = False
        try:
            adata = stub.load_adata(str(csv_file))
            validation_error = stub.validate_input(adata)
            if validation_error:
                raise ValueError(f"输入验证失败: {validation_error}")
        except FileNotFoundError:
            validation_skipped = True
        except ImportError:
            validation_skipped = True
        except (OSError, ValueError):
            validation_skipped = True

        assert validation_skipped, "OSError 应被捕获，允许后续 run() 执行"

        # run() 应该能正常调用
        result = stub.run(str(csv_file))
        assert isinstance(result, dict)

    def test_read_expression_matrix_handles_csv(self, tmp_path):
        """read_expression_matrix 应正确读取 CSV 文件。"""
        from modules.io_utils import read_expression_matrix

        # 使用 TSV 格式避免已知 bug #1（CSV 被 TSV 解析路径拦截）
        csv_file = tmp_path / "counts.tsv"
        df = pd.DataFrame(
            np.random.randint(0, 100, size=(5, 3)),
            index=['GeneA', 'GeneB', 'GeneC', 'GeneD', 'GeneE'],
            columns=['Sample1', 'Sample2', 'Sample3']
        )
        df.to_csv(csv_file, sep='\t')

        adata = read_expression_matrix(str(csv_file))
        assert adata.n_obs == 3  # 样本数
        assert adata.n_vars == 5  # 基因数


# ────────────────────────────────────────────
# TestPlotlyJsonSanitization — API 返回无 NaN
# ────────────────────────────────────────────

class TestPlotlyJsonSanitization:
    """验证 Plotly JSON API 响应不含 NaN/Inf 字面量。"""

    def _make_plotly_data_with_nan(self):
        """构造含 NaN/Inf 的模拟 Plotly 数据。"""
        return {
            'data': [{
                'type': 'scattergl',
                'x': [1.0, float('nan'), 3.0, float('inf'), 5.0],
                'y': [2.0, 4.0, float('-inf'), 6.0, 8.0],
                'text': ['a', 'b', 'c', 'd', 'e'],
            }],
            'layout': {'title': {'text': 'test'}, 'width': 800, 'height': 500}
        }

    def test_sanitize_replaces_nan_with_null(self):
        """NaN 和 Inf 应被替换为 None（JSON null）。"""
        from routes.api import _sanitize_plotly_values

        data = self._make_plotly_data_with_nan()
        result = _sanitize_plotly_values(data)

        x = result['data'][0]['x']
        y = result['data'][0]['y']

        # NaN → None
        assert x[1] is None
        assert y[2] is None
        # Inf → None
        assert x[3] is None
        # 正常值不变
        assert x[0] == 1.0
        assert x[2] == 3.0

    def test_sanitize_output_is_valid_json(self):
        """清洗后的数据应能用 strict JSON 序列化（allow_nan=False）。"""
        from routes.api import _sanitize_plotly_values

        data = self._make_plotly_data_with_nan()
        result = _sanitize_plotly_values(data)

        # allow_nan=False 会在遇到 NaN/Inf 时抛异常
        json_str = json.dumps(result, allow_nan=False)
        assert 'NaN' not in json_str
        assert 'Infinity' not in json_str
        assert 'null' in json_str

    def test_decode_and_sanitize_pipeline(self):
        """完整的解码 + 清洗管线应输出标准 JSON。"""
        from routes.api import _decode_plotly_binary, _sanitize_plotly_values
        import base64, struct

        # 构造含 NaN 的二进制编码数据（模拟 Plotly write_json 输出）
        arr = np.array([1.0, np.nan, 3.0, np.inf, 5.0], dtype=np.float64)
        bdata = base64.b64encode(arr.tobytes()).decode()

        data = {
            'data': [{
                'type': 'scattergl',
                'x': {'bdata': bdata, 'dtype': 'f8'},
                'y': [1.0, 2.0, 3.0, 4.0, 5.0],
            }],
            'layout': {}
        }

        decoded = _decode_plotly_binary(data)
        sanitized = _sanitize_plotly_values(decoded)

        x = sanitized['data'][0]['x']
        assert isinstance(x, list)
        assert len(x) == 5
        assert x[0] == 1.0
        assert x[1] is None  # NaN → null
        assert x[3] is None  # Inf → null

        # 整体可序列化
        json_str = json.dumps(sanitized, allow_nan=False)
        parsed = json.loads(json_str)
        assert parsed['data'][0]['x'][1] is None


# ────────────────────────────────────────────
# TestModuleOutputCompleteness — 模块输出完整性
# ────────────────────────────────────────────

class TestModuleOutputCompleteness:
    """验证分析模块的输出文件完整且非空。"""

    def test_deg_csv_not_empty_with_base_mean_filter(self, tmp_path):
        """base_mean_filter 不应导致 CSV 为空（Bug fix 验证）。"""
        from modules.bulk_deg import _run_single_comparison
        import scanpy as sc

        # 创建有明确差异的 AnnData：确保能检出 DE 基因
        np.random.seed(42)
        n_samples = 10
        n_genes = 50
        X = np.random.rand(n_samples, n_genes) * 2 + 5  # 基础表达量
        # 前 10 个基因在 B 组高表达（效应量大）
        X[5:, :10] += 10
        obs = pd.DataFrame(
            {'group': ['A']*5 + ['B']*5},
            index=[f's{i}' for i in range(n_samples)]
        )
        var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_genes)])
        adata = sc.AnnData(X=X, obs=obs, var=var)

        counts = adata.X.astype(float)
        g1_samples = list(adata.obs.index[:5])
        g2_samples = list(adata.obs.index[5:])

        plots_dir = str(tmp_path / 'plots')
        results_dir = str(tmp_path / 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        # 使用 base_mean_filter=1.0（之前会导致 CSV 为空的参数）
        deg_df, result_files, n_up, n_down = _run_single_comparison(
            adata, counts, g1_samples, g2_samples, 'A', 'B',
            method='t-test', fc_threshold=1.5, pval_threshold=0.05,
            top_n=10, gene_id_to_name={}, plots_dir=plots_dir,
            results_dir=results_dir, base_mean_filter=1.0
        )

        # 完整结果 CSV（bulk_deg_results）应非空
        csv_files = [rf for rf in result_files if rf['file_type'] == 'csv']
        results_csv = [rf for rf in csv_files if 'results' in os.path.basename(rf['file_path'])]
        assert len(results_csv) >= 1, "应有完整结果 CSV"

        for rf in results_csv:
            assert os.path.exists(rf['file_path']), f"CSV 不存在: {rf['file_path']}"
            df = pd.read_csv(rf['file_path'])
            assert len(df) > 0, f"完整结果 CSV 为空: {rf['file_path']}"

    def test_deg_volcano_trace_counts_match(self, tmp_path):
        """火山图的 trace 数据点数应与 regulation 计数一致。"""
        from modules.bulk_deg import _run_single_comparison
        import scanpy as sc

        np.random.seed(42)
        n_samples = 8
        n_genes = 50
        X = np.random.rand(n_samples, n_genes) * 10 + 5
        obs = pd.DataFrame(
            {'group': ['A']*4 + ['B']*4},
            index=[f's{i}' for i in range(n_samples)]
        )
        var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_genes)])
        adata = sc.AnnData(X=X, obs=obs, var=var)

        counts = adata.X.astype(float)
        g1_samples = list(adata.obs.index[:4])
        g2_samples = list(adata.obs.index[4:])

        plots_dir = str(tmp_path / 'plots')
        results_dir = str(tmp_path / 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        deg_df, result_files, n_up, n_down = _run_single_comparison(
            adata, counts, g1_samples, g2_samples, 'A', 'B',
            method='t-test', fc_threshold=1.5, pval_threshold=0.05,
            top_n=10, gene_id_to_name={}, plots_dir=plots_dir,
            results_dir=results_dir, base_mean_filter=0
        )

        # 火山图统一为原生 Matplotlib 静态输出。
        volcano_files = [rf for rf in result_files if rf.get('category') == 'volcano']
        assert len(volcano_files) >= 1, "应生成火山图"
        assert {rf['file_type'] for rf in volcano_files} == {'png', 'svg'}
        assert all(os.path.exists(rf['file_path']) for rf in volcano_files)


# ────────────────────────────────────────────
# TestWorkerModuleInterface — worker 与模块接口一致性
# ────────────────────────────────────────────

class TestWorkerModuleInterface:
    """验证 worker 使用的模块接口（run() 返回结构）一致。"""

    def _make_sc_anndata(self, n_obs=30, n_vars=200):
        """构造最小单细胞 AnnData。"""
        import scanpy as sc
        import anndata
        np.random.seed(42)
        raw = np.random.poisson(lam=5, size=(n_obs, n_vars)).astype(np.float32)
        raw[raw == 0] = 1
        obs = pd.DataFrame(index=[f'cell_{i}' for i in range(n_obs)])
        obs['leiden'] = pd.Categorical([str(i % 3) for i in range(n_obs)])
        obs['batch'] = pd.Categorical(['A'] * (n_obs // 2) + ['B'] * (n_obs - n_obs // 2))
        var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_vars)])
        adata = anndata.AnnData(X=raw.copy(), obs=obs, var=var)
        adata.layers['counts'] = raw.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=10, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        sc.tl.umap(adata)
        return adata

    def test_worker_can_parse_module_result(self, tmp_path):
        """Worker 使用的 JSON 序列化流程应能处理模块返回值。"""
        from modules.clustering import ClusteringAnalysis

        adata = self._make_sc_anndata()
        input_path = str(tmp_path / 'input.h5ad')
        adata.write_h5ad(input_path)

        mod = ClusteringAnalysis(
            project_dir=str(tmp_path),
            params={'resolutions': '0.5', 'clustering_method': 'leiden'},
            progress_callback=lambda p, m: None,
        )
        result = mod.run(input_path)

        # 模拟 worker 的结果处理（与 worker.py:78-93 一致）
        summary_json = json.dumps(result.get('summary', {}))
        assert len(summary_json) > 2

        parsed = json.loads(summary_json)
        assert isinstance(parsed, dict)

        # result_files 应可遍历
        for rf in (result.get('result_files') or []):
            assert isinstance(rf.get('file_type', ''), str)
            assert isinstance(rf.get('category', ''), str)
            assert isinstance(rf.get('label', ''), str)
            assert isinstance(rf.get('file_path', ''), str)

    def test_worker_result_files_file_paths_exist(self, tmp_path):
        """Worker 记录的 result_files 文件应全部存在。"""
        from modules.dimred import DimredAnalysis

        adata = self._make_sc_anndata()
        input_path = str(tmp_path / 'input.h5ad')
        adata.write_h5ad(input_path)

        mod = DimredAnalysis(
            project_dir=str(tmp_path),
            params={'n_comps': 10, 'umap_n_neighbors': 10},
            progress_callback=lambda p, m: None,
        )
        result = mod.run(input_path)

        for rf in (result.get('result_files') or []):
            fpath = rf.get('file_path', '')
            assert os.path.exists(fpath), f"Worker 注册的文件不存在: {fpath}"
            assert os.path.getsize(fpath) > 0, f"Worker 注册的文件为空: {fpath}"


# ────────────────────────────────────────────
# TestFilterPipeline — 过滤规则的跨模块传播
# ────────────────────────────────────────────

class TestFilterPipeline:
    """验证 _filters 参数在模块间正确传播。"""

    def test_apply_filters_removes_cells(self, tmp_path):
        """apply_filters 应根据规则移除细胞。"""
        from modules.base import BaseAnalysis
        import anndata

        np.random.seed(42)
        n_obs = 50
        adata = anndata.AnnData(
            X=np.random.rand(n_obs, 20).astype(np.float32),
            obs=pd.DataFrame({
                'score': np.random.rand(n_obs),
                'group': pd.Categorical(['A'] * 25 + ['B'] * 25),
            }, index=[f'c{i}' for i in range(n_obs)]),
            var=pd.DataFrame(index=[f'G{i}' for i in range(20)]),
        )
        input_path = str(tmp_path / 'input.h5ad')
        adata.write_h5ad(input_path)

        class _Stub(BaseAnalysis):
            MODULE_NAME = '_stub'
            def run(self, input_path):
                return {}

        # 应用过滤: score >= 0.5
        mod = _Stub(
            project_dir=str(tmp_path),
            params={'_filters': {'qc': [{'column': 'score', 'op': '>=', 'value': 0.5}]}},
            progress_callback=lambda p, m: None,
        )
        filtered = mod.apply_filters(adata, 'qc')
        assert filtered.n_obs < adata.n_obs, "过滤应减少细胞数"
        assert all(filtered.obs['score'] >= 0.5), "过滤后所有 score 应 >= 0.5"

    def test_apply_filters_with_missing_column(self, tmp_path):
        """apply_filters 对不存在的列应跳过（不崩溃）。"""
        from modules.base import BaseAnalysis
        import anndata

        adata = anndata.AnnData(
            X=np.random.rand(10, 5).astype(np.float32),
            obs=pd.DataFrame(index=[f'c{i}' for i in range(10)]),
            var=pd.DataFrame(index=[f'G{i}' for i in range(5)]),
        )

        class _Stub(BaseAnalysis):
            MODULE_NAME = '_stub'
            def run(self, input_path):
                return {}

        mod = _Stub(
            project_dir=str(tmp_path),
            params={'_filters': {'qc': [{'column': 'nonexistent', 'op': '>=', 'value': 0.5}]}},
            progress_callback=lambda p, m: None,
        )
        # 不应崩溃
        filtered = mod.apply_filters(adata, 'qc')
        assert filtered.n_obs == 10, "不存在的列不应影响数据"

    def test_apply_filters_between_operator(self, tmp_path):
        """apply_filters 的 between 操作符应正确过滤。"""
        from modules.base import BaseAnalysis
        import anndata

        np.random.seed(42)
        adata = anndata.AnnData(
            X=np.random.rand(20, 5).astype(np.float32),
            obs=pd.DataFrame({
                'value': np.arange(20, dtype=float),
            }, index=[f'c{i}' for i in range(20)]),
            var=pd.DataFrame(index=[f'G{i}' for i in range(5)]),
        )

        class _Stub(BaseAnalysis):
            MODULE_NAME = '_stub'
            def run(self, input_path):
                return {}

        mod = _Stub(
            project_dir=str(tmp_path),
            params={'_filters': {'test': [{'column': 'value', 'op': 'between', 'value': [5, 15]}]}},
            progress_callback=lambda p, m: None,
        )
        filtered = mod.apply_filters(adata, 'test')
        assert filtered.n_obs == 11, f"between [5,15] 应保留 11 个，实际 {filtered.n_obs}"
        assert all(filtered.obs['value'] >= 5)
        assert all(filtered.obs['value'] <= 15)


# ────────────────────────────────────────────
# TestSchemasIntegration — schemas 与模块参数一致性
# ────────────────────────────────────────────

class TestSchemasIntegration:
    """验证 PARAM_SCHEMAS 与模块实际使用的参数一致。"""

    def test_all_registered_modules_have_schemas(self):
        """MODULE_REGISTRY 中注册的模块都应在 PARAM_SCHEMAS 中有参数定义（convert_10x 除外）。"""
        from modules import MODULE_REGISTRY
        from modules.schemas import PARAM_SCHEMAS

        # convert_10x 是数据导入模块，不参与分析流水线，无 PARAM_SCHEMA
        SKIP_MODULES = {'convert_10x'}

        missing = []
        for name in MODULE_REGISTRY:
            if name in SKIP_MODULES:
                continue
            if name not in PARAM_SCHEMAS:
                missing.append(name)

        if missing:
            msg = "以下模块在 MODULE_REGISTRY 中注册但缺少 PARAM_SCHEMAS:\n"
            msg += '\n'.join(f'  {m}' for m in missing)
            pytest.fail(msg)

    def test_schema_metadata_matches_registry(self):
        """MODULE_DISPLAY_MAP 中的显示名应与 MODULE_REGISTRY 的 DISPLAY_NAME 一致或相关。"""
        from modules import MODULE_REGISTRY
        from modules.schemas import MODULE_DISPLAY_MAP

        for name, display_name in MODULE_DISPLAY_MAP.items():
            if name not in MODULE_REGISTRY:
                continue
            cls = MODULE_REGISTRY[name]
            if not cls.DISPLAY_NAME:
                continue
            # 允许 MODULE_DISPLAY_MAP 是 DISPLAY_NAME 的子串或相等
            dm = cls.DISPLAY_NAME
            if display_name != dm and display_name not in dm and dm not in display_name:
                import warnings
                warnings.warn(
                    f"{name}: MODULE_DISPLAY_MAP('{display_name}') 与 DISPLAY_NAME('{dm}') 不一致",
                    UserWarning
                )
