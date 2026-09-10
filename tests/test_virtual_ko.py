# tests/test_virtual_ko.py
"""虚拟敲除（CellOracle virtual KO）模块的注册、参数、路径校验与子进程协议测试。

本文件只测试平台侧（Python 3.12）逻辑，不要求 celloracle 环境存在；
真实 CellOracle 计算由 modules/celloracle_worker.py 在独立环境中完成；
端到端验证方式：准备 job JSON 后用 celloracle 环境执行
「<celloracle_env_python> modules/celloracle_worker.py --config <job.json>」。
"""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def _make_adata(n_obs=60, n_vars=200):
    import numpy as np
    import anndata as ad
    rng = np.random.default_rng(0)
    counts = rng.poisson(2.0, size=(n_obs, n_vars)).astype(np.float32)
    adata = ad.AnnData(X=counts)
    adata.layers['counts'] = counts.copy()
    adata.var_names = ['GENE%d' % i for i in range(n_vars)]
    adata.obs_names = ['cell%d' % i for i in range(n_obs)]
    adata.obs['celltype'] = ['typeA'] * 20 + ['typeB'] * 20 + ['typeC'] * 20
    adata.obs['leiden'] = ['0'] * 30 + ['1'] * 30
    adata.obsm['X_umap'] = rng.normal(size=(n_obs, 2)).astype(np.float32)
    adata.var_names = adata.var_names.astype(str)
    adata.var_names = adata.var_names[:2].tolist() + adata.var_names[2:].tolist()
    adata.var['gene_symbols'] = adata.var_names.values
    return adata


class TestVirtualKORegistration:
    def test_module_registered(self):
        from modules import MODULE_REGISTRY
        assert 'virtual_ko' in MODULE_REGISTRY
        cls = MODULE_REGISTRY['virtual_ko']
        assert cls.MODULE_NAME == 'virtual_ko'
        assert '虚拟' in cls.DISPLAY_NAME
        assert hasattr(cls, 'run')

    def test_pipeline_order_and_deps(self):
        from modules import PIPELINE_ORDER, PIPELINE_DEPS, validate_pipeline_order
        assert 'virtual_ko' in PIPELINE_ORDER
        assert PIPELINE_ORDER.index('clustering') < PIPELINE_ORDER.index('virtual_ko')
        assert 'clustering' in PIPELINE_DEPS['virtual_ko']
        valid, errors = validate_pipeline_order(PIPELINE_ORDER)
        assert valid, errors

    def test_schema_present(self):
        from modules.schemas import SC_MODULE_LIST, PARAM_SCHEMAS, SC_MODULE_NAMES
        assert any(m['name'] == 'virtual_ko' for m in SC_MODULE_LIST)
        assert 'virtual_ko' in SC_MODULE_NAMES
        keys = {p['key'] for p in PARAM_SCHEMAS['virtual_ko']}
        for required in ('base_grn_source', 'perturb_genes', 'alpha',
                         'n_propagation', 'cluster_key'):
            assert required in keys

    def test_worker_script_syntax(self):
        import py_compile
        worker = os.path.join(os.path.dirname(__file__), '..',
                              'modules', 'celloracle_worker.py')
        py_compile.compile(worker, doraise=True)


class TestVirtualKOParams:
    def test_parse_form_params_builtin_drops_file(self):
        from modules.schemas import PARAM_SCHEMAS, parse_form_params
        from werkzeug.datastructures import MultiDict
        schema = PARAM_SCHEMAS['virtual_ko']
        form = MultiDict([
            ('base_grn_source', 'builtin_human_hg38'),
            ('base_grn_file', '/some/file.parquet'),
            ('perturb_genes', 'SPI1'),
        ])
        params = parse_form_params(schema, form)
        assert 'base_grn_file' not in params
        assert params['base_grn_source'] == 'builtin_human_hg38'

    def test_parse_form_params_upload_keeps_file(self):
        from modules.schemas import PARAM_SCHEMAS, parse_form_params
        from werkzeug.datastructures import MultiDict
        schema = PARAM_SCHEMAS['virtual_ko']
        form = MultiDict([
            ('base_grn_source', 'upload'),
            ('base_grn_file', '/tmp/x.parquet'),
            ('perturb_genes', 'SPI1'),
        ])
        params = parse_form_params(schema, form)
        assert params['base_grn_file'] == '/tmp/x.parquet'


class TestVirtualKOPaths:
    def test_resolve_base_grn_accepts_upload_file(self, tmp_path):
        from modules.virtual_ko import VirtualKOAnalysis
        uploads = tmp_path / 'uploads'
        uploads.mkdir()
        grn = uploads / 'human_base_grn.parquet'
        grn.write_bytes(b'x')
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={
            'base_grn_source': 'upload', 'base_grn_file': 'human_base_grn.parquet',
        }, progress_callback=None)
        source, path = mod._resolve_base_grn()
        assert source == 'upload'
        assert path == str(grn.resolve()) or os.path.samefile(path, grn)

    def test_resolve_base_grn_rejects_outside_path(self, tmp_path):
        from modules.virtual_ko import VirtualKOAnalysis
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={
            'base_grn_source': 'upload', 'base_grn_file': '/etc/passwd',
        }, progress_callback=None)
        with pytest.raises(ValueError, match='当前项目的 uploads'):
            mod._resolve_base_grn()

    def test_resolve_base_grn_rejects_pickle_even_in_uploads(self, tmp_path):
        from modules.virtual_ko import VirtualKOAnalysis
        uploads = tmp_path / 'uploads'
        uploads.mkdir()
        unsafe = uploads / 'untrusted.pkl'
        unsafe.write_bytes(b'not a real pickle')
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={
            'base_grn_source': 'upload', 'base_grn_file': str(unsafe),
        }, progress_callback=None)
        with pytest.raises(ValueError, match='不接受 pickle'):
            mod._resolve_base_grn()

    def test_links_file_is_disabled_before_subprocess(self, tmp_path):
        from modules.virtual_ko import VirtualKOAnalysis
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={
            'links_file': 'unsafe.links',
        }, progress_callback=None)
        with pytest.raises(ValueError, match='不接受用户提供的 CellOracle Links'):
            mod._resolve_links_file()

    def test_resolve_base_grn_missing_file(self, tmp_path):
        from modules.virtual_ko import VirtualKOAnalysis
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={
            'base_grn_source': 'upload', 'base_grn_file': 'nope.parquet',
        }, progress_callback=None)
        with pytest.raises(ValueError, match='不存在'):
            mod._resolve_base_grn()

    def test_build_job_config_builtin(self, tmp_path):
        from modules.virtual_ko import VirtualKOAnalysis
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={
            'base_grn_source': 'builtin_human_hg38',
            'perturb_genes': 'SPI1, GATA1',
            'cluster_key': 'celltype',
        }, progress_callback=None)
        config = mod._build_job_config(str(tmp_path / 'in.h5ad'))
        assert config['base_grn_source'] == 'builtin_human_hg38'
        assert config['perturb_genes'] == 'SPI1, GATA1'
        assert config['alpha'] == 10
        assert config['output_dir'].endswith('virtual_ko')
        assert os.path.isdir(config['output_dir'])


class TestVirtualKOSubprocess:
    def test_subprocess_protocol_and_files(self, tmp_path):
        """用伪 worker 验证 stdout 事件协议与结果文件登记。"""
        from modules.virtual_ko import VirtualKOAnalysis
        fake_py = tmp_path / 'fake_python'
        fake_worker = tmp_path / 'worker.py'
        out_dir = tmp_path / 'results' / 'virtual_ko'
        plots_dir = tmp_path / 'plots'
        out_dir.mkdir(parents=True)
        plots_dir.mkdir()
        fake_worker.write_text(
            'import json, os, sys\n'
            'cfg = json.load(open(sys.argv[2], encoding="utf-8"))\n'
            'print(json.dumps({"type": "progress", "pct": 50, "message": "half"}))\n'
            'with open(os.path.join(cfg["output_dir"], "a.csv"), "w") as f:\n'
            '    f.write("a,b\\n1,2\\n")\n'
            'summary = {"status": "completed", "output_h5ad": "OUT.h5ad", "genes": {}, "environment": dict(os.environ)}\n'
            'print(json.dumps({"type": "result", "summary": summary, "files": ['
            '{"path": os.path.join(cfg["output_dir"], "a.csv"), "kind": "csv", '
            '"category": "table", "label": "t"}]}))\n',
            encoding='utf-8')
        fake_py.write_text('#!/bin/sh\nexec python3 "$@"\n', encoding='utf-8')
        fake_py.chmod(0o755)

        progress_events = []
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={},
                                progress_callback=lambda pct, msg: progress_events.append((pct, msg)))
        # 通过 monkeypatch 直接替换解析函数，避免依赖 config 默认路径
        mod._celloracle_python = lambda: str(fake_py)
        mod._celloracle_worker = lambda: str(fake_worker)
        config = {
            'input_h5ad': str(tmp_path / 'in.h5ad'),
            'output_dir': str(out_dir),
            'plots_dir': str(plots_dir),
            'intermediate_dir': str(tmp_path / 'intermediate'),
            'base_grn_source': 'builtin_human_hg38',
            'base_grn_file': '',
        }
        files, summary = mod._run_subprocess(config)
        assert summary['status'] == 'completed'
        assert any('50' == str(pct) for pct, _ in progress_events)
        assert any(f['path'].endswith('a.csv') for f in files)
        assert os.path.isfile(out_dir / 'a.csv')

    def test_subprocess_does_not_inherit_secrets(self, tmp_path, monkeypatch):
        """The isolated CellOracle process only receives an explicit env allowlist."""
        from modules.virtual_ko import VirtualKOAnalysis
        fake_py = tmp_path / 'fake_python'
        fake_worker = tmp_path / 'worker.py'
        fake_py.write_text('#!/bin/sh\nexec python3 "$@"\n', encoding='utf-8')
        fake_py.chmod(0o755)
        fake_worker.write_text(
            'import json, os, sys\n'
            'json.load(open(sys.argv[2], encoding="utf-8"))\n'
            'print(json.dumps({"type": "result", "summary": {"environment": dict(os.environ)}, "files": []}))\n',
            encoding='utf-8',
        )
        monkeypatch.setenv('AI_API_KEY', 'must-not-reach-subprocess')
        monkeypatch.setenv('PLATFORM_ACCESS_PASSWORD', 'must-not-reach-subprocess')
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={}, progress_callback=None)
        mod._celloracle_python = lambda: str(fake_py)
        mod._celloracle_worker = lambda: str(fake_worker)
        _files, summary = mod._run_subprocess({
            'input_h5ad': str(tmp_path / 'in.h5ad'),
            'output_dir': str(tmp_path / 'results'),
            'plots_dir': str(tmp_path / 'plots'),
            'intermediate_dir': str(tmp_path / 'intermediate'),
            'base_grn_source': 'builtin_human_hg38', 'base_grn_file': '',
        })
        assert 'AI_API_KEY' not in summary['environment']
        assert 'PLATFORM_ACCESS_PASSWORD' not in summary['environment']

    def test_run_unavailable_when_python_missing(self, tmp_path):
        from modules.virtual_ko import VirtualKOAnalysis
        adata = _make_adata()
        in_h5ad = tmp_path / 'in.h5ad'
        adata.write_h5ad(in_h5ad)
        mod = VirtualKOAnalysis(project_dir=str(tmp_path), params={
            'base_grn_source': 'builtin_human_hg38',
            'perturb_genes': 'GENE1',
        }, progress_callback=None)
        mod._celloracle_python = lambda: str(tmp_path / 'missing_python')
        mod._celloracle_worker = lambda: str(tmp_path / 'missing_worker.py')
        result = mod.run(str(in_h5ad))
        assert result['summary']['status'] == 'unavailable'
        assert 'celloracle' in result['summary']['missing_dependencies']

    def test_validate_input_fallback_columns(self):
        from modules.virtual_ko import VirtualKOAnalysis
        adata = _make_adata()
        mod = VirtualKOAnalysis(project_dir='.', params={'cluster_key': ''},
                                progress_callback=None)
        assert mod.validate_input(adata) is None
        adata2 = adata.copy()
        for col in ('celltype', 'final_annotation', 'leiden'):
            if col in adata2.obs.columns:
                del adata2.obs[col]
        assert mod.validate_input(adata2) is not None


class TestBaseGRNListing:
    def test_list_base_grn_files(self, tmp_path, monkeypatch):
        from config import Config
        from routes.analysis import _list_base_grn_files
        pid = 'proj123'
        uploads = tmp_path / 'projects' / pid / 'uploads'
        uploads.mkdir(parents=True)
        (uploads / 'base.parquet').write_bytes(b'1')
        (uploads / 'tf_info.csv').write_bytes(b'2')
        (uploads / 'other.h5ad').write_bytes(b'3')
        monkeypatch.setattr(Config, 'DATA_DIR', str(tmp_path))
        files = _list_base_grn_files(pid)
        names = {f['name'] for f in files}
        assert names == {'base.parquet', 'tf_info.csv'}
