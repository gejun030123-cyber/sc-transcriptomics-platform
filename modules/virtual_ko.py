import json
import logging
import os
import signal
import subprocess
import uuid

from modules.base import BaseAnalysis

logger = logging.getLogger(__name__)

# CellOracle 可用的 base GRN 来源（参数面板与 worker 共用）
BUILTIN_BASE_GRN_OPTIONS = {
    'builtin_human_hg38': '内置人类 promoter base GRN（hg38, gimmemotifs v5 fpr2）',
    'builtin_human_hg19': '内置人类 promoter base GRN（hg19, gimmemotifs v5 fpr2）',
}

CLUSTER_KEY_CANDIDATES = (
    'celltype', 'final_annotation', 'cell_type_l2', 'cell_type_l1',
    'leiden', 'louvain_annot',
)

# ``pickle`` is code, not data.  A project upload is user-controlled and must
# never be passed to CellOracle's pickle loader.  Tabular TF-info matrices are
# sufficient for user-provided base GRNs; the maintained CellOracle resources
# remain available through the builtin options above.
SAFE_BASE_GRN_EXTENSIONS = ('.parquet', '.pq', '.csv', '.tsv', '.txt', '.gz')


class VirtualKOAnalysis(BaseAnalysis):
    """CellOracle 虚拟敲除（in silico perturbation）分析。

    平台主环境为 Python 3.12，无法 import celloracle，因此本模块只负责
    校验输入、生成任务配置，并以 subprocess 调用独立 celloracle 环境
    （Python 3.9/3.10）中的 modules/celloracle_worker.py 完成
    GRN 推断（fit_GRN_for_simulation）与敲除模拟（simulate_shift /
    transition probability / p_mass），结果图表由 worker 落盘后登记
    到任务结果清单。
    """

    MODULE_NAME = "virtual_ko"
    DISPLAY_NAME = "虚拟敲除"
    DESCRIPTION = "基于 CellOracle 的 GRN 推断与 in silico 基因敲除扰动模拟"
    INPUT_REQUIRES = ['X_umap']

    # ---------- 工具方法 ----------

    def _celloracle_python(self):
        from config import Config
        return os.path.abspath(getattr(Config, 'CELLORACLE_PYTHON', ''))

    def _celloracle_worker(self):
        from config import Config
        return os.path.abspath(getattr(Config, 'CELLORACLE_WORKER_SCRIPT', ''))

    def _resolve_uploaded_base_grn(self, raw):
        """Resolve one non-symlink tabular GRN from this project's uploads."""
        upload_root = os.path.realpath(os.path.join(self.project_dir, 'uploads'))
        requested = str(raw or '').strip()
        if not requested:
            raise ValueError('请选择上传的 base GRN 文件。')
        candidate = requested if os.path.isabs(requested) else os.path.join(upload_root, requested)
        if os.path.islink(candidate):
            raise ValueError('base GRN 文件不能是符号链接。')
        path = os.path.realpath(candidate)
        try:
            contained = os.path.commonpath([upload_root, path]) == upload_root
        except ValueError:
            contained = False
        if not contained:
            raise ValueError('base GRN 文件必须位于当前项目的 uploads 目录。')
        if not os.path.isfile(path):
            raise ValueError('base GRN 文件不存在: ' + path)
        if not path.lower().endswith(SAFE_BASE_GRN_EXTENSIONS):
            raise ValueError(
                'base GRN 只支持非可执行的表格格式：.parquet、.csv、.tsv 或 .txt；'
                '不接受 pickle/Oracle/Links 文件。'
            )
        return path

    def _resolve_base_grn(self):
        """Return a builtin source or a safe tabular file from project uploads."""
        source = str(self.params.get('base_grn_source', 'upload') or 'upload')
        if source.startswith('builtin_'):
            return source, ''
        raw = str(self.params.get('base_grn_file', '') or '').strip()
        return 'upload', self._resolve_uploaded_base_grn(raw)

    def _resolve_links_file(self):
        raw = str(self.params.get('links_file', '') or '').strip()
        if not raw:
            return ''
        raise ValueError(
            '当前版本不接受用户提供的 CellOracle Links/pickle 文件，'
            '因为反序列化不可信 pickle 会执行任意代码。请使用内置或表格 base GRN。'
        )

    def _parse_genes(self):
        import re
        text = str(self.params.get('perturb_genes', '') or '')
        # 注意：字符类中的 \s 是空白符；过去误写为字母 s，会把 Hes1/Csf1r
        # 等含小写 s 的基因名拆坏。
        parts = re.split(r'[,\s;，；]+', text)
        genes = [p.strip() for p in parts if p.strip()]
        return list(dict.fromkeys(genes))

    def _build_job_config(self, input_path):
        results_dir = os.path.join(self.project_dir, 'results', 'virtual_ko')
        plots_dir = os.path.join(self.project_dir, 'plots')
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(intermediate_dir, exist_ok=True)

        source, base_grn_path = self._resolve_base_grn()
        links_path = self._resolve_links_file()
        cluster_key = str(self.params.get('cluster_key', '') or '').strip()
        n_propagation = int(self.params.get('n_propagation', 3) or 3)
        # CellOracle 只接受 1–5；在启动子进程（GRN 推断可能数十分钟）之前
        # 就拒绝非法值，而不是等模拟阶段才抛底层错误。
        if not 1 <= n_propagation <= 5:
            raise ValueError('n_propagation 必须在 1–5 之间（CellOracle 限制）。')
        config = {
            'input_h5ad': os.path.abspath(input_path),
            'output_dir': results_dir,
            'plots_dir': plots_dir,
            'intermediate_dir': intermediate_dir,
            'base_grn_source': source,
            'base_grn_file': base_grn_path,
            'base_grn_version': str(self.params.get('base_grn_version', '') or ''),
            'links_file': links_path,
            'perturb_genes': ', '.join(self._parse_genes()),
            'cluster_key': cluster_key,
            'embedding_name': str(self.params.get('embedding_name', 'X_umap') or 'X_umap'),
            'counts_layer': str(self.params.get('counts_layer', 'counts') or ''),
            'use_hvg': bool(self.params.get('use_hvg', True)),
            'max_cells': int(self.params.get('max_cells', 30000) or 0),
            'n_comps': int(self.params.get('n_comps', 0) or 0),
            'n_comps_max': int(self.params.get('n_comps_max', 50) or 50),
            'knn_k': int(self.params.get('knn_k', 0) or 0),
            'n_jobs': int(self.params.get('n_jobs', 4) or 1),
            'grn_unit': str(self.params.get('grn_unit', 'cluster') or 'cluster'),
            'alpha': float(self.params.get('alpha', 10) or 10),
            'edge_coef_cutoff': float(self.params.get('edge_coef_cutoff', 0) or 0),
            'top_edges_per_cluster': int(self.params.get('top_edges_per_cluster', 500) or 100),
            'n_propagation': n_propagation,
            'n_neighbors': int(self.params.get('n_neighbors', 200) or 50),
            'min_mass': float(self.params.get('min_mass', 0.01) or 0.01),
            'n_grid': int(self.params.get('n_grid', 40) or 20),
            'top_regulated_genes': int(self.params.get('top_regulated_genes', 50) or 10),
            'combine_perturbations': bool(self.params.get('combine_perturbations', False)),
            'save_oracle_object': bool(self.params.get('save_oracle_object', False)),
            # 默认拒绝把 log1p 归一化 X 当作 counts 建模；仅当用户在界面
            # 显式确认时才允许回退（结果不可靠，需在报告中标注）。
            'allow_non_count_fallback': bool(
                self.params.get('allow_non_count_fallback', False)
            ),
            # 随机对照（knn_random / 随机流场）的固定种子，保证可复现。
            'random_seed': int(self.params.get('random_seed', 0) or 0),
        }
        return config

    def _run_subprocess(self, config):
        """运行 celloracle worker，流式转发进度。返回 (files, summary)。"""
        python = self._celloracle_python()
        worker = self._celloracle_worker()
        if not python or not os.path.isfile(python):
            return None, {
                'status': 'unavailable',
                'error': ('CellOracle 环境未安装（找不到 ' + python + '）。'
                          '请在服务器上创建独立的 celloracle Python 3.9/3.10 '
                          '环境并在 config.py 中配置 CELLORACLE_PYTHON。'),
                'missing_dependencies': ['celloracle'],
            }
        if not worker or not os.path.isfile(worker):
            return None, {
                'status': 'unavailable',
                'error': '找不到 celloracle worker 脚本: ' + worker,
                'missing_dependencies': ['celloracle_worker'],
            }

        from config import Config
        runtime_dir = os.path.join(self.project_dir, 'runtime')
        os.makedirs(runtime_dir, exist_ok=True)
        job_path = os.path.join(runtime_dir,
                                'virtual_ko_job_' + uuid.uuid4().hex[:8] + '.json')
        stderr_path = os.path.join(runtime_dir,
                                   'virtual_ko_stderr_' + uuid.uuid4().hex[:8] + '.log')
        with open(job_path, 'w', encoding='utf-8') as fh:
            json.dump(config, fh, ensure_ascii=False, indent=2)

        home_dir = os.path.abspath(getattr(
            Config, 'CELLORACLE_HOME_DIR',
            os.path.join(Config.RUNTIME_TMP_DIR, 'celloracle_home')))
        os.makedirs(home_dir, exist_ok=True)
        # Do not inherit the Flask worker's full environment.  In particular,
        # .env secrets such as platform passwords and AI credentials must not
        # become readable by a process that handles user-controlled inputs.
        env_keys = (
            'PATH', 'LANG', 'LC_ALL', 'LC_CTYPE', 'TZ',
            'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_MAX_THREADS', 'CUDA_VISIBLE_DEVICES',
        )
        env = {key: os.environ[key] for key in env_keys if os.environ.get(key)}
        env.update({
            'HOME': home_dir,
            'XDG_CONFIG_HOME': os.path.join(home_dir, 'config'),
            'MPLCONFIGDIR': os.path.join(home_dir, 'mplconfig'),
            'NUMBA_CACHE_DIR': os.path.join(home_dir, 'numba_cache'),
            'MPLBACKEND': 'Agg',
        })
        for sub in ('config', 'mplconfig', 'numba_cache'):
            os.makedirs(os.path.join(home_dir, sub), exist_ok=True)

        cmd = [python, worker, '--config', job_path]
        self.progress(1, '启动 CellOracle 子进程: ' + python)
        logger.info('[virtual_ko] %s', ' '.join(cmd))
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=env, text=True, bufsize=1, start_new_session=True,
        )

        def terminate_process_group():
            """Terminate CellOracle and any children it spawned."""
            if proc.poll() is not None:
                return
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                try:
                    proc.kill()
                except ProcessLookupError:
                    pass

        files = []
        summary = {'status': 'completed'}
        error_message = ''
        stderr_tail = []
        err_fh = open(stderr_path, 'w', encoding='utf-8')

        def drain_stderr():
            # 在独立线程里持续排空 stderr（CellOracle 的 tqdm 进度条会写
            # stderr），避免管道缓冲写满导致父子进程互相等待。
            try:
                for err_line in proc.stderr:
                    err_fh.write(err_line)
                    stderr_tail.append(err_line)
                    if len(stderr_tail) > 60:
                        stderr_tail.pop(0)
            except Exception:
                pass

        import threading
        drain_thread = threading.Thread(target=drain_stderr, daemon=True)
        drain_thread.start()

        # 看门狗：CellOracle 推断可能长时间卡住（KNN 死锁、base GRN 首次
        # 下载网络挂起），若不设上限会永久占用线程池的 2 个工作线程。
        # 默认 6 小时，可用 CELLORACLE_TIMEOUT_SECONDS 覆盖。
        import os as _os
        timeout_seconds = float(_os.environ.get('CELLORACLE_TIMEOUT_SECONDS', '21600'))
        timed_out = []

        def _timeout_killer():
            timed_out.append(True)
            terminate_process_group()

        watchdog = threading.Timer(timeout_seconds, _timeout_killer)
        watchdog.daemon = True
        watchdog.start()
        try:
            for line in proc.stdout:
                line = line.strip()
                if not line:
                    continue
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    logger.debug('[virtual_ko] 非 JSON 输出: %s', line[:200])
                    continue
                etype = event.get('type')
                if etype == 'progress':
                    self.progress(int(event.get('pct', 0)),
                                  str(event.get('message', '')))
                elif etype == 'result':
                    files = event.get('files') or []
                    summary = event.get('summary') or {'status': 'completed'}
                elif etype == 'error':
                    error_message = (event.get('message') or '') + chr(10) + \
                                    (event.get('traceback') or '')
            ret = proc.wait()
            drain_thread.join(timeout=10)
        finally:
            watchdog.cancel()
            terminate_process_group()
            err_fh.close()
            try:
                os.remove(job_path)
            except OSError:
                pass

        if timed_out:
            raise RuntimeError(
                'CellOracle 虚拟敲除超过时限（%.0f 秒）已强制终止；'
                '可通过环境变量 CELLORACLE_TIMEOUT_SECONDS 调整。' % timeout_seconds)
        if error_message or ret != 0:
            tail = ''.join(stderr_tail[-30:])
            raise RuntimeError(
                'CellOracle 虚拟敲除失败（exit=' + str(ret) + '）:' + chr(10)
                + (error_message or '未知错误') + chr(10)
                + '--- stderr tail ---' + chr(10) + tail[-4000:])
        return files, summary

    # ---------- 框架接口 ----------

    def validate_input(self, adata):
        if 'X_umap' not in adata.obsm and 'X_draw_graph_fa' not in adata.obsm:
            if 'X_tsne' not in adata.obsm:
                return ('缺少 2D embedding（obsm 中没有 X_umap 等键）。'
                        '请先运行 dimred 生成 UMAP。')
        cluster_key = str(self.params.get('cluster_key', '') or '').strip()
        if cluster_key and cluster_key in adata.obs.columns:
            return None
        for cand in CLUSTER_KEY_CANDIDATES:
            if cand in adata.obs.columns:
                return None
        return ('缺少聚类/注释列（obs 中没有 celltype/final_annotation/leiden 等）。'
                '请先运行 clustering 或 annotation 模块。')

    def run(self, input_path):
        adata = self.load_adata(input_path)
        adata = self.apply_scope(adata)
        # counts 可用性检查（worker 内部还会再校验一次并给出 CellOracle 侧信息）
        has_counts = any(l in adata.layers for l in ('counts', 'raw_count'))
        if not has_counts:
            # worker 默认拒绝在非 counts 的 X 上推断 GRN，且该回退开关
            # 未在界面暴露；这里必须如实告知会失败，而不是承诺回退。
            logger.warning(
                '[virtual_ko] 输入 h5ad 没有 counts layer；'
                'worker 将拒绝运行（allow_non_count_fallback 未启用）。'
                '请先运行保留 counts 层的 QC/标准化流程。'
            )

        genes = self._parse_genes()
        if not genes:
            raise ValueError('请至少填写一个扰动基因（如 SPI1、GATA1）。')
        missing = [g for g in genes if g not in adata.var_names]
        if missing:
            raise ValueError('扰动基因不在表达矩阵中: %s。'
                             '请检查基因名（区分大小写）。' % ', '.join(missing))

        # CellOracle runs in a separate Python environment and therefore reads
        # the path from its job JSON.  Passing the original path here would
        # silently discard the in-memory scope applied above.  Materialise a
        # private scoped h5ad and remove it only after the worker exits.
        analysis_input_path = input_path
        temporary_input_path = None
        scope_key = str(self.params.get('scope_key', '') or '').strip()
        if scope_key:
            runtime_dir = os.path.join(self.project_dir, 'runtime')
            os.makedirs(runtime_dir, exist_ok=True)
            temporary_input_path = os.path.join(
                runtime_dir, f'virtual_ko_input_{uuid.uuid4().hex}.h5ad'
            )
            adata.write_h5ad(temporary_input_path)
            analysis_input_path = temporary_input_path

        try:
            config = self._build_job_config(analysis_input_path)
            self.progress(1, '准备 CellOracle 虚拟敲除任务...')
            files, summary = self._run_subprocess(config)
        finally:
            if temporary_input_path:
                try:
                    os.remove(temporary_input_path)
                except OSError:
                    pass

        if summary.get('status') == 'unavailable':
            return {
                'output_adata': input_path,
                'result_files': [],
                'summary': summary,
            }

        # 登记 worker 输出的所有文件
        result_files = []
        for item in files or []:
            path = item.get('path')
            if not path or not os.path.isfile(path):
                continue
            kind = item.get('kind', 'txt')
            file_type = kind if kind in ('csv', 'xlsx', 'png', 'svg', 'pdf',
                                         'tiff', 'jpg', 'jpeg', 'json', 'txt',
                                         'h5ad', 'oracle') else 'txt'
            result_files.append({
                'file_path': path,
                'file_type': file_type,
                'category': item.get('category', 'plot'),
                'label': item.get('label') or os.path.basename(path),
            })

        output_h5ad = os.path.join(self.project_dir, 'intermediate', 'virtual_ko_output.h5ad')
        if not os.path.isfile(output_h5ad):
            output_h5ad = summary.get('output_h5ad') or output_h5ad
        summary['status'] = 'completed'
        summary['module'] = self.MODULE_NAME
        summary['input_path'] = input_path
        summary['n_result_files'] = len(result_files)
        summary['scope_key'] = str(self.params.get('scope_key', '') or '').strip() or None
        summary['scope_values'] = ([v.strip() for v in str(self.params.get('scope_values', '') or '').split(',') if v.strip()] or None)
        self.progress(100, '虚拟敲除分析完成')
        return {
            'output_adata': output_h5ad,
            'result_files': result_files,
            'summary': summary,
        }
