import os
import re
import tempfile


def _load_dotenv(path=None):
    """Load simple KEY=VALUE pairs from .env without overriding real env vars."""
    env_path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), '.env')
    if not os.path.exists(env_path):
        return
    try:
        with open(env_path, encoding='utf-8') as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line or line.startswith('#') or '=' not in line:
                    continue
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                if not key or key in os.environ:
                    continue
                if (value.startswith('"') and value.endswith('"')) or (
                    value.startswith("'") and value.endswith("'")
                ):
                    value = value[1:-1]
                os.environ[key] = value
    except OSError:
        return


_load_dotenv()


class Config:
    _BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    @staticmethod
    def _load_secret():
        key = os.environ.get('SECRET_KEY')
        if key:
            return key
        secret_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'instance', '.secret_key')
        try:
            with open(secret_file) as f:
                return f.read().strip()
        except FileNotFoundError:
            import secrets
            key = secrets.token_hex(24)
            os.makedirs(os.path.dirname(secret_file), exist_ok=True)
            with open(secret_file, 'w') as f:
                f.write(key)
            os.chmod(secret_file, 0o600)
            return key

    SECRET_KEY = _load_secret.__func__()
    # Optional shared access gate for controlled lab deployments.  Keep this
    # empty for local development/tests; deployments exposed through a tunnel
    # should set a strong value outside the repository (for example in .env).
    PLATFORM_ACCESS_PASSWORD = os.environ.get('PLATFORM_ACCESS_PASSWORD', '')
    try:
        PLATFORM_ACCESS_SESSION_HOURS = max(
            1.0, float(os.environ.get('PLATFORM_ACCESS_SESSION_HOURS', '12'))
        )
    except ValueError:
        PLATFORM_ACCESS_SESSION_HOURS = 12.0
    SESSION_COOKIE_HTTPONLY = True
    SESSION_COOKIE_SAMESITE = 'Lax'
    SESSION_COOKIE_SECURE = os.environ.get(
        'PLATFORM_SESSION_COOKIE_SECURE', 'false'
    ).strip().lower() in {'1', 'true', 'yes', 'on'}
    DATA_DIR = os.environ.get('DATA_DIR', os.path.join(_BASE_DIR, 'data'))
    # Keep transient analysis files on the project data volume by default.
    # ``/tmp`` is commonly mounted on the root filesystem and can fill up when
    # 10x compatibility staging or plotting jobs handle large matrices.
    RUNTIME_TMP_DIR = os.path.abspath(os.environ.get(
        'RUNTIME_TMP_DIR', os.path.join(DATA_DIR, 'runtime_tmp')
    ))
    NUMBA_CACHE_DIR = os.path.abspath(os.environ.get(
        'NUMBA_CACHE_DIR', os.path.join(RUNTIME_TMP_DIR, 'numba_cache')
    ))
    MPLCONFIG_DIR = os.path.abspath(os.environ.get(
        'MPLCONFIGDIR', os.path.join(RUNTIME_TMP_DIR, 'mplconfig')
    ))
    CACHE_DIR = os.path.abspath(os.environ.get(
        'CACHE_DIR', os.path.join(DATA_DIR, 'cache')
    ))
    # Curated public resources used by optional single-cell modules live on
    # the controlled data volume rather than in user project uploads.  The
    # directory can be redirected by an administrator, while tests continue
    # to follow their isolated DATA_DIR.
    FUNCTIONAL_STATE_RESOURCE_DIR = os.environ.get(
        'FUNCTIONAL_STATE_RESOURCE_DIR', ''
    ).strip()
    # The repository carries a small, versioned public baseline so a clean
    # clone can run local pathway and TF-activity analyses without fetching
    # resources at runtime.  An administrator may override it with a
    # separately maintained controlled directory.
    SC_CELL_GO_GENE_SET_DIR = os.environ.get('SC_CELL_GO_GENE_SET_DIR', '').strip()
    DB_PATH = os.environ.get('DB_PATH', os.path.join(_BASE_DIR, 'instance', 'bioinfo.db'))
    CELLMARKER_PATH = os.environ.get('CELLMARKER_PATH', os.path.join(os.path.dirname(_BASE_DIR), 'CellMarker_Augmented_2021.txt'))

    # CellOracle 虚拟敲除（virtual KO）运行环境。
    # CellOracle 依赖 Python 3.8-3.10，而平台主环境为 Python 3.12，
    # 因此虚拟敲除模块通过 subprocess 调用独立的 celloracle 环境执行。
    # 可用环境变量 CELLORACLE_PYTHON / CELLORACLE_WORKER_SCRIPT /
    # CELLORACLE_HOME_DIR 覆盖以下默认值。
    CELLORACLE_PYTHON = os.path.abspath(os.environ.get(
        'CELLORACLE_PYTHON',
        '/home/oelab/AnaData/GJ/celloracle/celloracle_env/bin/python'
    ))
    CELLORACLE_WORKER_SCRIPT = os.path.abspath(os.environ.get(
        'CELLORACLE_WORKER_SCRIPT',
        os.path.join(_BASE_DIR, 'modules', 'celloracle_worker.py')
    ))
    # celloracle 会在 HOME/.config（genomepy）、HOME/celloracle_data 等位置
    # 写缓存，统一放到独立目录，避免污染用户主目录。
    CELLORACLE_HOME_DIR = os.path.abspath(os.environ.get(
        'CELLORACLE_HOME_DIR',
        '/home/oelab/AnaData/GJ/celloracle/celloracle_home'
    ))

    MAX_WORKERS = int(os.environ.get('MAX_WORKERS', '2'))
    CHUNK_SIZE_MB = 50
    PLOTLY_MAX_CELLS = 50000
    MIN_FREE_RAM_GB = float(os.environ.get('MIN_FREE_RAM_GB', '4'))
    CUDA_DEVICES = os.environ.get('CUDA_DEVICES', '0,1')

    # Server-resident 10x data are intentionally opt-in.  A web request may
    # only refer to a directory beneath one of these administrator configured
    # roots; arbitrary absolute paths and symlinks remain disallowed.  Example:
    # SC_BATCH_SOURCE_ROOTS=/home/oelab/data/GJ:/mnt/sc_data
    SC_BATCH_SOURCE_ROOTS = os.environ.get('SC_BATCH_SOURCE_ROOTS', '')

    # Server-resident WES inputs are also opt-in.  FASTQ/BAM/CRAM/VCF files
    # are usually too large for browser upload; only administrators may expose
    # explicitly configured, read-only roots to a project preflight.
    WES_SOURCE_ROOTS = os.environ.get('WES_SOURCE_ROOTS', '')
    # Small user-uploaded capture BEDs are stored in a platform-owned area.
    # The path can be moved to a lab data volume without changing the API.
    WES_UPLOAD_ROOT = os.environ.get('WES_UPLOAD_ROOT', '').strip()
    WES_MAX_CAPTURE_BED_MB = float(os.environ.get('WES_MAX_CAPTURE_BED_MB', '100'))
    # WES execution is opt-in.  Preparing a run and reviewing its bundle stay
    # available in the lab UI even when Nextflow is not installed on this host.
    WES_EXECUTOR_ENABLED = os.environ.get('WES_EXECUTOR_ENABLED', '').strip().lower() in {
        '1', 'true', 'yes', 'on'
    }
    WES_NEXTFLOW_BIN = os.environ.get('WES_NEXTFLOW_BIN', 'nextflow').strip() or 'nextflow'
    # Sarek's built-in Docker profile still uses the local Nextflow executor;
    # administrators may replace it with ``apptainer`` or an institutional
    # profile after validating the container runtime.
    WES_NEXTFLOW_PROFILE = os.environ.get('WES_NEXTFLOW_PROFILE', 'docker').strip() or 'docker'
    WES_NEXTFLOW_PIPELINE = os.environ.get(
        'WES_NEXTFLOW_PIPELINE', 'nf-core/sarek'
    ).strip() or 'nf-core/sarek'
    # This is the administrator's validated Sarek/iGenomes key for the
    # installed reference bundle.  A deployment using a custom FASTA should
    # set the corresponding Sarek genome configuration explicitly.
    WES_NEXTFLOW_GENOME = os.environ.get('WES_NEXTFLOW_GENOME', 'GATK.GRCh38').strip()
    # Optional local Sarek resources.  Keep them unset until the whole bundle
    # has been copied, checksummed and registered; an empty value lets staff
    # prepare/review a run but the production launch gate remains closed.
    WES_NEXTFLOW_IGENOMES_BASE = os.environ.get(
        'WES_NEXTFLOW_IGENOMES_BASE', ''
    ).strip()
    WES_NEXTFLOW_VEP_CACHE = os.environ.get('WES_NEXTFLOW_VEP_CACHE', '').strip()
    WES_NEXTFLOW_PON = os.environ.get('WES_NEXTFLOW_PON', '').strip()
    WES_NEXTFLOW_GERMLINE_RESOURCE = os.environ.get(
        'WES_NEXTFLOW_GERMLINE_RESOURCE', ''
    ).strip()
    # SRA is converted locally before entering the regular FASTQ manifest
    # path.  Keep this separate from Nextflow so an SRA conversion cannot
    # silently start a WES run.
    WES_SRA_FASTERQ_BIN = os.environ.get('WES_SRA_FASTERQ_BIN', 'fasterq-dump').strip() or 'fasterq-dump'
    try:
        WES_SRA_THREADS = max(1, int(os.environ.get('WES_SRA_THREADS', '4')))
    except ValueError:
        WES_SRA_THREADS = 4
    try:
        WES_SRA_MIN_FREE_GB = max(1.0, float(os.environ.get('WES_SRA_MIN_FREE_GB', '50')))
    except ValueError:
        WES_SRA_MIN_FREE_GB = 50.0
    try:
        WES_SRA_DISK_FACTOR = max(2.0, float(os.environ.get('WES_SRA_DISK_FACTOR', '3')))
    except ValueError:
        WES_SRA_DISK_FACTOR = 3.0
    # Real runs must use catalogued, checksum-verified resources.  Tests and
    # administrator-only engineering drills can opt out explicitly on the
    # executor object; the web API never disables this gate.
    WES_REQUIRE_VALIDATED_REFERENCES = os.environ.get(
        'WES_REQUIRE_VALIDATED_REFERENCES', 'true'
    ).strip().lower() in {'1', 'true', 'yes', 'on'}
    WES_MIN_FREE_GB = float(os.environ.get('WES_MIN_FREE_GB', '50'))
    WES_HAPPY_BIN = os.environ.get('WES_HAPPY_BIN', 'hap.py').strip() or 'hap.py'
    WES_SOMPY_BIN = os.environ.get('WES_SOMPY_BIN', 'som.py').strip() or 'som.py'
    WES_NEXTFLOW_POLL_SECONDS = float(os.environ.get('WES_NEXTFLOW_POLL_SECONDS', '5'))
    WES_CANCEL_GRACE_SECONDS = int(os.environ.get('WES_CANCEL_GRACE_SECONDS', '20'))

    # AI 对话配置（支持 OpenAI compatible 和 Anthropic messages compatible API）
    AI_API_KEY = os.environ.get('AI_API_KEY', '')
    AI_API_URL = os.environ.get('AI_API_URL', 'https://token-plan-cn.xiaomimimo.com/anthropic')
    AI_MODEL = os.environ.get('AI_MODEL', 'mimo-v2.5-pro')
    AI_API_TOKEN = os.environ.get('AI_API_TOKEN', '')

    _PID_RE = re.compile(r'^[a-zA-Z0-9_-]{1,64}$')

    @classmethod
    def _validate_pid(cls, pid):
        if not cls._PID_RE.match(pid):
            raise ValueError(f"Invalid project ID: {pid!r}")
        return pid

    @classmethod
    def project_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid)

    @classmethod
    def uploads_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'uploads')

    @classmethod
    def results_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'results')

    @classmethod
    def intermediate_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'intermediate')

    @classmethod
    def plots_dir(cls, pid):
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'plots')

    @classmethod
    def runtime_tmp_dir(cls):
        """Return a writable, project-owned directory for disposable files."""
        path = os.path.abspath(cls.RUNTIME_TMP_DIR)
        os.makedirs(path, exist_ok=True)
        return path

    @classmethod
    def functional_state_resource_dir(cls):
        """Return an administrator override, local snapshot, or bundled baseline."""
        configured = str(cls.FUNCTIONAL_STATE_RESOURCE_DIR or '').strip()
        if configured:
            return os.path.abspath(configured)
        local_root = os.path.abspath(os.path.join(cls.DATA_DIR, 'functional_state_resources'))
        if (
            os.path.isfile(os.path.join(local_root, 'gene_sets', 'gene_set_registry.json'))
            or os.path.isfile(os.path.join(local_root, 'collectri_human.tsv'))
        ):
            return local_root
        return os.path.abspath(os.path.join(cls._BASE_DIR, 'resources', 'functional_state_resources'))

    @classmethod
    def sc_cell_go_gene_set_dir(cls):
        """Return the local override or bundled GMT directory for ``sc_cell_go``."""
        configured = str(cls.SC_CELL_GO_GENE_SET_DIR or '').strip()
        if configured:
            return os.path.abspath(configured)
        local_root = os.path.abspath(os.path.join(cls.DATA_DIR, 'go_gene_sets'))
        if os.path.isdir(local_root) and any(
            name.lower().endswith(('.gmt', '.txt')) for name in os.listdir(local_root)
        ):
            return local_root
        functional_root = os.path.join(cls.functional_state_resource_dir(), 'gene_sets')
        if os.path.isfile(os.path.join(functional_root, 'gene_set_registry.json')):
            return functional_root
        return os.path.join(cls._BASE_DIR, 'resources', 'functional_state_resources', 'gene_sets')

    @classmethod
    def configure_runtime_tmpdir(cls):
        """Route Python and common Unix temporary-file lookups away from ``/tmp``.

        This runs while ``config`` is imported, before analysis libraries are
        loaded.  ``RUNTIME_TMP_DIR`` remains an explicit deployment override
        for installations that keep data on a separate mounted volume.
        """
        path = cls.runtime_tmp_dir()
        for variable in ('TMPDIR', 'TMP', 'TEMP'):
            os.environ[variable] = path
        for variable, cache_dir in (
            ('NUMBA_CACHE_DIR', cls.NUMBA_CACHE_DIR),
            ('MPLCONFIGDIR', cls.MPLCONFIG_DIR),
        ):
            os.makedirs(cache_dir, exist_ok=True)
            os.environ[variable] = cache_dir
        # ``tempfile`` caches its selected directory, so update it explicitly
        # as well as the environment variables.
        tempfile.tempdir = path
        return path

    @classmethod
    def branches_dir(cls, pid):
        """候选分支输出目录，隔离于主线 intermediate/."""
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'branches')

    @classmethod
    def branch_dir(cls, pid, branch_id):
        """单个候选分支目录."""
        cls._validate_pid(pid)
        return os.path.join(cls.DATA_DIR, 'projects', pid, 'branches', branch_id)

    @classmethod
    def _validate_path(cls, path, pid):
        """验证路径在项目目录内，拒绝符号链接。返回 (is_valid: bool, error: str|None)."""
        cls._validate_pid(pid)
        proj_dir = os.path.realpath(cls.project_dir(pid))
        try:
            real_path = os.path.realpath(path)
        except (OSError, ValueError) as e:
            return False, f"路径解析失败: {e}"
        if os.path.islink(path):
            return False, "不支持符号链接文件"
        if not real_path.startswith(proj_dir + os.sep) and real_path != proj_dir:
            return False, f"路径不在项目目录内: {path}"
        return True, None

    @classmethod
    def sc_batch_source_roots(cls):
        """Return canonical administrator-approved roots for batch 10x input.

        Read the environment on each call so a long-running development server
        can be configured without importing arbitrary filesystem paths into the
        application at module import time.  Empty by default is the safe mode.
        """
        raw = os.environ.get('SC_BATCH_SOURCE_ROOTS', cls.SC_BATCH_SOURCE_ROOTS)
        roots = []
        for item in str(raw or '').split(os.pathsep):
            item = item.strip()
            if not item:
                continue
            try:
                root = os.path.realpath(item)
            except (OSError, ValueError):
                continue
            if os.path.isdir(root) and root not in roots:
                roots.append(root)
        return tuple(roots)

    @classmethod
    def validate_sc_batch_source_path(cls, path, *, require_directory=True):
        """Validate a server-side batch source against ``SC_BATCH_SOURCE_ROOTS``.

        This is deliberately separate from ``_validate_path``: batch input is
        read-only and can live on a mounted data volume, while all platform
        outputs still remain beneath the project directory.
        """
        if not path:
            raise ValueError('缺少服务器数据目录')
        roots = cls.sc_batch_source_roots()
        if not roots:
            raise ValueError(
                '服务器目录批量导入尚未启用；管理员需在 .env 中设置 '
                'SC_BATCH_SOURCE_ROOTS（例如 /home/oelab/data/GJ）'
            )
        try:
            supplied = os.path.abspath(str(path))
            resolved = os.path.realpath(supplied)
        except (OSError, ValueError) as exc:
            raise ValueError(f'服务器目录解析失败: {exc}') from exc
        if os.path.islink(supplied):
            raise ValueError('服务器批量导入不接受符号链接路径')
        allowed = any(resolved == root or resolved.startswith(root + os.sep)
                      for root in roots)
        if not allowed:
            raise ValueError('服务器目录不在允许的 SC_BATCH_SOURCE_ROOTS 下')
        if require_directory and not os.path.isdir(resolved):
            raise ValueError(f'服务器目录不存在或不是目录: {path}')
        return resolved

    @classmethod
    def wes_source_roots(cls):
        """Return canonical administrator-approved roots for WES inputs."""
        raw = os.environ.get('WES_SOURCE_ROOTS', cls.WES_SOURCE_ROOTS)
        roots = []
        for item in str(raw or '').split(os.pathsep):
            item = item.strip()
            if not item:
                continue
            try:
                root = os.path.realpath(item)
            except (OSError, ValueError):
                continue
            if os.path.isdir(root) and root not in roots:
                roots.append(root)
        upload_root = cls.wes_upload_root()
        if os.path.isdir(upload_root) and upload_root not in roots:
            roots.append(upload_root)
        return tuple(roots)

    @classmethod
    def wes_upload_root(cls):
        """Return the writable, platform-owned root for capture BED uploads."""
        configured = os.environ.get('WES_UPLOAD_ROOT', cls.WES_UPLOAD_ROOT).strip()
        return os.path.abspath(configured or os.path.join(cls.DATA_DIR, 'wes_uploads'))

    @classmethod
    def validate_wes_source_path(cls, path, *, require_file=True):
        """Validate a read-only WES path under ``WES_SOURCE_ROOTS``."""
        if not path:
            raise ValueError('缺少 WES 服务器数据路径')
        roots = cls.wes_source_roots()
        if not roots:
            raise ValueError(
                'WES 服务器目录尚未启用；管理员需在 .env 中设置 WES_SOURCE_ROOTS'
            )
        try:
            supplied = os.path.abspath(str(path))
            resolved = os.path.realpath(supplied)
        except (OSError, ValueError) as exc:
            raise ValueError(f'WES 服务器路径解析失败: {exc}') from exc
        if os.path.islink(supplied):
            raise ValueError('WES 服务器数据不接受符号链接路径')
        if not any(resolved == root or resolved.startswith(root + os.sep) for root in roots):
            raise ValueError('WES 服务器路径不在允许的 WES_SOURCE_ROOTS 下')
        if require_file and not os.path.isfile(resolved):
            raise ValueError(f'WES 文件不存在或不是普通文件: {path}')
        return resolved

    @classmethod
    def validate_wes_reference_path(cls, path, *, require_directory=False):
        """Validate one configured reference file/directory without following links."""
        resolved = cls.validate_wes_source_path(path, require_file=False)
        if require_directory:
            if not os.path.isdir(resolved):
                raise ValueError(f'WES 参考资源不存在或不是目录: {path}')
        elif not os.path.isfile(resolved):
            raise ValueError(f'WES 参考资源不存在或不是普通文件: {path}')
        return resolved


# Apply the non-root temporary directory before Scanpy, Matplotlib, or worker
# modules can initialize their own caches.
Config.configure_runtime_tmpdir()
