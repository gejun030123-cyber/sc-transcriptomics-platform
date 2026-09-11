import csv
import json
import os
import shutil
import struct
import base64
import uuid
import logging
import psutil
from collections.abc import Mapping
from flask import Blueprint, jsonify, request, send_file
from models import AnalysisTask, Project, ResultFile
from config import Config
from modules.workflows.capture_uploads import store_capture_bed_upload

logger = logging.getLogger(__name__)

api_bp = Blueprint('api', __name__)


BUNDLED_PRESETS_DIR = os.path.join(Config._BASE_DIR, 'resources', 'presets')


@api_bp.route('/projects/<pid>/wes/capture-kits/upload', methods=['POST'])
def upload_wes_capture_kit(pid):
    """Register a user-provided calling BED as a non-production profile."""
    if not Project.get_by_id(pid):
        return jsonify({'error': 'Project not found'}), 404
    uploaded = request.files.get('calling_bed') or request.files.get('bed') or request.files.get('file')
    try:
        result = store_capture_bed_upload(
            project_id=pid,
            file_storage=uploaded,
            capture_kit_id=request.form.get('capture_kit_id', ''),
            name=request.form.get('name', ''),
            version=request.form.get('version', ''),
            assembly=request.form.get('assembly', 'GRCh38'),
        )
    except (OSError, ValueError) as exc:
        return jsonify({'error': str(exc)}), 400
    return jsonify({
        'message': 'BED 已上传并登记为 test_only，需管理员审核后才能用于生产分析',
        'profile': result['profile'],
        'checksum': result['checksum'],
        'size_bytes': result['size_bytes'],
    }), 201


def _validate_file_path(file_path):
    """校验文件路径在允许的目录内，防止路径遍历攻击"""
    if not file_path:
        return False
    if os.path.islink(file_path):
        return False
    try:
        real_path = os.path.realpath(file_path)
        data_dir = os.path.realpath(Config.DATA_DIR)
    except (OSError, ValueError):
        return False
    return real_path.startswith(data_dir + os.sep) or real_path == data_dir


def _validate_project_file_path(file_path, project_id):
    if not file_path or not project_id:
        return False
    is_valid, _ = Config._validate_path(file_path, project_id)
    return is_valid and os.path.isfile(file_path)


def _get_project_presets_dir(project_id):
    if '..' in project_id or '/' in project_id:
        raise ValueError('无效的项目 ID')
    return Config.project_dir(project_id) + '/presets'


def _get_global_presets_dir():
    """Resolve the global preset directory at call time.

    Tests and controlled deployments may change ``Config.DATA_DIR`` after this
    module is imported.  Resolving lazily prevents the bundled presets from
    being copied into a stale data root.
    """
    return os.path.join(Config.DATA_DIR, 'presets', '_global')


def _ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def _ensure_bundled_presets():
    """Make version-controlled default presets available in the data store.

    A scientist may customise or remove files in the global preset directory,
    so bundled defaults never overwrite an existing file or a preset with the
    same display name.  Copying them to the normal data-backed directory keeps
    the existing list/load/run APIs as the single preset contract.
    """
    if not os.path.isdir(BUNDLED_PRESETS_DIR):
        return
    target_dir = _get_global_presets_dir()
    _ensure_dir(target_dir)

    existing_names = set()
    for filename in os.listdir(target_dir):
        if not filename.endswith('.json'):
            continue
        existing = _load_preset(os.path.join(target_dir, filename))
        if existing and str(existing.get('name', '') or '').strip():
            existing_names.add(str(existing['name']).strip().casefold())

    for filename in sorted(os.listdir(BUNDLED_PRESETS_DIR)):
        if not filename.endswith('.json'):
            continue
        source = os.path.join(BUNDLED_PRESETS_DIR, filename)
        bundled = _load_preset(source)
        if not bundled:
            logger.warning('忽略无效的内置预设文件: %s', filename)
            continue
        name = str(bundled.get('name', '') or '').strip()
        target = os.path.join(target_dir, filename)
        if os.path.exists(target) or (name and name.casefold() in existing_names):
            continue
        temporary = f'{target}.{uuid.uuid4().hex}.tmp'
        try:
            shutil.copyfile(source, temporary)
            os.replace(temporary, target)
            if name:
                existing_names.add(name.casefold())
        except OSError as exc:
            logger.warning('无法初始化内置预设 %s: %s', filename, exc)
            try:
                if os.path.exists(temporary):
                    os.remove(temporary)
            except OSError:
                pass


def _load_preset(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        data['_filepath'] = filepath
        data['_id'] = os.path.splitext(os.path.basename(filepath))[0]
        return data
    except (json.JSONDecodeError, IOError):
        return None


def _load_accessible_preset(preset_id, project_id=''):
    """Load one project-local preset first, then a shared preset.

    Pipeline templates are server-owned JSON records.  The one-click endpoint
    intentionally resolves them here instead of trusting a browser copy of
    module names or parameter values.
    """
    if not preset_id or '..' in preset_id or '/' in preset_id:
        return None
    if project_id:
        preset = _load_preset(os.path.join(
            _get_project_presets_dir(project_id), f'{preset_id}.json'
        ))
        if preset:
            preset['_scope'] = 'project'
            return preset
    _ensure_bundled_presets()
    preset = _load_preset(os.path.join(_get_global_presets_dir(), f'{preset_id}.json'))
    if preset:
        preset['_scope'] = 'global'
    return preset


def _list_presets_in_dir(directory, scope):
    if not os.path.isdir(directory):
        return []
    presets = []
    for fname in sorted(os.listdir(directory)):
        if not fname.endswith('.json'):
            continue
        p = _load_preset(os.path.join(directory, fname))
        if p:
            p['_scope'] = scope
            presets.append(p)
    return presets


import math


def _sanitize_plotly_values(obj):
    """Replace NaN/Inf with None (JSON null) for standard-compliant JSON."""
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    if isinstance(obj, list):
        return [_sanitize_plotly_values(v) for v in obj]
    if isinstance(obj, dict):
        return {k: _sanitize_plotly_values(v) for k, v in obj.items()}
    return obj


def _decode_plotly_binary(obj):
    """Convert Plotly base64 binary fields to plain lists for Plotly.js compatibility."""
    if isinstance(obj, dict):
        if 'bdata' in obj and 'dtype' in obj:
            dtype_map = {'f4': 'f', 'f8': 'd', 'i1': 'b', 'i2': 'h', 'i4': 'i', 'i8': 'q', 'u1': 'B', 'u2': 'H', 'u4': 'I', 'u8': 'Q'}
            fmt = dtype_map.get(obj['dtype'], 'd')
            decoded = base64.b64decode(obj['bdata'])
            count = len(decoded) // struct.calcsize(fmt)
            return list(struct.unpack(f'<{count}{fmt}', decoded))
        return {k: _decode_plotly_binary(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_decode_plotly_binary(item) for item in obj]
    return obj

@api_bp.route('/tasks/<task_id>/status')
def task_status(task_id):
    t = AnalysisTask.get_by_id(task_id)
    if not t:
        return jsonify({'error': 'Task not found'}), 404
    log_entries = []
    if t.log_text:
        try:
            log_entries = json.loads(t.log_text)
        except Exception:
            pass
    result = {}
    if t.result_json:
        try:
            result = json.loads(t.result_json)
        except Exception:
            pass
    return jsonify({
        'status': t.status,
        'progress': t.progress,
        'progress_message': t.progress_message,
        'error': t.error_traceback,
        'started_at': t.started_at,
        'finished_at': t.finished_at,
        'result': result,
        'log': log_entries,
    })

@api_bp.route('/projects/<pid>/tasks')
def project_tasks(pid):
    tasks = AnalysisTask.get_by_project(pid)
    return jsonify([t.to_dict() for t in tasks])

@api_bp.route('/system/status')
def system_status():
    disk = shutil.disk_usage(Config.DATA_DIR)
    mem = psutil.virtual_memory()
    from worker import active_count
    return jsonify({
        'disk_free_gb': round(disk.free / (1024**3), 1),
        'ram_avail_gb': round(mem.available / (1024**3), 1),
        'ram_percent': mem.percent,
        'active_tasks': active_count(),
    })


@api_bp.route('/system/dependencies')
def system_dependencies():
    from modules.platform.system_health import dependency_status
    return jsonify(dependency_status())


@api_bp.route('/projects/<pid>/adata-info')
def adata_info(pid):
    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'error': 'Not found'}), 404
    adata_path = p.get_latest_adata_path()
    if not adata_path:
        uploads_dir = Config.uploads_dir(pid)
        if os.path.isdir(uploads_dir):
            for f in os.listdir(uploads_dir):
                if f.endswith('.h5ad'):
                    adata_path = os.path.join(uploads_dir, f)
                    break
    if not adata_path or not os.path.exists(adata_path):
        return jsonify({'error': 'No h5ad file found'}), 404
    try:
        import anndata
        from modules.inspect_utils import inspect_adata
        adata = anndata.read_h5ad(adata_path, backed='r')
        info = inspect_adata(adata, include_samples=True)
        info['n_obs'] = int(adata.n_obs)
        info['n_vars'] = int(adata.n_vars)
        info['file'] = os.path.basename(adata_path)
        adata.file.close()
        return jsonify(info)
    except Exception as e:
        logger.exception("API error")
        return jsonify({'error': str(e)}), 500


_PIPELINE_INPUT_EXTENSIONS = {
    'sc': ('.h5ad',),
    'bulk': ('.h5ad', '.csv', '.txt', '.tsv', '.xlsx', '.xls'),
}


def _list_pipeline_inputs(pid, analysis_type):
    """Return safe, project-owned files that can be chosen as pipeline input.

    The launcher must submit an absolute project path to the worker.  Do not
    derive that path in the browser from a display filename: aside from being
    unreliable, that would invite a user-controlled server path.  Files are
    deliberately limited to direct uploads plus registered completed outputs.
    """
    allowed_extensions = _PIPELINE_INPUT_EXTENSIONS[analysis_type]
    entries = []
    seen_paths = set()

    def add(path, label, source):
        path = os.path.abspath(str(path or ''))
        if not path or path in seen_paths:
            return
        if not path.lower().endswith(allowed_extensions):
            return
        valid, _ = Config._validate_path(path, pid)
        if not valid or os.path.islink(path) or not os.path.isfile(path):
            return
        seen_paths.add(path)
        entries.append({
            'path': path,
            'label': label,
            'source': source,
        })

    uploads_dir = Config.uploads_dir(pid)
    if os.path.isdir(uploads_dir):
        for filename in sorted(os.listdir(uploads_dir), key=str.lower):
            add(
                os.path.join(uploads_dir, filename),
                f'上传文件 · {filename}',
                'upload',
            )

    # A completed h5ad can be used when the scientist deliberately wants to
    # rerun a compatible template from an earlier in-project result.  It is
    # still checked as a regular, non-symlink project file above.
    for task in AnalysisTask.get_by_project(pid):
        if task.status != 'completed' or not task.output_adata_path:
            continue
        add(
            task.output_adata_path,
            f'已完成输出 · {os.path.basename(task.output_adata_path)}',
            'task_output',
        )
    return entries


@api_bp.route('/projects/<pid>/pipeline-inputs')
def list_pipeline_inputs(pid):
    """List selectable, server-validated input files for an SC/Bulk pipeline."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    analysis_type = str(request.args.get('type', '') or '').strip().lower()
    if analysis_type not in _PIPELINE_INPUT_EXTENSIONS:
        return jsonify({'error': 'type 必须是 sc 或 bulk'}), 400
    return jsonify({'files': _list_pipeline_inputs(pid, analysis_type)})


@api_bp.route('/projects/<pid>/result-file/<file_id>')
@api_bp.route('/result-file/<file_id>')
def get_result_file(file_id, pid=None):
    f = ResultFile.get_by_id(file_id)
    if not f:
        return jsonify({'error': 'Not found'}), 404
    pid = pid or request.args.get('project_id', '').strip()
    if not pid:
        return jsonify({'error': 'project_id required'}), 400
    if f.project_id != pid:
        return jsonify({'error': 'Result file 不属于该项目'}), 403
    if not _validate_project_file_path(f.file_path, f.project_id):
        return jsonify({'error': '文件路径不在所属项目内'}), 403
    if f.file_type == 'plotly_json':
        return jsonify({'error': 'Plotly interactive results have been retired; rerun the task for PNG/SVG output.'}), 410
    return send_file(f.file_path)


@api_bp.route('/projects/<pid>/result-file/<file_id>/table-preview')
def result_table_preview(pid, file_id):
    """Return a bounded, filterable preview of a project-owned CSV artifact."""
    result_file = ResultFile.get_by_id(file_id)
    if not result_file:
        return jsonify({'error': 'Not found'}), 404
    if result_file.project_id != pid:
        return jsonify({'error': 'Result file 不属于该项目'}), 403
    if result_file.file_type != 'csv':
        return jsonify({'error': '仅支持 CSV 结果表预览'}), 400
    if not _validate_project_file_path(result_file.file_path, pid):
        return jsonify({'error': '文件路径不在所属项目内'}), 403

    def bounded_int(raw, default, minimum, maximum):
        try:
            return max(minimum, min(maximum, int(raw)))
        except (TypeError, ValueError):
            return default

    max_scan = 5000
    limit = bounded_int(request.args.get('limit'), 200, 1, 1000)
    search = str(request.args.get('search', '') or '').strip()
    filter_column = str(request.args.get('column', '') or '').strip()
    min_value = request.args.get('min')
    max_value = request.args.get('max')
    try:
        import numpy as np
        import pandas as pd

        frame = pd.read_csv(result_file.file_path, nrows=max_scan)
        scanned_rows = int(len(frame))
        if search:
            mask = frame.astype(str).apply(
                lambda column: column.str.contains(search, case=False, na=False, regex=False)
            ).any(axis=1)
            frame = frame.loc[mask]
        if filter_column:
            if filter_column not in frame.columns:
                return jsonify({'error': f'列不存在: {filter_column}'}), 400
            numeric = pd.to_numeric(frame[filter_column], errors='coerce')
            if min_value not in (None, ''):
                try:
                    frame = frame.loc[numeric >= float(min_value)]
                    numeric = numeric.loc[frame.index]
                except ValueError:
                    return jsonify({'error': '最小值必须是数字'}), 400
            if max_value not in (None, ''):
                try:
                    frame = frame.loc[numeric <= float(max_value)]
                except ValueError:
                    return jsonify({'error': '最大值必须是数字'}), 400
        matched_rows = int(len(frame))
        preview = frame.head(limit).replace({np.nan: None})
        # JSON serializers do not consistently handle numpy scalar types.
        records = json.loads(preview.to_json(orient='records', force_ascii=False))
        return jsonify({
            'file_id': result_file.id,
            'label': result_file.label,
            'columns': [str(column) for column in frame.columns],
            'rows': records,
            'matched_rows': matched_rows,
            'returned_rows': int(len(records)),
            'scanned_rows': scanned_rows,
            'scan_limit': max_scan,
            'truncated': scanned_rows >= max_scan,
        })
    except UnicodeDecodeError:
        return jsonify({'error': 'CSV 编码无法读取，请直接下载文件查看'}), 400
    except Exception as exc:
        logger.exception('Result table preview failed')
        return jsonify({'error': f'表格预览失败: {exc}'}), 500


@api_bp.route('/projects/<pid>/enrichment-result/<task_id>')
@api_bp.route('/enrichment-result/<task_id>')
def enrichment_result(task_id, pid=None):
    """返回富集分析的静态结果文件元数据。"""
    from models import ResultFile
    task = AnalysisTask.get_by_id(task_id)
    if not task:
        return jsonify({'error': 'Task not found'}), 404
    pid = pid or request.args.get('project_id', '').strip()
    if not pid:
        return jsonify({'error': 'project_id required'}), 400
    if task.project_id != pid:
        return jsonify({'error': 'Task 不属于该项目'}), 403
    files = ResultFile.get_by_task(task_id)
    enrichment_files = [f for f in files if f.category == 'enrichment' and f.file_type in {'png', 'svg', 'pdf', 'tiff'}]
    result = []
    for f in enrichment_files:
        if f.project_id != task.project_id or not _validate_project_file_path(f.file_path, f.project_id):
            return jsonify({'error': '文件路径不在所属项目内'}), 403
        result.append({'id': f.id, 'label': f.label, 'file_type': f.file_type,
                       'url': f'/api/projects/{task.project_id}/result-file/{f.id}'})
    return jsonify(result)


@api_bp.route('/column-values')
def column_values():
    file_path = request.args.get('file_path', '')
    column = request.args.get('column', '')
    if not file_path or not column:
        return jsonify({'values': []})
    if not _validate_file_path(file_path):
        return jsonify({'error': '文件路径不在允许范围内'}), 403
    try:
        from modules.io_utils import read_expression_matrix, infer_expression_measurement
        adata = read_expression_matrix(file_path)
        if column in adata.obs.columns:
            values = sorted(adata.obs[column].astype(str).unique().tolist())
        else:
            values = []
        return jsonify({'values': values})
    except Exception:
        return jsonify({'values': []})


@api_bp.route('/projects/<pid>/design-preflight', methods=['POST'])
def design_preflight(pid):
    """Return a non-mutating experimental-design and contrast preview.

    This endpoint intentionally accepts only a file belonging to ``pid``.  It
    gives the form and embedded AI the same evidence without persisting a
    guessed metadata mapping.
    """
    project = Project.get_by_id(pid)
    if not project:
        return jsonify({'error': 'Not found'}), 404
    payload = request.get_json(silent=True) or request.form or {}
    file_path = str(payload.get('file_path', '') or '')
    module_name = str(payload.get('module_name', '') or '')
    params = payload.get('params', {}) or {}
    if not isinstance(params, dict):
        return jsonify({'error': 'params 必须是对象'}), 400
    if not file_path:
        return jsonify({'error': '缺少输入文件'}), 400
    if not _validate_project_file_path(file_path, pid):
        return jsonify({'error': '输入文件不属于当前项目或不可读取'}), 403
    try:
        from modules import MODULE_REGISTRY
        from modules.design_preflight import build_design_preflight
        from modules.io_utils import read_expression_matrix

        if module_name not in MODULE_REGISTRY:
            return jsonify({'error': f'未知分析模块: {module_name}'}), 400
        adata = read_expression_matrix(file_path)
        return jsonify(build_design_preflight(adata, module_name, params, file_path))
    except Exception as exc:
        logger.exception('Design preflight failed')
        return jsonify({'error': f'分析前检查失败: {exc}'}), 500


@api_bp.route('/projects/<pid>/deg-comparisons')
def deg_comparisons(pid):
    results_dir = Config.results_dir(pid)
    if not os.path.isdir(results_dir):
        return jsonify({'comparisons': []})

    csv_files = sorted([f for f in os.listdir(results_dir)
                        if f.startswith('bulk_deg_results') and f.endswith('.csv')
                        and 'merged' not in f and 'all_comparisons' not in f
                        and 'lrt' not in f and 'top_genes' not in f])

    # 从 ResultFile 表提取真实比较名（如 "moclel vs hmc3"）
    label_map = {}
    try:
        from models import AnalysisTask, ResultFile
        for t in AnalysisTask.get_by_project(pid):
            if t.module_name != 'bulk_deg':
                continue
            for rf in ResultFile.get_by_task(t.id):
                fname = os.path.basename(rf.file_path)
                if fname in csv_files and rf.label and 'vs' in rf.label:
                    # 从 "差异表达基因列表 (moclel vs hmc3)" 提取 "moclel vs hmc3"
                    lbl = rf.label
                    if '(' in lbl and ')' in lbl:
                        lbl = lbl.split('(')[-1].rstrip(')')
                    label_map[fname] = lbl
    except Exception:
        pass

    comparisons = []
    for f in csv_files:
        key = f.replace('bulk_deg_', '').replace('.csv', '')
        comparisons.append({'key': key, 'label': label_map.get(f, key)})

    return jsonify({'comparisons': comparisons})


@api_bp.route('/projects/<pid>/validate-filter-expression', methods=['POST'])
def validate_filter_expression(pid):
    """验证筛选表达式并返回匹配基因数"""
    from modules.expression_parser import validate, evaluate, ParseError
    import numpy as np

    data = request.get_json(silent=True) or {}
    expr = data.get('expression', '').strip()
    selected = data.get('selected_comparisons', [])

    if not expr:
        return jsonify({'valid': False, 'error': '表达式为空'})

    results_dir = Config.results_dir(pid)
    if not os.path.isdir(results_dir):
        return jsonify({'valid': False, 'error': '无 DEG 结果文件，请先运行 bulk_deg'})

    # Load comparison labels (same logic as deg_comparisons)
    csv_files = sorted([f for f in os.listdir(results_dir)
                        if f.startswith('bulk_deg_results') and f.endswith('.csv')
                        and 'merged' not in f and 'all_comparisons' not in f
                        and 'lrt' not in f and 'top_genes' not in f])
    if selected:
        csv_files = [f for f in csv_files
                     if f.replace('bulk_deg_', '').replace('.csv', '') in set(selected)]

    label_map = {}
    try:
        from models import AnalysisTask, ResultFile
        all_csv = set(f for f in os.listdir(results_dir)
                      if f.startswith('bulk_deg_results_') and f.endswith('.csv'))
        for t in AnalysisTask.get_by_project(pid):
            if t.module_name != 'bulk_deg':
                continue
            for rf in ResultFile.get_by_task(t.id):
                fname = os.path.basename(rf.file_path)
                if fname in all_csv and rf.label and 'vs' in rf.label:
                    lbl = rf.label
                    if '(' in lbl and ')' in lbl:
                        lbl = lbl.split('(')[-1].rstrip(')')
                    label_map[fname] = lbl.replace(' vs ', '-vs-').strip()
    except Exception:
        pass

    # Build matrices
    import pandas as pd
    comparisons = {}
    for f in csv_files:
        df = pd.read_csv(os.path.join(results_dir, f))
        if 'gene' not in df.columns or 'log2FC' not in df.columns:
            continue
        real_name = label_map.get(f, f.replace('bulk_deg_', '').replace('.csv', ''))
        comparisons[real_name] = df

    if len(comparisons) < 1:
        return jsonify({'valid': False, 'error': '未找到有效的比较结果'})

    all_genes = set()
    for df in comparisons.values():
        all_genes.update(df['gene'].tolist())
    all_genes = sorted(all_genes)
    comp_names = sorted(comparisons.keys())

    fc_threshold = float(data.get('fc_threshold', 2.0))
    pval_threshold = float(data.get('pval_threshold', 0.05))
    log2fc_thresh = np.log2(fc_threshold)

    padj_matrix = pd.DataFrame(1.0, index=all_genes, columns=comp_names)
    logfc_matrix = pd.DataFrame(0.0, index=all_genes, columns=comp_names)
    all_genes_idx = pd.Index(all_genes)
    for name, df in comparisons.items():
        df_idx = df.set_index('gene')
        common = all_genes_idx.intersection(df_idx.index)
        if len(common) > 0:
            logfc_matrix.loc[common, name] = df_idx.loc[common, 'log2FC']
            padj_matrix.loc[common, name] = df_idx.loc[common, 'padj']

    # Validate expression
    ast, err, warnings = validate(expr, comp_names)
    if err:
        return jsonify({'valid': False, 'error': err})

    # Build gene sets and evaluate
    gene_sets = {}
    for c in comp_names:
        sig_mask = (padj_matrix[c] < pval_threshold) & (abs(logfc_matrix[c]) >= log2fc_thresh)
        gene_sets[c] = set(logfc_matrix.index[sig_mask])

    try:
        filtered = evaluate(ast, comp_names, gene_sets, padj_matrix, logfc_matrix,
                            pval_threshold, log2fc_thresh)
    except ParseError as e:
        return jsonify({'valid': False, 'error': str(e)})

    return jsonify({
        'valid': True,
        'gene_count': len(filtered),
        'comparisons_used': comp_names,
        'warnings': warnings,
        'parse_tree': ast.to_dict(),
    })


@api_bp.route('/data-info')
def data_info():
    """返回输入文件的基本数据信息（样本数、基因数、样本名、注释列）"""
    file_path = request.args.get('file_path', '')
    if not file_path:
        return jsonify({'error': 'No file path'}), 400
    if not _validate_file_path(file_path):
        return jsonify({'error': '文件路径不在允许范围内'}), 403
    try:
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(file_path)
        qc_columns = {'total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'size_factor',
                       'total_counts_mt', 'total_counts_ribo', 'pct_counts_ribo', '_auto_group'}
        obs_cols = [c for c in adata.obs.columns if c not in qc_columns]
        return jsonify({
            'n_obs': adata.n_obs,
            'n_vars': adata.n_vars,
            'obs_columns': obs_cols,
            'sample_names': list(adata.obs.index),
            'input_measurement': infer_expression_measurement(adata, file_path),
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@api_bp.route('/obs-columns')
def obs_columns():
    """返回输入文件的 obs 列名和自动检测的分组信息"""
    file_path = request.args.get('file_path', '')
    module_name = request.args.get('module_name', '')
    if not file_path:
        return jsonify({'columns': [], 'grouping_candidates': [], 'sample_groups': {},
                        'time_candidates': [], 'sample_candidates': [],
                        'suggestions': {}})
    if not _validate_file_path(file_path):
        return jsonify({'error': '文件路径不在允许范围内'}), 403
    try:
        from modules.io_utils import (
            read_expression_matrix, infer_expression_measurement, infer_sample_group_candidates,
            rank_obs_grouping_candidates, QC_OBS_COLUMNS,
        )
        import pandas as pd
        adata = read_expression_matrix(file_path)
        qc_columns = QC_OBS_COLUMNS
        cols = [c for c in adata.obs.columns if c not in qc_columns]
        grouping_candidates = rank_obs_grouping_candidates(
            adata, module_name, purpose='groupby',
        )
        group_column_candidates = rank_obs_grouping_candidates(
            adata, module_name, purpose='group_column',
        )
        batch_candidates = rank_obs_grouping_candidates(
            adata, module_name, purpose='batch_key',
        )
        cluster_candidates = rank_obs_grouping_candidates(
            adata, module_name, purpose='cluster_key',
        )
        sample_key_candidates = rank_obs_grouping_candidates(
            adata, module_name, purpose='sample_key',
        )
        condition_candidates = rank_obs_grouping_candidates(
            adata, module_name, purpose='condition_key',
        )

        # Each parameter has a different meaning.  In particular, a technical
        # batch column must not be reused as the biological DEG groupby column.
        suggestions = {
            'groupby': grouping_candidates[0] if grouping_candidates else '',
            'group_column': group_column_candidates[0] if group_column_candidates else '',
            'color_by': grouping_candidates[0] if grouping_candidates else '',
            'batch_key': batch_candidates[0] if batch_candidates else '',
            'cluster_key': cluster_candidates[0] if cluster_candidates else '',
            'sample_key': sample_key_candidates[0] if sample_key_candidates else '',
            'condition_key': condition_candidates[0] if condition_candidates else '',
        }

        # 检测可能的时间列：除低基数数值外，也识别 D0/D3/D7、0h/12h
        # 等常见标签。前端仍要求用户确认顺序，不把普通分类列当成时间。
        from modules.sc_timecourse import _time_value
        time_candidates = []
        time_name_tokens = ('time', 'day', 'hour', 'week', 'stage', 'minute')
        technical_batch_names = {'batch', 'technical_batch', 'sequencing_batch', 'library_batch'}
        for c in cols:
            if c in qc_columns or str(c).lower() in technical_batch_names:
                continue
            col_data = adata.obs[c]
            unique_values = col_data.dropna().astype(str).unique().tolist()
            low_cardinality = 3 <= len(unique_values) <= 30
            is_numeric = pd.api.types.is_numeric_dtype(col_data)
            label_looks_temporal = bool(unique_values) and all(
                _time_value(value) is not None for value in unique_values
            )
            name_looks_temporal = any(token in str(c).lower() for token in time_name_tokens)
            if low_cardinality and (is_numeric or label_looks_temporal or name_looks_temporal):
                time_candidates.append(c)

        # 单细胞时序模块需要生物学样本 ID；它不能用普通细胞类型列替代。
        # ``batch`` and donor IDs are not automatically treated as independent
        # biological samples for temporal inference.  They remain selectable
        # manually only after the user verifies the experimental design.
        sample_priority = ('sample_id', 'sample', 'library_id', 'orig.ident')
        lower_to_original = {str(c).lower(): c for c in cols}
        sample_candidates = [lower_to_original[name] for name in sample_priority
                             if name in lower_to_original]
        if sample_candidates:
            # sample_id can be unique per row in bulk metadata and is still a
            # valid sample_key, even though it is intentionally rejected as a
            # generic grouping column.
            suggestions['sample_key'] = sample_candidates[0]

        # 无论是否存在技术 obs 列，都尝试从样本名推断分组。此前只在
        # ``cols`` 为空时执行，导致带 barcode/batch/sample 元数据的 h5ad
        # 永远得不到 _auto_group_ 候选。
        sample_groups = {}
        if adata.n_obs > 0:
            sample_names = adata.obs.index.tolist()
            candidates = infer_sample_group_candidates(sample_names)
            if candidates:
                sample_groups['auto_group'] = candidates[0]
                sample_groups['auto_group_candidates'] = candidates

        return jsonify({
            'columns': cols,
            'grouping_candidates': grouping_candidates,
            'sample_groups': sample_groups,
            'time_candidates': time_candidates,
            'sample_candidates': sample_candidates,
            'suggestions': suggestions,
            'input_measurement': infer_expression_measurement(adata, file_path),
        })
    except Exception:
        return jsonify({
            'columns': [], 'grouping_candidates': [],
            'sample_groups': {}, 'time_candidates': [], 'sample_candidates': [],
            'suggestions': {},
        })


@api_bp.route('/data-info-full')
def data_info_full():
    """返回 h5ad 文件的完整结构信息（obs/var 详情、layers、obsm、varm、obsp、uns）"""
    file_path = request.args.get('file_path', '')
    if not file_path:
        return jsonify({'error': 'No file path'}), 400
    if not _validate_file_path(file_path):
        return jsonify({'error': '文件路径不在允许范围内'}), 403
    if not file_path.endswith('.h5ad'):
        return jsonify({'error': 'Not an h5ad file'}), 400
    try:
        import anndata
        from modules.inspect_utils import inspect_adata
        adata = anndata.read_h5ad(file_path, backed='r')
        info = inspect_adata(adata, include_samples=True)
        adata.file.close()
        return jsonify(info)
    except Exception as e:
        logger.exception("API error")
        return jsonify({'error': str(e)}), 500


@api_bp.route('/presets')
def list_presets():
    project_id = request.args.get('project_id', '')
    filter_type = request.args.get('type', '')
    presets = []
    _ensure_bundled_presets()
    presets.extend(_list_presets_in_dir(_get_global_presets_dir(), 'global'))
    if project_id:
        proj_dir = _get_project_presets_dir(project_id)
        presets.extend(_list_presets_in_dir(proj_dir, 'project'))
    if filter_type:
        presets = [p for p in presets if p.get('analysis_type') == filter_type]
    result = []
    for p in presets:
        pipeline = p.get('pipeline')
        modules = pipeline.get('modules') if isinstance(pipeline, Mapping) else None
        has_pipeline = isinstance(modules, list) and bool(modules)
        result.append({
            'id': p['_id'],
            'name': p.get('name', ''),
            'description': p.get('description', ''),
            'analysis_type': p.get('analysis_type', ''),
            'scope': p['_scope'],
            # Do not expose parameter values in a list response.  The
            # one-click endpoint resolves the saved template server-side.
            'has_pipeline': has_pipeline,
            'pipeline': {'modules': list(modules)} if has_pipeline else None,
        })
    return jsonify({'presets': result})


@api_bp.route('/presets/<preset_id>')
def get_preset(preset_id):
    if '..' in preset_id or '/' in preset_id:
        return jsonify({'error': '无效的预设 ID'}), 400
    project_id = request.args.get('project_id', '')
    p = _load_accessible_preset(preset_id, project_id)
    if p:
        return jsonify({'preset': p})
    return jsonify({'error': '预设不存在'}), 404


@api_bp.route('/presets', methods=['POST'])
def save_preset():
    data = request.get_json(silent=True)
    if not isinstance(data, Mapping) or not str(data.get('name', '') or '').strip():
        return jsonify({'error': '预设名称不能为空'}), 400
    name = str(data.get('name') or '').strip()
    if len(name) > 160:
        return jsonify({'error': '预设名称不能超过 160 个字符'}), 400

    analysis_type = str(data.get('analysis_type', '') or '').strip()
    pipeline = data.get('pipeline', {})
    if not isinstance(pipeline, Mapping):
        return jsonify({'error': 'pipeline 必须是对象'}), 400
    pipeline_modules = pipeline.get('modules')
    if pipeline_modules is not None:
        if analysis_type not in {'sc', 'bulk'}:
            return jsonify({'error': '流程模板的 analysis_type 必须是 sc 或 bulk'}), 400
        if not isinstance(pipeline_modules, list) or not pipeline_modules:
            return jsonify({'error': '流程模板至少需要一个步骤'}), 400
        if any(not isinstance(module, str) or not module.strip()
               for module in pipeline_modules):
            return jsonify({'error': '流程模板包含无效步骤'}), 400
        if len(set(pipeline_modules)) != len(pipeline_modules):
            return jsonify({'error': '流程模板不能重复包含同一步骤'}), 400
        from modules import (
            BULK_MODULE_NAMES, MODULE_REGISTRY, SC_MODULE_NAMES,
            validate_pipeline_order,
        )
        allowed_modules = SC_MODULE_NAMES if analysis_type == 'sc' else BULK_MODULE_NAMES
        unknown = [module for module in pipeline_modules
                   if module not in MODULE_REGISTRY or module not in allowed_modules]
        if unknown:
            return jsonify({'error': f'流程模板包含不属于 {analysis_type} 的步骤: {unknown[0]}'}), 400
        valid_order, errors = validate_pipeline_order(pipeline_modules)
        if not valid_order:
            return jsonify({'error': '模块顺序或前置依赖不满足: ' + '；'.join(errors)}), 400

    preset_id = str(uuid.uuid4())[:8]
    preset = {
        'name': name,
        'description': data.get('description', ''),
        'analysis_type': analysis_type,
        'params': data.get('params', {}),
        'pipeline': dict(pipeline),
        'filters': data.get('filters', {}),
        'visualization': data.get('visualization', {}),
    }
    scope = data.get('scope', 'project')
    project_id = data.get('project_id', '')
    if scope == 'global':
        target_dir = _get_global_presets_dir()
    else:
        if not project_id:
            return jsonify({'error': '项目级预设需要 project_id'}), 400
        if not Project.get_by_id(project_id):
            return jsonify({'error': '项目不存在'}), 404
        target_dir = _get_project_presets_dir(project_id)
    _ensure_dir(target_dir)
    fpath = os.path.join(target_dir, f'{preset_id}.json')
    with open(fpath, 'w', encoding='utf-8') as f:
        json.dump(preset, f, ensure_ascii=False, indent=2)
    return jsonify({'id': preset_id, 'message': '预设已保存'})


@api_bp.route('/presets/<preset_id>', methods=['DELETE'])
def delete_preset(preset_id):
    if '..' in preset_id or '/' in preset_id:
        return jsonify({'error': '无效的预设 ID'}), 400
    project_id = request.args.get('project_id', '')
    if project_id:
        fpath = os.path.join(_get_project_presets_dir(project_id), f'{preset_id}.json')
        if os.path.isfile(fpath):
            os.remove(fpath)
            return jsonify({'message': '预设已删除'})
    fpath = os.path.join(_get_global_presets_dir(), f'{preset_id}.json')
    if os.path.isfile(fpath):
        os.remove(fpath)
        return jsonify({'message': '预设已删除'})
    return jsonify({'error': '预设不存在'}), 404


# ============ WES Workflow API ============

@api_bp.route('/projects/<pid>/wes/sra/upload', methods=['POST'])
def upload_wes_sra(pid):
    """Upload one local SRA archive and queue asynchronous FASTQ conversion."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    try:
        from modules.workflows.sra import store_sra_upload
        result = store_sra_upload(
            project_id=pid,
            file_storage=request.files.get('sra'),
            sample_id=request.form.get('sample_id', ''),
            patient_id=request.form.get('patient_id', ''),
            role=request.form.get('role', ''),
            capture_kit_id=request.form.get('capture_kit_id', ''),
        )
    except (OSError, ValueError) as exc:
        return jsonify({'error': str(exc)}), 422
    job = result['job']
    return jsonify({
        'message': 'SRA 已上传，FASTQ 转换任务已进入后台队列',
        'job': job,
        'source_asset': result['source_asset'],
    }), 202


@api_bp.route('/projects/<pid>/wes/sra/jobs', methods=['GET'])
def list_wes_sra_jobs(pid):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.sra import list_sra_jobs
    return jsonify({'jobs': list_sra_jobs(pid)})


@api_bp.route('/projects/<pid>/wes/sra/jobs/<job_id>', methods=['GET'])
def get_wes_sra_job(pid, job_id):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.sra import get_sra_job
    job = get_sra_job(job_id, pid)
    if not job:
        return jsonify({'error': 'SRA 转换任务不存在'}), 404
    return jsonify({'job': job})


@api_bp.route('/projects/<pid>/wes/data/upload', methods=['POST'])
def upload_wes_input_data(pid):
    """Upload one complete FASTQ/BAM/CRAM input set for the WES wizard."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    input_type = str(request.form.get('input_type', '') or '').strip().lower()
    files = {
        field: request.files.get(field)
        for field in ('fastq_1', 'fastq_2', 'bam', 'bai', 'cram', 'crai')
    }
    try:
        from modules.workflows.input_uploads import store_wes_input_upload
        result = store_wes_input_upload(
            project_id=pid,
            sample_id=request.form.get('sample_id', ''),
            patient_id=request.form.get('patient_id', ''),
            role=request.form.get('role', ''),
            input_type=input_type,
            files=files,
            capture_kit_id=request.form.get('capture_kit_id', ''),
        )
    except (OSError, ValueError) as exc:
        return jsonify({'error': str(exc)}), 422
    return jsonify({
        'message': 'WES 输入文件上传并登记成功',
        **result,
    }), 201


@api_bp.route('/projects/<pid>/wes/workflows', methods=['GET'])
def list_wes_workflows(pid):
    """List declarative WES workflows without exposing shell commands."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.registry import list_workflows
    return jsonify({'workflows': list_workflows('wes')})


@api_bp.route('/projects/<pid>/wes/manifests', methods=['GET'])
def list_wes_manifests(pid):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.storage import list_manifests
    return jsonify({'manifests': list_manifests(pid)})


@api_bp.route('/projects/<pid>/wes/references', methods=['GET'])
def list_wes_references(pid):
    """List the administrator-managed reference assets visible to WES."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.references import list_reference_assets
    return jsonify({'references': list_reference_assets(
        assembly=request.args.get('assembly', '').strip(),
        bundle_version=request.args.get('bundle_version', '').strip(),
        status=request.args.get('status', '').strip(),
    )})


@api_bp.route('/projects/<pid>/wes/capture-kits', methods=['GET'])
def list_wes_capture_kits(pid):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.references import list_capture_kit_profiles
    return jsonify({'capture_kits': list_capture_kit_profiles(
        assembly=request.args.get('assembly', '').strip(),
        include_retired=request.args.get('include_retired', '').lower() == 'true',
    )})


@api_bp.route('/projects/<pid>/wes/references/readiness', methods=['GET'])
def wes_reference_readiness(pid):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    workflow_key = request.args.get('workflow_key', '').strip()
    reference_bundle_id = request.args.get('reference_bundle_id', '').strip()
    capture_bed_id = request.args.get('capture_bed_id', '').strip()
    capture_bed_path = request.args.get('capture_bed_path', '').strip()
    if not workflow_key or not reference_bundle_id:
        return jsonify({'error': '需要 workflow_key 和 reference_bundle_id'}), 400
    from modules.workflows.references import launch_reference_readiness
    readiness = launch_reference_readiness(
        reference_bundle_id=reference_bundle_id,
        capture_bed_id=capture_bed_id,
        capture_bed_path=capture_bed_path,
        workflow_key=workflow_key,
    )
    return jsonify({'readiness': readiness}), 200 if readiness['valid'] else 422


@api_bp.route('/projects/<pid>/wes/runs', methods=['GET'])
def list_wes_runs(pid):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.runs import list_workflow_runs
    from modules.workflows.nextflow import NextflowExecutor
    runs = list_workflow_runs(pid)
    for item in runs:
        if item.get('status') in {'running', 'cancel_requested'}:
            NextflowExecutor.poll(item['id'], pid)
    return jsonify({'runs': list_workflow_runs(pid)})


@api_bp.route('/projects/<pid>/wes/runs/<run_id>', methods=['GET'])
def get_wes_run(pid, run_id):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.runs import get_workflow_run
    from modules.workflows.nextflow import NextflowExecutor
    run = get_workflow_run(run_id, pid)
    if not run:
        return jsonify({'error': 'WES run 不存在'}), 404
    if run.get('status') in {'running', 'cancel_requested'}:
        run = NextflowExecutor.poll(run_id, pid)
    return jsonify({'run': run})


@api_bp.route('/projects/<pid>/wes/runs/<run_id>/artifacts', methods=['GET'])
def wes_run_artifacts(pid, run_id):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.runs import get_workflow_run
    from modules.workflows.artifacts import collect_workflow_artifacts, list_workflow_artifacts
    run = get_workflow_run(run_id, pid)
    if not run:
        return jsonify({'error': 'WES run 不存在'}), 404
    collection = None
    if run.get('status') == 'completed' and (
            request.args.get('refresh', '').lower() == 'true' or
            not list_workflow_artifacts(run_id, pid)):
        collection = collect_workflow_artifacts(run_id, pid, strict=True)
        run = get_workflow_run(run_id, pid)
    return jsonify({
        'run': run,
        'artifacts': list_workflow_artifacts(run_id, pid),
        'collection': collection,
    })


def _controlled_stage_samplesheet(parent_run, stage):
    """Return a validated Sarek CSV produced by the immediately prior stage.

    Only a completed run's own results directory is trusted.  Every path in
    the generated sheet must be an absolute, non-symlink file under that same
    directory; this prevents a hand-edited CSV from becoming a general file
    reference mechanism for downstream Nextflow runs.
    """
    run_root = os.path.realpath(str(parent_run.get('run_dir') or ''))
    results_dir = os.path.realpath(os.path.join(run_root, 'results'))
    if not run_root or not results_dir.startswith(run_root + os.sep):
        raise ValueError('上游 run 结果目录无效')
    expected_name = str(stage.get('input_kind') or '')
    filename = {
        'recalibrated_csv': 'recalibrated.csv',
        'variantcalled_csv': 'variantcalled.csv',
    }.get(expected_name)
    if not filename:
        raise ValueError('该 WES 阶段不接受上游 samplesheet')
    candidates = (
        os.path.join(results_dir, 'preprocessing', 'csv', filename),
        os.path.join(results_dir, 'csv', filename),
    )
    source = next((path for path in candidates if os.path.isfile(path) and not os.path.islink(path)), '')
    real_source = os.path.realpath(source) if source else ''
    if not real_source or not real_source.startswith(results_dir + os.sep):
        raise ValueError(f'上游阶段未产生可用的 {filename}')
    try:
        with open(real_source, encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
            fields = set(reader.fieldnames or [])
    except (OSError, csv.Error) as exc:
        raise ValueError(f'无法读取上游 {filename}: {exc}') from exc
    if not rows:
        raise ValueError(f'上游 {filename} 没有样本行')
    required = {'patient', 'sample'}
    if expected_name == 'recalibrated_csv':
        if not required.issubset(fields) or not ({'cram', 'crai'}.issubset(fields) or {'bam', 'bai'}.issubset(fields)):
            raise ValueError('上游 recalibrated.csv 列不符合 Sarek 接口')
        path_fields = ('cram', 'crai') if {'cram', 'crai'}.issubset(fields) else ('bam', 'bai')
    else:
        if not required.issubset(fields) or 'vcf' not in fields:
            raise ValueError('上游 variantcalled.csv 列不符合 Sarek 接口')
        path_fields = ('vcf',)
    for row in rows:
        for field in path_fields:
            path = str(row.get(field) or '').strip()
            if not path or not os.path.isabs(path) or os.path.islink(path):
                raise ValueError(f'上游 {filename} 包含无效的 {field} 路径')
            real_path = os.path.realpath(path)
            if not os.path.isfile(real_path) or not real_path.startswith(results_dir + os.sep):
                raise ValueError(f'上游 {filename} 的 {field} 不属于已完成阶段结果')
    return real_source


def _prepare_wes_stage_run(pid, manifest, stage, *, parent_run=None, stage_options=None):
    """Create one reviewable stage run; this never starts Nextflow."""
    from modules.workflows.registry import get_workflow
    from modules.workflows.sarek import write_launch_bundle
    from modules.workflows.nextflow import NextflowExecutor, WorkflowNotConfiguredError
    from modules.workflows.runs import create_workflow_run

    workflow = get_workflow(stage['run_workflow_key'])
    if not workflow:
        raise ValueError(f"未登记 WES 阶段 workflow: {stage['run_workflow_key']}")
    run_id = 'wesrun_' + uuid.uuid4().hex[:16]
    run_root = os.path.join(Config.project_dir(pid), 'workflow_runs', run_id)
    work_dir = os.path.join(run_root, 'work')
    results_dir = os.path.join(run_root, 'results')
    launch_dir = os.path.join(run_root, 'launch')
    logs_dir = os.path.join(run_root, 'logs')
    for directory in (work_dir, results_dir, launch_dir, logs_dir):
        os.makedirs(directory, exist_ok=True)
    samplesheet_path = os.path.join(launch_dir, 'samplesheet.csv')
    source_sheet = _controlled_stage_samplesheet(parent_run, stage) if parent_run else ''
    manifest_body = dict(manifest.get('manifest') or {})
    if stage_options is not None:
        if not isinstance(stage_options, dict):
            raise ValueError('阶段参数必须是对象')
        allowed = set(stage.get('parameter_keys') or ())
        unexpected = sorted(set(stage_options) - allowed)
        if unexpected:
            raise ValueError('该 WES 阶段不支持参数: ' + ', '.join(unexpected))
        selected_options = dict(manifest_body.get('pipeline_options') or {})
        selected_options.update(stage_options)
        manifest_body['pipeline_options'] = selected_options
    manifest_samples = manifest_body.get('samples') or []
    raw_input_types = {
        str(sample.get('input_type') or '').strip().lower()
        for sample in manifest_samples
        if str(sample.get('input_type') or '').strip()
    }
    if not parent_run and len(raw_input_types) != 1:
        raise ValueError('WES manifest 必须只包含一种已校验的 input_type')
    raw_input_type = next(iter(raw_input_types), '')
    executor = NextflowExecutor()
    try:
        launch = executor.prepare(
            workflow, pid, run_id, manifest['id'], work_dir, results_dir,
            samplesheet_path=samplesheet_path,
            profile=Config.WES_NEXTFLOW_PROFILE,
            stdout_path=os.path.join(logs_dir, 'stdout.log'),
            stderr_path=os.path.join(logs_dir, 'stderr.log'),
            intervals_path=manifest_body.get('capture_bed_path', ''),
            input_type=raw_input_type if not parent_run else str(stage['input_kind']),
            reference_bundle_id=manifest_body.get('reference_bundle_id', ''),
            capture_bed_id=manifest_body.get('capture_bed_id', ''),
            pipeline_options=manifest_body.get('pipeline_options') or {},
            stage_key=stage['key'],
            analysis_workflow_key=stage['analysis_workflow_key'],
        )
    except (ValueError, WorkflowNotConfiguredError) as exc:
        raise ValueError(f'WES 阶段固定参数无效: {exc}') from exc
    try:
        bundle = write_launch_bundle(
            run_root, workflow=workflow.to_dict(), run_id=run_id, project_id=pid,
            # Preserve the originally registered manifest ID and checksum in
            # provenance, while storing the stage-specific parameter snapshot
            # that was actually reviewed for this independent run.
            manifest_record={**manifest, 'manifest': manifest_body},
            launch=launch.to_dict(), parameters=launch.parameters,
            profile=Config.WES_NEXTFLOW_PROFILE, samplesheet_source=source_sheet,
        )
    except (OSError, ValueError) as exc:
        raise ValueError(f'WES 阶段 launch bundle 生成失败: {exc}') from exc
    launch_payload = launch.to_dict()
    launch_payload['bundle'] = bundle
    launch_payload['stage'] = stage
    if parent_run:
        launch_payload['parent_run_id'] = parent_run['id']
    record = create_workflow_run(
        pid, workflow.key, manifest['id'], run_id=run_id, status='prepared',
        launch=launch_payload, executor=workflow.executor, workflow_release=workflow.release,
        profile=Config.WES_NEXTFLOW_PROFILE, run_dir=run_root,
        stdout_path=launch.stdout_path, stderr_path=launch.stderr_path,
        provenance={
            'workflow_release': workflow.release,
            'executor': workflow.executor,
            'profile': Config.WES_NEXTFLOW_PROFILE,
            'manifest_checksum': manifest.get('checksum', ''),
            'bundle_checksums': bundle.get('checksums', {}),
            'stage_key': stage['key'],
            'parent_run_id': parent_run['id'] if parent_run else '',
        },
    )
    return record, launch_payload


@api_bp.route('/projects/<pid>/wes/stages/prepare', methods=['POST'])
def prepare_first_wes_stage(pid):
    """Prepare the first explicit preprocessing/mapping stage."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    payload = request.get_json(silent=True) or {}
    workflow_key = str(payload.get('workflow_key', '') or '').strip()
    manifest_id = str(payload.get('manifest_id', '') or '').strip()
    from modules.workflows.stages import first_stage, get_stage
    from modules.workflows.storage import get_manifest
    manifest = get_manifest(manifest_id, pid) if manifest_id else None
    if not manifest:
        return jsonify({'error': '需要有效的胚系/体细胞 workflow_key 和 manifest_id'}), 422
    input_types = {
        str(sample.get('input_type') or '').strip().lower()
        for sample in (manifest.get('manifest') or {}).get('samples') or []
        if str(sample.get('input_type') or '').strip()
    }
    if len(input_types) != 1:
        return jsonify({'error': 'WES manifest 必须只包含一种已校验的 input_type'}), 422
    input_type = next(iter(input_types))
    # FASTQ requires preprocessing/mapping.  Existing BAM/CRAM submissions
    # retain their supported path and enter at caller stage, rather than being
    # forced through an invalid FASTQ-only mapping stage.
    stage = first_stage(workflow_key) if input_type == 'fastq' else get_stage('variant_calling', workflow_key)
    if input_type not in {'fastq', 'bam', 'cram'} or not stage:
        return jsonify({'error': '仅支持已校验的 FASTQ、BAM 或 CRAM WES manifest'}), 422
    try:
        record, launch = _prepare_wes_stage_run(pid, manifest, stage)
    except ValueError as exc:
        return jsonify({'error': str(exc)}), 422
    return jsonify({
        'run': record, 'launch': launch, 'stage': stage,
        'message': f"已准备“{stage['title']}”的 LaunchSpec；尚未启动 Nextflow。",
    }), 201


@api_bp.route('/projects/<pid>/wes/runs/<run_id>/next-stage/prepare', methods=['POST'])
def prepare_next_wes_stage(pid, run_id):
    """Prepare exactly the next stage from a completed, controlled Sarek output."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.runs import get_workflow_run
    from modules.workflows.storage import get_manifest
    from modules.workflows.stages import next_stage, stage_from_launch
    parent = get_workflow_run(run_id, pid)
    if not parent:
        return jsonify({'error': 'WES run 不存在'}), 404
    if parent.get('status') != 'completed':
        return jsonify({'error': '仅已完成的 WES 阶段可以进入下一步', 'run': parent}), 409
    current_stage = stage_from_launch(parent.get('launch') or {})
    if not current_stage:
        return jsonify({'error': '该 run 不是可串联的分步 WES run'}), 422
    stage = next_stage(current_stage['key'], current_stage['analysis_workflow_key'])
    if not stage:
        return jsonify({'error': '当前 WES 阶段已经是最终步骤'}), 409
    manifest = get_manifest(parent.get('manifest_id', ''), pid)
    if not manifest:
        return jsonify({'error': '上游 run 的 manifest 不存在'}), 422
    payload = request.get_json(silent=True) or {}
    try:
        record, launch = _prepare_wes_stage_run(
            pid, manifest, stage, parent_run=parent,
            stage_options=payload.get('pipeline_options') if 'pipeline_options' in payload else None,
        )
    except ValueError as exc:
        return jsonify({'error': str(exc)}), 422
    return jsonify({
        'run': record, 'launch': launch, 'stage': stage,
        'message': f"已准备“{stage['title']}”的 LaunchSpec；尚未启动 Nextflow。",
    }), 201


@api_bp.route('/projects/<pid>/wes/runs/prepare', methods=['POST'])
def prepare_wes_run(pid):
    """Create a reproducible, reviewable launch bundle without starting Nextflow."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    payload = request.get_json(silent=True) or {}
    workflow_key = str(payload.get('workflow_key', '') or '').strip()
    manifest_id = str(payload.get('manifest_id', '') or '').strip()
    if not workflow_key or not manifest_id:
        return jsonify({'error': '需要 workflow_key 和 manifest_id'}), 400

    from modules.workflows.registry import get_workflow
    from modules.workflows.storage import get_manifest
    from modules.workflows.sarek import write_launch_bundle
    from modules.workflows.nextflow import NextflowExecutor, WorkflowNotConfiguredError
    from modules.workflows.runs import create_workflow_run

    workflow = get_workflow(workflow_key)
    if not workflow:
        return jsonify({'error': f'未知 WES workflow: {workflow_key}'}), 422
    manifest = get_manifest(manifest_id, pid)
    if not manifest:
        return jsonify({'error': 'manifest_id 不存在或不属于当前项目'}), 422

    run_id = 'wesrun_' + uuid.uuid4().hex[:16]
    run_root = os.path.join(Config.project_dir(pid), 'workflow_runs', run_id)
    work_dir = os.path.join(run_root, 'work')
    results_dir = os.path.join(run_root, 'results')
    launch_dir = os.path.join(run_root, 'launch')
    logs_dir = os.path.join(run_root, 'logs')
    for directory in (work_dir, results_dir, launch_dir, logs_dir):
        os.makedirs(directory, exist_ok=True)
    samplesheet_path = os.path.join(launch_dir, 'samplesheet.csv')
    stdout_path = os.path.join(logs_dir, 'stdout.log')
    stderr_path = os.path.join(logs_dir, 'stderr.log')
    manifest_samples = (manifest.get('manifest') or {}).get('samples') or []
    input_types = {
        str(sample.get('input_type') or '').strip().lower()
        for sample in manifest_samples
        if str(sample.get('input_type') or '').strip()
    }
    if len(input_types) != 1:
        return jsonify({'error': 'WES manifest 必须只包含一种已校验的 input_type'}), 422
    input_type = next(iter(input_types))
    executor = NextflowExecutor()
    try:
        launch = executor.prepare(
            workflow, pid, run_id, manifest_id, work_dir, results_dir,
            samplesheet_path=samplesheet_path,
            profile=Config.WES_NEXTFLOW_PROFILE,
            stdout_path=stdout_path, stderr_path=stderr_path,
            intervals_path=(manifest.get('manifest') or {}).get('capture_bed_path', ''),
            input_type=input_type,
            reference_bundle_id=(manifest.get('manifest') or {}).get('reference_bundle_id', ''),
            capture_bed_id=(manifest.get('manifest') or {}).get('capture_bed_id', ''),
            pipeline_options=(manifest.get('manifest') or {}).get('pipeline_options') or {},
        )
    except (ValueError, WorkflowNotConfiguredError) as exc:
        return jsonify({'error': f'WES 固定运行参数无效: {exc}'}), 422
    try:
        bundle = write_launch_bundle(
            run_root, workflow=workflow.to_dict(), run_id=run_id, project_id=pid,
            manifest_record=manifest, launch=launch.to_dict(), parameters=launch.parameters,
            profile=Config.WES_NEXTFLOW_PROFILE,
        )
    except (OSError, ValueError) as exc:
        return jsonify({'error': f'WES launch bundle 生成失败: {exc}'}), 422
    launch_payload = launch.to_dict()
    launch_payload['bundle'] = bundle
    record = create_workflow_run(
        pid, workflow_key, manifest_id, run_id=run_id, status='prepared',
        launch=launch_payload,
        executor=workflow.executor, workflow_release=workflow.release,
        profile=Config.WES_NEXTFLOW_PROFILE, run_dir=run_root,
        stdout_path=stdout_path, stderr_path=stderr_path,
        provenance={'workflow_release': workflow.release, 'executor': workflow.executor,
                    'profile': Config.WES_NEXTFLOW_PROFILE,
                    'manifest_checksum': manifest.get('checksum', ''),
                    'bundle_checksums': bundle.get('checksums', {})},
    )
    return jsonify({
        'run': record,
        'launch': launch_payload,
        'message': '已生成可审阅的 WES LaunchSpec 和 samplesheet；尚未启动 Nextflow。',
    }), 201


@api_bp.route('/projects/<pid>/wes/runs/<run_id>/launch', methods=['POST'])
def launch_wes_run(pid, run_id):
    """Start one prepared run after an explicit human confirmation."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    payload = request.get_json(silent=True) or {}
    if payload.get('confirm') is not True:
        return jsonify({'error': '启动 WES run 需要 confirm=true'}), 400
    from modules.workflows.runs import get_workflow_run
    from modules.workflows.registry import get_workflow
    from modules.workflows.nextflow import NextflowExecutor, WorkflowNotConfiguredError
    run = get_workflow_run(run_id, pid)
    if not run:
        return jsonify({'error': 'WES run 不存在'}), 404
    if run.get('status') != 'prepared':
        return jsonify({'error': f"当前状态不可启动: {run.get('status')}", 'run': run}), 409
    workflow = get_workflow(run.get('workflow_key'))
    launch = run.get('launch') or {}
    if not workflow or not launch:
        return jsonify({'error': 'WES run 缺少可执行的 launch bundle'}), 422
    executor = NextflowExecutor()
    try:
        spec = executor.prepare(
            workflow, pid, run_id, run['manifest_id'], launch.get('work_dir', ''),
            launch.get('results_dir', ''), samplesheet_path=launch.get('samplesheet_path', ''),
            profile=launch.get('profile') or Config.WES_NEXTFLOW_PROFILE,
            stdout_path=launch.get('stdout_path', ''), stderr_path=launch.get('stderr_path', ''),
            intervals_path=launch.get('intervals_path', ''),
            params_file_path=launch.get('params_file_path', ''),
            input_type=launch.get('input_type', ''),
            reference_bundle_id=launch.get('reference_bundle_id', ''),
            capture_bed_id=launch.get('capture_bed_id', ''),
            pipeline_options=launch.get('pipeline_options') or {},
            stage_key=launch.get('stage_key', ''),
            analysis_workflow_key=launch.get('analysis_workflow_key', ''),
        )
        runtime = executor.launch(spec, project_id=pid)
    except (ValueError, WorkflowNotConfiguredError) as exc:
        return jsonify({'error': str(exc), 'run': run}), 409
    return jsonify({'run': runtime.get('run'), 'runtime': {
        key: value for key, value in runtime.items() if key != 'process'
    }, 'message': 'Nextflow 已启动；可通过 run 查询状态和日志。'}), 202


@api_bp.route('/projects/<pid>/wes/runs/<run_id>/cancel', methods=['POST'])
def cancel_wes_run(pid, run_id):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    payload = request.get_json(silent=True) or {}
    if payload.get('confirm') is not True:
        return jsonify({'error': '取消 WES run 需要 confirm=true'}), 400
    from modules.workflows.nextflow import NextflowExecutor
    run = NextflowExecutor.cancel(run_id, pid)
    if not run:
        return jsonify({'error': 'WES run 不存在'}), 404
    return jsonify({'run': run, 'message': '已请求停止 WES 进程；下一次状态查询会确认最终状态。'}), 202


@api_bp.route('/projects/<pid>/wes/runs/<run_id>/resume', methods=['POST'])
def resume_wes_run(pid, run_id):
    """Resume a failed/interrupted/cancelled local run with Nextflow -resume."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    payload = request.get_json(silent=True) or {}
    if payload.get('confirm') is not True:
        return jsonify({'error': '恢复 WES run 需要 confirm=true'}), 400
    from modules.workflows.runs import get_workflow_run
    from modules.workflows.registry import get_workflow
    from modules.workflows.nextflow import NextflowExecutor, WorkflowNotConfiguredError
    run = get_workflow_run(run_id, pid)
    if not run:
        return jsonify({'error': 'WES run 不存在'}), 404
    if run.get('status') not in {'failed', 'interrupted', 'cancelled'}:
        return jsonify({'error': f"当前状态不可恢复: {run.get('status')}", 'run': run}), 409
    workflow = get_workflow(run.get('workflow_key'))
    launch = run.get('launch') or {}
    if not workflow or not launch:
        return jsonify({'error': 'WES run 缺少可恢复的 launch bundle'}), 422
    executor = NextflowExecutor()
    try:
        spec = executor.prepare(
            workflow, pid, run_id, run['manifest_id'], launch.get('work_dir', ''),
            launch.get('results_dir', ''), samplesheet_path=launch.get('samplesheet_path', ''),
            profile=launch.get('profile') or Config.WES_NEXTFLOW_PROFILE,
            stdout_path=launch.get('stdout_path', ''), stderr_path=launch.get('stderr_path', ''),
            intervals_path=launch.get('intervals_path', ''),
            params_file_path=launch.get('params_file_path', ''),
            input_type=launch.get('input_type', ''),
            resume=True,
            reference_bundle_id=launch.get('reference_bundle_id', ''),
            capture_bed_id=launch.get('capture_bed_id', ''),
            pipeline_options=launch.get('pipeline_options') or {},
            stage_key=launch.get('stage_key', ''),
            analysis_workflow_key=launch.get('analysis_workflow_key', ''),
        )
        runtime = executor.launch(
            spec, project_id=pid,
            expected_statuses=(run.get('status'),),
        )
    except (ValueError, WorkflowNotConfiguredError) as exc:
        return jsonify({'error': str(exc), 'run': run}), 409
    return jsonify({'run': runtime.get('run'), 'runtime': {
        key: value for key, value in runtime.items() if key != 'process'
    }, 'message': '已使用 Nextflow -resume 启动；可继续查询状态和日志。'}), 202


@api_bp.route('/projects/<pid>/wes/runs/<run_id>/logs', methods=['GET'])
def wes_run_logs(pid, run_id):
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    from modules.workflows.runs import get_workflow_run
    run = get_workflow_run(run_id, pid)
    if not run:
        return jsonify({'error': 'WES run 不存在'}), 404
    stream = request.args.get('stream', 'stdout').strip().lower()
    if stream not in {'stdout', 'stderr'}:
        return jsonify({'error': 'stream 必须为 stdout 或 stderr'}), 400
    try:
        tail_bytes = max(1, min(int(request.args.get('tail_bytes', '65536')), 1024 * 1024))
    except ValueError:
        return jsonify({'error': 'tail_bytes 必须是整数'}), 400
    path = run.get('stdout_path') if stream == 'stdout' else run.get('stderr_path')
    run_dir = os.path.realpath(str(run.get('run_dir') or ''))
    try:
        real_log_path = os.path.realpath(str(path or ''))
    except (OSError, ValueError):
        return jsonify({'error': '日志路径无效'}), 400
    if not run_dir or not real_log_path.startswith(run_dir + os.sep):
        return jsonify({'error': '日志路径不在当前 WES run 目录内'}), 403
    path = real_log_path
    if not path or not os.path.isfile(path):
        return jsonify({'run_id': run_id, 'stream': stream, 'text': ''})
    try:
        with open(path, 'rb') as handle:
            handle.seek(0, os.SEEK_END)
            handle.seek(max(0, handle.tell() - tail_bytes))
            content = handle.read().decode('utf-8', errors='replace')
    except OSError as exc:
        return jsonify({'error': f'日志读取失败: {exc}'}), 500
    return jsonify({'run_id': run_id, 'stream': stream, 'text': content})


@api_bp.route('/projects/<pid>/wes/preflight', methods=['POST'])
def wes_preflight(pid):
    """Validate and version a WES manifest; does not launch a workflow."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    payload = request.get_json(silent=True) or {}
    workflow_key = str(payload.get('workflow_key', '') or '').strip()
    if not workflow_key:
        return jsonify({'error': '缺少 workflow_key'}), 400

    manifest = payload.get('manifest')
    if manifest is None:
        manifest_path = str(payload.get('manifest_path', '') or '').strip()
        if not manifest_path:
            return jsonify({'error': '需要 manifest 或 manifest_path'}), 400
        if not os.path.isabs(manifest_path):
            manifest_path = os.path.join(Config.project_dir(pid), manifest_path)
        in_project = _validate_project_file_path(manifest_path, pid)
        if not in_project:
            try:
                manifest_path = Config.validate_wes_source_path(manifest_path)
            except ValueError:
                return jsonify({'error': 'manifest_path 必须位于项目目录或管理员配置的 WES_SOURCE_ROOTS'}), 400
        try:
            from modules.workflows.preflight import load_manifest_file
            manifest = load_manifest_file(manifest_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            return jsonify({'error': f'manifest 读取失败: {exc}'}), 400

    # Lab users select a stable capture-kit profile; the server resolves its
    # immutable path.  Legacy/ad-hoc manifests still have to provide a path
    # explicitly and remain blocked at production launch unless catalogued.
    capture_profile = None
    if isinstance(manifest, dict):
        from modules.workflows.guided import normalize_guided_options
        try:
            manifest = dict(manifest)
            manifest['pipeline_options'] = normalize_guided_options(
                manifest.get('pipeline_options'), workflow_key
            )
        except ValueError as exc:
            return jsonify({'valid': False, 'errors': [str(exc)], 'warnings': []}), 422
        capture_id = str(manifest.get('capture_bed_id', '') or '').strip()
        if capture_id:
            from modules.workflows.references import get_capture_kit_profile
            capture_profile = get_capture_kit_profile(capture_id)
            calling = ((capture_profile or {}).get('assets') or {}).get('calling')
            if calling:
                supplied = str(manifest.get('capture_bed_path', '') or '').strip()
                catalog_path = calling.get('file_path', '')
                if supplied and os.path.realpath(supplied) != os.path.realpath(catalog_path):
                    return jsonify({'valid': False, 'errors': [
                        'capture_bed_path 与所选 capture kit profile 不一致'
                    ], 'warnings': []}), 422
                manifest = dict(manifest)
                manifest['capture_bed_path'] = catalog_path

    from modules.workflows.preflight import validate_manifest
    result = validate_manifest(
        manifest,
        project_dir=Config.project_dir(pid),
        workflow_key=workflow_key,
        require_files=bool(payload.get('check_files', True)),
        source_roots=Config.wes_source_roots(),
        content_checks=bool(payload.get('check_content', False)),
        full_fastq_integrity=bool(payload.get('full_fastq_integrity', False)),
    )
    # Reject a mistyped display label before it becomes an immutable manifest.
    # The full checksum/path readiness check remains at the actual launch gate.
    if result.valid and Config.WES_EXECUTOR_ENABLED:
        from modules.workflows.references import list_reference_assets
        bundle_id = str((result.normalized_manifest or {}).get('reference_bundle_id', '') or '')
        if not list_reference_assets(bundle_version=bundle_id, status='validated'):
            result.valid = False
            result.errors.append(f'reference_bundle_id 未匹配已认证的 bundle: {bundle_id}')
    if capture_profile and not capture_profile.get('production_allowed'):
        result.warnings.append(
            f"capture kit profile {capture_profile['capture_kit_id']} 仅为 "
            f"{capture_profile.get('status')}，可以预检但不能启动真实 run"
        )
    response = result.to_dict()
    if not result.valid:
        return jsonify(response), 422

    from modules.workflows.storage import register_manifest
    registered = register_manifest(pid, result.normalized_manifest, status='validated',
                                  validation=response)
    response['manifest_id'] = registered['id']
    response['manifest_version'] = registered['version']
    response['manifest_checksum'] = registered['checksum']
    response['message'] = 'WES manifest 已通过预检查并登记；尚未启动外部工作流。'
    return jsonify(response), 201


# ============ Pipeline Run API ============

def _build_pipeline_submission(pid, data):
    """Validate one constrained SC/Bulk pipeline submission.

    This is shared by the manual API and saved-template launcher.  A template
    may choose only registered module names and schema-recognised parameters;
    it cannot become an alternate route for arbitrary worker code or file
    paths.
    """
    from modules import MODULE_REGISTRY, SC_MODULE_NAMES, BULK_MODULE_NAMES, validate_pipeline_order
    from modules.schemas import PARAM_SCHEMAS, filter_active_params

    if not isinstance(data, Mapping):
        raise ValueError('请求体必须是对象')
    name = str(data.get('name', '') or '').strip()
    analysis_type = str(data.get('analysis_type', '') or '').strip()
    modules = data.get('modules', [])
    params = data.get('params', {})
    input_path = str(data.get('input_path', '') or '').strip()

    if not name:
        raise ValueError('缺少流程名称')
    if len(name) > 160:
        raise ValueError('流程名称不能超过 160 个字符')
    if analysis_type not in {'sc', 'bulk'}:
        raise ValueError('analysis_type 必须是 sc 或 bulk')
    if not isinstance(modules, list) or not modules:
        raise ValueError('模块列表为空')
    if any(not isinstance(module, str) or not module.strip() for module in modules):
        raise ValueError('模块列表包含无效项')
    if len(set(modules)) != len(modules):
        raise ValueError('流程模板不能重复包含同一模块')
    if not isinstance(params, Mapping):
        raise ValueError('流程参数必须是对象')
    if not input_path:
        raise ValueError('缺少输入文件路径')

    module_set = SC_MODULE_NAMES if analysis_type == 'sc' else BULK_MODULE_NAMES
    for module in modules:
        if module not in MODULE_REGISTRY:
            raise ValueError(f'未知模块: {module}')
        if module not in module_set:
            raise ValueError(f'模块 {module} 不属于 {analysis_type} 类型')
    is_valid, errors = validate_pipeline_order(modules)
    if not is_valid:
        raise ValueError('模块顺序不满足依赖约束: ' + '；'.join(errors))

    project_dir = os.path.abspath(Config.project_dir(pid))
    abs_input = os.path.abspath(input_path)
    valid_path, error = Config._validate_path(abs_input, pid)
    if not valid_path or not abs_input.startswith(project_dir + os.sep):
        raise ValueError(error or '输入文件不在项目目录内')
    if os.path.islink(abs_input):
        raise ValueError('输入文件不能是符号链接')
    if not os.path.isfile(abs_input):
        raise ValueError('输入文件不存在')

    # Templates saved by older page versions contain a flat object.  Apply
    # those values to every matching module *after* default construction;
    # applying them only when a default is absent silently ignored templates.
    params_by_module = {}
    module_keys = set(MODULE_REGISTRY)
    flat_params = {
        key: value for key, value in params.items()
        if key not in module_keys
    }
    for module in modules:
        schema = PARAM_SCHEMAS.get(module, [])
        schema_keys = {field['key'] for field in schema}
        default_params = {}
        for field in schema:
            default = field.get('default')
            if field.get('type') == 'select' and field.get('options'):
                options = field['options']
                if default not in options:
                    default = options[0]
            default_params[field['key']] = default
        module_params = params.get(module)
        if module_params is not None:
            if not isinstance(module_params, Mapping):
                raise ValueError(f'{module} 的参数必须是对象')
            default_params.update({
                key: value for key, value in module_params.items()
                if key in schema_keys
            })
        default_params.update({
            key: value for key, value in flat_params.items()
            if key in schema_keys
        })
        params_by_module[module] = filter_active_params(schema, default_params)
    return {
        'name': name,
        'analysis_type': analysis_type,
        'modules': modules,
        'params_by_module': params_by_module,
        'input_path': abs_input,
        'project_dir': project_dir,
    }


def _start_pipeline_submission(pid, data):
    """Persist and submit a validated pipeline to the existing background worker."""
    from models import PipelineRun
    from worker import submit_pipeline_run

    submission = _build_pipeline_submission(pid, data)
    pipeline_run = PipelineRun(
        project_id=pid,
        name=submission['name'],
        analysis_type=submission['analysis_type'],
        input_path=submission['input_path'],
        modules_json=json.dumps(submission['modules']),
        params_json=json.dumps(submission['params_by_module'], ensure_ascii=False),
    )
    pipeline_run.save()
    submit_pipeline_run(
        run_id=pipeline_run.id,
        project_id=pid,
        modules=submission['modules'],
        params_by_module=submission['params_by_module'],
        project_dir=submission['project_dir'],
        input_path=submission['input_path'],
    )
    return pipeline_run, submission

@api_bp.route('/projects/<pid>/pipeline-runs', methods=['POST'])
def create_pipeline_run(pid):
    """创建并启动 pipeline run。"""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    data = request.get_json(silent=True)
    if not data:
        return jsonify({'error': '请求体为空'}), 400
    try:
        pipeline_run, _ = _start_pipeline_submission(pid, data)
    except ValueError as exc:
        return jsonify({'error': str(exc)}), 400
    return jsonify({'id': pipeline_run.id, 'message': '流程已启动'}), 201


@api_bp.route('/projects/<pid>/pipeline-templates/<preset_id>/run', methods=['POST'])
def run_pipeline_template(pid, preset_id):
    """Start a saved SC/Bulk template without trusting browser-supplied steps.

    The user still explicitly clicks the launch control, while the server reads
    the modules and parameters from the saved template and applies the same
    path, dependency and schema validation as the normal pipeline endpoint.
    """
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目不存在'}), 404
    template = _load_accessible_preset(preset_id, pid)
    if not template:
        return jsonify({'error': '流程模板不存在'}), 404
    pipeline = template.get('pipeline')
    if not isinstance(pipeline, Mapping) or not pipeline.get('modules'):
        return jsonify({'error': '该预设不是可运行的流程模板'}), 422
    request_data = request.get_json(silent=True) or {}
    payload = {
        'name': str(template.get('name') or '').strip(),
        'analysis_type': template.get('analysis_type'),
        'modules': pipeline.get('modules'),
        'params': template.get('params') or {},
        'input_path': request_data.get('input_path', ''),
    }
    try:
        pipeline_run, submission = _start_pipeline_submission(pid, payload)
    except ValueError as exc:
        return jsonify({'error': str(exc)}), 422
    return jsonify({
        'id': pipeline_run.id,
        'template_id': preset_id,
        'template_scope': template.get('_scope'),
        'modules': submission['modules'],
        'message': '流程模板已在后台启动。',
    }), 201


@api_bp.route('/projects/<pid>/pipeline-runs')
def list_pipeline_runs(pid):
    """列出项目的所有 pipeline runs。"""
    from models import PipelineRun

    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'error': '项目不存在'}), 404

    runs = PipelineRun.get_by_project(pid)
    return jsonify([r.to_dict() for r in runs])


@api_bp.route('/pipeline-runs/<run_id>/status')
def pipeline_run_status(run_id):
    """查询 pipeline run 状态。"""
    from models import PipelineRun

    run = PipelineRun.get_by_id(run_id)
    if not run:
        return jsonify({'error': 'Pipeline run 不存在'}), 404

    return jsonify(run.to_dict())


@api_bp.route('/projects/<pid>/pipeline-runs/<run_id>/resume', methods=['POST'])
def resume_pipeline_run(pid, run_id):
    """Create a new PipelineRun starting from the failed module."""
    from models import PipelineRun, AnalysisTask
    from modules.reporting.pipeline_report import build_resume_plan
    from worker import submit_pipeline_run

    run = PipelineRun.get_by_id(run_id)
    if not run:
        return jsonify({'error': 'Pipeline run 不存在'}), 404
    if run.project_id != pid:
        return jsonify({'error': 'Pipeline run 不属于该项目'}), 403
    if run.status != 'failed':
        return jsonify({'error': '只有 failed 状态的 Pipeline run 可以续跑'}), 400

    try:
        task_ids = json.loads(run.task_ids_json) if run.task_ids_json else []
    except Exception:
        task_ids = []
    tasks = [t for tid in task_ids for t in [AnalysisTask.get_by_id(tid)] if t]
    resume_plan = build_resume_plan(run, tasks)
    if not resume_plan.get('can_resume'):
        return jsonify({'error': resume_plan.get('reason', '无法续跑'), 'resume_plan': resume_plan}), 400

    input_path = resume_plan['resume_input_path']
    is_valid, err = Config._validate_path(input_path, pid)
    if not is_valid:
        return jsonify({'error': err}), 400
    if not os.path.isfile(input_path):
        return jsonify({'error': '续跑输入文件不存在'}), 400

    modules = resume_plan['remaining_modules']
    params_by_module = resume_plan.get('params', {})
    new_run = PipelineRun(
        project_id=pid,
        name=f"{run.name} - resume",
        analysis_type=run.analysis_type,
        input_path=input_path,
        modules_json=json.dumps(modules),
        params_json=json.dumps(params_by_module, ensure_ascii=False),
    )
    new_run.save()
    submit_pipeline_run(
        run_id=new_run.id,
        project_id=pid,
        modules=modules,
        params_by_module=params_by_module,
        project_dir=Config.project_dir(pid),
        input_path=input_path,
    )
    return jsonify({
        'id': new_run.id,
        'message': '续跑流程已启动',
        'resume_plan': resume_plan,
    }), 201
