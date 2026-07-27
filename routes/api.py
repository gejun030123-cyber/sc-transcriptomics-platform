import json
import os
import shutil
import struct
import base64
import uuid
import logging
import psutil
from flask import Blueprint, jsonify, request, send_file
from models import AnalysisTask, Project, ResultFile
from config import Config

logger = logging.getLogger(__name__)

api_bp = Blueprint('api', __name__)


PRESETS_GLOBAL_DIR = os.path.join(Config.DATA_DIR, 'presets', '_global')


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


def _ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def _load_preset(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        data['_filepath'] = filepath
        data['_id'] = os.path.splitext(os.path.basename(filepath))[0]
        return data
    except (json.JSONDecodeError, IOError):
        return None


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
    enrichment_files = [f for f in files if f.category == 'enrichment' and f.file_type in {'png', 'svg'}]
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
        from modules.io_utils import read_expression_matrix
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
            read_expression_matrix, infer_sample_group_candidates,
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
    _ensure_dir(PRESETS_GLOBAL_DIR)
    presets.extend(_list_presets_in_dir(PRESETS_GLOBAL_DIR, 'global'))
    if project_id:
        proj_dir = _get_project_presets_dir(project_id)
        presets.extend(_list_presets_in_dir(proj_dir, 'project'))
    if filter_type:
        presets = [p for p in presets if p.get('analysis_type') == filter_type]
    result = []
    for p in presets:
        result.append({
            'id': p['_id'],
            'name': p.get('name', ''),
            'description': p.get('description', ''),
            'analysis_type': p.get('analysis_type', ''),
            'scope': p['_scope'],
        })
    return jsonify({'presets': result})


@api_bp.route('/presets/<preset_id>')
def get_preset(preset_id):
    if '..' in preset_id or '/' in preset_id:
        return jsonify({'error': '无效的预设 ID'}), 400
    project_id = request.args.get('project_id', '')
    if project_id:
        fpath = os.path.join(_get_project_presets_dir(project_id), f'{preset_id}.json')
        p = _load_preset(fpath)
        if p:
            return jsonify({'preset': p})
    fpath = os.path.join(PRESETS_GLOBAL_DIR, f'{preset_id}.json')
    p = _load_preset(fpath)
    if p:
        return jsonify({'preset': p})
    return jsonify({'error': '预设不存在'}), 404


@api_bp.route('/presets', methods=['POST'])
def save_preset():
    data = request.get_json()
    if not data or not data.get('name'):
        return jsonify({'error': '预设名称不能为空'}), 400
    preset_id = str(uuid.uuid4())[:8]
    preset = {
        'name': data['name'],
        'description': data.get('description', ''),
        'analysis_type': data.get('analysis_type', ''),
        'params': data.get('params', {}),
        'pipeline': data.get('pipeline', {}),
        'filters': data.get('filters', {}),
        'visualization': data.get('visualization', {}),
    }
    scope = data.get('scope', 'project')
    project_id = data.get('project_id', '')
    if scope == 'global':
        target_dir = PRESETS_GLOBAL_DIR
    else:
        if not project_id:
            return jsonify({'error': '项目级预设需要 project_id'}), 400
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
    fpath = os.path.join(PRESETS_GLOBAL_DIR, f'{preset_id}.json')
    if os.path.isfile(fpath):
        os.remove(fpath)
        return jsonify({'message': '预设已删除'})
    return jsonify({'error': '预设不存在'}), 404


# ============ Pipeline Run API ============

@api_bp.route('/projects/<pid>/pipeline-runs', methods=['POST'])
def create_pipeline_run(pid):
    """创建并启动 pipeline run。"""
    from modules import MODULE_REGISTRY, SC_MODULE_NAMES, BULK_MODULE_NAMES, validate_pipeline_order
    from models import PipelineRun
    from worker import submit_pipeline_run

    # 校验项目
    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'error': '项目不存在'}), 404

    data = request.get_json()
    if not data:
        return jsonify({'error': '请求体为空'}), 400

    # 解析参数
    name = data.get('name', '').strip()
    analysis_type = data.get('analysis_type', '').strip()
    modules = data.get('modules', [])
    params = data.get('params', {})
    input_path = data.get('input_path', '').strip()

    # 校验必填字段
    if not name:
        return jsonify({'error': '缺少流程名称'}), 400
    if not analysis_type:
        return jsonify({'error': '缺少分析类型'}), 400
    if not modules:
        return jsonify({'error': '模块列表为空'}), 400
    if not input_path:
        return jsonify({'error': '缺少输入文件路径'}), 400

    # 校验 analysis_type
    if analysis_type not in ('sc', 'bulk'):
        return jsonify({'error': 'analysis_type 必须是 sc 或 bulk'}), 400

    # 校验模块类型
    module_set = SC_MODULE_NAMES if analysis_type == 'sc' else BULK_MODULE_NAMES
    for mod in modules:
        if mod not in MODULE_REGISTRY:
            return jsonify({'error': f'未知模块: {mod}'}), 400
        if mod not in module_set:
            return jsonify({'error': f'模块 {mod} 不属于 {analysis_type} 类型'}), 400

    # 校验依赖顺序
    is_valid, errors = validate_pipeline_order(modules)
    if not is_valid:
        return jsonify({'error': '模块顺序不满足依赖约束', 'details': errors}), 400

    # 校验输入文件路径
    abs_input = os.path.abspath(input_path)
    project_dir = os.path.abspath(Config.project_dir(pid))
    if not abs_input.startswith(project_dir + os.sep):
        return jsonify({'error': '输入文件不在项目目录内'}), 400
    if os.path.islink(abs_input):
        return jsonify({'error': '输入文件不能是符号链接'}), 400
    if not os.path.isfile(abs_input):
        return jsonify({'error': '输入文件不存在'}), 400

    # 构建参数（从 schema 默认值 + 请求参数合并）
    from modules.schemas import PARAM_SCHEMAS, filter_active_params
    params_by_module = {}
    for mod in modules:
        # 从 schema 获取默认值
        schema = PARAM_SCHEMAS.get(mod, [])
        default_params = {}
        for field in schema:
            key = field['key']
            default = field.get('default')
            if field.get('type') == 'select' and 'options' in field:
                opts = field['options']
                if default not in opts and opts:
                    default = opts[0]
            default_params[key] = default
        # 合并请求参数（新形态：按模块名分组）
        if mod in params:
            module_params = params[mod]
            if isinstance(module_params, dict):
                default_params.update(module_params)
        # 兼容旧形态：扁平参数应用到所有匹配的模块
        for key, val in params.items():
            if key not in ('qc', 'normalize', 'hvg', 'dimred', 'batch_correct',
                          'clustering', 'subcluster', 'qc_reassess', 'annotation', 'deg',
                          'trajectory', 'sc_timecourse', 'proportion', 'cell_communication',
                          'bulk_qc', 'bulk_normalize', 'bulk_deg', 'bulk_pca',
                          'bulk_heatmap', 'bulk_enrichment', 'bulk_timecourse',
                          'bulk_deg_integration', 'convert_10x'):
                if isinstance(val, (str, int, float, bool)):
                    # 检查该模块的 schema 是否有这个 key
                    schema_keys = {f['key'] for f in PARAM_SCHEMAS.get(mod, [])}
                    if key in schema_keys and key not in default_params:
                        default_params[key] = val
        params_by_module[mod] = filter_active_params(schema, default_params)

    # 创建 PipelineRun
    pipeline_run = PipelineRun(
        project_id=pid,
        name=name,
        analysis_type=analysis_type,
        input_path=input_path,
        modules_json=json.dumps(modules),
        params_json=json.dumps(params_by_module, ensure_ascii=False)
    )
    pipeline_run.save()

    # 提交到线程池
    submit_pipeline_run(
        run_id=pipeline_run.id,
        project_id=pid,
        modules=modules,
        params_by_module=params_by_module,
        project_dir=project_dir,
        input_path=abs_input
    )

    return jsonify({'id': pipeline_run.id, 'message': '流程已启动'}), 201


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
