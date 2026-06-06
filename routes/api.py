import json
import os
import shutil
import struct
import base64
import psutil
from flask import Blueprint, jsonify, request, send_file
from models import AnalysisTask, Project, ResultFile
from config import Config

api_bp = Blueprint('api', __name__)


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

@api_bp.route('/projects/<pid>/adata-info')
def adata_info(pid):
    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'error': 'Not found'}), 404
    adata_path = p.get_latest_adata_path()
    if not adata_path:
        uploads_dir = os.path.join(Config.DATA_DIR, 'projects', pid, 'uploads')
        if os.path.isdir(uploads_dir):
            for f in os.listdir(uploads_dir):
                if f.endswith('.h5ad'):
                    adata_path = os.path.join(uploads_dir, f)
                    break
    if not adata_path or not os.path.exists(adata_path):
        return jsonify({'error': 'No h5ad file found'}), 404
    try:
        import anndata
        import numpy as np
        adata = anndata.read_h5ad(adata_path, backed='r')

        def _dtype_str(series):
            return str(series.dtype)

        def _dtype_str_arr(arr):
            return str(arr.dtype)

        # obs 信息
        obs_info = []
        for col in adata.obs.columns:
            s = adata.obs[col]
            entry = {'name': col, 'dtype': _dtype_str(s)}
            if s.dtype == object or s.dtype.name == 'category':
                nuniq = s.nunique()
                entry['n_unique'] = int(nuniq)
                if nuniq <= 20:
                    entry['values'] = [str(v) for v in sorted(s.dropna().unique().tolist())]
            obs_info.append(entry)

        # var 信息
        var_info = []
        for col in adata.var.columns:
            s = adata.var[col]
            entry = {'name': col, 'dtype': _dtype_str(s)}
            if s.dtype == object or s.dtype.name == 'category':
                nuniq = s.nunique()
                entry['n_unique'] = int(nuniq)
                if nuniq <= 20:
                    entry['values'] = [str(v) for v in sorted(s.dropna().unique().tolist())[:20]]
            var_info.append(entry)

        # layers
        layers_keys = list(adata.layers.keys())

        # obsm
        obsm_info = []
        for k in adata.obsm.keys():
            arr = adata.obsm[k]
            obsm_info.append({'key': k, 'shape': list(arr.shape), 'dtype': _dtype_str_arr(arr)})

        # varm
        varm_info = []
        for k in adata.varm.keys():
            arr = adata.varm[k]
            varm_info.append({'key': k, 'shape': list(arr.shape), 'dtype': _dtype_str_arr(arr)})

        # obsp
        obsp_keys = list(adata.obsp.keys())

        # uns
        uns_info = []
        for k in adata.uns.keys():
            v = adata.uns[k]
            if isinstance(v, np.ndarray):
                uns_info.append({'key': k, 'type': 'ndarray', 'shape': list(v.shape), 'dtype': str(v.dtype)})
            elif isinstance(v, dict):
                uns_info.append({'key': k, 'type': 'dict', 'keys': list(v.keys())[:20]})
            else:
                val_str = str(v)[:200]
                uns_info.append({'key': k, 'type': type(v).__name__, 'preview': val_str})

        # var_names 样本
        var_names_sample = list(adata.var_names[:min(20, adata.n_vars)])
        obs_names_sample = list(adata.obs_names[:min(10, adata.n_obs)])

        info = {
            'n_obs': int(adata.n_obs),
            'n_vars': int(adata.n_vars),
            'file': os.path.basename(adata_path),
            'obs_columns': obs_info,
            'var_columns': var_info,
            'layers': layers_keys,
            'obsm': obsm_info,
            'varm': varm_info,
            'obsp': obsp_keys,
            'uns': uns_info,
            'var_names_sample': [str(v) for v in var_names_sample],
            'obs_names_sample': [str(v) for v in obs_names_sample],
        }
        adata.file.close()
        return jsonify(info)
    except Exception as e:
        import traceback
        return jsonify({'error': str(e), 'traceback': traceback.format_exc()}), 500

@api_bp.route('/result-file/<file_id>')
def get_result_file(file_id):
    f = ResultFile.get_by_id(file_id)
    if not f:
        return jsonify({'error': 'Not found'}), 404
    if f.file_type == 'plotly_json':
        with open(f.file_path, 'r') as fh:
            data = json.load(fh)
        data = _decode_plotly_binary(data)
        return jsonify(data)
    return send_file(f.file_path)


@api_bp.route('/enrichment-result/<task_id>')
def enrichment_result(task_id):
    """返回富集分析的 Plotly JSON 结果"""
    from models import ResultFile
    files = ResultFile.get_by_task(task_id)
    enrichment_files = [f for f in files if f.category == 'enrichment']
    result = []
    for f in enrichment_files:
        with open(f.file_path, 'r') as fh:
            data = json.load(fh)
        result.append({'id': f.id, 'label': f.label, 'data': data})
    return jsonify(result)


@api_bp.route('/column-values')
def column_values():
    file_path = request.args.get('file_path', '')
    column = request.args.get('column', '')
    if not file_path or not column:
        return jsonify({'values': []})
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


@api_bp.route('/data-info')
def data_info():
    """返回输入文件的基本数据信息（样本数、基因数、样本名、注释列）"""
    file_path = request.args.get('file_path', '')
    if not file_path:
        return jsonify({'error': 'No file path'}), 400
    try:
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(file_path)
        qc_columns = {'total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'size_factor', 'total_counts_mt'}
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
    if not file_path:
        return jsonify({'columns': [], 'sample_groups': {}})
    try:
        from modules.io_utils import read_expression_matrix
        import re
        adata = read_expression_matrix(file_path)
        qc_columns = {'total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'size_factor', 'total_counts_mt'}
        cols = [c for c in adata.obs.columns if c not in qc_columns]

        # 如果 obs 有注释列，直接返回
        if cols:
            return jsonify({'columns': cols, 'sample_groups': {}})

        # 如果 obs 没有注释列，尝试从样本名中提取分组前缀
        sample_groups = {}
        if adata.n_obs > 0:
            sample_names = adata.obs.index.tolist()
            prefixes = []
            valid_indices = []
            for i, name in enumerate(sample_names):
                name_str = str(name)
                name_clean = re.sub(r'_(count|FPKM|TPM|fpkm|tpm|Counts|normalized)$', '', name_str)
                prefix = re.sub(r'[-_]\d+.*$', '', name_clean)
                # 只保留匹配 prefix-number 模式的样本名（排除注释列）
                if prefix and re.match(r'^[a-zA-Z][a-zA-Z0-9]*[-_]\d', name_clean):
                    prefixes.append(prefix)
                    valid_indices.append(i)
                else:
                    prefixes.append(None)
            unique_prefixes = sorted(set(p for p in prefixes if p is not None))
            if len(unique_prefixes) > 1 and len(valid_indices) > len(unique_prefixes):
                mapping = {str(sample_names[i]): prefixes[i] for i in valid_indices}
                sample_groups['auto_group'] = {
                    'values': unique_prefixes,
                    'mapping': mapping
                }

        return jsonify({'columns': cols, 'sample_groups': sample_groups})
    except Exception:
        return jsonify({'columns': [], 'sample_groups': {}})


@api_bp.route('/data-info-full')
def data_info_full():
    """返回 h5ad 文件的完整结构信息（obs/var 详情、layers、obsm、varm、obsp、uns）"""
    file_path = request.args.get('file_path', '')
    if not file_path:
        return jsonify({'error': 'No file path'}), 400
    if not file_path.endswith('.h5ad'):
        return jsonify({'error': 'Not an h5ad file'}), 400
    try:
        import anndata
        import numpy as np
        adata = anndata.read_h5ad(file_path, backed='r')

        def _dtype_str(series):
            return str(series.dtype)

        obs_info = []
        for col in adata.obs.columns:
            s = adata.obs[col]
            entry = {'name': col, 'dtype': _dtype_str(s)}
            if s.dtype == object or s.dtype.name == 'category':
                nuniq = s.nunique()
                entry['n_unique'] = int(nuniq)
                if nuniq <= 20:
                    entry['values'] = [str(v) for v in sorted(s.dropna().unique().tolist())]
            obs_info.append(entry)

        var_info = []
        for col in adata.var.columns:
            s = adata.var[col]
            entry = {'name': col, 'dtype': _dtype_str(s)}
            if s.dtype == object or s.dtype.name == 'category':
                nuniq = s.nunique()
                entry['n_unique'] = int(nuniq)
                if nuniq <= 20:
                    entry['values'] = [str(v) for v in sorted(s.dropna().unique().tolist())[:20]]
            var_info.append(entry)

        layers_keys = list(adata.layers.keys())

        obsm_info = []
        for k in adata.obsm.keys():
            arr = adata.obsm[k]
            obsm_info.append({'key': k, 'shape': list(arr.shape), 'dtype': str(arr.dtype)})

        varm_info = []
        for k in adata.varm.keys():
            arr = adata.varm[k]
            varm_info.append({'key': k, 'shape': list(arr.shape), 'dtype': str(arr.dtype)})

        obsp_keys = list(adata.obsp.keys())

        uns_info = []
        for k in adata.uns.keys():
            v = adata.uns[k]
            if isinstance(v, np.ndarray):
                uns_info.append({'key': k, 'type': 'ndarray', 'shape': list(v.shape), 'dtype': str(v.dtype)})
            elif isinstance(v, dict):
                uns_info.append({'key': k, 'type': 'dict', 'keys': list(v.keys())[:20]})
            else:
                uns_info.append({'key': k, 'type': type(v).__name__, 'preview': str(v)[:200]})

        var_names_sample = [str(v) for v in adata.var_names[:min(20, adata.n_vars)]]

        info = {
            'obs_columns': obs_info,
            'var_columns': var_info,
            'layers': layers_keys,
            'obsm': obsm_info,
            'varm': varm_info,
            'obsp': obsp_keys,
            'uns': uns_info,
            'var_names_sample': var_names_sample,
        }
        adata.file.close()
        return jsonify(info)
    except Exception as e:
        import traceback
        return jsonify({'error': str(e), 'traceback': traceback.format_exc()}), 500
