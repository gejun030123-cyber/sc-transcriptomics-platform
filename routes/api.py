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
    return jsonify({
        'status': t.status,
        'progress': t.progress,
        'progress_message': t.progress_message,
        'error': t.error_traceback,
        'started_at': t.started_at,
        'finished_at': t.finished_at,
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
        adata = anndata.read_h5ad(adata_path, backed='r')
        info = {
            'n_obs': int(adata.n_obs),
            'n_vars': int(adata.n_vars),
            'obs_columns': list(adata.obs.columns),
            'obsm_keys': list(adata.obsm.keys()),
            'file': os.path.basename(adata_path),
        }
        adata.file.close()
        return jsonify(info)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

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
