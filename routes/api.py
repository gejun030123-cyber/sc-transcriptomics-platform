import json
import os
import shutil
import psutil
from flask import Blueprint, jsonify, request, send_file
from models import AnalysisTask, Project, ResultFile
from config import Config

api_bp = Blueprint('api', __name__)

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
            return jsonify(json.load(fh))
    return send_file(f.file_path)
