import os
from flask import Blueprint, render_template, request, redirect, url_for, flash, jsonify
from werkzeug.utils import secure_filename
from models import Project
from config import Config

upload_bp = Blueprint('upload', __name__)

ALLOWED_EXT = {'.h5ad', '.h5', '.csv', '.txt', '.mtx', '.gz', '.xlsx', '.xls', '.tsv'}

# 10x 文件名匹配模式
_10X_FILE_PATTERNS = {
    'barcodes': ['barcodes.tsv', 'barcodes.tsv.gz'],
    'genes_v2': ['genes.tsv', 'genes.tsv.gz'],
    'genes_v3': ['features.tsv', 'features.tsv.gz'],
    'matrix': ['matrix.mtx', 'matrix.mtx.gz'],
}

def check_10x_files(uploads_dir):
    """检测 uploads 目录中是否包含完整的 10x 三文件组合。"""
    if not os.path.isdir(uploads_dir):
        return {'has_10x': False}

    existing = set(os.listdir(uploads_dir))
    result = {'has_10x': False, 'files': {}, 'version': None}

    # 检测 barcodes
    barcodes = None
    for name in _10X_FILE_PATTERNS['barcodes']:
        if name in existing:
            barcodes = name
            break
    if not barcodes:
        return result

    # 检测 matrix
    matrix = None
    for name in _10X_FILE_PATTERNS['matrix']:
        if name in existing:
            matrix = name
            break
    if not matrix:
        return result

    # 检测 genes（优先 v3）
    genes = None
    version = None
    for name in _10X_FILE_PATTERNS['genes_v3']:
        if name in existing:
            genes = name
            version = 'v3'
            break
    if not genes:
        for name in _10X_FILE_PATTERNS['genes_v2']:
            if name in existing:
                genes = name
                version = 'v2'
                break
    if not genes:
        return result

    return {
        'has_10x': True,
        'files': {'barcodes': barcodes, 'genes': genes, 'matrix': matrix},
        'version': version,
    }

@upload_bp.route('/<pid>/upload', methods=['GET', 'POST'])
def upload(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('项目未找到', 'danger')
        return redirect(url_for('main.index'))
    if request.method == 'POST':
        if 'file' not in request.files:
            flash('请选择文件', 'danger')
            return redirect(url_for('upload.upload', pid=pid))
        file = request.files['file']
        if file.filename == '':
            flash('请选择文件', 'danger')
            return redirect(url_for('upload.upload', pid=pid))
        fname = secure_filename(file.filename)
        ext = os.path.splitext(fname)[1].lower()
        if ext not in ALLOWED_EXT and not fname.endswith('.mtx.gz') and not fname.endswith('.tsv.gz'):
            flash(f'不支持的文件格式: {ext}', 'danger')
            return redirect(url_for('upload.upload', pid=pid))
        uploads_dir = Config.uploads_dir(pid)
        os.makedirs(uploads_dir, exist_ok=True)
        fpath = os.path.join(uploads_dir, fname)
        file.save(fpath)
        p.status = 'data_ready'
        p.save()
        flash(f'文件 "{fname}" 上传成功', 'success')
        return redirect(url_for('projects.detail', pid=pid))
    return render_template('upload.html', project=p)


@upload_bp.route('/<pid>/upload/check-10x')
def check_10x(pid):
    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'has_10x': False, 'error': '项目未找到'}), 404
    uploads_dir = Config.uploads_dir(pid)
    return jsonify(check_10x_files(uploads_dir))


@upload_bp.route('/<pid>/upload/convert-10x', methods=['POST'])
def convert_10x(pid):
    import json as _json
    from models import AnalysisTask
    from worker import submit_task

    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'error': '项目未找到'}), 404

    uploads_dir = Config.uploads_dir(pid)
    check = check_10x_files(uploads_dir)
    if not check['has_10x']:
        return jsonify({'error': '未检测到完整的 10x 数据文件'}), 400

    species = request.form.get('species', '').strip() or None
    genome = request.form.get('genome', '').strip() or None

    task = AnalysisTask(
        project_id=pid,
        module_name='convert_10x',
        status='pending',
        params_json=_json.dumps({
            'mtx_dir': uploads_dir,
            'species': species,
            'genome': genome,
        }, ensure_ascii=False),
    )
    task.save()

    proj_dir = Config.project_dir(pid)
    submit_task(task.id, pid, 'convert_10x',
                {'mtx_dir': uploads_dir, 'species': species, 'genome': genome},
                proj_dir, uploads_dir)

    return jsonify({'task_id': task.id})
