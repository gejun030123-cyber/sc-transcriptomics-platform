import os
from flask import Blueprint, render_template, request, redirect, url_for, flash, jsonify
from werkzeug.utils import secure_filename
from models import Project
from config import Config

upload_bp = Blueprint('upload', __name__)

ALLOWED_EXT = {
    '.h5ad', '.h5', '.hdf5', '.loom', '.zarr',
    '.csv', '.txt', '.mtx', '.gz', '.xlsx', '.xls', '.tsv',
    '.zip',
}

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


def check_importable_sc_files(uploads_dir):
    """检测 uploads 目录中可统一导入为 h5ad 的单细胞文件。"""
    if not os.path.isdir(uploads_dir):
        return {'has_importable': False, 'files': []}

    from modules.io_utils import infer_sc_data_format

    supported = []
    skip_names = {'converted_10x.h5ad'}
    for name in sorted(os.listdir(uploads_dir)):
        if name in skip_names or name.endswith('_imported.h5ad'):
            continue
        path = os.path.join(uploads_dir, name)
        fmt = infer_sc_data_format(path)
        if fmt in {'h5ad', '10x_h5', 'loom', 'zarr'}:
            try:
                size_mb = round(os.path.getsize(path) / (1024 * 1024), 1) if os.path.isfile(path) else None
            except OSError:
                size_mb = None
            supported.append({
                'name': name,
                'format': fmt,
                'size_mb': size_mb,
            })

    return {'has_importable': bool(supported), 'files': supported}

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


@upload_bp.route('/<pid>/upload/check-sc-import')
def check_sc_import(pid):
    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'has_importable': False, 'error': '项目未找到'}), 404
    uploads_dir = Config.uploads_dir(pid)
    return jsonify(check_importable_sc_files(uploads_dir))


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


@upload_bp.route('/<pid>/upload/import-sc', methods=['POST'])
def import_sc(pid):
    import json as _json
    from models import AnalysisTask
    from worker import submit_task
    from modules.io_utils import infer_sc_data_format

    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'error': '项目未找到'}), 404

    uploads_dir = Config.uploads_dir(pid)
    source_file = request.form.get('source_file', '').strip()
    if not source_file:
        return jsonify({'error': '请选择要导入的单细胞文件'}), 400
    if os.path.basename(source_file) != source_file:
        return jsonify({'error': '文件名不合法'}), 400

    source_path = os.path.join(uploads_dir, source_file)
    if not os.path.exists(source_path):
        return jsonify({'error': f'文件不存在: {source_file}'}), 404

    input_format = request.form.get('input_format', '').strip() or infer_sc_data_format(source_path)
    if input_format not in {'h5ad', '10x_h5', 'loom', 'zarr', 'auto'}:
        return jsonify({'error': f'不支持的单细胞导入格式: {input_format}'}), 400

    species = request.form.get('species', '').strip() or None
    genome = request.form.get('genome', '').strip() or None

    params = {
        'source_path': source_path,
        'input_format': input_format,
        'species': species,
        'genome': genome,
    }
    task = AnalysisTask(
        project_id=pid,
        module_name='convert_10x',
        status='pending',
        params_json=_json.dumps(params, ensure_ascii=False),
    )
    task.save()

    proj_dir = Config.project_dir(pid)
    submit_task(task.id, pid, 'convert_10x', params, proj_dir, source_path)

    return jsonify({'task_id': task.id})


@upload_bp.route('/<pid>/upload/import-10x-batches', methods=['POST'])
def import_10x_batches(pid):
    """Upload and merge two zipped 10x batches with an explicit batch label."""
    import json as _json
    from models import AnalysisTask
    from worker import submit_task

    p = Project.get_by_id(pid)
    if not p:
        return jsonify({'error': '项目未找到'}), 404

    files = request.files.getlist('batch_zip')
    if not files:
        files = [item for item in (
            request.files.get('batch_a_zip') or request.files.get('batch_a'),
            request.files.get('batch_b_zip') or request.files.get('batch_b'),
        ) if item]
    if len(files) < 2:
        return jsonify({'error': '请上传两组 ZIP 文件'}), 400
    if len(files) > 2:
        return jsonify({'error': '当前接口只支持两组 ZIP 文件'}), 400

    names = request.form.getlist('batch_name')
    if not names:
        names = [request.form.get('batch_a_name', ''), request.form.get('batch_b_name', '')]
    while len(names) < 2:
        names.append('')
    default_names = []
    for index, file in enumerate(files, start=1):
        original = secure_filename(file.filename or '')
        if not original.lower().endswith('.zip'):
            return jsonify({'error': f'第 {index} 个文件不是 ZIP'}), 400
        stem = os.path.splitext(original)[0] or f'batch_{index}'
        name = names[index - 1].strip() or stem
        name = ''.join(ch if ch.isalnum() or ch in '._-' else '_' for ch in name).strip('._-')
        if not name:
            name = f'batch_{index}'
        default_names.append(name)
    if len(set(default_names)) != len(default_names):
        return jsonify({'error': '两组批次名称不能相同'}), 400

    species = request.form.get('species', '').strip() or None
    genome = request.form.get('genome', '').strip() or None
    uploads_dir = Config.uploads_dir(pid)
    zip_dir = os.path.join(uploads_dir, 'batch_zips')
    os.makedirs(zip_dir, exist_ok=True)
    batch_sources = []
    for index, (file, batch_name) in enumerate(zip(files, default_names), start=1):
        original = secure_filename(file.filename or '')
        stored_name = f'{index}_{batch_name}_{original}'
        zip_path = os.path.join(zip_dir, stored_name)
        file.save(zip_path)
        batch_sources.append({'zip_path': zip_path, 'batch_name': batch_name})

    params = {
        'batch_sources': batch_sources,
        'species': species,
        'genome': genome,
    }
    task = AnalysisTask(
        project_id=pid,
        module_name='convert_10x',
        status='pending',
        params_json=_json.dumps(params, ensure_ascii=False),
    )
    task.save()
    p.status = 'processing'
    p.save()
    submit_task(task.id, pid, 'convert_10x', params,
                Config.project_dir(pid), uploads_dir)
    return jsonify({'task_id': task.id, 'batch_names': default_names})
