import os
import uuid
from flask import Blueprint, render_template, request, redirect, url_for, flash, jsonify, send_file
from werkzeug.utils import secure_filename
from models import Project
from config import Config

upload_bp = Blueprint('upload', __name__)

ALLOWED_EXT = {
    '.h5ad', '.h5', '.hdf5', '.loom', '.zarr',
    '.csv', '.txt', '.mtx', '.gz', '.xlsx', '.xls', '.tsv',
    '.zip',
    # CellOracle virtual KO accepts only tabular TF-info matrices.  Never
    # accept project-uploaded pickle/Oracle/Links files for deserialization.
    '.parquet', '.pq',
}


def _unique_upload_path(directory, stem, suffix):
    """Allocate an import artifact name without trusting a user-supplied path."""
    os.makedirs(directory, exist_ok=True)
    safe_stem = secure_filename(stem) or 'bulk_counts'
    token = uuid.uuid4().hex[:10]
    return os.path.join(directory, f'{safe_stem}_{token}{suffix}')


def _upload_suffix(filename):
    """Keep compound annotation suffixes such as ``.gtf.gz`` intact."""
    name = str(filename or '')
    stem, suffix = os.path.splitext(name)
    if suffix.lower() == '.gz':
        _, annotation_suffix = os.path.splitext(stem)
        if annotation_suffix.lower() in {'.gtf', '.gff', '.gff3'}:
            return annotation_suffix.lower() + '.gz'
    return suffix.lower()

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
    # A 10x Matrix Market dataset must be imported as its complete three-file
    # set.  Do not mistakenly offer its barcodes/features files as standalone
    # expression matrices merely because they are TSV files.
    skip_names = {
        'converted_10x.h5ad',
        'matrix.mtx', 'matrix.mtx.gz',
        'barcodes.tsv', 'barcodes.tsv.gz',
        'features.tsv', 'features.tsv.gz',
        'genes.tsv', 'genes.tsv.gz',
    }
    for name in sorted(os.listdir(uploads_dir)):
        if name in skip_names or name.endswith('_imported.h5ad'):
            continue
        path = os.path.join(uploads_dir, name)
        fmt = infer_sc_data_format(path)
        if fmt in {'h5ad', '10x_h5', 'loom', 'zarr', 'expression_matrix'}:
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
        ajax_request = request.headers.get('X-Requested-With') == 'XMLHttpRequest'

        def reject(message, status=400):
            if ajax_request:
                return jsonify({'error': message}), status
            flash(message, 'danger')
            return redirect(url_for('upload.upload', pid=pid))

        if 'file' not in request.files:
            return reject('请选择文件')
        file = request.files['file']
        if file.filename == '':
            return reject('请选择文件')
        fname = secure_filename(file.filename)
        if not fname:
            return reject('文件名不合法')
        ext = os.path.splitext(fname)[1].lower()
        if ext not in ALLOWED_EXT and not fname.endswith('.mtx.gz') and not fname.endswith('.tsv.gz'):
            return reject(f'不支持的文件格式: {ext or "无扩展名"}')
        uploads_dir = Config.uploads_dir(pid)
        os.makedirs(uploads_dir, exist_ok=True)
        fpath = os.path.join(uploads_dir, fname)
        file.save(fpath)
        p.status = 'data_ready'
        p.save()
        if ajax_request:
            return jsonify({
                'ok': True,
                'filename': fname,
                'size_bytes': os.path.getsize(fpath),
            }), 201
        flash(f'文件 "{fname}" 上传成功', 'success')
        return redirect(url_for('projects.detail', pid=pid))
    return render_template('upload.html', project=p)


@upload_bp.route('/<pid>/upload/import-bulk-counts', methods=['POST'])
def import_bulk_counts(pid):
    """Import a raw count matrix and explicit sample metadata as one h5ad input."""
    project = Project.get_by_id(pid)
    if not project:
        return jsonify({'error': '项目未找到'}), 404

    count_file = request.files.get('count_matrix')
    metadata_file = request.files.get('sample_metadata')
    annotation_file = request.files.get('gene_annotation')
    if not count_file or not count_file.filename:
        return jsonify({'error': '请选择原始 counts 矩阵。'}), 400
    if not metadata_file or not metadata_file.filename:
        return jsonify({'error': '请选择样本信息表。'}), 400

    count_name = secure_filename(count_file.filename)
    metadata_name = secure_filename(metadata_file.filename)
    annotation_name = secure_filename(annotation_file.filename) if annotation_file and annotation_file.filename else ''
    if not count_name or not metadata_name or (annotation_file and annotation_file.filename and not annotation_name):
        return jsonify({'error': '文件名不合法。'}), 400

    from modules.bulk_import import (
        build_bulk_counts_adata, gene_annotation_extension_allowed,
        tabular_extension_allowed,
    )
    if not tabular_extension_allowed(count_name) or not tabular_extension_allowed(metadata_name):
        return jsonify({'error': 'counts 矩阵和样本信息表仅支持 CSV、TSV、TXT 或 Excel。'}), 400
    if annotation_name and not gene_annotation_extension_allowed(annotation_name):
        return jsonify({
            'error': '基因注释仅支持 CSV、TSV、TXT、Excel、GTF/GFF（可为 .gz）。'
        }), 400

    uploads_dir = Config.uploads_dir(pid)
    source_dir = os.path.join(uploads_dir, 'bulk_count_imports')
    os.makedirs(source_dir, exist_ok=True)
    count_path = _unique_upload_path(
        source_dir, os.path.splitext(count_name)[0], os.path.splitext(count_name)[1].lower(),
    )
    metadata_path = _unique_upload_path(
        source_dir, os.path.splitext(metadata_name)[0], _upload_suffix(metadata_name),
    )
    count_file.save(count_path)
    metadata_file.save(metadata_path)
    annotation_path = None
    if annotation_name:
        annotation_path = _unique_upload_path(
            source_dir, os.path.splitext(annotation_name)[0], _upload_suffix(annotation_name),
        )
        annotation_file.save(annotation_path)

    try:
        adata = build_bulk_counts_adata(count_path, metadata_path, annotation_path)
        output_path = _unique_upload_path(
            uploads_dir, f'bulk_raw_counts_{os.path.splitext(count_name)[0]}', '.h5ad',
        )
        adata.write_h5ad(output_path)
    except (OSError, ValueError) as exc:
        return jsonify({'error': str(exc)}), 400

    project.status = 'data_ready'
    project.save()
    import_info = dict(adata.uns.get('bulk_import', {}) or {})
    return jsonify({
        'ok': True,
        'output_file': os.path.basename(output_path),
        'n_samples': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'n_gene_names': int(import_info.get('n_gene_names', 0)),
        'gene_name_source': import_info.get('gene_name_source', 'none'),
        'conditions': import_info.get('condition_sizes', {}),
        # The importer establishes a valid count/design contract, but it does
        # not replace QC or low-expression filtering.  Start the guided Bulk
        # chain at QC so an uploaded raw table cannot jump straight to DEG.
        'next_url': url_for('analysis.analyze', pid=pid, module_name='bulk_qc'),
    }), 201


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
    if input_format not in {
        'h5ad', '10x_h5', 'loom', 'zarr', 'expression_matrix', 'auto',
    }:
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
    """Upload and merge two or more zipped 10x batches with explicit labels."""
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
        return jsonify({'error': '请至少上传两组 ZIP 文件'}), 400

    names = request.form.getlist('batch_name')
    if not names:
        names = [request.form.get('batch_a_name', ''), request.form.get('batch_b_name', '')]
    while len(names) < len(files):
        names.append('')
    sample_ids = request.form.getlist('sample_id')
    while len(sample_ids) < len(files):
        sample_ids.append('')
    conditions = request.form.getlist('condition')
    while len(conditions) < len(files):
        conditions.append('')
    default_names = []
    resolved_sample_ids = []
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
        sample_id = sample_ids[index - 1].strip() or name
        sample_id = ''.join(ch if ch.isalnum() or ch in '._-' else '_' for ch in sample_id).strip('._-')
        if not sample_id:
            sample_id = name
        resolved_sample_ids.append(sample_id)
    if len(set(default_names)) != len(default_names):
        return jsonify({'error': '批次名称不能重复'}), 400
    if len(set(resolved_sample_ids)) != len(resolved_sample_ids):
        return jsonify({'error': '样本 ID（sample_id）不能重复'}), 400

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
        batch_sources.append({
            'zip_path': zip_path,
            'batch_name': batch_name,
            'sample_id': resolved_sample_ids[index - 1],
            'condition': conditions[index - 1].strip(),
        })

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


@upload_bp.route('/<pid>/upload/discover-10x-directory', methods=['POST'])
def discover_10x_directory(pid):
    """Create an editable manifest template from an approved server data root."""
    from modules.sc_batch import _available_path, discovered_manifest_frame

    project = Project.get_by_id(pid)
    if not project:
        return jsonify({'error': '项目未找到'}), 404
    source_root = (request.form.get('source_root') or '').strip()
    try:
        source_root = Config.validate_sc_batch_source_path(source_root)
        manifest = discovered_manifest_frame(source_root)
    except (OSError, ValueError) as exc:
        return jsonify({'error': str(exc)}), 400
    if manifest.empty:
        return jsonify({'error': '该目录下未发现完整的 10x matrix.mtx + barcodes + features/genes 文件'}), 400

    manifest_dir = os.path.join(Config.uploads_dir(pid), 'sc_batch_manifests')
    path = _available_path(manifest_dir, 'sc_batch_manifest_template', '.csv')
    manifest.to_csv(path, index=False)
    return jsonify({
        'manifest_file': os.path.basename(path),
        'n_samples': int(len(manifest)),
        'sample_ids': manifest['sample_id'].tolist(),
        'download_url': url_for('upload.download_sc_batch_manifest', pid=pid,
                                filename=os.path.basename(path)),
        'message': '已生成模板。请填写 condition 和 replicate 后再上传并导入；目录名不会被自动当作分组。',
    })


@upload_bp.route('/<pid>/upload/sc-batch-manifest/<filename>')
def download_sc_batch_manifest(pid, filename):
    """Download only a platform-generated manifest from the project upload area."""
    if not Project.get_by_id(pid):
        return jsonify({'error': '项目未找到'}), 404
    if secure_filename(filename) != filename:
        return jsonify({'error': '文件名不合法'}), 400
    path = os.path.join(Config.uploads_dir(pid), 'sc_batch_manifests', filename)
    valid, error = Config._validate_path(path, pid)
    if not valid or not os.path.isfile(path):
        return jsonify({'error': error or '模板不存在'}), 404
    return send_file(path, as_attachment=True, download_name=filename)


@upload_bp.route('/<pid>/upload/import-10x-manifest', methods=['POST'])
def import_10x_manifest(pid):
    """Submit a manifest-driven, arbitrary-size 10x import task."""
    import json as _json
    from models import AnalysisTask
    from worker import submit_task
    from modules.sc_batch import _available_path, read_sample_manifest

    project = Project.get_by_id(pid)
    if not project:
        return jsonify({'error': '项目未找到'}), 404
    source_root = (request.form.get('source_root') or '').strip()
    dataset_name = (request.form.get('dataset_name') or 'sc_batch').strip()
    species = request.form.get('species', '').strip() or None
    genome = request.form.get('genome', '').strip() or None
    manifest_dir = os.path.join(Config.uploads_dir(pid), 'sc_batch_manifests')
    os.makedirs(manifest_dir, exist_ok=True)

    uploaded = request.files.get('manifest_file')
    if uploaded and uploaded.filename:
        original = secure_filename(uploaded.filename)
        ext = os.path.splitext(original)[1].lower()
        if ext not in {'.csv', '.tsv', '.txt', '.xlsx', '.xls'}:
            return jsonify({'error': 'manifest 仅支持 CSV、TSV 或 Excel'}), 400
        manifest_path = _available_path(manifest_dir, os.path.splitext(original)[0] or 'sample_manifest', ext)
        uploaded.save(manifest_path)
    else:
        selected = (request.form.get('manifest_name') or '').strip()
        if not selected or secure_filename(selected) != selected:
            return jsonify({'error': '请上传已填写的样本 manifest'}), 400
        manifest_path = os.path.join(manifest_dir, selected)

    try:
        root = Config.validate_sc_batch_source_path(source_root)
        manifest = read_sample_manifest(manifest_path, root)
    except (OSError, ValueError) as exc:
        return jsonify({'error': str(exc)}), 400

    params = {
        'manifest_path': manifest_path,
        'source_root': root,
        'dataset_name': dataset_name,
        'species': species,
        'genome': genome,
    }
    task = AnalysisTask(
        project_id=pid,
        module_name='sc_batch_import',
        status='pending',
        params_json=_json.dumps(params, ensure_ascii=False),
    )
    task.save()
    project.status = 'processing'
    project.save()
    submit_task(task.id, pid, 'sc_batch_import', params,
                Config.project_dir(pid), manifest_path)
    return jsonify({
        'task_id': task.id,
        'n_samples': int(len(manifest)),
        'conditions': sorted(manifest['condition'].astype(str).unique().tolist()),
    })
