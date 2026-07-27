"""Routes for non-destructive figure selection, preview and versioning."""

from pathlib import Path
from uuid import uuid4

from flask import Blueprint, flash, jsonify, redirect, render_template, request, send_file, url_for
from werkzeug.utils import secure_filename

from config import Config
from models import FigureAsset, FigureVersion, Project
from modules.figure_studio import (
    FigureStudioError,
    list_sources,
    render_preview_data_uri,
    resolve_source,
    save_figure_version,
    validate_uploaded_image,
)


figure_studio_bp = Blueprint('figure_studio', __name__)


def _project_or_redirect(pid):
    project = Project.get_by_id(pid)
    if not project:
        flash('项目未找到。', 'danger')
        return None
    return project


def _version_payload(pid, version):
    return {
        'id': version.id,
        'label': version.label,
        'edit_mode': version.edit_mode,
        'png_url': (f'/projects/{pid}/figure-studio/versions/{version.id}/png'
                    if version.png_path else ''),
        'svg_url': (f'/projects/{pid}/figure-studio/versions/{version.id}/svg'
                    if version.svg_path else ''),
    }


@figure_studio_bp.route('/<pid>/figure-studio')
def studio(pid):
    project = _project_or_redirect(pid)
    if not project:
        return redirect(url_for('main.index'))
    return render_template(
        'figure_studio.html', project=project, sources=list_sources(pid),
        versions=[_version_payload(pid, version) for version in FigureVersion.get_by_project(pid)],
        selected_kind=request.args.get('source_kind', ''),
        selected_id=request.args.get('source_id', ''),
    )


@figure_studio_bp.route('/<pid>/figure-studio/sources')
def sources(pid):
    if not _project_or_redirect(pid):
        return jsonify({'error': '项目未找到。'}), 404
    return jsonify({'sources': list_sources(pid)})


@figure_studio_bp.route('/<pid>/figure-studio/upload', methods=['POST'])
def upload(pid):
    if not _project_or_redirect(pid):
        return jsonify({'error': '项目未找到。'}), 404
    image = request.files.get('image')
    if image is None or not image.filename:
        return jsonify({'error': '请选择需要上传的图片。'}), 400
    try:
        extension = validate_uploaded_image(image)
        safe_name = secure_filename(image.filename)
        label = (Path(safe_name).stem or '上传图片')[:120]
        output_dir = Path(Config.project_dir(pid)) / 'figure_studio' / 'uploads'
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / f'{uuid4().hex[:12]}{extension}'
        image.save(destination)
        asset = FigureAsset.create(pid, label, extension.lstrip('.'), str(destination))
    except FigureStudioError as exc:
        return jsonify({'error': str(exc)}), 400
    return jsonify({
        'source': {
            'kind': 'asset', 'id': asset.id, 'label': asset.label,
            'file_type': asset.file_type, 'edit_mode': 'style_only',
            'preview_url': f'/projects/{pid}/figure-studio/assets/{asset.id}',
        }
    }), 201


@figure_studio_bp.route('/<pid>/figure-studio/preview', methods=['POST'])
def preview(pid):
    if not _project_or_redirect(pid):
        return jsonify({'error': '项目未找到。'}), 404
    payload = request.get_json(silent=True) or {}
    try:
        source = resolve_source(pid, payload.get('source_kind', ''), payload.get('source_id', ''))
        rendered = render_preview_data_uri(source, payload.get('style') or {})
    except FigureStudioError as exc:
        return jsonify({'error': str(exc)}), 400
    except Exception:
        return jsonify({'error': '预览生成失败，请检查原图和可编辑数据是否完整。'}), 500
    return jsonify(rendered)


@figure_studio_bp.route('/<pid>/figure-studio/save', methods=['POST'])
def save(pid):
    if not _project_or_redirect(pid):
        return jsonify({'error': '项目未找到。'}), 404
    payload = request.get_json(silent=True) or {}
    try:
        source = resolve_source(pid, payload.get('source_kind', ''), payload.get('source_id', ''))
        version = save_figure_version(pid, source, payload.get('style') or {}, payload.get('label', ''))
    except FigureStudioError as exc:
        return jsonify({'error': str(exc)}), 400
    except Exception:
        return jsonify({'error': '保存图形版本失败，请稍后重试。'}), 500
    return jsonify({'version': _version_payload(pid, version)}), 201


@figure_studio_bp.route('/<pid>/figure-studio/assets/<asset_id>')
def view_asset(pid, asset_id):
    asset = FigureAsset.get_by_id(asset_id)
    if not asset or asset.project_id != pid:
        flash('上传图片未找到。', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    valid, _ = Config._validate_path(asset.file_path, pid)
    if not valid or not Path(asset.file_path).is_file():
        flash('上传图片路径无效。', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    return send_file(asset.file_path)


@figure_studio_bp.route('/<pid>/figure-studio/versions/<version_id>/<file_type>')
def view_version(pid, version_id, file_type):
    if file_type not in {'png', 'svg'}:
        return jsonify({'error': '不支持的版本格式。'}), 404
    version = FigureVersion.get_by_id(version_id)
    if not version or version.project_id != pid:
        flash('图形版本未找到。', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    file_path = version.png_path if file_type == 'png' else version.svg_path
    valid, _ = Config._validate_path(file_path, pid)
    if not file_path or not valid or not Path(file_path).is_file():
        return jsonify({'error': '该图形版本没有此格式。'}), 404
    return send_file(file_path)
