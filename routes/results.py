import json
import os
import io
import re
import zipfile
from collections.abc import Mapping, Sequence
from numbers import Number
from flask import Blueprint, render_template, send_file, flash, redirect, url_for
from models import Project, AnalysisTask, ResultFile, PipelineRun
from modules.schemas import MODULE_DISPLAY_MAP, STATUS_MAP
from modules import PIPELINE_ORDER
from config import Config


# ``result_json`` is an audit record shared by the worker, reports and APIs.
# It may therefore contain paths, checksums, per-batch dictionaries or long
# lists that are useful for reproducibility but unsuitable for a result-page
# summary. Keep the stored result intact and apply this presentation policy
# only at the HTML boundary.
_RESULT_SUMMARY_MAX_ITEMS = 16
_RESULT_VALUE_MAX_CHARS = 96
_RESULT_VALUE_MAX_COLLECTION_ITEMS = 5
_TECHNICAL_RESULT_KEY_TOKENS = {
    'path', 'paths', 'file', 'files', 'dir', 'directory', 'checksum',
    'sha256', 'hash', 'fingerprint', 'manifest', 'provenance', 'request_id',
    'request_hash',
}


def _result_key_is_technical(key):
    """Whether a result key is an implementation detail, not a page metric."""
    tokens = set(re.findall(r'[a-z0-9]+', str(key or '').lower()))
    return bool(tokens & _TECHNICAL_RESULT_KEY_TOKENS) or str(key or '').startswith('_')


def _looks_like_absolute_path(value):
    value = str(value or '').strip()
    return (
        value.startswith(('/', '\\\\'))
        or bool(re.match(r'^[A-Za-z]:[\\\\/]', value))
    )


def _compact_result_value(value, *, max_chars=_RESULT_VALUE_MAX_CHARS):
    """Return a safe, single-line result preview, or ``None`` when omitted.

    A compact scalar and a very small scalar collection are meaningful in a
    summary card. Nested/large data belongs in registered result files or the
    manifest, where it can be inspected without making the page unreadable.
    """
    if value is None:
        return None
    if isinstance(value, bool):
        return '是' if value else '否'
    if isinstance(value, Number):
        return str(value)
    if isinstance(value, str):
        text = ' '.join(value.split())
        if not text or len(text) > max_chars or _looks_like_absolute_path(text):
            return None
        return text

    if isinstance(value, Mapping):
        if not value or len(value) > _RESULT_VALUE_MAX_COLLECTION_ITEMS:
            return None
        entries = []
        child_limit = max(24, max_chars // _RESULT_VALUE_MAX_COLLECTION_ITEMS)
        for key, child in value.items():
            is_nested = (
                isinstance(child, (Mapping, Sequence))
                and not isinstance(child, (str, bytes, bytearray))
            )
            if _result_key_is_technical(key) or is_nested:
                return None
            compact_child = _compact_result_value(child, max_chars=child_limit)
            if compact_child is None:
                return None
            entries.append(f'{str(key).replace("_", " ")}: {compact_child}')
        text = '；'.join(entries)
        return text if len(text) <= max_chars else None

    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        if not value or len(value) > _RESULT_VALUE_MAX_COLLECTION_ITEMS:
            return None
        child_limit = max(24, max_chars // _RESULT_VALUE_MAX_COLLECTION_ITEMS)
        entries = [
            _compact_result_value(item, max_chars=child_limit)
            for item in value
        ]
        if any(item is None for item in entries):
            return None
        text = '；'.join(entries)
        return text if len(text) <= max_chars else None
    return None


def _result_summary_items(result_data):
    """Select compact, user-facing metrics from a task result summary."""
    if not isinstance(result_data, Mapping):
        return [], 0

    items = []
    omitted = 0
    for key, value in result_data.items():
        if _result_key_is_technical(key):
            omitted += 1
            continue
        compact_value = _compact_result_value(value)
        if compact_value is None or len(items) >= _RESULT_SUMMARY_MAX_ITEMS:
            omitted += 1
            continue
        items.append({'key': str(key), 'value': compact_value})
    return items, omitted


def _review_evidence_for_display(evidence):
    """Copy review evidence with only compact values available to the template."""
    if not isinstance(evidence, Mapping):
        return evidence
    display = dict(evidence)
    checks = []
    for check in evidence.get('checks') or []:
        item = dict(check)
        if 'value' in item:
            compact_value = _compact_result_value(item['value'])
            if compact_value is None:
                item.pop('value', None)
            else:
                item['display_value'] = compact_value
        checks.append(item)
    display['checks'] = checks
    return display


def _validate_path(file_path, project_id):
    if not file_path:
        return False
    is_valid, _ = Config._validate_path(file_path, project_id)
    return is_valid and os.path.isfile(file_path)


def _materialize_static_files(task, project_id, files):
    """Convert legacy Plotly payloads so old tasks also use the static contract."""
    legacy = [f for f in files if f.file_type == 'plotly_json' and _validate_path(f.file_path, project_id)]
    if not legacy:
        return files
    from modules.reporting.static_rendering import render_plotly_json

    known_paths = {f.file_path for f in files}
    for rf in legacy:
        try:
            rendered = render_plotly_json(rf.file_path, label=rf.label or rf.category or '分析图')
            for item in rendered:
                item['category'] = rf.category or 'plot'
                if item['file_path'] in known_paths:
                    continue
                files.append(ResultFile.create(
                    task_id=task.id,
                    project_id=project_id,
                    file_type=item['file_type'],
                    category=item['category'],
                    label=item['label'],
                    file_path=item['file_path'],
                ))
                known_paths.add(item['file_path'])
            try:
                os.remove(rf.file_path)
            except OSError:
                pass
        except Exception:
            continue
    return files


results_bp = Blueprint('results', __name__)

@results_bp.route('/<pid>/task/<task_id>')
def task_detail(pid, task_id):
    p = Project.get_by_id(pid)
    t = AnalysisTask.get_by_id(task_id)
    if not p or not t or t.project_id != pid:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    files = _materialize_static_files(t, pid, ResultFile.get_by_task(task_id))
    for f in files:
        try:
            f._file_size = os.path.getsize(f.file_path) if f.file_path and os.path.exists(f.file_path) else 0
        except Exception:
            f._file_size = 0
    # Plotly JSON is a legacy implementation detail and is never exposed.
    plotly_files = []
    image_files = [f for f in files if f.file_type in ('png', 'svg', 'pdf', 'tiff', 'jpg', 'jpeg')]
    # A Matplotlib figure is emitted as both 300 dpi PNG and SVG.  Show PNG once
    # in the browser (predictable cross-browser rendering), while preserving the
    # SVG companion as a publication-ready download.
    image_groups = {}
    for image in image_files:
        key = (image.category, image.label)
        current = image_groups.get(key)
        if current is None or (image.file_type == 'png' and current.file_type != 'png'):
            image_groups[key] = image
    preferred_image_files = list(image_groups.values())
    for image in preferred_image_files:
        image.svg_variant = next(
            (candidate for candidate in image_files
             if candidate.category == image.category
             and candidate.label == image.label
             and candidate.file_type == 'svg'),
            None,
        )
        image.pdf_variant = next(
            (candidate for candidate in files
             if candidate.category == image.category
             and candidate.label == image.label
             and candidate.file_type == 'pdf'),
            None,
        )
        image.tiff_variant = next(
            (candidate for candidate in files
             if candidate.category == image.category
             and candidate.label == image.label
             and candidate.file_type == 'tiff'),
            None,
        )
    csv_files = [f for f in files if f.file_type == 'csv']
    json_files = [f for f in files if f.file_type == 'json']
    result_data = {}
    try:
        result_data = json.loads(t.result_json) if t.result_json else {}
    except Exception:
        pass
    review_evidence = None
    result_interpretation = None
    try:
        from modules.reporting.review_evidence import build_review_evidence, build_result_interpretation
        review_evidence = build_review_evidence(t.module_name, result_data, files)
        result_interpretation = build_result_interpretation(t.module_name, result_data, files)
    except Exception:
        review_evidence = None

    result_summary_items, omitted_result_fields = _result_summary_items(result_data)
    review_evidence = _review_evidence_for_display(review_evidence)

    # 计算下一步模块
    next_module = None
    current_module = t.module_name
    try:
        idx = PIPELINE_ORDER.index(current_module)
        if idx < len(PIPELINE_ORDER) - 1:
            next_module = PIPELINE_ORDER[idx + 1]
    except ValueError:
        pass

    return render_template('analysis_result.html', project=p, task=t,
                          plotly_files=plotly_files, image_files=preferred_image_files, csv_files=csv_files,
                          json_files=json_files,
                          review_evidence=review_evidence,
                          result_interpretation=result_interpretation,
                          result_summary_items=result_summary_items,
                          omitted_result_fields=omitted_result_fields,
                          module_display=MODULE_DISPLAY_MAP.get(t.module_name, t.module_name),
                          status_cn=STATUS_MAP.get(t.status, t.status),
                          next_module=next_module,
                          module_display_map=MODULE_DISPLAY_MAP)

@results_bp.route('/<pid>/results')
def results_gallery(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    pipeline_runs = PipelineRun.get_by_project(pid)
    files_by_task = {t.id: _materialize_static_files(t, pid, ResultFile.get_by_task(t.id)) for t in tasks}
    report_info = {}
    try:
        from modules.reporting.project_report import ensure_project_report
        report_info = ensure_project_report(
            project=p,
            project_dir=Config.project_dir(pid),
            tasks=tasks,
            files_by_task=files_by_task,
            module_display_map=MODULE_DISPLAY_MAP,
        )
    except Exception:
        report_info = {}
    return render_template('results_gallery.html', project=p, tasks=tasks,
                          pipeline_runs=pipeline_runs,
                          module_display_map=MODULE_DISPLAY_MAP, status_map=STATUS_MAP,
                          report_info=report_info)


@results_bp.route('/<pid>/results/project-report')
def project_report(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    files_by_task = {t.id: _materialize_static_files(t, pid, ResultFile.get_by_task(t.id)) for t in tasks}
    from modules.reporting.project_report import ensure_project_report
    report_info = ensure_project_report(
        project=p,
        project_dir=Config.project_dir(pid),
        tasks=tasks,
        files_by_task=files_by_task,
        module_display_map=MODULE_DISPLAY_MAP,
    )
    report_path = report_info.get('report_path')
    if not _validate_path(report_path, pid):
        flash('报告路径不合法', 'danger')
        return redirect(url_for('main.index'))
    return send_file(report_path, as_attachment=False)


@results_bp.route('/<pid>/results/plot-gallery')
def plot_gallery(pid):
    p = Project.get_by_id(pid)
    if not p:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    tasks = AnalysisTask.get_by_project(pid)
    files_by_task = {t.id: _materialize_static_files(t, pid, ResultFile.get_by_task(t.id)) for t in tasks}
    from modules.reporting.project_report import ensure_project_report
    report_info = ensure_project_report(
        project=p,
        project_dir=Config.project_dir(pid),
        tasks=tasks,
        files_by_task=files_by_task,
        module_display_map=MODULE_DISPLAY_MAP,
    )
    gallery_path = report_info.get('gallery_path')
    if not _validate_path(gallery_path, pid):
        flash('图库路径不合法', 'danger')
        return redirect(url_for('main.index'))
    return send_file(gallery_path, as_attachment=False)


@results_bp.route('/<pid>/results/plots-archive')
def plots_archive(pid):
    """Download the offline gallery and static publication figures as one ZIP."""
    p = Project.get_by_id(pid)
    if not p:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))

    tasks = AnalysisTask.get_by_project(pid)
    files_by_task = {t.id: _materialize_static_files(t, pid, ResultFile.get_by_task(t.id)) for t in tasks}
    from modules.reporting.project_report import ensure_project_report
    report_info = ensure_project_report(
        project=p,
        project_dir=Config.project_dir(pid),
        tasks=tasks,
        files_by_task=files_by_task,
        module_display_map=MODULE_DISPLAY_MAP,
    )

    buffer = io.BytesIO()
    manifest = []
    with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        gallery_path = report_info.get('gallery_path')
        if _validate_path(gallery_path, pid):
            archive.write(gallery_path, 'plot_gallery.html')

        for task in tasks:
            for rf in files_by_task.get(task.id, []):
                if rf.file_type not in {'png', 'svg', 'pdf', 'tiff', 'jpg', 'jpeg'} or not _validate_path(rf.file_path, pid):
                    continue
                safe_label = re.sub(r'[^a-zA-Z0-9._-]+', '_', rf.label or rf.category or 'plot').strip('_') or 'plot'
                extension = os.path.splitext(rf.file_path)[1].lower() or f'.{rf.file_type}'
                arcname = f'static_images/{task.module_name}_{task.id}_{safe_label}_{rf.id}{extension}'
                archive.write(rf.file_path, arcname)
                manifest.append({
                    'task_id': task.id,
                    'module': task.module_name,
                    'label': rf.label,
                    'archive_path': arcname,
                })
        archive.writestr('manifest.json', json.dumps(manifest, ensure_ascii=False, indent=2))

    buffer.seek(0)
    return send_file(
        buffer,
        mimetype='application/zip',
        as_attachment=True,
        download_name=f'{pid}_all_plots.zip',
    )


def _pipeline_artifacts(pid, run):
    task_ids = []
    try:
        task_ids = json.loads(run.task_ids_json) if run.task_ids_json else []
    except Exception:
        task_ids = []
    tasks = []
    for tid in task_ids:
        task = AnalysisTask.get_by_id(tid)
        if task and task.project_id == pid:
            tasks.append(task)
    files_by_task = {t.id: _materialize_static_files(t, pid, ResultFile.get_by_task(t.id)) for t in tasks}
    from modules.reporting.pipeline_report import write_pipeline_report
    return tasks, files_by_task, write_pipeline_report(
        Config.project_dir(pid),
        run,
        tasks,
        files_by_task,
        MODULE_DISPLAY_MAP,
    )


@results_bp.route('/<pid>/pipeline-runs/<run_id>')
def pipeline_run_detail(pid, run_id):
    p = Project.get_by_id(pid)
    run = PipelineRun.get_by_id(run_id)
    if not p or not run:
        flash('未找到', 'danger')
        return redirect(url_for('main.index'))
    if run.project_id != pid:
        flash('流程不属于该项目', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    tasks, files_by_task, artifacts = _pipeline_artifacts(pid, run)
    return render_template(
        'pipeline_run_detail.html',
        project=p,
        run=run,
        run_data=run.to_dict(),
        tasks=tasks,
        files_by_task=files_by_task,
        artifacts=artifacts,
        module_display_map=MODULE_DISPLAY_MAP,
        status_map=STATUS_MAP,
    )


@results_bp.route('/<pid>/pipeline-runs/<run_id>/manifest')
def pipeline_run_manifest(pid, run_id):
    run = PipelineRun.get_by_id(run_id)
    if not run or run.project_id != pid:
        flash('流程未找到', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    _, _, artifacts = _pipeline_artifacts(pid, run)
    manifest_path = artifacts.get('manifest_path')
    if not _validate_path(manifest_path, pid):
        flash('Manifest 路径不合法', 'danger')
        return redirect(url_for('main.index'))
    return send_file(manifest_path, as_attachment=False)


@results_bp.route('/<pid>/pipeline-runs/<run_id>/report')
def pipeline_run_report(pid, run_id):
    run = PipelineRun.get_by_id(run_id)
    if not run or run.project_id != pid:
        flash('流程未找到', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    _, _, artifacts = _pipeline_artifacts(pid, run)
    report_path = artifacts.get('report_path')
    if not _validate_path(report_path, pid):
        flash('报告路径不合法', 'danger')
        return redirect(url_for('main.index'))
    return send_file(report_path, as_attachment=False)

@results_bp.route('/<pid>/results/file/<file_id>')
def view_file(pid, file_id):
    f = ResultFile.get_by_id(file_id)
    if not f or f.project_id != pid:
        flash('文件未找到', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    if f.file_type == 'plotly_json':
        flash('Plotly 交互图已停用，请重新运行任务生成 PNG/SVG。', 'warning')
        return redirect(url_for('projects.detail', pid=pid))
    if not _validate_path(f.file_path, pid):
        flash('文件路径不合法', 'danger')
        return redirect(url_for('main.index'))
    if f.file_type == 'csv':
        return send_file(f.file_path, as_attachment=True)
    return send_file(f.file_path)

@results_bp.route('/<pid>/download/<task_id>')
def download_adata(pid, task_id):
    t = AnalysisTask.get_by_id(task_id)
    if not t or t.project_id != pid or not t.output_adata_path:
        flash('文件未找到', 'danger')
        return redirect(url_for('projects.detail', pid=pid))
    if not _validate_path(t.output_adata_path, pid):
        flash('文件路径不合法', 'danger')
        return redirect(url_for('main.index'))
    return send_file(t.output_adata_path, as_attachment=True)
