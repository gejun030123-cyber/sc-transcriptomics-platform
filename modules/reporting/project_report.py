"""Project-level report and offline static figure gallery generation."""

import base64
import html as html_lib
import json
import os
from datetime import datetime, timezone

from modules.reporting.review_evidence import build_review_evidence


def _safe_json_loads(value, default):
    if not value:
        return default
    try:
        return json.loads(value)
    except (TypeError, ValueError, json.JSONDecodeError):
        return default


def _safe_relpath(path, root):
    if not path:
        return ""
    try:
        return os.path.relpath(path, root)
    except ValueError:
        return path


def _data_uri(path):
    if not path or not os.path.isfile(path):
        return ""
    try:
        with open(path, "rb") as fh:
            encoded = base64.b64encode(fh.read()).decode("ascii")
        extension = os.path.splitext(path)[1].lower()
        mime = "image/svg+xml" if extension == ".svg" else "image/png"
        return f"data:{mime};base64,{encoded}"
    except OSError:
        return ""


def _task_label(task, module_display_map=None):
    module_display_map = module_display_map or {}
    return module_display_map.get(task.module_name, task.module_name)


def write_plot_gallery(project_dir, tasks, files_by_task, module_display_map=None):
    """Write a self-contained HTML gallery for static PNG/SVG result files."""
    html_dir = os.path.join(project_dir, "results", "html_plots")
    os.makedirs(html_dir, exist_ok=True)
    output_path = os.path.join(html_dir, "plot_gallery.html")
    charts = []

    for task in tasks:
        for rf in files_by_task.get(task.id, []):
            if rf.file_type != "png":
                continue
            image_uri = _data_uri(rf.file_path)
            if not image_uri:
                continue
            svg = next((candidate for candidate in files_by_task.get(task.id, [])
                        if candidate.category == rf.category
                        and candidate.label == rf.label
                        and candidate.file_type == "svg"), None)
            charts.append({
                "task": task,
                "file": rf,
                "image_uri": image_uri,
                "svg_uri": _data_uri(svg.file_path) if svg else "",
            })

    parts = [
        "<!doctype html>",
        "<html lang=\"zh-CN\">",
        "<head>",
        "<meta charset=\"utf-8\">",
        "<title>Static Figure Gallery</title>",
        "<style>",
        "body{font-family:Arial,sans-serif;margin:24px;background:#f7f8fb;color:#1f2933}",
        ".chart{background:white;border:1px solid #d9e2ec;border-radius:8px;padding:16px;margin:0 0 20px}",
        ".meta{color:#6b7280;font-size:12px;margin:4px 0 12px}",
        ".plot{width:100%;text-align:center}",
        ".plot img{max-width:100%;height:auto;border-radius:4px}",
        ".chart-actions{display:flex;gap:8px;margin:8px 0 12px}",
        ".chart-actions a{border:1px solid #2563eb;background:white;color:#1d4ed8;border-radius:6px;padding:7px 11px;text-decoration:none}",
        "h1{font-size:24px;margin-bottom:4px} h2{font-size:16px;margin:0}",
        "</style>",
        "</head>",
        "<body>",
        "<h1>项目图表图库</h1>",
        f"<div class=\"meta\">生成时间：{html_lib.escape(datetime.now(timezone.utc).isoformat())}</div>",
    ]
    if not charts:
        parts.append("<p>暂无可展示的静态科研图。</p>")

    for idx, chart in enumerate(charts):
        task = chart["task"]
        rf = chart["file"]
        title = html_lib.escape(rf.label or rf.category or rf.file_type)
        module_name = html_lib.escape(_task_label(task, module_display_map))
        relpath = html_lib.escape(_safe_relpath(rf.file_path, project_dir))
        svg_link = ""
        if chart["svg_uri"]:
            svg_link = f'<a href="{chart["svg_uri"]}" download="{html_lib.escape(rf.label or "plot")}.svg">下载 SVG</a>'
        parts.extend([
            "<section class=\"chart\">",
            f"<h2>{title}</h2>",
            f"<div class=\"meta\">模块：{module_name} | Task：{html_lib.escape(task.id)} | 文件：{relpath}</div>",
            f"<div class=\"chart-actions\"><a href=\"{chart['image_uri']}\" download=\"{html_lib.escape(rf.label or 'plot')}.png\">下载 PNG</a>{svg_link}</div>",
            f"<div class=\"plot\"><img src=\"{chart['image_uri']}\" alt=\"{title}\"></div>",
            "</section>",
        ])

    parts.extend(["</body>", "</html>"])
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(parts))
    return output_path


def _manifest_summary(files):
    for rf in files:
        if rf.file_type == "json" and rf.category == "manifest" and os.path.isfile(rf.file_path):
            try:
                with open(rf.file_path, "r", encoding="utf-8") as fh:
                    return json.load(fh)
            except (OSError, ValueError, json.JSONDecodeError):
                return {}
    return {}


def _review_evidence_from_task(task, summary, files, manifest):
    evidence = manifest.get("review_evidence") if isinstance(manifest, dict) else None
    if evidence:
        return evidence
    return build_review_evidence(task.module_name, summary, files)


def _status_cn(status):
    return {
        "pass": "通过",
        "review": "需复核",
        "warning": "警告",
    }.get(status, status)


def write_project_report(project, project_dir, tasks, files_by_task, module_display_map=None,
                         gallery_path=None):
    """Write a Markdown project report summarizing task results and artifacts."""
    report_dir = os.path.join(project_dir, "results", "reports")
    os.makedirs(report_dir, exist_ok=True)
    output_path = os.path.join(report_dir, "project_report.md")
    module_display_map = module_display_map or {}

    completed = [t for t in tasks if t.status == "completed"]
    failed = [t for t in tasks if t.status == "failed"]
    lines = [
        f"# 项目结果报告：{project.name}",
        "",
        f"- 项目 ID：`{project.id}`",
        f"- 生成时间：{datetime.now(timezone.utc).isoformat()}",
        f"- 任务总数：{len(tasks)}",
        f"- 已完成任务：{len(completed)}",
        f"- 失败任务：{len(failed)}",
    ]
    if gallery_path:
        lines.append(f"- 图表图库：[{_safe_relpath(gallery_path, report_dir)}]({_safe_relpath(gallery_path, report_dir)})")
    lines.extend(["", "## 任务概览", ""])

    if not tasks:
        lines.append("暂无分析任务。")
    for task in tasks:
        files = files_by_task.get(task.id, [])
        summary = _safe_json_loads(task.result_json, {})
        manifest = _manifest_summary(files)
        evidence = _review_evidence_from_task(task, summary, files, manifest)
        display = _task_label(task, module_display_map)
        lines.extend([
            f"### {display}",
            "",
            f"- Task ID：`{task.id}`",
            f"- 状态：`{task.status}`",
            f"- 完成时间：{task.finished_at or ''}",
            f"- 输出 h5ad：`{_safe_relpath(task.output_adata_path, project_dir)}`" if task.output_adata_path else "- 输出 h5ad：无",
            f"- 结果文件数：{len(files)}",
        ])
        if manifest:
            lines.append(f"- Manifest：`{_safe_relpath(next((f.file_path for f in files if f.category == 'manifest'), ''), project_dir)}`")
        if summary:
            lines.append("")
            lines.append("关键 summary：")
            for key in sorted(summary.keys())[:12]:
                val = summary[key]
                if isinstance(val, (dict, list)):
                    val = json.dumps(val, ensure_ascii=False)
                lines.append(f"- `{key}`：{val}")
        if evidence:
            lines.append("")
            lines.append(f"审批证据：{_status_cn(evidence.get('status'))}")
            for check in evidence.get("checks", []):
                value = f"（{check['value']}）" if "value" in check else ""
                lines.append(
                    f"- {_status_cn(check.get('status'))}：{check.get('name')} {value} - {check.get('message')}"
                )
        plot_files = [f for f in files if f.file_type in {"png", "svg", "pdf", "tiff"}]
        table_files = [f for f in files if f.file_type == "csv"]
        if plot_files or table_files:
            lines.append("")
            lines.append("主要产物：")
            for rf in plot_files[:8] + table_files[:8]:
                lines.append(f"- {rf.label or rf.category}：`{_safe_relpath(rf.file_path, project_dir)}`")
        lines.append("")

    lines.extend([
        "## 建议下一步",
        "",
        "- 优先复核 QC、聚类、注释和 DEG 的关键图表与 summary。",
        "- 如果目标细胞类型或候选分支仍不明确，进入 Agent 候选分支对比阶段。",
        "- 对需要交付的项目，保留本报告、manifest、图表图库和 h5ad 输出作为可复现结果包。",
        "",
    ])

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return output_path


def ensure_project_report(project, project_dir, tasks, files_by_task, module_display_map=None):
    gallery_path = write_plot_gallery(project_dir, tasks, files_by_task, module_display_map)
    report_path = write_project_report(
        project=project,
        project_dir=project_dir,
        tasks=tasks,
        files_by_task=files_by_task,
        module_display_map=module_display_map,
        gallery_path=gallery_path,
    )
    return {
        "report_path": report_path,
        "gallery_path": gallery_path,
    }
