"""Pipeline run aggregation reports and resume planning."""

import json
import os
from datetime import datetime, timezone


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


def _load_manifest(files):
    for rf in files or []:
        if getattr(rf, "file_type", "") == "json" and getattr(rf, "category", "") == "manifest":
            path = getattr(rf, "file_path", "")
            if os.path.isfile(path):
                try:
                    with open(path, "r", encoding="utf-8") as fh:
                        return json.load(fh)
                except (OSError, ValueError, json.JSONDecodeError):
                    return {}
    return {}


def build_resume_plan(run, tasks):
    """Return a best-effort plan for resuming a failed pipeline run."""
    modules = _safe_json_loads(run.modules_json, [])
    params = _safe_json_loads(run.params_json, {})
    task_by_module = {task.module_name: task for task in tasks}
    failed_module = run.current_module or ""

    if not failed_module:
        failed_task = next((t for t in tasks if t.status == "failed"), None)
        failed_module = failed_task.module_name if failed_task else ""

    if failed_module and failed_module in modules:
        failed_index = modules.index(failed_module)
    else:
        completed_count = len([t for t in tasks if t.status == "completed"])
        failed_index = min(completed_count, len(modules))
        failed_module = modules[failed_index] if failed_index < len(modules) else ""

    if not failed_module:
        return {
            "can_resume": False,
            "reason": "未找到失败模块或流程已完成。",
            "failed_module": "",
            "remaining_modules": [],
            "resume_input_path": "",
            "params": {},
        }

    if failed_index > 0:
        previous_module = modules[failed_index - 1]
        previous_task = task_by_module.get(previous_module)
        resume_input = previous_task.output_adata_path if previous_task and previous_task.status == "completed" else run.input_path
    else:
        resume_input = run.input_path

    remaining = modules[failed_index:]
    return {
        "can_resume": bool(remaining and resume_input),
        "reason": "从失败模块重新执行后续流程。",
        "failed_module": failed_module,
        "remaining_modules": remaining,
        "resume_input_path": resume_input,
        "params": {m: params.get(m, {}) for m in remaining},
    }


def build_pipeline_manifest(run, tasks, files_by_task):
    modules = _safe_json_loads(run.modules_json, [])
    params = _safe_json_loads(run.params_json, {})
    task_rows = []
    for idx, task in enumerate(tasks):
        files = files_by_task.get(task.id, [])
        manifest = _load_manifest(files)
        task_rows.append({
            "step": idx + 1,
            "task_id": task.id,
            "module_name": task.module_name,
            "status": task.status,
            "started_at": task.started_at,
            "finished_at": task.finished_at,
            "output_adata_path": task.output_adata_path,
            "result_files": [rf.to_dict() for rf in files],
            "summary": _safe_json_loads(task.result_json, {}),
            "manifest": manifest,
            "error_traceback": task.error_traceback,
        })

    return {
        "manifest_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "pipeline_run": run.to_dict(),
        "modules": modules,
        "params": params,
        "tasks": task_rows,
        "resume_plan": build_resume_plan(run, tasks),
    }


def _paths(project_dir, run_id):
    out_dir = os.path.join(project_dir, "results", "pipeline_runs")
    os.makedirs(out_dir, exist_ok=True)
    return (
        os.path.join(out_dir, f"pipeline_run_{run_id}_manifest.json"),
        os.path.join(out_dir, f"pipeline_run_{run_id}_report.md"),
    )


def write_pipeline_report(project_dir, run, tasks, files_by_task, module_display_map=None):
    module_display_map = module_display_map or {}
    manifest_path, report_path = _paths(project_dir, run.id)
    manifest = build_pipeline_manifest(run, tasks, files_by_task)
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=2)

    lines = [
        f"# Pipeline Run 报告：{run.name}",
        "",
        f"- Run ID：`{run.id}`",
        f"- 项目 ID：`{run.project_id}`",
        f"- 分析类型：`{run.analysis_type}`",
        f"- 状态：`{run.status}`",
        f"- 进度：{run.progress}%",
        f"- 输入文件：`{_safe_relpath(run.input_path, project_dir)}`",
        f"- 生成时间：{manifest['generated_at']}",
        "",
        "## 步骤",
        "",
        "| Step | 模块 | 状态 | 输出 | 结果文件数 |",
        "| --- | --- | --- | --- | --- |",
    ]
    for task in manifest["tasks"]:
        module_label = module_display_map.get(task["module_name"], task["module_name"])
        lines.append(
            f"| {task['step']} | {module_label} | {task['status']} | "
            f"`{_safe_relpath(task.get('output_adata_path'), project_dir)}` | "
            f"{len(task.get('result_files') or [])} |"
        )

    resume = manifest["resume_plan"]
    lines.extend(["", "## 失败与续跑建议", ""])
    if run.status == "failed":
        lines.append(f"- 失败模块：`{resume.get('failed_module') or run.current_module}`")
        lines.append(f"- 可续跑：{'是' if resume.get('can_resume') else '否'}")
        lines.append(f"- 续跑输入：`{_safe_relpath(resume.get('resume_input_path'), project_dir)}`")
        lines.append(f"- 续跑模块：`{', '.join(resume.get('remaining_modules') or [])}`")
    else:
        lines.append("- 当前流程未失败，无需续跑。")

    lines.extend(["", "## Manifest", "", f"- `{_safe_relpath(manifest_path, project_dir)}`", ""])
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return {"manifest_path": manifest_path, "report_path": report_path, "manifest": manifest}
