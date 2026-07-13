"""Utilities for writing structured analysis result manifests."""

import json
import os
import platform
from datetime import datetime, timezone

from modules.reporting.review_evidence import build_review_evidence


MANIFEST_VERSION = 1


def _safe_json_loads(value, default):
    if not value:
        return default
    try:
        return json.loads(value)
    except (TypeError, ValueError, json.JSONDecodeError):
        return default


def _file_entry(path):
    entry = {
        "path": path,
        "exists": False,
        "size_bytes": None,
    }
    if path and os.path.exists(path):
        entry["exists"] = True
        try:
            entry["size_bytes"] = os.path.getsize(path)
        except OSError:
            entry["size_bytes"] = None
    return entry


def build_task_manifest(task, result, result_files=None, pipeline_run_id=None):
    """Build a serializable manifest for a completed analysis task."""
    params = _safe_json_loads(getattr(task, "params_json", "{}"), {})
    summary = result.get("summary", {}) if isinstance(result, dict) else {}
    output_adata = result.get("output_adata") if isinstance(result, dict) else None
    raw_files = result_files if result_files is not None else result.get("result_files", [])
    result_file_entries = []

    for rf in raw_files or []:
        if hasattr(rf, "to_dict"):
            item = rf.to_dict()
        else:
            item = dict(rf)
        path = item.get("file_path", "")
        item["file"] = _file_entry(path)
        result_file_entries.append(item)

    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "project_id": task.project_id,
        "task_id": task.id,
        "pipeline_run_id": pipeline_run_id,
        "branch_id": getattr(task, "branch_id", None),
        "module_name": task.module_name,
        "status": "completed",
        "params": params,
        "summary": summary,
        "output_adata": _file_entry(output_adata),
        "result_files": result_file_entries,
        "runtime": {
            "python_version": platform.python_version(),
            "platform": platform.platform(),
        },
    }
    review_evidence = build_review_evidence(task.module_name, summary, raw_files)
    if review_evidence:
        manifest["review_evidence"] = review_evidence
    return manifest


def write_task_manifest(project_dir, task, result, result_files=None, pipeline_run_id=None):
    """Write a task manifest under results/manifests and return a ResultFile-like dict."""
    manifests_dir = os.path.join(project_dir, "results", "manifests")
    os.makedirs(manifests_dir, exist_ok=True)
    manifest = build_task_manifest(
        task=task,
        result=result,
        result_files=result_files,
        pipeline_run_id=pipeline_run_id,
    )
    manifest_path = os.path.join(
        manifests_dir,
        f"{task.module_name}_{task.id}_analysis_manifest.json",
    )
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=2)
    return {
        "file_path": manifest_path,
        "file_type": "json",
        "category": "manifest",
        "label": f"{task.module_name} analysis manifest",
    }
