"""Constrained collection of user-facing WES workflow artifacts."""

import hashlib
import json
import os
import uuid
from typing import Any, Dict, Iterable

from database import get_conn

from .runs import get_workflow_run, transition_workflow_run


def _sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _decode(row):
    if not row:
        return None
    item = dict(row)
    try:
        item["metadata"] = json.loads(item.pop("metadata_json") or "{}")
    except (TypeError, ValueError, json.JSONDecodeError):
        item["metadata"] = {}
        item.pop("metadata_json", None)
    item["required"] = bool(item.get("required"))
    return item


def list_workflow_artifacts(run_id: str, project_id: str = ""):
    conn = get_conn()
    try:
        query = (
            "SELECT a.* FROM workflow_artifacts a JOIN workflow_runs r ON r.id=a.run_id "
            "WHERE a.run_id=?"
        )
        params = [run_id]
        if project_id:
            query += " AND r.project_id=?"
            params.append(project_id)
        query += " ORDER BY a.required DESC, a.artifact_kind, a.file_path"
        rows = conn.execute(query, tuple(params)).fetchall()
    finally:
        conn.close()
    return [_decode(row) for row in rows]


def get_workflow_artifact(artifact_id: str, project_id: str = ""):
    conn = get_conn()
    try:
        query = (
            "SELECT a.* FROM workflow_artifacts a JOIN workflow_runs r ON r.id=a.run_id "
            "WHERE a.id=?"
        )
        params = [artifact_id]
        if project_id:
            query += " AND r.project_id=?"
            params.append(project_id)
        row = conn.execute(query, tuple(params)).fetchone()
    finally:
        conn.close()
    return _decode(row)


def _kind(relative_path: str) -> str:
    name = relative_path.lower()
    if name.endswith("multiqc_report.html"):
        return "multiqc_report"
    if name.endswith(".vcf.gz.tbi"):
        return "vcf_index"
    if name.endswith(".ann.vcf.gz"):
        return "annotated_vcf"
    if name.endswith(".filtered.vcf.gz"):
        return "filtered_vcf"
    if name.endswith(".vcf.gz"):
        return "variant_vcf"
    if "/pipeline_info/" in "/" + name or name.startswith("pipeline_info/"):
        return "pipeline_provenance"
    if name.endswith("variantcalled.csv"):
        return "variant_summary"
    return ""


def _candidate_paths(results_dir: str, *, max_files: int = 10000) -> Iterable[tuple[str, str]]:
    seen = 0
    for directory, dirnames, filenames in os.walk(results_dir, followlinks=False):
        dirnames[:] = [name for name in dirnames
                       if not os.path.islink(os.path.join(directory, name))]
        for filename in filenames:
            seen += 1
            if seen > max_files:
                raise ValueError("WES results 文件数超过收集上限")
            path = os.path.join(directory, filename)
            if os.path.islink(path) or not os.path.isfile(path):
                continue
            relative = os.path.relpath(path, results_dir)
            kind = _kind(relative)
            if kind:
                yield path, kind


def _register(run_id: str, path: str, kind: str, required: bool,
              metadata: Dict[str, Any]):
    checksum = _sha256(path)
    conn = get_conn()
    try:
        existing = conn.execute(
            "SELECT * FROM workflow_artifacts WHERE run_id=? AND file_path=? "
            "AND artifact_kind=?",
            (run_id, path, kind),
        ).fetchone()
        if existing:
            if existing["checksum"] != checksum:
                raise ValueError(f"已登记工件 checksum 发生变化: {path}")
            return _decode(existing)
        artifact_id = "wesart_" + uuid.uuid4().hex[:16]
        conn.execute(
            "INSERT INTO workflow_artifacts "
            "(id, run_id, artifact_kind, file_path, checksum, required, metadata_json) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (artifact_id, run_id, kind, path, checksum, int(required),
             json.dumps(metadata, ensure_ascii=False, sort_keys=True)),
        )
        conn.commit()
        row = conn.execute(
            "SELECT * FROM workflow_artifacts WHERE id=?", (artifact_id,)
        ).fetchone()
    finally:
        conn.close()
    return _decode(row)


def collect_workflow_artifacts(run_id: str, project_id: str, *, strict: bool = True):
    """Index known Sarek outputs and enforce the small v1 result contract."""
    run = get_workflow_run(run_id, project_id)
    if not run:
        raise ValueError("WES run 不存在")
    run_root = os.path.realpath(str(run.get("run_dir") or ""))
    results_dir = os.path.realpath(os.path.join(run_root, "results"))
    if not run_root or not results_dir.startswith(run_root + os.sep):
        raise ValueError("WES results 路径越过 run 目录边界")
    if not os.path.isdir(results_dir) or os.path.islink(results_dir):
        errors = ["WES results 目录不存在或无效"]
        if strict and run.get("status") == "completed":
            transition_workflow_run(
                run_id, project_id, "failed", expected_statuses=("completed",),
                error_text="; ".join(errors),
            )
        return {"valid": False, "artifacts": [], "errors": errors}

    workflow_key = run.get("workflow_key")
    primary_kinds = {
        "wes_germline": {"filtered_vcf"},
        "wes_somatic": {"variant_vcf"},
        "wes_annotate_only": {"annotated_vcf"},
    }.get(workflow_key, set())
    artifacts = []
    for path, kind in _candidate_paths(results_dir):
        relative = os.path.relpath(path, results_dir)
        required = kind in primary_kinds or kind in {"vcf_index", "multiqc_report"}
        artifacts.append(_register(
            run_id, path, kind, required,
            {"relative_path": relative, "size_bytes": os.path.getsize(path)},
        ))

    kinds = {item["artifact_kind"] for item in artifacts}
    errors = []
    if primary_kinds and not kinds.intersection(primary_kinds):
        errors.append("缺少 workflow 主 VCF 工件")
    if "vcf_index" not in kinds:
        errors.append("缺少 VCF tabix index")
    if "multiqc_report" not in kinds:
        errors.append("缺少 MultiQC HTML")
    if strict and errors and run.get("status") == "completed":
        transition_workflow_run(
            run_id, project_id, "failed", expected_statuses=("completed",),
            error_text="WES 工件收集失败: " + "; ".join(errors),
        )
    return {"valid": not errors, "artifacts": artifacts, "errors": errors}
