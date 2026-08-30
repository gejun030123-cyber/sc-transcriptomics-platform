"""Persistence for WES workflow run records and launch provenance."""

import json
import uuid
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, Optional

from database import get_conn


RUN_STATUSES = {
    "pending", "prepared", "running", "completed", "failed", "cancel_requested",
    "cancelled", "interrupted",
}

RUN_TRANSITIONS = {
    "pending": {"prepared", "cancelled"},
    "prepared": {"running", "cancel_requested", "cancelled"},
    "running": {"completed", "failed", "cancel_requested", "interrupted"},
    "cancel_requested": {"cancelled", "failed"},
    "interrupted": {"running", "failed", "cancelled"},
    # A zero exit code is provisional until required result artifacts and
    # checksums have been collected.
    "completed": {"failed"},
    "failed": {"running", "cancelled"},
    "cancelled": {"running"},
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _decode(row):
    if not row:
        return None
    data = dict(row)
    for field in ("launch_json", "provenance_json"):
        key = field.removesuffix("_json")
        try:
            data[key] = json.loads(data.pop(field) or "{}")
        except (TypeError, ValueError, json.JSONDecodeError):
            data[key] = {}
            data.pop(field, None)
    return data


def create_workflow_run(project_id: str, workflow_key: str, manifest_id: str,
                        *, status: str = "prepared", launch: Optional[Dict[str, Any]] = None,
                        provenance: Optional[Dict[str, Any]] = None,
                        parameter_signature: str = "", run_id: str = "",
                        executor: str = "nextflow", workflow_release: str = "",
                        profile: str = "", run_dir: str = "",
                        stdout_path: str = "", stderr_path: str = "") -> Dict[str, Any]:
    if status not in RUN_STATUSES:
        raise ValueError(f"未知 workflow run 状态: {status}")
    run_id = run_id or ("wesrun_" + uuid.uuid4().hex[:16])
    conn = get_conn()
    try:
        conn.execute(
            "INSERT INTO workflow_runs "
            "(id, project_id, workflow_key, manifest_id, status, parameter_signature, "
            "launch_json, provenance_json, executor, workflow_release, profile, run_dir, "
            "stdout_path, stderr_path, updated_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (run_id, project_id, workflow_key, manifest_id, status, parameter_signature,
             json.dumps(launch or {}, ensure_ascii=False),
             json.dumps(provenance or {}, ensure_ascii=False), executor, workflow_release,
             profile, run_dir, stdout_path, stderr_path, _utc_now()),
        )
        conn.commit()
    finally:
        conn.close()
    return get_workflow_run(run_id, project_id)


def get_workflow_run(run_id: str, project_id: Optional[str] = None):
    conn = get_conn()
    try:
        query = "SELECT * FROM workflow_runs WHERE id=?"
        params = [run_id]
        if project_id:
            query += " AND project_id=?"
            params.append(project_id)
        row = conn.execute(query, tuple(params)).fetchone()
    finally:
        conn.close()
    return _decode(row)


def list_workflow_runs(project_id: str):
    conn = get_conn()
    try:
        rows = conn.execute(
            "SELECT * FROM workflow_runs WHERE project_id=? ORDER BY created_at DESC, id DESC",
            (project_id,),
        ).fetchall()
    finally:
        conn.close()
    return [_decode(row) for row in rows]


def update_workflow_run(run_id: str, project_id: str, *, status: Optional[str] = None,
                        external_run_id: Optional[str] = None,
                        error_text: Optional[str] = None,
                        **fields: Any) -> Optional[Dict[str, Any]]:
    """Update a run while keeping lifecycle transitions explicit.

    The additional fields are a small fixed runtime surface used by the local
    Nextflow executor (PID, log paths and exit code).  Arbitrary SQL columns
    are not accepted.
    """
    if status is not None:
        current = get_workflow_run(run_id, project_id)
        if not current:
            return None
        if status != current.get("status") and status not in RUN_TRANSITIONS.get(current.get("status"), set()):
            raise ValueError(f"不允许的 workflow run 状态转换: {current.get('status')} -> {status}")
    updates = []
    values = []
    if status is not None:
        updates.append("status=?")
        values.append(status)
    if external_run_id is not None:
        updates.append("external_run_id=?")
        values.append(external_run_id)
    if error_text is not None:
        updates.append("error_text=?")
        values.append(error_text)
    allowed_fields = {
        "executor", "workflow_release", "profile", "run_dir", "pid",
        "process_group_id", "stdout_path", "stderr_path", "exit_code",
        "cancel_requested_at", "started_at", "finished_at",
    }
    for field, value in fields.items():
        if field not in allowed_fields:
            raise ValueError(f"不支持的 workflow run 字段: {field}")
        updates.append(f"{field}=?")
        values.append(value)
    if updates:
        updates.append("updated_at=?")
        values.append(_utc_now())
        values.extend([run_id, project_id])
        conn = get_conn()
        try:
            conn.execute(
                f"UPDATE workflow_runs SET {', '.join(updates)} WHERE id=? AND project_id=?",
                tuple(values),
            )
            conn.commit()
        finally:
            conn.close()
    return get_workflow_run(run_id, project_id)


def transition_workflow_run(run_id: str, project_id: str, target_status: str,
                            *, expected_statuses: Optional[Iterable[str]] = None,
                            error_text: Optional[str] = None,
                            **fields: Any) -> Optional[Dict[str, Any]]:
    """Atomically move a run through the small lab lifecycle.

    ``expected_statuses`` protects launch/cancel requests from racing with a
    second browser request.  A failed compare-and-swap returns the current
    record so the API can report the actual state without guessing.
    """
    if target_status not in RUN_STATUSES:
        raise ValueError(f"未知 workflow run 状态: {target_status}")
    current = get_workflow_run(run_id, project_id)
    if not current:
        return None
    current_status = current.get("status")
    expected = set(expected_statuses or (current_status,))
    if current_status not in expected:
        return current
    if target_status != current_status and target_status not in RUN_TRANSITIONS.get(current_status, set()):
        raise ValueError(f"不允许的 workflow run 状态转换: {current_status} -> {target_status}")

    updates = ["status=?", "updated_at=?"]
    values = [target_status, _utc_now()]
    if error_text is not None:
        updates.append("error_text=?")
        values.append(error_text)
    allowed_fields = {
        "executor", "workflow_release", "profile", "run_dir", "pid",
        "process_group_id", "stdout_path", "stderr_path", "exit_code",
        "cancel_requested_at", "started_at", "finished_at",
    }
    for field, value in fields.items():
        if field not in allowed_fields:
            raise ValueError(f"不支持的 workflow run 字段: {field}")
        updates.append(f"{field}=?")
        values.append(value)
    placeholders = ",".join("?" for _ in expected)
    values.extend([run_id, project_id, *sorted(expected)])
    conn = get_conn()
    try:
        cursor = conn.execute(
            f"UPDATE workflow_runs SET {', '.join(updates)} "
            f"WHERE id=? AND project_id=? AND status IN ({placeholders})",
            tuple(values),
        )
        conn.commit()
    finally:
        conn.close()
    return get_workflow_run(run_id, project_id)
