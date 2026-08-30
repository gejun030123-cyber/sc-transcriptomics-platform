"""Persistence helpers for versioned WES manifests and workflow assets."""

import hashlib
import json
import uuid
from typing import Any, Dict, Optional

from database import get_conn


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _manifest_id() -> str:
    return "manifest_" + uuid.uuid4().hex[:16]


def register_manifest(project_id: str, manifest: Dict[str, Any],
                      status: str = "validated", validation: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    content = _canonical_json(manifest)
    checksum = hashlib.sha256(content.encode("utf-8")).hexdigest()
    manifest_id = _manifest_id()
    validation_json = _canonical_json(validation or {})
    conn = get_conn()
    existing_id = None
    try:
        existing = conn.execute(
            "SELECT id FROM sample_manifests WHERE project_id=? AND checksum=?",
            (project_id, checksum),
        ).fetchone()
        if existing:
            existing_id = existing["id"]
            version = None
        else:
            row = conn.execute(
                "SELECT COALESCE(MAX(version), 0) AS version FROM sample_manifests WHERE project_id=?",
                (project_id,),
            ).fetchone()
            version = int(row["version"] or 0) + 1
            conn.execute(
                "INSERT INTO sample_manifests "
                "(id, project_id, version, content_json, checksum, status, validation_json) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (manifest_id, project_id, version, content, checksum, status, validation_json),
            )
            conn.commit()
    finally:
        conn.close()
    if existing_id:
        registered = get_manifest(existing_id, project_id)
        if registered:
            return {
                "id": registered["id"],
                "project_id": registered["project_id"],
                "version": registered["version"],
                "checksum": registered["checksum"],
                "status": registered["status"],
                "manifest": registered["manifest"],
                "validation": registered.get("validation", {}),
            }
    return {
        "id": manifest_id,
        "project_id": project_id,
        "version": version,
        "checksum": checksum,
        "status": status,
        "manifest": manifest,
        "validation": validation or {},
    }


def get_manifest(manifest_id: str, project_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
    conn = get_conn()
    try:
        query = "SELECT * FROM sample_manifests WHERE id=?"
        params = [manifest_id]
        if project_id:
            query += " AND project_id=?"
            params.append(project_id)
        row = conn.execute(query, tuple(params)).fetchone()
    finally:
        conn.close()
    if not row:
        return None
    data = dict(row)
    try:
        data["manifest"] = json.loads(data.pop("content_json") or "{}")
        data["validation"] = json.loads(data.get("validation_json") or "{}")
    except (TypeError, ValueError, json.JSONDecodeError):
        data["manifest"] = {}
        data["validation"] = {}
    data.pop("validation_json", None)
    return data


def list_manifests(project_id: str):
    conn = get_conn()
    try:
        rows = conn.execute(
            "SELECT id, project_id, version, checksum, status, created_at "
            "FROM sample_manifests WHERE project_id=? ORDER BY version DESC",
            (project_id,),
        ).fetchall()
    finally:
        conn.close()
    return [dict(row) for row in rows]
