"""Typed, checksum-aware storage helpers for WES input/output assets."""

import hashlib
import json
import os
import uuid
from typing import Any, Dict, Iterable, Optional

from config import Config
from database import get_conn


ARTIFACT_KINDS = {
    "sra", "fastq", "bam", "bai", "cram", "crai", "vcf", "gvcf", "tbi",
    "bed", "interval_list", "table", "figure", "report", "manifest", "log",
}


def sha256_file(path: str, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_asset_path(project_id: str, file_path: str,
                        source_roots: Iterable[str] = ()) -> str:
    if not file_path:
        raise ValueError("缺少资产文件路径")
    supplied = str(file_path)
    if not os.path.isabs(supplied):
        supplied = os.path.join(Config.project_dir(project_id), supplied)
    if os.path.islink(supplied):
        raise ValueError("WES 资产不允许符号链接")
    resolved = os.path.realpath(supplied)
    project_root = os.path.realpath(Config.project_dir(project_id))
    roots = [project_root]
    roots.extend(os.path.realpath(root) for root in (source_roots or Config.wes_source_roots()))
    if not any(resolved == root or resolved.startswith(root + os.sep) for root in roots):
        raise ValueError("WES 资产路径不在项目目录或管理员配置的数据目录内")
    if not os.path.isfile(resolved):
        raise ValueError("WES 资产文件不存在或不是普通文件")
    return resolved


def _decode_asset(row):
    data = dict(row)
    for field in ("metadata_json", "parent_asset_ids_json"):
        raw = data.pop(field, "{}" if field == "metadata_json" else "[]")
        try:
            data[field.removesuffix("_json")] = json.loads(raw or ("{}" if field == "metadata_json" else "[]"))
        except (TypeError, ValueError, json.JSONDecodeError):
            data[field.removesuffix("_json")] = {} if field == "metadata_json" else []
    return data


def register_data_asset(project_id: str, *, artifact_kind: str, file_path: str,
                        sample_id: str = "", assay_type: str = "wes", role: str = "",
                        format: str = "", checksum: str = "", size_bytes: Optional[int] = None,
                        reference_build: str = "", metadata: Optional[Dict[str, Any]] = None,
                        parent_asset_ids: Optional[Iterable[str]] = None,
                        created_by_run_id: str = "", source_roots: Iterable[str] = ()) -> Dict[str, Any]:
    kind = str(artifact_kind or "").strip().lower()
    if kind not in ARTIFACT_KINDS:
        raise ValueError(f"不支持的 WES asset artifact_kind: {artifact_kind}")
    resolved = _resolve_asset_path(project_id, file_path, source_roots)
    if size_bytes is None:
        size_bytes = os.path.getsize(resolved)
    checksum = checksum or sha256_file(resolved)
    asset_id = "asset_" + uuid.uuid4().hex[:16]
    metadata = metadata or {}
    parent_asset_ids = list(parent_asset_ids or [])
    conn = get_conn()
    try:
        conn.execute(
            "INSERT INTO data_assets "
            "(id, project_id, sample_id, assay_type, role, artifact_kind, format, file_path, "
            "size_bytes, checksum, reference_build, metadata_json, parent_asset_ids_json, created_by_run_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (asset_id, project_id, sample_id, assay_type, role, kind, format, resolved,
             int(size_bytes), checksum, reference_build, json.dumps(metadata, ensure_ascii=False),
             json.dumps(parent_asset_ids, ensure_ascii=False), created_by_run_id),
        )
        conn.commit()
    finally:
        conn.close()
    return get_data_asset(asset_id, project_id)


def get_data_asset(asset_id: str, project_id: Optional[str] = None):
    conn = get_conn()
    try:
        query = "SELECT * FROM data_assets WHERE id=?"
        params = [asset_id]
        if project_id:
            query += " AND project_id=?"
            params.append(project_id)
        row = conn.execute(query, tuple(params)).fetchone()
    finally:
        conn.close()
    return _decode_asset(row) if row else None


def list_data_assets(project_id: str, *, sample_id: str = "", artifact_kind: str = ""):
    conn = get_conn()
    try:
        query = "SELECT * FROM data_assets WHERE project_id=?"
        params = [project_id]
        if sample_id:
            query += " AND sample_id=?"
            params.append(sample_id)
        if artifact_kind:
            query += " AND artifact_kind=?"
            params.append(artifact_kind)
        query += " ORDER BY created_at, id"
        rows = conn.execute(query, tuple(params)).fetchall()
    finally:
        conn.close()
    return [_decode_asset(row) for row in rows]
