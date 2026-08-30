"""Controlled capture BED uploads for the internal WES catalog."""

import hashlib
import os
import re
import uuid

from werkzeug.utils import secure_filename

from config import Config
from modules.workflows.references import (
    get_capture_kit_profile,
    register_capture_kit_profile,
)


_KIT_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
_ALLOWED_SUFFIXES = (".bed", ".bed.gz")


def _checksum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def store_capture_bed_upload(*, project_id, file_storage, capture_kit_id,
                             name, version, assembly="GRCh38"):
    """Save and register one user-provided calling BED as ``test_only``.

    The upload is deliberately not promoted to ``validated`` and never starts
    a workflow.  An administrator must review the profile and promote its
    immutable asset later.
    """
    Config._validate_pid(project_id)
    kit_id = str(capture_kit_id or "").strip()
    display_name = str(name or "").strip()
    kit_version = str(version or "").strip()
    assembly = str(assembly or "GRCh38").strip()
    if not _KIT_ID_RE.fullmatch(kit_id):
        raise ValueError("capture_kit_id 只能包含字母、数字、下划线和短横线（最多 64 个字符）")
    if not display_name or len(display_name) > 200:
        raise ValueError("capture kit 名称不能为空且不能超过 200 个字符")
    if not kit_version or len(kit_version) > 100:
        raise ValueError("capture kit 版本不能为空且不能超过 100 个字符")
    if not assembly or len(assembly) > 50:
        raise ValueError("assembly 无效")
    if not file_storage or not getattr(file_storage, "filename", ""):
        raise ValueError("请选择 BED 或 BED.GZ 文件")

    original_name = str(file_storage.filename)
    safe_name = secure_filename(original_name)
    lower_name = safe_name.lower()
    if not safe_name or not lower_name.endswith(_ALLOWED_SUFFIXES):
        raise ValueError("只允许上传 .bed 或 .bed.gz 文件")

    existing = get_capture_kit_profile(kit_id)
    if existing and str(existing.get("version") or "") != kit_version:
        raise ValueError("capture_kit_id 已存在但版本不同，请为新版本使用新的稳定 ID")

    upload_root = Config.wes_upload_root()
    target_dir = os.path.join(upload_root, "projects", project_id, "capture", kit_id)
    os.makedirs(target_dir, mode=0o750, exist_ok=True)
    target_path = os.path.join(target_dir, f"{uuid.uuid4().hex[:16]}_{safe_name}")
    file_storage.save(target_path)
    try:
        size = os.path.getsize(target_path)
        max_bytes = int(max(1.0, Config.WES_MAX_CAPTURE_BED_MB) * 1024 * 1024)
        if size > max_bytes:
            raise ValueError(f"BED 文件超过 {Config.WES_MAX_CAPTURE_BED_MB:g} MB 限制")
        profile = register_capture_kit_profile(
            capture_kit_id=kit_id,
            name=display_name,
            version=kit_version,
            assembly=assembly,
            calling_bed_path=target_path,
            status="test_only",
            source=f"user_upload:project:{project_id}",
            metadata={
                "upload_kind": "user_capture_bed",
                "uploaded_by_project": project_id,
                "original_filename": original_name[:255],
            },
            source_roots=(upload_root,),
        )
        return {
            "profile": profile,
            "path": target_path,
            "checksum": _checksum(target_path),
            "size_bytes": size,
        }
    except Exception:
        try:
            os.remove(target_path)
        except OSError:
            pass
        raise
