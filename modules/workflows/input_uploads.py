"""Project-scoped WES input uploads used by the browser wizard.

The endpoint deliberately accepts only the files needed to make one typed
sample input complete: FASTQ R1/R2, BAM+BAI, or CRAM+CRAI.  SRA archives and
unpaired alignment files are not accepted as workflow inputs.
"""

import os
import re
import uuid
from typing import Dict, Iterable

from werkzeug.utils import secure_filename

from config import Config
from modules.workflows.assets import register_data_asset


_SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_INPUT_FIELDS = {
    "fastq": ("fastq_1", "fastq_2"),
    "bam": ("bam", "bai"),
    "cram": ("cram", "crai"),
}
_ALLOWED_SUFFIXES = {
    "fastq_1": (".fastq", ".fq", ".fastq.gz", ".fq.gz"),
    "fastq_2": (".fastq", ".fq", ".fastq.gz", ".fq.gz"),
    "bam": (".bam",),
    "bai": (".bai",),
    "cram": (".cram",),
    "crai": (".crai",),
}


def _validate_file(field: str, file_storage):
    if not file_storage or not getattr(file_storage, "filename", ""):
        raise ValueError(f"{field} 文件不能为空")
    original = str(file_storage.filename)
    safe = secure_filename(original)
    if not safe:
        raise ValueError(f"{field} 文件名不合法")
    lower = safe.lower()
    if not any(lower.endswith(suffix) for suffix in _ALLOWED_SUFFIXES[field]):
        allowed = ", ".join(_ALLOWED_SUFFIXES[field])
        raise ValueError(f"{field} 只允许以下后缀: {allowed}")
    return original, safe


def _format_for(field: str, filename: str) -> str:
    lower = filename.lower()
    if field.startswith("fastq"):
        return "fastq.gz" if lower.endswith(".gz") else "fastq"
    return field


def store_wes_input_upload(*, project_id: str, sample_id: str, patient_id: str,
                           role: str, input_type: str, files: Dict[str, object],
                           capture_kit_id: str = "") -> Dict[str, object]:
    """Store and register one complete typed WES sample input."""
    Config._validate_pid(project_id)
    sample_id = str(sample_id or "").strip()
    patient_id = str(patient_id or "").strip()
    role = str(role or "").strip().lower()
    input_type = str(input_type or "").strip().lower()
    if not _SAMPLE_ID_RE.fullmatch(sample_id):
        raise ValueError("sample_id 只能包含字母、数字、点、下划线和短横线（最多 128 个字符）")
    if not patient_id or len(patient_id) > 128:
        raise ValueError("patient_id 不能为空且不能超过 128 个字符")
    if role not in {"germline", "normal", "tumor"}:
        raise ValueError("role 必须是 germline、normal 或 tumor")
    if input_type not in _INPUT_FIELDS:
        raise ValueError("input_type 只支持 fastq、bam 或 cram；SRA 不支持直接上传")

    selected = {}
    for field in _INPUT_FIELDS[input_type]:
        original, safe = _validate_file(field, files.get(field))
        selected[field] = {
            "upload": files[field], "original": original, "safe": safe,
        }

    root = os.path.realpath(Config.project_dir(project_id))
    target_dir = os.path.join(root, "wes_inputs", secure_filename(sample_id), uuid.uuid4().hex[:16])
    os.makedirs(target_dir, mode=0o750, exist_ok=True)
    saved = {}
    try:
        for field, item in selected.items():
            destination = os.path.join(target_dir, item["safe"])
            item["upload"].save(destination)
            if not os.path.isfile(destination) or os.path.getsize(destination) <= 0:
                raise ValueError(f"{field} 上传后为空")
            saved[field] = {
                "path": destination,
                "original": item["original"],
                "size_bytes": os.path.getsize(destination),
            }

        assets = []
        paths = {}
        for field, item in saved.items():
            asset = register_data_asset(
                project_id,
                artifact_kind=field if field in {"bam", "bai", "cram", "crai"} else "fastq",
                file_path=item["path"],
                sample_id=sample_id,
                role=role,
                format=_format_for(field, item["original"]),
                metadata={
                    "upload_kind": "wes_browser_input",
                    "original_filename": item["original"][:255],
                    "input_type": input_type,
                    "input_field": field,
                    "mate": "R1" if field == "fastq_1" else "R2" if field == "fastq_2" else "",
                    "patient_id": patient_id,
                    "capture_kit_id": str(capture_kit_id or "").strip(),
                },
            )
            assets.append(asset)
            paths[field] = item["path"]
    except Exception:
        for item in saved.values():
            try:
                os.remove(item["path"])
            except OSError:
                pass
        try:
            os.rmdir(target_dir)
        except OSError:
            pass
        raise

    return {
        "sample_id": sample_id,
        "patient_id": patient_id,
        "role": role,
        "input_type": input_type,
        "paths": paths,
        "assets": assets,
    }
