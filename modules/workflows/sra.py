"""Asynchronous local SRA-to-FASTQ conversion for the WES browser workflow."""

import gzip
import json
import os
import re
import shutil
import subprocess
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

from werkzeug.utils import secure_filename

from config import Config
from database import get_conn
from modules.workflows.assets import register_data_asset
from modules.workflows.content import inspect_fastq


_SRA_EXECUTOR = ThreadPoolExecutor(max_workers=1)
_FUTURES: Dict[str, Any] = {}
_FUTURES_LOCK = threading.Lock()
_SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_JOB_FIELDS = {
    "status", "progress", "progress_message", "output_paths_json",
    "output_asset_ids_json", "log_text", "error_text", "started_at", "finished_at",
}


def _now():
    return datetime.now(timezone.utc).isoformat()


def _decode(row):
    if not row:
        return None
    item = dict(row)
    for field, default in (
        ("output_paths_json", {}), ("output_asset_ids_json", []),
    ):
        raw = item.pop(field, "")
        try:
            item[field.removesuffix("_json")] = json.loads(raw or json.dumps(default))
        except (TypeError, ValueError, json.JSONDecodeError):
            item[field.removesuffix("_json")] = default
    return item


def get_sra_job(job_id: str, project_id: Optional[str] = None):
    conn = get_conn()
    try:
        query = "SELECT * FROM wes_sra_jobs WHERE id=?"
        params = [job_id]
        if project_id:
            query += " AND project_id=?"
            params.append(project_id)
        row = conn.execute(query, tuple(params)).fetchone()
    finally:
        conn.close()
    return _decode(row)


def list_sra_jobs(project_id: str):
    conn = get_conn()
    try:
        rows = conn.execute(
            "SELECT * FROM wes_sra_jobs WHERE project_id=? ORDER BY created_at DESC, id DESC",
            (project_id,),
        ).fetchall()
    finally:
        conn.close()
    return [_decode(row) for row in rows]


def _update_sra_job(job_id: str, **updates):
    values = {}
    for key, value in updates.items():
        if key not in _JOB_FIELDS:
            raise ValueError(f"不支持的 SRA job 字段: {key}")
        if key in {"output_paths_json", "output_asset_ids_json"}:
            value = json.dumps(value, ensure_ascii=False)
        values[key] = value
    if not values:
        return get_sra_job(job_id)
    assignments = ", ".join(f"{key}=?" for key in values)
    conn = get_conn()
    try:
        conn.execute(
            f"UPDATE wes_sra_jobs SET {assignments} WHERE id=?",
            (*values.values(), job_id),
        )
        conn.commit()
    finally:
        conn.close()
    return get_sra_job(job_id)


def _disk_check(path: str, input_size: int) -> Dict[str, int | bool]:
    free = shutil.disk_usage(path).free
    minimum = int(max(1.0, float(getattr(Config, "WES_SRA_MIN_FREE_GB", Config.WES_MIN_FREE_GB))) * 1024 ** 3)
    factor = max(2.0, float(getattr(Config, "WES_SRA_DISK_FACTOR", 3.0)))
    required = max(minimum, int(max(1, input_size) * factor))
    return {"valid": free >= required, "free_bytes": free, "required_bytes": required}


def create_sra_job(*, project_id: str, sample_id: str, patient_id: str, role: str,
                   capture_kit_id: str,
                   source_path: str, source_asset_id: str, output_dir: str,
                   input_size_bytes: int) -> Dict[str, Any]:
    job_id = "sra_" + uuid.uuid4().hex[:16]
    conn = get_conn()
    try:
        conn.execute(
            "INSERT INTO wes_sra_jobs "
            "(id, project_id, sample_id, patient_id, role, capture_kit_id, source_path, "
            "source_asset_id, output_dir, input_size_bytes) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (job_id, project_id, sample_id, patient_id, role, capture_kit_id, source_path,
             source_asset_id, output_dir, int(input_size_bytes)),
        )
        conn.commit()
    finally:
        conn.close()
    return get_sra_job(job_id, project_id)


def _append_log(job_id: str, lines, message: str, progress: Optional[int] = None):
    lines.append(message.rstrip())
    lines[:] = lines[-200:]
    updates = {"log_text": "\n".join(lines), "progress_message": message.rstrip()[:500]}
    if progress is not None:
        updates["progress"] = progress
    _update_sra_job(job_id, **updates)


def _find_fastq_pair(output_dir: str):
    files = sorted(Path(output_dir).glob("*.fastq"))
    pairs = []
    for path in files:
        name = path.name
        suffix = None
        if name.endswith("_1.fastq"):
            suffix = "_1.fastq"
            mate_name = name[:-len(suffix)] + "_2.fastq"
        elif name.endswith("_R1.fastq"):
            suffix = "_R1.fastq"
            mate_name = name[:-len(suffix)] + "_R2.fastq"
        else:
            continue
        mate = path.with_name(mate_name)
        if mate.is_file():
            pairs.append((path, mate))
    if len(files) != 2 or len(pairs) != 1:
        raise ValueError(
            "SRA 转换后必须恰好得到一对 FASTQ（R1/R2）；"
            f"当前检测到 {len(files)} 个 FASTQ、{len(pairs)} 对"
        )
    return pairs[0]


def _gzip_fastq(path: Path) -> Path:
    target = Path(str(path) + ".gz")
    subprocess.run(["gzip", "-n", "-f", str(path)], check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if not target.is_file() or target.stat().st_size <= 0:
        raise ValueError(f"FASTQ 压缩失败: {path.name}")
    return target


def _fastq_record_count(path: Path) -> int:
    """Count complete FASTQ records without loading the file into memory."""
    records = 0
    line_count = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line_count += 1
            if line_count % 4 == 0:
                records += 1
    if line_count == 0 or line_count % 4:
        raise ValueError(f"FASTQ 记录不完整: {path.name}")
    return records


def _validate_fastq_pair(r1: Path, r2: Path):
    """Validate gzip content and ensure paired files contain equal read counts."""
    checks = []
    for path in (r1, r2):
        check = inspect_fastq(str(path))
        if not check.get("valid"):
            raise ValueError(f"FASTQ 内容检查失败: {path.name}: {check.get('error', '未知错误')}")
        checks.append(_fastq_record_count(path))
    if checks[0] != checks[1]:
        raise ValueError(f"R1/R2 read 数不一致: R1={checks[0]}, R2={checks[1]}")
    return checks[0]


def _run_sra_job(job_id: str):
    job = get_sra_job(job_id)
    if not job:
        return
    lines = []
    _update_sra_job(job_id, status="running", progress=5,
                    progress_message="开始检查磁盘和 SRA 工具", started_at=_now())
    try:
        output_dir = os.path.realpath(job["output_dir"])
        os.makedirs(output_dir, mode=0o750, exist_ok=True)
        disk = _disk_check(output_dir, int(job["input_size_bytes"]))
        if not disk["valid"]:
            raise ValueError(
                f"磁盘空间不足：至少需要 {disk['required_bytes'] / 1024 ** 3:.1f} GB 可用，"
                f"当前 {disk['free_bytes'] / 1024 ** 3:.1f} GB"
            )

        binary_name = str(getattr(Config, "WES_SRA_FASTERQ_BIN", "fasterq-dump") or "fasterq-dump")
        binary = binary_name if os.path.isfile(binary_name) else shutil.which(binary_name)
        if not binary:
            raise ValueError("服务器未找到 fasterq-dump，请安装 SRA Toolkit 或配置 WES_SRA_FASTERQ_BIN")

        temp_dir = os.path.join(output_dir, "tmp")
        os.makedirs(temp_dir, mode=0o750, exist_ok=True)
        threads = max(1, int(getattr(Config, "WES_SRA_THREADS", 4)))
        command = [
            binary, "--split-files", "--threads", str(threads),
            "--temp", temp_dir, "-O", output_dir, job["source_path"],
        ]
        _append_log(job_id, lines, "执行: " + " ".join(command), 10)
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        if process.stdout:
            for line in process.stdout:
                _append_log(job_id, lines, line, 20)
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(f"fasterq-dump 退出码 {return_code}")
        _append_log(job_id, lines, "SRA 转换完成，检查 paired-end FASTQ", 65)

        r1_raw, r2_raw = _find_fastq_pair(output_dir)
        r1_gz = _gzip_fastq(r1_raw)
        r2_gz = _gzip_fastq(r2_raw)
        final_r1 = Path(output_dir) / f"{job['sample_id']}_R1.fastq.gz"
        final_r2 = Path(output_dir) / f"{job['sample_id']}_R2.fastq.gz"
        os.replace(r1_gz, final_r1)
        os.replace(r2_gz, final_r2)
        read_count = _validate_fastq_pair(final_r1, final_r2)
        _append_log(
            job_id, lines,
            f"R1/R2 配对和 gzip 内容检查通过（{read_count} reads），登记 FASTQ 数据资产",
            85,
        )

        assets = []
        paths = {"fastq_1": str(final_r1), "fastq_2": str(final_r2)}
        for mate, path in (("R1", final_r1), ("R2", final_r2)):
            asset = register_data_asset(
                job["project_id"], artifact_kind="fastq", file_path=str(path),
                sample_id=job["sample_id"], role=job["role"], format="fastq.gz",
                parent_asset_ids=[job["source_asset_id"]] if job["source_asset_id"] else [],
                metadata={
                    "upload_kind": "sra_conversion",
                    "sra_job_id": job_id,
                    "mate": mate,
                    "patient_id": job["patient_id"],
                    "capture_kit_id": job["capture_kit_id"],
                    "source_sra": job["source_path"],
                },
            )
            assets.append(asset)
        # The temp directory can contain large uncompressed intermediates.  Keep
        # the original SRA and final gzip files for provenance, but reclaim the
        # disposable conversion workspace after successful registration.
        shutil.rmtree(temp_dir, ignore_errors=True)

        _update_sra_job(
            job_id, status="completed", progress=100, progress_message="转换完成",
            output_paths_json=paths, output_asset_ids_json=[item["id"] for item in assets],
            log_text="\n".join(lines), finished_at=_now(),
        )
    except Exception as exc:
        _append_log(job_id, lines, f"失败: {exc}", 100)
        _update_sra_job(job_id, status="failed", error_text=str(exc), finished_at=_now())
    finally:
        with _FUTURES_LOCK:
            _FUTURES.pop(job_id, None)


def submit_sra_job(job_id: str) -> bool:
    with _FUTURES_LOCK:
        if job_id in _FUTURES and not _FUTURES[job_id].done():
            return False
        future = _SRA_EXECUTOR.submit(_run_sra_job, job_id)
        _FUTURES[job_id] = future
    return True


def store_sra_upload(*, project_id: str, file_storage, sample_id: str,
                     patient_id: str, role: str, capture_kit_id: str = "") -> Dict[str, Any]:
    Config._validate_pid(project_id)
    sample_id = str(sample_id or "").strip()
    patient_id = str(patient_id or "").strip()
    role = str(role or "").strip().lower()
    if not _SAMPLE_ID_RE.fullmatch(sample_id):
        raise ValueError("sample_id 只能包含字母、数字、点、下划线和短横线（最多 128 个字符）")
    if not patient_id or len(patient_id) > 128:
        raise ValueError("patient_id 不能为空且不能超过 128 个字符")
    if role not in {"germline", "normal", "tumor"}:
        raise ValueError("role 必须是 germline、normal 或 tumor")
    if not file_storage or not getattr(file_storage, "filename", ""):
        raise ValueError("请选择 SRA 文件")
    original = str(file_storage.filename)
    safe = os.path.basename(original)
    if not safe.lower().endswith(".sra"):
        raise ValueError("SRA 上传只允许 .sra 文件；SRR accession 请先下载为本地 SRA")
    safe_name = secure_filename(safe)
    if not safe_name:
        raise ValueError("SRA 文件名无效")

    job_id = "sra_" + uuid.uuid4().hex[:16]
    project_root = os.path.realpath(Config.project_dir(project_id))
    job_dir = os.path.join(project_root, "wes_sra", job_id)
    output_dir = os.path.join(job_dir, "fastq")
    os.makedirs(output_dir, mode=0o750, exist_ok=True)
    source_path = os.path.join(job_dir, safe_name)
    file_storage.save(source_path)
    try:
        size = os.path.getsize(source_path)
        if size <= 0:
            raise ValueError("SRA 文件为空")
        disk = _disk_check(job_dir, size)
        if not disk["valid"]:
            raise ValueError(
                f"磁盘空间不足：转换前至少需要 {disk['required_bytes'] / 1024 ** 3:.1f} GB 可用，"
                f"当前 {disk['free_bytes'] / 1024 ** 3:.1f} GB"
            )
        source_asset = register_data_asset(
            project_id, artifact_kind="sra", file_path=source_path,
            sample_id=sample_id, role=role, format="sra",
            metadata={
                "upload_kind": "wes_browser_sra",
                "original_filename": original[:255],
                "patient_id": patient_id,
                "capture_kit_id": str(capture_kit_id or "").strip(),
            },
        )
        job = create_sra_job(
            project_id=project_id, sample_id=sample_id, patient_id=patient_id,
            role=role, capture_kit_id=str(capture_kit_id or "").strip(),
            source_path=source_path, source_asset_id=source_asset["id"],
            output_dir=output_dir, input_size_bytes=size,
        )
        submit_sra_job(job["id"])
        return {"job": get_sra_job(job["id"], project_id), "source_asset": source_asset}
    except Exception:
        try:
            shutil.rmtree(job_dir)
        except OSError:
            pass
        raise
