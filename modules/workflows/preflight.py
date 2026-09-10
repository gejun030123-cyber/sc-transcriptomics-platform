"""WES manifest validation without reading large sequence files into memory."""

import csv
import json
import os
from typing import Any, Dict, Iterable, List, Optional

from .contracts import PreflightResult
from .registry import get_workflow


ROLE_VALUES = {"germline", "normal", "tumor"}
INPUT_VALUES = {"fastq", "bam", "cram", "vcf"}


_PATH_FIELDS = (
    "fastq_1", "fastq_2", "bam", "bai", "cram", "crai", "vcf", "tbi",
    "pedigree_path", "capture_bed_path",
)


def _resolve_relative_paths(manifest: Dict[str, Any], base_dir: str) -> Dict[str, Any]:
    """Resolve file fields relative to the manifest/project directory."""
    result = dict(manifest)
    normalized = []
    for raw_sample in manifest.get("samples", []):
        if not isinstance(raw_sample, dict):
            normalized.append(raw_sample)
            continue
        sample = dict(raw_sample)
        for field in _PATH_FIELDS:
            value = _as_text(sample.get(field))
            if value and not os.path.isabs(value):
                sample[field] = os.path.abspath(os.path.join(base_dir, value))
        normalized.append(sample)
    result["samples"] = normalized
    return result


def load_manifest_file(path: str) -> Dict[str, Any]:
    """Load a JSON or CSV/TSV manifest from a caller-validated project path."""
    suffix = os.path.basename(path).lower()
    if suffix.endswith(".json"):
        with open(path, "r", encoding="utf-8") as handle:
            value = json.load(handle)
        if isinstance(value, list):
            return _resolve_relative_paths({"samples": value}, os.path.dirname(os.path.abspath(path)))
        if not isinstance(value, dict):
            raise ValueError("manifest JSON 必须是对象或样本数组")
        return _resolve_relative_paths(value, os.path.dirname(os.path.abspath(path)))

    delimiter = "\t" if suffix.endswith((".tsv", ".tab")) else ","
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    if not rows:
        raise ValueError("manifest 文件没有样本行")
    return _resolve_relative_paths({"samples": rows}, os.path.dirname(os.path.abspath(path)))


def _as_text(value: Any) -> str:
    return str(value or "").strip()


def _path_error(path: str, project_dir: Optional[str], label: str,
                require_files: bool, source_roots: Iterable[str] = ()) -> Optional[str]:
    if not path:
        return f"{label} 缺少路径"
    if os.path.islink(path):
        return f"{label} 不允许符号链接"
    if project_dir:
        try:
            real_path = os.path.realpath(path)
            real_project = os.path.realpath(project_dir)
        except (OSError, ValueError) as exc:
            return f"{label} 路径解析失败: {exc}"
        in_project = real_path == real_project or real_path.startswith(real_project + os.sep)
        allowed_roots = []
        for root in source_roots or ():
            try:
                allowed_roots.append(os.path.realpath(root))
            except (OSError, ValueError):
                continue
        in_source_root = any(real_path == root or real_path.startswith(root + os.sep)
                             for root in allowed_roots)
        if not in_project and not in_source_root:
            return f"{label} 不在项目目录或管理员配置的 WES 数据目录内"
    if require_files and not os.path.isfile(path):
        return f"{label} 文件不存在"
    return None


def _detect_input_type(sample: Dict[str, Any]) -> str:
    declared = _as_text(sample.get("input_type")).lower()
    present = []
    if sample.get("fastq_1") or sample.get("fastq_2"):
        present.append("fastq")
    if sample.get("bam"):
        present.append("bam")
    if sample.get("cram"):
        present.append("cram")
    if sample.get("vcf"):
        present.append("vcf")
    if len(present) > 1:
        return "mixed"
    if declared:
        if present and declared != present[0]:
            return "mismatch"
        return declared
    if present:
        return present[0]
    return ""


def validate_manifest(manifest: Any, project_dir: Optional[str] = None,
                      workflow_key: str = "", require_files: bool = True,
                      source_roots: Iterable[str] = (), content_checks: bool = False,
                      reference_path: Optional[str] = None,
                      full_fastq_integrity: bool = False) -> PreflightResult:
    """Validate the platform-facing WES manifest contract.

    The function is deterministic and side-effect free.  Registration and
    workflow launch are intentionally separate operations.
    """
    errors: List[str] = []
    warnings: List[str] = []
    checks: Dict[str, Any] = {}

    if isinstance(manifest, list):
        manifest = {"samples": manifest}
    if not isinstance(manifest, dict):
        return PreflightResult(False, ["manifest 必须是对象或样本数组"])

    samples = manifest.get("samples")
    if not isinstance(samples, list) or not samples:
        return PreflightResult(False, ["manifest.samples 必须是非空数组"])

    spec = get_workflow(workflow_key) if workflow_key else None
    if workflow_key and not spec:
        errors.append(f"未知 WES workflow: {workflow_key}")

    reference_bundle_id = _as_text(manifest.get("reference_bundle_id"))
    capture_bed_id = _as_text(manifest.get("capture_bed_id"))
    if not reference_bundle_id:
        errors.append("缺少 reference_bundle_id")
    if not capture_bed_id:
        errors.append("缺少 capture_bed_id")
    capture_bed_path = _as_text(manifest.get("capture_bed_path"))
    if spec and spec.key != "wes_annotate_only" and not capture_bed_path:
        errors.append("WES calling workflow 缺少 capture_bed_path")
    if capture_bed_path:
        if project_dir and not os.path.isabs(capture_bed_path):
            capture_bed_path = os.path.abspath(os.path.join(project_dir, capture_bed_path))
        path_issue = _path_error(
            capture_bed_path, project_dir, "capture_bed_path", require_files, source_roots
        )
        if path_issue:
            errors.append(path_issue)

    normalized_samples: List[Dict[str, Any]] = []
    sample_ids = set()
    normal_ids = set()
    roles = set()
    content_summaries: List[Dict[str, Any]] = []

    for index, raw_sample in enumerate(samples, start=1):
        if not isinstance(raw_sample, dict):
            errors.append(f"samples[{index}] 必须是对象")
            continue
        sample = {str(key): value for key, value in raw_sample.items()}
        if project_dir:
            sample = _resolve_relative_paths({"samples": [sample]}, project_dir)["samples"][0]
        prefix = f"samples[{index}]"
        sample_id = _as_text(sample.get("sample_id"))
        patient_id = _as_text(sample.get("patient_id"))
        role = _as_text(sample.get("role")).lower()
        input_type = _detect_input_type(sample)
        if not sample_id:
            errors.append(f"{prefix}.sample_id 缺失")
        elif sample_id in sample_ids:
            errors.append(f"sample_id 重复: {sample_id}")
        else:
            sample_ids.add(sample_id)
        if not patient_id:
            errors.append(f"{prefix}.patient_id 缺失")
        if role not in ROLE_VALUES:
            errors.append(f"{prefix}.role 必须是 germline、normal 或 tumor")
        else:
            roles.add(role)
            if role == "normal":
                normal_ids.add(sample_id)
        if input_type not in INPUT_VALUES:
            errors.append(f"{prefix}.input_type 无法识别或不支持: {input_type or '<empty>'}")
        elif spec and input_type not in spec.input_types:
            errors.append(f"{prefix} 的 {input_type} 不适用于 {workflow_key}")

        required_paths = {
            "fastq": (("fastq_1", "R1"), ("fastq_2", "R2")),
            "bam": (("bam", "BAM"), ("bai", "BAI")),
            "cram": (("cram", "CRAM"), ("crai", "CRAI")),
            "vcf": (("vcf", "VCF"), ("tbi", "TBI")),
        }
        for field, label in required_paths.get(input_type, ()):
            path = _as_text(sample.get(field))
            if not path:
                errors.append(f"{prefix}.{field} 缺失")
                continue
            path_issue = _path_error(
                path, project_dir, f"{prefix}.{field}", require_files, source_roots
            )
            if path_issue:
                errors.append(path_issue)

        if role == "tumor":
            matched_normal = _as_text(sample.get("matched_normal_id"))
            if matched_normal:
                checks.setdefault("matched_normals", []).append({
                    "tumor": sample_id,
                    "normal": matched_normal,
                })
            elif workflow_key == "wes_somatic":
                warnings.append(f"{sample_id or prefix} 为 tumor-only；结果证据等级降低")
        if not _as_text(sample.get("capture_kit")):
            warnings.append(f"{sample_id or prefix} 未登记 capture_kit")
        sample_reference = _as_text(sample.get("reference_bundle_id", reference_bundle_id))
        if not sample_reference:
            errors.append(f"{prefix} 缺少 reference_bundle_id")
        elif reference_bundle_id and sample_reference != reference_bundle_id:
            errors.append(f"{prefix}.reference_bundle_id 与 manifest 顶层不一致")
        sample_capture = _as_text(sample.get("capture_bed_id", capture_bed_id))
        if capture_bed_id and sample_capture and sample_capture != capture_bed_id:
            errors.append(f"{prefix}.capture_bed_id 与 manifest 顶层不一致")

        sample["sample_id"] = sample_id
        sample["patient_id"] = patient_id
        sample["role"] = role
        sample["input_type"] = input_type
        if content_checks and input_type in INPUT_VALUES:
            from .content import inspect_sample_files
            content = inspect_sample_files(
                sample, reference_path=reference_path,
                full_fastq_integrity=full_fastq_integrity,
            )
            content_summaries.append(content)
            if not content.get("valid", False):
                errors.append(f"{prefix} 内容检查失败")
            for item in content.get("checks", []):
                if item.get("warning"):
                    warnings.append(f"{prefix}: {item['warning']}")
        normalized_samples.append(sample)

    if workflow_key == "wes_somatic":
        if "germline" in roles:
            errors.append("wes_somatic 只接受 role=normal 或 role=tumor")
        tumor_samples = [row for row in normalized_samples if row.get("role") == "tumor"]
        if not tumor_samples:
            errors.append("somatic workflow 至少需要一个 tumor 样本")
        tumor_only_samples = [
            row for row in tumor_samples if not _as_text(row.get("matched_normal_id"))
        ]
        if tumor_only_samples:
            if manifest.get("tumor_only_confirmed") is not True:
                errors.append("tumor-only workflow 必须显式设置 tumor_only_confirmed=true")
            if not _as_text(manifest.get("tumor_only_reason")):
                errors.append("tumor-only workflow 必须记录 tumor_only_reason")
        for row in tumor_samples:
            normal_id = _as_text(row.get("matched_normal_id"))
            if normal_id and normal_id not in normal_ids:
                errors.append(f"tumor {row.get('sample_id')} 的 matched_normal_id 不存在: {normal_id}")
            if normal_id:
                normal = next((item for item in normalized_samples
                               if item.get("sample_id") == normal_id), None)
                if normal and normal.get("patient_id") != row.get("patient_id"):
                    errors.append(
                        f"tumor {row.get('sample_id')} 与 matched_normal_id={normal_id} 的 patient_id 不一致"
                    )
    elif workflow_key == "wes_germline" and roles - {"germline"}:
        errors.append("wes_germline 只接受 role=germline")

    checks.update({
        "sample_count": len(normalized_samples),
        "roles": sorted(roles),
        "input_types": sorted({row.get("input_type") for row in normalized_samples if row.get("input_type")}),
        "reference_bundle_id": reference_bundle_id,
        "capture_bed_id": capture_bed_id,
    })
    if content_checks:
        checks["content_checks"] = content_summaries
    normalized = dict(manifest)
    normalized["samples"] = normalized_samples
    normalized["reference_bundle_id"] = reference_bundle_id
    normalized["capture_bed_id"] = capture_bed_id
    if capture_bed_path:
        normalized["capture_bed_path"] = capture_bed_path
    return PreflightResult(not errors, errors, warnings, normalized, checks)
