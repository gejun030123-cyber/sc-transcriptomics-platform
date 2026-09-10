"""Optional content-level checks for WES files.

The platform can run manifest/path checks without bioinformatics packages. If
``pysam`` is installed, this module performs bounded header/index checks; it
never replaces samtools, bcftools or the variant caller.
"""

import gzip
import os
from typing import Any, Dict, List, Optional


def _pysam():
    try:
        import pysam  # type: ignore
        return pysam, None
    except ImportError as exc:
        return None, str(exc)


def inspect_fastq(path: str, *, full_integrity: bool = False) -> Dict[str, Any]:
    result = {"path": path, "available": True, "valid": True, "checks": {}}
    try:
        # Both uploaded plain FASTQ and the platform's .fastq.gz conversion
        # output are accepted.  Do not infer compression from a user-supplied
        # filename; the magic bytes are the actual file contract.
        with open(path, "rb") as probe:
            is_gzip = probe.read(2) == b"\x1f\x8b"
        opener = gzip.open if is_gzip else open
        with opener(path, "rb") as handle:
            payload = handle.read(16 * 1024)
            if full_integrity:
                while handle.read(8 * 1024 * 1024):
                    pass
        lines = payload.splitlines()
        result["checks"]["gzip"] = is_gzip
        result["checks"]["full_integrity"] = bool(full_integrity)
        result["checks"]["sample_lines"] = len(lines)
        if len(lines) < 4:
            result["valid"] = False
            result["error"] = "FASTQ gzip 内容不足一个 read"
        elif not lines[0].startswith(b"@") or not lines[2].startswith(b"+"):
            result["valid"] = False
            result["error"] = "FASTQ 首个 read 格式不正确"
    except (OSError, EOFError) as exc:
        result["valid"] = False
        result["error"] = f"FASTQ gzip 无法读取: {exc}"
    return result


def inspect_alignment(path: str, input_type: str, sample_id: str = "",
                      reference_path: Optional[str] = None) -> Dict[str, Any]:
    pysam, import_error = _pysam()
    result: Dict[str, Any] = {
        "path": path, "input_type": input_type, "available": bool(pysam),
        "valid": True, "checks": {},
    }
    if not pysam:
        result["warning"] = "未安装 pysam，跳过 BAM/CRAM header/index 内容检查"
        result["dependency_error"] = import_error
        return result
    try:
        kwargs = {"reference_filename": reference_path} if reference_path else {}
        with pysam.AlignmentFile(path, "rb", **kwargs) as handle:
            header = handle.header.to_dict()
            sort_order = str(header.get("HD", {}).get("SO", "unknown"))
            read_groups = header.get("RG", []) or []
            sample_names = sorted({str(item.get("SM")) for item in read_groups if item.get("SM")})
            result["checks"].update({
                "sort_order": sort_order,
                "has_sequence_dictionary": bool(handle.references),
                "has_read_groups": bool(read_groups),
                "sample_names": sample_names,
                "has_index": bool(handle.has_index()),
            })
            if sort_order not in {"coordinate", "unknown"}:
                result["valid"] = False
                result["error"] = "alignment 未标记为 coordinate sort"
            if not handle.references:
                result["valid"] = False
                result["error"] = "alignment header 缺少参考序列字典"
            if sample_id and sample_names and sample_id not in sample_names:
                result["valid"] = False
                result["error"] = f"manifest sample_id={sample_id} 不在 BAM/CRAM read group 中"
            if input_type in {"bam", "cram"} and not handle.has_index():
                result["valid"] = False
                result["error"] = f"{input_type.upper()} 缺少可用索引"
    except (OSError, ValueError, RuntimeError) as exc:
        result["valid"] = False
        result["error"] = f"{input_type.upper()} 无法读取: {exc}"
    return result


def inspect_vcf(path: str, sample_id: str = "") -> Dict[str, Any]:
    pysam, import_error = _pysam()
    result: Dict[str, Any] = {
        "path": path, "input_type": "vcf", "available": bool(pysam),
        "valid": True, "checks": {},
    }
    if not pysam:
        result["warning"] = "未安装 pysam，跳过 VCF header/index 内容检查"
        result["dependency_error"] = import_error
        return result
    try:
        with pysam.VariantFile(path) as handle:
            samples = list(handle.header.samples)
            result["checks"].update({
                "contig_count": len(handle.header.contigs),
                "sample_names": samples,
                "has_format": bool(handle.header.formats),
            })
            if sample_id and samples and sample_id not in samples:
                result["valid"] = False
                result["error"] = f"manifest sample_id={sample_id} 不在 VCF sample 列中"
    except (OSError, ValueError, RuntimeError) as exc:
        result["valid"] = False
        result["error"] = f"VCF 无法读取: {exc}"
    return result


def inspect_sample_files(sample: Dict[str, Any], *, reference_path: Optional[str] = None,
                         full_fastq_integrity: bool = False) -> Dict[str, Any]:
    """Run optional checks for one sample and return a structured summary."""
    input_type = str(sample.get("input_type", "") or "").lower()
    sample_id = str(sample.get("sample_id", "") or "")
    checks: List[Dict[str, Any]] = []
    if input_type == "fastq":
        for key in ("fastq_1", "fastq_2"):
            if sample.get(key):
                checks.append(inspect_fastq(
                    str(sample[key]), full_integrity=full_fastq_integrity
                ))
    elif input_type in {"bam", "cram"}:
        checks.append(inspect_alignment(str(sample[input_type]), input_type, sample_id, reference_path))
    elif input_type == "vcf" and sample.get("vcf"):
        checks.append(inspect_vcf(str(sample["vcf"]), sample_id))
    return {
        "sample_id": sample_id,
        "valid": all(item.get("valid", False) for item in checks) if checks else True,
        "checks": checks,
    }
