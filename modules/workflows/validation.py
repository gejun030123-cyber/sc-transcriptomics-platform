"""P3 small-variant benchmark helpers for GIAB and HCC1395.

The platform uses hap.py/som.py for scientific comparison.  Python here only
prepares the evaluation BED, fixes the command contract, parses metrics and
writes provenance; it does not replace haplotype-aware benchmarking.
"""

import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

from config import Config


THRESHOLDS = {
    "germline": {
        "SNP": {"precision": 0.99, "recall": 0.97},
        "INDEL": {"precision": 0.98, "recall": 0.90},
    },
    "somatic": {
        "SNP": {"precision": 0.90, "recall": 0.90},
        "INDEL": {"precision": 0.80, "recall": 0.80},
    },
}


def _sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_bed(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if path.lower().endswith(".gz") else open(
        path, "rt", encoding="utf-8"
    )


def _read_merged_bed(path: str) -> tuple[Dict[str, List[Tuple[int, int]]], List[str]]:
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    contig_order = []
    with _open_bed(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number} BED 少于 3 列")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number} BED 坐标不是整数") from exc
            contig = fields[0]
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number} BED 区间无效")
            if contig not in intervals:
                intervals[contig] = []
                contig_order.append(contig)
            intervals[contig].append((start, end))
    if not intervals:
        raise ValueError(f"BED 不包含有效区间: {path}")
    for contig, values in intervals.items():
        merged = []
        for start, end in sorted(values):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        intervals[contig] = merged
    return intervals, contig_order


def _intersect_two(left: Mapping[str, Sequence[Tuple[int, int]]],
                   right: Mapping[str, Sequence[Tuple[int, int]]]):
    result: Dict[str, List[Tuple[int, int]]] = {}
    for contig, left_values in left.items():
        right_values = right.get(contig, ())
        i = j = 0
        overlaps = []
        while i < len(left_values) and j < len(right_values):
            left_start, left_end = left_values[i]
            right_start, right_end = right_values[j]
            start, end = max(left_start, right_start), min(left_end, right_end)
            if start < end:
                overlaps.append((start, end))
            if left_end <= right_end:
                i += 1
            else:
                j += 1
        if overlaps:
            result[contig] = overlaps
    return result


def write_evaluation_bed(bed_paths: Sequence[str], output_path: str) -> Dict[str, Any]:
    """Write truth-confidence ∩ capture ∩ callable intervals."""
    if len(bed_paths) < 2:
        raise ValueError("科学 benchmark 至少需要 truth confident BED 和 capture BED")
    merged, order = _read_merged_bed(os.path.abspath(bed_paths[0]))
    input_summaries = []
    for path in bed_paths:
        path = os.path.abspath(path)
        if not os.path.isfile(path) or os.path.islink(path):
            raise ValueError(f"benchmark BED 不存在或为符号链接: {path}")
        values, _ = _read_merged_bed(path)
        styles = {"chr" if contig.startswith("chr") else "bare" for contig in values}
        if len(styles) > 1:
            raise ValueError(f"benchmark BED 混用 contig 风格: {path}")
        input_summaries.append({"path": path, "checksum": _sha256(path)})
    for path in bed_paths[1:]:
        right, _ = _read_merged_bed(os.path.abspath(path))
        merged = _intersect_two(merged, right)
    interval_count = sum(len(values) for values in merged.values())
    bases = sum(end - start for values in merged.values() for start, end in values)
    if not interval_count:
        raise ValueError("truth/capture/callable BED 交集为空；可能是 assembly 或 contig 命名不一致")
    output_path = os.path.abspath(output_path)
    if os.path.exists(output_path):
        raise ValueError("evaluation BED 已存在；请使用新的版本化输出路径")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        for contig in order:
            for start, end in merged.get(contig, ()):
                handle.write(f"{contig}\t{start}\t{end}\n")
    return {
        "path": output_path,
        "checksum": _sha256(output_path),
        "interval_count": interval_count,
        "bases": bases,
        "inputs": input_summaries,
    }


def build_benchmark_command(*, mode: str, truth_vcf: str, query_vcf: str,
                            evaluation_bed: str, reference_fasta: str,
                            output_prefix: str, benchmark_bin: str = "") -> List[str]:
    if mode not in {"germline", "somatic"}:
        raise ValueError("benchmark mode 必须是 germline 或 somatic")
    binary = benchmark_bin or (
        Config.WES_HAPPY_BIN if mode == "germline" else Config.WES_SOMPY_BIN
    )
    return [
        binary, os.path.abspath(truth_vcf), os.path.abspath(query_vcf),
        "-f", os.path.abspath(evaluation_bed),
        "-r", os.path.abspath(reference_fasta),
        "-o", os.path.abspath(output_prefix),
        "--pass-only",
    ]


def parse_benchmark_summary(summary_path: str) -> Dict[str, Dict[str, float]]:
    metrics: Dict[str, Dict[str, float]] = {}
    with open(summary_path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            variant_type = str(row.get("Type") or row.get("TYPE") or "").upper()
            filter_name = str(row.get("Filter") or row.get("FILTER") or "").upper()
            if variant_type not in {"SNP", "INDEL"} or filter_name not in {"PASS", "ALL"}:
                continue
            # PASS takes precedence; ALL is retained only when a tool does not
            # emit a separate PASS row despite --pass-only.
            if variant_type in metrics and filter_name != "PASS":
                continue
            def number(*names):
                for name in names:
                    value = row.get(name)
                    if value not in (None, ""):
                        return float(value)
                return 0.0
            metrics[variant_type] = {
                "precision": number("METRIC.Precision", "Precision", "precision"),
                "recall": number("METRIC.Recall", "Recall", "recall"),
                "f1": number("METRIC.F1_Score", "F1", "f1"),
                "truth_tp": number("TRUTH.TP"),
                "truth_fn": number("TRUTH.FN"),
                "query_fp": number("QUERY.FP"),
                "filter": filter_name,
            }
    if not metrics:
        raise ValueError("benchmark summary 未找到 SNP/INDEL 指标")
    return metrics


def evaluate_metrics(metrics: Mapping[str, Mapping[str, float]], *, mode: str,
                     stratum: str = "", stratum_evidence_valid: bool = False) -> Dict[str, Any]:
    if mode not in THRESHOLDS:
        raise ValueError("未知 benchmark mode")
    checks = []
    for variant_type, limits in THRESHOLDS[mode].items():
        observed = metrics.get(variant_type)
        if not observed:
            checks.append({"variant_type": variant_type, "passed": False,
                           "reason": "缺少指标", "thresholds": limits})
            continue
        passed = (observed.get("precision", 0) >= limits["precision"] and
                  observed.get("recall", 0) >= limits["recall"])
        checks.append({"variant_type": variant_type, "passed": passed,
                       "observed": dict(observed), "thresholds": limits})
    scientific_scope_ok = (mode != "somatic" or
                           (stratum == "vaf_ge_0.10" and stratum_evidence_valid))
    return {
        "passed": all(item["passed"] for item in checks) and scientific_scope_ok,
        "scientific_scope_ok": scientific_scope_ok,
        "status": ("validated" if all(item["passed"] for item in checks) and scientific_scope_ok
                   else "not_validated"),
        "mode": mode,
        "stratum": stratum,
        "checks": checks,
        "limitations": ([] if scientific_scope_ok else [
            "somatic 指标缺少与当前 truth/query checksum 绑定的 VAF≥0.10 分层证据，不能用于生产验收"
        ]),
    }


def run_benchmark(*, mode: str, truth_vcf: str, query_vcf: str,
                  evaluation_bed: str, reference_fasta: str,
                  output_prefix: str, benchmark_bin: str = "",
                  stratum: str = "", stratum_evidence: str = "") -> Dict[str, Any]:
    inputs = [truth_vcf, query_vcf, evaluation_bed, reference_fasta]
    for path in inputs:
        if not os.path.isfile(path) or os.path.islink(path):
            raise ValueError(f"benchmark 输入不存在或为符号链接: {path}")
    command = build_benchmark_command(
        mode=mode, truth_vcf=truth_vcf, query_vcf=query_vcf,
        evaluation_bed=evaluation_bed, reference_fasta=reference_fasta,
        output_prefix=output_prefix, benchmark_bin=benchmark_bin,
    )
    binary = command[0]
    resolved = os.path.abspath(binary) if os.path.sep in binary else shutil.which(binary)
    if not resolved or not os.path.isfile(resolved) or not os.access(resolved, os.X_OK):
        raise ValueError(f"benchmark 工具不可用: {binary}")
    command[0] = resolved
    output_prefix = os.path.abspath(output_prefix)
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    for suffix in (".validation.json", ".summary.csv"):
        if os.path.exists(output_prefix + suffix):
            raise ValueError("benchmark 输出已存在；请使用新的版本化 output_prefix")
    evidence_record = None
    evidence_valid = False
    if stratum_evidence:
        evidence_path = os.path.abspath(stratum_evidence)
        if not os.path.isfile(evidence_path) or os.path.islink(evidence_path):
            raise ValueError("somatic stratum evidence 不存在或为符号链接")
        with open(evidence_path, encoding="utf-8") as handle:
            evidence_record = json.load(handle)
        evidence_valid = bool(
            isinstance(evidence_record, dict) and
            evidence_record.get("stratum") == stratum and
            evidence_record.get("truth_vcf_checksum") == _sha256(truth_vcf) and
            evidence_record.get("query_vcf_checksum") == _sha256(query_vcf) and
            str(evidence_record.get("filter_expression") or "").strip()
        )
    completed = subprocess.run(
        command, cwd=os.path.dirname(output_prefix), stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=False,
    )
    summary_path = output_prefix + ".summary.csv"
    metrics = parse_benchmark_summary(summary_path) if completed.returncode == 0 else {}
    evaluation = (evaluate_metrics(
        metrics, mode=mode, stratum=stratum, stratum_evidence_valid=evidence_valid
    )
                  if metrics else {"passed": False, "status": "not_validated",
                                   "checks": [], "limitations": ["benchmark 执行失败"]})
    record = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "mode": mode,
        "command": command,
        "exit_code": completed.returncode,
        "inputs": {os.path.basename(path): {"path": os.path.abspath(path), "checksum": _sha256(path)}
                   for path in inputs},
        "summary_path": summary_path if os.path.isfile(summary_path) else "",
        "metrics": metrics,
        "evaluation": evaluation,
        "stratum_evidence": ({
            "path": os.path.abspath(stratum_evidence),
            "checksum": _sha256(stratum_evidence),
            "valid": evidence_valid,
            "record": evidence_record,
        } if stratum_evidence else None),
        "stdout_tail": completed.stdout[-8000:],
        "stderr_tail": completed.stderr[-8000:],
    }
    record_path = output_prefix + ".validation.json"
    with open(record_path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, ensure_ascii=False, sort_keys=True, indent=2)
        handle.write("\n")
    record["record_path"] = record_path
    return record
