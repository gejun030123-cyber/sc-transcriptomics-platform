"""Small, reviewable launch-bundle helpers for nf-core/sarek.

The platform keeps the user-facing manifest richer than a Nextflow
samplesheet (roles, capture kit, matched normal, reference bundle).  This
module renders only the stable Sarek columns and writes an immutable bundle
next to each prepared run so a lab member can inspect exactly what would be
executed.
"""

import csv
import hashlib
import json
import os
from typing import Any, Dict, Iterable, List, Mapping


SAREK_INPUT_HEADER = [
    "patient", "sex", "status", "sample", "lane",
    "fastq_1", "fastq_2", "bam", "bai", "cram", "crai",
]
SAREK_VCF_HEADER = ["patient", "sample", "vcf"]


def _text(value: Any, default: str = "") -> str:
    value = str(value or "").strip()
    return value or default


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _rows_for_sarek(manifest: Mapping[str, Any], workflow_key: str) -> tuple[List[str], List[Dict[str, str]]]:
    samples = manifest.get("samples") or []
    if workflow_key == "wes_annotate_only":
        header = SAREK_VCF_HEADER
        rows = []
        for sample in samples:
            if _text(sample.get("input_type")).lower() != "vcf":
                raise ValueError("wes_annotate_only 只接受 VCF 样本")
            rows.append({
                "patient": _text(sample.get("patient_id"), _text(sample.get("sample_id"))),
                "sample": _text(sample.get("sample_id")),
                "vcf": _text(sample.get("vcf")),
            })
        return header, rows

    header = SAREK_INPUT_HEADER
    rows = []
    for sample in samples:
        input_type = _text(sample.get("input_type")).lower()
        if input_type not in {"fastq", "bam", "cram"}:
            raise ValueError(f"{sample.get('sample_id', '<unknown>')} 的输入类型不适用于 Sarek: {input_type}")
        role = _text(sample.get("role")).lower()
        rows.append({
            "patient": _text(sample.get("patient_id")),
            "sex": _text(sample.get("sex"), "NA"),
            # Sarek status follows the normal=0 / tumor=1 convention.  A
            # germline sample is represented as 0 and does not become tumor
            # merely because the caller selected the somatic workflow.
            "status": "1" if role == "tumor" else "0",
            "sample": _text(sample.get("sample_id")),
            "lane": _text(sample.get("lane"), "lane_1"),
            "fastq_1": _text(sample.get("fastq_1")) if input_type == "fastq" else "",
            "fastq_2": _text(sample.get("fastq_2")) if input_type == "fastq" else "",
            "bam": _text(sample.get("bam")) if input_type == "bam" else "",
            "bai": _text(sample.get("bai")) if input_type == "bam" else "",
            "cram": _text(sample.get("cram")) if input_type == "cram" else "",
            "crai": _text(sample.get("crai")) if input_type == "cram" else "",
        })
    return header, rows


def render_samplesheet(manifest: Mapping[str, Any], workflow_key: str, output_path: str) -> Dict[str, Any]:
    """Render a deterministic CSV samplesheet and return its metadata."""
    header, rows = _rows_for_sarek(manifest, workflow_key)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return {
        "path": os.path.abspath(output_path),
        "header": header,
        "row_count": len(rows),
        "checksum": _sha256(output_path),
    }


def _write_json(path: str, value: Any) -> str:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(_canonical_json(value))
        handle.write("\n")
    return os.path.abspath(path)


def write_launch_bundle(run_root: str, *, workflow: Mapping[str, Any], run_id: str,
                        project_id: str, manifest_record: Mapping[str, Any],
                        launch: Mapping[str, Any], parameters: Mapping[str, Any] | None = None,
                        profile: str = "docker") -> Dict[str, Any]:
    """Write all reviewable inputs for one prepared WES run.

    The returned checksums cover the generated manifest, samplesheet,
    parameters and provenance files.  Large sequence files are referenced by
    path and are not copied into the bundle.
    """
    launch_dir = os.path.join(os.path.abspath(run_root), "launch")
    os.makedirs(launch_dir, exist_ok=True)
    manifest = manifest_record.get("manifest") or {}
    manifest_path = _write_json(os.path.join(launch_dir, "manifest.json"), manifest)
    samplesheet = render_samplesheet(
        manifest, str(workflow.get("key") or ""), os.path.join(launch_dir, "samplesheet.csv")
    )
    parameters_path = _write_json(os.path.join(launch_dir, "parameters.json"), parameters or {})
    provenance = {
        "project_id": project_id,
        "run_id": run_id,
        "manifest_id": manifest_record.get("id", ""),
        "manifest_checksum": manifest_record.get("checksum", ""),
        "workflow": dict(workflow),
        "profile": profile,
        "launch": dict(launch),
    }
    provenance_path = _write_json(os.path.join(launch_dir, "provenance.json"), provenance)
    files = {
        "manifest": manifest_path,
        "samplesheet": samplesheet["path"],
        "parameters": parameters_path,
        "provenance": provenance_path,
    }
    checksums = {name: _sha256(path) for name, path in files.items()}
    checksums_path = _write_json(os.path.join(launch_dir, "checksums.json"), checksums)
    return {
        "launch_dir": launch_dir,
        "files": files,
        "checksums": checksums,
        "checksums_path": checksums_path,
        "samplesheet": samplesheet,
    }
