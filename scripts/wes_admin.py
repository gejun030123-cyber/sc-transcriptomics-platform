#!/usr/bin/env python3
"""Minimal administrator CLI for WES resources and P3 benchmarks."""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from database import init_db
from modules.workflows.references import (
    list_capture_kit_profiles,
    list_reference_assets,
    register_capture_kit_profile,
    register_reference_asset,
    set_reference_asset_status,
)
from modules.workflows.validation import run_benchmark, write_evaluation_bed


def _print(value):
    print(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True))


def _write_vep_manifest(*, cache_root, output, release, archive="", archive_sha256=""):
    cache_root = os.path.realpath(cache_root)
    output = os.path.abspath(output)
    if not os.path.isdir(cache_root):
        raise ValueError("VEP cache 目录不存在")
    if os.path.islink(output):
        raise ValueError("manifest 不允许符号链接")
    if archive:
        archive = os.path.realpath(archive)
        if not os.path.isfile(archive):
            raise ValueError("VEP cache archive 不存在")
        if not archive_sha256:
            digest = hashlib.sha256()
            with open(archive, "rb") as handle:
                for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                    digest.update(chunk)
            archive_sha256 = digest.hexdigest()
    payload = {
        "schema_version": 1,
        "cache_root": cache_root,
        "cache_release": str(release or "").strip(),
        "archive": archive,
        "archive_sha256": archive_sha256,
    }
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return payload | {"manifest_path": output}


def main():
    parser = argparse.ArgumentParser(description="WES 实验室管理员工具")
    sub = parser.add_subparsers(dest="command", required=True)

    capture = sub.add_parser("register-capture", help="登记一个 capture kit profile")
    capture.add_argument("--id", required=True)
    capture.add_argument("--name", required=True)
    capture.add_argument("--version", required=True)
    capture.add_argument("--assembly", default="GRCh38")
    capture.add_argument("--calling-bed", required=True)
    capture.add_argument("--vendor-bed", default="")
    capture.add_argument("--qc-bed", default="")
    capture.add_argument("--padding-bp", type=int, default=0)
    capture.add_argument("--status", choices=("test_only", "validated"), default="test_only")
    capture.add_argument("--source", default="")

    reference = sub.add_parser("register-reference", help="登记一个 reference bundle 文件")
    reference.add_argument("--assembly", default="GRCh38")
    reference.add_argument("--bundle", required=True)
    reference.add_argument("--type", required=True)
    reference.add_argument("--path", required=True)
    reference.add_argument("--status", choices=("registered", "test_only", "validated"), default="registered")
    reference.add_argument("--source", default="")
    reference.add_argument("--metadata-json", default="", help="可选 JSON 文件，写入 asset metadata")

    status = sub.add_parser("set-status", help="提升或停用一个 immutable asset")
    status.add_argument("asset_id")
    status.add_argument("status", choices=("registered", "test_only", "validated", "retired"))

    listing = sub.add_parser("list", help="列出参考资源")
    listing.add_argument("--capture-kits", action="store_true")

    intersect = sub.add_parser("evaluation-bed", help="生成 truth/capture/callable 交集 BED")
    intersect.add_argument("--bed", action="append", required=True)
    intersect.add_argument("--output", required=True)

    benchmark = sub.add_parser("benchmark", help="运行 hap.py/som.py 并记录验收指标")
    benchmark.add_argument("--mode", choices=("germline", "somatic"), required=True)
    benchmark.add_argument("--truth-vcf", required=True)
    benchmark.add_argument("--query-vcf", required=True)
    benchmark.add_argument("--evaluation-bed", required=True)
    benchmark.add_argument("--reference-fasta", required=True)
    benchmark.add_argument("--output-prefix", required=True)
    benchmark.add_argument("--bin", default="")
    benchmark.add_argument("--stratum", default="")
    benchmark.add_argument("--stratum-evidence", default="")

    vep_manifest = sub.add_parser("vep-manifest", help="生成 VEP cache provenance manifest")
    vep_manifest.add_argument("--cache-root", required=True)
    vep_manifest.add_argument("--output", required=True)
    vep_manifest.add_argument("--release", default="")
    vep_manifest.add_argument("--archive", default="")
    vep_manifest.add_argument("--archive-sha256", default="")

    args = parser.parse_args()
    init_db()
    if args.command == "register-capture":
        value = register_capture_kit_profile(
            capture_kit_id=args.id, name=args.name, version=args.version,
            assembly=args.assembly, calling_bed_path=args.calling_bed,
            vendor_bed_path=args.vendor_bed, qc_bed_path=args.qc_bed,
            padding_bp=args.padding_bp, status=args.status, source=args.source,
        )
    elif args.command == "register-reference":
        metadata = {}
        if args.metadata_json:
            metadata = json.loads(Path(args.metadata_json).read_text(encoding="utf-8"))
            if not isinstance(metadata, dict):
                raise ValueError("--metadata-json 必须是 JSON object")
        value = register_reference_asset(
            assembly=args.assembly, bundle_version=args.bundle,
            asset_type=args.type, file_path=args.path, status=args.status,
            source=args.source, metadata=metadata,
        )
    elif args.command == "set-status":
        value = set_reference_asset_status(args.asset_id, args.status)
    elif args.command == "list":
        value = (list_capture_kit_profiles(include_retired=True)
                 if args.capture_kits else list_reference_assets())
    elif args.command == "evaluation-bed":
        value = write_evaluation_bed(args.bed, args.output)
    elif args.command == "vep-manifest":
        value = _write_vep_manifest(
            cache_root=args.cache_root, output=args.output, release=args.release,
            archive=args.archive, archive_sha256=args.archive_sha256,
        )
    else:
        value = run_benchmark(
            mode=args.mode, truth_vcf=args.truth_vcf, query_vcf=args.query_vcf,
            evaluation_bed=args.evaluation_bed, reference_fasta=args.reference_fasta,
            output_prefix=args.output_prefix, benchmark_bin=args.bin,
            stratum=args.stratum, stratum_evidence=args.stratum_evidence,
        )
    _print(value)


if __name__ == "__main__":
    main()
