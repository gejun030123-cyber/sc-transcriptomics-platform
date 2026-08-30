#!/usr/bin/env python3
"""Backfill sample_id / condition metadata into a project's h5ad files.

Multi-batch ZIP imports (10x_mtx_zip_batches) historically wrote only
obs["batch"] / obs["batch_source"].  This script adds:

* obs["sample_id"]  — taken from the existing batch column (or, for the
  per-batch files under uploads/batch_imports/N_<name>/, from the folder
  name), so every sample is addressable as an independent biological sample.
* obs["condition"]  — from a user-supplied batch -> condition mapping
  (e.g. {"batch_1": "疾病组", "NEC_251015": "对照组"}).

Existing non-empty values are never overwritten; a JSON audit of every changed
file is written next to the data.  The script must be run where the data
directory is writable (the platform server), e.g.:

    python scripts/backfill_sample_conditions.py \
        --project-dir /home/oelab/AnaData/GJ/02_platform_runtime/data/projects/9af1c57a-0da \
        --mapping '{"batch_1": "疾病组", "batch_2": "疾病组", "IP3": "疾病组", "NEC_251015": "对照组"}'
"""

import argparse
import json
import os
import re
import sys
import time

import anndata as ad
import pandas as pd


def _safe_name(value, fallback):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("._-")
    return name or fallback


def _batch_from_directory(path):
    """uploads/batch_imports/2_batch_2/batch_imported.h5ad -> batch_2"""
    parent = os.path.basename(os.path.dirname(path))
    parts = parent.split("_", 1)
    return parts[1] if len(parts) > 1 and parts[1] else parent


def _batch_series(adata, path):
    """Return the per-cell batch label and how it was derived."""
    for column in ("batch", "batch_source"):
        if column in adata.obs.columns:
            return adata.obs[column].astype(str), column
    inferred = _batch_from_directory(path)
    return pd.Series(inferred, index=adata.obs.index), "directory"


def annotate_h5ad(path, mapping, dry_run=False):
    """Add sample_id/condition to one h5ad; returns a report dict."""
    path = os.path.abspath(path)
    a = ad.read_h5ad(path)
    report = {
        "path": path,
        "n_obs": int(a.n_obs),
        "n_vars": int(a.n_vars),
        "batch_source_column": None,
        "sample_id_added": False,
        "condition_added": False,
        "condition_filled": 0,
        "batch_values": {},
        "mapping_applied": {},
    }

    batch, source_column = _batch_series(a, path)
    report["batch_source_column"] = source_column
    report["batch_values"] = {str(k): int(v) for k, v in batch.value_counts().items()}
    report["mapping_applied"] = {
        str(k): str(mapping.get(k, ""))
        for k in sorted(set(batch.unique().tolist()) | set(mapping))
    }

    changed = False
    if "sample_id" not in a.obs.columns:
        a.obs["sample_id"] = batch.map(lambda value: _safe_name(value, "sample"))
        report["sample_id_added"] = True
        changed = True

    if "condition" not in a.obs.columns:
        a.obs["condition"] = batch.map(lambda value: str(mapping.get(value, "")))
        report["condition_added"] = True
        changed = True
    else:
        existing = a.obs["condition"].astype(str)
        blank = existing.str.strip().eq("") & batch.isin(mapping)
        if bool(blank.any()):
            filled = existing.copy()
            filled.loc[blank] = batch.loc[blank].map(lambda value: str(mapping.get(value, "")))
            a.obs["condition"] = filled
            report["condition_filled"] = int(blank.sum())
            changed = True

    if not changed:
        report["status"] = "unchanged"
        return report

    a.obs["sample_id"] = a.obs["sample_id"].astype(str)
    a.obs["condition"] = a.obs["condition"].astype(str)
    a.uns["sc_batch_import"] = {
        "sample_key": "sample_id",
        "condition_key": "condition",
        "source": "backfill_sample_conditions",
        "mapping": dict(mapping),
    }
    if dry_run:
        report["status"] = "dry_run_would_write"
        return report

    tmp_path = path + ".tmp_backfill"
    a.write_h5ad(tmp_path)
    os.replace(tmp_path, path)
    report["status"] = "written"
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", required=True,
                        help="Project directory containing uploads/ and intermediate/")
    parser.add_argument("--mapping", required=True,
                        help='JSON object batch->condition, e.g. {"IP3": "疾病组", "NEC_251015": "对照组"}')
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would change without writing files")
    parser.add_argument("--limit", type=int, default=0,
                        help="Process at most N files (testing)")
    parser.add_argument("--audit-dir", default="",
                        help="Directory for the JSON audit (default: <project>/results)")
    args = parser.parse_args()

    mapping = json.loads(args.mapping)
    if not isinstance(mapping, dict):
        sys.exit("--mapping 必须是 JSON 对象")

    project_dir = os.path.abspath(args.project_dir)
    candidates = []
    for sub in ("uploads", "intermediate"):
        root = os.path.join(project_dir, sub)
        if not os.path.isdir(root):
            continue
        for current, _dirs, files in os.walk(root):
            for name in files:
                if name.endswith(".h5ad"):
                    candidates.append(os.path.join(current, name))
    candidates.sort()

    report = {
        "project_dir": project_dir,
        "mapping": mapping,
        "dry_run": bool(args.dry_run),
        "files": [],
        "written": 0,
        "unchanged": 0,
        "errors": [],
        "finished_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    for path in candidates[:args.limit] if args.limit else candidates:
        print(f"[backfill] {path}", flush=True)
        try:
            entry = annotate_h5ad(path, mapping, dry_run=args.dry_run)
        except Exception as exc:  # noqa: BLE001 - keep auditing other files
            report["errors"].append({"path": path, "error": str(exc)})
            print(f"  ERROR: {exc}", flush=True)
            continue
        report["files"].append(entry)
        if entry["status"] == "written":
            report["written"] += 1
        else:
            report["unchanged"] += 1
        print(f"  -> {entry['status']} "
              f"(sample_id_added={entry['sample_id_added']}, "
              f"condition_added={entry['condition_added']}, "
              f"condition_filled={entry['condition_filled']})", flush=True)

    results_dir = os.path.abspath(args.audit_dir) if args.audit_dir else os.path.join(project_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    report_path = os.path.join(
        results_dir,
        "metadata_backfill_" + time.strftime("%Y%m%d_%H%M%S") + ".json",
    )
    with open(report_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    print(f"[backfill] audit written to {report_path}")
    print(f"[backfill] done: {report['written']} written, "
          f"{report['unchanged']} unchanged, {len(report['errors'])} errors")


if __name__ == "__main__":
    main()