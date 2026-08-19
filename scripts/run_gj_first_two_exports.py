#!/usr/bin/env python3
"""Run four independent GJ scRNA projects and export only package folders 01/02.

The runner intentionally treats each top-level source directory as a separate
project.  It never infers biological conditions from similar-looking sample
names: until a reviewed design is supplied, ``condition`` equals ``sample_id``.
Intermediate h5ad files are retained under ``90_analysis_work`` so an
interrupted run can resume without repeating completed stages.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import json
import os
import shutil
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

DEFAULT_SOURCE = Path("/home/oelab/AnaData/GJ/00_source_matrices")
DEFAULT_OUTPUT = Path("/home/oelab/AnaData/GJ/01_analysis_results")
PROJECT_ORDER = ("CXY_megtacolon", "WYQ_VLO", "LHS", "WYQ_rett")

PIPELINE = (
    ("qc", {
        "mito_perc": 0.20,
        "nUMIs": 500,
        "detected_genes": 250,
        "max_detected_genes": 0,
        "ribo_perc": 0,
        "hb_perc": 0,
        "batch_key": "sample_id",
        "batch_adaptive_qc": False,
        "save_counts_layer": True,
        "score_cell_cycle": False,
        "show_qc_filter_summary": True,
        "show_doublet_histogram": True,
    }),
    ("normalize", {
        "method": "log1p",
        "target_sum": 10000,
        "show_expression_distribution": True,
    }),
    ("hvg", {
        "n_top_genes": 2000,
        "hvg_flavor": "seurat_v3",
        "batch_key": "",
        "exclude_mt_genes": True,
        "exclude_cc_genes": False,
        "show_hvg_rank_plot": True,
    }),
    ("dimred", {
        "n_comps": 50,
        "pca_hvg_only": True,
        "auto_n_comps": "none",
        "umap_n_neighbors": 20,
        "umap_min_dist": 0.30,
        "umap_metric": "euclidean",
        "umap_spread": 1.0,
        "enable_tsne": False,
        "show_pca_scatter": True,
    }),
    ("clustering", {
        "resolutions": "0.8",
        "primary_resolution": "0.8",
        "n_neighbors": 20,
        "clustering_method": "leiden",
        "n_iterations": 2,
        "distance_metric": "euclidean",
        "use_corrected": False,
        "batch_key": "sample_id",
        "compute_marker_preview": False,
        "show_cluster_size_bar": True,
        "show_cluster_batch_composition": True,
        "show_labeled_umap": True,
        "show_resolution_sankey": False,
    }),
    ("sc_csv_export", {
        "export_scope": "proportions_expression",
        "sample_key": "sample_id",
        "condition_key": "condition",
        "cluster_key": "leiden",
        "export_folder": "sc_batch_results",
        "export_prefix": "sc_batch",
    }),
)


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def json_default(value):
    if hasattr(value, "item"):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=json_default),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def load_state(path: Path) -> dict:
    if not path.is_file():
        return {"created_at": now(), "stages": {}}
    return json.loads(path.read_text(encoding="utf-8"))


def barcode_count(matrix_dir: Path) -> int:
    with gzip.open(matrix_dir / "barcodes.tsv.gz", "rt", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def discover_project(source_dir: Path, project: str) -> pd.DataFrame:
    from modules.sc_batch import find_10x_matrix_dirs

    project_source = source_dir / project
    matrix_dirs = [Path(item) for item in find_10x_matrix_dirs(str(project_source))]
    if not matrix_dirs:
        raise FileNotFoundError(f"No complete 10x matrix directories found: {project_source}")

    rows = []
    for matrix_dir in matrix_dirs:
        sample_id = matrix_dir.name
        rows.append({
            "sample_id": sample_id,
            "matrix_dir": str(matrix_dir.relative_to(project_source)),
            "condition": sample_id,
            "replicate": "1",
            "batch": sample_id,
            "experiment_id": project,
            "raw_filtered_cells": barcode_count(matrix_dir),
            "notes": "Condition is intentionally sample_id; no biological grouping inferred.",
        })
    frame = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    if frame["sample_id"].duplicated().any():
        raise ValueError(f"Duplicate sample_id in {project}")
    return frame


def write_manifest(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "sample_id", "matrix_dir", "condition", "replicate", "batch",
        "experiment_id", "notes",
    ]
    frame[columns].to_csv(path, index=False)


def completed_stage(state: dict, stage: str) -> str | None:
    record = state.get("stages", {}).get(stage, {})
    output = str(record.get("output_adata") or "")
    if record.get("status") == "completed" and output and Path(output).is_file():
        return output
    return None


def run_module(project: str, work_dir: Path, stage: str, params: dict,
               input_path: str | None, state: dict, state_path: Path) -> tuple[str, dict]:
    from modules import MODULE_REGISTRY

    previous = completed_stage(state, stage)
    if previous:
        log(f"{project}: resume — skip completed stage {stage}")
        return previous, state["stages"][stage].get("summary", {})

    log(f"{project}: start {stage}")

    def progress(percent: int, message: str) -> None:
        label = "warning" if percent < 0 else f"{percent}%"
        log(f"{project}/{stage} [{label}] {message}")

    module = MODULE_REGISTRY[stage](
        project_dir=str(work_dir), params=dict(params), progress_callback=progress,
    )
    started = now()
    try:
        result = module.run(input_path)
        output_adata = str(result.get("output_adata") or "")
        if not output_adata or not Path(output_adata).is_file():
            raise RuntimeError(f"{stage} did not produce a readable output_adata")
        record = {
            "status": "completed",
            "started_at": started,
            "finished_at": now(),
            "output_adata": output_adata,
            "summary": result.get("summary", {}),
            "result_files": result.get("result_files", []),
        }
        state.setdefault("stages", {})[stage] = record
        state["updated_at"] = now()
        atomic_json(state_path, state)
        log(f"{project}: completed {stage}")
        del module
        gc.collect()
        return output_adata, record["summary"]
    except Exception as exc:
        state.setdefault("stages", {})[stage] = {
            "status": "failed",
            "started_at": started,
            "finished_at": now(),
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        state["updated_at"] = now()
        atomic_json(state_path, state)
        raise


def run_import(project: str, source_dir: Path, work_dir: Path, manifest_path: Path,
               state: dict, state_path: Path) -> tuple[str, dict]:
    from modules.sc_batch import SCBatchImport

    previous = completed_stage(state, "sc_batch_import")
    if previous:
        log(f"{project}: resume — skip completed stage sc_batch_import")
        return previous, state["stages"]["sc_batch_import"].get("summary", {})

    log(f"{project}: start sc_batch_import")

    def progress(percent: int, message: str) -> None:
        label = "warning" if percent < 0 else f"{percent}%"
        log(f"{project}/sc_batch_import [{label}] {message}")

    module = SCBatchImport(
        project_dir=str(work_dir),
        params={
            "manifest_path": str(manifest_path),
            "source_root": str(source_dir / project),
            "dataset_name": project,
            "species": "Human",
            "genome": "GRCh38",
        },
        progress_callback=progress,
    )
    started = now()
    try:
        result = module.run(None)
        output_adata = str(result.get("output_adata") or "")
        if not output_adata or not Path(output_adata).is_file():
            raise RuntimeError("sc_batch_import did not produce a readable h5ad")
        state.setdefault("stages", {})["sc_batch_import"] = {
            "status": "completed",
            "started_at": started,
            "finished_at": now(),
            "output_adata": output_adata,
            "summary": result.get("summary", {}),
            "result_files": result.get("result_files", []),
        }
        state["updated_at"] = now()
        atomic_json(state_path, state)
        log(f"{project}: completed sc_batch_import")
        del module
        gc.collect()
        return output_adata, result.get("summary", {})
    except Exception as exc:
        state.setdefault("stages", {})["sc_batch_import"] = {
            "status": "failed",
            "started_at": started,
            "finished_at": now(),
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        state["updated_at"] = now()
        atomic_json(state_path, state)
        raise


def publish_outputs(project_root: Path, work_dir: Path, frame: pd.DataFrame,
                    state: dict) -> None:
    package = work_dir / "results" / "sc_batch_results"
    for folder in ("01_group_proportions", "02_gene_expression"):
        source = package / folder
        if not source.is_dir():
            raise FileNotFoundError(f"Missing expected result folder: {source}")
        shutil.copytree(source, project_root / folder, dirs_exist_ok=True)

    frame.to_csv(project_root / "sample_inventory.csv", index=False)
    export_summary = state.get("stages", {}).get("sc_csv_export", {}).get("summary", {})
    readme = f"""# {project_root.name} single-cell batch results

- Generated: {now()}
- Source matrices: `{DEFAULT_SOURCE / project_root.name}`
- Raw filtered cells: {int(frame['raw_filtered_cells'].sum())}
- Samples: {len(frame)}
- Cells after QC: {export_summary.get('n_cells', 'see run_state.json')}
- Export scope: `proportions_expression`
- Biological grouping: not inferred; `condition` currently equals `sample_id`.

Only `01_group_proportions` and `02_gene_expression` are published.  Reusable
h5ad intermediates, plots, module summaries, and the reviewed manifest are in
`90_analysis_work`.
"""
    (project_root / "README.md").write_text(readme, encoding="utf-8")


def run_project(source_dir: Path, output_dir: Path, project: str) -> dict:
    project_root = output_dir / project
    work_dir = project_root / "90_analysis_work"
    for child in ("uploads", "intermediate", "results", "plots"):
        (work_dir / child).mkdir(parents=True, exist_ok=True)

    state_path = work_dir / "run_state.json"
    state = load_state(state_path)
    state.update({
        "project": project,
        "source_dir": str(source_dir / project),
        "project_root": str(project_root),
        "pipeline": [name for name, _ in PIPELINE],
    })

    frame = discover_project(source_dir, project)
    manifest_path = work_dir / "uploads" / f"{project}_manifest.csv"
    write_manifest(frame, manifest_path)
    state["n_samples"] = int(len(frame))
    state["raw_filtered_cells"] = int(frame["raw_filtered_cells"].sum())
    atomic_json(state_path, state)

    current_input, import_summary = run_import(
        project, source_dir, work_dir, manifest_path, state, state_path,
    )
    final_summary = import_summary
    for stage, params in PIPELINE:
        current_input, final_summary = run_module(
            project, work_dir, stage, params, current_input, state, state_path,
        )

    publish_outputs(project_root, work_dir, frame, state)
    return {
        "project": project,
        "status": "completed",
        "n_samples": int(len(frame)),
        "raw_filtered_cells": int(frame["raw_filtered_cells"].sum()),
        "cells_after_qc": int(final_summary.get("n_cells", 0)),
        "result_root": str(project_root),
    }


def write_quality_summaries(output_dir: Path, projects: tuple[str, ...]) -> None:
    """Create compact project/sample QC indexes from the published CSV files."""
    project_rows = []
    sample_frames = []
    for project in projects:
        project_root = output_dir / project
        inventory = pd.read_csv(project_root / "sample_inventory.csv")
        proportions = pd.read_csv(
            project_root / "01_group_proportions"
            / "sc_batch_sample_cluster_proportions.csv"
        )
        retained = (
            proportions.groupby("sample_id", observed=True)["total_cells"]
            .first().astype(int)
        )
        sample_qc = inventory[["sample_id", "raw_filtered_cells"]].copy()
        sample_qc.insert(0, "project", project)
        sample_qc["cells_after_qc"] = (
            sample_qc["sample_id"].map(retained).fillna(0).astype(int)
        )
        sample_qc["cells_removed"] = (
            sample_qc["raw_filtered_cells"] - sample_qc["cells_after_qc"]
        )
        sample_qc["retention_rate"] = (
            sample_qc["cells_after_qc"]
            / sample_qc["raw_filtered_cells"].where(
                sample_qc["raw_filtered_cells"] > 0, pd.NA
            )
        )
        sample_qc["review_flag"] = ""
        sample_qc.loc[sample_qc["cells_after_qc"] < 500, "review_flag"] = (
            "low_post_qc_cells"
        )
        low_retention = sample_qc["retention_rate"] < 0.85
        sample_qc.loc[low_retention, "review_flag"] = sample_qc.loc[
            low_retention, "review_flag"
        ].map(lambda value: f"{value};low_retention".strip(";"))
        sample_frames.append(sample_qc)

        counts_path = (
            project_root / "02_gene_expression"
            / "sc_batch_pseudobulk_counts.csv"
        )
        with counts_path.open("r", encoding="utf-8") as handle:
            n_genes = max(sum(1 for _ in handle) - 1, 0)
        project_rows.append({
            "project": project,
            "status": "completed",
            "n_samples": int(len(sample_qc)),
            "raw_filtered_cells": int(sample_qc["raw_filtered_cells"].sum()),
            "cells_after_qc": int(sample_qc["cells_after_qc"].sum()),
            "cells_removed": int(sample_qc["cells_removed"].sum()),
            "retention_rate": float(
                sample_qc["cells_after_qc"].sum()
                / sample_qc["raw_filtered_cells"].sum()
            ),
            "n_clusters": int(proportions["cluster"].astype(str).nunique()),
            "n_genes_exported": n_genes,
            "samples_flagged_for_review": int(
                sample_qc["review_flag"].ne("").sum()
            ),
            "result_root": str(project_root),
        })

    pd.DataFrame(project_rows).to_csv(
        output_dir / "batch_quality_summary.csv", index=False,
    )
    pd.concat(sample_frames, ignore_index=True).to_csv(
        output_dir / "sample_quality_summary.csv", index=False,
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--projects", default=",".join(PROJECT_ORDER),
        help="Comma-separated top-level project directories.",
    )
    parser.add_argument("--list-only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    projects = tuple(item.strip() for item in args.projects.split(",") if item.strip())
    unknown = sorted(set(projects) - set(PROJECT_ORDER))
    if unknown:
        raise ValueError(f"Unknown projects: {', '.join(unknown)}")
    if not args.source.is_dir():
        raise FileNotFoundError(args.source)

    if args.list_only:
        for project in projects:
            frame = discover_project(args.source, project)
            log(
                f"{project}: {len(frame)} samples, "
                f"{int(frame['raw_filtered_cells'].sum())} raw filtered cells"
            )
        return 0

    args.output.mkdir(parents=True, exist_ok=True)
    os.environ["SC_BATCH_SOURCE_ROOTS"] = str(args.source)
    runtime_tmp = args.output / "00_runtime_tmp"
    os.environ["RUNTIME_TMP_DIR"] = str(runtime_tmp)
    os.environ["NUMBA_CACHE_DIR"] = str(runtime_tmp / "numba_cache")
    os.environ["MPLCONFIGDIR"] = str(runtime_tmp / "mplconfig")

    summaries = []
    for project in projects:
        log(f"===== {project} =====")
        try:
            summaries.append(run_project(args.source, args.output, project))
        except Exception as exc:
            summaries.append({
                "project": project,
                "status": "failed",
                "error": str(exc),
                "result_root": str(args.output / project),
            })
            pd.DataFrame(summaries).to_csv(args.output / "batch_run_summary.csv", index=False)
            log(traceback.format_exc())
            return 1
        pd.DataFrame(summaries).to_csv(args.output / "batch_run_summary.csv", index=False)

    write_quality_summaries(args.output, projects)
    log("All requested projects completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
