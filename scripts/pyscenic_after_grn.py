#!/usr/bin/env python3
"""Continue one locally managed pySCENIC run after its GRN stage exits.

This helper is deliberately parameterised only with trusted server paths.  It
records durable status and logs, runs the pinned Docker image with argument
arrays, and never transmits single-cell data off the server.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np

# The locally managed pySCENIC environment currently supplies LoomPy 3.0.x
# alongside NumPy 2.  LoomPy still refers to these removed dtype aliases when
# opening the AUCell Loom, so restore their previous meanings before import.
# This is a compatibility shim only; it does not transform any matrix values.
if not hasattr(np, "string_"):
    np.string_ = np.bytes_
if not hasattr(np, "object_"):
    np.object_ = np.dtype("O").type

import loompy
import pandas as pd

# The controller can be started with a minimal PYTHONPATH for its isolated
# Loom writer dependency.  Add the repository root explicitly so it can reuse
# the platform's audited SCENIC display module as well.
REPOSITORY_ROOT = str(Path(__file__).resolve().parents[1])
if REPOSITORY_ROOT not in sys.path:
    sys.path.insert(0, REPOSITORY_ROOT)

from modules.scenic import ScenicAnalysis


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--container-name", required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--resource-dir", type=Path, required=True)
    parser.add_argument("--source-h5ad", type=Path, required=True)
    parser.add_argument("--poll-seconds", type=int, default=60)
    return parser.parse_args()


def read_state(container_name: str) -> dict:
    command = ["docker", "inspect", container_name, "--format", "{{json .State}}"]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Docker container unavailable: {completed.stderr.strip()}")
    return json.loads(completed.stdout)


def snapshot_docker_logs(container_name: str, destination: Path) -> None:
    with destination.open("w", encoding="utf-8") as handle:
        subprocess.run(
            ["docker", "logs", "--timestamps", container_name],
            stdout=handle, stderr=subprocess.STDOUT, text=True, check=False,
        )


def run_docker(command: list[str], log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Docker stage failed ({completed.returncode}); see {log_path.name}.")


def extract_targets(regulons_csv: Path) -> dict[str, list[str]]:
    frame = pd.read_csv(regulons_csv)
    lower_columns = {str(column).lower(): str(column) for column in frame.columns}
    regulon_column = next((lower_columns[key] for key in ("regulon", "tf") if key in lower_columns), None)
    targets_column = next((lower_columns[key] for key in ("targetgenes", "target_genes", "targets") if key in lower_columns), None)
    if not regulon_column or not targets_column:
        return {}
    targets: dict[str, list[str]] = {}
    for row in frame[[regulon_column, targets_column]].dropna().itertuples(index=False):
        regulon = str(row[0]).strip()
        genes = [gene.strip() for gene in re.split(r"[;,]", str(row[1])) if gene.strip()]
        if regulon and genes:
            targets[regulon] = list(dict.fromkeys(genes))
    return targets


def auc_dataframe(aucell_loom: Path) -> pd.DataFrame:
    with loompy.connect(str(aucell_loom), mode="r", validate=False) as dataset:
        if "RegulonsAUC" not in dataset.ca:
            raise ValueError("AUCell Loom 缺少 RegulonsAUC 列属性。")
        raw = np.asarray(dataset.ca["RegulonsAUC"])
        cell_ids = np.asarray(dataset.ca["CellID"]).astype(str)
    if raw.dtype.names:
        frame = pd.DataFrame({name: raw[name] for name in raw.dtype.names}, index=cell_ids)
    elif raw.ndim == 2:
        raise ValueError("AUCell 矩阵没有 regulon 列名，不能安全导入 AnnData。")
    else:
        raise ValueError("无法识别 AUCell Loom 中的 RegulonsAUC 格式。")
    if frame.index.has_duplicates or frame.columns.has_duplicates or not np.isfinite(frame.to_numpy(dtype=float)).all():
        raise ValueError("AUCell 输出包含重复 ID/调控子或非有限值。")
    return frame


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    resource_dir = args.resource_dir.resolve()
    source_h5ad = args.source_h5ad.resolve()
    status_path = run_dir / "run_status.json"

    def update_status(**updates) -> None:
        current = json.loads(status_path.read_text(encoding="utf-8")) if status_path.exists() else {}
        current.update(updates)
        current["updated_utc"] = utc_now()
        status_path.write_text(json.dumps(current, ensure_ascii=False, indent=2), encoding="utf-8")

    try:
        while True:
            state = read_state(args.container_name)
            status = str(state.get("Status", "unknown"))
            update_status(state="running" if status == "running" else status, stage="grnboost2_network_inference",
                          container_pid=state.get("Pid"), container_status=status)
            if status == "running":
                time.sleep(max(10, args.poll_seconds))
                continue
            if status != "exited" or int(state.get("ExitCode", 1)) != 0:
                snapshot_docker_logs(args.container_name, run_dir / "grn_docker.log")
                update_status(state="failed", stage="grnboost2_network_inference",
                              exit_code=state.get("ExitCode"), error=state.get("Error") or "GRN container did not exit successfully")
                return 1
            break

        snapshot_docker_logs(args.container_name, run_dir / "grn_docker.log")
        adjacency = run_dir / "adjacencies.tsv"
        if not adjacency.exists() or adjacency.stat().st_size == 0:
            raise ValueError("GRN 阶段未生成非空 adjacencies.tsv。")

        image = "aertslab/pyscenic:0.12.1"
        mounts = [
            "-v", f"{run_dir}:/work",
            "-v", f"{resource_dir}:/resources:ro",
        ]
        update_status(state="running", stage="cistarget_motif_pruning")
        run_docker([
            "docker", "run", "--rm", "--user", "1000:1000", "--cpus=16", "--memory=96g", *mounts,
            "--entrypoint", "pyscenic", image, "ctx", "/work/adjacencies.tsv",
            "/resources/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
            "/resources/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
            "--annotations_fname", "/resources/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
            "--expression_mtx_fname", "/work/input_counts_hvg10k_plus_tfs.loom",
            "--cell_id_attribute", "CellID", "--gene_attribute", "Gene",
            "--num_workers", "16", "--output", "/work/regulons.csv",
        ], run_dir / "ctx.log")
        regulons = run_dir / "regulons.csv"
        if not regulons.exists() or regulons.stat().st_size == 0:
            raise ValueError("cisTarget 阶段未生成非空 regulons.csv。")

        update_status(state="running", stage="aucell_scoring")
        run_docker([
            "docker", "run", "--rm", "--user", "1000:1000", "--cpus=16", "--memory=96g", *mounts,
            "--entrypoint", "pyscenic", image, "aucell", "/work/input_counts_hvg10k_plus_tfs.loom",
            "/work/regulons.csv", "--cell_id_attribute", "CellID", "--gene_attribute", "Gene",
            "--num_workers", "16", "--seed", "1337", "--output", "/work/aucell.loom",
        ], run_dir / "aucell.log")

        update_status(state="running", stage="import_aucell_and_generate_scenic_panels")
        auc = auc_dataframe(run_dir / "aucell.loom")
        adata = ad.read_h5ad(source_h5ad)
        missing_cells = adata.obs_names.astype(str).difference(auc.index)
        if len(missing_cells):
            raise ValueError(f"AUCell 输出缺少 {len(missing_cells)} 个源细胞。")
        adata.obsm["X_aucell"] = auc.loc[adata.obs_names.astype(str)].copy()
        targets = extract_targets(regulons)
        adata.uns["scenic_regulon_targets"] = targets
        adata.uns["scenic"] = {
            "method": "pySCENIC 0.12.1 / GRNBoost2 + cisTarget v10 + AUCell",
            "input_layer": "counts",
            "run_manifest": "run_manifest.json",
            "n_regulons": int(auc.shape[1]),
            "n_regulon_target_sets": int(len(targets)),
        }
        augmented_h5ad = run_dir / "pyscenic_aucell_annotated.h5ad"
        adata.write_h5ad(augmented_h5ad)

        def analysis_progress(percent: int, message: str) -> None:
            update_status(state="running", stage="scenic_aucell_rss_visualization", display_progress=int(percent), display_message=str(message))

        analysis = ScenicAnalysis(str(run_dir), {
            "_analysis_id": "full_pyscenic",
            "groupby": "final_annotation",
            "auc_obsm_key": "X_aucell",
            "regulon_targets_uns_key": "scenic_regulon_targets",
            "top_n_regulons": 30,
            "max_heatmap_cells": 1000,
            "top_n_rss": 15,
            "show_network": True,
            "network_max_targets_per_regulon": 12,
            "show_correlation": True,
            "top_n_correlation_regulons": 30,
        }, analysis_progress)
        display = analysis.run(str(augmented_h5ad))
        update_status(state="completed", stage="completed", exit_code=0,
                      n_regulons=int(auc.shape[1]), n_target_sets=int(len(targets)),
                      augmented_h5ad=augmented_h5ad.name, display_summary=display.get("summary", {}))
        return 0
    except Exception as exc:  # status file is the durable handoff for a long run.
        update_status(state="failed", error=f"{type(exc).__name__}: {exc}")
        raise


if __name__ == "__main__":
    sys.exit(main())
