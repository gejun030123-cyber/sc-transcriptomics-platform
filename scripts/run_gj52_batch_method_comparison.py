#!/usr/bin/env python3
"""Compare batch-integration methods on the GJ52 organoid/tissue dataset.

The goal is diagnostic visualization, not choosing a final biological truth.
Each method gets its own neighbors/UMAP/Leiden result and comparable metrics.
"""

from __future__ import annotations

import argparse
import gc
import html
import importlib.util
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd


ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

PROJECT_ID = "gj52_organoid_tissue_batch"
PROJECT_DIR = ROOT_DIR / "data" / "projects" / PROJECT_ID
DEFAULT_INPUT = PROJECT_DIR / "intermediate" / "hvg_output.h5ad"
DEFAULT_FINAL = PROJECT_DIR / "intermediate" / "proportion_output.h5ad"
DEFAULT_OUTDIR = PROJECT_DIR / "results" / "batch_method_comparison_20260707"
REPORT_PATH = ROOT_DIR / "GJ52_BATCH_METHOD_COMPARISON_20260707.md"


def has_module(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return ""

    def fmt(value) -> str:
        if isinstance(value, (dict, list, tuple)):
            text = json.dumps(value, ensure_ascii=False)
        else:
            if pd.isna(value):
                return ""
            text = str(value)
        return text.replace("|", "\\|").replace("\n", " ")

    cols = list(df.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(fmt(row[col]) for col in cols) + " |")
    return "\n".join(lines)


def read_hvg(input_path: Path, n_top_genes: int | None = None):
    import scanpy as sc

    adata = sc.read_h5ad(input_path)
    if "highly_variable" in adata.var.columns:
        hv = adata.var["highly_variable"].astype(bool).values
        if n_top_genes and "highly_variable_rank" in adata.var.columns:
            ranks = pd.Series(adata.var["highly_variable_rank"].values, index=adata.var_names)
            keep_names = ranks.dropna().sort_values().head(n_top_genes).index
            hv = adata.var_names.isin(keep_names)
        adata = adata[:, hv].copy()
    if adata.n_vars < 2:
        raise ValueError("Not enough HVGs for comparison")
    adata.obs["batch"] = adata.obs["batch"].astype(str)
    return adata


def add_reference_labels(adata, final_path: Path):
    import scanpy as sc

    if not final_path.exists():
        return adata
    final = sc.read_h5ad(final_path, backed="r")
    for key in ["celltype", "leiden"]:
        if key in final.obs.columns:
            series = final.obs[key].astype(str)
            adata.obs[f"reference_{key}"] = series.reindex(adata.obs_names).fillna("NA").values
    final.file.close()
    return adata


def prepare_reference(input_path: Path, final_path: Path, n_top_genes: int, n_pcs: int):
    import scanpy as sc

    log("Loading HVG data for PCA reference")
    adata = read_hvg(input_path, n_top_genes=n_top_genes)
    adata = add_reference_labels(adata, final_path)
    log(f"Reference shape: {adata.n_obs} cells x {adata.n_vars} HVGs")
    sc.pp.scale(adata, max_value=10)
    actual_pcs = max(2, min(n_pcs, adata.n_obs - 1, adata.n_vars - 1))
    sc.tl.pca(adata, n_comps=actual_pcs, svd_solver="arpack")
    return {
        "obs": adata.obs.copy(),
        "var_names": adata.var_names.astype(str).tolist(),
        "x_pca": adata.obsm["X_pca"].copy(),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_pcs": int(actual_pcs),
    }


def skeleton(ref: dict, rep_key: str, embedding: np.ndarray):
    import anndata as ad

    x = np.zeros((ref["n_cells"], 1), dtype=np.float32)
    work = ad.AnnData(X=x, obs=ref["obs"].copy())
    work.obsm[rep_key] = np.asarray(embedding, dtype=np.float32)
    return work


def run_neighbors_umap_leiden(work, rep_key: str, resolution: float, n_neighbors: int, method: str):
    import scanpy as sc

    sc.pp.neighbors(work, n_neighbors=n_neighbors, use_rep=rep_key)
    sc.tl.umap(work, min_dist=0.4)
    sc.tl.leiden(work, resolution=resolution, key_added="leiden", flavor="igraph", n_iterations=2)
    work.uns["comparison_method"] = method
    return work


def run_baseline(ref: dict, args):
    work = skeleton(ref, "X_pca", ref["x_pca"])
    return run_neighbors_umap_leiden(work, "X_pca", args.resolution, args.n_neighbors, "baseline")


def run_harmony(ref: dict, args):
    import harmonypy as hm

    theta = args.harmony_theta
    result = hm.run_harmony(
        ref["x_pca"],
        ref["obs"],
        vars_use=[args.batch_key],
        theta=theta,
        max_iter_harmony=args.harmony_max_iter,
    )
    z = result.Z_corr
    if z.shape[0] == ref["x_pca"].shape[1]:
        z = z.T
    work = skeleton(ref, "X_pca_harmony", z)
    return run_neighbors_umap_leiden(work, "X_pca_harmony", args.resolution, args.n_neighbors, "harmony")


def run_combat(input_path: Path, final_path: Path, args):
    import scanpy as sc

    adata = read_hvg(input_path, n_top_genes=args.n_top_genes)
    adata = add_reference_labels(adata, final_path)
    sc.pp.combat(adata, key=args.batch_key)
    sc.pp.scale(adata, max_value=10)
    actual_pcs = max(2, min(args.n_pcs, adata.n_obs - 1, adata.n_vars - 1))
    sc.tl.pca(adata, n_comps=actual_pcs, svd_solver="arpack")
    adata.obsm["X_pca_combat"] = adata.obsm["X_pca"].copy()
    return run_neighbors_umap_leiden(adata, "X_pca_combat", args.resolution, args.n_neighbors, "combat")


def run_bbknn(ref: dict, args):
    import bbknn

    work = skeleton(ref, "X_pca", ref["x_pca"])
    bbknn.bbknn(
        work,
        batch_key=args.batch_key,
        neighbors_within_batch=args.bbknn_neighbors_within_batch,
        use_rep="X_pca",
    )
    import scanpy as sc

    sc.tl.umap(work, min_dist=0.4)
    sc.tl.leiden(work, resolution=args.resolution, key_added="leiden", flavor="igraph", n_iterations=2)
    return work


def run_scanorama(input_path: Path, final_path: Path, args):
    import anndata as ad
    import scanorama

    adata = read_hvg(input_path, n_top_genes=args.n_top_genes)
    adata = add_reference_labels(adata, final_path)
    batch_values = sorted(adata.obs[args.batch_key].astype(str).unique())
    batches = [adata[adata.obs[args.batch_key].astype(str) == b].copy() for b in batch_values]
    scanorama.integrate_scanpy(batches, dimred=args.n_pcs)
    embedding = np.zeros((adata.n_obs, args.n_pcs), dtype=np.float32)
    for batch_adata in batches:
        if "X_scanorama" not in batch_adata.obsm:
            raise ValueError("scanorama did not produce X_scanorama")
        corr = batch_adata.obsm["X_scanorama"]
        positions = adata.obs_names.get_indexer(batch_adata.obs_names)
        embedding[positions, : corr.shape[1]] = corr
    work = ad.AnnData(X=np.zeros((adata.n_obs, 1), dtype=np.float32), obs=adata.obs.copy())
    work.obsm["X_scanorama"] = embedding
    return run_neighbors_umap_leiden(work, "X_scanorama", args.resolution, args.n_neighbors, "scanorama")


def run_scvi(input_path: Path, final_path: Path, args):
    import scanpy as sc
    import scvi

    adata = read_hvg(input_path, n_top_genes=args.n_top_genes)
    adata = add_reference_labels(adata, final_path)
    if "counts" not in adata.layers:
        raise ValueError("scVI requires a counts layer")
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=args.batch_key)
    model = scvi.model.SCVI(
        adata,
        n_latent=args.scvi_n_latent,
        n_hidden=args.scvi_n_hidden,
        n_layers=args.scvi_n_layers,
        dropout_rate=args.scvi_dropout_rate,
    )
    train_kwargs = {"max_epochs": args.scvi_epochs}
    try:
        train_kwargs.update({"accelerator": "auto", "devices": "auto"})
        model.train(**train_kwargs)
    except TypeError:
        train_kwargs.pop("accelerator", None)
        train_kwargs.pop("devices", None)
        model.train(**train_kwargs)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    return run_neighbors_umap_leiden(adata, "X_scVI", args.resolution, args.n_neighbors, "scvi")


def run_sysvi(input_path: Path, final_path: Path, args):
    import scanpy as sc
    import scvi
    from scvi.external import SysVI

    adata = read_hvg(input_path, n_top_genes=args.n_top_genes)
    adata = add_reference_labels(adata, final_path)
    scvi.settings.seed = args.sysvi_seed

    # SysVI expects normalized/log-transformed HVG features in X for scRNA integration.
    # This differs from scVI, which models raw counts through a counts layer.
    SysVI.setup_anndata(adata=adata, batch_key=args.batch_key)
    model = SysVI(
        adata=adata,
        prior=args.sysvi_prior,
        n_prior_components=args.sysvi_n_prior_components,
        n_latent=args.scvi_n_latent,
        n_hidden=args.scvi_n_hidden,
        n_layers=args.scvi_n_layers,
        dropout_rate=args.scvi_dropout_rate,
        embed_categorical_covariates=args.sysvi_embed_categorical_covariates,
    )
    train_kwargs = {
        "max_epochs": args.sysvi_epochs,
        "check_val_every_n_epoch": 1,
        "plan_kwargs": {
            "kl_weight": args.sysvi_kl_weight,
            "z_distance_cycle_weight": args.sysvi_cycle_weight,
        },
    }
    try:
        train_kwargs.update({"accelerator": "auto", "devices": "auto"})
        model.train(**train_kwargs)
    except TypeError:
        train_kwargs.pop("accelerator", None)
        train_kwargs.pop("devices", None)
        model.train(**train_kwargs)

    adata.obsm["X_sysvi"] = model.get_latent_representation(adata=adata)
    return run_neighbors_umap_leiden(adata, "X_sysvi", args.resolution, args.n_neighbors, "sysvi")


def normalized_entropy(counts: np.ndarray) -> float:
    total = counts.sum()
    if total <= 0:
        return 0.0
    p = counts[counts > 0] / total
    if len(p) <= 1:
        return 0.0
    return float(-(p * np.log(p)).sum() / np.log(len(counts)))


def compute_metrics(work, rep_key: str, batch_key: str) -> dict:
    from sklearn.metrics import silhouette_score

    metrics = {
        "n_cells": int(work.n_obs),
        "n_clusters": int(work.obs["leiden"].nunique()) if "leiden" in work.obs else None,
    }
    labels = work.obs[batch_key].astype("category")
    sample_size = min(10000, work.n_obs)
    try:
        metrics["asw_batch"] = round(
            float(silhouette_score(work.obsm[rep_key], labels.cat.codes.values, sample_size=sample_size, random_state=7)),
            4,
        )
    except Exception as exc:
        metrics["asw_batch_error"] = str(exc)

    if "leiden" in work.obs:
        table = pd.crosstab(work.obs["leiden"].astype(str), work.obs[batch_key].astype(str))
        cluster_sizes = table.sum(axis=1)
        entropies = table.apply(lambda row: normalized_entropy(row.values.astype(float)), axis=1)
        max_fracs = table.div(cluster_sizes, axis=0).max(axis=1)
        weights = cluster_sizes / cluster_sizes.sum()
        metrics.update({
            "weighted_batch_entropy": round(float((entropies * weights).sum()), 4),
            "weighted_max_batch_fraction": round(float((max_fracs * weights).sum()), 4),
            "min_cluster_size": int(cluster_sizes.min()),
            "median_cluster_size": float(cluster_sizes.median()),
            "small_clusters_lt_50": int((cluster_sizes < 50).sum()),
        })

    if "reference_celltype" in work.obs:
        ct = work.obs["reference_celltype"].astype("category")
        if ct.nunique() > 1:
            try:
                metrics["asw_reference_celltype"] = round(
                    float(silhouette_score(work.obsm[rep_key], ct.cat.codes.values, sample_size=sample_size, random_state=7)),
                    4,
                )
            except Exception as exc:
                metrics["asw_reference_celltype_error"] = str(exc)
    return metrics


def make_scatter(work, color_key: str, title: str):
    import plotly.graph_objects as go

    coords = work.obsm["X_umap"][:, :2]
    labels = work.obs[color_key].astype(str).fillna("NA")
    fig = go.Figure()
    for value in sorted(labels.unique(), key=lambda x: (len(x), x)):
        mask = labels.values == value
        fig.add_trace(go.Scattergl(
            x=coords[mask, 0],
            y=coords[mask, 1],
            mode="markers",
            name=str(value),
            marker={"size": 3, "opacity": 0.72},
            hovertemplate=f"{html.escape(color_key)}={html.escape(str(value))}<extra></extra>",
        ))
    fig.update_layout(
        title=title,
        xaxis_title="UMAP-1",
        yaxis_title="UMAP-2",
        plot_bgcolor="white",
        width=760,
        height=560,
        legend={"itemsizing": "constant"},
    )
    return fig


def save_method_outputs(work, method: str, rep_key: str, outdir: Path, batch_key: str) -> list[dict]:
    import anndata as ad
    import plotly.io as pio

    outdir.mkdir(parents=True, exist_ok=True)
    files = []
    for color_key, suffix in [
        (batch_key, "batch"),
        ("leiden", "cluster"),
        ("reference_celltype", "reference_celltype"),
    ]:
        if color_key not in work.obs.columns:
            continue
        fig = make_scatter(work, color_key, f"{method}: UMAP by {color_key}")
        path = outdir / f"{method}_umap_{suffix}.json"
        path.write_text(pio.to_json(fig), encoding="utf-8")
        files.append({"method": method, "kind": suffix, "path": str(path)})

    h5ad_path = outdir / f"{method}_comparison.h5ad"
    keep = ad.AnnData(X=np.zeros((work.n_obs, 1), dtype=np.float32), obs=work.obs.copy())
    for key, value in work.obsm.items():
        if key in {"X_umap", rep_key}:
            keep.obsm[key] = np.asarray(value, dtype=np.float32)
    keep.uns["comparison_method"] = method
    keep.uns["comparison_embedding_key"] = rep_key
    keep.write_h5ad(h5ad_path)
    files.append({"method": method, "kind": "h5ad", "path": str(h5ad_path)})
    return files


def make_gallery(plot_files: list[dict], outdir: Path) -> Path:
    import plotly.io as pio

    sections = []
    first = True
    for item in plot_files:
        if item["kind"] == "h5ad":
            continue
        path = Path(item["path"])
        fig = pio.from_json(path.read_text(encoding="utf-8"))
        body = pio.to_html(fig, include_plotlyjs=("include" if first else False), full_html=False)
        first = False
        sections.append(f"<section><h2>{html.escape(item['method'])} - {html.escape(item['kind'])}</h2>{body}</section>")
    doc = """<!doctype html>
<html><head><meta charset="utf-8"><title>GJ52 batch method comparison</title>
<style>
body{font-family:Arial,sans-serif;margin:24px;background:#f7f8fb;color:#1f2933}
section{background:#fff;border:1px solid #dde3ea;border-radius:8px;margin:0 0 20px;padding:16px}
h1{margin:0 0 8px} h2{font-size:16px;margin:0 0 12px}
.meta{color:#667085;margin-bottom:20px}
</style></head><body>
<h1>GJ52 Batch Method Comparison</h1>
<div class="meta">Baseline, ComBat, Harmony, BBKNN, Scanorama, and scVI where dependencies are available.</div>
""" + "\n".join(sections or ["<p>No plots were produced.</p>"]) + "\n</body></html>\n"
    gallery = outdir / "gj52_batch_method_comparison.html"
    gallery.write_text(doc, encoding="utf-8")
    return gallery


def write_report(summary: list[dict], skipped: list[dict], files: list[dict], gallery: Path, args) -> None:
    df = pd.DataFrame(summary)
    csv_path = args.outdir / "batch_method_metrics.csv"
    df.to_csv(csv_path, index=False)

    lines = [
        "# GJ52 batch method comparison",
        "",
        f"- Input: `{args.input}`",
        f"- Batch key: `{args.batch_key}`",
        f"- Resolution: `{args.resolution}`",
        f"- Gallery: `{gallery}`",
        f"- Metrics CSV: `{csv_path}`",
        "",
        "## Completed methods",
        "",
    ]
    if not df.empty:
        lines.append(dataframe_to_markdown(df))
    else:
        lines.append("No methods completed.")

    lines.extend(["", "## Skipped or failed methods", ""])
    if skipped:
        lines.append(dataframe_to_markdown(pd.DataFrame(skipped)))
    else:
        lines.append("None.")

    lines.extend([
        "",
        "## Interpretation notes",
        "",
        "- Lower absolute batch ASW suggests stronger batch mixing, but for organoid vs tissue this can also mean biological differences were compressed.",
        "- Higher weighted batch entropy means clusters contain both organoid and tissue more evenly.",
        "- Lower weighted max batch fraction means fewer clusters are dominated by a single source.",
        "- Reference celltype labels come from the previous pipeline annotation and are used only as a visual guide.",
        "",
        "## Output files",
        "",
    ])
    for item in files:
        lines.append(f"- `{item['method']}` `{item['kind']}`: `{item['path']}`")
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--final", type=Path, default=DEFAULT_FINAL)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--methods", default="baseline,combat,harmony,bbknn,scanorama,scvi")
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--resolution", type=float, default=0.8)
    parser.add_argument("--harmony-theta", type=float, default=1.0)
    parser.add_argument("--harmony-max-iter", type=int, default=20)
    parser.add_argument("--bbknn-neighbors-within-batch", type=int, default=5)
    parser.add_argument("--scvi-epochs", type=int, default=80)
    parser.add_argument("--scvi-n-latent", type=int, default=30)
    parser.add_argument("--scvi-n-hidden", type=int, default=128)
    parser.add_argument("--scvi-n-layers", type=int, default=1)
    parser.add_argument("--scvi-dropout-rate", type=float, default=0.1)
    parser.add_argument("--sysvi-epochs", type=int, default=60)
    parser.add_argument("--sysvi-cycle-weight", type=float, default=5.0)
    parser.add_argument("--sysvi-kl-weight", type=float, default=1.0)
    parser.add_argument("--sysvi-prior", default="vamp", choices=["vamp", "standard_normal"])
    parser.add_argument("--sysvi-n-prior-components", type=int, default=5)
    parser.add_argument("--sysvi-seed", type=int, default=0)
    parser.add_argument("--sysvi-embed-categorical-covariates", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    methods = [m.strip().lower() for m in args.methods.split(",") if m.strip()]
    log(f"Requested methods: {methods}")
    log(f"Output dir: {args.outdir}")

    ref = prepare_reference(args.input, args.final, args.n_top_genes, args.n_pcs)
    summary = []
    skipped = []
    plot_files = []

    runners = {
        "baseline": (None, lambda: ("X_pca", run_baseline(ref, args))),
        "harmony": ("harmonypy", lambda: ("X_pca_harmony", run_harmony(ref, args))),
        "combat": (None, lambda: ("X_pca_combat", run_combat(args.input, args.final, args))),
        "bbknn": ("bbknn", lambda: ("X_pca", run_bbknn(ref, args))),
        "scanorama": ("scanorama", lambda: ("X_scanorama", run_scanorama(args.input, args.final, args))),
        "scvi": ("scvi", lambda: ("X_scVI", run_scvi(args.input, args.final, args))),
        "sysvi": ("scvi", lambda: ("X_sysvi", run_sysvi(args.input, args.final, args))),
    }

    for method in methods:
        if method not in runners:
            skipped.append({"method": method, "reason": "unknown method"})
            continue
        module_name, runner = runners[method]
        if module_name and not has_module(module_name):
            skipped.append({"method": method, "reason": f"missing Python package: {module_name}"})
            continue
        log(f"Running {method}")
        start = time.time()
        try:
            rep_key, work = runner()
            metrics = compute_metrics(work, rep_key, args.batch_key)
            metrics.update({
                "method": method,
                "embedding_key": rep_key,
                "runtime_sec": round(time.time() - start, 1),
            })
            summary.append(metrics)
            plot_files.extend(save_method_outputs(work, method, rep_key, args.outdir, args.batch_key))
            del work
            gc.collect()
        except Exception as exc:
            skipped.append({
                "method": method,
                "reason": str(exc),
                "traceback": traceback.format_exc(limit=4),
            })
            log(f"{method} failed: {exc}")

    gallery = make_gallery(plot_files, args.outdir)
    write_report(summary, skipped, plot_files, gallery, args)
    print(json.dumps({
        "report": str(REPORT_PATH),
        "gallery": str(gallery),
        "metrics_csv": str(args.outdir / "batch_method_metrics.csv"),
        "completed": [row["method"] for row in summary],
        "skipped": skipped,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
