#!/usr/bin/env python3
"""Generate reproducible Phase 1 before/after figures from repository Bulk data.

The script deliberately uses the checked-in example project instead of random
demo data for PCA, volcano, heatmap and ORA.  No completed GSEA result is
available in that project, so the GSEA renderer is exercised with a clearly
labelled schema fixture whose pathway names and FDR values come from the real
ORA table.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figure_engine import NatureFigureDirector, export_figure


DEFAULT_PROJECT = ROOT / "data" / "projects" / "699134b0-e1c"
DEFAULT_OUTPUT = ROOT / "artifacts" / "nature_figure_phase1"


def _sample_metadata(sample_names):
    metadata = {"Cell type": [], "Genotype": [], "Treatment": []}
    for sample in sample_names:
        parts = str(sample).split("_")
        metadata["Cell type"].append(parts[0] if parts else "Unknown")
        metadata["Genotype"].append(parts[1] if len(parts) > 1 else "Unknown")
        treatment = re.sub(r"\d+$", "", parts[2]) if len(parts) > 2 else "Unknown"
        metadata["Treatment"].append(treatment)
    return metadata


def _export(director, name, spec, data, output_dir):
    figure = director.render(spec, data)
    paths, report = export_figure(
        figure,
        output_dir / name,
        spec,
        report_path=output_dir / f"{name}_nature_readiness.json",
    )
    plt.close(figure)
    return paths, report


def _build_overview(project_dir, output_dir, after_paths):
    pairs = [
        ("PCA", project_dir / "plots" / "bulk_pca.png", Path(after_paths["pca"]["png"])),
        ("Volcano", project_dir / "plots" / "bulk_deg_volcano.png", Path(after_paths["volcano"]["png"])),
        ("Heatmap", project_dir / "plots" / "bulk_heatmap.png", Path(after_paths["heatmap"]["png"])),
        ("GO ORA", project_dir / "plots" / "enrichment_ora_go_bp.png", Path(after_paths["ora"]["png"])),
    ]
    figure, axes = plt.subplots(len(pairs), 2, figsize=(11.0, 16.0), dpi=120)
    for row, (label, before_path, after_path) in enumerate(pairs):
        for column, (stage, path) in enumerate((("Before", before_path), ("After", after_path))):
            axis = axes[row, column]
            axis.imshow(plt.imread(path))
            axis.set_title(f"{stage} · {label}", loc="left", fontsize=10, fontweight="semibold")
            axis.set_axis_off()
    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02, hspace=0.16, wspace=0.05)
    overview = output_dir / "before_after_overview.png"
    figure.savefig(overview, dpi=180, facecolor="white")
    plt.close(figure)
    return overview


def validate(project_dir, output_dir):
    project_dir = Path(project_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = project_dir / "plots"
    results_dir = project_dir / "results"
    intermediate_dir = project_dir / "intermediate"

    pca_adata = anndata.read_h5ad(intermediate_dir / "bulk_pca_output.h5ad")
    heatmap_adata = anndata.read_h5ad(intermediate_dir / "bulk_heatmap_output.h5ad")
    deg = pd.read_csv(results_dir / "bulk_deg_results.csv")
    ora = pd.read_csv(results_dir / "enrichment_ora_go_bp_results.csv")
    samples = [str(value) for value in pca_adata.obs_names]
    metadata = _sample_metadata(samples)
    director = NatureFigureDirector()
    formats = ("svg", "pdf", "png")

    pca_spec = director.create_spec(
        "pca", width="single", title="Principal component analysis",
        show_legend=True, confidence_ellipse=True, formats=formats,
    )
    pca_paths, pca_report = _export(
        director, "after_pca_89mm", pca_spec,
        {
            "coordinates": np.asarray(pca_adata.obsm["X_pca"])[:, :2],
            "groups": [f"{cell} · {treatment}" for cell, treatment in zip(
                metadata["Cell type"], metadata["Treatment"]
            )],
            "batches": metadata["Genotype"],
            "samples": samples,
            "explained_variance": np.asarray(pca_adata.uns["pca"]["variance_ratio"])[:2],
        },
        output_dir,
    )

    volcano_spec = director.create_spec(
        "volcano", width="single", title="Differential expression",
        fc_threshold=1.0, fdr_threshold=0.05, label_n=8,
        label_genes=("AQP1", "COL11A1"), formats=formats,
    )
    volcano_paths, volcano_report = _export(
        director, "after_volcano_89mm", volcano_spec, deg, output_dir,
    )

    significant = deg[deg["regulation"].isin(["Up", "Down"])].copy()
    selected = pd.concat([
        significant[significant["regulation"] == "Up"].nsmallest(12, "padj"),
        significant[significant["regulation"] == "Down"].nsmallest(12, "padj"),
    ]).drop_duplicates("gene")
    gene_lookup = {str(gene): index for index, gene in enumerate(heatmap_adata.var_names)}
    selected = selected[selected["gene"].astype(str).isin(gene_lookup)].copy()
    gene_indices = [gene_lookup[str(gene)] for gene in selected["gene"]]
    heatmap_matrix = np.asarray(heatmap_adata.X)[:, gene_indices].T
    heatmap_spec = director.create_spec(
        "heatmap", width="single", title="Top differentially expressed genes",
        zscore="row", row_cluster=True, col_cluster=True,
        distance_metric="correlation", max_row_labels=24,
        max_col_labels=12, formats=formats,
    )
    heatmap_paths, heatmap_report = _export(
        director, "after_heatmap_89mm", heatmap_spec,
        {
            "matrix": heatmap_matrix,
            "gene_labels": selected["gene"].astype(str).tolist(),
            "sample_labels": samples,
            "annotations": metadata,
        },
        output_dir,
    )

    ora_spec = director.create_spec(
        "enrichment", width="single", title="GO biological process enrichment",
        top_n=min(10, len(ora)), formats=formats,
    )
    ora_paths, ora_report = _export(
        director, "after_ora_dotplot_89mm", ora_spec, ora, output_dir,
    )

    fixture_rows = ora.head(min(8, len(ora))).copy()
    gsea_fixture = pd.DataFrame({
        "Term": fixture_rows["Term"].astype(str).tolist(),
        "NES": np.array([2.35, 2.02, 1.78, 1.52, -1.45, -1.71, -1.96, -2.24])[:len(fixture_rows)],
        "FDR": pd.to_numeric(fixture_rows["Adjusted P-value"], errors="coerce").to_numpy(),
        "setSize": np.linspace(42, 126, len(fixture_rows)).round().astype(int),
    })
    gsea_spec = director.create_spec(
        "gsea", width="single", title="Gene set enrichment analysis",
        top_n=len(gsea_fixture), formats=formats,
    )
    gsea_paths, gsea_report = _export(
        director, "after_gsea_schema_fixture_89mm", gsea_spec, gsea_fixture, output_dir,
    )

    after_paths = {
        "pca": pca_paths,
        "volcano": volcano_paths,
        "heatmap": heatmap_paths,
        "ora": ora_paths,
        "gsea": gsea_paths,
    }
    overview = _build_overview(project_dir, output_dir, after_paths)
    reports = {
        "pca": pca_report,
        "volcano": volcano_report,
        "heatmap": heatmap_report,
        "ora": ora_report,
        "gsea": gsea_report,
    }
    manifest = {
        "source_project": str(project_dir.relative_to(ROOT)),
        "data_provenance": {
            "pca": "real bulk_pca_output.h5ad",
            "volcano": "real bulk_deg_results.csv",
            "heatmap": "real bulk_heatmap_output.h5ad + DEG-selected genes",
            "ora": "real enrichment_ora_go_bp_results.csv",
            "gsea": "schema fixture; real ORA terms/FDR with deterministic NES/setSize",
        },
        "outputs": {
            key: {fmt: str(Path(path).relative_to(ROOT)) for fmt, path in paths.items()}
            for key, paths in after_paths.items()
        },
        "before_after_overview": str(overview.relative_to(ROOT)),
        "readiness": {
            key: report.to_dict() for key, report in reports.items()
        },
    }
    (output_dir / "validation_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    lines = [
        "# Nature Figure Engine Phase 1 validation",
        "",
        "All templates were rendered at a physical width of 89 mm.",
        "",
        "| Template | Source | Score | Ready |",
        "|---|---|---:|:---:|",
    ]
    for key, source in manifest["data_provenance"].items():
        report = reports[key]
        lines.append(f"| {key} | {source} | {report.score}/100 | {'yes' if report.ready else 'no'} |")
    lines.extend([
        "",
        "GSEA is explicitly a schema fixture because this repository contains no completed real GSEA result.",
        "See each `*_nature_readiness.json` for issue-level QA.",
    ])
    (output_dir / "VALIDATION_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = validate(args.project_dir, args.output_dir)
    scores = {
        key: payload["nature_readiness_score"]
        for key, payload in manifest["readiness"].items()
    }
    print(json.dumps(scores, ensure_ascii=False))


if __name__ == "__main__":
    main()
