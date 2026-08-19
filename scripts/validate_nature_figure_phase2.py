#!/usr/bin/env python3
"""Render and audit Phase 2 Nature templates from repository Bulk RNA data."""

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

from figure_engine import (
    FigurePanel, FigureSpec, NatureFigureComposer, NatureFigureDirector,
    export_figure, figure_signature, simulate_image,
)


DEFAULT_PROJECT = ROOT / 'data' / 'projects' / '699134b0-e1c'
DEFAULT_OUTPUT = ROOT / 'artifacts' / 'nature_figure_phase2'


def _dense(matrix):
    return matrix.toarray() if hasattr(matrix, 'toarray') else np.asarray(matrix)


def _metadata(sample_names):
    cell_type, treatment, genotype = [], [], []
    for name in sample_names:
        parts = str(name).split('_')
        cell_type.append(parts[0] if parts else 'Unknown')
        genotype.append(parts[1] if len(parts) > 1 else 'Unknown')
        treatment.append(re.sub(r'\d+$', '', parts[2]) if len(parts) > 2 else 'Unknown')
    return {'Cell type': cell_type, 'Treatment': treatment, 'Genotype': genotype}


def _export(director, name, spec, data, output_dir):
    figure = director.render(spec, data)
    paths, report = export_figure(
        figure, output_dir / name, spec,
        report_path=output_dir / f'{name}_nature_readiness.json',
    )
    signature = figure_signature(figure)
    plt.close(figure)
    return paths, report, signature


def _pathway_scores(expression):
    variances = np.nanvar(expression, axis=0)
    top = np.argsort(variances)[-64:]
    return np.vstack([
        np.nanmean(expression[:, indices], axis=1)
        for indices in np.array_split(top, 8)
    ])


def _module_trait(expression, metadata):
    variances = np.nanvar(expression, axis=0)
    top = np.argsort(variances)[-72:]
    eigengenes = np.vstack([
        np.nanmean(expression[:, indices], axis=1)
        for indices in np.array_split(top, 6)
    ])
    trait_names, trait_values = [], []
    for name in ('Cell type', 'Treatment'):
        values = metadata[name]
        for level in list(dict.fromkeys(values))[:3]:
            trait_names.append(f'{name}: {level}')
            trait_values.append(np.asarray([value == level for value in values], dtype=float))
    traits = np.vstack(trait_values)
    correlations = np.corrcoef(np.vstack([eigengenes, traits]))[:len(eigengenes), len(eigengenes):]
    return correlations, trait_names


def _upset_sets(results_dir):
    paths = sorted(results_dir.glob('bulk_deg_results*.csv'))[:4]
    sets = {}
    for path in paths:
        frame = pd.read_csv(path)
        if {'gene', 'padj', 'log2FC'}.issubset(frame.columns):
            selected = frame[(frame['padj'] < 0.05) & (frame['log2FC'].abs() >= 1.0)]
            sets[path.stem.replace('bulk_deg_results', 'Contrast') or 'Contrast'] = set(
                selected['gene'].dropna().astype(str)
            )
    return sets


def _before_after(project_dir, output_dir, after):
    pairs = [
        ('MA', project_dir / 'plots' / 'bulk_deg_ma.png', Path(after['ma']['png'])),
        ('Correlation', project_dir / 'plots' / 'bulk_corr_heatmap.png',
         Path(after['correlation']['png'])),
    ]
    available = [pair for pair in pairs if pair[1].is_file() and pair[2].is_file()]
    if not available:
        return None
    figure, axes = plt.subplots(len(available), 2, figsize=(10, 4.6 * len(available)), dpi=120)
    axes = np.asarray(axes).reshape(len(available), 2)
    for row, (label, before_path, after_path) in enumerate(available):
        for column, (stage, path) in enumerate((('Before', before_path), ('After', after_path))):
            axes[row, column].imshow(plt.imread(path))
            axes[row, column].set_title(f'{stage} · {label}', loc='left', fontsize=10,
                                        fontweight='semibold')
            axes[row, column].set_axis_off()
    figure.patch.set_facecolor('white')
    figure.tight_layout()
    path = output_dir / 'phase2_before_after.png'
    figure.savefig(path, dpi=180, facecolor='white')
    plt.close(figure)
    return path


def validate(project_dir, output_dir):
    project_dir = Path(project_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    intermediate = project_dir / 'intermediate'
    results = project_dir / 'results'
    pca_adata = anndata.read_h5ad(intermediate / 'bulk_pca_output.h5ad')
    heatmap_adata = anndata.read_h5ad(intermediate / 'bulk_heatmap_output.h5ad')
    deg = pd.read_csv(results / 'bulk_deg_results.csv')
    ora = pd.read_csv(results / 'enrichment_ora_go_bp_results.csv')
    expression = _dense(heatmap_adata.X).astype(float)
    samples = [str(value) for value in heatmap_adata.obs_names]
    metadata = _metadata(samples)
    director = NatureFigureDirector()
    formats = ('svg', 'pdf', 'png', 'tiff')

    deg['mean_expression'] = (
        pd.to_numeric(deg['mean_group1'], errors='coerce')
        + pd.to_numeric(deg['mean_group2'], errors='coerce')
    ) / 2.0
    ma_spec = director.create_spec(
        'ma', title='Mean expression and differential signal', label_n=6,
        fc_threshold=1.0, fdr_threshold=0.05, formats=formats,
    )
    corr_spec = director.create_spec(
        'correlation', width='double', title='Sample correlation',
        max_col_labels=18, formats=formats,
    )
    corr = np.corrcoef(expression)
    scores = _pathway_scores(expression)
    score_payload = {
        'scores': scores,
        'pathway_labels': [f'Expression-derived pathway block {index + 1}'
                           for index in range(scores.shape[0])],
        'sample_labels': samples,
        'annotations': metadata,
    }
    gsva_spec = director.create_spec(
        'gsva', width='double', title='GSVA score matrix', formats=formats,
        max_row_labels=12, max_col_labels=18,
    )
    ssgsea_spec = director.create_spec(
        'ssgsea', width='double', title='ssGSEA score matrix', formats=formats,
        max_row_labels=12, max_col_labels=18,
    )
    upset_payload = _upset_sets(results)
    upset_spec = director.create_spec(
        'upset', width='double', title='DEG intersections', formats=formats,
        top_intersections=16,
    )
    module_trait, trait_names = _module_trait(expression, metadata)
    wgcna_spec = director.create_spec(
        'wgcna', width='double', title='Module–trait relationships', formats=formats,
        annotate_cells=True,
    )

    jobs = {
        'ma': (ma_spec, deg),
        'correlation': (corr_spec, {
            'correlation_matrix': corr, 'sample_labels': samples,
            'groups': metadata['Treatment'],
        }),
        'gsva': (gsva_spec, score_payload),
        'ssgsea': (ssgsea_spec, score_payload),
        'upset': (upset_spec, upset_payload),
        'wgcna': (wgcna_spec, {
            'correlation_matrix': module_trait,
            'module_labels': [f'ME{index + 1}' for index in range(module_trait.shape[0])],
            'trait_labels': trait_names,
        }),
    }
    outputs, reports, signatures = {}, {}, {}
    for key, (spec, data) in jobs.items():
        paths, report, signature = _export(
            director, f'after_{key}', spec, data, output_dir)
        outputs[key], reports[key], signatures[key] = paths, report, signature

    pca_data = {
        'coordinates': np.asarray(pca_adata.obsm['X_pca'])[:, :2],
        'groups': metadata['Treatment'], 'batches': metadata['Genotype'],
        'samples': samples,
        'explained_variance': np.asarray(pca_adata.uns['pca']['variance_ratio'])[:2],
    }
    ora_for_composer = ora.head(8).copy()
    composite_spec = FigureSpec(
        plot_type='composite', width='double', height_mm=172,
        title='Bulk RNA-seq evidence summary', formats=formats,
    )
    panels = [
        FigurePanel('pca', pca_data, director.create_spec(
            'pca', title='Principal component analysis', show_legend=True), row=0, column=0),
        FigurePanel('ma', deg, director.create_spec(
            'ma', title='Differential expression', label_n=4), row=0, column=1),
        FigurePanel('correlation', jobs['correlation'][1], director.create_spec(
            'correlation', title='Sample correlation', max_col_labels=10), row=1, column=0),
        FigurePanel('enrichment', ora_for_composer, director.create_spec(
            'enrichment', title='GO enrichment', top_n=8), row=1, column=1),
    ]
    composite = NatureFigureComposer(director).compose(
        panels, composite_spec, nrows=2, ncols=2, shared_legend=True)
    composite_paths, composite_report = export_figure(
        composite, output_dir / 'bulk_rna_composite', composite_spec,
        report_path=output_dir / 'bulk_rna_composite_nature_readiness.json',
    )
    outputs['composite'] = composite_paths
    reports['composite'] = composite_report
    signatures['composite'] = figure_signature(composite)
    plt.close(composite)

    simulations = {
        deficiency: simulate_image(
            composite_paths['png'], output_dir / f'bulk_rna_composite_{deficiency}.png',
            deficiency,
        )
        for deficiency in ('protanopia', 'deuteranopia', 'tritanopia')
    }
    before_after = _before_after(project_dir, output_dir, outputs)
    manifest = {
        'source_project': str(project_dir.relative_to(ROOT)),
        'data_provenance': {
            'ma': 'real bulk_deg_results.csv',
            'correlation': 'real bulk_heatmap_output.h5ad',
            'gsva_ssgsea': 'deterministic pathway-score fixture derived from real expression',
            'upset': 'real multi-contrast bulk_deg_results*.csv',
            'wgcna': 'module–trait schema fixture derived from real expression and sample metadata',
            'composite': 'real PCA/DEG/correlation/ORA results',
        },
        'outputs': {
            key: {fmt: str(Path(path).relative_to(ROOT)) for fmt, path in paths.items()}
            for key, paths in outputs.items()
        },
        'readiness': {key: report.to_dict() for key, report in reports.items()},
        'visual_snapshots': signatures,
        'color_vision_simulations': {
            key: str(Path(path).relative_to(ROOT)) for key, path in simulations.items()
        },
        'before_after': str(before_after.relative_to(ROOT)) if before_after else None,
    }
    (output_dir / 'validation_manifest.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    lines = [
        '# Nature Figure Engine Phase 2 validation', '',
        '| Template | Source | Score | Gate |', '|---|---|---:|:---:|',
    ]
    for key, report in reports.items():
        source = manifest['data_provenance'].get(key, manifest['data_provenance'].get('gsva_ssgsea', ''))
        lines.append(f"| {key} | {source} | {report.score}/100 | {'pass' if report.ready else 'fail'} |")
    lines.extend([
        '',
        'SVG/PDF are vector-first outputs; PNG/TIFF use the final physical canvas and declared DPI.',
        'Colour-vision previews cover protanopia, deuteranopia and tritanopia.',
        'GSVA/ssGSEA and WGCNA are clearly marked data-derived schema fixtures because the project has no completed result tables for those methods.',
    ])
    (output_dir / 'VALIDATION_REPORT.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project-dir', type=Path, default=DEFAULT_PROJECT)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = validate(args.project_dir, args.output_dir)
    print(json.dumps({key: value['nature_readiness_score']
                      for key, value in manifest['readiness'].items()}, ensure_ascii=False))


if __name__ == '__main__':
    main()
