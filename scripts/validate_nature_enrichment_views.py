#!/usr/bin/env python
"""Render and validate the pathway-enrichment forms on project data.

The ORA panels use the stored GO-BP result table from project ``bec48c50-bc9``.
The running-score panel uses the project's real ranked DEG statistics and the
stored pathway membership to exercise the lossless GSEA running-data contract;
because no completed GSEA task is present in this project, NES/FDR are left
blank rather than invented.
"""

from __future__ import annotations

import html
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mpl-nature-enrichment')
BASE = Path(__file__).resolve().parents[1]
PROJECT = BASE / 'data/projects/bec48c50-bc9'
OUT = BASE / 'artifacts/nature_enrichment_views'
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(BASE))

from figure_engine import NatureFigureDirector, export_figure


def main():
    ora_path = PROJECT / 'results/enrichment_ora_go_bp_hep_i_ffa_vs_hep_ct_ffa_results.csv'
    deg_path = PROJECT / 'results/bulk_deg_results.csv'
    ora = pd.read_csv(ora_path)
    deg = pd.read_csv(deg_path).dropna(subset=['gene', 'log2FC']).copy()
    deg = deg.sort_values('log2FC', ascending=False).drop_duplicates('gene').reset_index(drop=True)
    # The stored ORA table predates the current DEG CSV and its manifest records
    # the exact input-list size used for that task.  Prefer that provenance over
    # inferring a denominator from a possibly newer DEG result.
    manifest = PROJECT / 'results/manifests/bulk_enrichment_16a805cd-7c7_analysis_manifest.json'
    manifest_count = 0
    if manifest.is_file():
        manifest_count = int(json.loads(manifest.read_text(encoding='utf-8')).get('summary', {}).get('n_input_genes', 0) or 0)
    input_count = manifest_count or max(1, int((deg.get('regulation', pd.Series(dtype=str)).astype(str) != 'NS').sum()))
    counts = pd.to_numeric(ora['Overlap'].astype(str).str.split('/', n=1).str[0], errors='coerce')
    ora['GeneRatio'] = counts / float(input_count)
    ora['InputGeneCount'] = input_count
    ora = ora.sort_values('Adjusted P-value').head(30).reset_index(drop=True)
    ora.to_csv(OUT / 'source_table_with_gene_ratio.csv', index=False)

    director = NatureFigureDirector()
    formats = ('svg', 'pdf', 'png', 'tiff')
    rendered = []
    reports = {}

    views = [
        ('01_dotplot', 'enrichment', ora, 'ORA · Nature enrichment dotplot'),
        ('02_barplot', 'barplot', ora, 'ORA · Nature enrichment barplot'),
        ('03_chord', 'chord', ora, 'ORA · Nature chord-style membership'),
        ('04_cnetplot', 'cnetplot', ora, 'ORA · Nature pathway–gene cnetplot'),
        ('05_emapplot', 'emapplot', ora, 'ORA · Nature pathway similarity map'),
    ]
    for stem, plot_type, data, title in views:
        spec = director.create_spec(
            plot_type, width='single', top_n=12, max_genes=30,
            formats=formats, title=title,
        )
        figure = director.render(spec, data)
        paths, report = export_figure(figure, OUT / stem, spec)
        report.write_json(OUT / f'{stem}_nature_readiness.json')
        rendered.append((stem, title, paths, report))
        reports[stem] = report.to_dict()
        import matplotlib.pyplot as plt
        plt.close(figure)

    gene_to_rank = {str(gene): index for index, gene in enumerate(deg['gene'].astype(str))}
    curves = []
    for _, row in ora.head(2).iterrows():
        genes = [gene for gene in str(row.get('Genes', '')).split(';') if gene in gene_to_rank]
        if len(genes) < 2:
            continue
        curves.append({
            'term': str(row['Term']),
            'ranking': deg['log2FC'].to_numpy(float),
            'hit_indices': np.array(sorted({gene_to_rank[gene] for gene in genes}), dtype=int),
            # No completed GSEA task exists in this project, so leave these
            # inferential values absent and let the renderer show n/a.
        })
    if curves:
        spec = director.create_spec(
            'gsea_running', width='single', running_term_n=len(curves),
            formats=formats, title='GSEA running score · ranking contract validation',
        )
        figure = director.render(spec, {'curves': curves})
        paths, report = export_figure(figure, OUT / '06_gsea_running', spec)
        report.write_json(OUT / '06_gsea_running_nature_readiness.json')
        rendered.append(('06_gsea_running', 'GSEA · Nature running score', paths, report))
        reports['06_gsea_running'] = report.to_dict()
        import matplotlib.pyplot as plt
        plt.close(figure)

    cards = []
    for stem, title, paths, report in rendered:
        cards.append(
            '<article><h2>{}</h2><img src="{}" alt="{}"><p>Nature readiness: {} / 100 · {} · '
            '<a href="{}">SVG</a> · <a href="{}">PDF</a> · <a href="{}">TIFF</a> · '
            '<a href="{}_nature_readiness.json">QA JSON</a></p></article>'.format(
                html.escape(title), Path(paths['png']).name, html.escape(title), report.score,
                'pass' if report.ready else 'warning', Path(paths['svg']).name,
                Path(paths['pdf']).name, Path(paths['tiff']).name, stem,
            )
        )
    note = (
        'ORA panels use the stored project result table. The running-score panel uses the real project '
        'ranked DEG statistics and stored pathway membership to validate the rendering contract; this '
        'project has no completed GSEA task, so NES/FDR are intentionally shown as n/a. '
        f'GeneRatio uses Count / InputGeneCount with InputGeneCount={input_count}.'
    )
    page = '<!doctype html><meta charset="utf-8"><title>Nature enrichment views</title><style>'
    page += 'body{font-family:Arial,sans-serif;color:#20262e;max-width:1100px;margin:24px auto} '
    page += 'article{margin:28px 0;border-top:1px solid #d8dce2;padding-top:16px} '
    page += 'img{max-width:100%;height:auto;border:1px solid #e6e9ed} p{color:#626c78}'
    page += '</style><h1>Nature Figure Engine · pathway enrichment</h1>'
    page += f'<p>{html.escape(note)} <a href="source_table_with_gene_ratio.csv">derived QA table</a>.</p>' + ''.join(cards)
    (OUT / 'index.html').write_text(page, encoding='utf-8')
    (OUT / 'validation_summary.json').write_text(json.dumps(reports, ensure_ascii=False, indent=2), encoding='utf-8')
    for stem, _, _, report in rendered:
        print(stem, report.score, report.status)
    print('output:', OUT)


if __name__ == '__main__':
    main()
