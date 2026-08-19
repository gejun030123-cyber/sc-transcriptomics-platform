"""多数据库富集概览：GO/KEGG 等数据库分面而不混合语义。"""

from __future__ import annotations

import re
import textwrap

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import attach_contract, compress_redundant_terms, requested_term_rows, wrap_term


def _first_column(frame, candidates):
    lookup = {str(column).strip().lower(): column for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def _ratio_value(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    if isinstance(value, str) and '/' in value:
        left, right = value.split('/', 1)
        try:
            denominator = float(right)
            return float(left) / denominator if denominator else np.nan
        except ValueError:
            return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _count_value(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return 1.0
    if isinstance(value, str) and '/' in value:
        value = value.split('/', 1)[0]
    try:
        return max(1.0, float(value))
    except (TypeError, ValueError):
        return 1.0


def _term_id(value):
    match = re.search(r'\((?:GO|KEGG|REACTOME|WP|WIKIPATHWAYS)\s*:[^)]+\)', str(value), re.I)
    if match:
        return match.group(0).strip('()')
    match = re.search(r'\b(?:GO|KEGG|REACTOME|WP|WIKIPATHWAYS)\s*:[A-Za-z0-9_.-]+\b', str(value), re.I)
    return match.group(0) if match else ''


def _database_label(value):
    text = str(value or '').upper().replace('-', '_')
    labels = {
        'GO_BP': 'GO-BP', 'GO_BIOLOGICAL_PROCESS': 'GO-BP', 'BP': 'GO-BP',
        'GO_CC': 'GO-CC', 'GO_CELLULAR_COMPONENT': 'GO-CC', 'CC': 'GO-CC',
        'GO_MF': 'GO-MF', 'GO_MOLECULAR_FUNCTION': 'GO-MF', 'MF': 'GO-MF',
        'KEGG': 'KEGG', 'REACTOME': 'Reactome',
        'WIKIPATHWAYS': 'WikiPathways', 'WIKIPATHWAY': 'WikiPathways',
    }
    return labels.get(text, str(value or 'Other'))


def _bubble_sizes(values, low=24.0, high=92.0):
    values = np.asarray(values, dtype=float)
    if len(values) == 0 or np.isclose(np.nanmin(values), np.nanmax(values)):
        return np.full(len(values), (low + high) / 2.0)
    transformed = np.sqrt(np.clip(values, 0, None))
    scaled = (transformed - transformed.min()) / (transformed.max() - transformed.min())
    return low + scaled * (high - low)


def _fdr_norm(values):
    from matplotlib.colors import LogNorm

    values = np.clip(np.asarray(values, dtype=float), np.finfo(float).tiny, 1.0)
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    if np.isclose(vmin, vmax):
        vmin = max(np.finfo(float).tiny, vmin / 2.0)
        vmax = min(1.0, vmax * 2.0)
    return LogNorm(vmin=vmin, vmax=vmax, clip=True)


class NatureEnrichmentOverview:
    """双栏多数据库概览；每个数据库保留独立的 term 语义。"""

    plot_type = 'enrichment_overview'

    def render(self, data, spec: FigureSpec, container=None):
        if container is not None:
            raise ValueError('enrichment_overview 需要独立 Figure，不支持嵌入单一 panel')
        frame = pd.DataFrame(data).copy()
        database_column = _first_column(frame, ('Database', 'database', 'Gene_set', 'gene_set'))
        term_column = _first_column(frame, ('Term', 'Description', 'pathway', 'term'))
        fdr_column = _first_column(frame, (
            'Enrichment FDR', 'Adjusted P-value', 'Adjusted p-value', 'FDR', 'fdr', 'padj',
        ))
        method_column = _first_column(frame, ('Method', 'method'))
        if not database_column or not term_column:
            raise ValueError('enrichment_overview 需要 Database 和 Term 字段')
        if not fdr_column:
            raise ValueError('enrichment_overview 需要 FDR/Adjusted P-value 字段')

        method = str(frame[method_column].dropna().iloc[0]).upper() if method_column and frame[method_column].notna().any() else 'ORA'
        is_gsea = method == 'GSEA'
        score_column = _first_column(frame, ('NES', 'nes', 'normalized enrichment score')) if is_gsea else None
        ratio_column = _first_column(frame, ('GeneRatio', 'gene_ratio', 'fraction', 'Overlap'))
        count_column = _first_column(frame, ('Count', 'num', 'Gene Count', 'gene_count', 'Overlap', 'setSize', 'size'))
        if is_gsea and not score_column:
            raise ValueError('GSEA overview 需要 NES 字段')
        if not is_gsea and not ratio_column:
            raise ValueError('ORA overview 需要 GeneRatio 或 Overlap 字段')

        frame['_database'] = frame[database_column].map(_database_label)
        frame['_term'] = frame[term_column].map(str)
        frame['_display_term'] = frame['_term'].map(lambda value: re.sub(r'\s+', ' ', value).strip())
        frame['_fdr'] = pd.to_numeric(frame[fdr_column], errors='coerce').clip(lower=np.finfo(float).tiny, upper=1.0)
        frame['_x'] = (pd.to_numeric(frame[score_column], errors='coerce') if is_gsea
                       else frame[ratio_column].map(_ratio_value))
        frame['_count'] = frame[count_column].map(_count_value) if count_column else 1.0
        frame = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=['_fdr', '_x']).copy()
        if frame.empty:
            raise ValueError('没有可绘制的多数据库富集结果')

        preferred = ['GO-BP', 'GO-CC', 'GO-MF', 'KEGG', 'Reactome', 'WikiPathways', 'Other']
        available = [item for item in preferred if item in set(frame['_database'])]
        available.extend(item for item in frame['_database'].drop_duplicates() if item not in available)
        if spec.database_scope:
            requested = {_database_label(value) for value in spec.database_scope}
            available = [item for item in available if item in requested]
        if not available:
            raise ValueError('没有匹配到可展示的数据库')

        selected_frames = []
        warnings = []
        compressed_databases = []
        target_mode = bool(spec.pathway_terms) and spec.pathway_selection in {'selected', 'selected_plus_top'}
        for database in available:
            subset = frame[frame['_database'] == database].sort_values(['_fdr', '_x'], ascending=[True, False]).copy()
            target_rows = subset.iloc[0:0].copy()
            if target_mode:
                target_rows, unmatched, _ = requested_term_rows(subset, spec.pathway_terms)
                if unmatched:
                    warnings.append(f'{database} 未匹配目标通路：' + '、'.join(unmatched))
                if spec.pathway_selection == 'selected':
                    subset = target_rows
                else:
                    remainder = subset.loc[~subset.index.isin(target_rows.index)]
                    subset = pd.concat([target_rows, remainder], axis=0)
            gene_column = _first_column(subset, ('Genes', 'genes', 'geneID', 'matched_genes', 'lead_genes'))
            if gene_column:
                subset, removed = compress_redundant_terms(
                    subset, gene_column, threshold=spec.redundancy_threshold, maximum=spec.top_n,
                )
                if removed:
                    compressed_databases.append(f'{database} ({removed})')
            else:
                subset = subset.head(spec.top_n)
            if len(subset) > spec.top_n:
                subset = subset.head(spec.top_n)
            if subset.empty:
                continue
            selected_frames.append(subset)

        if not selected_frames:
            raise ValueError('目标通路选择后没有可绘制结果')
        if compressed_databases:
            warnings.append('按数据库内命中基因 Jaccard 压缩冗余通路：' + '、'.join(compressed_databases))
        display = pd.concat(selected_frames, ignore_index=True)
        profile_style = get_style(spec.style)
        profile = profile_style.profile(spec)
        panel_count = len(selected_frames)
        ncols = 2 if panel_count > 1 else 1
        nrows = int(np.ceil(panel_count / ncols))
        with profile_style.context(spec):
            import matplotlib.pyplot as plt
            from matplotlib.lines import Line2D

            fig, axes = plt.subplots(
                nrows, ncols, figsize=profile.figsize, dpi=profile.dpi,
                squeeze=False, sharex=True if not is_gsea else False,
            )
            axes_flat = axes.ravel()
            norm = _fdr_norm(display['_fdr'].to_numpy(float))
            colors = profile_style.fdr_cmap()
            x_label = 'Normalized enrichment score (NES)' if is_gsea else 'Gene ratio'
            letters = 'abcdefghijklmnopqrstuvwxyz'
            max_x = float(np.nanmax(display['_x']))
            min_x = float(np.nanmin(display['_x']))
            for index, (axis, subset) in enumerate(zip(axes_flat, selected_frames)):
                local = subset.copy().reset_index(drop=True)
                labels = [wrap_term(value, width=25 if ncols == 2 else 34, max_chars=72) for value in local['_term']]
                positions = np.arange(len(local))
                point_sizes = _bubble_sizes(local['_count'].to_numpy(float), low=28 if spec.width == 'double' else 20, high=100 if spec.width == 'double' else 72)
                artist = axis.scatter(
                    local['_x'], positions, s=point_sizes, c=local['_fdr'], cmap=colors,
                    norm=norm, edgecolors='white', linewidths=.45, alpha=.96, zorder=3,
                )
                axis.set_yticks(positions, labels)
                axis.invert_yaxis()
                axis.set_title(_database_label(local['_database'].iloc[0]), loc='left', pad=6)
                axis.text(-0.16, 1.04, letters[index], transform=axis.transAxes,
                          ha='left', va='bottom', fontweight='bold', color=profile_style.text)
                axis.tick_params(axis='y', length=0, pad=3, labelsize=max(5.2, profile.tick_font_pt - .2))
                profile_style.apply_axis(axis, profile)
                axis.spines['left'].set_visible(False)
                axis.tick_params(axis='y', length=0, pad=3)
                axis.margins(x=.10, y=.10)
                if is_gsea:
                    axis.axvline(0, color=profile_style.neutral_dark, linestyle=(0, (3, 2)), linewidth=.55, zorder=0)
            for axis in axes_flat[panel_count:]:
                axis.set_visible(False)
            for axis in axes_flat[:panel_count]:
                axis.set_xlabel(x_label)
            fig.subplots_adjust(left=.28 if ncols == 2 else .42, right=.97, top=.88, bottom=.19,
                                wspace=.52, hspace=.62)
            cbar_ax = fig.add_axes([.38, .075, .28, .018])
            cbar_ax._nature_auxiliary = True
            colorbar = fig.colorbar(artist, cax=cbar_ax, orientation='horizontal')
            colorbar.outline.set_visible(False)
            colorbar.ax.tick_params(labelsize=max(5.2, profile.tick_font_pt - .7), length=1.8,
                                    width=profile.tick_width_pt, colors=profile_style.axis)
            tick_values = np.geomspace(norm.vmin, norm.vmax, 3)
            colorbar.set_ticks(tick_values)
            colorbar.minorticks_off()
            colorbar.ax.set_xticklabels([f'{value:.1g}' for value in tick_values])
            colorbar.set_label('FDR', fontsize=profile.legend_font_pt, labelpad=2, color=profile_style.text)
            if np.unique(display['_count']).size > 1:
                legend_values = np.unique(np.quantile(display['_count'], [0, .5, 1]).round())
                legend_values = legend_values[legend_values > 0]
                legend_sizes = _bubble_sizes(legend_values, low=28, high=100)
                handles = [Line2D([], [], linestyle='None', marker='o',
                                   markersize=np.sqrt(size) / 1.55,
                                   markerfacecolor=profile_style.signal_teal,
                                   markeredgecolor='white', markeredgewidth=.4,
                                   label=str(int(count)))
                           for count, size in zip(legend_values, legend_sizes)]
                fig.legend(handles=handles, title='Count', loc='lower right',
                           bbox_to_anchor=(.97, .025), frameon=False,
                           ncol=min(3, len(handles)), columnspacing=.55,
                           handletextpad=.25, fontsize=profile.legend_font_pt,
                           title_fontsize=profile.legend_font_pt)
            title = spec.title or 'Multi-database pathway enrichment overview'
            fig.suptitle(textwrap.fill(str(title), width=80 if spec.width == 'double' else 42),
                         x=.02, ha='left', y=.965, fontsize=profile.title_font_pt,
                         fontweight='semibold', color=profile_style.text)
            fig._nature_panel_grid = {'nrows': nrows, 'ncols': ncols, 'panels': panel_count}
            return attach_contract(
                fig, spec, profile_style,
                semantic_warnings=warnings,
                encodings={'x': x_label, 'size': 'Count', 'color': 'FDR', 'facet': 'Database'},
            )
