"""Nature GSEA 与 GO/KEGG enrichment dotplot 模板。"""

import re
import textwrap

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import (
    attach_contract,
    compress_redundant_terms,
    finite_numeric,
    requested_term_rows,
    wrap_term,
    strip_term_id,
)


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


def _bubble_sizes(counts, single=True):
    counts = np.asarray(counts, dtype=float)
    low, high = (16.0, 58.0) if single else (22.0, 82.0)
    if len(counts) == 0 or np.isclose(np.nanmin(counts), np.nanmax(counts)):
        return np.full(len(counts), (low + high) / 2.0)
    transformed = np.sqrt(np.clip(counts, 0, None))
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
        if np.isclose(vmin, vmax):
            vmax = min(1.0, vmin * 10.0)
    return LogNorm(vmin=vmin, vmax=vmax, clip=True)


def _draw_dotplot(frame, spec, *, x_column, x_label, count_column, fdr_column,
                  term_column, zero_line=False, plot_type='enrichment', container=None):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    style = get_style(spec.style)
    profile = style.profile(spec.with_updates(plot_type=plot_type))
    subset = frame.copy()
    subset['_x'] = finite_numeric(subset[x_column])
    subset['_fdr'] = np.clip(finite_numeric(subset[fdr_column], fill=1.0),
                             np.finfo(float).tiny, 1.0)
    subset['_count'] = [_count_value(value) for value in subset[count_column]]
    subset = subset[np.isfinite(subset['_x']) & np.isfinite(subset['_fdr'])].copy()
    if subset.empty:
        raise ValueError('富集 dotplot 没有可绘制的 x/FDR 数值')
    subset['_rank'] = subset['_fdr'].rank(method='first')
    ordered = subset.sort_values(['_fdr', '_x'], ascending=[True, False])
    redundant_removed = 0
    unmatched_targets = []
    selected_targets = []
    target_truncated = 0
    target_mode = bool(spec.pathway_terms) and spec.pathway_selection in {'selected', 'selected_plus_top'}
    if zero_line and target_mode:
        matching = ordered.copy()
        matching['_term'] = matching[term_column].map(str)
        matching['_display_term'] = matching['_term'].map(strip_term_id)
        targets, unmatched_targets, selected_targets = requested_term_rows(
            matching, spec.pathway_terms,
        )
        if len(targets) > spec.top_n:
            target_truncated = len(targets) - spec.top_n
            targets = targets.head(spec.top_n)
            selected_targets = selected_targets[:spec.top_n]
        if spec.pathway_selection == 'selected':
            subset = targets.sort_values('_x', ascending=True)
        else:
            remainder = ordered.loc[~ordered.index.isin(targets.index)]
            positive = remainder[remainder['_x'] >= 0].sort_values(['_fdr', '_x']).head(int(np.ceil(max(0, spec.top_n - len(targets)) / 2)))
            negative = remainder[remainder['_x'] < 0].sort_values(['_fdr', '_x'], ascending=[True, False]).head(int(np.floor(max(0, spec.top_n - len(targets)) / 2)))
            subset = pd.concat([targets, positive, negative]).drop_duplicates()
            subset = subset.sort_values('_x', ascending=True)
        if subset.empty:
            raise ValueError('未匹配到目标通路；请使用结果表中的完整通路名称或 GO/KEGG ID。')
    elif zero_line:
        # Balance positive/negative NES before applying Top N; otherwise a
        # strong positive tail can hide every negatively enriched pathway.
        positive = ordered[ordered['_x'] >= 0].sort_values(['_fdr', '_x']).head(int(np.ceil(spec.top_n / 2)))
        negative = ordered[ordered['_x'] < 0].sort_values(['_fdr', '_x']).head(int(np.floor(spec.top_n / 2)))
        balanced = pd.concat([positive, negative])
        if not balanced.empty:
            subset = balanced.sort_values('_x', ascending=True)
        else:
            subset = ordered.head(spec.top_n)
    else:
        gene_column = _first_column(ordered, ('Genes', 'genes', 'geneID', 'matched_genes', 'lead_genes'))
        if target_mode:
            # Add helper display columns only for matching; source labels and
            # IDs remain untouched in the result table.
            matching = ordered.copy()
            matching['_term'] = matching[term_column].map(str)
            matching['_display_term'] = matching['_term'].map(strip_term_id)
            targets, unmatched_targets, selected_targets = requested_term_rows(
                matching, spec.pathway_terms,
            )
            if len(targets) > spec.top_n:
                target_truncated = len(targets) - spec.top_n
                targets = targets.head(spec.top_n)
                selected_targets = selected_targets[:spec.top_n]
            if spec.pathway_selection == 'selected':
                subset = targets
            else:
                remainder = ordered.loc[~ordered.index.isin(targets.index)]
                remainder = remainder.sort_values(['_fdr', '_x'], ascending=[True, False])
                if gene_column:
                    remainder, redundant_removed = compress_redundant_terms(
                        remainder, gene_column, threshold=spec.redundancy_threshold,
                        maximum=max(0, spec.top_n - len(targets)),
                    )
                else:
                    remainder = remainder.head(max(0, spec.top_n - len(targets)))
                subset = pd.concat([targets, remainder], axis=0)
            subset = subset.sort_values(['_fdr', '_x'], ascending=[False, True])
            if subset.empty:
                raise ValueError('未匹配到目标通路；请使用结果表中的完整通路名称或 GO/KEGG ID。')
        else:
            if gene_column:
                ordered, redundant_removed = compress_redundant_terms(
                    ordered, gene_column, threshold=spec.redundancy_threshold,
                    maximum=spec.top_n,
                )
            else:
                ordered = ordered.head(spec.top_n)
            subset = ordered.sort_values(['_fdr', '_x'], ascending=[False, True])

    terms = [
        wrap_term(value, width=spec.term_wrap, max_chars=spec.term_max_chars)
        for value in subset[term_column]
    ]
    x_values = subset['_x'].to_numpy(dtype=float)
    fdr_values = subset['_fdr'].to_numpy(dtype=float)
    counts = subset['_count'].to_numpy(dtype=float)
    row_steps = 1.0 + 0.28 * np.asarray([str(term).count('\n') for term in terms], dtype=float)
    y = np.cumsum(row_steps) - row_steps / 2.0
    norm = _fdr_norm(fdr_values)
    warnings = []
    semantic_errors = []
    if np.all(counts <= 1):
        warnings.append(
            '所有展示通路的 overlap/count 均为 1；图形可用于探索，但不应作为强富集证据。'
        )
    if np.unique(np.round(x_values, 8)).size == 1:
        warnings.append('所有展示通路的横轴数值相同，dotplot 缺少区分度。')
    if not zero_line and 'Gene ratio' in x_label and np.any(x_values > 1.0 + 1e-9):
        semantic_errors.append(
            'GeneRatio 大于 1；请检查 ORA 输入基因总数与结果表是否来自同一次分析，'
            '禁止将该图作为论文图使用。'
        )
    if len(frame) > spec.top_n:
        if target_mode:
            warnings.append('图中优先展示了用户指定的目标通路；未选择的完整结果保留在来源表。')
        else:
            warnings.append(f'仅展示排序前 {spec.top_n} 个通路；完整结果保留在来源表。')
    if unmatched_targets:
        warnings.append('未匹配目标通路：' + '、'.join(unmatched_targets))
    if target_mode and selected_targets:
        warnings.append(f'已按用户选择展示 {len(selected_targets)} 个目标通路。')
    if target_truncated:
        warnings.append(f'目标通路超过当前图型上限，隐藏 {target_truncated} 个；可提高展示通路数或分批绘图。')
    if redundant_removed:
        warnings.append(f'按命中基因 Jaccard ≥ {spec.redundancy_threshold:.2f} 压缩了 {redundant_removed} 个冗余通路。')

    with style.context(spec):
        if container is None:
            fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
        else:
            fig = container
            ax = container.subplots()
        artist = ax.scatter(
            x_values, y,
            s=_bubble_sizes(counts, single=spec.width == 'single'),
            c=fdr_values, cmap=style.fdr_cmap(), norm=norm,
            edgecolors='white', linewidths=0.45, alpha=0.95, zorder=3,
        )
        if zero_line:
            ax.axvline(0, color=style.neutral_dark, linestyle=(0, (3, 2)),
                       linewidth=0.55, zorder=0)
        ax.set_yticks(y, terms)
        ax.tick_params(axis='y', length=0, pad=3, labelsize=max(5.5, profile.tick_font_pt - 0.2))
        ax.set_xlabel(x_label)
        ax.set_ylabel('')
        title = spec.title or ('GSEA' if zero_line else 'Pathway enrichment')
        # The single-column dotplot reserves roughly half of the canvas for
        # pathway labels, so its title starts at the plot-column boundary.
        # Wrapping at 30 characters keeps the title inside the 89 mm canvas;
        # the previous 40-character wrap visibly clipped long comparison
        # labels such as ``NH4Cl_B vs NH4Cl_En · GO_BP ORA``.
        title_width = 30 if spec.width == 'single' else 64
        ax.set_title(textwrap.fill(str(title), width=title_width, break_long_words=False),
                     loc='left', pad=5)
        style.apply_axis(ax, profile)
        ax.spines['left'].set_visible(False)
        # The y-axis is a label column, not a quantitative axis.  Applying the
        # shared style resets tick lengths, so explicitly remove the small
        # dashes that otherwise appear after every pathway name.
        ax.tick_params(axis='y', length=0, pad=3)
        ax.margins(x=0.09, y=0.04)

        left = 0.50 if spec.width == 'single' else 0.40
        panel_mode = container is not None
        # Reserve two independent rows below the x-axis: one for the FDR
        # colourbar and one for the Count size key.  A figure-level legend
        # anchored at the bottom-right is tempting here, but its title and
        # marker row can be pulled into the colourbar by ``bbox_inches='tight'``
        # (this was the source of the overlapping ``Count``/``FDR`` labels in
        # Figure Studio previews).  Dedicated axes keep the two encodings
        # geometrically separated in PNG, SVG and PDF exports.
        fig.subplots_adjust(left=left, right=0.96, top=0.90,
                            bottom=0.34 if panel_mode else 0.31)
        cbar_ax = fig.add_axes([
            left, 0.235 if panel_mode else 0.190,
            min(0.30, 0.94 - left), 0.018,
        ])
        cbar_ax._nature_auxiliary = True
        colorbar = fig.colorbar(artist, cax=cbar_ax, orientation='horizontal')
        colorbar.outline.set_visible(False)
        colorbar.ax.tick_params(
            labelsize=max(5.2, profile.tick_font_pt - 0.7),
            length=1.8, width=profile.tick_width_pt, colors=style.axis,
        )
        tick_values = np.geomspace(norm.vmin, norm.vmax, 3)
        colorbar.set_ticks(tick_values)
        colorbar.minorticks_off()
        colorbar.ax.set_xticklabels([f'{value:.1g}' for value in tick_values])
        colorbar.set_label('FDR', fontsize=profile.legend_font_pt, labelpad=2, color=style.text)

        if np.unique(counts).size > 1:
            size_values = np.unique(np.quantile(counts, [0.0, 0.5, 1.0]).round())
            size_values = size_values[size_values > 0]
            legend_sizes = _bubble_sizes(size_values, spec.width == 'single')
            size_handles = [
                Line2D([], [], linestyle='None', marker='o',
                       markersize=np.sqrt(marker_size) / 1.5,
                       markerfacecolor=style.signal_teal, markeredgecolor='white',
                       markeredgewidth=0.4, label=f'{int(value)}')
                for value, marker_size in zip(size_values, legend_sizes)
            ]
            count_ax = fig.add_axes([
                left, 0.025 if panel_mode else 0.012,
                min(0.38, 0.94 - left), 0.065,
            ])
            count_ax._nature_auxiliary = True
            count_ax.set_axis_off()
            count_ax.legend(
                handles=size_handles, title='Count', loc='center',
                frameon=False, borderaxespad=0.0, handletextpad=0.25,
                ncol=min(3, len(size_handles)), columnspacing=0.55,
                fontsize=profile.legend_font_pt, title_fontsize=profile.legend_font_pt,
                labelcolor=style.text,
            )
            if count_ax.legend_ is not None:
                count_ax.legend_.get_title().set_color(style.text)

    if container is not None:
        return container
    return attach_contract(
        fig, spec, style,
        semantic_warnings=warnings,
        semantic_errors=semantic_errors,
        encodings={'x': x_label, 'size': 'Count', 'color': 'FDR'},
    )


class NatureEnrichmentDotplot:
    """GO/KEGG ORA 默认 dotplot：GeneRatio、Count、FDR。"""

    plot_type = 'enrichment_dotplot'

    def render(self, data, spec: FigureSpec, container=None):
        frame = pd.DataFrame(data).copy()
        term_column = _first_column(frame, ('Term', 'Description', 'pathway'))
        fdr_column = _first_column(frame, (
            'Adjusted P-value', 'Adjusted p-value', 'FDR', 'padj', 'qvalue', 'p.adjust',
        ))
        ratio_column = _first_column(frame, ('GeneRatio',))
        overlap_fallback = False
        if ratio_column is None:
            ratio_column = _first_column(frame, ('fraction', 'Overlap'))
            overlap_fallback = ratio_column is not None
        count_column = _first_column(frame, ('Count', 'num', 'Gene Count', 'Overlap', 'setSize'))
        if not all((term_column, fdr_column, ratio_column, count_column)):
            raise ValueError('NatureEnrichmentDotplot 需要 term、FDR、GeneRatio/Overlap 和 Count 字段')
        frame['_gene_ratio'] = [_ratio_value(value) for value in frame[ratio_column]]
        return _draw_dotplot(
            frame, spec.with_updates(plot_type='enrichment_dotplot'),
            x_column='_gene_ratio',
            x_label='Overlap fraction' if overlap_fallback else 'Gene ratio',
            count_column=count_column,
            fdr_column=fdr_column, term_column=term_column, zero_line=False,
            plot_type='enrichment_dotplot', container=container,
        )


class NatureGSEA:
    """GSEA dotplot：NES、gene-set size/Count、FDR。"""

    plot_type = 'gsea'

    def render(self, data, spec: FigureSpec, container=None):
        frame = pd.DataFrame(data).copy()
        term_column = _first_column(frame, ('Term', 'Description', 'pathway'))
        nes_column = _first_column(frame, ('NES', 'normalized enrichment score'))
        fdr_column = _first_column(frame, ('fdr', 'FDR q-val', 'FDR', 'qvalue', 'padj'))
        count_column = _first_column(frame, ('matched_size', 'setSize', 'size', 'Count', 'Tag %'))
        if count_column and frame[count_column].dtype == object:
            frame['_gsea_count'] = [
                _count_value(re.match(r'\s*([0-9.]+)', str(value)).group(1)
                             if re.match(r'\s*([0-9.]+)', str(value)) else value)
                for value in frame[count_column]
            ]
            count_column = '_gsea_count'
        if count_column is None:
            frame['_gsea_count'] = 1
            count_column = '_gsea_count'
        if not all((term_column, nes_column, fdr_column)):
            raise ValueError('NatureGSEA 需要 Term、NES 和 FDR 字段')
        return _draw_dotplot(
            frame, spec.with_updates(plot_type='gsea'),
            x_column=nes_column, x_label='Normalized enrichment score (NES)',
            count_column=count_column, fdr_column=fdr_column,
            term_column=term_column, zero_line=True, plot_type='gsea', container=container,
        )
