"""Publication-grade pathway enrichment views.

The templates in this module deliberately accept a small, explicit data
contract.  They never infer gene membership from a summary statistic: Chord,
Cnet and Emap require a real ``Genes``/``geneID`` field, while the GSEA running
score requires ``ranking`` plus ``hit_indices``/``RES``.  All five views use the
same Nature style profile and return the normal figure contract consumed by
the validator and exporter.
"""

from __future__ import annotations

import ast
import math
import re
from collections import Counter

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import (
    attach_contract,
    compress_redundant_terms,
    gene_label_order,
    requested_term_rows,
    strip_term_id,
    wrap_term,
)


def _first_column(frame, candidates):
    lookup = {str(column).strip().lower(): column for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def _number(value, default=np.nan):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return default
    if isinstance(value, str) and '/' in value:
        value = value.split('/', 1)[0]
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _genes(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)):
        values = list(value)
    else:
        text = str(value).strip()
        if not text or text.lower() in {'nan', 'none'}:
            return []
        if text.startswith('[') and text.endswith(']'):
            try:
                values = ast.literal_eval(text)
            except (SyntaxError, ValueError):
                values = re.split(r'[;,/|]', text.strip('[]'))
        else:
            values = re.split(r'[;,/|]', text)
    return list(dict.fromkeys(str(item).strip() for item in values if str(item).strip()))


def _enrichment_frame(data, spec: FigureSpec, *, require_genes=False):
    frame = pd.DataFrame(data).copy()
    if frame.empty:
        raise ValueError('富集结果为空，无法生成通路图。')
    term_column = _first_column(frame, ('Term', 'Description', 'pathway', 'term'))
    fdr_column = _first_column(frame, (
        'Adjusted P-value', 'Adjusted p-value', 'FDR q-val', 'FDR',
        'fdr', 'qvalue', 'padj', 'p.adjust',
    ))
    count_column = _first_column(frame, (
        'Count', 'count', 'num', 'Gene Count', 'gene_count',
        'matched_size', 'setSize', 'geneset_size', 'Overlap',
    ))
    gene_column = _first_column(frame, (
        'Genes', 'genes', 'geneID', 'gene_id', 'matched_genes',
        'lead_genes', 'Lead_genes',
    ))
    score_column = _first_column(frame, ('NES', 'nes', 'normalized enrichment score'))
    if term_column is None or fdr_column is None:
        raise ValueError('富集模板需要 Term/Description 和 FDR/Adjusted P-value 字段。')
    if require_genes and gene_column is None:
        raise ValueError('该图型需要真实基因成员字段（Genes/geneID/matched_genes）。')
    frame['_term'] = frame[term_column].map(str)
    # Keep the database-qualified term in the source table, but never put GO/
    # KEGG identifiers into the primary figure labels.
    frame['_display_term'] = frame['_term'].map(strip_term_id)
    frame['_fdr'] = np.clip(pd.to_numeric(frame[fdr_column], errors='coerce'),
                            np.finfo(float).tiny, 1.0)
    frame['_count'] = frame[count_column].map(_number) if count_column else pd.Series(1.0, index=frame.index)
    frame['_count'] = pd.to_numeric(frame['_count'], errors='coerce').fillna(1).clip(lower=1)
    frame['_genes'] = frame[gene_column].map(_genes) if gene_column else [[] for _ in range(len(frame))]
    if score_column:
        frame['_nes'] = pd.to_numeric(frame[score_column], errors='coerce')
    return frame.dropna(subset=['_fdr']).copy()


def _new_figure(spec, style, container=None):
    import matplotlib.pyplot as plt

    profile = style.profile(spec)
    if container is None:
        return plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
    return container, container.subplots()


def _term_limits(spec):
    return min(spec.top_n, 8 if spec.width == 'single' else 12)


def _network_colors():
    # This fixed subset remains distinguishable under the engine's colour
    # vision audit; pathway names and node shapes provide redundant labels.
    return ('#355F8A', '#90AEC5', '#7A6FA8', '#D08A5B', '#6F9D72', '#4C9099')


def _select_terms(frame, spec, *, direction=False):
    limit = _term_limits(spec)
    target_mode = bool(spec.pathway_terms) and spec.pathway_selection in {'selected', 'selected_plus_top'}
    unmatched = []
    requested = []
    if target_mode:
        requested_frame, unmatched, requested = requested_term_rows(frame, spec.pathway_terms)
        if requested_frame.empty:
            raise ValueError('未匹配到目标通路；请使用结果表中的完整通路名称或 GO/KEGG ID。')
        target_truncated = max(0, len(requested_frame) - limit)
        if target_truncated:
            requested_frame = requested_frame.head(limit)
            requested = requested[:limit]
        if spec.pathway_selection == 'selected':
            selected = requested_frame
            selected.attrs['redundant_removed'] = 0
            selected.attrs['target_selected'] = len(requested)
            selected.attrs['target_unmatched'] = unmatched
            selected.attrs['target_truncated'] = target_truncated
            return selected.reset_index(drop=True)
    else:
        requested_frame = frame.iloc[0:0].copy()
    if direction and '_nes' in frame:
        source = frame.loc[~frame.index.isin(requested_frame.index)] if target_mode else frame
        remaining = max(0, limit - len(requested_frame))
        positive = source[source['_nes'] >= 0].sort_values(['_fdr', '_nes']).head(math.ceil(remaining / 2))
        negative = source[source['_nes'] < 0].sort_values(['_fdr', '_nes'], ascending=[True, False]).head(math.floor(remaining / 2))
        selected = pd.concat([requested_frame, negative, positive], ignore_index=True)
        selected = selected.sort_values('_nes').reset_index(drop=True)
        selected.attrs['redundant_removed'] = 0
        selected.attrs['target_selected'] = len(requested)
        selected.attrs['target_unmatched'] = unmatched
        selected.attrs['target_truncated'] = target_truncated
        return selected
    candidates = frame.loc[~frame.index.isin(requested_frame.index)] if target_mode else frame
    candidates = candidates.sort_values(['_fdr', '_count'], ascending=[True, False])
    selected, removed = compress_redundant_terms(
        candidates, '_genes', threshold=spec.redundancy_threshold,
        maximum=max(0, limit - len(requested_frame)),
    )
    if target_mode:
        selected = pd.concat([requested_frame, selected], ignore_index=False)
    # Keep the most significant representative terms first.  Network views
    # apply their physical term cap after this selection; descending FDR here
    # would silently retain the least significant terms in a single-column
    # Chord/Cnet figure.
    selected = selected.sort_values('_fdr', ascending=True).reset_index(drop=True)
    selected.attrs['redundant_removed'] = removed
    selected.attrs['target_selected'] = len(requested)
    selected.attrs['target_unmatched'] = unmatched
    selected.attrs['target_truncated'] = target_truncated if target_mode else 0
    return selected


class NatureEnrichmentBarplot:
    """Horizontal ORA/GSEA bars with explicit score and gene-count labels."""

    plot_type = 'enrichment_barplot'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt

        style = get_style(spec.style)
        frame = _enrichment_frame(data, spec)
        is_gsea = '_nes' in frame and frame['_nes'].notna().any()
        subset = _select_terms(frame, spec, direction=is_gsea)
        values = subset['_nes'].to_numpy(float) if is_gsea else -np.log10(subset['_fdr'].to_numpy(float))
        colors = [style.signal_red if value >= 0 else style.signal_blue for value in values] if is_gsea else [style.signal_blue] * len(values)
        terms = [wrap_term(value, width=spec.term_wrap, max_chars=spec.term_max_chars) for value in subset['_display_term']]
        row_steps = 1.0 + 0.28 * np.asarray([str(term).count('\n') for term in terms], dtype=float)
        y = np.cumsum(row_steps) - row_steps / 2.0
        with style.context(spec):
            fig, ax = _new_figure(spec.with_updates(plot_type=self.plot_type), style, container)
            ax.barh(y, values, color=colors, height=row_steps * .64, edgecolor='none', alpha=.96, zorder=2)
            ax.axvline(0, color=style.neutral_dark, linewidth=.55, zorder=1) if is_gsea else None
            ax.set_yticks(y, terms)
            ax.invert_yaxis()
            ax.set_xlabel('Normalized enrichment score (NES)' if is_gsea else r'$-\log_{10}$(adjusted P-value)')
            ax.set_ylabel('')
            ax.set_title(spec.title or ('GSEA enrichment' if is_gsea else 'ORA enrichment'), loc='left', pad=5)
            style.apply_axis(ax, style.profile(spec.with_updates(plot_type=self.plot_type)))
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis='y', length=0, pad=3)
            limit = max(float(np.nanmax(np.abs(values))), 1.0)
            # Reserve deliberate space for the audit labels instead of letting
            # the text collide with the right edge at a single-column width.
            ax.set_xlim((-limit * 1.48, limit * 1.48) if is_gsea else (0, limit * 1.48))
            for yi, value, count in zip(y, values, subset['_count']):
                offset = limit * .035
                x = value + (offset if value >= 0 else -offset)
                value_label = f'{value:+.2f}' if is_gsea else f'{value:.2f}'
                ax.text(x, yi, f'{value_label}  ·  n={int(round(count))}',
                        ha='left' if value >= 0 else 'right',
                        va='center', color=style.muted_text, fontsize=style.profile(spec).legend_font_pt)
            fig.subplots_adjust(left=.51 if spec.width == 'single' else .40, right=.97, top=.88, bottom=.20)
        warnings = []
        if len(frame) > len(subset):
            warnings.append('图中仅显示当前选择的通路子集；完整结果保留在结果表。')
        if subset.attrs.get('target_selected'):
            warnings.append(f"已按用户选择展示 {subset.attrs['target_selected']} 个目标通路。")
        if subset.attrs.get('target_unmatched'):
            warnings.append('未匹配目标通路：' + '、'.join(subset.attrs['target_unmatched']))
        if subset.attrs.get('target_truncated'):
            warnings.append(f"目标通路超过当前图型上限，隐藏 {subset.attrs['target_truncated']} 个；可提高展示通路数或分批绘图。")
        if subset.attrs.get('redundant_removed', 0):
            warnings.append(
                f"按命中基因 Jaccard ≥ {spec.redundancy_threshold:.2f} "
                f"压缩了 {subset.attrs['redundant_removed']} 个冗余通路。"
            )
        return attach_contract(fig, spec.with_updates(plot_type=self.plot_type), style,
                               semantic_warnings=warnings,
                               encodings={'x': 'NES' if is_gsea else '-log10(FDR)', 'color': 'NES direction' if is_gsea else 'single significance colour', 'size': 'Count annotation'})


def _membership_subset(frame, spec, term_cap, gene_cap):
    selected_all = _select_terms(frame, spec, direction=False)
    if spec.gene_labels:
        wanted = set(spec.gene_labels)
        hits = selected_all[selected_all['_genes'].map(lambda values: bool(wanted.intersection(values)))]
        remainder = selected_all.loc[~selected_all.index.isin(hits.index)]
        if not hits.empty:
            selected_all = pd.concat([hits, remainder], ignore_index=False)
    selected = selected_all.head(term_cap).copy()
    selected.attrs.update(selected_all.attrs)
    selected = selected.reset_index(drop=True)
    if len(selected_all) > term_cap and spec.pathway_terms and spec.pathway_selection in {'selected', 'selected_plus_top'}:
        selected.attrs['target_truncated'] = len(selected_all) - term_cap
    frequency = Counter(gene for genes in selected['_genes'] for gene in genes)
    ranked_genes = sorted(frequency, key=lambda item: (-frequency[item], item))
    requested_genes = [gene for gene in spec.gene_labels if gene in frequency]
    genes = requested_genes + [gene for gene in ranked_genes if gene not in requested_genes]
    genes = genes[:gene_cap]
    selected.attrs['target_gene_truncated'] = max(0, len(requested_genes) - len(genes))
    selected['_genes_kept'] = selected['_genes'].map(lambda values: [gene for gene in values if gene in genes])
    return selected, genes


class NatureEnrichmentChord:
    """Readable chord-style membership view with explicit truncation caps."""

    plot_type = 'enrichment_chord'

    def render(self, data, spec: FigureSpec, container=None):
        from matplotlib.patches import FancyArrowPatch

        style = get_style(spec.style)
        term_cap = 4 if spec.width == 'single' else 6
        gene_cap = min(spec.max_genes, 14 if spec.width == 'single' else 30)
        frame = _enrichment_frame(data, spec, require_genes=True)
        selected, genes = _membership_subset(frame, spec, term_cap, gene_cap)
        if not genes:
            raise ValueError('Chord 图没有可用的基因成员。')
        with style.context(spec):
            fig, ax = _new_figure(spec.with_updates(plot_type=self.plot_type), style, container)
            term_y = np.linspace(.82, -.82, len(selected))
            gene_y = np.linspace(.92, -.92, len(genes))
            term_pos = {}
            gene_pos = {}
            palette = _network_colors()
            for index, y in enumerate(term_y):
                term_pos[index] = (-.72, y)
                ax.scatter([-.72], [y], s=58, marker='s', color=palette[index % len(palette)], edgecolor='white', linewidth=.6, zorder=4)
                ax.text(-.84, y, wrap_term(selected.loc[index, '_display_term'], width=22 if spec.width == 'single' else 31, max_chars=spec.term_max_chars), ha='right', va='center', color=style.text, fontsize=style.profile(spec).tick_font_pt)
            for index, gene in enumerate(genes):
                gene_pos[gene] = (.72, gene_y[index])
                ax.scatter([.72], [gene_y[index]], s=18, color='#B8C1CB', edgecolor='white', linewidth=.45, zorder=4)
                # Gene nodes and labels are separate controls: all membership
                # nodes remain connected, while the requested label strategy
                # determines which symbols are printed.
            label_genes, omitted_gene_labels = gene_label_order(
                genes, spec.gene_labels, strategy=spec.gene_label_strategy,
                membership=selected['_genes_kept'], maximum=spec.max_gene_labels,
            )
            for gene in label_genes:
                ax.text(.84, gene_pos[gene][1], gene, ha='left', va='center', color=style.muted_text, fontsize=max(5.0, style.profile(spec).tick_font_pt - .6))
            for term_index, row in selected.iterrows():
                for gene in row['_genes_kept']:
                    patch = FancyArrowPatch(term_pos[term_index], gene_pos[gene], arrowstyle='-',
                                            connectionstyle='arc3,rad=0.18', color='#B7C0C9',
                                            linewidth=.55, alpha=.78, zorder=1)
                    ax.add_patch(patch)
            ax.text(-.72, 1.07, 'Pathways', ha='center', va='bottom', color=style.muted_text, fontsize=style.profile(spec).legend_font_pt)
            ax.text(.72, 1.07, 'Genes', ha='center', va='bottom', color=style.muted_text, fontsize=style.profile(spec).legend_font_pt)
            if len(frame) > len(selected) or sum(len(item) for item in frame['_genes']) > len(genes):
                ax.text(0.0, -1.055,
                        f'Displayed subset: {len(selected)} pathways · {len(genes)} genes',
                        ha='center', va='top', color=style.muted_text,
                        fontsize=max(4.8, style.profile(spec).legend_font_pt - .6))
            ax.set_title(spec.title or 'Pathway–gene membership', loc='left', pad=5)
            ax.set_xlim(-1.92, 1.62)
            ax.set_ylim(-1.10, 1.18)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.axis('off')
            fig.subplots_adjust(left=.02, right=.98, top=.88, bottom=.04)
        warnings = []
        if len(frame) > len(selected) or sum(len(item) for item in frame['_genes']) > len(genes):
            warnings.append(f'Chord 图为可读性截断视图：展示 {len(selected)} 个通路和 {len(genes)} 个基因。')
        if omitted_gene_labels:
            warnings.append(f'仅标注 {len(label_genes)}/{len(genes)} 个基因；可提高“最大基因标签数”或指定目标基因。')
        if selected.attrs.get('target_selected'):
            warnings.append(f"已按用户选择展示 {selected.attrs['target_selected']} 个目标通路。")
        if selected.attrs.get('target_unmatched'):
            warnings.append('未匹配目标通路：' + '、'.join(selected.attrs['target_unmatched']))
        if selected.attrs.get('target_truncated'):
            warnings.append(f"目标通路超过单栏网络图上限，隐藏 {selected.attrs['target_truncated']} 个；可切换双栏或减少目标通路。")
        if selected.attrs.get('target_gene_truncated'):
            warnings.append(f"指定基因超过网络图节点上限，隐藏 {selected.attrs['target_gene_truncated']} 个；可切换双栏或减少目标基因。")
        if selected.attrs.get('redundant_removed', 0):
            warnings.append(f"按命中基因 Jaccard ≥ {spec.redundancy_threshold:.2f} 压缩了 {selected.attrs['redundant_removed']} 个冗余通路。")
        return attach_contract(fig, spec.with_updates(plot_type=self.plot_type), style,
                               semantic_warnings=warnings,
                               encodings={'edge': 'gene membership', 'color': 'directly labelled pathway', 'shape': 'pathway square / gene circle'})


class NatureEnrichmentCnetplot:
    """Compact cnetplot-like network with a numbered pathway key."""

    plot_type = 'enrichment_cnetplot'

    def render(self, data, spec: FigureSpec, container=None):
        style = get_style(spec.style)
        term_cap = 4 if spec.width == 'single' else 6
        gene_cap = min(spec.max_genes, 14 if spec.width == 'single' else 28)
        frame = _enrichment_frame(data, spec, require_genes=True)
        selected, genes = _membership_subset(frame, spec, term_cap, gene_cap)
        if not genes:
            raise ValueError('Cnetplot 没有可用的基因成员。')
        with style.context(spec):
            fig, ax = _new_figure(spec.with_updates(plot_type=self.plot_type), style, container)
            palette = _network_colors()
            # Put the network to the right and reserve a stable key on the left.
            term_angles = np.linspace(np.pi * .15, np.pi * 1.85, len(selected), endpoint=False)
            gene_angles = np.linspace(0, 2 * np.pi, len(genes), endpoint=False)
            network_center = (.32, .02)
            term_pos = {index: (network_center[0] + .22 * np.cos(angle), network_center[1] + .22 * np.sin(angle)) for index, angle in enumerate(term_angles)}
            gene_pos = {gene: (network_center[0] + .76 * np.cos(angle), network_center[1] + .76 * np.sin(angle)) for gene, angle in zip(genes, gene_angles)}
            for term_index, row in selected.iterrows():
                for gene in row['_genes_kept']:
                    ax.plot([term_pos[term_index][0], gene_pos[gene][0]], [term_pos[term_index][1], gene_pos[gene][1]], color=palette[term_index % len(palette)], linewidth=.55, alpha=.28, zorder=1)
            for gene, (x, y) in gene_pos.items():
                ax.scatter([x], [y], s=22, color='#C4CCD5', edgecolor='white', linewidth=.5, zorder=3)
            for term_index, (x, y) in term_pos.items():
                ax.scatter([x], [y], s=100, marker='s', color=palette[term_index % len(palette)], edgecolor='white', linewidth=.8, zorder=4)
                ax.text(x, y, f'P{term_index + 1}', ha='center', va='center', fontsize=max(5.2, style.profile(spec).tick_font_pt - .4), color='white', fontweight='bold', zorder=5)
            label_genes, omitted_gene_labels = gene_label_order(
                genes, spec.gene_labels, strategy=spec.gene_label_strategy,
                membership=selected['_genes_kept'], maximum=spec.max_gene_labels,
            )
            for gene in label_genes:
                x, y = gene_pos[gene]
                angle = math.atan2(y - network_center[1], x - network_center[0])
                ha = 'left' if np.cos(angle) >= 0 else 'right'
                ax.text(x + .06 * np.cos(angle), y + .06 * np.sin(angle), gene, ha=ha, va='center', color=style.muted_text, fontsize=max(5.0, style.profile(spec).tick_font_pt - .7), zorder=5)
            key_y = .82
            ax.text(-1.55, 1.03, 'Pathway key', ha='left', va='bottom', color=style.muted_text, fontsize=style.profile(spec).legend_font_pt)
            for term_index, row in selected.iterrows():
                label = wrap_term(row['_display_term'], width=21 if spec.width == 'single' else 34, max_chars=spec.term_max_chars)
                label = f'P{term_index + 1}  {label}'
                ax.text(-1.55, key_y, label, ha='left', va='top', color=style.text, fontsize=max(5.0, style.profile(spec).tick_font_pt - .4))
                key_y -= .17 + .075 * label.count('\n')
            if len(frame) > len(selected) or sum(len(item) for item in frame['_genes']) > len(genes):
                ax.text(-.15, -1.055,
                        f'Displayed subset: {len(selected)} pathways · {len(genes)} genes',
                        ha='center', va='top', color=style.muted_text,
                        fontsize=max(4.8, style.profile(spec).legend_font_pt - .6))
            ax.set_title(spec.title or 'Pathway–gene cnetplot', loc='left', pad=5)
            ax.set_xlim(-1.64, 1.33)
            ax.set_ylim(-1.12, 1.14)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.axis('off')
            fig.subplots_adjust(left=.02, right=.98, top=.88, bottom=.03)
        warnings = []
        if len(frame) > len(selected) or sum(len(item) for item in frame['_genes']) > len(genes):
            warnings.append(f'Cnetplot 为可读性截断视图：展示 {len(selected)} 个通路和 {len(genes)} 个基因。')
        if omitted_gene_labels:
            warnings.append(f'仅标注 {len(label_genes)}/{len(genes)} 个基因；可提高“最大基因标签数”或指定目标基因。')
        if selected.attrs.get('target_selected'):
            warnings.append(f"已按用户选择展示 {selected.attrs['target_selected']} 个目标通路。")
        if selected.attrs.get('target_unmatched'):
            warnings.append('未匹配目标通路：' + '、'.join(selected.attrs['target_unmatched']))
        if selected.attrs.get('target_truncated'):
            warnings.append(f"目标通路超过单栏网络图上限，隐藏 {selected.attrs['target_truncated']} 个；可切换双栏或减少目标通路。")
        if selected.attrs.get('target_gene_truncated'):
            warnings.append(f"指定基因超过网络图节点上限，隐藏 {selected.attrs['target_gene_truncated']} 个；可切换双栏或减少目标基因。")
        if selected.attrs.get('redundant_removed', 0):
            warnings.append(f"按命中基因 Jaccard ≥ {spec.redundancy_threshold:.2f} 压缩了 {selected.attrs['redundant_removed']} 个冗余通路。")
        return attach_contract(fig, spec.with_updates(plot_type=self.plot_type), style,
                               semantic_warnings=warnings,
                               encodings={'edge': 'pathway–gene membership', 'color': 'directly labelled pathway', 'shape': 'pathway square / gene circle'})


class NatureEnrichmentEmapplot:
    """Pathway similarity map using Jaccard edges and a compact term key."""

    plot_type = 'enrichment_emapplot'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize

        style = get_style(spec.style)
        term_cap = 8 if spec.width == 'single' else 10
        candidates = _enrichment_frame(data, spec, require_genes=True).sort_values(['_fdr', '_count'])
        selected = _select_terms(candidates, spec, direction=False)
        frame = selected.head(term_cap).copy()
        selection_attrs = dict(selected.attrs)
        if len(selected) > term_cap:
            selection_attrs['target_truncated'] = len(selected) - term_cap
        redundant_removed = int(selected.attrs.get('redundant_removed', 0))
        frame = frame.reset_index(drop=True)
        sets = [set(values) for values in frame['_genes']]
        n = len(frame)
        if n < 1:
            raise ValueError('Emapplot 没有含基因成员的目标通路。')
        similarity = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i + 1, n):
                union = sets[i] | sets[j]
                similarity[i, j] = similarity[j, i] = len(sets[i] & sets[j]) / len(union) if union else 0
        with style.context(spec):
            fig, ax = _new_figure(spec.with_updates(plot_type=self.plot_type), style, container)
            center = np.array([-.35, .02])
            radius = .68
            angles = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, n, endpoint=False)
            pos = np.column_stack([center[0] + radius * np.cos(angles), center[1] + radius * np.sin(angles)])
            for i in range(n):
                for j in range(i + 1, n):
                    if similarity[i, j] >= spec.similarity_threshold:
                        ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], color=style.neutral_dark, linewidth=.45 + 2.4 * similarity[i, j], alpha=.48, zorder=1)
            counts = frame['_count'].to_numpy(float)
            sqrt_counts = np.sqrt(counts)
            sizes = 90 + 320 * (sqrt_counts - sqrt_counts.min()) / max(np.ptp(sqrt_counts), 1e-9)
            log_fdr = -np.log10(frame['_fdr'].to_numpy(float))
            norm = Normalize(vmin=float(log_fdr.min()), vmax=float(log_fdr.max()) if log_fdr.max() > log_fdr.min() else float(log_fdr.min() + 1))
            nodes = ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=log_fdr, cmap=style.fdr_cmap(), norm=norm, edgecolor='white', linewidth=.7, zorder=4)
            for i, (x, y) in enumerate(pos):
                ax.text(x, y, f'E{i + 1}', ha='center', va='center', color=style.text, fontsize=max(5.2, style.profile(spec).tick_font_pt - .5), fontweight='semibold', zorder=5)
            ax.text(.63, 1.03, 'Term key', ha='left', va='bottom', color=style.muted_text, fontsize=style.profile(spec).legend_font_pt)
            key_y = .84
            for i, term in enumerate(frame['_display_term']):
                text = wrap_term(term, width=20 if spec.width == 'single' else 30, max_chars=spec.term_max_chars)
                text = f'E{i + 1}  {text}'
                ax.text(.63, key_y, text, ha='left', va='top', color=style.text, fontsize=max(5.0, style.profile(spec).tick_font_pt - .5))
                key_y -= .16 + .07 * text.count('\n')
            fig.subplots_adjust(left=.02, right=.98, top=.88, bottom=.06)
            cbar_ax = fig.add_axes([.57, .075, .28, .025])
            cbar = fig.colorbar(nodes, cax=cbar_ax, orientation='horizontal')
            cbar.outline.set_visible(False)
            cbar.set_label(r'$-\log_{10}$(FDR)', fontsize=style.profile(spec).legend_font_pt, labelpad=2)
            cbar.ax.tick_params(length=1.8, labelsize=max(5.0, style.profile(spec).tick_font_pt - .5), width=.45)
            ax.set_title(spec.title or 'Pathway similarity map', loc='left', pad=5)
            ax.text(-.35, -1.03, f'Node size = Count · edge = Jaccard ≥ {spec.similarity_threshold:.2f}', ha='center', va='top', color=style.muted_text, fontsize=max(5.0, style.profile(spec).legend_font_pt - .4))
            ax.set_xlim(-1.26, 1.72)
            ax.set_ylim(-1.16, 1.14)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.axis('off')
        warnings = []
        edge_count = int(np.sum(np.triu(similarity >= spec.similarity_threshold, 1)))
        if edge_count == 0:
            warnings.append('当前 Jaccard 阈值下没有通路相似性边；请将其作为独立通路解读。')
        if len(frame) >= term_cap:
            warnings.append(f'仅展示相似性图前 {term_cap} 个显著通路。')
        if selection_attrs.get('target_selected'):
            warnings.append(f"已按用户选择展示 {selection_attrs['target_selected']} 个目标通路。")
        if selection_attrs.get('target_unmatched'):
            warnings.append('未匹配目标通路：' + '、'.join(selection_attrs['target_unmatched']))
        if selection_attrs.get('target_truncated'):
            warnings.append(f"目标通路超过图型上限，隐藏 {selection_attrs['target_truncated']} 个；可减少目标通路或切换双栏。")
        if redundant_removed:
            warnings.append(f'按命中基因 Jaccard ≥ {spec.redundancy_threshold:.2f} 压缩了 {redundant_removed} 个冗余通路。')
        return attach_contract(fig, spec.with_updates(plot_type=self.plot_type), style,
                               semantic_warnings=warnings,
                               encodings={'edge': 'Jaccard shared genes', 'color': '-log10(FDR)', 'size': 'Count'})


def _curve_list(data):
    if isinstance(data, dict) and 'curves' in data:
        curves = list(data['curves'])
        common = {key: value for key, value in data.items() if key != 'curves'}
        return [{**common, **dict(curve)} for curve in curves]
    if isinstance(data, (list, tuple)):
        return [dict(item) for item in data]
    if isinstance(data, dict):
        return [dict(data)]
    raise ValueError('GSEA running 模板需要 mapping 或 curves 列表。')


def _running_curve(curve):
    ranking = curve.get('ranking', curve.get('rank_metric'))
    ranking = np.asarray(ranking, dtype=float).reshape(-1) if ranking is not None else None
    hits = curve.get('hit_indices', curve.get('hits'))
    if isinstance(hits, str):
        hits = _genes(hits)
    hits = np.asarray(hits if hits is not None else [], dtype=int).reshape(-1)
    res = curve.get('running_score', curve.get('RES', curve.get('res')))
    if isinstance(res, str):
        try:
            res = ast.literal_eval(res)
        except (SyntaxError, ValueError):
            res = None
    res = np.asarray(res, dtype=float).reshape(-1) if res is not None else None
    if res is not None and (not np.isfinite(res).all() or np.nanmax(np.abs(res)) > 1.5):
        # The preview script used an unnormalised cumulative count.  Recompute
        # from ranking/hits when possible rather than presenting a misleading ES.
        res = None
    if res is None:
        if ranking is None or ranking.size < 2 or hits.size == 0:
            raise ValueError('GSEA running 图需要 ranking + hit_indices，或标准化 RES。')
        hits = np.unique(hits[(hits >= 0) & (hits < ranking.size)])
        if hits.size == 0 or hits.size >= ranking.size:
            raise ValueError('GSEA hit_indices 不在 ranking 范围内。')
        weights = np.abs(ranking[hits])
        if weights.sum() <= 0:
            weights = np.ones_like(weights)
        increments = np.full(ranking.size, -1.0 / (ranking.size - hits.size))
        increments[hits] = weights / weights.sum()
        res = np.cumsum(increments)
    if ranking is None:
        ranking = np.zeros(res.size, dtype=float)
    if res.size != ranking.size:
        # OmicVerse/gseapy should always provide a full-length RES.  A shorter
        # curve has no unambiguous rank coordinate, so reject it explicitly.
        raise ValueError('GSEA running RES 长度必须与 ranking 一致。')
    return ranking, hits, res


class NatureGSEARunning:
    """Classic GSEA running ES + hit rug, using real ranking details."""

    plot_type = 'gsea_running'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt

        style = get_style(spec.style)
        all_curves = _curve_list(data)
        target_mode = bool(spec.pathway_terms) and spec.pathway_selection in {'selected', 'selected_plus_top'}
        unmatched_targets = []
        selected_target_count = 0
        if target_mode:
            curve_frame = pd.DataFrame({
                '_term': [str(curve.get('term', curve.get('Term', curve.get('Description', '')))) for curve in all_curves],
            })
            curve_frame['_display_term'] = curve_frame['_term'].map(strip_term_id)
            requested, unmatched_targets, indices = requested_term_rows(
                curve_frame, spec.pathway_terms,
            )
            selected_target_count = len(indices)
            selected_indices = list(requested.index)
            if spec.pathway_selection == 'selected':
                all_curves = [all_curves[index] for index in selected_indices]
            else:
                remainder = [curve for index, curve in enumerate(all_curves) if index not in selected_indices]
                all_curves = [all_curves[index] for index in selected_indices] + remainder
        curves = all_curves[:spec.running_term_n]
        if not curves:
            raise ValueError('GSEA running 图没有可绘制的通路。')
        with style.context(spec):
            if container is None:
                profile = style.profile(spec.with_updates(plot_type=self.plot_type))
                fig = plt.figure(figsize=profile.figsize, dpi=profile.dpi)
            else:
                fig = container
            gs = fig.add_gridspec(len(curves) * 2, 1, height_ratios=sum(([3.2, .58] for _ in curves), []), hspace=.10)
            for index, curve in enumerate(curves):
                ranking, hits, res = _running_curve(curve)
                direction_score = _number(curve.get('NES', curve.get('nes')), default=float(res[int(np.argmax(np.abs(res)))]))
                color = style.signal_red if direction_score >= 0 else style.signal_blue
                ax = fig.add_subplot(gs[index * 2, 0])
                rug = fig.add_subplot(gs[index * 2 + 1, 0], sharex=ax)
                x = np.arange(res.size)
                peak = int(np.argmax(np.abs(res)))
                ax.plot(x, res, color=color, linewidth=1.1, zorder=3)
                ax.axhline(0, color=style.neutral, linewidth=.55, zorder=0)
                ax.scatter([peak], [res[peak]], s=22, color=color, edgecolor='white', linewidth=.55, zorder=4)
                ax.text(.015, .86, wrap_term(curve.get('term', curve.get('Term', 'GSEA pathway')), width=30, max_chars=spec.term_max_chars), transform=ax.transAxes, ha='left', va='top', color=style.text, fontsize=style.profile(spec).legend_font_pt, fontweight='semibold')
                nes = _number(curve.get('NES', curve.get('nes')), default=np.nan)
                fdr = _number(curve.get('FDR', curve.get('fdr', curve.get('FDR q-val'))), default=np.nan)
                metrics = f'NES = {nes:+.2f}' if np.isfinite(nes) else 'NES = n/a'
                metrics += f'   FDR = {fdr:.3g}' if np.isfinite(fdr) else '   FDR = n/a'
                ax.text(.985, .86, metrics, transform=ax.transAxes, ha='right', va='top', color=style.text, fontsize=style.profile(spec).legend_font_pt, bbox={'boxstyle': 'round,pad=.22', 'facecolor': 'white', 'edgecolor': style.subtle, 'linewidth': .45})
                style.apply_axis(ax, style.profile(spec.with_updates(plot_type=self.plot_type)))
                ax.spines['left'].set_visible(True)
                ax.set_ylabel('ES')
                ax.set_xlim(0, max(1, res.size - 1))
                pad = max(.05, .12 * max(abs(float(res.min())), abs(float(res.max()))))
                ax.set_ylim(float(res.min()) - pad, float(res.max()) + pad)
                # Keep one opaque direction colour across the curve and rug;
                # alpha variants would be interpreted as duplicate palette
                # entries by the colour-vision validator.
                rug.vlines(hits, 0, 1, color=color, linewidth=.55, alpha=1.0)
                rug.set_ylim(0, 1)
                rug.set_yticks([])
                rug.spines['top'].set_visible(False)
                rug.spines['right'].set_visible(False)
                rug.spines['left'].set_visible(False)
                rug.spines['bottom'].set_visible(index == len(curves) - 1)
                rug.tick_params(axis='x', labelsize=max(5.2, style.profile(spec).tick_font_pt - .4), length=1.8)
                if index < len(curves) - 1:
                    plt.setp(ax.get_xticklabels(), visible=False)
                    plt.setp(rug.get_xticklabels(), visible=False)
            rug.set_xlabel('Rank in ordered dataset')
            if container is None:
                fig.suptitle(spec.title or 'GSEA running enrichment score', x=.02, ha='left', y=.985, fontsize=style.profile(spec).title_font_pt, fontweight='semibold', color=style.text)
                fig.subplots_adjust(left=.14 if spec.width == 'single' else .10, right=.92, top=.88, bottom=.14)
        warnings = []
        if target_mode and selected_target_count:
            warnings.append(f'已按用户选择展示 {selected_target_count} 个目标通路。')
        if unmatched_targets:
            warnings.append('未匹配目标通路：' + '、'.join(unmatched_targets))
        if target_mode and len(all_curves) > len(curves):
            warnings.append(f'运行曲线最多展示 {spec.running_term_n} 条；完整曲线保留在结果表。')
        return attach_contract(fig, spec.with_updates(plot_type=self.plot_type), style,
                               semantic_warnings=warnings,
                               encodings={'x': 'ranked gene position', 'y': 'running enrichment score (ES)', 'rug': 'gene-set hits', 'color': 'NES direction'})
