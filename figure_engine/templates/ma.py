"""NatureMA：平均表达量与 log2 fold-change 的固定投稿模板。"""

import textwrap

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import adjust_labels, attach_contract, finite_numeric, prune_overlapping_labels


def _first_column(frame, candidates):
    lookup = {str(column).strip().lower(): column for column in frame.columns}
    return next((lookup[name.lower()] for name in candidates if name.lower() in lookup), None)


def _wrap_title(value, default, width):
    """Wrap each supplied title line independently to retain context grouping."""
    text = str(value or default).replace('_', ' ')
    lines = [line.strip() for line in text.splitlines() if line.strip()] or [default]
    return '\n'.join(
        textwrap.fill(line, width=width, break_long_words=False)
        for line in lines
    )


def _select_labels(frame, spec, regulation):
    """Prioritise requested genes, then fill the remaining label budget fairly."""
    selected = []
    genes = frame['_gene'].astype(str)
    # Preserve the scientist's requested ordering rather than the incidental
    # order of the DEG table. A requested gene may be non-significant; it is
    # still a legitimate point to inspect on an MA plot.
    for gene in spec.label_genes:
        for index in frame.index[genes.eq(str(gene))]:
            if index not in selected:
                selected.append(index)
    remaining = max(0, int(spec.label_n) - len(selected))
    if remaining and spec.label_strategy != 'none':
        significant = frame.loc[regulation != 'NS'].copy()
        significant['_abs_fc'] = significant['_log2fc'].abs()
        significant = significant.sort_values(['_fdr', '_abs_fc'], ascending=[True, False])
        per_side = max(1, int(np.ceil(remaining / 2)))
        candidates = pd.concat([
            significant.loc[regulation.loc[significant.index] == 'Up'].head(per_side),
            significant.loc[regulation.loc[significant.index] == 'Down'].head(per_side),
        ]).sort_values(['_fdr', '_abs_fc'], ascending=[True, False])
        for index in candidates.index:
            if index not in selected:
                selected.append(index)
            if len(selected) >= int(spec.label_n):
                break
    return selected


class NatureMA:
    plot_type = 'ma'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        style = get_style(spec.style)
        profile = style.profile(spec)
        frame = pd.DataFrame(data).copy()
        mean_column = _first_column(frame, (
            'baseMean', 'mean_expression', 'meanExpression', 'AveExpr',
            'average_expression', 'mean',
        ))
        fc_column = _first_column(frame, ('log2FC', 'log2FoldChange', 'logFC'))
        fdr_column = _first_column(frame, ('padj', 'FDR', 'fdr', 'qvalue'))
        gene_column = _first_column(frame, ('gene', 'gene_name', 'symbol', 'Gene'))
        if mean_column is None or fc_column is None:
            raise ValueError('NatureMA 需要平均表达量和 log2FC 字段')

        frame['_mean'] = finite_numeric(frame[mean_column])
        frame['_log2fc'] = finite_numeric(frame[fc_column])
        frame['_fdr'] = (np.clip(finite_numeric(frame[fdr_column], fill=1.0),
                                 np.finfo(float).tiny, 1.0)
                         if fdr_column else np.ones(len(frame)))
        frame['_gene'] = (frame[gene_column].astype(str) if gene_column
                          else frame.index.astype(str))
        frame = frame[np.isfinite(frame['_mean']) & np.isfinite(frame['_log2fc'])].copy()
        if frame.empty:
            raise ValueError('NatureMA 没有可绘制的有限数值')

        x = np.log2(np.clip(frame['_mean'].to_numpy(dtype=float), 0, None) + 1.0)
        y = frame['_log2fc'].to_numpy(dtype=float)
        significant = frame['_fdr'].to_numpy(dtype=float) < spec.fdr_threshold
        regulation = np.full(len(frame), 'NS', dtype=object)
        regulation[significant & (y >= spec.fc_threshold)] = 'Up'
        regulation[significant & (y <= -spec.fc_threshold)] = 'Down'
        # Do not squeeze legitimate extreme fold changes into artificial rows
        # against the upper/lower MA boundary.  Unlike Volcano FDR tails, the
        # y-axis is itself the effect size, so clipping/compressing it would
        # alter the visual evidence.  The full observed range with 10% headroom
        # keeps points and their selected labels inside the axes.
        max_abs_fc = float(np.nanmax(np.abs(y)))
        y_limit = max(spec.fc_threshold * 1.6, max_abs_fc * 1.10, 1.5)
        y_limit = float(np.ceil(y_limit * 2) / 2)
        y_display = y

        with style.context(spec):
            if container is None:
                fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
            else:
                fig = container
                ax = container.subplots()
            for label, size, alpha, zorder in (
                ('NS', 5.5, 0.38, 1), ('Down', 8.0, 0.82, 2), ('Up', 8.0, 0.82, 2),
            ):
                mask = regulation == label
                if mask.any():
                    ax.scatter(x[mask], y_display[mask], s=size, color=style.deg_palette[label],
                               alpha=alpha, linewidths=0, rasterized=True, zorder=zorder)
            ax.axhline(0, color=style.neutral_dark, linewidth=0.55, zorder=0)
            ax.axhline(spec.fc_threshold, color=style.neutral_dark,
                       linestyle=(0, (3, 2)), linewidth=0.5, zorder=0)
            ax.axhline(-spec.fc_threshold, color=style.neutral_dark,
                       linestyle=(0, (3, 2)), linewidth=0.5, zorder=0)
            ax.set_ylim(-y_limit, y_limit)
            ax.set_xlabel(r'log$_2$(mean expression + 1)')
            ax.set_ylabel(r'log$_2$(fold change)')
            title_width = 38 if spec.width == 'single' else 72
            wrapped_title = _wrap_title(spec.title, 'MA plot', title_width)
            ax.set_title(
                wrapped_title,
                loc='left', pad=5,
            )
            style.apply_axis(ax, profile)

            texts = []
            candidates = _select_labels(frame, spec, pd.Series(regulation, index=frame.index))
            max_labels = 8 if spec.width == 'single' else 12
            for index in candidates[:max_labels]:
                position = frame.index.get_loc(index)
                texts.append(ax.text(
                    x[position], y_display[position], frame.at[index, '_gene'],
                    fontsize=profile.tick_font_pt, color=style.text,
                    ha='left' if y_display[position] >= 0 else 'right', va='bottom',
                    bbox={'boxstyle': 'round,pad=0.12', 'facecolor': 'white',
                          'edgecolor': style.subtle, 'linewidth': 0.35}, zorder=5,
                ))
            adjust_labels(texts, ax, arrow_color=style.neutral_dark)
            semantic_warnings = []
            hidden_labels = prune_overlapping_labels(texts, ax)
            if hidden_labels:
                semantic_warnings.append(
                    f'为避免最终尺寸文字重叠，已隐藏 {len(hidden_labels)} 个低优先级基因标签；完整基因列表保留在结果表。'
                )
            if spec.show_legend:
                handles = [
                    Line2D([], [], linestyle='None', marker='o', markersize=3.8,
                           markerfacecolor=style.deg_palette[label], markeredgewidth=0,
                           label=label)
                    for label in ('Up', 'Down', 'NS')
                ]
                ax.legend(handles=handles, loc='upper right', frameon=False)
            title_lines = wrapped_title.count('\n') + 1
            top = max(0.76, 0.90 - 0.055 * (title_lines - 1))
            fig.subplots_adjust(left=0.19, right=0.97, bottom=0.17, top=top)

        if container is not None:
            return container
        return attach_contract(
            fig, spec, style,
            semantic_warnings=semantic_warnings,
            encodings={'x': 'mean expression', 'y': 'log2 fold change', 'color': 'DEG direction'},
            validation_texts=texts,
        )
