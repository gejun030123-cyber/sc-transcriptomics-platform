"""NatureMA：平均表达量与 log2 fold-change 的固定投稿模板。"""

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import adjust_labels, attach_contract, finite_numeric, prune_overlapping_labels


def _first_column(frame, candidates):
    lookup = {str(column).strip().lower(): column for column in frame.columns}
    return next((lookup[name.lower()] for name in candidates if name.lower() in lookup), None)


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
        y_limit = max(spec.fc_threshold * 1.6, float(np.nanpercentile(np.abs(y), 99.5)) * 1.08, 1.5)
        y_limit = min(12.0, float(np.ceil(y_limit * 2) / 2))
        y_display = np.clip(y, -y_limit * 0.985, y_limit * 0.985)

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
            ax.set_title(spec.title or 'MA plot', loc='left', pad=5)
            style.apply_axis(ax, profile)

            texts = []
            candidates = list(frame.index[frame['_gene'].isin(spec.label_genes)])
            if spec.label_strategy != 'none' and len(candidates) < spec.label_n:
                ranked = frame.assign(_abs_fc=np.abs(y)).sort_values(
                    ['_fdr', '_abs_fc'], ascending=[True, False])
                ranked = ranked[ranked.index.isin(frame.index[regulation != 'NS'])]
                for index in ranked.index:
                    if index not in candidates:
                        candidates.append(index)
                    if len(candidates) >= spec.label_n:
                        break
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
            fig.subplots_adjust(left=0.19, right=0.97, bottom=0.17, top=0.90)

        if container is not None:
            return container
        return attach_contract(
            fig, spec, style,
            semantic_warnings=semantic_warnings,
            encodings={'x': 'mean expression', 'y': 'log2 fold change', 'color': 'DEG direction'},
            validation_texts=texts,
        )
