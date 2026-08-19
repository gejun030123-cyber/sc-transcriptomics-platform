"""NatureVolcano：固定阈值语义和受控基因标注。"""

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import adjust_labels, attach_contract, finite_numeric, prune_overlapping_labels


def _select_labels(frame, spec, regulation):
    genes = frame['_gene'].astype(str)
    selected = []
    custom = set(spec.label_genes)
    if custom:
        selected.extend(frame.index[genes.isin(custom)].tolist())
    remaining = max(0, int(spec.label_n) - len(selected))
    if remaining and spec.label_strategy != 'none':
        significant = frame[regulation.isin(['Up', 'Down'])].copy()
        significant['_abs_fc'] = significant['_log2fc'].abs()
        significant = significant.sort_values(['_fdr', '_abs_fc'], ascending=[True, False])
        per_side = max(1, int(np.ceil(remaining / 2)))
        candidates = pd.concat([
            significant[regulation.loc[significant.index] == 'Up'].head(per_side),
            significant[regulation.loc[significant.index] == 'Down'].head(per_side),
        ]).sort_values(['_fdr', '_abs_fc'], ascending=[True, False])
        for index in candidates.index:
            if index not in selected:
                selected.append(index)
            if len(selected) >= int(spec.label_n):
                break
    return selected


class NatureVolcano:
    """固定的 Bulk RNA-seq Volcano publication 模板。"""

    plot_type = 'volcano'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        style = get_style(spec.style)
        profile = style.profile(spec)
        frame = pd.DataFrame(data).copy()
        required = {'gene', 'log2FC'}
        if not required.issubset(frame.columns):
            raise ValueError(f'NatureVolcano 缺少字段: {sorted(required - set(frame.columns))}')
        fdr_column = next((name for name in ('padj', 'fdr', 'FDR', 'qvalue') if name in frame.columns), None)
        if fdr_column is None:
            raise ValueError('NatureVolcano 需要 padj/fdr/qvalue 字段')
        frame['_gene'] = frame['gene'].astype(str)
        frame['_log2fc'] = finite_numeric(frame['log2FC'])
        frame['_fdr'] = np.clip(finite_numeric(frame[fdr_column], fill=1.0),
                                np.finfo(float).tiny, 1.0)
        frame = frame[np.isfinite(frame['_log2fc'])].copy()
        if frame.empty:
            raise ValueError('NatureVolcano 没有有限的 log2FC 数据')
        regulation = pd.Series('NS', index=frame.index, dtype=object)
        regulation[(frame['_fdr'] < spec.fdr_threshold)
                   & (frame['_log2fc'] >= spec.fc_threshold)] = 'Up'
        regulation[(frame['_fdr'] < spec.fdr_threshold)
                   & (frame['_log2fc'] <= -spec.fc_threshold)] = 'Down'
        y_raw = -np.log10(frame['_fdr'].to_numpy(dtype=float))
        robust_y = float(np.nanpercentile(y_raw, 99.5)) if len(y_raw) else 0.0
        threshold_y = -np.log10(spec.fdr_threshold)
        y_limit = max(4.0, threshold_y + 1.5, np.ceil(robust_y + 0.5))
        y_limit = min(40.0, float(y_limit))
        clipped = y_raw > y_limit
        y_display = np.minimum(y_raw, y_limit * 0.985)
        x_abs = np.abs(frame['_log2fc'].to_numpy(dtype=float))
        x_limit = max(spec.fc_threshold * 1.6, float(np.nanpercentile(x_abs, 99.8)) * 1.08, 1.5)
        x_limit = float(np.ceil(x_limit * 2) / 2)
        x_display = np.clip(frame['_log2fc'].to_numpy(dtype=float), -x_limit * 0.99, x_limit * 0.99)

        warnings = []
        if len(frame) > 100000:
            warnings.append('Volcano 点数超过 100,000；SVG 中的数据层将栅格化以控制文件体积。')
        max_labels = 8 if spec.width == 'single' else 12
        if spec.label_n > max_labels:
            warnings.append(f'{spec.width} 画布最多建议标注 {max_labels} 个基因；已自动限制。')
            label_spec = spec.with_updates(label_n=max_labels)
        else:
            label_spec = spec

        with style.context(spec):
            if container is None:
                fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
            else:
                fig = container
                ax = container.subplots()
            for label, point_size, alpha, zorder in (
                ('NS', 5.5 if spec.width == 'single' else 7.0, 0.38, 1),
                ('Down', 8.0 if spec.width == 'single' else 10.0, 0.82, 2),
                ('Up', 8.0 if spec.width == 'single' else 10.0, 0.82, 2),
            ):
                mask = regulation.to_numpy() == label
                if mask.any():
                    ax.scatter(
                        x_display[mask], y_display[mask],
                        s=point_size, color=style.deg_palette[label],
                        alpha=alpha, linewidths=0, rasterized=True, zorder=zorder,
                    )
            if clipped.any():
                ax.scatter(
                    x_display[clipped], y_display[clipped],
                    s=11, marker='^',
                    c=[style.deg_palette[value] for value in regulation.to_numpy()[clipped]],
                    linewidths=0, alpha=0.9, rasterized=True, zorder=3,
                )

            threshold_color = style.neutral_dark
            ax.axhline(threshold_y, color=threshold_color, linestyle=(0, (3, 2)),
                       linewidth=0.55, zorder=0)
            ax.axvline(spec.fc_threshold, color=threshold_color, linestyle=(0, (3, 2)),
                       linewidth=0.55, zorder=0)
            ax.axvline(-spec.fc_threshold, color=threshold_color, linestyle=(0, (3, 2)),
                       linewidth=0.55, zorder=0)
            ax.set_xlim(-x_limit, x_limit)
            ax.set_ylim(0, y_limit)
            ax.set_xlabel(r'log$_2$(fold change)')
            ax.set_ylabel(r'$-\log_{10}$(FDR)')
            ax.set_title(spec.title or 'Differential expression', loc='left', pad=5)
            style.apply_axis(ax, profile)

            counts = {label: int((regulation == label).sum()) for label in ('Up', 'Down', 'NS')}
            selected = _select_labels(frame, label_spec, regulation)
            texts = []
            for index in selected:
                row_position = frame.index.get_loc(index)
                x_value = x_display[row_position]
                y_value = min(y_display[row_position], y_limit * 0.88)
                texts.append(ax.text(
                    x_value, y_value, frame.at[index, '_gene'],
                    ha='left' if x_value >= 0 else 'right', va='bottom',
                    fontsize=profile.tick_font_pt, color=style.text,
                    bbox={'boxstyle': 'round,pad=0.12', 'facecolor': 'white',
                          'edgecolor': style.subtle, 'linewidth': 0.35, 'alpha': 0.92},
                    zorder=5,
                ))
            adjust_labels(texts, ax, arrow_color=style.neutral_dark)
            hidden_labels = prune_overlapping_labels(texts, ax)
            semantic_warnings = []
            if hidden_labels:
                semantic_warnings.append(
                    f'为避免最终尺寸文字重叠，已隐藏 {len(hidden_labels)} 个低优先级基因标签；完整基因列表保留在结果表。'
                )

            if spec.show_legend:
                handles = [
                    Line2D([], [], linestyle='None', marker='o', markersize=3.8,
                           markerfacecolor=style.deg_palette[label], markeredgewidth=0,
                           label=f'{label} (n={counts[label]})')
                    for label in ('Up', 'Down', 'NS')
                ]
                ax.legend(handles=handles, loc='lower right', frameon=False,
                          handletextpad=0.35, borderaxespad=0.2)
            fig.subplots_adjust(left=0.19, right=0.97, bottom=0.17, top=0.90)

        if container is not None:
            return container
        return attach_contract(
            fig, spec, style,
            semantic_warnings=warnings + semantic_warnings,
            encodings={'color': 'DEG direction', 'shape': 'clipped FDR'},
            validation_texts=texts,
        )
