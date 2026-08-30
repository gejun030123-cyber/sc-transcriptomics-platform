"""NatureVolcano：固定阈值语义和受控基因标注。"""

import textwrap

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import adjust_labels, attach_contract, finite_numeric, prune_overlapping_labels


def _compress_fdr_tail(y_raw, threshold_y):
    """Keep an extreme FDR tail legible without flattening it at a ceiling.

    DESeq2 can legitimately report adjusted p-values far below the useful
    display range.  A hard y-axis cap turns all of those points into one dense
    horizontal row, which hides their relative strength and makes the Volcano
    plot look broken.  Above a clearly marked break we therefore use a
    monotonic log compression.  The plotted order is preserved and tick labels
    remain in the original ``-log10(FDR)`` units.
    """
    values = np.asarray(y_raw, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return values, None

    # Values below 1e-12 retain the ordinary, linear volcano scale.  The
    # threshold line (normally 1.30 for FDR 0.05) is consequently never
    # affected by compression.
    break_y = max(12.0, float(np.ceil(float(threshold_y) + 2.0)))
    max_y = float(np.nanmax(finite))
    if max_y <= break_y + 1.0:
        return values, None

    scale = max(4.0, break_y * 0.5)
    displayed = values.copy()
    tail = displayed > break_y
    displayed[tail] = break_y + scale * np.log1p((displayed[tail] - break_y) / scale)
    return displayed, {"break_y": break_y, "scale": scale, "max_y": max_y}


def _tail_tick_values(tail_info):
    """Return sparse, truthful tick locations for a compressed FDR tail."""
    break_y = float(tail_info["break_y"])
    scale = float(tail_info["scale"])
    max_y = float(tail_info["max_y"])

    # Four linear-range ticks leave enough space for three sparse tail ticks
    # on a single-column figure.  The tail ticks have familiar 10^x values,
    # not arbitrary display coordinates.
    raw_ticks = list(np.linspace(0.0, break_y, num=4))
    candidates = np.asarray((15, 20, 30, 50, 75, 100, 150, 200, 300), dtype=float)
    candidates = candidates[(candidates > break_y) & (candidates <= max_y)]
    if len(candidates) > 3:
        candidates = candidates[np.linspace(0, len(candidates) - 1, num=3, dtype=int)]
    raw_ticks.extend(candidates.tolist())
    raw_ticks = np.asarray(list(dict.fromkeys(float(value) for value in raw_ticks)))
    displayed = raw_ticks.copy()
    tail = displayed > break_y
    displayed[tail] = break_y + scale * np.log1p((displayed[tail] - break_y) / scale)
    labels = [f"{value:g}" for value in raw_ticks]
    return displayed, labels


def _wrap_title(value, default, width):
    """Wrap a title while preserving an intentional scientific context line."""
    text = str(value or default).replace('_', ' ')
    lines = [line.strip() for line in text.splitlines() if line.strip()] or [default]
    return '\n'.join(
        textwrap.fill(line, width=width, break_long_words=False)
        for line in lines
    )


def _select_labels(frame, spec, regulation):
    genes = frame['_gene'].astype(str)
    selected = []
    # Follow explicit user order; a custom label is a review request, not an
    # automatic rank, and may intentionally be non-significant.
    for gene in spec.label_genes:
        for index in frame.index[genes.eq(str(gene))]:
            if index not in selected:
                selected.append(index)
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
    """固定的 transcriptomics Volcano publication 模板。"""

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
        threshold_y = -np.log10(spec.fdr_threshold)
        y_display, tail_info = _compress_fdr_tail(y_raw, threshold_y)
        y_limit = max(4.0, threshold_y + 1.5, float(np.ceil(np.nanmax(y_display) + 0.6)))
        x_abs = np.abs(frame['_log2fc'].to_numpy(dtype=float))
        # Use the full observed range so extreme-fold-change points are drawn at
        # their true positions instead of being clipped to a vertical strip on
        # both axes edges (the “两边又跑出画面” volcano artifact).  This mirrors
        # the bulk DEG volcano contract.
        finite_x = x_abs[np.isfinite(x_abs)]
        max_abs_x = float(np.nanmax(finite_x)) if finite_x.size else 0.0
        x_limit = max(spec.fc_threshold * 1.8, max_abs_x * 1.08, 2.5)
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
            ax.set_ylabel(
                r'$-\log_{10}$(FDR) · compressed tail'
                if tail_info else r'$-\log_{10}$(FDR)'
            )
            if tail_info:
                tick_locations, tick_labels = _tail_tick_values(tail_info)
                ax.set_yticks(tick_locations, tick_labels)
            title_width = 38 if spec.width == 'single' else 72
            wrapped_title = _wrap_title(
                spec.title, 'Differential expression', title_width,
            )
            ax.set_title(
                wrapped_title,
                loc='left', pad=5,
            )
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
            title_lines = wrapped_title.count('\n') + 1
            top = max(0.76, 0.90 - 0.055 * (title_lines - 1))
            fig.subplots_adjust(left=0.19, right=0.97, bottom=0.17, top=top)

        if container is not None:
            return container
        return attach_contract(
            fig, spec, style,
            semantic_warnings=warnings + semantic_warnings,
            encodings={
                'color': 'DEG direction',
                'y': ('-log10(FDR), monotonic tail compression' if tail_info
                      else '-log10(FDR)'),
            },
            validation_texts=texts,
        )
