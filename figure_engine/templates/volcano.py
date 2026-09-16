"""NatureVolcano：固定阈值语义和受控基因标注。"""

import re
import textwrap

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import (
    adjust_labels, attach_contract, ensure_angled_annotation_leaders,
    finite_numeric, prune_overlapping_labels,
)


DEFAULT_VOLCANO_X_LIMIT = 2.5
DEFAULT_VOLCANO_Y_LIMIT = 18.0


def _volcano_display_limits(y_raw, threshold_y, fc_threshold=0.0):
    """Return publication-scale axes without changing the reported statistics.

    The main cloud in a Bulk RNA-seq volcano plot is ordinarily most legible
    within ±2.5 log2FC.  Larger observed effects remain in the figure as edge
    triangles rather than expanding the panel or being silently discarded.
    Likewise, the FDR axis stays linear in its displayed range and uses a
    triangular cap above 18 instead of a nonlinear tail transformation.
    """
    # A user may deliberately choose an unusually stringent fold-change gate;
    # keep both threshold lines visible even then, while retaining the compact
    # ±2.5 default for normal Bulk RNA-seq contrasts.
    x_limit = max(
        DEFAULT_VOLCANO_X_LIMIT,
        float(np.ceil((float(fc_threshold) + 0.25) * 2.0) / 2.0),
    )

    finite_y = np.asarray(y_raw, dtype=float)
    finite_y = finite_y[np.isfinite(finite_y)]
    if finite_y.size == 0:
        y_limit = max(4.0, float(threshold_y) + 1.5)
    elif float(np.nanmax(finite_y)) > DEFAULT_VOLCANO_Y_LIMIT:
        y_limit = DEFAULT_VOLCANO_Y_LIMIT
    else:
        y_limit = max(
            4.0,
            float(threshold_y) + 1.5,
            float(np.ceil(float(np.nanmax(finite_y)) + 0.6)),
        )
    return float(x_limit), float(y_limit)


def _wrap_title(value, default, width):
    """Wrap a title while preserving an intentional scientific context line."""
    text = str(value or default).translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789'))
    # Common treatment spelling benefits from a real chemical subscript while
    # retaining the grouping suffix literally (for example ``NH4Cl_B``).
    text = re.sub(r'(?<![A-Za-z0-9])NH4Cl(?![A-Za-z0-9])', r'NH$_4$Cl', text)
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
        # Feature IDs are retained in the results table, but a bare Ensembl
        # identifier is seldom an informative publication label.  Explicitly
        # requested identifiers above are still honoured.
        readable = significant[~significant['_gene'].str.match(
            r'^(?:ENS(?:G|MUSG|DARG|RNOG)\d+|[A-Za-z]+\d{8,})$', na=False,
        )]
        if readable.empty:
            readable = significant

        # FDR alone would repeatedly label a very narrow set of effects.  Keep
        # the strongest statistical evidence from each direction, then use the
        # remaining capacity for the largest effects.  Explicitly requested
        # genes above always remain first-priority labels.
        # With the four-label platform default, show two genes per direction.
        # At six or more labels, reserve the additional slots for the largest
        # observed effects after that balanced FDR-first core.
        fdr_per_side = 2 if remaining >= 4 else 1
        fdr_candidates = pd.concat([
            readable[regulation.loc[readable.index] == 'Up']
            .sort_values(['_fdr', '_abs_fc'], ascending=[True, False]).head(fdr_per_side),
            readable[regulation.loc[readable.index] == 'Down']
            .sort_values(['_fdr', '_abs_fc'], ascending=[True, False]).head(fdr_per_side),
        ])
        effect_candidates = readable.sort_values(
            ['_abs_fc', '_fdr'], ascending=[False, True],
        )
        for candidates in (fdr_candidates, effect_candidates):
            for index in candidates.index:
                if index not in selected:
                    selected.append(index)
                if len(selected) >= int(spec.label_n):
                    break
            if len(selected) >= int(spec.label_n):
                break
    return selected


def _label_offsets(label_rows, y_limit, x_limit):
    """Start labels with short angled leaders; ``adjustText`` handles overlap."""
    if not label_rows:
        return []
    placed = []
    for side in (-1, 1):
        rows = [row for row in label_rows if row['side'] == side]
        for rank, row in enumerate(sorted(rows, key=lambda item: item['y_point'])):
            # Labels near an x/y edge point inwards so their initial short
            # leader stays on-canvas.  The small rank offset breaks exact ties
            # before the final repel pass, without becoming a horizontal rail.
            x_side = side
            if abs(float(row['x_point'])) >= float(x_limit) * 0.78:
                x_side *= -1
            y_side = -1 if float(row['y_point']) >= float(y_limit) * 0.84 else 1
            placed.append({
                **row,
                # Data-coordinate starts let adjustText enforce the axes
                # boundary correctly; Annotation + offset-points positions do
                # not share its movement transform reliably.
                'x_text': float(row['x_point']) + x_side * float(x_limit) * 0.08,
                'y_text': float(row['y_point']) + y_side * float(y_limit) * (0.05 + 0.02 * (rank // 2)),
                'label_side': x_side,
                'y_side': y_side,
            })
    return placed


def _place_legend_in_clearest_corner(ax, fig, legend, x_values, y_values, labels):
    """Choose the least occupied corner for a compact in-panel legend.

    Bulk contrasts can put their strongest signal on either side of the
    volcano.  A fixed corner either covers that signal or, worse, covers a
    label rail.  Evaluate the four conventional corners against the rendered
    point cloud and label boxes so the legend remains informative but does not
    compete with the data.
    """
    finite = np.isfinite(x_values) & np.isfinite(y_values)
    points = ax.transData.transform(np.column_stack((x_values[finite], y_values[finite])))
    candidates = ('upper left', 'upper right', 'lower right', 'lower left')
    best_location = candidates[0]
    best_score = None
    for location in candidates:
        try:
            legend.set_loc(location)
        except AttributeError:  # Matplotlib < 3.8 compatibility
            legend._loc = location
        fig.canvas.draw()
        legend_box = legend.get_window_extent(fig.canvas.get_renderer())
        point_hits = int(sum(legend_box.contains(x, y) for x, y in points))
        label_hits = 0
        for label in labels:
            patch = label.get_bbox_patch() if hasattr(label, 'get_bbox_patch') else None
            if patch is not None and patch.get_visible():
                label_box = patch.get_window_extent(fig.canvas.get_renderer())
            else:
                from matplotlib.text import Text
                label_box = Text.get_window_extent(label, fig.canvas.get_renderer())
            label_hits += int(legend_box.overlaps(label_box))
        # A readable selected gene is more valuable than an unlabelled
        # background point, so it carries a stronger placement penalty.
        score = point_hits + 20 * label_hits
        if best_score is None or score < best_score:
            best_location, best_score = location, score
    try:
        legend.set_loc(best_location)
    except AttributeError:  # Matplotlib < 3.8 compatibility
        legend._loc = best_location
    # FigureValidator can distinguish this measured in-panel placement from a
    # fixed legend that may silently cover data.
    legend._nature_data_aware_placement = True


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
        x_raw = frame['_log2fc'].to_numpy(dtype=float)
        y_raw = -np.log10(frame['_fdr'].to_numpy(dtype=float))
        threshold_y = -np.log10(spec.fdr_threshold)
        x_limit, y_limit = _volcano_display_limits(y_raw, threshold_y, spec.fc_threshold)
        x_clipped = np.abs(x_raw) > x_limit
        y_clipped = y_raw > y_limit
        x_display = np.clip(x_raw, -x_limit * 0.985, x_limit * 0.985)
        y_display = np.minimum(y_raw, y_limit * 0.985)

        warnings = []
        if len(frame) > 100000:
            warnings.append('Volcano 点数超过 100,000；SVG 中的数据层将栅格化以控制文件体积。')
        # Six automatic labels (three per direction) remain readable in an
        # 89-mm Bulk RNA-seq panel.  More labels belong in the full DEG table
        # or can be selected deliberately in Figure Studio.
        max_labels = 6 if spec.width == 'single' else 10
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
                ('NS', 3.5 if spec.width == 'single' else 4.5, 0.26, 1),
                ('Down', 9.0 if spec.width == 'single' else 11.0, 0.90, 2),
                ('Up', 9.0 if spec.width == 'single' else 11.0, 0.90, 2),
            ):
                mask = (regulation.to_numpy() == label) & ~(x_clipped | y_clipped)
                if mask.any():
                    ax.scatter(
                        x_display[mask], y_display[mask],
                        s=point_size, color=style.deg_palette[label],
                        alpha=alpha, linewidths=0, rasterized=True, zorder=zorder,
                    )
            # Edge triangles preserve the direction of truncated points while
            # leaving the compact central cloud readable.
            for label in ('NS', 'Down', 'Up'):
                label_mask = regulation.to_numpy() == label
                vertical = label_mask & y_clipped
                if vertical.any():
                    ax.scatter(
                        x_display[vertical], y_display[vertical], s=17,
                        marker='^', color=style.deg_palette[label], alpha=0.92,
                        linewidths=0, rasterized=True, zorder=4,
                    )
                left = label_mask & x_clipped & ~y_clipped & (x_raw < 0)
                right = label_mask & x_clipped & ~y_clipped & (x_raw > 0)
                if left.any():
                    ax.scatter(
                        x_display[left], y_display[left], s=17,
                        marker='<', color=style.deg_palette[label], alpha=0.92,
                        linewidths=0, rasterized=True, zorder=4,
                    )
                if right.any():
                    ax.scatter(
                        x_display[right], y_display[right], s=17,
                        marker='>', color=style.deg_palette[label], alpha=0.92,
                        linewidths=0, rasterized=True, zorder=4,
                    )
            threshold_color = style.subtle
            ax.axhline(threshold_y, color=threshold_color, linestyle=(0, (3, 2)),
                       linewidth=0.55, alpha=0.72, zorder=0)
            ax.axvline(spec.fc_threshold, color=threshold_color, linestyle=(0, (3, 2)),
                       linewidth=0.55, alpha=0.72, zorder=0)
            ax.axvline(-spec.fc_threshold, color=threshold_color, linestyle=(0, (3, 2)),
                       linewidth=0.55, alpha=0.72, zorder=0)
            ax.set_xlim(-x_limit, x_limit)
            ax.set_ylim(0, y_limit)
            ax.set_xlabel(r'log$_2$(fold change)')
            ax.set_ylabel(r'$-\log_{10}(\mathrm{FDR})$')
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
            selected = _select_labels(frame, label_spec, regulation)[:max_labels]
            label_rows = []
            for index in selected:
                row_position = frame.index.get_loc(index)
                x_value = x_display[row_position]
                label_rows.append({
                    'gene': frame.at[index, '_gene'],
                    'x_point': float(x_value),
                    'y_point': min(float(y_display[row_position]), y_limit * 0.88),
                    'side': 1 if x_value >= 0 else -1,
                })
            texts = []
            for label in _label_offsets(label_rows, y_limit, x_limit):
                texts.append(ax.annotate(
                    label['gene'],
                    xy=(label['x_point'], label['y_point']),
                    xytext=(label['x_text'], label['y_text']), textcoords='data',
                    ha='left' if label['label_side'] > 0 else 'right',
                    va='bottom' if label['y_side'] > 0 else 'top',
                    fontsize=profile.tick_font_pt, color=style.text,
                    arrowprops={'arrowstyle': '-', 'color': style.neutral_dark,
                                'linewidth': 0.45, 'shrinkA': 1, 'shrinkB': 1},
                    zorder=5,
                ))
            adjust_labels(texts, ax, arrow_color=style.neutral_dark)
            ensure_angled_annotation_leaders(texts, y_limit)
            hidden_labels = prune_overlapping_labels(texts, ax)
            semantic_warnings = []
            if hidden_labels:
                semantic_warnings.append(
                    f'为避免最终尺寸文字重叠，已隐藏 {len(hidden_labels)} 个低优先级基因标签；完整基因列表保留在结果表。'
                )

            if spec.show_legend:
                show_legend_counts = bool(spec.extra.get('legend_counts', True))
                handles = [
                    Line2D([], [], linestyle='None', marker='o', markersize=3.8,
                           markerfacecolor=style.deg_palette[label], markeredgewidth=0,
                           label=(f'{label}  {counts[label]:,}'
                                  if show_legend_counts else label))
                    for label in ('Down', 'NS', 'Up')
                ]
                legend = ax.legend(
                    handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.18),
                    frameon=False, ncol=3, handletextpad=0.30,
                    columnspacing=0.85, borderaxespad=0.0,
                )
            title_lines = wrapped_title.count('\n') + 1
            top = max(0.76, 0.90 - 0.055 * (title_lines - 1))
            fig.subplots_adjust(left=0.19, right=0.97,
                                bottom=0.27 if spec.show_legend else 0.17, top=top)

        if container is not None:
            return container
        return attach_contract(
            fig, spec, style,
            semantic_warnings=warnings + semantic_warnings,
            encodings={
                'color': 'DEG direction',
                'x': 'log2(fold change); edge triangles mark values beyond displayed range',
                'y': '-log10(FDR); top triangles mark values above displayed range',
            },
            validation_texts=texts,
        )
