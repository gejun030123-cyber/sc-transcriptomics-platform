"""Native Matplotlib builders for publication/display single-cell figures.

These helpers deliberately contain no Plotly conversion code.  They are used by
the display-facing analysis paths so the figure's source, preview, and SVG
export are all the same Matplotlib canvas.
"""

import numpy as np

from modules.figure_style import (
    NATURE_AXIS,
    NATURE_GRID,
    NATURE_PALETTE,
    NATURE_TEXT,
    nature_continuous_cmap,
    nature_rcparams,
    register_nature_cjk_font,
    stable_category_colors,
)
from figure_engine.templates.common import adjust_labels

# The categorical bar palette is intentionally softer than the general-purpose
# palette. Paired with a dark outline, it stays readable without making the
# bars visually heavier than their labels.
BAR_PALETTE = (
    '#7BA9C8', '#F07B73', '#B7B2CF', '#79B4A4', '#D8AD62', '#9AA6B2',
)

# Apply the shared font contract before any native canvas is constructed.  This
# prevents tight_layout from resolving CJK labels against DejaVu Sans before
# BaseAnalysis.save_matplotlib_figure has a chance to normalize the figure.
try:
    import matplotlib as _mpl
    register_nature_cjk_font()
    _mpl.rcParams.update(nature_rcparams())
except Exception:
    pass


def _as_categories(values):
    if hasattr(values, 'cat'):
        return values.astype('category')
    import pandas as pd

    return pd.Series(values).astype('category')


def _style_axis(ax):
    ax.grid(False)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, colors=NATURE_AXIS, width=0.7, length=3)
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name in ('left', 'bottom'))
        spine.set_color('#98A2B3')
        spine.set_linewidth(0.7)
    return ax


def _style_bar_axis(ax, *, show_grid=False):
    """Apply the clean, outline-led bar-chart treatment used in result views."""
    _style_axis(ax)
    for spine_name in ('left', 'bottom'):
        spine = ax.spines[spine_name]
        spine.set_color('#15191E')
        spine.set_linewidth(1.25)
    ax.tick_params(colors='#15191E', width=1.05, length=4)
    if show_grid:
        ax.grid(axis='both', color=NATURE_GRID, linewidth=0.55, alpha=0.88, zorder=0)
    return ax


def sample_tick_plan(labels, max_labels=18):
    """Choose readable, transparent tick density for sample-labelled panels.

    Dense sample IDs need not be drawn on every tick to preserve the data
    matrix.  The returned positions always include the first and last sample;
    the complete ordered sample list remains available in the source data and
    correlation metadata.
    """
    labels = [str(label) for label in labels]
    n_labels = len(labels)
    if n_labels == 0:
        return {'indices': [], 'labels': [], 'stride': 1, 'font_size': 8.0}
    if max_labels is None or n_labels <= max(1, int(max_labels)):
        indices = list(range(n_labels))
        stride = 1
    else:
        stride = int(np.ceil(n_labels / max(1, int(max_labels))))
        indices = list(range(0, n_labels, stride))
        if indices[-1] != n_labels - 1:
            indices.append(n_labels - 1)
    longest = max(len(label) for label in labels)
    font_size = 8.0 if len(indices) <= 12 else (7.0 if len(indices) <= 18 else 6.5)
    if longest > 14:
        font_size = min(font_size, 6.5)
    return {
        'indices': indices,
        'labels': [labels[index] for index in indices],
        'stride': stride,
        'font_size': font_size,
    }


def apply_sample_tick_labels(ax, labels, axis='x', max_labels=18, positions=None,
                             rotation=45, font_size=None):
    """Apply a sample tick plan to an x/y axis and return its provenance."""
    plan = sample_tick_plan(labels, max_labels=max_labels)
    base_positions = np.arange(len(labels)) if positions is None else np.asarray(positions)
    selected_positions = [base_positions[index] for index in plan['indices']]
    effective_font_size = min(float(font_size), plan['font_size']) if font_size else plan['font_size']
    if axis == 'x':
        ax.set_xticks(selected_positions, plan['labels'], rotation=rotation, ha='right')
        ax.tick_params(axis='x', labelsize=effective_font_size)
    elif axis == 'y':
        ax.set_yticks(selected_positions, plan['labels'])
        ax.tick_params(axis='y', labelsize=effective_font_size)
    else:
        raise ValueError("axis must be 'x' or 'y'.")
    return {**plan, 'font_size': effective_font_size}


def plot_umap_axis(ax, adata, color_key, title='', basis='X_umap',
                   point_size=7, opacity=0.78, show_legend=True,
                   label_categories=False, hide_axes=True):
    """Draw one Nature-style UMAP panel on an existing axis."""
    coords = np.asarray(adata.obsm[basis])[:, :2]
    values = adata.obs[color_key] if color_key in adata.obs.columns else None
    import pandas as pd
    label_texts = []

    if values is None:
        ax.scatter(coords[:, 0], coords[:, 1], s=point_size,
                   color=NATURE_PALETTE[0], alpha=opacity,
                   linewidths=0, rasterized=True)
    elif pd.api.types.is_numeric_dtype(values):
        artist = ax.scatter(
            coords[:, 0], coords[:, 1], s=point_size,
            c=np.asarray(values, dtype=float), cmap=nature_continuous_cmap(),
            alpha=opacity, linewidths=0, rasterized=True,
        )
        colorbar = ax.figure.colorbar(artist, ax=ax, fraction=0.035,
                                      pad=0.025, aspect=32)
        colorbar.outline.set_visible(False)
        colorbar.ax.tick_params(labelsize=7, width=0.6, colors=NATURE_AXIS)
        colorbar.set_label(str(color_key), fontsize=8, labelpad=5)
    else:
        categorical = _as_categories(values)
        color_map = stable_category_colors(
            [str(category) for category in categorical.cat.categories],
            existing=adata.uns.get(f'{color_key}_colors', []),
        )
        for index, category in enumerate(categorical.cat.categories):
            mask = np.asarray(categorical == category)
            if not mask.any():
                continue
            ax.scatter(
                coords[mask, 0], coords[mask, 1], s=point_size,
                color=color_map.get(str(category), NATURE_PALETTE[index % len(NATURE_PALETTE)]),
                alpha=opacity, linewidths=0, rasterized=True,
                label=str(category),
            )
            if label_categories:
                label_texts.append(ax.text(
                    float(np.median(coords[mask, 0])),
                    float(np.median(coords[mask, 1])),
                    str(category), ha='center', va='center', fontsize=8,
                    color=NATURE_TEXT,
                    bbox={'boxstyle': 'round,pad=0.18', 'facecolor': 'white',
                          'edgecolor': '#D0D5DD', 'alpha': 0.82, 'linewidth': 0.5},
                    zorder=5,
                ))
        if show_legend and len(categorical.cat.categories) and not label_categories:
            legend = ax.legend(
                loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False,
                fontsize=7, markerscale=1.15, borderaxespad=0,
                title=str(color_key),
            )
            legend.get_title().set_fontsize(8)

    if label_texts:
        adjust_labels(label_texts, ax, arrow_color='#98A2B3')
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel('UMAP 1', fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel('UMAP 2', fontsize=9, color=NATURE_TEXT)
    ax.margins(0.035)
    _style_axis(ax)
    if hide_axes:
        ax.set_xticks([])
        ax.set_yticks([])
    return ax


def umap_figure(adata, color_key, title='', basis='X_umap',
                point_size=7, opacity=0.78, label_categories=False,
                hide_axes=True):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=150)
    plot_umap_axis(
        ax, adata, color_key, title=title, basis=basis,
        point_size=point_size, opacity=opacity,
        label_categories=label_categories, hide_axes=hide_axes,
    )
    import pandas as pd
    is_numeric = pd.api.types.is_numeric_dtype(adata.obs[color_key])
    fig.subplots_adjust(left=0.10, right=0.80 if not is_numeric and not label_categories else 0.93,
        bottom=0.11, top=0.88)
    return fig


def umap_panel_figure(adata, color_keys, titles=None, basis='X_umap',
                      point_size=5, opacity=0.72, ncols=3, hide_axes=True):
    import matplotlib.pyplot as plt

    color_keys = [key for key in color_keys if key in adata.obs.columns]
    if not color_keys:
        return None
    ncols = max(1, min(int(ncols), len(color_keys)))
    nrows = int(np.ceil(len(color_keys) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.4 * ncols, 4.0 * nrows),
        dpi=150, squeeze=False,
    )
    axes_flat = axes.ravel()
    titles = titles or [str(key) for key in color_keys]
    for index, key in enumerate(color_keys):
        plot_umap_axis(
            axes_flat[index], adata, key,
            title=titles[index] if index < len(titles) else str(key),
            basis=basis, point_size=point_size, opacity=opacity,
            show_legend=False, hide_axes=hide_axes,
        )
    for ax in axes_flat[len(color_keys):]:
        ax.set_visible(False)
    fig.tight_layout(pad=1.1)
    return fig


def heatmap_figure(matrix, x_labels=None, y_labels=None, title='',
                   x_label='', y_label='', colorbar_label='z-score',
                   vmin=-3, vmax=3, x_group_labels=None,
                   x_label_rotation=45):
    """Render a compact heatmap with an optional second-level x-axis header.

    ``x_group_labels`` is useful when adjacent columns share a long parent
    label (for example, cell type) while the x ticks carry a short child label
    (for example, condition).  Keeping the two levels separate prevents the
    diagonal, overlapping labels that make dense biological heatmaps hard to
    read.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    matrix = np.asarray(matrix, dtype=float)
    cmap = LinearSegmentedColormap.from_list(
        'nature_diverging', ['#0F4D92', '#DCEAF0', '#FFFFFF', '#F6CFCB', '#B64342']
    )
    grouped_columns = x_group_labels is not None
    if grouped_columns and len(x_group_labels) != matrix.shape[1]:
        raise ValueError('x_group_labels must match the heatmap column count')
    width_per_column = 0.78 if grouped_columns else 0.22
    width = max(7.0, min(16.0, 4.8 + width_per_column * max(1, matrix.shape[1])))
    height = max(4.8, min(12.0, 3.6 + 0.18 * max(1, matrix.shape[0])))
    fig, ax = plt.subplots(figsize=(width, height), dpi=150)
    image = ax.imshow(matrix, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    if x_labels is not None:
        ax.set_xticks(np.arange(len(x_labels)), [str(item) for item in x_labels])
        ax.tick_params(axis='x', labelrotation=x_label_rotation,
                       labelsize=8, pad=5)
    if y_labels is not None:
        ax.set_yticks(np.arange(len(y_labels)), [str(item) for item in y_labels])
        ax.tick_params(axis='y', labelsize=8)
    if grouped_columns:
        # Place the title at figure level, leaving a dedicated row for the
        # grouped headers rather than letting both strings occupy the same
        # narrow top margin.
        fig.suptitle(title, x=0.01, y=0.99, ha='left', fontsize=10,
                     fontweight='semibold', color=NATURE_TEXT)
        groups = [str(item) for item in x_group_labels]
        start = 0
        for index in range(1, len(groups) + 1):
            if index != len(groups) and groups[index] == groups[start]:
                continue
            end = index - 1
            centre = (start + end) / 2
            ax.text(centre, 1.015, groups[start], transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=7.4, fontweight='semibold',
                    color=NATURE_TEXT, clip_on=False, linespacing=0.92)
            if start:
                ax.axvline(start - 0.5, color='#98A2B3', linewidth=0.75, alpha=0.85)
            start = index
    else:
        ax.set_title(title, loc='left', pad=10, fontsize=10,
                     fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025, aspect=32)
    colorbar.outline.set_visible(False)
    colorbar.set_label(colorbar_label, fontsize=8, labelpad=5)
    _style_axis(ax)
    fig.tight_layout(pad=1.1, rect=(0.0, 0.0, 1.0, 0.94) if grouped_columns else None)
    return fig


def marker_violin_figure(values_by_category, categories, genes, title='',
                         y_label='Expression (display scale)'):
    """Render compact per-gene violin panels for one cluster-vs-rest contrast.

    ``values_by_category`` is ordered as ``[category][gene]``.  Keeping the
    input as arrays rather than an AnnData object makes the renderer reusable
    and lets callers deterministically subsample large cell collections before
    drawing a static figure.
    """
    import matplotlib.pyplot as plt

    categories = [str(category) for category in categories]
    genes = [str(gene) for gene in genes]
    if not categories or not genes or len(values_by_category) != len(categories):
        raise ValueError('Marker violin requires matching categories and genes.')

    n_panels = len(genes)
    ncols = min(3, max(1, n_panels))
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(max(7.2, 3.35 * ncols), max(3.9, 3.25 * nrows + 0.55)),
        dpi=150,
        squeeze=False,
    )
    colors = [NATURE_PALETTE[index % len(NATURE_PALETTE)]
              for index in range(len(categories))]

    for gene_index, gene in enumerate(genes):
        axis = axes.flat[gene_index]
        series = []
        for category_index in range(len(categories)):
            raw_values = np.asarray(values_by_category[category_index][gene_index], dtype=float)
            finite_values = raw_values[np.isfinite(raw_values)]
            # matplotlib's KDE estimator needs at least two observations.  A
            # single observed cell is still informative, so draw it as a thin
            # distribution instead of discarding the whole marker panel.
            if finite_values.size == 1:
                finite_values = np.repeat(finite_values, 2)
            elif finite_values.size == 0:
                finite_values = np.zeros(2, dtype=float)
            series.append(finite_values)
        parts = axis.violinplot(
            series, positions=np.arange(1, len(categories) + 1),
            showmedians=True, showextrema=False, widths=0.78,
        )
        for body_index, body in enumerate(parts['bodies']):
            body.set_facecolor(colors[body_index])
            body.set_edgecolor(colors[body_index])
            body.set_alpha(0.78)
        for key in ('cmedians',):
            if key in parts:
                parts[key].set_color('#15191E')
                parts[key].set_linewidth(1.0)
        axis.set_xticks(np.arange(1, len(categories) + 1), categories,
                         rotation=28, ha='right')
        axis.set_title(gene, loc='left', pad=7, fontsize=9,
                       fontweight='semibold', color=NATURE_TEXT)
        axis.set_ylabel(y_label, fontsize=8, color=NATURE_TEXT)
        axis.grid(axis='y', color=NATURE_GRID, linewidth=0.55, alpha=0.72)
        _style_axis(axis)

    for axis in axes.flat[n_panels:]:
        axis.set_visible(False)
    fig.suptitle(title, x=0.01, ha='left', y=0.995, fontsize=11,
                 fontweight='semibold', color=NATURE_TEXT)
    layout_rect = (0.0, 0.0, 1.0, 0.93)
    fig._native_layout_rect = layout_rect
    fig.tight_layout(rect=layout_rect, pad=1.1)
    return fig


def _correlation_display_limits(correlation_matrix):
    """Choose an interpretable display scale from non-diagonal correlations.

    The returned limits only control colour mapping.  They never alter the
    underlying correlation matrix or the pairwise values exported alongside it.
    """
    matrix = np.asarray(correlation_matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError('Correlation matrix must be square.')
    off_diagonal = matrix[~np.eye(matrix.shape[0], dtype=bool)]
    finite = off_diagonal[np.isfinite(off_diagonal)]
    if finite.size == 0:
        return 0.0, 1.0, {'min': np.nan, 'median': np.nan, 'max': np.nan}

    actual_min = float(np.nanmin(finite))
    actual_max = float(np.nanmax(finite))
    stats = {
        'min': actual_min,
        'median': float(np.nanmedian(finite)),
        'max': actual_max,
    }
    if actual_min < 0:
        # Negative correlations require an honest diverging scale around zero.
        vmin = max(-1.0, np.floor(actual_min * 20.0) / 20.0)
    else:
        # Keep an absolute 0.95 lower bound when samples are almost identical,
        # so tiny numerical differences are not visually overstated.  Otherwise
        # round down to a 0.05 step to expose meaningful QC differences.
        vmin = max(0.0, min(0.95, np.floor(actual_min * 20.0) / 20.0))
    if vmin >= 1.0:
        vmin = 0.95
    return float(vmin), 1.0, stats


def _correlation_cluster_order(correlation_matrix):
    """Return a stable average-linkage leaf order based on 1 - correlation."""
    matrix = np.asarray(correlation_matrix, dtype=float)
    n_samples = matrix.shape[0]
    if n_samples <= 2:
        return list(range(n_samples))

    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    finite = matrix[np.isfinite(matrix) & ~np.eye(n_samples, dtype=bool)]
    fallback = float(np.nanmin(finite)) if finite.size else 0.0
    safe_matrix = np.nan_to_num(matrix, nan=fallback, posinf=fallback, neginf=fallback)
    distance = np.clip(1.0 - safe_matrix, 0.0, 2.0)
    distance = (distance + distance.T) / 2.0
    np.fill_diagonal(distance, 0.0)
    try:
        linkage_matrix = linkage(
            squareform(distance, checks=False), method='average', optimal_ordering=True,
        )
        return leaves_list(linkage_matrix).astype(int).tolist()
    except ValueError:
        return list(range(n_samples))


def _correlation_cmap(name, includes_negative=False):
    """Return a named, colourblind-safe-enough correlation colour scale."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if includes_negative or name == 'RdBu_r':
        cmap = plt.get_cmap('RdBu_r').copy()
    elif name == 'viridis':
        cmap = plt.get_cmap('viridis').copy()
    elif name == 'YlOrRd':
        cmap = plt.get_cmap('YlOrRd').copy()
    else:
        cmap = LinearSegmentedColormap.from_list(
            'nature_correlation_blues', ['#F7FAFC', '#DCEAF0', '#8FA9C9', '#0F4D92']
        )
    cmap.set_bad('#F2F4F7')
    return cmap


def correlation_heatmap_figure(correlation_matrix, sample_labels, title='',
                               method='Pearson', group_labels=None,
                               colorscale='Blues', cluster=True,
                               mask_diagonal=True):
    """Draw an ordered sample-correlation heatmap with transparent scaling.

    The diagonal can be greyed out because self-correlation is always one and
    otherwise dominates a high-correlation QC panel.  Sample order and display
    limits are returned for reproducible reporting.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Rectangle

    matrix = np.asarray(correlation_matrix, dtype=float)
    labels = [str(label) for label in sample_labels]
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError('Correlation matrix must be square.')
    if matrix.shape[0] != len(labels):
        raise ValueError('The number of sample labels must match the correlation matrix.')

    n_samples = matrix.shape[0]
    valid_groups = group_labels is not None and len(group_labels) == n_samples
    groups = [str(group) for group in group_labels] if valid_groups else None
    order = _correlation_cluster_order(matrix) if cluster else list(range(n_samples))
    ordered_matrix = matrix[np.ix_(order, order)]
    ordered_labels = [labels[index] for index in order]
    ordered_groups = [groups[index] for index in order] if groups is not None else None
    display_vmin, display_vmax, stats = _correlation_display_limits(ordered_matrix)

    width = max(7.4, min(14.0, 5.0 + 0.24 * max(1, n_samples)))
    height = max(5.2, min(12.0, 4.0 + 0.20 * max(1, n_samples)))
    fig, ax = plt.subplots(figsize=(width, height), dpi=150)

    plot_matrix = ordered_matrix.copy()
    if mask_diagonal:
        np.fill_diagonal(plot_matrix, np.nan)
    image_kwargs = {
        'aspect': 'auto', 'interpolation': 'nearest',
        'cmap': _correlation_cmap(colorscale, includes_negative=display_vmin < 0),
    }
    if display_vmin < 0:
        image_kwargs['norm'] = TwoSlopeNorm(
            vmin=display_vmin, vcenter=0.0, vmax=display_vmax)
    else:
        image_kwargs.update(vmin=display_vmin, vmax=display_vmax)
    image = ax.imshow(plot_matrix, **image_kwargs)
    sample_ticks_x = apply_sample_tick_labels(
        ax, ordered_labels, axis='x', max_labels=18, rotation=55,
    )
    sample_ticks_y = apply_sample_tick_labels(
        ax, ordered_labels, axis='y', max_labels=18,
    )
    ax.set_title(title, loc='left', pad=38 if ordered_groups is not None else 22, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    if np.isfinite(stats['min']):
        ax.text(
            0, 1.065 if ordered_groups is not None else 1.015,
            f'Non-diagonal {method} r: {stats["min"]:.3f}\u2013{stats["max"]:.3f}; '
            f'median {stats["median"]:.3f}',
            transform=ax.transAxes, ha='left', va='bottom', fontsize=7.5,
            color='#667085',
        )
    x_axis_label = 'Sample'
    if sample_ticks_x['stride'] > 1:
        x_axis_label += f" (labels every {sample_ticks_x['stride']} samples)"
    ax.set_xlabel(x_axis_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel('Sample', fontsize=9, color=NATURE_TEXT)
    _style_axis(ax)

    colorbar = fig.colorbar(image, ax=ax, fraction=0.038, pad=0.025, aspect=30)
    colorbar.outline.set_visible(False)
    colorbar.set_label(
        f'{method} r\n(display {display_vmin:.2f}\u2013{display_vmax:.2f})',
        fontsize=8, labelpad=6,
    )
    colorbar.ax.tick_params(labelsize=7, width=0.6, colors=NATURE_AXIS)

    group_color_map = {}
    if ordered_groups is not None:
        unique_groups = list(dict.fromkeys(ordered_groups))
        group_color_map = {
            group: NATURE_PALETTE[index % len(NATURE_PALETTE)]
            for index, group in enumerate(unique_groups)
        }
        # Keep the group annotation on the same axes: this survives layout
        # normalization during export and stays perfectly aligned to samples.
        for index, group in enumerate(ordered_groups):
            ax.add_patch(Rectangle(
                (index - 0.5, 1.008), 1.0, 0.03,
                transform=ax.get_xaxis_transform(), clip_on=False,
                facecolor=group_color_map[group], edgecolor='none', zorder=6,
            ))
        run_start = 0
        for index in range(1, len(ordered_groups) + 1):
            if index == len(ordered_groups) or ordered_groups[index] != ordered_groups[run_start]:
                if index - run_start >= 2:
                    ax.text(
                        (run_start + index - 1) / 2, 1.023, ordered_groups[run_start],
                        transform=ax.get_xaxis_transform(), ha='center', va='center',
                        fontsize=6.5, color='white', clip_on=False, zorder=7,
                    )
                run_start = index
        ax.text(-0.015, 1.023, 'Group', transform=ax.transAxes,
                ha='right', va='center', fontsize=7, color=NATURE_TEXT,
                clip_on=False)

    # Reserve a header band for the title, the numeric correlation range, and
    # (when present) the aligned group strip.  BaseAnalysis reuses this rect
    # during export, so SVG/PNG cannot crop those annotations.
    layout_rect = (0.0, 0.0, 1.0, 0.82 if ordered_groups is not None else 0.90)
    fig._native_layout_rect = layout_rect
    fig.tight_layout(rect=layout_rect, pad=1.1)
    return fig, {
        'sample_order': order,
        'ordered_labels': ordered_labels,
        'sample_label_stride': sample_ticks_x['stride'],
        'sample_labels_shown': len(sample_ticks_x['indices']),
        'ordered_groups': ordered_groups,
        'group_color_map': group_color_map,
        'display_vmin': display_vmin,
        'display_vmax': display_vmax,
        'off_diagonal': stats,
    }


def correlation_pairwise_table(correlation_matrix, sample_labels, group_labels=None,
                               method='pearson'):
    """Return all unique sample pairs and their raw correlation values."""
    import pandas as pd

    matrix = np.asarray(correlation_matrix, dtype=float)
    labels = [str(label) for label in sample_labels]
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] != len(labels):
        raise ValueError('Correlation matrix and sample labels must have matching square dimensions.')
    groups = [str(group) for group in group_labels] if group_labels is not None else None
    if groups is not None and len(groups) != len(labels):
        raise ValueError('The number of group labels must match the sample labels.')

    rows = []
    for first in range(len(labels)):
        for second in range(first + 1, len(labels)):
            group_first = groups[first] if groups is not None else ''
            group_second = groups[second] if groups is not None else ''
            rows.append({
                'sample_1': labels[first],
                'sample_2': labels[second],
                'group_1': group_first,
                'group_2': group_second,
                'pair_type': 'within_group' if groups is not None and group_first == group_second else
                             ('between_group' if groups is not None else 'all_samples'),
                f'{method}_r': float(matrix[first, second]),
            })
    return pd.DataFrame(rows)


def summarize_correlation_pairs(pairwise_table, method='pearson'):
    """Summarize raw pairwise correlations overall and by group relationship."""
    import pandas as pd

    value_column = f'{method}_r'
    if value_column not in pairwise_table.columns:
        raise ValueError(f'Missing correlation column: {value_column}')
    rows = []
    grouped_rows = []
    if not pairwise_table.empty and not (pairwise_table['pair_type'] == 'all_samples').all():
        grouped_rows = [
            (str(pair_type), frame)
            for pair_type, frame in pairwise_table.groupby('pair_type', dropna=False)
        ]
    for label, subset in [('all_pairs', pairwise_table), *grouped_rows]:
        values = pd.to_numeric(subset[value_column], errors='coerce').dropna()
        rows.append({
            'pair_set': label,
            'n_pairs': int(values.shape[0]),
            'mean_r': float(values.mean()) if not values.empty else np.nan,
            'median_r': float(values.median()) if not values.empty else np.nan,
            'min_r': float(values.min()) if not values.empty else np.nan,
            'max_r': float(values.max()) if not values.empty else np.nan,
        })
    return pd.DataFrame(rows)


def bar_figure(labels, values, title='', x_label='', y_label='',
               colors=None, annotations=None, rotation=35):
    """Draw a compact categorical bar chart with legible value callouts.

    Category order is intentionally supplied by the caller: sample, cluster,
    and time-course order are scientific context and must not be silently
    reordered merely for visual effect. The chart instead uses restrained
    spacing, a dark outline, and optional value labels to make the comparisons
    easier to scan.
    """
    import matplotlib.pyplot as plt

    labels = [str(label) for label in labels]
    values = np.asarray(values, dtype=float)
    colors = colors or [NATURE_PALETTE[0]] * len(labels)
    longest_label = max((len(label) for label in labels), default=0)
    label_rotation = max(rotation, 55) if longest_label > 24 else rotation
    fig, ax = plt.subplots(
        figsize=(max(7.0, 0.80 * len(labels) + 4.0),
                 max(5.0, 4.2 + 0.025 * longest_label)),
        dpi=150,
    )
    bars = ax.bar(np.arange(len(labels)), values, width=0.66, color=colors, alpha=0.96,
                  edgecolor='#15191E', linewidth=0.85, zorder=3)
    if annotations is not None:
        finite = np.abs(values[np.isfinite(values)])
        value_range = max(float(finite.max()) if finite.size else 0.0, 1.0)
        offset = value_range * 0.018
        for bar, annotation in zip(bars, annotations):
            height = bar.get_height()
            direction = 1 if height >= 0 else -1
            ax.text(bar.get_x() + bar.get_width() / 2, height + direction * offset,
                    str(annotation), ha='center', va='bottom' if direction > 0 else 'top', fontsize=8,
                    color=NATURE_TEXT)
    ax.set_xticks(np.arange(len(labels)), labels, rotation=label_rotation, ha='right')
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    _style_bar_axis(ax)
    ax.margins(x=0.045, y=0.10)
    fig.tight_layout(pad=1.1)
    return fig


def scatter_figure(x, y, title='', x_label='', y_label='', colors=None,
                   labels=None, size=8, alpha=0.78, annotate_top=None):
    """Draw a compact publication scatter plot without an interactive layer."""
    import matplotlib.pyplot as plt

    x = np.asarray(x)
    y = np.asarray(y, dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 5.0), dpi=150)
    if colors is None:
        colors = NATURE_PALETTE[0]
    ax.scatter(x, y, s=size, c=colors, alpha=alpha, linewidths=0,
               rasterized=True)
    if labels is not None and annotate_top:
        order = np.argsort(np.nan_to_num(y, nan=-np.inf))[-int(annotate_top):]
        for index in order:
            if np.isfinite(y[index]):
                ax.annotate(str(labels[index]), (x[index], y[index]),
                            xytext=(4, 4), textcoords='offset points', fontsize=7,
                            color=NATURE_TEXT)
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    _style_axis(ax)
    fig.tight_layout(pad=1.1)
    return fig


def histogram_figure(series, labels=None, title='', x_label='', y_label='Frequency',
                     bins=50, alpha=0.55, colors=None):
    """Draw one or more overlaid histograms for static QC/normalization output."""
    import matplotlib.pyplot as plt

    if isinstance(series, np.ndarray) or not isinstance(series, (list, tuple)):
        series = [series]
    labels = list(labels or [f'Series {i + 1}' for i in range(len(series))])
    colors = list(colors or [NATURE_PALETTE[i % len(NATURE_PALETTE)]
                             for i in range(len(series))])
    fig, ax = plt.subplots(figsize=(7.5, 5.0), dpi=150)
    for index, values in enumerate(series):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size:
            ax.hist(values, bins=bins, alpha=alpha, color=colors[index % len(colors)],
                    label=labels[index] if index < len(labels) else f'Series {index + 1}',
                    edgecolor='white', linewidth=0.25)
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    if len(series) > 1:
        ax.legend(frameon=False, fontsize=8)
    ax.grid(axis='y', color=NATURE_GRID, linewidth=0.55, alpha=0.72)
    _style_axis(ax)
    fig.tight_layout(pad=1.1)
    return fig


def line_figure(x, series, title='', x_label='', y_label='', labels=None,
                colors=None, markers=None):
    """Draw one or more static lines for trend/variance summaries."""
    import matplotlib.pyplot as plt

    labels = list(labels or [f'Series {i + 1}' for i in range(len(series))])
    colors = list(colors or [NATURE_PALETTE[i % len(NATURE_PALETTE)]
                             for i in range(len(series))])
    fig, ax = plt.subplots(figsize=(7.5, 5.0), dpi=150)
    for index, values in enumerate(series):
        ax.plot(np.asarray(x), np.asarray(values, dtype=float),
                color=colors[index % len(colors)], linewidth=1.8,
                marker=(markers[index] if markers else None), markersize=3.5,
                label=labels[index] if index < len(labels) else f'Series {index + 1}')
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    if len(series) > 1:
        ax.legend(frameon=False, fontsize=8)
    ax.grid(axis='y', color=NATURE_GRID, linewidth=0.55, alpha=0.72)
    _style_axis(ax)
    fig.tight_layout(pad=1.1)
    return fig


def pie_figure(labels, values, title='', donut=False):
    """Draw a compact static composition chart."""
    import matplotlib.pyplot as plt

    labels = [str(label) for label in labels]
    values = np.asarray(values, dtype=float)
    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=150)
    ax.pie(values, labels=labels, colors=[NATURE_PALETTE[i % len(NATURE_PALETTE)]
                                          for i in range(len(labels))],
           startangle=90, counterclock=False,
           autopct='%1.1f%%', pctdistance=0.78 if donut else 0.68,
           wedgeprops={'width': 0.42 if donut else 1.0, 'edgecolor': 'white',
                       'linewidth': 0.7}, textprops={'fontsize': 8})
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    fig.tight_layout(pad=1.1)
    return fig


def grouped_bar_figure(labels, series, title='', x_label='', y_label='',
                       rotation=35, stacked=False, max_series=12):
    import matplotlib.pyplot as plt

    labels = [str(label) for label in labels]
    series = [(str(name), np.asarray(values, dtype=float)) for name, values in series]
    # A batch column can contain hundreds of technical labels.  Drawing every
    # one creates a vertical legend taller than the figure itself.  Preserve
    # the dominant categories and pool the remainder into an explicit Other
    # series so the composition remains interpretable and bounded.
    if max_series and len(series) > int(max_series):
        ranked = sorted(series, key=lambda item: float(np.nansum(item[1])), reverse=True)
        keep_n = max(1, int(max_series) - 1)
        kept = ranked[:keep_n]
        other_values = np.zeros(len(labels), dtype=float)
        for _, values in ranked[keep_n:]:
            other_values += np.nan_to_num(values, nan=0.0)
        series = kept + [('Other', other_values)]
    longest_label = max((len(label) for label in labels), default=0)
    label_rotation = max(rotation, 55) if longest_label > 24 else rotation
    fig, ax = plt.subplots(
        figsize=(max(7.0, 0.80 * len(labels) + 4.0),
                 max(5.0, 4.2 + 0.025 * longest_label)),
        dpi=150,
    )
    x = np.arange(len(labels))
    width = min(0.74 / max(1, len(series)), 0.30)
    bottoms = np.zeros(len(labels), dtype=float)
    for index, (name, values) in enumerate(series):
        offset = 0 if stacked else (index - (len(series) - 1) / 2) * width
        ax.bar(x + offset, values, width=width if not stacked else 0.70,
               color=BAR_PALETTE[index % len(BAR_PALETTE)],
               alpha=0.96, edgecolor='#15191E', linewidth=0.75, label=name,
               bottom=bottoms if stacked else None)
        if stacked:
            bottoms += values
    ax.set_xticks(x, labels, rotation=label_rotation, ha='right')
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    if series:
        ax.legend(frameon=False, fontsize=7, ncol=min(4, len(series)),
                  loc='upper center', bbox_to_anchor=(0.5, 1.08),
                  borderaxespad=0.0, handlelength=1.2, columnspacing=1.0)
    _style_bar_axis(ax)
    ax.margins(x=0.045, y=0.08)
    fig.tight_layout(pad=1.1, rect=(0, 0, 1, 0.94) if series else None)
    return fig


def diverging_bar_figure(labels, up_values, down_values, title='', x_label='',
                          y_label='Number of DEGs', rotation=35):
    """Render up/down counts as an explicit zero-centred diverging bar chart.

    Down-regulated counts are plotted below the zero baseline solely for the
    display. Callers keep their original non-negative counts in tables and
    summaries, avoiding ambiguity about the underlying statistic.
    """
    import matplotlib.pyplot as plt

    labels = [str(label) for label in labels]
    up_values = np.asarray(up_values, dtype=float)
    down_values = np.asarray(down_values, dtype=float)
    if len(labels) != len(up_values) or len(labels) != len(down_values):
        raise ValueError('labels、up_values 和 down_values 的长度必须一致。')
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(max(7.0, 0.60 * len(labels) + 4.0), 5.4), dpi=150)
    ax.bar(x, -down_values, width=0.72, color='#7299C0', alpha=0.96,
           edgecolor='#15191E', linewidth=0.75, label='Down-regulated', zorder=3)
    ax.bar(x, up_values, width=0.72, color='#F07B73', alpha=0.96,
           edgecolor='#15191E', linewidth=0.75, label='Up-regulated', zorder=3)
    finite = np.abs(np.concatenate((up_values, down_values)))
    limit = max(float(np.nanmax(finite)) if finite.size else 0.0, 1.0)
    ax.set_ylim(-limit * 1.14, limit * 1.14)
    ax.axhline(0, color='#15191E', linewidth=0.85, zorder=4)
    ax.set_xticks(x, labels, rotation=rotation, ha='right')
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    ax.legend(frameon=False, fontsize=8, loc='upper right', handlelength=1.15)
    _style_bar_axis(ax, show_grid=True)
    ax.margins(x=0.045)
    fig.tight_layout(pad=1.1)
    return fig


def grouped_box_figure(values, categories, series_labels, title='',
                       x_label='', y_label='', rotation=35):
    """Draw grouped boxplots for a matrix of observations x measurements."""
    import matplotlib.pyplot as plt

    values = np.asarray(values, dtype=float)
    categories = _as_categories(categories)
    category_labels = [str(item) for item in categories.cat.categories]
    series_labels = [str(item) for item in series_labels]
    n_categories = len(category_labels)
    n_series = max(1, len(series_labels))
    fig, ax = plt.subplots(
        figsize=(max(8.0, 0.8 * n_categories + 4.5), 5.2), dpi=150
    )
    positions = np.arange(n_categories)
    width = min(0.72 / n_series, 0.18)
    for series_index, label in enumerate(series_labels):
        data = []
        valid_positions = []
        for category_index, category in enumerate(category_labels):
            mask = np.asarray(categories == category)
            column = values[mask, series_index]
            column = column[np.isfinite(column)]
            if column.size:
                data.append(column)
                valid_positions.append(
                    positions[category_index]
                    + (series_index - (n_series - 1) / 2) * width
                )
        if not data:
            continue
        color = NATURE_PALETTE[series_index % len(NATURE_PALETTE)]
        boxes = ax.boxplot(
            data, positions=valid_positions, widths=width * 0.9,
            patch_artist=True, showfliers=False,
            boxprops={'facecolor': color, 'alpha': 0.62, 'edgecolor': color},
            medianprops={'color': NATURE_TEXT, 'linewidth': 1.0},
            whiskerprops={'color': color, 'linewidth': 0.8},
            capprops={'color': color, 'linewidth': 0.8},
        )
        for patch in boxes['boxes']:
            patch.set_label(label)
    ax.set_xticks(positions, category_labels, rotation=rotation, ha='right')
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, frameon=False, fontsize=8)
    ax.grid(axis='y', color=NATURE_GRID, linewidth=0.55, alpha=0.72)
    _style_axis(ax)
    fig.tight_layout(pad=1.1)
    return fig


def marker_dotplot_figure(mean_expression, detection_fraction, categories,
                          gene_labels, title='', x_label='Cell type',
                          y_label='Marker gene', value_label='Mean log1p expression',
                          detection_label='% cells expressing', gene_sections=None):
    """Draw an annotation-validation dot plot.

    Columns are annotated groups and rows are marker genes.  Dot colour carries
    the mean displayed expression while dot area carries the fraction of cells
    with detectable expression.  Keeping these two signals separate is more
    informative for zero-inflated single-cell data than a collection of
    side-by-side boxplots, and avoids a legend entry for every individual box.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.colors import Normalize

    mean_expression = np.asarray(mean_expression, dtype=float)
    detection_fraction = np.asarray(detection_fraction, dtype=float)
    if mean_expression.ndim != 2 or detection_fraction.shape != mean_expression.shape:
        raise ValueError('mean_expression and detection_fraction must be 2D arrays with the same shape')
    categories = [str(item) for item in categories]
    gene_labels = [str(item) for item in gene_labels]
    n_categories, n_genes = mean_expression.shape
    if len(categories) != n_categories or len(gene_labels) != n_genes:
        raise ValueError('category/gene labels do not match marker matrix dimensions')
    if not n_categories or not n_genes:
        raise ValueError('marker dotplot requires at least one category and one gene')

    finite_values = mean_expression[np.isfinite(mean_expression)]
    if finite_values.size:
        vmin = float(np.nanpercentile(finite_values, 2))
        vmax = float(np.nanpercentile(finite_values, 98))
        if np.isclose(vmin, vmax):
            vmin = float(np.nanmin(finite_values))
            vmax = float(np.nanmax(finite_values))
        if np.isclose(vmin, vmax):
            vmax = vmin + 1.0
    else:
        vmin, vmax = 0.0, 1.0
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)

    longest_category = max((len(category) for category in categories), default=0)
    category_rotation = 55 if longest_category > 24 else 35
    width = max(8.0, 0.95 * n_categories + 3.8)
    height = max(5.2, 0.38 * n_genes + 2.7, 4.2 + 0.025 * longest_category)
    fig, ax = plt.subplots(figsize=(width, height), dpi=150)
    x, y = np.meshgrid(np.arange(n_categories), np.arange(n_genes))
    values = np.nan_to_num(mean_expression.T, nan=vmin).ravel()
    fractions = np.clip(np.nan_to_num(detection_fraction.T, nan=0.0), 0.0, 1.0).ravel()
    sizes = 12.0 + 155.0 * fractions
    artist = ax.scatter(
        x.ravel(), y.ravel(), s=sizes, c=values, cmap=nature_continuous_cmap(),
        norm=norm, edgecolors='white', linewidths=0.35, alpha=0.95,
    )
    ax.set_xticks(np.arange(n_categories), categories, rotation=category_rotation, ha='right')
    section_values = []
    if gene_sections is not None:
        if isinstance(gene_sections, dict):
            section_values = [str(gene_sections.get(label, '')) for label in gene_labels]
        else:
            section_values = [str(item or '') for item in gene_sections]
        if len(section_values) != n_genes:
            raise ValueError('gene_sections must match gene_labels when provided')
    display_gene_labels = list(gene_labels)
    if section_values:
        previous = None
        for index, section in enumerate(section_values):
            if section and section != previous:
                display_gene_labels[index] = f'{section}\n{gene_labels[index]}'
                if index:
                    ax.axhline(index - 0.5, color='#98A2B3', linewidth=0.85, zorder=0)
            previous = section
    ax.set_yticks(np.arange(n_genes), display_gene_labels)
    ax.set_xlim(-0.55, n_categories - 0.45)
    ax.set_ylim(-0.55, n_genes - 0.45)
    ax.invert_yaxis()
    ax.set_title(title, loc='left', pad=10, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    ax.set_axisbelow(True)
    ax.grid(color=NATURE_GRID, linewidth=0.55, alpha=0.72)
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name in ('left', 'bottom'))
        spine.set_color('#98A2B3')
        spine.set_linewidth(0.7)

    colorbar = fig.colorbar(artist, ax=ax, fraction=0.035, pad=0.025, aspect=30)
    colorbar.outline.set_visible(False)
    colorbar.set_label(value_label, fontsize=8, labelpad=5)
    colorbar.ax.tick_params(labelsize=7, width=0.6, colors=NATURE_AXIS)

    size_values = (0.10, 0.50, 0.90)
    size_handles = [
        Line2D([], [], linestyle='', marker='o', markerfacecolor=NATURE_PALETTE[1],
               markeredgecolor='white', markersize=np.sqrt(12.0 + 155.0 * value) / 1.5,
               label=f'{int(value * 100)}%')
        for value in size_values
    ]
    fig.legend(
        handles=size_handles, title=detection_label, loc='upper right',
        bbox_to_anchor=(0.995, 0.92), frameon=False, fontsize=8,
        title_fontsize=8, borderaxespad=0.0,
    )
    fig.tight_layout(pad=1.1, rect=(0.0, 0.0, 0.78, 0.96))
    return fig


def annotation_strip_figure(values, title='', color_map=None, orientation='vertical'):
    """Draw a compact categorical annotation strip for a heatmap."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    values = [str(value) for value in values]
    unique = list(dict.fromkeys(values))
    color_map = color_map or {
        value: NATURE_PALETTE[index % len(NATURE_PALETTE)]
        for index, value in enumerate(unique)
    }
    indices = [unique.index(value) for value in values]
    if orientation == 'horizontal':
        matrix = np.asarray(indices, dtype=float)[None, :]
        figsize = (max(7.0, 0.22 * len(values) + 3), 1.8)
    else:
        matrix = np.asarray(indices, dtype=float)[:, None]
        figsize = (4.0, max(5.0, 0.20 * len(values) + 2.5))
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    ax.imshow(matrix, aspect='auto', cmap=ListedColormap([color_map[value] for value in unique]),
              vmin=-0.5, vmax=max(0.5, len(unique) - 0.5))
    if orientation == 'horizontal':
        ax.set_yticks([0], [title])
        ax.set_xticks([])
    else:
        ax.set_xticks([0], [title])
        ax.set_yticks([])
    handles = [Patch(facecolor=color_map[value], edgecolor='none', label=value)
               for value in unique]
    if handles:
        ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5),
                  frameon=False, fontsize=7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout(pad=0.8)
    return fig
