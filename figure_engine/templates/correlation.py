"""NatureCorrelationHeatmap：样本相关性的固定、诚实显示模板。"""

import numpy as np

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import attach_contract, label_tick_indices


def _cluster_order(matrix):
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    if matrix.shape[0] <= 2:
        return np.arange(matrix.shape[0])
    distance = np.clip(1.0 - matrix, 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    return leaves_list(linkage(squareform(distance, checks=False), method='average'))


class NatureCorrelationHeatmap:
    plot_type = 'correlation_heatmap'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

        style = get_style(spec.style)
        profile = style.profile(spec)
        if isinstance(data, dict):
            labels = [str(value) for value in data.get('sample_labels', ())]
            groups = data.get('groups')
            if data.get('correlation_matrix') is not None:
                matrix = np.asarray(data['correlation_matrix'], dtype=float)
            else:
                values = np.asarray(data.get('matrix'), dtype=float)
                if values.ndim != 2:
                    raise ValueError('NatureCorrelationHeatmap 需要二维 matrix')
                matrix = np.corrcoef(
                    values.T if data.get('samples_axis', 1) == 1 else values
                )
        else:
            matrix = np.asarray(data, dtype=float)
            labels, groups = [], None
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError('NatureCorrelationHeatmap 需要方阵 correlation_matrix')
        n = matrix.shape[0]
        if not labels:
            labels = [f'S{index + 1}' for index in range(n)]
        if len(labels) != n:
            raise ValueError('sample_labels 数量与相关矩阵不一致')
        order = _cluster_order(matrix) if spec.col_cluster else np.arange(n)
        ordered = matrix[np.ix_(order, order)]
        labels = [labels[index] for index in order]
        groups = [str(groups[index]) for index in order] if groups is not None else None
        off_diag = ordered[~np.eye(n, dtype=bool)] if n > 1 else ordered.ravel()
        finite = off_diag[np.isfinite(off_diag)]
        includes_negative = bool(finite.size and finite.min() < 0)
        if includes_negative:
            cmap = style.heatmap_cmap()
            limit = max(0.25, float(np.nanpercentile(np.abs(finite), 99)))
            norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
            image_kwargs = {'cmap': cmap, 'norm': norm}
        else:
            low = max(0.0, float(np.nanpercentile(finite, 2)) - 0.03) if finite.size else 0.0
            cmap = LinearSegmentedColormap.from_list(
                'nature_correlation', ['#F7FAFC', '#DCEAF0', '#8FA9C9', style.signal_blue])
            cmap.set_bad('#F2F3F5')
            image_kwargs = {'cmap': cmap, 'vmin': low, 'vmax': 1.0}
        display = ordered.copy()
        if spec.mask_diagonal:
            np.fill_diagonal(display, np.nan)

        with style.context(spec):
            if container is None:
                fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
            else:
                fig = container
                ax = container.subplots()
            image = ax.imshow(display, aspect='equal', interpolation='nearest',
                              rasterized=True, **image_kwargs)
            # Long Bulk sample IDs need substantially more horizontal room
            # than compact gene symbols at final submission size.
            physical_cap = 5 if spec.width == 'single' else 8
            maximum = min(spec.max_col_labels or physical_cap, physical_cap)
            ticks = label_tick_indices(n, maximum)
            ax.set_xticks(ticks, [labels[index] for index in ticks], rotation=55,
                          ha='right', rotation_mode='anchor')
            ax.set_yticks(ticks, [labels[index] for index in ticks])
            ax.tick_params(length=0, labelsize=max(5.5, profile.tick_font_pt - 0.3))
            ax.set_xlabel('Sample')
            ax.set_ylabel('Sample')
            # Reserve a figure-level title band above the annotation strip;
            # putting an Axes title at the top edge lets the strip visually
            # cut through long titles in the final export.
            title = spec.title or 'Sample correlation'
            fig.text(0.19, 0.965, title, ha='left', va='top',
                     fontsize=profile.title_font_pt, fontweight='semibold',
                     color=style.text)
            for spine in ax.spines.values():
                spine.set_visible(False)
            if spec.annotate_cells and n <= 12:
                for row in range(n):
                    for column in range(n):
                        if row == column and spec.mask_diagonal:
                            continue
                        value = ordered[row, column]
                        ax.text(column, row, f'{value:.2f}', ha='center', va='center',
                                fontsize=max(5.0, profile.tick_font_pt - 1.0),
                                color=style.text)
            colorbar = fig.colorbar(image, ax=ax, fraction=0.038, pad=0.025, aspect=28)
            colorbar.outline.set_visible(False)
            colorbar.set_label(f'{spec.correlation_method.title()} r',
                               fontsize=profile.legend_font_pt)
            colorbar.ax.tick_params(labelsize=max(5.3, profile.tick_font_pt - 0.6),
                                    length=1.8, width=profile.tick_width_pt)
            if groups:
                from matplotlib.patches import Rectangle
                color_map = style.group_colors(groups)
                for index, group in enumerate(groups):
                    ax.add_patch(Rectangle(
                        (index - 0.5, 1.035), 1.0, 0.025,
                        transform=ax.get_xaxis_transform(), clip_on=False,
                        facecolor=color_map[group], edgecolor='none', zorder=5,
                    ))
            fig.subplots_adjust(left=0.19, right=0.88, bottom=0.25, top=0.88)

        metadata = {
            'sample_order': [int(value) for value in order],
            'ordered_labels': labels,
            'sample_label_stride': max(1, int(np.ceil(n / max(1, len(ticks))))),
            'sample_labels_shown': int(len(ticks)),
            'ordered_groups': groups,
            'display_vmin': float(-limit if includes_negative else low),
            'display_vmax': float(limit if includes_negative else 1.0),
            'off_diagonal': {
                'min': float(np.nanmin(finite)) if finite.size else np.nan,
                'max': float(np.nanmax(finite)) if finite.size else np.nan,
                'median': float(np.nanmedian(finite)) if finite.size else np.nan,
            },
        }
        if container is not None:
            container._nature_correlation_metadata = metadata
            return container
        fig = attach_contract(
            fig, spec, style,
            encodings={'color': 'correlation', 'annotation': 'sample group' if groups else ''},
        )
        fig._nature_correlation_metadata = metadata
        return fig
