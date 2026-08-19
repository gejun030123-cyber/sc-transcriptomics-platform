"""NaturePCA：颜色表达 group，marker 表达 batch。"""

import numpy as np

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import adjust_labels, as_1d, attach_contract


def _confidence_ellipse(ax, x, y, color, level=0.95):
    """绘制轻量置信椭圆；少于 3 个点时不绘制。"""
    from matplotlib.patches import Ellipse
    from scipy.stats import chi2

    points = np.column_stack([x, y]).astype(float)
    points = points[np.isfinite(points).all(axis=1)]
    if len(points) < 3:
        return None
    covariance = np.cov(points, rowvar=False)
    if covariance.shape != (2, 2) or not np.isfinite(covariance).all():
        return None
    values, vectors = np.linalg.eigh(covariance)
    if np.any(values <= 0):
        return None
    order = values.argsort()[::-1]
    values, vectors = values[order], vectors[:, order]
    angle = float(np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0])))
    scale = float(np.sqrt(chi2.ppf(level, df=2)))
    ellipse = Ellipse(
        xy=np.nanmean(points, axis=0),
        width=2 * scale * np.sqrt(values[0]),
        height=2 * scale * np.sqrt(values[1]),
        angle=angle,
        facecolor=color,
        edgecolor=color,
        linewidth=0.65,
        alpha=0.10,
        zorder=1,
    )
    ax.add_patch(ellipse)
    return ellipse


class NaturePCA:
    """固定的 PCA publication 模板。"""

    plot_type = 'pca'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        style = get_style(spec.style)
        profile = style.profile(spec)
        coordinates = np.asarray(data.get('coordinates'), dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError('NaturePCA 需要 n_samples x 2 的 coordinates')
        coordinates = coordinates[:, :2]
        n_samples = coordinates.shape[0]
        groups = as_1d(data.get('groups'), n_samples, 'All')
        batches = as_1d(data.get('batches'), n_samples, 'All')
        samples = as_1d(data.get('samples'), n_samples, '')
        variance = np.asarray(data.get('explained_variance', [np.nan, np.nan]), dtype=float)
        variance = variance[:2]
        finite_variance = np.abs(variance[np.isfinite(variance)])
        if finite_variance.size and np.nanmax(finite_variance) <= 1.0:
            variance = variance * 100.0
        group_labels = list(dict.fromkeys(str(value) for value in groups))
        batch_labels = list(dict.fromkeys(str(value) for value in batches))
        color_map = dict(data.get('group_colors') or style.group_colors(group_labels))
        marker_map = dict(data.get('batch_markers') or style.batch_markers(batch_labels))

        warnings = []
        if len(group_labels) > 8:
            warnings.append(
                f'PCA 包含 {len(group_labels)} 个颜色组；建议将复合分组拆为 color + marker 或 facet。'
            )
        if len(batch_labels) > len(style.marker_cycle):
            warnings.append(f'PCA batch 数量 {len(batch_labels)} 超过可稳定区分的 marker 数量。')

        with style.context(spec):
            if container is None:
                fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
            else:
                fig = container
                ax = container.subplots()
            for group in group_labels:
                group_mask = groups.astype(str) == group
                if spec.confidence_ellipse:
                    _confidence_ellipse(
                        ax, coordinates[group_mask, 0], coordinates[group_mask, 1],
                        color_map[group], spec.ellipse_level,
                    )
                for batch in batch_labels:
                    mask = group_mask & (batches.astype(str) == batch)
                    if not mask.any():
                        continue
                    ax.scatter(
                        coordinates[mask, 0], coordinates[mask, 1],
                        s=profile.marker_size_pt2,
                        marker=marker_map[batch],
                        color=color_map[group],
                        edgecolor='white', linewidth=0.45,
                        alpha=0.92, zorder=3,
                    )

            label_mask = np.zeros(n_samples, dtype=bool)
            if spec.show_sample_labels:
                label_mask[:] = True
                if n_samples > 16:
                    warnings.append('样本数超过 16，全部 sample label 可能降低 PCA 可读性。')
            if spec.outlier_labels:
                label_mask |= np.isin(samples.astype(str), list(spec.outlier_labels))
            texts = []
            x_span = max(float(np.ptp(coordinates[:, 0])), 1.0)
            y_span = max(float(np.ptp(coordinates[:, 1])), 1.0)
            for index in np.flatnonzero(label_mask):
                label = str(samples[index])
                if not label:
                    continue
                texts.append(ax.text(
                    coordinates[index, 0] + x_span * 0.012,
                    coordinates[index, 1] + y_span * 0.012,
                    label, fontsize=profile.tick_font_pt,
                    color=style.text, ha='left', va='bottom', zorder=5,
                ))
            adjust_labels(texts, ax, arrow_color=style.neutral_dark)

            if np.nanmin(coordinates[:, 0]) < 0 < np.nanmax(coordinates[:, 0]):
                ax.axvline(0, color=style.subtle, linewidth=0.45, zorder=0)
            if np.nanmin(coordinates[:, 1]) < 0 < np.nanmax(coordinates[:, 1]):
                ax.axhline(0, color=style.subtle, linewidth=0.45, zorder=0)
            x_suffix = f' ({variance[0]:.1f}%)' if variance.size > 0 and np.isfinite(variance[0]) else ''
            y_suffix = f' ({variance[1]:.1f}%)' if variance.size > 1 and np.isfinite(variance[1]) else ''
            x_label = str(data.get('x_label') or 'PC1')
            y_label = str(data.get('y_label') or 'PC2')
            ax.set_xlabel(f'{x_label}{x_suffix}')
            ax.set_ylabel(f'{y_label}{y_suffix}')
            ax.set_title(spec.title or 'Principal component analysis', loc='left', pad=5)
            ax.margins(0.10)
            style.apply_axis(ax, profile)

            handles = [
                Line2D([], [], linestyle='None', marker='o', markersize=4.2,
                       markerfacecolor=color_map[group], markeredgecolor='white',
                       markeredgewidth=0.4,
                       label=f'Group · {group}' if len(batch_labels) > 1 else group)
                for group in group_labels
            ]
            if len(batch_labels) > 1:
                handles.extend([
                    Line2D([], [], linestyle='None', marker=marker_map[batch], markersize=4.2,
                           markerfacecolor=style.neutral_dark, markeredgecolor='white',
                           markeredgewidth=0.4, label=f'Batch · {batch}')
                    for batch in batch_labels
                ])
            if handles and spec.show_legend:
                columns = min(4 if spec.width == 'double' else 3, max(1, len(handles)))
                rows = int(np.ceil(len(handles) / columns))
                bottom = min(0.42, 0.16 + rows * 0.055)
                fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=bottom)
                fig.legend(
                    handles=handles, loc='lower center', ncol=columns,
                    bbox_to_anchor=(0.5, 0.015), frameon=False,
                    handletextpad=0.35, columnspacing=0.8,
                )
            else:
                fig.subplots_adjust(left=0.18, right=0.96, top=0.90, bottom=0.17)

        if container is not None:
            return container
        return attach_contract(
            fig, spec, style,
            semantic_warnings=warnings,
            encodings={'color': 'group', 'marker': 'batch' if len(batch_labels) > 1 else ''},
            validation_texts=texts,
        )
