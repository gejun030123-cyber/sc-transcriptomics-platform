"""NatureEmbedding: publication-safe single-cell UMAP/t-SNE/PCA point clouds."""

from __future__ import annotations

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import adjust_labels, attach_contract, prune_overlapping_labels


def _marker_size(n_points, width):
    """Derive point size from density and final physical width."""
    base = 7.5 if width == 'single' else 9.0
    scale = np.sqrt(2500.0 / max(1, int(n_points)))
    return float(np.clip(base * scale, 0.45, 10.0))


def _continuous_limits(values):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0, 1.0, False
    lower, upper = np.nanpercentile(finite, [1.0, 99.0])
    central_span = max(float(upper - lower), np.finfo(float).eps)
    clipped = bool(
        np.any(finite < lower - central_span * 0.10)
        or np.any(finite > upper + central_span * 0.10)
    )
    if not np.isfinite(lower) or not np.isfinite(upper) or np.isclose(lower, upper):
        lower, upper = float(np.nanmin(finite)), float(np.nanmax(finite))
    if np.isclose(lower, upper):
        padding = max(abs(float(lower)) * 0.05, 0.5)
        lower, upper = float(lower - padding), float(upper + padding)
    return float(lower), float(upper), clipped


class NatureEmbedding:
    """Render dense cell embeddings without applying sample-level PCA semantics."""

    plot_type = 'embedding'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        style = get_style(spec.style)
        profile = style.profile(spec)
        coordinates = np.asarray(data.get('coordinates'), dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError('NatureEmbedding 需要 n_cells x 2 的 coordinates')
        coordinates = coordinates[:, :2]
        finite_coordinates = np.isfinite(coordinates).all(axis=1)
        if not finite_coordinates.any():
            raise ValueError('NatureEmbedding 没有有限坐标')

        values = data.get('values')
        if values is None:
            values = np.repeat('All cells', len(coordinates))
        values = np.asarray(values).reshape(-1)
        if len(values) != len(coordinates):
            raise ValueError('NatureEmbedding values 与 coordinates 长度不一致')
        coordinates = coordinates[finite_coordinates]
        values = values[finite_coordinates]

        requested_type = str(data.get('value_type') or 'auto').lower()
        numeric_values = pd.to_numeric(pd.Series(values), errors='coerce').to_numpy(dtype=float)
        continuous = requested_type == 'continuous' or (
            requested_type == 'auto' and pd.api.types.is_numeric_dtype(values)
        )
        value_label = str(data.get('value_label') or 'Cells')
        x_label = str(data.get('x_label') or 'UMAP 1')
        y_label = str(data.get('y_label') or 'UMAP 2')
        hide_axes = bool(data.get('hide_axes', True))
        direct_labels = bool(data.get('direct_labels', False))
        point_size = _marker_size(len(coordinates), spec.width)
        warnings = []
        texts = []
        use_direct_labels = False

        with style.context(spec):
            if container is None:
                fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
            else:
                fig = container
                ax = container.subplots()

            if continuous:
                finite_values = np.isfinite(numeric_values)
                if not finite_values.any():
                    raise ValueError('NatureEmbedding 连续颜色变量没有有限数值')
                if (~finite_values).any():
                    ax.scatter(
                        coordinates[~finite_values, 0], coordinates[~finite_values, 1],
                        s=point_size, color=style.neutral, alpha=0.42,
                        linewidths=0, rasterized=True, zorder=1,
                    )
                lower, upper, clipped = _continuous_limits(numeric_values[finite_values])
                order = np.argsort(numeric_values[finite_values], kind='stable')
                numeric_coordinates = coordinates[finite_values][order]
                shown_values = numeric_values[finite_values][order]
                artist = ax.scatter(
                    numeric_coordinates[:, 0], numeric_coordinates[:, 1],
                    s=point_size, c=shown_values, cmap=style.expression_cmap(),
                    vmin=lower, vmax=upper, alpha=0.78, linewidths=0,
                    rasterized=True, zorder=2,
                )
                colorbar = fig.colorbar(artist, ax=ax, fraction=0.035, pad=0.025, aspect=28)
                colorbar.outline.set_visible(False)
                colorbar.ax.tick_params(
                    labelsize=profile.tick_font_pt, width=profile.tick_width_pt,
                    length=2, colors=style.axis,
                )
                colorbar.set_label(value_label, fontsize=profile.axis_font_pt, color=style.text)
                if clipped:
                    warnings.append('连续颜色按第 1–99 百分位显示；完整数值保留在 AnnData。')
                encodings = {'color': f'continuous {value_label}', 'point': 'cell'}
                # Reserve a full right-side band for a vertical colorbar label;
                # long labels such as "log-normalized expression" otherwise
                # extend beyond an 89 mm single-column canvas.
                right = 0.84
            else:
                labels = np.asarray([str(value) for value in values], dtype=object)
                requested_order = [str(value) for value in data.get('category_order', [])]
                categories = requested_order or list(dict.fromkeys(labels.tolist()))
                categories = [value for value in categories if np.any(labels == value)]
                colors = style.group_colors(categories)
                counts = {category: int(np.sum(labels == category)) for category in categories}
                for category in categories:
                    mask = labels == category
                    ax.scatter(
                        coordinates[mask, 0], coordinates[mask, 1],
                        s=point_size, color=colors[category], alpha=0.74,
                        linewidths=0, rasterized=True, zorder=2,
                    )
                if len(categories) > len(style.categorical_palette):
                    warnings.append(
                        f'嵌入图包含 {len(categories)} 个类别，超过稳定色板容量；'
                        '建议按主要细胞谱系拆图或合并低频类别。'
                    )
                use_direct_labels = direct_labels or len(categories) > 12
                if use_direct_labels:
                    label_cap = min(20, len(categories))
                    labelled = sorted(categories, key=lambda value: (-counts[value], value))[:label_cap]
                    for category in labelled:
                        mask = labels == category
                        center = np.nanmedian(coordinates[mask], axis=0)
                        texts.append(ax.text(
                            center[0], center[1], category, ha='center', va='center',
                            fontsize=max(5.3, profile.tick_font_pt - 0.25), color=style.text,
                            bbox={'boxstyle': 'round,pad=0.14', 'facecolor': 'white',
                                  'edgecolor': style.subtle, 'alpha': 0.88, 'linewidth': 0.35},
                            zorder=5,
                        ))
                    if len(categories) > label_cap:
                        warnings.append(
                            f'仅直接标注细胞数最多的 {label_cap} 个类别；完整类别保留在 AnnData。'
                        )
                    right = 0.97
                else:
                    handles = [
                        Line2D([], [], linestyle='None', marker='o', markersize=3.8,
                               markerfacecolor=colors[category], markeredgewidth=0,
                               label=category)
                        for category in categories
                    ]
                    ax.legend(
                        handles=handles, title=value_label, loc='center left',
                        bbox_to_anchor=(1.01, 0.5), frameon=False,
                        fontsize=profile.legend_font_pt,
                        title_fontsize=profile.legend_font_pt,
                        borderaxespad=0, handletextpad=0.35,
                    )
                    right = 0.77
                encodings = {'color': f'categorical {value_label}', 'point': 'cell'}

            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(spec.title or 'Single-cell embedding', loc='left', pad=5)
            ax.margins(0.035)
            style.apply_axis(ax, profile, despine=True)
            ax.grid(False)
            for spine in ax.spines.values():
                spine.set_visible(False)
            if hide_axes:
                ax.set_xticks([])
                ax.set_yticks([])
            fig.subplots_adjust(left=0.08, right=right, bottom=0.10, top=0.90)

            if use_direct_labels:
                # Label positions must be adjusted after the final axes
                # rectangle is established.  The previous call happened
                # before ``subplots_adjust`` and therefore used a different
                # pixel geometry from the exported PNG/SVG; dense annotation
                # UMAPs could collide again at the final size.
                adjust_labels(texts, ax, arrow_color=style.neutral_dark)
                hidden_labels = prune_overlapping_labels(texts, ax)
                if hidden_labels:
                    warnings.append(
                        '为避免最终尺寸文字重叠，已隐藏 '
                        f'{len(hidden_labels)} 个低优先级细胞类型标签；完整类别保留在 AnnData。'
                    )

        if len(coordinates) > 100000:
            warnings.append('细胞点层已栅格化以控制 SVG/PDF 文件体积；文字和坐标轴保持矢量。')
        if container is not None:
            return container
        return attach_contract(
            fig, spec, style, semantic_warnings=warnings,
            encodings=encodings, validation_texts=texts,
        )
