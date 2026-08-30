"""NatureHeatmap：基因在行、样本在列并集成注释与聚类。"""

import numpy as np

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import attach_contract, label_tick_indices, robust_symmetric_limit


def _zscore_rows(matrix):
    matrix = np.asarray(matrix, dtype=float)
    means = np.nanmean(matrix, axis=1, keepdims=True)
    scales = np.nanstd(matrix, axis=1, ddof=0, keepdims=True)
    scales[~np.isfinite(scales) | (scales == 0)] = 1.0
    result = (matrix - means) / scales
    return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)


def _safe_matrix(matrix):
    result = np.asarray(matrix, dtype=float).copy()
    if not np.isnan(result).any() and np.isfinite(result).all():
        return result
    row_means = np.nanmean(np.where(np.isfinite(result), result, np.nan), axis=1)
    row_means[~np.isfinite(row_means)] = 0.0
    invalid = ~np.isfinite(result)
    result[invalid] = np.take(row_means, np.where(invalid)[0])
    return result


def _linkage_and_order(matrix, axis, method, metric):
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist

    values = matrix if axis == 0 else matrix.T
    if values.shape[0] <= 2:
        return None, np.arange(values.shape[0], dtype=int)
    safe_values = _safe_matrix(values)
    use_metric = 'euclidean' if method == 'ward' else str(metric).lower()
    if use_metric == 'pearson':
        use_metric = 'correlation'
    elif use_metric == 'spearman':
        from scipy.stats import rankdata

        safe_values = np.apply_along_axis(rankdata, 1, safe_values)
        use_metric = 'correlation'
    distances = pdist(safe_values, metric=use_metric)
    if not np.isfinite(distances).all() or not np.any(distances > 0):
        return None, np.arange(values.shape[0], dtype=int)
    tree = linkage(distances, method=method, optimal_ordering=True)
    return tree, leaves_list(tree).astype(int)


class NatureHeatmap:
    """固定的 transcriptomics publication heatmap 模板。"""

    plot_type = 'heatmap'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt
        from matplotlib.colors import to_rgba
        from matplotlib.patches import Patch
        from scipy.cluster.hierarchy import dendrogram

        style = get_style(spec.style)
        profile = style.profile(spec)
        matrix = np.asarray(data.get('matrix'), dtype=float)
        if matrix.ndim != 2 or min(matrix.shape) == 0:
            raise ValueError('NatureHeatmap 需要非空的 genes x samples 二维 matrix')
        gene_labels = [str(value) for value in data.get('gene_labels', [])]
        sample_labels = [str(value) for value in data.get('sample_labels', [])]
        if len(gene_labels) != matrix.shape[0] or len(sample_labels) != matrix.shape[1]:
            raise ValueError('Heatmap gene/sample labels 与 matrix 维度不一致')
        annotations = {
            str(name): [str(value) for value in values]
            for name, values in dict(data.get('annotations') or {}).items()
        }
        for name, values in annotations.items():
            if len(values) != matrix.shape[1]:
                raise ValueError(f'Heatmap annotation {name} 与样本数不一致')

        colorbar_label = str(data.get('colorbar_label') or (
            'Gene-wise z-score' if spec.zscore == 'row' else 'Expression'
        ))

        display = _zscore_rows(matrix) if spec.zscore == 'row' else _safe_matrix(matrix)
        row_tree, row_order = (
            _linkage_and_order(display, 0, spec.cluster_method, spec.distance_metric)
            if spec.row_cluster else (None, np.arange(display.shape[0]))
        )
        col_tree, col_order = (
            _linkage_and_order(display, 1, spec.cluster_method, spec.distance_metric)
            if spec.col_cluster else (None, np.arange(display.shape[1]))
        )
        display = display[np.ix_(row_order, col_order)]
        ordered_genes = [gene_labels[index] for index in row_order]
        ordered_samples = [sample_labels[index] for index in col_order]
        ordered_annotations = {
            name: [values[index] for index in col_order]
            for name, values in annotations.items()
        }
        color_limit = float(spec.color_limit or robust_symmetric_limit(display))
        symmetric_color = bool(data.get('symmetric_color', True))
        explicit_vmin = data.get('vmin')
        explicit_vmax = data.get('vmax')
        warnings = []
        if matrix.shape[0] > 50 and spec.width == 'single':
            warnings.append(
                f'单栏热图包含 {matrix.shape[0]} 个基因；建议主文限制为 20-30 个或改用双栏。'
            )
        if matrix.shape[1] > 60 and spec.width == 'single':
            warnings.append(f'单栏热图包含 {matrix.shape[1]} 个样本，样本标签将按密度抽样。')

        with style.context(spec):
            fig = (plt.figure(figsize=profile.figsize, dpi=profile.dpi)
                   if container is None else container)
            has_row_tree = row_tree is not None
            has_col_tree = col_tree is not None
            n_annotations = len(ordered_annotations)
            row_heights = [0.13 if has_col_tree else 0.005,
                           max(0.005, 0.052 * n_annotations), 0.865]
            col_widths = [0.13 if has_row_tree else 0.005, 0.87]
            right = 0.76 if spec.width == 'single' else 0.84
            grid = fig.add_gridspec(
                3, 2, height_ratios=row_heights, width_ratios=col_widths,
                left=0.08, right=right, bottom=0.35, top=0.89,
                wspace=0.025, hspace=0.025,
            )
            ax_col = fig.add_subplot(grid[0, 1])
            ax_annotation = fig.add_subplot(grid[1, 1])
            ax_row = fig.add_subplot(grid[2, 0])
            ax_heat = fig.add_subplot(grid[2, 1])

            if has_col_tree:
                dendrogram(
                    col_tree, ax=ax_col, no_labels=True,
                    color_threshold=0, above_threshold_color=style.neutral_dark,
                    link_color_func=lambda _: style.neutral_dark,
                )
                for collection in ax_col.collections:
                    collection.set_linewidth(0.45)
                ax_col.set_axis_off()
            else:
                ax_col.set_axis_off()

            if has_row_tree:
                dendrogram(
                    row_tree, ax=ax_row, orientation='left', no_labels=True,
                    color_threshold=0, above_threshold_color=style.neutral_dark,
                    link_color_func=lambda _: style.neutral_dark,
                )
                for collection in ax_row.collections:
                    collection.set_linewidth(0.45)
                ax_row.invert_yaxis()
                ax_row.set_axis_off()
            else:
                ax_row.set_axis_off()

            annotation_handles = []
            if n_annotations:
                rgba = np.ones((n_annotations, display.shape[1], 4), dtype=float)
                palette_offset = 0
                for row_index, (name, values) in enumerate(ordered_annotations.items()):
                    unique = list(dict.fromkeys(values))
                    color_map = {
                        value: style.categorical_palette[(palette_offset + index)
                                                         % len(style.categorical_palette)]
                        for index, value in enumerate(unique)
                    }
                    palette_offset += max(1, len(unique))
                    rgba[row_index] = np.asarray([to_rgba(color_map[value]) for value in values])
                    annotation_handles.extend([
                        Patch(facecolor=color_map[value], edgecolor='none', label=f'{name}: {value}')
                        for value in unique
                    ])
                ax_annotation.imshow(rgba, aspect='auto', interpolation='nearest')
                ax_annotation.set_yticks(
                    np.arange(n_annotations), list(ordered_annotations),
                    fontsize=max(5.5, profile.tick_font_pt - 0.4),
                )
                ax_annotation.set_xticks([])
                ax_annotation.tick_params(axis='y', length=0, pad=2, colors=style.axis)
                for spine in ax_annotation.spines.values():
                    spine.set_visible(False)
            else:
                ax_annotation.set_axis_off()

            image_kwargs = {
                'cmap': style.heatmap_cmap(),
                'vmin': (-color_limit if symmetric_color else (0.0 if explicit_vmin is None else float(explicit_vmin))),
                'vmax': (color_limit if symmetric_color else (color_limit if explicit_vmax is None else float(explicit_vmax))),
                'rasterized': True,
            }
            image = ax_heat.imshow(
                display, aspect='auto', interpolation='nearest', **image_kwargs,
            )
            # The final-size row capacity is a physical typography constraint,
            # not an arbitrary UI preference.  At 89 mm, more than 21 labels
            # in this layout touch even at the allowed minimum font size.
            physical_row_cap = 21 if spec.width == 'single' else 45
            row_ticks = label_tick_indices(
                len(ordered_genes), min(spec.max_row_labels, physical_row_cap)
            )
            col_ticks = label_tick_indices(len(ordered_samples), spec.max_col_labels)
            ax_heat.set_yticks(row_ticks, [ordered_genes[index] for index in row_ticks])
            ax_heat.yaxis.tick_right()
            ax_heat.tick_params(axis='y', labelsize=profile.tick_font_pt, length=0, pad=2)
            ax_heat.set_xticks(col_ticks, [ordered_samples[index] for index in col_ticks],
                               rotation=90, ha='center', va='top')
            ax_heat.tick_params(axis='x', labelsize=max(5.5, profile.tick_font_pt - 0.4),
                                length=0, pad=2)
            for spine in ax_heat.spines.values():
                spine.set_visible(False)
            ax_heat.set_xlabel('')
            ax_heat.set_ylabel('')

            title = fig.text(
                0.08, 0.96, spec.title or 'Expression heatmap',
                ha='left', va='top', fontsize=profile.title_font_pt,
                fontweight='semibold', color=style.text,
            )
            cbar_left = 0.31 if spec.width == 'single' else 0.38
            cbar_width = 0.30 if spec.width == 'single' else 0.25
            cbar_ax = fig.add_axes([cbar_left, 0.145, cbar_width, 0.018])
            colorbar = fig.colorbar(image, cax=cbar_ax, orientation='horizontal')
            colorbar.outline.set_visible(False)
            colorbar.ax.tick_params(
                labelsize=max(5.5, profile.tick_font_pt - 0.5),
                width=profile.tick_width_pt, length=1.8, colors=style.axis,
            )
            colorbar.set_label(
                colorbar_label,
                fontsize=profile.legend_font_pt, color=style.text, labelpad=2,
            )
            if annotation_handles:
                if len(annotation_handles) > 10:
                    warnings.append(
                        f'Heatmap annotation legend 包含 {len(annotation_handles)} 个类别；已省略图例以避免单栏拥挤，建议拆图或合并低频类别。'
                    )
                else:
                    fig.legend(
                        handles=annotation_handles, loc='lower center',
                        bbox_to_anchor=(0.5, 0.005), ncol=min(4 if spec.width == 'single' else 6,
                                                            len(annotation_handles)),
                        frameon=False, handlelength=0.8, handleheight=0.65,
                        handletextpad=0.3, columnspacing=0.65,
                        fontsize=max(5.3, profile.legend_font_pt - 0.5),
                    )

        if container is not None:
            return container
        return attach_contract(
            fig, spec, style,
            semantic_warnings=warnings,
            encodings={
                'color': colorbar_label.lower(),
                'annotation': ', '.join(ordered_annotations),
            },
            validation_texts=[title],
        )
