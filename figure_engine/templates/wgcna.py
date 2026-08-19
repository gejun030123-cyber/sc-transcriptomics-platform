"""NatureWGCNA：module–trait correlation heatmap."""

import numpy as np

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import attach_contract, label_tick_indices


class NatureWGCNA:
    plot_type = 'wgcna'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt

        style = get_style(spec.style)
        profile = style.profile(spec)
        payload = dict(data or {})
        matrix = np.asarray(payload.get('correlation_matrix', payload.get('matrix')), dtype=float)
        if matrix.ndim != 2:
            raise ValueError('NatureWGCNA 需要 module x trait correlation_matrix')
        modules = [str(value) for value in payload.get(
            'module_labels', [f'Module {index + 1}' for index in range(matrix.shape[0])])]
        traits = [str(value) for value in payload.get(
            'trait_labels', [f'Trait {index + 1}' for index in range(matrix.shape[1])])]
        pvalues = payload.get('pvalues')
        pvalues = np.asarray(pvalues, dtype=float) if pvalues is not None else None
        if len(modules) != matrix.shape[0] or len(traits) != matrix.shape[1]:
            raise ValueError('WGCNA module/trait 标签与矩阵维度不一致')
        if pvalues is not None and pvalues.shape != matrix.shape:
            raise ValueError('WGCNA pvalues 与 correlation_matrix 维度不一致')
        limit = max(0.25, min(1.0, float(np.nanmax(np.abs(matrix)))))

        with style.context(spec):
            if container is None:
                fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
            else:
                fig = container
                ax = container.subplots()
            image = ax.imshow(matrix, aspect='auto', interpolation='nearest',
                              cmap=style.heatmap_cmap(), vmin=-limit, vmax=limit,
                              rasterized=True)
            row_ticks = label_tick_indices(len(modules), spec.max_row_labels)
            col_ticks = label_tick_indices(len(traits), spec.max_col_labels)
            ax.set_yticks(row_ticks, [modules[index] for index in row_ticks])
            ax.set_xticks(col_ticks, [traits[index] for index in col_ticks],
                          rotation=45, ha='right', rotation_mode='anchor')
            ax.tick_params(length=0, labelsize=max(5.5, profile.tick_font_pt - 0.3))
            ax.set_xlabel('Trait')
            ax.set_ylabel('Co-expression module')
            ax.set_title(spec.title or 'Module–trait relationships', loc='left', pad=5)
            for spine in ax.spines.values():
                spine.set_visible(False)
            annotate = spec.annotate_cells or matrix.size <= 80
            if annotate:
                for row in range(matrix.shape[0]):
                    for column in range(matrix.shape[1]):
                        value = matrix[row, column]
                        suffix = ''
                        if pvalues is not None and np.isfinite(pvalues[row, column]):
                            pvalue = pvalues[row, column]
                            suffix = '***' if pvalue < 0.001 else ('**' if pvalue < 0.01 else ('*' if pvalue < 0.05 else ''))
                        ax.text(column, row, f'{value:.2f}{suffix}', ha='center', va='center',
                                fontsize=max(5.0, profile.tick_font_pt - 0.8),
                                color='white' if abs(value) > limit * 0.62 else style.text)
            colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025, aspect=30)
            colorbar.outline.set_visible(False)
            colorbar.set_label('Module–trait correlation', fontsize=profile.legend_font_pt)
            colorbar.ax.tick_params(labelsize=max(5.2, profile.tick_font_pt - 0.7),
                                    length=1.8, width=profile.tick_width_pt)
            fig.subplots_adjust(left=0.25 if spec.width == 'single' else 0.18,
                                right=0.88, bottom=0.23, top=0.88)

        if container is not None:
            return container
        return attach_contract(
            fig, spec, style,
            encodings={'color': 'correlation', 'text': 'correlation and significance'},
        )
