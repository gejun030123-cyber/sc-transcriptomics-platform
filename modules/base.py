from abc import ABC, abstractmethod
from typing import Callable, Optional
from modules.figure_style import (
    NATURE_BG,
    NATURE_FONT_FAMILY,
    NATURE_PALETTE,
    nature_continuous_cmap,
)

# Result files are persisted in the database and exposed by the result routes.
# Keep the file-type contract in one place so modules, tests and API consumers
# agree on optional spreadsheet exports as well as image/table artifacts.
VALID_RESULT_FILE_TYPES = frozenset({
    'csv', 'xlsx', 'plotly_json', 'png', 'svg', 'jpg', 'jpeg',
    'info', 'json', 'txt', 'h5ad',
})

VISUALIZATION_THEMES = {
    'default': {
        'bg_color': NATURE_BG,
        'color_palette': list(NATURE_PALETTE),
        'font_family': NATURE_FONT_FAMILY,
    },
    'nature': {
        'bg_color': NATURE_BG,
        'color_palette': list(NATURE_PALETTE),
        'font_family': 'Helvetica',
    },
    'dark': {
        'bg_color': '#1a1a2e',
        'color_palette': ['#e94560', '#0f3460', '#16213e', '#533483',
                          '#00b4d8', '#48cae4', '#90e0ef', '#f8f9fa'],
        'font_family': 'Arial',
    },
}


class BaseAnalysis(ABC):
    MODULE_NAME = ""
    DISPLAY_NAME = ""
    DESCRIPTION = ""
    INPUT_REQUIRES = []
    PARAM_SCHEMA = {}

    def __init__(self, project_dir: str, params: dict,
                 progress_callback: Callable[[int, str], None]):
        self.project_dir = project_dir
        self.params = params
        self._progress = progress_callback

    def progress(self, pct: int, message: str):
        if not self._progress:
            return
        self._progress(min(pct, 100), message)

    def apply_filters(self, adata, module_name):
        """根据 self.params['_filters'][module_name] 中的声明式规则过滤 adata。"""
        import logging
        logger = logging.getLogger(__name__)
        filters = self.params.get('_filters', {}).get(module_name, [])
        for rule in filters:
            col = rule.get('column', '')
            op = rule.get('op', '')
            val = rule.get('value')
            if col not in adata.obs.columns:
                self.progress(-1, f"警告: 列 '{col}' 不存在于 obs 中，跳过该过滤规则")
                continue
            if op == '>=':
                mask = adata.obs[col] >= val
            elif op == '<=':
                mask = adata.obs[col] <= val
            elif op == '==':
                mask = adata.obs[col] == val
            elif op == '!=':
                mask = adata.obs[col] != val
            elif op == 'in':
                if not hasattr(val, '__iter__') or isinstance(val, str):
                    logger.warning(f"过滤规则 'in' 需要列表类型的 value，收到: {type(val).__name__}")
                    continue
                mask = adata.obs[col].isin(val)
            elif op == 'not_in':
                if not hasattr(val, '__iter__') or isinstance(val, str):
                    logger.warning(f"过滤规则 'not_in' 需要列表类型的 value，收到: {type(val).__name__}")
                    continue
                mask = ~adata.obs[col].isin(val)
            elif op == 'between':
                if not isinstance(val, (list, tuple)) or len(val) < 2:
                    logger.warning(f"过滤规则 'between' 需要 [min, max]，收到: {val}")
                    continue
                mask = adata.obs[col].between(val[0], val[1])
            else:
                logger.warning(f"不支持的过滤操作符: {op!r}，规则: {rule}")
                continue
            n_before = adata.n_obs
            adata = adata[mask].copy()
            n_filtered = n_before - adata.n_obs
            if n_filtered > 0:
                self.progress(-1, f"过滤 {col} {op} {val}: 移除 {n_filtered} 个样本")
        return adata

    def get_plotly_layout(self, title='', **overrides):
        """根据 self.params['_visualization'] 生成统一的 Plotly layout dict。"""
        viz = self.params.get('_visualization', {})
        theme_name = viz.get('theme', 'default')
        theme = VISUALIZATION_THEMES.get(theme_name, VISUALIZATION_THEMES['default'])
        layout = {
            'title': title,
            'width': viz.get('figure_width', 800),
            'height': viz.get('figure_height', 500),
            'plot_bgcolor': viz.get('bg_color', theme['bg_color']),
            'colorway': viz.get('color_palette', theme['color_palette']),
            'font': {
                'family': viz.get('font_family', theme['font_family']),
                'size': viz.get('font_size', 12),
            },
        }
        layout.update(overrides)
        return layout

    def get_viz_params(self):
        """返回传递给 umap_scatter 的可视化参数 dict。"""
        viz = self.params.get('_visualization', {})
        return {
            'umap_point_size': viz.get('umap_point_size', 5),
            'umap_opacity': viz.get('umap_opacity', 0.7),
            'umap_legend_fontsize': viz.get('umap_legend_fontsize', 10),
        }

    def export_static_fig(self, fig, plots_dir, filename, export_formats=None):
        """将 Plotly figure 导出为静态图片（SVG/PNG）。依赖 kaleido。"""
        if not export_formats:
            return []
        import os
        from modules.figure_style import style_plotly_figure
        style_plotly_figure(fig, self.params.get('_visualization', {}))
        exported = []
        for fmt in export_formats:
            if fmt in ('svg', 'png'):
                try:
                    import plotly.io as pio
                    out_path = os.path.join(plots_dir, filename.replace('.json', f'.{fmt}'))
                    pio.write_image(fig, out_path, format=fmt, engine='kaleido')
                    exported.append(out_path)
                except ImportError:
                    self.progress(-1, "kaleido 未安装，跳过静态图片导出")
                except Exception as e:
                    self.progress(-1, f"静态导出 {fmt} 失败: {e}")
        return exported

    def save_matplotlib_figure(self, fig, plots_dir, filename, category, label,
                               formats=None, dpi=300, preserve_aspect=False):
        """Persist a publication-quality Matplotlib/Scanpy figure and register it.

        Static figures are intentionally generated independently of Plotly/Kaleido:
        Scanpy and OmicVerse native plots can therefore retain their typography and
        vector geometry even when the optional Kaleido renderer is unavailable.
        """
        import os
        from modules.figure_style import apply_matplotlib_style, _font_for_text

        viz = self.params.get('_visualization', {})
        if formats is None:
            formats = viz.get(
                'static_formats', ('png', 'svg')
            )
        formats = [fmt.lower() for fmt in formats if fmt.lower() in ('png', 'svg')]
        if not formats:
            return []

        # Normalize the native Scanpy/Matplotlib canvas for publication output.
        # This changes the figure itself before both PNG and SVG are written, so
        # the browser preview and the vector download share the same composition.
        try:
            # Native builders choose their own aspect ratio (for example a
            # multi-resolution UMAP grid or a long cluster composition panel).
            # Callers can preserve that aspect explicitly; otherwise the
            # requested visualization dimensions are applied to the canvas.
            if ('figure_width' in viz or 'figure_height' in viz) and not preserve_aspect:
                current_width, current_height = fig.get_size_inches()
                width = max(6.5, float(viz.get('figure_width', current_width * 100)) / 100)
                height = max(4.8, float(viz.get('figure_height', current_height * 100)) / 100)
                fig.set_size_inches(width, height, forward=True)
            font_size = max(9, float(viz.get('font_size', 12)))
            font_family = viz.get('font_family', NATURE_FONT_FAMILY)
            apply_matplotlib_style(fig, viz)
            for axis in getattr(fig, 'axes', []):
                axis.set_facecolor(viz.get('bg_color', 'white'))
                axis.tick_params(labelsize=max(8, font_size - 2), width=0.7,
                                 colors='#374151')
                axis.xaxis.label.set_size(font_size)
                axis.yaxis.label.set_size(font_size)
                axis.xaxis.label.set_color('#1f2937')
                axis.yaxis.label.set_color('#1f2937')
                title = axis.title
                title.set_fontsize(font_size + 1)
                title.set_fontweight('semibold')
                title.set_color('#111827')
                title.set_fontfamily(_font_for_text(title.get_text(), font_family))
                for spine in axis.spines.values():
                    spine.set_linewidth(0.65)
                    spine.set_color('#c7cdd6')
                legend = axis.get_legend()
                if legend is not None:
                    legend.set_frame_on(False)
                    for text in legend.get_texts():
                        text.set_fontsize(max(8, font_size - 2))
                        text.set_fontfamily(_font_for_text(text.get_text(), font_family))
            # Figure-level legends (used by enrichment charts) need an explicit
            # right-side reservation; otherwise tight_layout expands the axes
            # underneath the legend and clips/overlaps the ontology key.
            if getattr(fig, 'legends', None):
                fig.tight_layout(rect=(0.0, 0.0, 0.76, 0.93), pad=1.1)
            else:
                layout_rect = getattr(fig, '_native_layout_rect', None)
                if layout_rect is not None:
                    fig.tight_layout(rect=layout_rect, pad=1.1)
                else:
                    fig.tight_layout(pad=1.1)

            # Some dense native figures provide an opt-in post-layout audit.
            # Run it after all shared typography and layout adjustments so the
            # report reflects the exact canvas that will be exported.
            post_layout_audit = getattr(fig, '_post_layout_audit', None)
            if callable(post_layout_audit):
                fig._post_layout_audit_result = post_layout_audit(fig)
        except Exception:
            # A third-party figure may expose a non-standard Axes object; export
            # it unchanged rather than making a successful analysis fail.
            pass

        stem, _ = os.path.splitext(filename)
        result_files = []
        for fmt in formats:
            output_path = os.path.join(plots_dir, f'{stem}.{fmt}')
            save_kwargs = {
                'format': fmt, 'bbox_inches': 'tight', 'pad_inches': 0.15,
                'facecolor': 'white',
            }
            if fmt == 'png':
                save_kwargs['dpi'] = dpi
            fig.savefig(output_path, **save_kwargs)
            result_files.append({
                'file_path': output_path,
                'file_type': fmt,
                'category': category,
                'label': label,
            })
        # Analysis runs may produce dozens of panels; release the canvas once
        # both the raster preview and editable vector artifact are persisted.
        try:
            import matplotlib.pyplot as plt
            plt.close(fig)
        except Exception:
            pass
        return result_files

    def build_publication_umap(self, adata, color_key, title='', basis='X_umap'):
        """Build a compact, publication-oriented UMAP Matplotlib figure.

        Scanpy's default figure is excellent for exploration, but its point size
        and legend/colorbar scale poorly when embedded in a web result card. This
        renderer keeps the shared Nature palette while adapting marker size and
        layout to the number of cells.
        """
        import math
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt

        coords = np.asarray(adata.obsm[basis])[:, :2]
        viz = self.params.get('_visualization', {})
        n_obs = max(1, int(coords.shape[0]))
        point_size = float(viz.get('umap_point_size', 5))
        marker_size = max(7, min(30, point_size * math.sqrt(10000 / n_obs)))
        opacity = float(viz.get('umap_opacity', 0.78))
        font_size = max(9, float(viz.get('font_size', 12)))
        font_family = viz.get('font_family', 'Arial')
        fig, ax = plt.subplots(figsize=(9, 6.8), dpi=150)

        values = adata.obs[color_key] if color_key in adata.obs.columns else None
        is_numeric = values is not None and pd.api.types.is_numeric_dtype(values)
        if values is None:
            ax.scatter(coords[:, 0], coords[:, 1], s=marker_size,
                       c=NATURE_PALETTE[0],
                       alpha=opacity, linewidths=0, rasterized=True)
        elif is_numeric:
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1], s=marker_size,
                c=np.asarray(values, dtype=float), cmap=nature_continuous_cmap(),
                alpha=opacity,
                linewidths=0, rasterized=True,
            )
            colorbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.025,
                                    aspect=32)
            colorbar.ax.tick_params(labelsize=max(8, font_size - 2), width=0.6)
            colorbar.set_label(str(color_key), fontsize=font_size - 1,
                               labelpad=6)
        else:
            categorical = values.astype('category')
            color_map = {
                str(category): NATURE_PALETTE[index % len(NATURE_PALETTE)]
                for index, category in enumerate(categorical.cat.categories)
            }
            for category in categorical.cat.categories:
                mask = np.asarray(categorical == category)
                if not mask.any():
                    continue
                ax.scatter(
                    coords[mask, 0], coords[mask, 1], s=marker_size,
                    color=color_map.get(str(category), '#9aa3b2'),
                    alpha=opacity, linewidths=0, rasterized=True,
                    label=str(category),
                )
            n_categories = len(categorical.cat.categories)
            if n_categories:
                legend = ax.legend(
                    loc='center left', bbox_to_anchor=(1.01, 0.5),
                    frameon=False, fontsize=max(8, font_size - 2),
                    markerscale=1.25, borderaxespad=0,
                    ncol=2 if n_categories > 16 else 1,
                )
                legend.set_title(str(color_key), prop={'size': font_size - 1})

        ax.set_title(title, fontsize=font_size + 2, fontweight='semibold',
                     color='#172033', pad=12, family=font_family)
        ax.set_xlabel('UMAP 1', fontsize=font_size, color='#374151', family=font_family)
        ax.set_ylabel('UMAP 2', fontsize=font_size, color='#374151', family=font_family)
        ax.tick_params(labelsize=max(8, font_size - 2), colors='#4b5563', width=0.6)
        ax.set_facecolor('white')
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.margins(0.035)
        fig.patch.set_facecolor('white')
        fig.subplots_adjust(left=0.09, right=0.82 if values is not None and not is_numeric else 0.94,
                            bottom=0.1, top=0.88)
        return fig

    def export_results(self, adata, output_dir, export_format='h5ad',
                       include_layers=None, include_obsm=None, compression='gzip'):
        """将 adata 导出为指定格式。返回导出文件路径列表。"""
        import os, pandas as pd
        os.makedirs(output_dir, exist_ok=True)
        exported = []

        if export_format == 'h5ad':
            out_path = os.path.join(output_dir, 'result.h5ad')
            adata.write_h5ad(out_path, compression=compression if compression != 'none' else None)
            exported.append(out_path)
        elif export_format == 'csv':
            adata.obs.to_csv(os.path.join(output_dir, 'obs.csv'))
            exported.append(os.path.join(output_dir, 'obs.csv'))
            adata.var.to_csv(os.path.join(output_dir, 'var.csv'))
            exported.append(os.path.join(output_dir, 'var.csv'))
            x_path = os.path.join(output_dir, 'expression.csv')
            if hasattr(adata.X, 'toarray'):
                df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
            else:
                df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
            df.to_csv(x_path, compression='gzip' if compression == 'gzip' else None)
            exported.append(x_path)
        elif export_format == 'loom':
            out_path = os.path.join(output_dir, 'result.loom')
            adata.write_loom(out_path)
            exported.append(out_path)

        if include_obsm:
            for key in include_obsm:
                if key in adata.obsm:
                    obsm_path = os.path.join(output_dir, f'obsm_{key}.csv')
                    pd.DataFrame(adata.obsm[key], index=adata.obs_names).to_csv(obsm_path)
                    exported.append(obsm_path)
        return exported

    def validate_input(self, adata) -> Optional[str]:
        return None

    def load_adata(self, input_path):
        """加载 h5ad 并 remap 基因名。"""
        import scanpy as sc
        from modules.io_utils import remap_var_names
        adata = sc.read_h5ad(input_path)
        return remap_var_names(adata)

    def save_output(self, adata, module_name):
        """保存中间结果到 intermediate/，返回 output_path。"""
        import os
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, f'{module_name}_output.h5ad')
        adata.write_h5ad(output_path)
        return output_path

    def ensure_plots_dir(self):
        """确保 plots/ 目录存在，返回路径。"""
        import os
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        return plots_dir

    def save_plotly_json(self, fig, plots_dir, filename, category, label):
        """将 Plotly figure 保存为 JSON 并返回 result_file dict。"""
        import os, re
        from modules.figure_style import style_plotly_figure
        safe_filename = re.sub(r'[^a-zA-Z0-9_.\-]', '_', filename)
        fpath = os.path.join(plots_dir, safe_filename)
        style_plotly_figure(fig, self.params.get('_visualization', {}))
        fig.write_json(fpath)
        return {'file_path': fpath, 'file_type': 'plotly_json', 'category': category, 'label': label}

    @abstractmethod
    def run(self, input_path: str) -> dict:
        pass
