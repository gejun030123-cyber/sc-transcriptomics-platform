import re
from abc import ABC, abstractmethod
from typing import Callable, Optional

from modules.figure_style import (
    NATURE_BG,
    NATURE_FONT_FAMILY,
    NATURE_PALETTE,
    normalize_visualization_params,
)

# Result files are persisted in the database and exposed by the result routes.
# Keep the file-type contract in one place so modules, tests and API consumers
# agree on optional spreadsheet exports as well as image/table artifacts.
VALID_RESULT_FILE_TYPES = frozenset({
    'csv', 'xlsx', 'plotly_json', 'png', 'svg', 'pdf', 'tiff', 'jpg', 'jpeg',
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

    def apply_scope(self, adata, module_name=None):
        """将分析限制在指定 obs 列（如 sample_id）的选定值对应的细胞上。

        这是跨模块统一的“在哪个样品/条件下比较”入口：``scope_key`` 指定
        obs 列（例如 sample_id、condition、timepoint），``scope_values`` 为
        逗号/分号/换行分隔的取值列表。留空 scope_key 表示分析全部细胞。
        模块应尽早调用本方法，使后续分组、比较和统计都只基于该子集。
        """
        scope_key = str(self.params.get('scope_key', '') or '').strip()
        raw_values = str(self.params.get('scope_values', '') or '').strip()
        if not scope_key:
            return adata
        if scope_key not in adata.obs.columns:
            raise ValueError(
                f"分析范围列 '{scope_key}' 不在 adata.obs 中；"
                "请选择 sample_id/condition 等真实 obs 列，或留空分析全部细胞。"
            )
        values = [
            item.strip() for item in re.split(r'[,;\n]+', raw_values) if item.strip()
        ]
        if not values:
            raise ValueError(
                f"已填写分析范围列 '{scope_key}'，但未选择具体取值"
                "（scope_values 至少需要一个值）。"
            )
        obs_values = adata.obs[scope_key].astype(str)
        present = [value for value in values if value in set(obs_values)]
        missing = [value for value in values if value not in present]
        if missing:
            self.progress(
                -1,
                f"分析范围 '{scope_key}' 中不存在取值: {missing}；已忽略。",
            )
        if not present:
            raise ValueError(
                f"scope_values 中的取值都不存在于 '{scope_key}'：{values}"
            )
        n_before = adata.n_obs
        adata = adata[obs_values.isin(present)].copy()
        self.progress(
            -1,
            f"分析范围限制: {scope_key} ∈ {present}，"
            f"保留 {adata.n_obs}/{n_before} 个细胞。",
        )
        return adata
    def get_plotly_layout(self, title='', **overrides):
        """根据 self.params['_visualization'] 生成统一的 Plotly layout dict。"""
        viz = normalize_visualization_params(self.params.get('_visualization', {}))
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
        viz = normalize_visualization_params(self.params.get('_visualization', {}))
        return {
            'umap_point_size': viz.get('umap_point_size', 5),
            'umap_opacity': viz.get('umap_opacity', 0.7),
            'umap_legend_fontsize': viz.get('umap_legend_fontsize', 10),
            'umap_label_categories': viz.get('umap_label_categories', False),
            'umap_hide_axes': viz.get('umap_hide_axes', True),
        }

    def export_static_fig(self, fig, plots_dir, filename, export_formats=None):
        """将 Plotly figure 导出为静态图片（SVG/PNG/PDF/TIFF）。依赖 kaleido。"""
        if not export_formats:
            return []
        import os
        from modules.figure_style import style_plotly_figure
        style_plotly_figure(fig, self.params.get('_visualization', {}))
        exported = []
        for fmt in export_formats:
            if fmt in ('svg', 'png', 'pdf', 'tiff'):
                try:
                    import plotly.io as pio
                    out_path = os.path.join(plots_dir, filename.replace('.json', f'.{fmt}'))
                    if fmt == 'tiff':
                        from PIL import Image
                        png_bytes = pio.to_image(fig, format='png', engine='kaleido', scale=2)
                        import io
                        with Image.open(io.BytesIO(png_bytes)) as image:
                            image.convert('RGB').save(out_path, format='TIFF',
                                                      compression='tiff_lzw', dpi=(300, 300))
                    else:
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

        raw_viz = self.params.get('_visualization', {})
        viz = normalize_visualization_params(raw_viz)
        nature_spec = getattr(fig, '_nature_spec_object', None)
        if raw_viz.get('static_formats') or raw_viz.get('export_formats'):
            formats = viz.get('static_formats', ('png', 'svg'))
        elif formats is None:
            formats = (
                tuple(nature_spec.formats)
                if nature_spec is not None
                else viz.get('static_formats', ('png', 'svg'))
            )
        formats = [fmt.lower() for fmt in formats
                   if fmt.lower() in ('png', 'svg', 'pdf', 'tiff')]
        if not formats:
            # 用户可以在界面取消全部静态图格式。此时没有任何文件需要写出，
            # 但画布已经创建；不关闭就会在长任务里持续累积（matplotlib 的
            # pyplot 注册表持有强引用）。
            try:
                import matplotlib.pyplot as plt
                plt.close(fig)
            except Exception:
                pass
            return []

        # Figures created by the deterministic engine already carry an exact
        # physical-size/style contract.  Reapplying the legacy post-processor
        # or tight bounding boxes would invalidate that contract, so route them
        # through the same exporter and readiness audit used by Bulk RNA-seq.
        if nature_spec is not None:
            from figure_engine import export_registered_figure

            stem, _ = os.path.splitext(filename)
            results_dir = os.path.join(self.project_dir, 'results')
            os.makedirs(results_dir, exist_ok=True)
            qa_path = os.path.join(results_dir, f'{stem}_nature_readiness.json')
            result_files, _ = export_registered_figure(
                fig, os.path.join(plots_dir, stem),
                nature_spec.with_updates(formats=tuple(formats)),
                category=category, label=label, formats=tuple(formats),
                qa_path=qa_path,
            )
            try:
                import matplotlib.pyplot as plt
                plt.close(fig)
            except Exception:
                pass
            return result_files

        # Normalize the native Scanpy/Matplotlib canvas for publication output.
        # This changes the figure itself before both PNG and SVG are written, so
        # the browser preview and the vector download share the same composition.
        try:
            # Native builders choose their own aspect ratio (for example a
            # multi-resolution UMAP grid or a long cluster composition panel).
            # Callers can preserve that aspect explicitly; content-aware native
            # canvases are also kept intact when the UI defaults would shrink
            # their label area.
            if ('figure_width' in raw_viz or 'figure_height' in raw_viz) and not preserve_aspect:
                current_width, current_height = fig.get_size_inches()
                width = max(6.5, float(viz.get('figure_width', current_width * 100)) / 100)
                height = max(4.8, float(viz.get('figure_height', current_height * 100)) / 100)
                # Several native single-cell figures size their canvas from the
                # number/length of labels (heatmaps, composition bars and
                # marker dot-plots).  Shrinking those canvases to the generic
                # 800x500 UI default leaves the font size unchanged while
                # removing the space allocated to each label, which produces
                # the visible text collisions on the result page.  The user
                # dimensions remain the lower bound for compact figures, but a
                # content-aware canvas is never compressed below its layout.
                if current_width > width or current_height > height:
                    # Preserve both dimensions together; independently taking
                    # max(width) and max(height) would subtly distort a wide
                    # panel and move labels relative to the plotted data.
                    width, height = float(current_width), float(current_height)
                fig.set_size_inches(width, height, forward=True)
            font_size = max(9, float(viz.get('font_size', 12)))
            font_family = viz.get('font_family', NATURE_FONT_FAMILY)
            is_dark = str(viz.get('bg_color', 'white')).lower() in {
                '#1a1a2e', '#111827', '#0f172a'
            }
            text_color = '#F8FAFC' if is_dark else '#1f2937'
            axis_color = '#CBD5E1' if is_dark else '#374151'
            spine_color = '#64748B' if is_dark else '#c7cdd6'
            apply_matplotlib_style(fig, viz)
            for axis in getattr(fig, 'axes', []):
                axis.set_facecolor(viz.get('bg_color', 'white'))
                axis.tick_params(labelsize=max(8, font_size - 2), width=0.7,
                                 colors=axis_color)
                axis.xaxis.label.set_size(font_size)
                axis.yaxis.label.set_size(font_size)
                axis.xaxis.label.set_color(text_color)
                axis.yaxis.label.set_color(text_color)
                title = axis.title
                title.set_fontsize(font_size + 1)
                title.set_fontweight('semibold')
                title.set_color(text_color)
                title.set_fontfamily(_font_for_text(title.get_text(), font_family))
                for spine_name, spine in axis.spines.items():
                    spine.set_linewidth(0.65)
                    spine.set_color(spine_color)
                    spine.set_visible(spine_name in ('left', 'bottom'))
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
                'facecolor': viz.get('bg_color', 'white'),
            }
            if fmt in {'png', 'tiff'}:
                save_kwargs['dpi'] = dpi
            if fmt == 'tiff':
                save_kwargs['pil_kwargs'] = {'compression': 'tiff_lzw'}
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

    def build_publication_embedding(self, adata, color_key, title='', basis='X_umap',
                                    x_label=None, y_label=None):
        """Build a deterministic cell-level embedding with Nature Figure Engine.

        Categorical ``obs`` values use the fixed cell-type palette. Numeric
        ``obs`` values and genes use the shared continuous expression scale.
        The dense cell layer is rasterized while text remains editable.
        """
        import numpy as np
        import pandas as pd
        from figure_engine import NatureFigureDirector

        if basis not in adata.obsm:
            raise KeyError(f"AnnData.obsm 缺少嵌入坐标: {basis}")
        coordinates = np.asarray(adata.obsm[basis])[:, :2]
        values = None
        value_type = 'categorical'
        category_order = []
        if color_key in adata.obs.columns:
            series = adata.obs[color_key]
            values = series.to_numpy()
            key_hint = str(color_key).strip().lower()
            discrete_hint = (
                key_hint in {'batch', 'phase', 'leiden', 'louvain', 'celltype', 'cell_type'}
                or 'cluster' in key_hint or key_hint.endswith('_group')
            )
            if pd.api.types.is_numeric_dtype(series) and not discrete_hint:
                value_type = 'continuous'
            else:
                values = series.astype(str).to_numpy()
                if isinstance(series.dtype, pd.CategoricalDtype):
                    category_order = [str(value) for value in series.cat.categories]
                else:
                    category_order = list(dict.fromkeys(values.tolist()))
        elif color_key in adata.var_names:
            expression = adata[:, [color_key]].X
            if hasattr(expression, 'toarray'):
                expression = expression.toarray()
            values = np.asarray(expression, dtype=float).reshape(-1)
            value_type = 'continuous'

        axis_prefix = basis.removeprefix('X_').upper()
        x_label = x_label or f'{axis_prefix} 1'
        y_label = y_label or f'{axis_prefix} 2'
        viz = normalize_visualization_params(self.params.get('_visualization', {}))
        n_categories = len(category_order) if value_type == 'categorical' else 0
        width = 'double' if n_categories > 12 else 'single'
        director = NatureFigureDirector()
        spec = director.spec_from_params(
            'embedding', self.params, width=width,
            title=title or f'{axis_prefix} by {color_key}',
            evidence_role='discovery',
        )
        return director.render(spec, {
            'coordinates': coordinates,
            'values': values,
            'value_type': value_type,
            'value_label': str(color_key),
            'category_order': category_order,
            'x_label': x_label,
            'y_label': y_label,
            'direct_labels': viz.get('umap_label_categories', False),
            'hide_axes': viz.get('umap_hide_axes', True),
        })

    def build_publication_umap(self, adata, color_key, title='', basis='X_umap'):
        """Backward-compatible wrapper for the shared embedding renderer."""
        return self.build_publication_embedding(
            adata, color_key, title=title, basis=basis,
            x_label='UMAP 1', y_label='UMAP 2',
        )

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
        """加载 h5ad 并 remap 基因名。

        重映射可能产生重复基因名（例如数据中已存在 ``G_1``，新后缀又生成了
        同名条目），而 ``Index.get_indexer``/``get_loc`` 要求索引唯一；在
        入口统一保证 var_names 唯一，避免分析跑到末尾才抛 InvalidIndexError。
        """
        import scanpy as sc
        from modules.io_utils import remap_var_names
        adata = sc.read_h5ad(input_path)
        adata = remap_var_names(adata)
        if not adata.var_names.is_unique:
            adata.var_names_make_unique()
        return adata

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

    def remove_plot_artifacts(self, *filenames):
        """Remove stale plot artifacts through the shared output boundary.

        Modules should not need to know how the project output directory is
        laid out.  Cleanup remains deliberately explicit and limited to the
        filenames supplied by the caller.
        """
        import os

        plots_dir = self.ensure_plots_dir()
        removed = []
        for filename in filenames:
            if not filename:
                continue
            path = os.path.join(plots_dir, os.path.basename(str(filename)))
            if not os.path.exists(path):
                continue
            try:
                os.remove(path)
                removed.append(path)
            except OSError as exc:
                raise OSError(f"无法清理旧图 {path}: {exc}") from exc
        return removed

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
