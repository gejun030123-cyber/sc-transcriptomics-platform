from abc import ABC, abstractmethod
from typing import Callable, Optional

VISUALIZATION_THEMES = {
    'default': {
        'bg_color': 'white',
        'color_palette': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'],
        'font_family': 'Arial',
    },
    'nature': {
        'bg_color': 'white',
        'color_palette': ['#E64B35', '#4DBBD5', '#00A087', '#3C5488',
                          '#F39B7F', '#8491B4', '#91D1C2', '#DC0000'],
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

    @abstractmethod
    def validate_input(self, adata) -> Optional[str]:
        pass

    @abstractmethod
    def run(self, input_path: str) -> dict:
        pass
