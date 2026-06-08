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

    @abstractmethod
    def validate_input(self, adata) -> Optional[str]:
        pass

    @abstractmethod
    def run(self, input_path: str) -> dict:
        pass
