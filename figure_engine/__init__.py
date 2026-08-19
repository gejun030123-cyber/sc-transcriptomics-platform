"""确定性的 Nature Portfolio 科研绘图引擎。"""

from .director import NatureFigureDirector
from .composer import FigurePanel, NatureFigureComposer
from .export import export_figure, export_registered_figure
from .spec import FigureSpec
from .validator import FigureReadinessError, FigureValidator, ValidationReport
from .accessibility import audit_palette, simulate_image, simulate_rgb
from .visual_regression import compare_signatures, figure_signature, hash_distance

__all__ = [
    'FigureSpec', 'FigureValidator', 'NatureFigureDirector',
    'ValidationReport', 'export_figure', 'export_registered_figure',
    'FigurePanel', 'NatureFigureComposer',
    'FigureReadinessError', 'audit_palette', 'simulate_image', 'simulate_rgb',
    'compare_signatures', 'figure_signature', 'hash_distance',
]
