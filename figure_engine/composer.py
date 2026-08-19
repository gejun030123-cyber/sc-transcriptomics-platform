"""Multi-panel Nature figure composition with one physical-size contract."""

from dataclasses import dataclass, field
from typing import Any, Callable, Optional

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from figure_engine.templates.common import attach_contract


@dataclass
class FigurePanel:
    """One panel rendered by the Director or by an axes-native callback."""

    plot_type: str = ''
    data: Any = None
    spec: Optional[FigureSpec] = None
    label: str = ''
    row: int = 0
    column: int = 0
    rowspan: int = 1
    colspan: int = 1
    draw: Optional[Callable] = None
    metadata: dict = field(default_factory=dict)


class NatureFigureComposer:
    """Compose editable panels without rasterising child figures."""

    def __init__(self, director=None):
        if director is None:
            from figure_engine.director import NatureFigureDirector
            director = NatureFigureDirector()
        self.director = director

    def compose(self, panels, spec: FigureSpec, *, nrows=None, ncols=None,
                shared_legend=True, width_ratios=None, height_ratios=None):
        import matplotlib.pyplot as plt

        panels = list(panels)
        if not panels:
            raise ValueError('Figure Composer 至少需要一个 panel')
        if spec.width != 'double':
            raise ValueError('Multi-panel 主图必须使用 183 mm double-column 画布')
        nrows = int(nrows or max(panel.row + panel.rowspan for panel in panels))
        ncols = int(ncols or max(panel.column + panel.colspan for panel in panels))
        if nrows < 1 or ncols < 1:
            raise ValueError('Composer layout 必须至少为 1 x 1')
        style = get_style(spec.style)
        profile = style.profile(spec.with_updates(plot_type='composite'))

        with style.context(spec):
            fig = plt.figure(figsize=profile.figsize, dpi=profile.dpi)
            grid = fig.add_gridspec(
                nrows, ncols,
                width_ratios=width_ratios, height_ratios=height_ratios,
                left=0.035, right=0.985, bottom=0.14 if shared_legend else 0.055,
                top=0.97, wspace=0.12, hspace=0.14,
            )
            occupied = set()
            subfigures = []
            panel_letters = []
            for index, panel in enumerate(panels):
                cells = {
                    (row, column)
                    for row in range(panel.row, panel.row + panel.rowspan)
                    for column in range(panel.column, panel.column + panel.colspan)
                }
                if occupied.intersection(cells):
                    raise ValueError(f'Panel {index + 1} 与其他 panel 占用相同网格')
                occupied.update(cells)
                subfigure = fig.add_subfigure(
                    grid[panel.row:panel.row + panel.rowspan,
                         panel.column:panel.column + panel.colspan]
                )
                label = panel.label or chr(ord('a') + index)
                subfigure.text(
                    0.008, 0.995, label, transform=subfigure.transSubfigure,
                    ha='left', va='top', fontsize=profile.title_font_pt + 1.0,
                    fontweight='bold', color=style.text, zorder=100,
                )
                if panel.draw is not None:
                    panel.draw(subfigure, style, profile)
                else:
                    panel_spec = panel.spec or self.director.create_spec(panel.plot_type)
                    panel_spec = panel_spec.with_updates(
                        mode=spec.mode, style=spec.style,
                        width='single' if panel.colspan == 1 and ncols > 1 else 'double',
                    )
                    self.director.render_into(panel_spec, panel.data, subfigure)
                # A square correlation panel already identifies samples on
                # the y axis; repeating the same x-axis title competes with
                # the shared legend at final double-column size.
                for axis in subfigure.axes:
                    if axis.get_xlabel() == axis.get_ylabel() == 'Sample':
                        axis.set_xlabel('')
                subfigures.append(subfigure)
                panel_letters.append(label)

            shared_labels = []
            if shared_legend:
                handle_by_label = {}
                for subfigure in subfigures:
                    for axis in subfigure.axes:
                        handles, labels = axis.get_legend_handles_labels()
                        for handle, label in zip(handles, labels):
                            if label and not label.startswith('_'):
                                handle_by_label.setdefault(label, handle)
                        legend = axis.get_legend()
                        if legend is not None:
                            legend.remove()
                    for legend in list(subfigure.legends):
                        labels = [text.get_text() for text in legend.get_texts()]
                        handles = getattr(legend, 'legend_handles', ())
                        for handle, label in zip(handles, labels):
                            if label and not label.startswith('_'):
                                handle_by_label.setdefault(label, handle)
                        legend.remove()
                if handle_by_label:
                    shared_labels = list(handle_by_label)
                    fig.legend(
                        list(handle_by_label.values()), shared_labels,
                        loc='lower center', bbox_to_anchor=(0.5, 0.012),
                        ncol=min(6, len(shared_labels)), frameon=False,
                        handletextpad=0.35, columnspacing=0.8,
                        fontsize=profile.legend_font_pt,
                    )

        fig._nature_panel_labels = panel_letters
        fig._nature_panel_grid = {'rows': nrows, 'columns': ncols}
        return attach_contract(
            fig, spec.with_updates(plot_type='composite'), style,
            encodings={'color': 'multi-panel group, DEG direction, correlation and FDR',
                       'panels': ', '.join(panel_letters),
                       'shared_legend': ', '.join(shared_labels)},
        )
