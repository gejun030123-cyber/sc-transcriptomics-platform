"""NatureUpSet：交集大小和集合隶属关系的完整 UpSet 模板。"""

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import attach_contract


def _intersection_frame(data):
    if isinstance(data, dict) and 'sets' not in data and all(
            isinstance(value, (set, list, tuple, np.ndarray)) for value in data.values()):
        names = list(data)
        memberships = {name: set(map(str, data[name])) for name in names}
        universe = set().union(*memberships.values()) if memberships else set()
        counts = {}
        for member in universe:
            signature = tuple(name for name in names if member in memberships[name])
            if signature:
                counts[signature] = counts.get(signature, 0) + 1
        return pd.DataFrame([
            {'sets': signature, 'count': count} for signature, count in counts.items()
        ])
    frame = pd.DataFrame(data).copy()
    if not {'sets', 'count'}.issubset(frame.columns):
        raise ValueError('NatureUpSet 需要 {sets, count} 交集表或 set-name -> members 映射')
    frame['sets'] = frame['sets'].map(
        lambda value: tuple(str(item).strip() for item in value)
        if isinstance(value, (list, tuple, set))
        else tuple(item.strip() for item in str(value).replace('∩', '|').split('|') if item.strip())
    )
    frame['count'] = pd.to_numeric(frame['count'], errors='coerce')
    return frame.dropna(subset=['count'])


class NatureUpSet:
    plot_type = 'upset'

    def render(self, data, spec: FigureSpec, container=None):
        import matplotlib.pyplot as plt

        style = get_style(spec.style)
        profile = style.profile(spec)
        frame = _intersection_frame(data)
        frame = frame[frame['count'] > 0].sort_values(
            ['count'], ascending=False).head(spec.top_intersections).reset_index(drop=True)
        if frame.empty:
            raise ValueError('NatureUpSet 没有非零交集')
        set_names = list(dict.fromkeys(name for values in frame['sets'] for name in values))
        x = np.arange(len(frame))

        with style.context(spec):
            if container is None:
                fig = plt.figure(figsize=profile.figsize, dpi=profile.dpi)
            else:
                fig = container
            grid = fig.add_gridspec(2, 1, height_ratios=(0.62, 0.38),
                                    left=0.13, right=0.98, bottom=0.18, top=0.89,
                                    hspace=0.05)
            ax_bar = fig.add_subplot(grid[0])
            ax_matrix = fig.add_subplot(grid[1], sharex=ax_bar)
            bars = ax_bar.bar(x, frame['count'].to_numpy(dtype=float), width=0.68,
                              color=style.signal_blue, edgecolor='none')
            ax_bar.bar_label(bars, labels=[f'{int(value):,}' for value in frame['count']],
                             padding=2, fontsize=max(5.2, profile.tick_font_pt - 0.6),
                             color=style.text)
            positive_counts = frame['count'].to_numpy(dtype=float)
            use_symlog = (positive_counts.max() / max(positive_counts.min(), 1.0)) >= 50
            if use_symlog:
                ax_bar.set_yscale('symlog', linthresh=1.0, linscale=0.8)
            ax_bar.set_ylabel('Intersection size' + (' (symlog)' if use_symlog else ''))
            ax_bar.set_title(spec.title or 'Set intersections', loc='left', pad=5)
            ax_bar.tick_params(axis='x', labelbottom=False, length=0)
            style.apply_axis(ax_bar, profile)
            ax_bar.spines['bottom'].set_visible(False)
            ax_bar.margins(y=0.14)

            for column, membership in enumerate(frame['sets']):
                rows = [set_names.index(name) for name in membership]
                ax_matrix.scatter(np.full(len(set_names), column), np.arange(len(set_names)),
                                  s=10, color=style.neutral, linewidths=0, zorder=1)
                if rows:
                    ax_matrix.scatter(np.full(len(rows), column), rows, s=16,
                                      color=style.signal_blue, linewidths=0, zorder=3)
                    if len(rows) > 1:
                        ax_matrix.plot([column, column], [min(rows), max(rows)],
                                       color=style.signal_blue, linewidth=0.8, zorder=2)
            ax_matrix.set_yticks(np.arange(len(set_names)), set_names)
            ax_matrix.set_xticks(x, [str(index + 1) for index in x], rotation=0)
            ax_matrix.set_xlabel('Intersection rank')
            ax_matrix.tick_params(axis='y', length=0)
            ax_matrix.set_ylim(len(set_names) - 0.5, -0.5)
            style.apply_axis(ax_matrix, profile)
            for name in ('top', 'right', 'left'):
                ax_matrix.spines[name].set_visible(False)

        if container is not None:
            return container
        warnings = []
        if len(set_names) > 10:
            warnings.append(f'UpSet 包含 {len(set_names)} 个集合，主文建议拆分或筛选。')
        return attach_contract(
            fig, spec, style, semantic_warnings=warnings,
            encodings={'bar': 'intersection size', 'matrix': 'set membership'},
        )
