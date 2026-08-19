"""Fixed publication templates for secondary Bulk RNA-seq diagnostics.

The main scientific plots have dedicated renderers.  This module deliberately
keeps the smaller diagnostic family in one deterministic renderer so that QC,
normalisation checks and time-course summaries cannot silently fall back to
Matplotlib defaults.  ``data['kind']`` selects a scientific view; visual
constants remain owned by :class:`NatureStyle` and ``FigureSpec``.
"""

from __future__ import annotations

import math

import numpy as np

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style

from .common import as_1d, attach_contract, label_tick_indices, robust_symmetric_limit


def _finite(values):
    return np.asarray(values, dtype=float).reshape(-1)


def _labels(values, length=None):
    values = as_1d(values, length, '')
    return np.asarray([str(value) for value in values], dtype=object)


def _short_labels(values, maximum=10):
    values = [str(value) for value in values]
    if len(values) <= maximum:
        return values
    return [f'S{i + 1}' for i in range(len(values))]


class NatureDiagnostic:
    """Render normalisation, QC, exploratory and time-course summaries."""

    plot_type = 'diagnostic'

    def render(self, data, spec: FigureSpec, container=None):
        kind = str((data or {}).get('kind') or spec.extra.get('kind') or spec.plot_type).lower()
        methods = {
            'normalization_library': self._library,
            'library_size': self._library,
            'normalization_boxplot': self._distribution,
            'distribution': self._distribution,
            'normalization_pca': self._pca_compare,
            'pca_compare': self._pca_compare,
            'qc_overview': self._qc_overview,
            'variance': self._variance,
            'bar': self._bar,
            'boxplot': self._boxplot,
            'pairs': self._pairs,
            'violin': self._violin,
            'qq': self._qq,
            'trajectory': self._trajectory,
        }
        if kind in {'heatmap', 'correlation_heatmap'}:
            from .heatmap import NatureHeatmap

            return NatureHeatmap().render(data, spec.with_updates(plot_type='heatmap'), container=container)
        try:
            return methods[kind](data, spec, container)
        except KeyError as exc:
            raise ValueError(f'NatureDiagnostic 不支持 kind={kind!r}') from exc

    def _canvas(self, spec, container, nrows=1, ncols=1, *, height=None):
        import matplotlib.pyplot as plt

        style = get_style(spec.style)
        # The export layer always restores ``spec.height_mm``.  Build the
        # canvas at that same physical size when it is explicit; otherwise a
        # template may be laid out at one height and then silently compressed
        # during export, which is a common source of title/label collisions.
        effective_height = spec.height_mm if spec.height_mm is not None else height
        profile = style.profile(
            spec.with_updates(height_mm=effective_height)
            if effective_height is not None else spec
        )
        with style.context(spec):
            if container is None:
                fig, axes = plt.subplots(
                    nrows, ncols, figsize=profile.figsize, dpi=profile.dpi,
                    squeeze=False,
                )
            else:
                fig = container
                axes = np.asarray(fig.subplots(nrows, ncols, squeeze=False))
        return fig, axes, style, profile

    @staticmethod
    def _finish(fig, axes, style, profile, spec, title, *, warnings=(), encodings=None,
                left=0.16, right=0.96, bottom=0.17, top=0.88, hspace=0.34):
        for ax in np.asarray(axes).reshape(-1):
            if ax.get_visible():
                style.apply_axis(ax, profile)
        fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                            wspace=0.28, hspace=hspace)
        if title:
            # Keep a dedicated title band above all axes titles.  A figure
            # title at 0.965 is close enough to a first-row ``Axes.title`` to
            # overlap after rasterisation at the final submission size.
            fig.text(left, 0.988, title, ha='left', va='top',
                     fontsize=profile.title_font_pt, fontweight='semibold',
                     color=style.text)
        return attach_contract(fig, spec, style, semantic_warnings=warnings,
                               encodings=encodings or {})

    def _library(self, data, spec, container):
        labels = _labels(data.get('sample_labels'))
        raw = _finite(data.get('raw_values'))
        normalized = _finite(data.get('normalized_values'))
        if len(labels) != len(raw) or len(labels) != len(normalized):
            raise ValueError('normalization_library 的样本和值长度不一致')
        fig, axes, style, profile = self._canvas(spec, container, 1, 2, height=60)
        idx = np.arange(len(labels))
        shown = _short_labels(labels, 10)
        shown = [label if len(label) <= 12 else label[:11] + '…' for label in shown]
        for ax, values, subtitle, color, ylabel in zip(
            axes.ravel(), (raw, normalized),
            ('Before normalization', 'After normalization'),
            (style.signal_orange, style.signal_blue),
            (str(data.get('raw_ylabel') or 'Library size'),
             str(data.get('normalized_ylabel') or 'Library size')),
        ):
            ax.bar(idx, values, color=color, alpha=0.86, edgecolor='white', linewidth=0.25)
            ax.set_title(subtitle, loc='left', pad=3)
            ax.set_ylabel(ylabel)
            ticks = label_tick_indices(len(labels), 8)
            ax.set_xticks(ticks, [shown[i] for i in ticks], rotation=35, ha='right')
            ax.set_xlabel('Sample')
        return self._finish(fig, axes, style, profile, spec, spec.title or 'Normalisation check',
                            left=0.11, bottom=0.28, top=0.88)

    def _distribution(self, data, spec, container):
        labels = _labels(data.get('sample_labels'))
        raw = np.asarray(data.get('raw_matrix'), dtype=float)
        normalized = np.asarray(data.get('normalized_matrix'), dtype=float)
        if raw.ndim != 2 or normalized.ndim != 2 or raw.shape[0] != len(labels):
            raise ValueError('normalization_boxplot 需要样本 x 基因矩阵')
        fig, axes, style, profile = self._canvas(spec, container, 1, 2, height=76)
        shown = _short_labels(labels, 10)
        ticks = label_tick_indices(len(labels), 8)
        for ax, matrix, subtitle, color in zip(
            axes.ravel(), (raw, normalized),
            ('Before transformation', 'After transformation'),
            (style.signal_orange, style.signal_blue),
        ):
            box = ax.boxplot(matrix.T, patch_artist=True, showfliers=False,
                             widths=0.72,
                             medianprops={'color': style.text, 'linewidth': 0.7},
                             whiskerprops={'color': style.neutral_dark, 'linewidth': 0.5},
                             capprops={'color': style.neutral_dark, 'linewidth': 0.5})
            for patch in box['boxes']:
                patch.set_facecolor(color)
                patch.set_edgecolor(color)
                patch.set_alpha(0.52)
            ax.set_title(subtitle, loc='left', pad=3)
            ax.set_ylabel(str(data.get('ylabel') or 'Expression'))
            ax.set_xticks(ticks + 1, [shown[i] for i in ticks], rotation=35, ha='right')
            ax.set_xlabel('Sample')
        note = str(data.get('note') or '')
        if note:
            fig.text(0.11, 0.02, note, ha='left', va='bottom', fontsize=max(5.5, profile.tick_font_pt - 0.5),
                     color=style.muted_text)
        return self._finish(fig, axes, style, profile, spec, spec.title or 'Sample expression distributions',
                            left=0.11, bottom=0.31, top=0.88)

    def _pca_compare(self, data, spec, container):
        coords = [np.asarray(item, dtype=float)[:, :2] for item in data.get('coordinates', [])]
        if len(coords) != 2:
            raise ValueError('pca_compare 需要两个 n_samples x 2 坐标矩阵')
        n = coords[0].shape[0]
        groups = _labels(data.get('groups'), n)
        colors = dict(data.get('group_colors') or get_style(spec.style).group_colors(groups))
        labels = _labels(data.get('sample_labels'), n)
        explained = data.get('explained_variance') or [None, None, None, None]
        fig, axes, style, profile = self._canvas(spec, container, 1, 2, height=70)
        for index, (ax, xy) in enumerate(zip(axes.ravel(), coords)):
            for group in dict.fromkeys(groups):
                mask = groups == group
                ax.scatter(xy[mask, 0], xy[mask, 1], s=profile.marker_size_pt2,
                           color=colors[str(group)], edgecolor='white', linewidth=0.4,
                           alpha=0.9, label=str(group))
            xv = explained[index * 2] if len(explained) > index * 2 else None
            yv = explained[index * 2 + 1] if len(explained) > index * 2 + 1 else None
            def _variance_suffix(value):
                try:
                    value = float(value)
                except (TypeError, ValueError):
                    return ''
                if not np.isfinite(value):
                    return ''
                return f' ({value * 100:.1f}%)' if value <= 1 else f' ({value:.1f}%)'

            suffix_x = _variance_suffix(xv)
            suffix_y = _variance_suffix(yv)
            ax.set_xlabel(f'PC1{suffix_x}')
            ax.set_ylabel(f'PC2{suffix_y}')
            ax.set_title(str(data.get('panel_titles', ['Before', 'After'])[index]), loc='left', pad=3)
        warning = []
        if len(dict.fromkeys(groups)) > 8:
            warning.append('PCA 对比包含超过 8 个颜色组，建议拆分分组或使用 batch marker。')
        return self._finish(fig, axes, style, profile, spec, spec.title or 'PCA before/after comparison',
                            warnings=warning, encodings={'color': 'group'}, left=0.14, right=0.96,
                            bottom=0.22, top=0.88)

    def _qc_overview(self, data, spec, container):
        metrics = []
        for item in data.get('metrics', []):
            values = _finite(item.get('values'))
            if values.size == 0 or not np.isfinite(values).any():
                continue
            finite = values[np.isfinite(values)]
            if finite.size == 0 or np.nanmax(finite) - np.nanmin(finite) <= 1e-12:
                continue
            metrics.append((values, str(item.get('title') or ''), str(item.get('ylabel') or '')))
        lib = _finite(data.get('library_size'))
        genes = _finite(data.get('detected_genes'))
        n_panels = len(metrics) + (1 if len(lib) and len(genes) else 0)
        n_panels = max(1, n_panels)
        ncols = 2
        nrows = int(math.ceil(n_panels / ncols))
        height = max(82.0, min(170.0, 43.0 * nrows))
        fig, axes, style, profile = self._canvas(spec, container, nrows, ncols, height=height)
        flat = axes.ravel()
        sample_idx = np.arange(len(lib) if len(lib) else len(metrics[0][0]))
        colors = data.get('pass_colors')
        if colors is None:
            colors = [style.signal_blue] * len(sample_idx)
        for ax, (values, title, ylabel) in zip(flat, metrics):
            ax.bar(sample_idx[:len(values)], values, color=style.signal_blue, alpha=0.88,
                   edgecolor='white', linewidth=0.2)
            ax.set_title(title, loc='left', pad=3)
            ax.set_ylabel(ylabel)
            ticks = label_tick_indices(len(values), 8)
            labels = _labels(data.get('sample_labels'), len(values))
            ax.set_xticks(ticks, [_short_labels(labels, 10)[i] for i in ticks], rotation=35, ha='right')
        if len(lib) and len(genes):
            ax = flat[len(metrics)]
            ax.scatter(lib, genes, c=colors, s=profile.marker_size_pt2,
                       alpha=0.88, edgecolor='white', linewidth=0.35)
            ax.set_title('Library size vs detected genes', loc='left', pad=3)
            ax.set_xlabel('Library size')
            ax.set_ylabel('Detected genes')
        for ax in flat[n_panels:]:
            ax.set_visible(False)
        warning = []
        omitted = str(data.get('omitted_metrics') or '')
        if omitted:
            warning.append(f'QC overview omitted invariant metric(s): {omitted}')
        return self._finish(fig, axes, style, profile, spec, spec.title or 'Bulk RNA-seq QC overview',
                            warnings=warning, encodings={'color': 'pass/fail'}, left=0.14, right=0.96,
                            # Three rows need a visibly larger title band than
                            # a single-row diagnostic; this also keeps the
                            # first-row panel title away from the figure title
                            # at 89/183 mm export sizes.
                            bottom=0.28, top=0.84, hspace=0.68)

    def _variance(self, data, spec, container):
        values = _finite(data.get('variance'))
        fig, axes, style, profile = self._canvas(spec, container, 1, 1, height=70)
        ax = axes.ravel()[0]
        x = np.arange(len(values))
        ax.bar(x, values, color=style.signal_blue, alpha=0.86, label='Variance ratio')
        cumulative = np.cumsum(values)
        ax2 = ax.twinx()
        ax2.plot(x, cumulative, color=style.signal_red, marker='o', markersize=3.2,
                 linewidth=1.0, label='Cumulative')
        labels = _labels(data.get('labels'), len(values))
        ticks = label_tick_indices(len(labels), 8 if spec.width == 'single' else 14)
        ax.set_xticks(ticks, [labels[index] for index in ticks], rotation=32, ha='right')
        ax.set_xlabel('Principal component')
        ax.set_ylabel('Variance ratio')
        ax2.set_ylabel('Cumulative ratio')
        style.apply_axis(ax2, profile)
        ax2.spines['left'].set_visible(False)
        handles, names = ax.get_legend_handles_labels()
        handles2, names2 = ax2.get_legend_handles_labels()
        # Keep the two-entry legend inside the physical canvas.  An external
        # bbox made the second label clip at the right edge of 89 mm exports.
        fig.legend(handles + handles2, names + names2, frameon=False,
                   loc='upper right', bbox_to_anchor=(0.96, 0.90), ncol=2,
                   borderaxespad=0, handlelength=1.2, columnspacing=0.8,
                   fontsize=profile.legend_font_pt)
        return self._finish(fig, axes, style, profile, spec, spec.title or 'PCA variance explained',
                            left=0.16, right=0.88, bottom=0.30, top=0.80)

    def _boxplot(self, data, spec, container):
        groups = [str(value) for value in data.get('groups', [])]
        values = [np.asarray(value, dtype=float) for value in data.get('values', [])]
        if len(groups) != len(values):
            raise ValueError('boxplot groups/values 长度不一致')
        fig, axes, style, profile = self._canvas(spec, container, 1, 1, height=68)
        ax = axes.ravel()[0]
        box = ax.boxplot(values, tick_labels=groups, patch_artist=True, showfliers=False)
        for index, patch in enumerate(box['boxes']):
            patch.set_facecolor(style.categorical_palette[index % len(style.categorical_palette)])
            patch.set_edgecolor(style.categorical_palette[index % len(style.categorical_palette)])
            patch.set_alpha(0.58)
        for index, vals in enumerate(values, 1):
            vals = vals[np.isfinite(vals)]
            jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) else np.asarray([])
            ax.scatter(np.full(len(vals), index) + jitter, vals, s=profile.marker_size_pt2 * 0.42,
                       color=style.neutral_dark, alpha=0.6, linewidths=0, zorder=3)
        ax.set_ylabel(str(data.get('ylabel') or 'Expression'))
        ax.set_title(spec.title or str(data.get('title') or ''), loc='left', pad=4)
        return self._finish(fig, axes, style, profile, spec, '', left=0.18, bottom=0.20, top=0.88)

    def _bar(self, data, spec, container):
        labels = [str(value) for value in data.get('labels', [])]
        values = _finite(data.get('values'))
        if len(labels) != len(values):
            raise ValueError('bar labels/values 长度不一致')
        fig, axes, style, profile = self._canvas(spec, container, 1, 1, height=78)
        ax = axes.ravel()[0]
        colors = [style.signal_red if value >= 0 else style.signal_blue for value in values]
        x = np.arange(len(values))
        ax.bar(x, values, color=colors, alpha=0.82, edgecolor='white', linewidth=0.25)
        ax.axhline(0, color=style.neutral_dark, linewidth=0.55)
        # Loading gene names are long and diagonal labels are especially easy
        # to collide at the final physical size.  Keep a sparse, deterministic
        # subset for publication panels; the complete ranked list remains in
        # the result table.
        tick_cap = 6 if spec.mode == 'nature_portfolio' else (8 if spec.width == 'single' else 8)
        ticks = label_tick_indices(len(labels), tick_cap)
        shown = [label if len(label) <= 9 else label[:8] + '…' for label in labels]
        ax.set_xticks(ticks, [shown[i] for i in ticks], rotation=50, ha='right')
        ax.set_xlabel(str(data.get('xlabel') or 'Gene'))
        ax.set_ylabel(str(data.get('ylabel') or 'Value'))
        return self._finish(fig, axes, style, profile, spec, spec.title or '',
                            left=0.16, bottom=0.30, top=0.88)

    def _pairs(self, data, spec, container):
        labels = [str(value) for value in data.get('labels', [])]
        values = [np.asarray(value, dtype=float) for value in data.get('values', [])]
        selected = []
        for label, value in zip(labels, values):
            finite = value[np.isfinite(value)]
            if finite.size and np.nanmax(finite) - np.nanmin(finite) > 1e-12:
                selected.append((label, value))
        if not selected:
            selected = list(zip(labels, values))[:1]
        labels, values = zip(*selected)
        n = len(values)
        fig, axes, style, profile = self._canvas(spec, container, n, n, height=max(78, min(170, n * 27)))
        groups = _labels(data.get('groups'), len(values[0]))
        colors = dict(data.get('group_colors') or style.group_colors(groups))
        markers = dict(data.get('group_markers') or style.batch_markers(groups))
        for row in range(n):
            for col in range(n):
                ax = axes[row, col]
                if row == col:
                    ax.hist(values[row], bins=18, color=style.signal_blue, alpha=0.72,
                            edgecolor='white', linewidth=0.2)
                elif row > col:
                    for group in dict.fromkeys(groups):
                        mask = groups == group
                        ax.scatter(values[col][mask], values[row][mask], s=profile.marker_size_pt2 * 0.35,
                                   color=colors[str(group)], marker=markers[str(group)],
                                   alpha=0.72, linewidths=0)
                else:
                    ax.set_visible(False)
                if row == n - 1:
                    ax.set_xlabel(labels[col])
                if col == 0 and row > col:
                    ax.set_ylabel(labels[row])
        return self._finish(fig, axes, style, profile, spec, spec.title or 'QC metric pairs',
                            encodings={'color': 'group', 'marker': 'group'},
                            left=0.16, right=0.98, bottom=0.14, top=0.90)

    def _violin(self, data, spec, container):
        groups = [str(value) for value in data.get('groups', [])]
        metrics = [(str(item.get('label')), [np.asarray(v, dtype=float) for v in item.get('values', [])])
                   for item in data.get('metrics', [])]
        metrics = [(label, vals) for label, vals in metrics if any(np.isfinite(v).any() for v in vals)]
        n = max(1, len(metrics))
        ncols = 2
        nrows = int(math.ceil(n / ncols))
        fig, axes, style, profile = self._canvas(spec, container, nrows, ncols, height=max(70, 40 * nrows))
        shown_groups = _short_labels(groups, 8)
        group_ticks = label_tick_indices(len(groups), 8) + 1
        for ax, (label, vals) in zip(axes.ravel(), metrics):
            parts = ax.violinplot(vals, positions=np.arange(1, len(vals) + 1),
                                  showmeans=False, showextrema=False)
            for index, body in enumerate(parts['bodies']):
                color = (style.signal_blue if len(groups) > 3
                         else style.categorical_palette[index % len(style.categorical_palette)])
                body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(0.56)
            ax.set_title(label, loc='left', pad=3)
            ax.set_xticks(group_ticks, [shown_groups[i - 1] for i in group_ticks],
                          rotation=35, ha='right')
        for ax in axes.ravel()[len(metrics):]:
            ax.set_visible(False)
        return self._finish(fig, axes, style, profile, spec, spec.title or 'QC metrics by group',
                            left=0.14, right=0.98, bottom=0.26, top=0.84, hspace=0.58)

    def _qq(self, data, spec, container):
        x = _finite(data.get('theoretical'))
        y = _finite(data.get('observed'))
        fig, axes, style, profile = self._canvas(spec, container, 1, 1, height=68)
        ax = axes.ravel()[0]
        ax.scatter(x, y, s=profile.marker_size_pt2 * 0.42, color=style.signal_blue,
                   alpha=0.72, linewidths=0, rasterized=True)
        max_val = max(float(np.nanmax(x)), float(np.nanmax(y))) * 1.08
        ax.plot([0, max_val], [0, max_val], color=style.signal_red, linestyle=(0, (3, 2)), linewidth=0.65)
        ax.set_xlim(0, max_val); ax.set_ylim(0, max_val)
        ax.set_xlabel(str(data.get('xlabel') or 'Theoretical quantiles'))
        ax.set_ylabel(str(data.get('ylabel') or 'Observed statistics'))
        return self._finish(fig, axes, style, profile, spec, spec.title or 'Q-Q plot', left=0.17, bottom=0.21)

    def _trajectory(self, data, spec, container):
        time = _finite(data.get('time'))
        centers = np.asarray(data.get('centers'), dtype=float)
        if centers.ndim != 2 or centers.shape[1] != len(time):
            raise ValueError('trajectory centers 需要 cluster x time 矩阵')
        fig, axes, style, profile = self._canvas(spec, container, 1, 1, height=72)
        ax = axes.ravel()[0]
        labels = [str(v) for v in data.get('labels', [f'C{i + 1}' for i in range(len(centers))])]
        for index, (row, label) in enumerate(zip(centers, labels)):
            ax.plot(time, row, color=style.categorical_palette[index % len(style.categorical_palette)],
                    linewidth=1.05, marker='o', markersize=2.8, label=label)
        ax.set_xlabel(str(data.get('xlabel') or 'Time'))
        ax.set_ylabel(str(data.get('ylabel') or 'Z-score'))
        if len(labels) <= 8:
            ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
        else:
            fig.text(0.11, 0.02, f'{len(labels)} clusters; legend omitted for readability',
                     fontsize=max(5.5, profile.tick_font_pt - 0.5), color=style.muted_text)
        return self._finish(fig, axes, style, profile, spec, spec.title or 'Cluster trajectories',
                            left=0.14, right=0.78 if len(labels) <= 8 else 0.96, bottom=0.22)
