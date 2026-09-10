"""Nature Portfolio 的真实尺寸、字体和语义色板。"""

from dataclasses import dataclass
from typing import Dict, Tuple

from figure_engine.spec import FigureSpec


MM_PER_INCH = 25.4


@dataclass(frozen=True)
class StyleProfile:
    width_mm: float
    height_mm: float
    base_font_pt: float
    title_font_pt: float
    axis_font_pt: float
    tick_font_pt: float
    legend_font_pt: float
    axis_linewidth_pt: float
    tick_width_pt: float
    tick_length_pt: float
    marker_size_pt2: float
    dpi: int

    @property
    def figsize(self) -> Tuple[float, float]:
        return self.width_mm / MM_PER_INCH, self.height_mm / MM_PER_INCH


class NatureStyle:
    """Nature Portfolio 基础风格；模板只能读取，不应在运行时修改。"""

    name = 'nature'
    background = '#FFFFFF'
    text = '#20262E'
    muted_text = '#626C78'
    axis = '#4A5563'
    subtle = '#CDD3DA'
    neutral = '#D8DCE2'
    neutral_dark = '#77808C'
    signal_blue = '#4C78A8'
    signal_teal = '#4C9099'
    signal_violet = '#7A6FA8'
    signal_red = '#B65C5C'
    signal_orange = '#D08A5B'

    deg_palette = {
        'NS': '#D8DCE2',
        'Up': '#B65C5C',
        'Down': '#4C78A8',
    }
    categorical_palette = (
        '#4C78A8', '#D08A5B', '#6F9D72', '#8A73A8', '#4C9099',
        '#B65C5C', '#8B8278', '#7C8FB3', '#B88AA7', '#8096A0',
        '#A98262', '#5F8B62', '#A66B7D', '#6A8298', '#B59A5B',
        '#6F7886', '#9A7D68', '#78908B', '#8C6F89', '#A86F62',
    )
    marker_cycle = ('o', 's', '^', 'D', 'P', 'X', 'v', '<', '>')
    heatmap_colors = ('#355F8A', '#90AEC5', '#E7EEF2', '#FFFFFF',
                      '#F1E4DF', '#D59B91', '#A94F50')
    fdr_colors = ('#183B66', '#356B91', '#73A5B5', '#C3D9D7', '#EEF1EC')

    _single_heights: Dict[str, float] = {
        'pca': 78.0,
        'embedding': 78.0,
        'volcano': 79.0,
        'heatmap': 116.0,
        'gsea': 102.0,
        'enrichment': 102.0,
        'enrichment_dotplot': 102.0,
        'enrichment_overview': 126.0,
        'go_priority_dotplot': 150.0,
        'go_priority_barplot': 150.0,
        'pathway_triptych_dotplot': 150.0,
        'pathway_triptych_barplot': 150.0,
        'enrichment_barplot': 102.0,
        'enrichment_chord': 112.0,
        'enrichment_cnetplot': 108.0,
        'enrichment_emapplot': 108.0,
        'gsea_running': 116.0,
        'ma': 79.0,
        'correlation_heatmap': 104.0,
        'gsva': 112.0,
        'ssgsea': 112.0,
        'upset': 105.0,
        'wgcna': 108.0,
        'diagnostic': 82.0,
        'composite': 120.0,
    }
    _double_heights: Dict[str, float] = {
        'pca': 108.0,
        'embedding': 108.0,
        'volcano': 112.0,
        'heatmap': 128.0,
        'gsea': 116.0,
        'enrichment': 116.0,
        'enrichment_dotplot': 116.0,
        'enrichment_overview': 132.0,
        'go_priority_dotplot': 164.0,
        'go_priority_barplot': 164.0,
        'pathway_triptych_dotplot': 164.0,
        'pathway_triptych_barplot': 164.0,
        'enrichment_barplot': 116.0,
        'enrichment_chord': 156.0,
        'enrichment_cnetplot': 128.0,
        'enrichment_emapplot': 128.0,
        'gsea_running': 128.0,
        'ma': 112.0,
        'correlation_heatmap': 128.0,
        'gsva': 128.0,
        'ssgsea': 128.0,
        'upset': 122.0,
        'wgcna': 128.0,
        'diagnostic': 110.0,
        'composite': 160.0,
    }

    def profile(self, spec: FigureSpec) -> StyleProfile:
        heights = self._single_heights if spec.width == 'single' else self._double_heights
        height_mm = float(spec.height_mm or heights.get(spec.plot_type, 90.0))
        if spec.mode == 'standard':
            base, title, axis, tick, legend, dpi = 8.5, 10.0, 8.5, 7.5, 7.5, min(spec.dpi, 300)
            marker = 30.0 if spec.width == 'single' else 38.0
        elif spec.mode == 'publication':
            base, title, axis, tick, legend, dpi = 7.5, 8.5, 7.5, 6.5, 6.5, min(spec.dpi, 600)
            marker = 24.0 if spec.width == 'single' else 30.0
        else:
            base, title, axis, tick, legend, dpi = 7.0, 8.0, 7.0, 6.2, 6.2, max(300, spec.dpi)
            marker = 21.0 if spec.width == 'single' else 28.0
        return StyleProfile(
            width_mm=spec.width_mm,
            height_mm=height_mm,
            base_font_pt=base,
            title_font_pt=title,
            axis_font_pt=axis,
            tick_font_pt=tick,
            legend_font_pt=legend,
            axis_linewidth_pt=0.65,
            tick_width_pt=0.55,
            tick_length_pt=2.4,
            marker_size_pt2=marker,
            dpi=dpi,
        )

    def rc_params(self, spec: FigureSpec) -> dict:
        profile = self.profile(spec)
        return {
            # Use an explicit CJK-capable face instead of the generic
            # ``sans-serif`` alias. Registration is idempotent in the shared
            # platform style helper, so this remains stable across render order.
            'font.family': 'Noto Sans CJK JP',
            'font.sans-serif': [
                'Noto Sans CJK JP', 'Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif',
            ],
            'font.size': profile.base_font_pt,
            'axes.titlesize': profile.title_font_pt,
            'axes.labelsize': profile.axis_font_pt,
            'xtick.labelsize': profile.tick_font_pt,
            'ytick.labelsize': profile.tick_font_pt,
            'legend.fontsize': profile.legend_font_pt,
            'axes.linewidth': profile.axis_linewidth_pt,
            'axes.edgecolor': self.axis,
            'axes.labelcolor': self.text,
            'axes.titlecolor': self.text,
            'xtick.color': self.axis,
            'ytick.color': self.axis,
            'xtick.major.width': profile.tick_width_pt,
            'ytick.major.width': profile.tick_width_pt,
            'xtick.major.size': profile.tick_length_pt,
            'ytick.major.size': profile.tick_length_pt,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.grid': False,
            'legend.frameon': False,
            'figure.facecolor': self.background,
            'axes.facecolor': self.background,
            'savefig.facecolor': self.background,
            'svg.fonttype': 'none',
            'pdf.fonttype': 42,
            'ps.fonttype': 42,
            'axes.unicode_minus': False,
        }

    def context(self, spec: FigureSpec):
        import matplotlib as mpl

        return mpl.rc_context(self.rc_params(spec))

    def apply_axis(self, ax, profile: StyleProfile, *, despine=True):
        ax.set_facecolor(self.background)
        ax.grid(False)
        ax.tick_params(
            colors=self.axis,
            labelsize=profile.tick_font_pt,
            width=profile.tick_width_pt,
            length=profile.tick_length_pt,
        )
        ax.xaxis.label.set_color(self.text)
        ax.yaxis.label.set_color(self.text)
        ax.xaxis.label.set_size(profile.axis_font_pt)
        ax.yaxis.label.set_size(profile.axis_font_pt)
        ax.title.set_color(self.text)
        ax.title.set_fontsize(profile.title_font_pt)
        ax.title.set_fontweight('semibold')
        for name, spine in ax.spines.items():
            spine.set_linewidth(profile.axis_linewidth_pt)
            spine.set_color(self.axis)
            if despine:
                spine.set_visible(name in {'left', 'bottom'})

    def group_colors(self, values):
        labels = list(dict.fromkeys(str(value) for value in values))
        return {
            label: self.categorical_palette[index % len(self.categorical_palette)]
            for index, label in enumerate(labels)
        }

    def batch_markers(self, values):
        labels = list(dict.fromkeys(str(value) for value in values))
        return {
            label: self.marker_cycle[index % len(self.marker_cycle)]
            for index, label in enumerate(labels)
        }

    def heatmap_cmap(self):
        from matplotlib.colors import LinearSegmentedColormap

        cmap = LinearSegmentedColormap.from_list('nature_muted_diverging', self.heatmap_colors)
        cmap.set_bad('#F2F3F5')
        return cmap

    def fdr_cmap(self):
        from matplotlib.colors import LinearSegmentedColormap

        return LinearSegmentedColormap.from_list('nature_fdr', self.fdr_colors)

    def expression_cmap(self):
        """Sequential map for expression, QC and pseudotime embeddings."""
        from matplotlib.colors import LinearSegmentedColormap

        cmap = LinearSegmentedColormap.from_list(
            'nature_expression',
            ('#E7EFF2', '#9BC4C9', self.signal_teal, self.signal_blue, self.signal_red),
        )
        cmap.set_bad('#D8DCE2')
        return cmap
