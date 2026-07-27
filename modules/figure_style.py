"""Shared Nature-style figure contract for every project plotting backend.

Plotly is retained as an internal figure description format, but all persisted
project figures pass through this module before they become user-visible static
artifacts.  Keeping the contract here means new modules inherit the same style
by using the standard BaseAnalysis save helpers.
"""

NATURE_PALETTE = [
    '#0F4D92',  # deep blue
    '#42949E',  # teal
    '#7E6BA8',  # restrained violet
    '#B64342',  # rose/red accent
    '#D98C5F',  # warm secondary accent
    '#667085',  # neutral dark
    '#8FA9C9',  # soft blue
    '#6BAE75',  # muted green
    '#C7A7C9',  # soft lilac
    '#B7B9C8',  # neutral soft
]

# A restrained sequential scale for continuous measurements (QC metrics,
# expression summaries, pseudotime, etc.).  Keeping this beside the discrete
# palette prevents native Matplotlib figures and legacy Plotly conversions
# from silently falling back to Viridis.
NATURE_CONTINUOUS_COLORS = ['#DCEAF0', '#42949E', '#0F4D92', '#B64342']
NATURE_PLOTLY_CONTINUOUS_SCALE = [
    [0.0, NATURE_CONTINUOUS_COLORS[0]],
    [0.35, NATURE_CONTINUOUS_COLORS[1]],
    [0.7, NATURE_CONTINUOUS_COLORS[2]],
    [1.0, NATURE_CONTINUOUS_COLORS[3]],
]

NATURE_FONT_FAMILY = 'Arial'
NATURE_BG = 'white'
NATURE_TEXT = '#1F2937'
NATURE_MUTED = '#667085'
NATURE_AXIS = '#475467'
NATURE_GRID = '#E4E7EC'
NATURE_CJK_FONT = 'Noto Sans CJK JP'
NATURE_CJK_FONT_PATH = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'


def _font_for_text(text, preferred):
    """Use the installed CJK fallback only where a label needs it."""
    if text and any('\u3400' <= char <= '\u9fff' for char in str(text)):
        import os

        if os.path.isfile(NATURE_CJK_FONT_PATH):
            return NATURE_CJK_FONT
    return preferred


def nature_rcparams():
    """Return Matplotlib rcParams shared by native and converted figures."""
    return {
        'font.family': 'sans-serif',
        # Put the CJK-capable face first so Chinese result labels do not turn
        # into missing-glyph boxes; it remains visually close to Arial for
        # Latin text and falls back gracefully when the font is unavailable.
        'font.sans-serif': [NATURE_CJK_FONT, 'Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
        'font.size': 9,
        'svg.fonttype': 'none',
        'pdf.fonttype': 42,
        'axes.spines.right': False,
        'axes.spines.top': False,
        'axes.linewidth': 0.8,
        'axes.edgecolor': NATURE_AXIS,
        'axes.labelcolor': NATURE_TEXT,
        'xtick.color': NATURE_AXIS,
        'ytick.color': NATURE_AXIS,
        'legend.frameon': False,
        'figure.facecolor': NATURE_BG,
        'savefig.facecolor': NATURE_BG,
    }


def nature_continuous_cmap():
    """Return the shared Matplotlib sequential colormap."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        'nature_continuous', list(NATURE_CONTINUOUS_COLORS)
    )


def apply_matplotlib_style(fig, viz=None):
    """Normalize a Matplotlib figure before PNG/SVG persistence."""
    import os
    import matplotlib as mpl

    try:
        from matplotlib import font_manager

        if os.path.isfile(NATURE_CJK_FONT_PATH):
            font_manager.fontManager.addfont(NATURE_CJK_FONT_PATH)
    except Exception:
        pass

    viz = viz or {}
    mpl.rcParams.update(nature_rcparams())
    font_size = max(8.0, float(viz.get('font_size', 9)))
    font_family = viz.get('font_family', NATURE_FONT_FAMILY)
    bg_color = viz.get('bg_color', NATURE_BG)
    for axis in getattr(fig, 'axes', []):
        axis.set_facecolor(bg_color)
        axis.tick_params(labelsize=max(7, font_size - 1), width=0.7,
                         length=3, colors=NATURE_AXIS)
        axis.xaxis.label.set_size(font_size)
        axis.yaxis.label.set_size(font_size)
        axis.xaxis.label.set_color(NATURE_TEXT)
        axis.yaxis.label.set_color(NATURE_TEXT)
        axis.xaxis.label.set_family(_font_for_text(axis.get_xlabel(), font_family))
        axis.yaxis.label.set_family(_font_for_text(axis.get_ylabel(), font_family))
        title = axis.title
        title.set_fontsize(font_size + 1)
        title.set_fontweight('semibold')
        title.set_color(NATURE_TEXT)
        title.set_family(_font_for_text(title.get_text(), font_family))
        for spine_name, spine in axis.spines.items():
            spine.set_linewidth(0.7)
            spine.set_color('#98A2B3')
            spine.set_visible(spine_name in ('left', 'bottom'))
        axis.grid(False)
        legend = axis.get_legend()
        if legend is not None:
            legend.set_frame_on(False)
            for text in legend.get_texts():
                text.set_fontsize(max(7, font_size - 1))
                text.set_family(_font_for_text(text.get_text(), font_family))
    fig.patch.set_facecolor(bg_color)
    return fig


def style_plotly_figure(fig, viz=None):
    """Apply the shared Plotly-side contract before JSON persistence."""
    viz = viz or {}
    width = int(viz.get('figure_width', 900))
    height = int(viz.get('figure_height', 560))
    font_family = viz.get('font_family', NATURE_FONT_FAMILY)
    font_size = max(9, int(float(viz.get('font_size', 11))))
    bg_color = viz.get('bg_color', NATURE_BG)
    fig.update_layout(
        template='simple_white',
        width=width,
        height=height,
        paper_bgcolor=bg_color,
        plot_bgcolor=bg_color,
        colorway=list(NATURE_PALETTE),
        font={'family': font_family, 'size': font_size, 'color': NATURE_TEXT},
        title_font={'family': font_family, 'size': font_size + 3,
                    'color': NATURE_TEXT},
        margin={'l': 72, 'r': 34, 't': 68, 'b': 64},
        legend={'bgcolor': 'rgba(255,255,255,0)', 'borderwidth': 0,
                'font': {'family': font_family, 'size': max(9, font_size - 1)}},
    )
    fig.update_xaxes(
        showline=True, linewidth=0.8, linecolor='#98A2B3', mirror=False,
        showgrid=False, zeroline=False,
        title_font={'family': font_family, 'size': font_size, 'color': NATURE_TEXT},
        tickfont={'family': font_family, 'size': max(8, font_size - 1), 'color': NATURE_AXIS},
    )
    fig.update_yaxes(
        showline=True, linewidth=0.8, linecolor='#98A2B3', mirror=False,
        showgrid=True, gridcolor=NATURE_GRID, gridwidth=0.6, zeroline=False,
        title_font={'family': font_family, 'size': font_size, 'color': NATURE_TEXT},
        tickfont={'family': font_family, 'size': max(8, font_size - 1), 'color': NATURE_AXIS},
    )
    fig.update_layout(coloraxis={'colorscale': list(NATURE_PLOTLY_CONTINUOUS_SCALE)})
    # Modules may provide an explicit Viridis/Plasma string.  Normalize those
    # trace-level scales as well, so the interactive JSON and the static
    # renderer cannot drift apart.
    for trace in getattr(fig, 'data', []):
        for channel in ('marker', 'line'):
            style = getattr(trace, channel, None)
            if style is not None and getattr(style, 'colorscale', None) is not None:
                style.colorscale = list(NATURE_PLOTLY_CONTINUOUS_SCALE)
        if getattr(trace, 'colorscale', None) is not None:
            trace.colorscale = list(NATURE_PLOTLY_CONTINUOUS_SCALE)
    return fig
