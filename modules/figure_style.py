"""Shared Nature-style figure contract for every project plotting backend.

Plotly is retained as an internal figure description format, but all persisted
project figures pass through this module before they become user-visible static
artifacts.  Keeping the contract here means new modules inherit the same style
by using the standard BaseAnalysis save helpers.
"""

import math

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

VISUALIZATION_THEME_DEFAULTS = {
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
        'bg_color': '#1A1A2E',
        'color_palette': ['#E94560', '#00B4D8', '#90E0EF', '#F8F9FA',
                          '#F6BD60', '#6BAE75', '#C7A7C9', '#B7B9C8'],
        'font_family': 'Arial',
    },
}


def _as_bool(value, default=False):
    """Coerce browser/form values without treating the string ``'false'`` as true."""
    if isinstance(value, str):
        if value.strip().lower() in {'false', '0', 'no', 'off', ''}:
            return False
        if value.strip().lower() in {'true', '1', 'yes', 'on'}:
            return True
    return default if value is None else bool(value)


def _bounded_number(value, default, lower, upper, integer=False):
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = float(default)
    if not math.isfinite(number):
        number = float(default)
    number = max(lower, min(upper, number))
    return int(number) if integer else number


def normalize_visualization_params(viz=None):
    """Return bounded, shared visualization settings for every figure backend."""
    raw = dict(viz or {})
    theme_name = raw.get('theme', 'default')
    theme = VISUALIZATION_THEME_DEFAULTS.get(
        theme_name, VISUALIZATION_THEME_DEFAULTS['default']
    )
    params = dict(raw)
    params['theme'] = theme_name if theme_name in VISUALIZATION_THEME_DEFAULTS else 'default'
    params['bg_color'] = raw.get('bg_color') or theme['bg_color']
    raw_palette = raw.get('color_palette')
    params['color_palette'] = ([raw_palette] if isinstance(raw_palette, str)
                               else list(raw_palette or theme['color_palette']))
    params['font_family'] = raw.get('font_family') or theme['font_family']
    params['font_size'] = _bounded_number(raw.get('font_size'), 12, 8, 28, integer=True)
    params['figure_width'] = _bounded_number(raw.get('figure_width'), 800, 320, 2400, integer=True)
    params['figure_height'] = _bounded_number(raw.get('figure_height'), 500, 240, 1800, integer=True)
    params['umap_point_size'] = _bounded_number(raw.get('umap_point_size'), 5, 1, 20)
    params['umap_opacity'] = _bounded_number(raw.get('umap_opacity'), 0.7, 0.1, 1.0)
    params['umap_legend_fontsize'] = _bounded_number(raw.get('umap_legend_fontsize'), 10, 7, 24, integer=True)
    params['umap_label_categories'] = _as_bool(raw.get('umap_label_categories'), False)
    params['umap_hide_axes'] = _as_bool(raw.get('umap_hide_axes'), True)
    mode = str(raw.get('figure_mode', 'publication')).strip().lower()
    params['figure_mode'] = mode if mode in {
        'standard', 'publication', 'nature_portfolio'
    } else 'publication'
    width_profile = str(raw.get('width_profile', 'single')).strip().lower()
    params['width_profile'] = width_profile if width_profile in {'single', 'double'} else 'single'
    nature_style = str(raw.get('nature_style', 'nature')).strip().lower()
    params['nature_style'] = nature_style if nature_style in {
        'nature', 'nature_communications', 'nature_aging'
    } else 'nature'
    raw_formats = raw.get('static_formats', raw.get('export_formats', ('png', 'svg')))
    if isinstance(raw_formats, str):
        raw_formats = raw_formats.replace(',', ' ').split()
    params['static_formats'] = tuple(dict.fromkeys(
        str(value).strip().lower() for value in (raw_formats or ())
        if str(value).strip().lower() in {'png', 'svg', 'pdf', 'tiff'}
    )) or ('png', 'svg')
    params['export_formats'] = params['static_formats']
    if params['figure_mode'] == 'nature_portfolio':
        # Publication visual constants are owned by Figure Engine, not the UI.
        params['theme'] = 'nature'
        params['bg_color'] = 'white'
    return params


def stable_category_colors(categories, existing=None, palette=None):
    """Assign repeatable colors while honoring AnnData's category color metadata."""
    labels = [str(category) for category in categories]
    supplied = list(existing or [])
    colors = list(palette or NATURE_PALETTE)
    if len(supplied) >= len(labels):
        return dict(zip(labels, supplied[:len(labels)]))
    return {
        label: colors[index % len(colors)]
        for index, label in enumerate(labels)
    }


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


def register_nature_cjk_font():
    """Register the local CJK face once and return whether it is available.

    Matplotlib's ``addfont`` appends to a process-global font list. Repeating
    it for the same TTC can change generic font resolution and make otherwise
    deterministic canvases depend on render order.
    """
    import os

    if not os.path.isfile(NATURE_CJK_FONT_PATH):
        return False
    try:
        from matplotlib import font_manager

        expected = os.path.realpath(NATURE_CJK_FONT_PATH)
        registered = any(
            os.path.realpath(getattr(entry, 'fname', '')) == expected
            for entry in font_manager.fontManager.ttflist
        )
        if not registered:
            font_manager.fontManager.addfont(NATURE_CJK_FONT_PATH)
        return True
    except Exception:
        return False


def apply_matplotlib_style(fig, viz=None):
    """Normalize a Matplotlib figure before PNG/SVG persistence."""
    import matplotlib as mpl

    register_nature_cjk_font()

    viz = normalize_visualization_params(viz)
    mpl.rcParams.update(nature_rcparams())
    font_size = max(8.0, float(viz.get('font_size', 9)))
    font_family = viz.get('font_family', NATURE_FONT_FAMILY)
    bg_color = viz.get('bg_color', NATURE_BG)
    is_dark = str(bg_color).lower() in {'#1a1a2e', '#111827', '#0f172a'}
    text_color = '#F8FAFC' if is_dark else NATURE_TEXT
    axis_color = '#CBD5E1' if is_dark else NATURE_AXIS
    spine_color = '#64748B' if is_dark else '#98A2B3'
    for axis in getattr(fig, 'axes', []):
        axis.set_facecolor(bg_color)
        axis.tick_params(labelsize=max(7, font_size - 1), width=0.7,
                         length=3, colors=axis_color)
        axis.xaxis.label.set_size(font_size)
        axis.yaxis.label.set_size(font_size)
        axis.xaxis.label.set_color(text_color)
        axis.yaxis.label.set_color(text_color)
        axis.xaxis.label.set_family(_font_for_text(axis.get_xlabel(), font_family))
        axis.yaxis.label.set_family(_font_for_text(axis.get_ylabel(), font_family))
        title = axis.title
        title.set_fontsize(font_size + 1)
        title.set_fontweight('semibold')
        title.set_color(text_color)
        title.set_family(_font_for_text(title.get_text(), font_family))
        for spine_name, spine in axis.spines.items():
            spine.set_linewidth(0.7)
            spine.set_color(spine_color)
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
    viz = normalize_visualization_params(viz)
    width = int(viz.get('figure_width', 900))
    height = int(viz.get('figure_height', 560))
    font_family = viz.get('font_family', NATURE_FONT_FAMILY)
    font_size = max(9, int(float(viz.get('font_size', 11))))
    bg_color = viz.get('bg_color', NATURE_BG)
    palette = list(viz.get('color_palette') or NATURE_PALETTE)
    is_dark = str(bg_color).lower() in {'#1a1a2e', '#111827', '#0f172a'}
    text_color = '#F8FAFC' if is_dark else NATURE_TEXT
    axis_color = '#CBD5E1' if is_dark else NATURE_AXIS
    grid_color = '#334155' if is_dark else NATURE_GRID
    fig.update_layout(
        template='plotly_dark' if is_dark else 'simple_white',
        width=width,
        height=height,
        paper_bgcolor=bg_color,
        plot_bgcolor=bg_color,
        colorway=palette,
        font={'family': font_family, 'size': font_size, 'color': text_color},
        title_font={'family': font_family, 'size': font_size + 3,
                    'color': text_color},
        margin={'l': 72, 'r': 34, 't': 68, 'b': 64},
        legend={'bgcolor': 'rgba(255,255,255,0)', 'borderwidth': 0,
                'font': {'family': font_family, 'size': max(9, font_size - 1),
                         'color': text_color}},
    )
    fig.update_xaxes(
        showline=True, linewidth=0.8, linecolor=axis_color, mirror=False,
        showgrid=False, zeroline=False,
        title_font={'family': font_family, 'size': font_size, 'color': text_color},
        tickfont={'family': font_family, 'size': max(8, font_size - 1), 'color': axis_color},
    )
    fig.update_yaxes(
        showline=True, linewidth=0.8, linecolor=axis_color, mirror=False,
        showgrid=not is_dark, gridcolor=grid_color, gridwidth=0.6, zeroline=False,
        title_font={'family': font_family, 'size': font_size, 'color': text_color},
        tickfont={'family': font_family, 'size': max(8, font_size - 1), 'color': axis_color},
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
