"""Matplotlib rendering for legacy Plotly result payloads.

Analysis modules historically produced Plotly JSON.  The web result contract is
now static-only, so this adapter converts those payloads at the worker boundary
and keeps the rendering policy in one place instead of duplicating plot code in
every analysis module.
"""

import base64
import io
import json
import math
import os
from typing import Any


from modules.figure_style import (
    NATURE_AXIS,
    NATURE_BG,
    NATURE_FONT_FAMILY,
    NATURE_GRID,
    NATURE_MUTED,
    NATURE_PALETTE,
    nature_continuous_cmap,
    NATURE_TEXT,
    apply_matplotlib_style,
    register_nature_cjk_font,
)

PALETTE = list(NATURE_PALETTE)


def _nature_cmap():
    return nature_continuous_cmap()


def _decode(value: Any):
    """Decode Plotly's compact dtype/bdata arrays into ordinary values."""
    if isinstance(value, dict) and "bdata" in value and "dtype" in value:
        import numpy as np

        raw = base64.b64decode(value["bdata"])
        return np.frombuffer(raw, dtype=np.dtype(value["dtype"])).tolist()
    if isinstance(value, dict):
        return {key: _decode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_decode(item) for item in value]
    return value


def _title(layout, fallback):
    value = layout.get("title", fallback)
    if isinstance(value, dict):
        value = value.get("text", fallback)
    return str(value or fallback)


def _axis_title(layout, axis, fallback):
    value = layout.get(axis, {}) or {}
    title = value.get("title", fallback) if isinstance(value, dict) else fallback
    if isinstance(title, dict):
        title = title.get("text", fallback)
    return str(title or fallback)


def _numeric(values):
    import numpy as np

    if values is None:
        return None
    try:
        return np.asarray(values, dtype=float)
    except (TypeError, ValueError):
        return None


def _sample_tick_labels(labels, max_labels=18):
    """Keep the first/last category while thinning dense axis labels."""
    labels = [str(value) for value in labels]
    if len(labels) <= max_labels:
        indices = list(range(len(labels)))
    else:
        stride = max(1, int(math.ceil(len(labels) / max_labels)))
        indices = list(range(0, len(labels), stride))
        if indices[-1] != len(labels) - 1:
            indices.append(len(labels) - 1)
    return indices, [labels[index] for index in indices]


def _trace_color(trace, index):
    marker = trace.get("marker", {}) or {}
    color = marker.get("color", trace.get("line", {}).get("color"))
    if isinstance(color, list):
        return color
    return color or PALETTE[index % len(PALETTE)]


def _render_image_trace(ax, trace):
    source = trace.get("source", "")
    if not isinstance(source, str) or "," not in source:
        return False
    try:
        from PIL import Image

        payload = source.split(",", 1)[1]
        image = Image.open(io.BytesIO(base64.b64decode(payload)))
        ax.imshow(image, aspect="auto")
        ax.set_axis_off()
        return True
    except Exception:
        return False


def _render_sankey(ax, trace):
    """Render a readable static Sankey summary when vector flows are unavailable."""
    import numpy as np

    node = trace.get("node", {}) or {}
    labels = [str(item) for item in node.get("label", [])]
    values = trace.get("link", {}).get("value", []) or []
    if not labels:
        ax.text(0.5, 0.5, "Sankey 数据为空", ha="center", va="center")
        return
    count = len(labels)
    y = np.linspace(0.95, 0.05, count)
    ax.barh(y, [1] * count, height=max(0.012, 0.7 / count), color=PALETTE[0], alpha=0.85)
    for position, label in zip(y, labels):
        ax.text(1.02, position, label, va="center", fontsize=8)
    ax.set_xlim(0, 1.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.text(0.02, 0.98, f"{len(values)} 条流向", transform=ax.transAxes,
            va="top", fontsize=9, color="#4b5563")


def render_plotly_payload(payload: dict, output_stem: str, label: str = ""):
    """Render a Plotly JSON payload as publication-style PNG and SVG."""
    import matplotlib.pyplot as plt
    import numpy as np

    register_nature_cjk_font()
    from matplotlib import rcParams

    rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Noto Sans CJK JP", "Noto Sans CJK SC",
                             "Arial", "Helvetica", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })
    payload = _decode(payload)
    layout = payload.get("layout", {}) or {}
    traces = payload.get("data", []) or []
    width = min(10.5, max(7.0, float(layout.get('width', 900)) / 100))
    height = min(8.0, max(4.8, float(layout.get('height', 560)) / 100))
    # Legacy Plotly payloads do not carry the content-aware canvas sizing used
    # by native figures.  Reserve additional horizontal/vertical space before
    # drawing categorical axes so long cluster/sample labels do not collide.
    vertical_labels = [
        str(value)
        for trace in traces
        if (trace or {}).get('type', 'scatter') == 'bar'
        and (trace or {}).get('orientation', 'v') != 'h'
        for value in ((trace or {}).get('x', []) or [])
    ]
    horizontal_labels = [
        str(value)
        for trace in traces
        if (trace or {}).get('type', 'scatter') == 'bar'
        and (trace or {}).get('orientation', 'v') == 'h'
        for value in ((trace or {}).get('x', []) or [])
    ]
    if vertical_labels:
        width = max(width, min(16.0, 6.8 + 0.42 * len(set(vertical_labels))))
        if max(map(len, vertical_labels), default=0) > 16:
            height = max(height, 5.8)
    if horizontal_labels:
        height = max(height, min(16.0, 4.8 + 0.28 * len(set(horizontal_labels))))
    fig, ax = plt.subplots(figsize=(width, height), dpi=150)
    fig.patch.set_facecolor(NATURE_BG)
    ax.set_facecolor(NATURE_BG)
    rendered_any = False
    legend_items = []
    has_colorbar = False
    violin_positions = []
    violin_labels = []
    horizontal_tick_labels = []
    bar_traces = [trace for trace in traces
                  if (trace or {}).get('type', 'scatter') == 'bar'
                  and (trace or {}).get('orientation', 'v') != 'h']
    bar_index = 0
    bar_tick_labels = []

    for index, raw_trace in enumerate(traces):
        trace = raw_trace or {}
        trace_type = trace.get("type", "scatter")
        name = str(trace.get("name", ""))
        if trace_type == "image":
            rendered_any = _render_image_trace(ax, trace) or rendered_any
            continue
        if trace_type == "sankey":
            _render_sankey(ax, trace)
            rendered_any = True
            continue

        x = trace.get("x", []) or []
        y = trace.get("y", []) or []
        color = _trace_color(trace, index)
        if trace_type in {"scatter", "scattergl"}:
            mode = str(trace.get("mode", "lines"))
            numeric_color = _numeric(color)
            if numeric_color is not None and len(numeric_color) == len(x):
                artist = ax.scatter(x, y, c=numeric_color, cmap=_nature_cmap(), s=20,
                                    alpha=0.68, linewidths=0, rasterized=True,
                                    label=name or None)
                colorbar = fig.colorbar(artist, ax=ax, fraction=0.035, pad=0.025, aspect=32)
                colorbar.outline.set_visible(False)
                colorbar.ax.tick_params(labelsize=8, width=0.6, colors=NATURE_AXIS)
                has_colorbar = True
            elif "markers" in mode:
                artist = ax.scatter(x, y, color=color, s=20, alpha=0.68,
                                    linewidths=0, rasterized=True,
                                    label=name or None)
            else:
                artist, = ax.plot(x, y, color=color, linewidth=1.6,
                                  alpha=0.9, label=name or None)
            rendered_any = True
            if name:
                legend_items.append(name)
        elif trace_type == "bar":
            numeric_y = _numeric(y)
            orientation = (trace.get("orientation") or "v").lower()
            if orientation == "h":
                horizontal_tick_labels = [str(item) for item in x]
                ax.barh(horizontal_tick_labels, numeric_y, color=color, alpha=0.84,
                        label=name or None)
            else:
                positions = np.arange(len(x))
                n_series = max(1, len(bar_traces))
                bar_width = min(0.78 / n_series, 0.32)
                offset = (bar_index - (n_series - 1) / 2) * bar_width
                ax.bar(positions + offset, numeric_y, width=bar_width,
                       color=color, alpha=0.86, edgecolor='white', linewidth=0.35,
                       label=name or None)
                bar_tick_labels = [str(item) for item in x]
                bar_index += 1
                rendered_any = True
        elif trace_type in {"histogram"}:
            values = _numeric(x if x else y)
            if values is not None:
                bins = int(trace.get('nbinsx') or trace.get('nbinsy') or 30)
                ax.hist(values, bins=max(10, min(80, bins)), color=color,
                        alpha=0.72, edgecolor='white', linewidth=0.3,
                        label=name or None)
                rendered_any = True
        elif trace_type in {"box", "violin"}:
            values = _numeric(y)
            if values is not None:
                values = values[np.isfinite(values)]
                if not len(values):
                    continue
                if trace_type == "violin":
                    position = len(violin_positions) + 1
                    parts = ax.violinplot(
                        values, positions=[position], widths=0.78,
                        showmeans=False, showmedians=True, showextrema=False,
                    )
                    for body in parts.get('bodies', []):
                        body.set_facecolor(color)
                        body.set_edgecolor(color)
                        body.set_alpha(0.62)
                        body.set_linewidth(0.7)
                    if 'cmedians' in parts:
                        parts['cmedians'].set_color(NATURE_TEXT)
                        parts['cmedians'].set_linewidth(1.1)
                    violin_positions.append(position)
                    violin_labels.append(name or f'Metric {position}')
                else:
                    position = len(violin_positions) + 1
                    box = ax.boxplot(values, positions=[position], widths=0.58,
                                     patch_artist=True, showfliers=False,
                                     boxprops={"facecolor": color, "alpha": 0.72,
                                               'edgecolor': color},
                                     medianprops={'color': NATURE_TEXT, 'linewidth': 1.1},
                                     whiskerprops={'color': color},
                                     capprops={'color': color})
                    violin_positions.append(position)
                    violin_labels.append(name or f'Metric {position}')
                rendered_any = True
        elif trace_type in {"heatmap", "contour"}:
            matrix = np.asarray(trace.get("z", []), dtype=float)
            if matrix.size:
                image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r")
                colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025, aspect=32)
                colorbar.outline.set_visible(False)
                colorbar.ax.tick_params(labelsize=8, width=0.6, colors=NATURE_AXIS)
                has_colorbar = True
                rendered_any = True
        elif trace_type == "pie":
            labels = [str(item) for item in trace.get("labels", [])]
            values = _numeric(trace.get("values", []))
            if values is not None and values.size:
                ax.pie(values, labels=labels, autopct="%1.1f%%", startangle=90,
                       colors=PALETTE[:len(values)])
                rendered_any = True
        elif trace_type == "table":
            header = trace.get("header", {}).get("values", [])
            cells = trace.get("cells", {}).get("values", [])
            table = ax.table(cellText=list(zip(*cells)) if cells else [],
                             colLabels=header, loc="center", cellLoc="center")
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            ax.set_axis_off()
            rendered_any = True

    if not rendered_any:
        ax.text(0.5, 0.5, "该图暂无可用静态数据", ha="center", va="center",
                fontsize=13, color="#6b7280")
    title = _title(layout, label or "分析图")
    if violin_positions:
        tick_cap = 12 if max(map(len, violin_labels), default=0) > 16 else 18
        selected, shown = _sample_tick_labels(violin_labels, max_labels=tick_cap)
        ax.set_xticks([violin_positions[index] for index in selected], shown)
        ax.tick_params(axis='x', labelrotation=28, labelsize=8)
    elif bar_tick_labels:
        tick_cap = 12 if max(map(len, bar_tick_labels), default=0) > 16 else 18
        selected, shown = _sample_tick_labels(bar_tick_labels, max_labels=tick_cap)
        centers = np.arange(len(bar_tick_labels))
        ax.set_xticks(centers[selected], shown, rotation=30, ha='right')
    if horizontal_tick_labels:
        tick_cap = 12 if max(map(len, horizontal_tick_labels), default=0) > 16 else 18
        selected, shown = _sample_tick_labels(horizontal_tick_labels, max_labels=tick_cap)
        ax.set_yticks(selected, shown)
        ax.tick_params(axis='y', labelsize=8)
    ax.set_title(title, fontsize=14, fontweight="semibold", color=NATURE_TEXT,
                 loc='left', pad=12)
    if not has_colorbar:
        ax.set_xlabel(_axis_title(layout, "xaxis", ""), fontsize=11, color="#374151")
        ax.set_ylabel(_axis_title(layout, "yaxis", ""), fontsize=11, color="#374151")
    ax.tick_params(labelsize=9, colors=NATURE_AXIS, width=0.6)
    ax.grid(axis='y', color=NATURE_GRID, linewidth=0.55, alpha=0.72)
    ax.set_axisbelow(True)
    apply_matplotlib_style(fig, {'font_size': 9, 'font_family': NATURE_FONT_FAMILY})
    if legend_items:
        legend = ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
                           frameon=False, fontsize=8)
        legend.set_title('Group', prop={'size': 8})
    fig.subplots_adjust(left=0.10, right=0.82 if legend_items or has_colorbar else 0.95,
                        bottom=0.16, top=0.88)
    png_path = f"{output_stem}.png"
    svg_path = f"{output_stem}.svg"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.15,
                facecolor="white")
    fig.savefig(svg_path, format="svg", bbox_inches="tight", pad_inches=0.15,
                facecolor="white")
    plt.close(fig)
    return [
        {"file_path": png_path, "file_type": "png", "label": label},
        {"file_path": svg_path, "file_type": "svg", "label": label},
    ]


def render_plotly_json(path: str, label: str = ""):
    """Load a Plotly JSON file and write static siblings beside it."""
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    stem, _ = os.path.splitext(path)
    return render_plotly_payload(payload, stem, label=label)
