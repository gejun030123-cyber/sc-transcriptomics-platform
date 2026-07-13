"""Matplotlib rendering for legacy Plotly result payloads.

Analysis modules historically produced Plotly JSON.  The web result contract is
now static-only, so this adapter converts those payloads at the worker boundary
and keeps the rendering policy in one place instead of duplicating plot code in
every analysis module.
"""

import base64
import io
import json
import os
from typing import Any


PALETTE = [
    "#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
]


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

    try:
        from matplotlib import font_manager

        cjk_font = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
        if os.path.isfile(cjk_font):
            font_manager.fontManager.addfont(cjk_font)
    except Exception:
        pass
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Noto Sans CJK JP", "Noto Sans CJK SC", "DejaVu Sans"],
        "axes.unicode_minus": False,
    })
    payload = _decode(payload)
    layout = payload.get("layout", {}) or {}
    traces = payload.get("data", []) or []
    fig, ax = plt.subplots(figsize=(9, 6.5), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    rendered_any = False
    legend_items = []
    has_colorbar = False

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
                artist = ax.scatter(x, y, c=numeric_color, cmap="viridis", s=18,
                                    alpha=0.78, linewidths=0, rasterized=True,
                                    label=name or None)
                fig.colorbar(artist, ax=ax, fraction=0.035, pad=0.025, aspect=32)
                has_colorbar = True
            elif "markers" in mode:
                artist = ax.scatter(x, y, color=color, s=18, alpha=0.78,
                                    linewidths=0, rasterized=True,
                                    label=name or None)
            else:
                artist, = ax.plot(x, y, color=color, linewidth=1.8,
                                  alpha=0.9, label=name or None)
            rendered_any = True
            if name:
                legend_items.append(name)
        elif trace_type == "bar":
            numeric_y = _numeric(y)
            orientation = (trace.get("orientation") or "v").lower()
            if orientation == "h":
                ax.barh([str(item) for item in x], numeric_y, color=color, alpha=0.88,
                        label=name or None)
            else:
                positions = np.arange(len(x))
                ax.bar(positions, numeric_y, color=color, alpha=0.88, label=name or None)
                ax.set_xticks(positions, [str(item) for item in x], rotation=35,
                              ha="right")
            rendered_any = True
        elif trace_type in {"histogram"}:
            values = _numeric(x if x else y)
            if values is not None:
                ax.hist(values, bins=30, color=color, alpha=0.72, label=name or None)
                rendered_any = True
        elif trace_type in {"box", "violin"}:
            values = _numeric(y)
            if values is not None:
                if trace_type == "violin":
                    ax.violinplot(values, showmeans=True)
                else:
                    ax.boxplot(values, patch_artist=True,
                               boxprops={"facecolor": color, "alpha": 0.75})
                if name:
                    ax.set_xticks([1], [name])
                rendered_any = True
        elif trace_type in {"heatmap", "contour"}:
            matrix = np.asarray(trace.get("z", []), dtype=float)
            if matrix.size:
                image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r")
                fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025, aspect=32)
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
    ax.set_title(title, fontsize=15, fontweight="semibold", color="#172033", pad=14)
    if not has_colorbar:
        ax.set_xlabel(_axis_title(layout, "xaxis", ""), fontsize=11, color="#374151")
        ax.set_ylabel(_axis_title(layout, "yaxis", ""), fontsize=11, color="#374151")
    ax.tick_params(labelsize=9, colors="#4b5563", width=0.6)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if legend_items:
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
                  fontsize=9)
    fig.subplots_adjust(left=0.09, right=0.82 if legend_items or has_colorbar else 0.95,
                        bottom=0.12, top=0.88)
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
