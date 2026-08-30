"""Nature Figure Validator：在最终物理尺寸上执行自动 QA。"""

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Optional

import numpy as np

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from figure_engine.accessibility import audit_palette


@dataclass
class ValidationIssue:
    code: str
    severity: str
    message: str
    deduction: int
    details: dict = field(default_factory=dict)


@dataclass
class ValidationReport:
    score: int
    status: str
    ready: bool
    issues: list
    metrics: dict
    spec: dict

    def to_dict(self):
        return {
            # Keep the historical field for consumers that already read it,
            # while making the scope of this score explicit.
            'figure_quality_score': self.score,
            'nature_readiness_score': self.score,
            'score_type': 'figure_quality',
            'status': self.status,
            'ready': self.ready,
            'issues': [asdict(issue) if isinstance(issue, ValidationIssue) else dict(issue)
                       for issue in self.issues],
            'metrics': self.metrics,
            'figure_spec': self.spec,
        }

    def write_json(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=2), encoding='utf-8')
        return str(path)


class FigureReadinessError(RuntimeError):
    """Raised only when a caller explicitly enforces the publication gate."""

    def __init__(self, report):
        messages = '; '.join(issue.message for issue in report.issues[:4])
        super().__init__(f'Nature readiness {report.score}/100 未通过：{messages}')
        self.report = report


def _rgba_tuple(color):
    from matplotlib.colors import to_rgba

    try:
        rgba = to_rgba(color)
        return tuple(round(float(value), 3) for value in rgba)
    except (TypeError, ValueError):
        return None


def _figure_colors(fig):
    colors = set()
    for ax in fig.axes:
        for line in ax.lines:
            rgba = _rgba_tuple(line.get_color())
            if rgba:
                colors.add(rgba)
        for patch in ax.patches:
            for value in (patch.get_facecolor(), patch.get_edgecolor()):
                rgba = _rgba_tuple(value)
                if rgba and rgba[3] > 0.05:
                    colors.add(rgba)
        for collection in ax.collections:
            for values in (collection.get_facecolors(), collection.get_edgecolors()):
                for value in np.asarray(values).reshape(-1, 4)[:100]:
                    rgba = tuple(round(float(item), 3) for item in value)
                    if rgba[3] > 0.05:
                        colors.add(rgba)
    return colors


def _is_chromatic(rgba):
    from matplotlib.colors import rgb_to_hsv

    if rgba is None or rgba[3] <= 0.05:
        return False
    hsv = rgb_to_hsv(np.asarray(rgba[:3], dtype=float))
    return bool(hsv[1] >= 0.16 and 0.12 <= hsv[2] <= 0.98)


def _red_green_present(colors):
    from matplotlib.colors import rgb_to_hsv

    red = green = False
    for rgba in colors:
        if not _is_chromatic(rgba):
            continue
        hue, saturation, _ = rgb_to_hsv(np.asarray(rgba[:3], dtype=float))
        red |= bool(saturation > 0.25 and (hue < 0.06 or hue > 0.94))
        green |= bool(saturation > 0.25 and 0.22 < hue < 0.48)
    return red and green


def _bbox_overlap(first, second, tolerance=6.0):
    width = min(first.x1, second.x1) - max(first.x0, second.x0)
    height = min(first.y1, second.y1) - max(first.y0, second.y0)
    return width > 1 and height > 1 and width * height > tolerance


class FigureValidator:
    """对 Matplotlib Figure 生成可执行的 Nature readiness 报告。"""

    def __init__(self, ready_threshold=90):
        self.ready_threshold = int(ready_threshold)

    def validate(self, fig, spec: FigureSpec, export_paths: Optional[Mapping[str, str]] = None):
        from matplotlib.colors import to_rgba
        from matplotlib.legend import Legend
        from matplotlib.text import Text
        from PIL import Image

        style = get_style(spec.style)
        profile = style.profile(spec)
        issues = []

        def add(code, severity, message, deduction, **details):
            issues.append(ValidationIssue(
                code=code, severity=severity, message=message,
                deduction=int(deduction), details=details,
            ))

        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        width_in, height_in = fig.get_size_inches()
        actual_mm = (float(width_in * 25.4), float(height_in * 25.4))
        expected_mm = (profile.width_mm, profile.height_mm)
        size_error = max(abs(actual_mm[0] - expected_mm[0]), abs(actual_mm[1] - expected_mm[1]))
        if size_error > 0.6:
            add('canvas_size', 'error', '画布不符合 FigureSpec 的最终物理尺寸。', 14,
                expected_mm=expected_mm, actual_mm=actual_mm)

        texts = [
            text for text in fig.findobj(match=lambda item: isinstance(item, Text))
            if text.get_visible() and str(text.get_text()).strip()
        ]
        font_sizes = [float(text.get_fontsize()) for text in texts]
        min_allowed = 5.0 if spec.mode == 'nature_portfolio' else 5.5
        minimum_font = min(font_sizes) if font_sizes else np.nan
        if font_sizes and minimum_font < min_allowed - 0.05:
            add('minimum_font_size', 'error',
                f'最小字号 {minimum_font:.1f} pt 低于 {min_allowed:.1f} pt。', 14,
                minimum_font_pt=minimum_font, required_pt=min_allowed)

        font_families = {
            tuple(str(value) for value in text.get_fontfamily())
            for text in texts
        }
        if len(font_families) > 2:
            add('font_consistency', 'warning', '图中存在过多字体 family。', 6,
                families=[list(value) for value in sorted(font_families)])

        colors = _figure_colors(fig)
        chromatic_colors = {value for value in colors if _is_chromatic(value)}
        encodings = dict(getattr(fig, '_nature_encodings', {}) or {})
        color_semantics = str(encodings.get('color', '')).lower()
        continuous_color = any(token in color_semantics for token in (
            'continuous', 'fdr', 'z-score', 'zscore', 'expression', 'correlation', 'nes',
            'jaccard', 'log2fc', 'direction', 'loading', 'variance',
            'distance', 'ratio', 'qvalue', '-log10',
        ))
        if len(chromatic_colors) > 12 and not continuous_color:
            add('excessive_colors', 'warning',
                f'检测到 {len(chromatic_colors)} 种非中性色，超过 Nature Portfolio 建议。', 8,
                count=len(chromatic_colors))
        if (_red_green_present(colors) and not continuous_color
                and not any(encodings.get(name) for name in ('marker', 'shape', 'hatch'))):
            add('red_green_conflict', 'warning',
                '红色和绿色可能是唯一的类别编码；应增加 marker/shape/hatch。', 7)
        color_vision = audit_palette([value[:3] for value in chromatic_colors])
        if not continuous_color and not color_vision['pass']:
            add('color_vision_collision', 'warning',
                '分类色在色觉模拟下难以区分，应同时使用 marker/shape 或调整色板。', 8,
                **color_vision)

        figure_bounds = fig.bbox
        text_boxes = []
        out_of_bounds = []
        clipped_tick_texts = set()
        for ax in fig.axes:
            x_low, x_high = sorted(ax.get_xlim())
            y_low, y_high = sorted(ax.get_ylim())
            for minor in (False, True):
                for label in ax.get_xticklabels(minor=minor):
                    location = float(label.get_position()[0])
                    if location < x_low - 1e-9 or location > x_high + 1e-9:
                        clipped_tick_texts.add(id(label))
                for label in ax.get_yticklabels(minor=minor):
                    location = float(label.get_position()[1])
                    if location < y_low - 1e-9 or location > y_high + 1e-9:
                        clipped_tick_texts.add(id(label))
        for text in texts:
            if id(text) in clipped_tick_texts:
                continue
            try:
                box = text.get_window_extent(renderer=renderer)
            except Exception:
                continue
            if box.width <= 0 or box.height <= 0:
                continue
            # Matplotlib creates tick Text objects just beyond the selected
            # limits.  They are clipped by their axes and never appear in the
            # export, so they must not become false clipping/overlap failures.
            owner = getattr(text, 'axes', None)
            if text.get_clip_on() and owner is not None:
                axes_box = owner.get_window_extent(renderer=renderer)
                if (box.x1 < axes_box.x0 or box.x0 > axes_box.x1
                        or box.y1 < axes_box.y0 or box.y0 > axes_box.y1):
                    continue
            label = str(text.get_text())[:100]
            text_boxes.append((label, box, text))
            if (box.x0 < figure_bounds.x0 - 1 or box.x1 > figure_bounds.x1 + 1
                    or box.y0 < figure_bounds.y0 - 1 or box.y1 > figure_bounds.y1 + 1):
                out_of_bounds.append(label)
        overlaps = []
        for index, (left_label, left_box, left_text) in enumerate(text_boxes):
            for right_label, right_box, right_text in text_boxes[index + 1:]:
                # 同一个字符串的轴标题/图例副本仍需检查；仅跳过同一 Text 对象。
                if left_text is right_text:
                    continue
                if _bbox_overlap(left_box, right_box):
                    overlaps.append((left_label, right_label))
                    if len(overlaps) >= 20:
                        break
            if len(overlaps) >= 20:
                break
        if overlaps:
            # In the strict Nature Portfolio mode a single real text overlap
            # is enough to fail the publication gate.  Warnings remain useful
            # for Standard/Publication exploratory output.
            overlap_is_error = spec.mode == 'nature_portfolio'
            add('label_overlap', 'error' if overlap_is_error else 'warning',
                f'检测到 {len(overlaps)} 组文字重叠（最多报告 20 组）。',
                min(20, 5 + len(overlaps)), pairs=overlaps)
        if out_of_bounds:
            add('axis_clipping', 'error',
                f'检测到 {len(out_of_bounds)} 个文字元素越出画布。', 12,
                labels=out_of_bounds[:20])

        legend_overlaps = []
        legends = list(fig.findobj(match=lambda item: isinstance(item, Legend)))
        for legend in legends:
            try:
                legend_box = legend.get_window_extent(renderer=renderer)
            except Exception:
                continue
            for ax in fig.axes:
                # Colourbars and dedicated size-key axes are intentionally
                # occupied by their own legend; they are not data regions.
                # Treating those auxiliary axes as overlap failures made a
                # valid dotplot fail QA even though the Count key was clearly
                # separated from the plotting panel.
                if (legend.axes is ax and not getattr(ax, '_nature_auxiliary', False)
                        and _bbox_overlap(legend_box, ax.get_window_extent(renderer), 30.0)):
                    legend_overlaps.append(ax.get_title() or ax.get_xlabel() or 'axis')
        if legend_overlaps:
            add('legend_overlap', 'warning',
                '图例位于数据区域内，可能遮挡数据。', 5,
                axes=legend_overlaps)

        background = tuple(round(value, 3) for value in fig.get_facecolor())
        white = tuple(round(value, 3) for value in to_rgba('white'))
        if background != white:
            add('white_background', 'error', 'Nature Portfolio publication 图必须使用白色背景。', 10,
                background=background)

        for warning in getattr(fig, '_nature_semantic_warnings', []) or []:
            add('semantic_warning', 'warning', str(warning), 4)
        for error in getattr(fig, '_nature_semantic_errors', []) or []:
            add('semantic_error', 'error', str(error), 20)

        export_paths = dict(export_paths or {})
        for fmt in spec.formats:
            path = Path(export_paths.get(fmt, '')) if export_paths.get(fmt) else None
            if path is None or not path.is_file() or path.stat().st_size == 0:
                add(f'{fmt}_export', 'error', f'缺少有效的 {fmt.upper()} 导出。', 10)
        svg_path = Path(export_paths['svg']) if export_paths.get('svg') else None
        if svg_path and svg_path.is_file():
            svg_text = svg_path.read_text(encoding='utf-8', errors='ignore')
            if '<text' not in svg_text:
                add('svg_editable_text', 'error', 'SVG 中没有可编辑 <text> 元素。', 14)
        pdf_path = Path(export_paths['pdf']) if export_paths.get('pdf') else None
        if pdf_path and pdf_path.is_file():
            if not pdf_path.read_bytes()[:5] == b'%PDF-':
                add('pdf_export', 'error', 'PDF 文件头无效。', 12)
        png_path = Path(export_paths['png']) if export_paths.get('png') else None
        png_dpi = None
        if png_path and png_path.is_file():
            with Image.open(png_path) as image:
                dpi = image.info.get('dpi')
                if dpi:
                    png_dpi = float(dpi[0])
                    if abs(png_dpi - profile.dpi) / max(profile.dpi, 1) > 0.06:
                        add('png_dpi', 'warning',
                            f'PNG DPI {png_dpi:.0f} 与合同 {profile.dpi} 不一致。', 5,
                            actual_dpi=png_dpi, expected_dpi=profile.dpi)
        tiff_path = Path(export_paths['tiff']) if export_paths.get('tiff') else None
        tiff_dpi = None
        if tiff_path and tiff_path.is_file():
            with Image.open(tiff_path) as image:
                dpi = image.info.get('dpi')
                if dpi:
                    tiff_dpi = float(dpi[0])
                    if abs(tiff_dpi - profile.dpi) / max(profile.dpi, 1) > 0.06:
                        add('tiff_dpi', 'error',
                            f'TIFF DPI {tiff_dpi:.0f} 与合同 {profile.dpi} 不一致。', 8,
                            actual_dpi=tiff_dpi, expected_dpi=profile.dpi)

        score = max(0, 100 - sum(issue.deduction for issue in issues))
        has_error = any(issue.severity == 'error' for issue in issues)
        ready = score >= self.ready_threshold and not has_error
        status = 'ready' if ready and score >= 90 else ('warning' if ready else 'fail')
        metrics = {
            'canvas_mm': [round(value, 3) for value in actual_mm],
            'expected_canvas_mm': [round(value, 3) for value in expected_mm],
            'minimum_font_pt': None if not font_sizes else round(float(minimum_font), 2),
            'font_families': [list(value) for value in sorted(font_families)],
            'chromatic_color_count': len(chromatic_colors),
            'text_element_count': len(texts),
            'label_overlap_count': len(overlaps),
            'out_of_bounds_count': len(out_of_bounds),
            'legend_overlap_count': len(legend_overlaps),
            'png_dpi': png_dpi,
            'tiff_dpi': tiff_dpi,
            'color_vision': color_vision,
            'panel_grid': getattr(fig, '_nature_panel_grid', None),
        }
        return ValidationReport(
            score=score, status=status, ready=ready,
            issues=issues, metrics=metrics, spec=spec.to_dict(),
        )

    def require(self, report):
        """Enforce the automatic submission gate for an explicit release step."""

        if not report.ready:
            raise FigureReadinessError(report)
        return report
