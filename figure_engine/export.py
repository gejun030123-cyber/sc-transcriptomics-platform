"""保持最终物理尺寸的 SVG/PDF/PNG 导出。"""

from pathlib import Path

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from figure_engine.validator import FigureValidator


def export_figure(fig, output_stem, spec: FigureSpec, *, formats=None,
                  validator=None, report_path=None):
    """导出 FigureSpec 指定的格式，不使用 tight bbox 改变物理画布。"""
    style = get_style(spec.style)
    profile = style.profile(spec)
    export_formats = tuple(formats or spec.formats)
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(*profile.figsize, forward=True)
    paths = {}
    with style.context(spec):
        for fmt in export_formats:
            fmt = str(fmt).lower()
            if fmt not in {'svg', 'pdf', 'png', 'tiff'}:
                continue
            path = output_stem.with_suffix(f'.{fmt}')
            kwargs = {
                'format': fmt,
                'facecolor': style.background,
                'edgecolor': 'none',
                'bbox_inches': None,
            }
            if fmt in {'png', 'tiff'}:
                kwargs['dpi'] = profile.dpi
            if fmt == 'tiff':
                kwargs['pil_kwargs'] = {'compression': 'tiff_lzw'}
            fig.savefig(path, **kwargs)
            paths[fmt] = str(path)
    validation_spec = spec.with_updates(formats=tuple(export_formats))
    report = (validator or FigureValidator()).validate(
        fig, validation_spec, export_paths=paths,
    )
    if report_path:
        report.write_json(report_path)
    return paths, report


def export_registered_figure(fig, output_stem, spec: FigureSpec, *, category, label,
                             formats=None, qa_path=None, enforce_gate=False):
    """Register requested PNG/SVG/PDF/TIFF artifacts and readiness JSON."""
    formats = tuple(formats or spec.formats)
    validator = FigureValidator()
    paths, report = export_figure(
        fig, output_stem, spec.with_updates(formats=tuple(formats)),
        formats=formats, report_path=qa_path, validator=validator,
    )
    result_files = [
        {
            'file_path': path,
            'file_type': fmt,
            'category': category,
            'label': label,
        }
        for fmt, path in paths.items()
    ]
    if qa_path:
        result_files.append({
            'file_path': str(qa_path),
            'file_type': 'json',
            'category': 'info',
            'label': f'{label} · Figure Quality / Nature readiness {report.score}/100',
        })
    if enforce_gate:
        validator.require(report)
    return result_files, report
