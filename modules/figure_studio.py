"""Versioned figure-studio rendering for platform plots and uploaded images."""

import base64
import json
import io
import os
import re
import tempfile
import xml.etree.ElementTree as ET
from html import escape as xml_escape
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd

from config import Config
from models import AnalysisTask, FigureAsset, FigureVersion, ResultFile


IMAGE_TYPES = {'png', 'jpg', 'jpeg', 'svg'}
RASTER_EXTENSIONS = {'.png', '.jpg', '.jpeg'}
MAX_UPLOAD_BYTES = 15 * 1024 * 1024
MAX_HEATMAP_PICKER_SAMPLES = 300
MAX_HEATMAP_PICKER_GROUPS = 60
COLOR_MAPS = {'RdBu_r', 'coolwarm', 'viridis', 'magma', 'Blues', 'Greens', 'YlOrRd'}
CORRELATION_COLOR_MAPS = {'Blues', 'RdBu_r', 'viridis', 'YlOrRd'}
LEGEND_POSITIONS = {
    'upper right', 'upper left', 'lower right', 'lower left', 'center right',
    'center left', 'lower center', 'upper center', 'best', 'none',
}
_COLOR_RE = re.compile(r'^#[0-9a-fA-F]{6}$')
_SAFE_SVG_RE = re.compile(r'<\s*(script|foreignObject)\b|\bon[a-z]+\s*=|(?:href|xlink:href)\s*=\s*["\']\s*(?:https?:|file:|javascript:)', re.I)


class FigureStudioError(ValueError):
    """Raised for a safe, user-facing figure-studio validation failure."""


def _clean_text(value, limit=180):
    return str(value or '').strip()[:limit]


def _number(value, default=None, lower=None, upper=None):
    if value in (None, ''):
        return default
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if not np.isfinite(number):
        return default
    if lower is not None:
        number = max(lower, number)
    if upper is not None:
        number = min(upper, number)
    return number


def _color(value, default):
    candidate = str(value or '').strip()
    return candidate if _COLOR_RE.match(candidate) else default


def normalize_style(payload):
    """Accept a small, typed, render-safe styling payload from the browser."""
    payload = payload or {}
    styles = {
        'title': _clean_text(payload.get('title')),
        'x_label': _clean_text(payload.get('x_label')),
        'y_label': _clean_text(payload.get('y_label')),
        'font_family': _clean_text(payload.get('font_family') or 'Arial', 80),
        'font_size': _number(payload.get('font_size'), 11, 7, 28),
        'background': _color(payload.get('background'), '#FFFFFF'),
        'grid': bool(payload.get('grid', False)),
        'legend_position': str(payload.get('legend_position') or 'upper right'),
        'legend_text': _clean_text(payload.get('legend_text'), 300),
        'x_min': _number(payload.get('x_min')),
        'x_max': _number(payload.get('x_max')),
        'y_min': _number(payload.get('y_min')),
        'y_max': _number(payload.get('y_max')),
        'pvalue_threshold': _number(payload.get('pvalue_threshold'), 0.05, 1e-300, 1.0),
        'fc_threshold': _number(payload.get('fc_threshold'), 1.5, 1e-6, 1e6),
        'up_color': _color(payload.get('up_color'), '#B64342'),
        'down_color': _color(payload.get('down_color'), '#0F4D92'),
        'ns_color': _color(payload.get('ns_color'), '#C7CDD6'),
        'label_genes': _clean_text(payload.get('label_genes'), 500),
        'heatmap_cmap': str(payload.get('heatmap_cmap') or 'RdBu_r'),
        'heatmap_top_n': int(_number(payload.get('heatmap_top_n'), 50, 5, 200)),
        'heatmap_sample_labels': str(payload.get('heatmap_sample_labels') or 'auto'),
        'heatmap_sample_scope': str(payload.get('heatmap_sample_scope') or 'all'),
        'heatmap_selected_samples': _clean_text(payload.get('heatmap_selected_samples'), 12000),
        'heatmap_group_column': _clean_text(payload.get('heatmap_group_column'), 80),
        'heatmap_selected_groups': _clean_text(payload.get('heatmap_selected_groups'), 4000),
        'row_cluster': bool(payload.get('row_cluster', True)),
        'col_cluster': bool(payload.get('col_cluster', True)),
        'annotation_column': _clean_text(payload.get('annotation_column'), 80),
        'corr_colorscale': str(payload.get('corr_colorscale') or 'Blues'),
        'enrichment_top_n': int(_number(payload.get('enrichment_top_n'), 20, 3, 50)),
        'enrichment_overview_top_n': int(_number(payload.get('enrichment_overview_top_n'), 8, 3, 20)),
        'enrichment_database_scope': _clean_text(payload.get('enrichment_database_scope'), 1000),
        'enrichment_pathway_selection': str(payload.get('enrichment_pathway_selection') or 'top'),
        'enrichment_target_pathways': _clean_text(payload.get('enrichment_target_pathways'), 4000),
        'enrichment_gene_label_strategy': str(payload.get('enrichment_gene_label_strategy') or 'all'),
        'enrichment_target_genes': _clean_text(payload.get('enrichment_target_genes'), 4000),
        'enrichment_max_gene_labels': int(_number(payload.get('enrichment_max_gene_labels'), 30, 0, 80)),
        'enrichment_bp_color': _color(payload.get('enrichment_bp_color'), '#FDBE85'),
        'enrichment_cc_color': _color(payload.get('enrichment_cc_color'), '#7A6FA8'),
        'enrichment_mf_color': _color(payload.get('enrichment_mf_color'), '#B8A9D1'),
        'enrichment_kegg_color': _color(payload.get('enrichment_kegg_color'), '#7BC77B'),
        'enrichment_other_color': _color(payload.get('enrichment_other_color'), '#9CB8D8'),
    }
    if styles['legend_position'] not in LEGEND_POSITIONS:
        styles['legend_position'] = 'upper right'
    if styles['heatmap_cmap'] not in COLOR_MAPS:
        styles['heatmap_cmap'] = 'RdBu_r'
    if styles['heatmap_sample_labels'] not in {'auto', 'all', 'none'}:
        styles['heatmap_sample_labels'] = 'auto'
    if styles['heatmap_sample_scope'] not in {'all', 'deg_groups', 'selected_groups', 'selected_samples'}:
        styles['heatmap_sample_scope'] = 'all'
    if styles['corr_colorscale'] not in CORRELATION_COLOR_MAPS:
        styles['corr_colorscale'] = 'Blues'
    if styles['enrichment_pathway_selection'] not in {'top', 'selected', 'selected_plus_top'}:
        styles['enrichment_pathway_selection'] = 'top'
    if styles['enrichment_gene_label_strategy'] not in {'all', 'shared', 'selected', 'none'}:
        styles['enrichment_gene_label_strategy'] = 'all'
    if styles['x_min'] is not None and styles['x_max'] is not None and styles['x_min'] >= styles['x_max']:
        styles['x_min'] = styles['x_max'] = None
    if styles['y_min'] is not None and styles['y_max'] is not None and styles['y_min'] >= styles['y_max']:
        styles['y_min'] = styles['y_max'] = None
    return styles


def _validate_project_path(path, project_id):
    valid, _ = Config._validate_path(path, project_id)
    return valid and os.path.isfile(path)


def _heatmap_picker_options(data_path):
    """Return bounded sample/group choices for heatmap-only figure controls."""
    try:
        import anndata as ad

        adata = ad.read_h5ad(data_path, backed='r')
        try:
            samples = [str(value) for value in adata.obs_names.tolist()]
            obs = adata.obs.copy()
        finally:
            file_handle = getattr(adata, 'file', None)
            if file_handle is not None:
                file_handle.close()
    except Exception:
        return {}
    group_columns = {}
    for column in obs.columns:
        values = obs[column]
        if pd.api.types.is_numeric_dtype(values):
            continue
        categories = list(dict.fromkeys(str(value) for value in values.dropna().tolist()))
        if 1 < len(categories) <= MAX_HEATMAP_PICKER_GROUPS:
            group_columns[str(column)] = categories
    return {
        'samples': samples[:MAX_HEATMAP_PICKER_SAMPLES],
        'sample_count': len(samples),
        'sample_picker_truncated': len(samples) > MAX_HEATMAP_PICKER_SAMPLES,
        'group_columns': group_columns,
    }


def _heatmap_context(task, data_path):
    """Expose only display-relevant task metadata for semantic redraws."""
    try:
        params = json.loads((task.params_json if task else '') or '{}')
    except (TypeError, ValueError, json.JSONDecodeError):
        params = {}
    return {
        'heatmap_options': _heatmap_picker_options(data_path),
        'heatmap_comparison_label': _clean_text(params.get('deg_comparison_label'), 200),
        'heatmap_groupby': _clean_text(params.get('groupby'), 80),
    }


def _source_mode_for_result(result_file):
    filename = os.path.basename(result_file.file_path).lower()
    task = AnalysisTask.get_by_id(result_file.task_id)
    if (
        (filename.startswith('bulk_deg_volcano') or filename.startswith('bulk_deg_ma'))
        and task and task.module_name == 'bulk_deg'
    ):
        plot_kind = 'ma' if filename.startswith('bulk_deg_ma') else 'volcano'
        stem = f'bulk_deg_{plot_kind}'
        suffix = filename[len(stem):].rsplit('.', 1)[0]
        companion = os.path.join(Config.results_dir(result_file.project_id), f'bulk_deg_results{suffix}.csv')
        if _validate_project_path(companion, result_file.project_id):
            return f'bulk_{plot_kind}', companion, {}
    if filename == 'bulk_heatmap.png' and task and task.module_name == 'bulk_heatmap':
        if _validate_project_path(task.output_adata_path, result_file.project_id):
            return 'bulk_heatmap', task.output_adata_path, {}
    if filename == 'bulk_corr_heatmap.png' and task and task.module_name == 'bulk_heatmap':
        companion = os.path.join(Config.results_dir(result_file.project_id), 'bulk_correlation_pairs.csv')
        if _validate_project_path(companion, result_file.project_id):
            return 'bulk_correlation', companion, {}
    if filename.startswith('enrichment_') and filename.endswith('.png') and task and task.module_name == 'bulk_enrichment':
        output_key = filename.rsplit('.', 1)[0]
        if output_key.startswith('enrichment_overview_'):
            integrated = os.path.join(
                Config.results_dir(result_file.project_id),
                'enrichment_integrated_results.csv',
            )
            if _validate_project_path(integrated, result_file.project_id):
                return 'bulk_enrichment_overview', integrated, {}
        companion = os.path.join(Config.results_dir(result_file.project_id), f'{output_key}_results.csv')
        if not _validate_project_path(companion, result_file.project_id):
            # Extended Nature views use a suffixed image stem (…_chord,
            # …_cnetplot, …_emapplot), while all views intentionally share
            # the unsuffixed enrichment result table.  Resolve that table so
            # the result-page “美化” button opens the semantic redraw controls
            # instead of silently falling back to style-only editing.
            for view_suffix in (
                'dotplot', 'barplot', 'chord', 'cnetplot', 'emapplot',
                'gsea_running', 'gsea',
            ):
                marker = f'_{view_suffix}'
                if output_key.endswith(marker):
                    base_key = output_key[:-len(marker)]
                    candidate = os.path.join(
                        Config.results_dir(result_file.project_id),
                        f'{base_key}_results.csv',
                    )
                    if _validate_project_path(candidate, result_file.project_id):
                        companion = candidate
                    break
        if _validate_project_path(companion, result_file.project_id):
            return 'bulk_enrichment', companion, {}
    # Single-cell DEG figures (cell-level or pseudobulk) are data redraws
    # with manual gene labelling.
    if task and task.module_name in {'sc_cell_deg', 'sc_pseudobulk_deg'}:
        for plot_kind in ('volcano', 'ma'):
            sc_deg = _sc_deg_figure_context(result_file, task, plot_kind)
            if sc_deg:
                data_path, context = sc_deg
                return f'sc_{plot_kind}', data_path, context
    # Single-cell GO enrichment figures use the same pathway-selection controls
    # as Bulk RNA-seq; keep the running curve non-redrawable because its score
    # depends on the full ranking object that a CSV cannot rebuild.
    if task and task.module_name == 'sc_cell_go':
        sc_enrichment = _sc_enrichment_context(result_file, task)
        if sc_enrichment:
            data_path, context = sc_enrichment
            return 'sc_enrichment', data_path, context
    return 'style_only', '', {}


def _sc_safe_name(value, fallback='value'):
    """Recreate the SC module name-safe fragment used in figure stems."""
    text = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value or '').strip()).strip('._-')
    return text or fallback


def _task_params(task_id):
    """Return the analysis-task parameters as a plain dict."""
    task = AnalysisTask.get_by_id(task_id)
    if not task:
        return {}
    try:
        params = json.loads(task.params_json or '{}')
    except (TypeError, ValueError, json.JSONDecodeError):
        params = {}
    return params if isinstance(params, dict) else {}


def _sc_package_csv_paths(project_id, folder_key):
    """Return project-local SC CSV package paths for one numeric subfolder.

    folder_key is 'deg' or 'go'.  Both registered user-visible tables and the
    protected '.internal' copies are considered so a figure can be redrawn even
    when 'export_full_tables' was left off.
    """
    folder_names = {
        'deg': '03_differential_expression',
        'go': '04_go_enrichment',
    }
    dirname = folder_names.get(folder_key)
    if not dirname:
        return []
    results_root = os.path.join(Config.project_dir(project_id), 'results')
    if not os.path.isdir(results_root):
        return []
    paths = []
    seen = set()
    for package in sorted(os.listdir(results_root)):
        package_dir = os.path.join(results_root, package)
        target_dir = os.path.join(package_dir, dirname)
        if not os.path.isdir(target_dir):
            continue
        for dirpath, _dirnames, filenames in os.walk(target_dir):
            for name in filenames:
                if not name.lower().endswith('.csv'):
                    continue
                candidate = os.path.join(dirpath, name)
                if candidate in seen or not _validate_project_path(candidate, project_id):
                    continue
                seen.add(candidate)
                paths.append(candidate)
    return paths


def _sc_task_csv_paths(result_file, folder_key):
    """Find task-artifact tables first, then legacy SC result packages."""
    candidates = []
    seen = set()
    artifact_dir = os.path.dirname(result_file.file_path)
    if os.path.isdir(artifact_dir):
        for name in sorted(os.listdir(artifact_dir)):
            if not name.lower().endswith('.csv'):
                continue
            candidate = os.path.join(artifact_dir, name)
            if candidate in seen or not _validate_project_path(candidate, result_file.project_id):
                continue
            seen.add(candidate)
            candidates.append(candidate)
    for candidate in _sc_package_csv_paths(result_file.project_id, folder_key):
        if candidate not in seen:
            seen.add(candidate)
            candidates.append(candidate)
    return candidates


def _read_csv_identity(path):
    """Read only the small provenance columns needed to bind an SC figure."""
    try:
        frame = pd.read_csv(path)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError, ValueError):
        return None

    def uniq(column):
        if column not in frame.columns:
            return []
        seen = []
        for value in frame[column].dropna().tolist():
            text = str(value).strip()
            if text and text not in seen:
                seen.append(text)
        return seen

    return {
        'path': path,
        'columns': set(frame.columns),
        'comparison_ids': uniq('comparison_id'),
        'deg_scopes': uniq('deg_scope'),
        'clusters': uniq('cluster'),
        'directions': uniq('direction'),
        'gene_sets': uniq('gene_set'),
        'methods': uniq('method'),
    }


def _sc_deg_figure_context(result_file, task, plot_kind):
    """Bind one SC DEG Volcano/MA figure to its exact comparison table."""
    filename = os.path.basename(result_file.file_path)
    cell_marker = f'sc_cell_deg_{plot_kind}_'
    pseudobulk_marker = f'sc_pseudobulk_{plot_kind}_'
    if cell_marker in filename:
        token = cell_marker
        default_prefix = 'sc_cell_level'
        module_name = 'sc_cell_deg'
    elif pseudobulk_marker in filename:
        token = pseudobulk_marker
        default_prefix = 'sc_pseudobulk'
        module_name = 'sc_pseudobulk_deg'
    else:
        return None
    token_index = filename.find(token)
    stem = filename.rsplit('.', 1)[0][token_index + len(token):]
    params = _task_params(task.id)
    prefix = _sc_safe_name(params.get('export_prefix', default_prefix), default_prefix)
    for candidate in _sc_task_csv_paths(result_file, 'deg'):
        identity = _read_csv_identity(candidate)
        if not identity:
            continue
        if 'gene' not in identity['columns']:
            continue
        if not any(column in identity['columns'] for column in ('log2FC', 'avg_log2FC', 'logFC')):
            continue
        if not any(column in identity['columns'] for column in ('padj', 'p.adjust', 'p_val_adj')):
            continue
        if module_name == 'sc_cell_deg':
            for comparison_id in identity['comparison_ids']:
                for deg_scope in identity['deg_scopes']:
                    for cluster in identity['clusters']:
                        expected = (
                            f"{prefix}_{_sc_safe_name(comparison_id)}_"
                            f"{_sc_safe_name(deg_scope)}_{_sc_safe_name(cluster)}"
                        )
                        if expected == stem:
                            return candidate, {
                                'module': module_name,
                                'comparison_id': comparison_id,
                                'deg_scope': deg_scope,
                                'cluster': cluster,
                                'evidence_role': 'discovery',
                            }
        else:
            for comparison_id in identity['comparison_ids']:
                for cluster in identity['clusters']:
                    expected = (
                        f"{prefix}_{_sc_safe_name(comparison_id)}_"
                        f"{_sc_safe_name(cluster)}"
                    )
                    if expected == stem:
                        return candidate, {
                            'module': module_name,
                            'comparison_id': comparison_id,
                            'deg_scope': '',
                            'cluster': cluster,
                            'evidence_role': 'comparison',
                        }
    return None


_ENRICHMENT_SUFFIX_BY_PLOT_TYPE = {
    'enrichment_dotplot': 'dotplot',
    'enrichment_barplot': 'barplot',
    'enrichment_chord': 'chord',
    'enrichment_cnetplot': 'cnetplot',
    'enrichment_emapplot': 'emapplot',
    'gsea': 'gsea',
}


def _sc_enrichment_context(result_file, task):
    """Bind an SC GO enrichment figure to its pathway table and cluster/direction."""
    filename = os.path.basename(result_file.file_path)
    marker = 'sc_cell_go_'
    marker_index = filename.find(marker)
    if marker_index < 0:
        return None
    plot_type = _enrichment_plot_type_from_filename(result_file.file_path, '')
    if plot_type not in _ENRICHMENT_SUFFIX_BY_PLOT_TYPE:
        return None
    suffix = _ENRICHMENT_SUFFIX_BY_PLOT_TYPE[plot_type]
    stem = filename.rsplit('.', 1)[0]
    # Task-artifact registration prefixes files with a stable sequence number
    # (for example ``004_sc_cell_go_…``); historical package files start
    # directly with the marker.  Both name forms identify the same analysis
    # contract.
    base = stem[marker_index + len(marker):]
    for candidate in _sc_task_csv_paths(result_file, 'go'):
        identity = _read_csv_identity(candidate)
        if not identity:
            continue
        if not {'Term', 'Adjusted P-value', 'Overlap'}.issubset(identity['columns']):
            continue
        for gene_set in identity['gene_sets']:
            for deg_scope in identity['deg_scopes']:
                for comparison_id in identity['comparison_ids']:
                    prefix = (
                        f"{_sc_safe_name(gene_set)}_"
                        f"{_sc_safe_name(deg_scope)}_{_sc_safe_name(comparison_id)}_"
                    )
                    if not base.startswith(prefix):
                        continue
                    remainder = base[len(prefix):]
                    if remainder.endswith('_' + suffix):
                        remainder = remainder[:-(len(suffix) + 1)]
                    for direction in ('Up', 'Down', 'All'):
                        if remainder.endswith('_' + direction):
                            cluster = remainder[:-(len(direction) + 1)]
                            if cluster in identity['clusters'] and direction in identity['directions']:
                                return candidate, {
                                    'module': 'sc_cell_go',
                                    'cluster': cluster,
                                    'direction': direction,
                                    'method': identity['methods'][0] if identity['methods'] else '',
                                    'plot_type': plot_type,
                                }
    return None


def _enrichment_plot_type_from_filename(file_path, data_path=''):
    """Recover the semantic Nature renderer represented by an enrichment PNG.

    Older analysis runs used an unsuffixed ``enrichment_*.png`` name for a
    legacy bar renderer, while current runs use that same stem for the primary
    ORA dotplot and put alternate views behind explicit suffixes.  For an
    unsuffixed image, inspect the companion table's semantic columns so the
    Figure Studio button opens the same dotplot the user is looking at.
    """
    stem = Path(str(file_path)).stem.lower()
    if stem.startswith('enrichment_overview_'):
        return 'enrichment_overview'
    suffixes = {
        '_dotplot': 'enrichment_dotplot',
        '_barplot': 'enrichment_barplot',
        '_chord': 'enrichment_chord',
        '_cnetplot': 'enrichment_cnetplot',
        '_emapplot': 'enrichment_emapplot',
        '_gsea_running': 'gsea_running',
        '_gsea': 'gsea',
    }
    for suffix, plot_type in suffixes.items():
        if stem.endswith(suffix):
            return plot_type
    if data_path and os.path.isfile(data_path):
        try:
            columns = {str(column).strip().lower() for column in pd.read_csv(data_path, nrows=0).columns}
        except (OSError, pd.errors.EmptyDataError, UnicodeDecodeError):
            columns = set()
        # GeneRatio + an explicit overlap/count column are the lossless ORA
        # dotplot contract.  Older tables without GeneRatio remain on the
        # legacy compatibility renderer until regenerated by the task.
        ratio_columns = {'generatio', 'gene ratio'}
        count_columns = {'count', 'num', 'gene count', 'overlap', 'setsize'}
        if columns.intersection(ratio_columns) and columns.intersection(count_columns):
            return 'enrichment_dotplot'
    # A legacy unsuffixed enrichment PNG may have been rendered by the older
    # ontology-colour bar renderer; keep that compatibility path distinct from
    # an explicit Nature ``_dotplot`` output.
    return 'legacy_enrichment'


def resolve_source(project_id, source_kind, source_id):
    """Resolve a project-scoped source and its available editing capabilities."""
    if source_kind == 'result_file':
        source = ResultFile.get_by_id(source_id)
        if not source or source.project_id != project_id or source.file_type not in IMAGE_TYPES:
            raise FigureStudioError('未找到可编辑的平台图。')
        if not _validate_project_path(source.file_path, project_id):
            raise FigureStudioError('原图路径无效。')
        edit_mode, data_path, sc_context = _source_mode_for_result(source)
        payload = {
            'kind': source_kind, 'id': source.id, 'label': source.label,
            'file_type': source.file_type, 'file_path': source.file_path,
            'edit_mode': edit_mode, 'data_path': data_path,
        }
        if edit_mode in {'bulk_enrichment', 'bulk_enrichment_overview', 'sc_enrichment'}:
            payload['plot_type'] = (
                sc_context.get('plot_type', '')
                if edit_mode == 'sc_enrichment' and sc_context
                else _enrichment_plot_type_from_filename(source.file_path, data_path)
            )
        if edit_mode in {'sc_volcano', 'sc_ma', 'sc_enrichment'} and sc_context:
            payload['sc_context'] = sc_context
        if edit_mode == 'bulk_heatmap':
            payload.update(_heatmap_context(AnalysisTask.get_by_id(source.task_id), data_path))
        return payload
    if source_kind == 'asset':
        source = FigureAsset.get_by_id(source_id)
        if not source or source.project_id != project_id or source.file_type not in IMAGE_TYPES:
            raise FigureStudioError('未找到上传图片。')
        if not _validate_project_path(source.file_path, project_id):
            raise FigureStudioError('上传图片路径无效。')
        return {
            'kind': source_kind, 'id': source.id, 'label': source.label,
            'file_type': source.file_type, 'file_path': source.file_path,
            'edit_mode': 'style_only', 'data_path': '',
        }
    if source_kind == 'version':
        version = FigureVersion.get_by_id(source_id)
        if not version or version.project_id != project_id:
            raise FigureStudioError('未找到已保存的图形版本。')
        preview_path = version.png_path or version.svg_path
        if not _validate_project_path(preview_path, project_id):
            raise FigureStudioError('已保存图形版本的路径无效。')
        try:
            saved_style = json.loads(version.style_json or '{}')
        except (TypeError, ValueError):
            saved_style = {}
        # A saved data-backed version is intentionally re-rendered from the
        # original analysis table/data.  This keeps later edits semantic rather
        # than applying raster edits on top of a prior preview.
        try:
            base = resolve_source(project_id, version.source_kind, version.source_id)
        except FigureStudioError:
            base = {
                'file_path': preview_path,
                'file_type': Path(preview_path).suffix.lstrip('.').lower(),
                'edit_mode': 'style_only', 'data_path': '',
            }
        base.update({
            'kind': 'version', 'id': version.id, 'label': version.label,
            'saved_style': saved_style,
        })
        return base
    raise FigureStudioError('不支持的图形来源。')


def list_sources(project_id):
    """List platform images and user uploads, preferring PNG over paired SVG."""
    sources = []
    seen = set()
    files = [item for item in ResultFile.get_by_project(project_id) if item.file_type in IMAGE_TYPES]
    files.sort(key=lambda item: (item.task_id, item.category, item.label, item.file_type != 'png'))
    for item in files:
        if not _validate_project_path(item.file_path, project_id):
            continue
        key = (item.task_id, item.category, item.label)
        if key in seen:
            continue
        seen.add(key)
        mode, data_path, sc_context = _source_mode_for_result(item)
        source_payload = {
            'kind': 'result_file', 'id': item.id, 'label': item.label or '平台分析图',
            'file_type': item.file_type, 'edit_mode': mode,
            'preview_url': f'/projects/{project_id}/results/file/{item.id}',
        }
        if mode in {'bulk_enrichment', 'bulk_enrichment_overview', 'sc_enrichment'}:
            source_payload['plot_type'] = (
                sc_context.get('plot_type', '')
                if mode == 'sc_enrichment' and sc_context
                else _enrichment_plot_type_from_filename(item.file_path, data_path)
            )
        if mode in {'sc_volcano', 'sc_ma', 'sc_enrichment'} and sc_context:
            source_payload['sc_context'] = sc_context
        if mode == 'bulk_heatmap':
            task = AnalysisTask.get_by_id(item.task_id)
            if task:
                source_payload.update(_heatmap_context(task, data_path))
        sources.append(source_payload)
    for item in FigureAsset.get_by_project(project_id):
        if _validate_project_path(item.file_path, project_id):
            sources.append({
                'kind': 'asset', 'id': item.id, 'label': item.label,
                'file_type': item.file_type, 'edit_mode': 'style_only',
                'preview_url': f'/projects/{project_id}/figure-studio/assets/{item.id}',
            })
    for item in FigureVersion.get_by_project(project_id):
        preview_path = item.png_path or item.svg_path
        if not _validate_project_path(preview_path, project_id):
            continue
        try:
            saved_style = json.loads(item.style_json or '{}')
        except (TypeError, ValueError):
            saved_style = {}
        source_payload = {
            'kind': 'version', 'id': item.id, 'label': f'已保存：{item.label}',
            'file_type': Path(preview_path).suffix.lstrip('.').lower(),
            'edit_mode': item.edit_mode or 'style_only',
            'preview_url': f'/projects/{project_id}/figure-studio/versions/{item.id}/'
                           f'{"png" if item.png_path else "svg"}',
            'style': saved_style,
        }
        if source_payload['edit_mode'] == 'bulk_heatmap':
            try:
                resolved = resolve_source(project_id, 'version', item.id)
                source_payload.update({
                    key: resolved[key] for key in ('heatmap_options', 'heatmap_comparison_label', 'heatmap_groupby')
                    if key in resolved
                })
            except FigureStudioError:
                pass
        sources.append(source_payload)
    return sources


def _apply_common_style(fig, ax, style):
    from modules.figure_style import NATURE_AXIS, NATURE_GRID, NATURE_TEXT, _font_for_text

    fig.patch.set_facecolor(style['background'])
    ax.set_facecolor(style['background'])
    if style['title']:
        ax.set_title(style['title'], loc='left', fontweight='semibold')
    if style['x_label']:
        ax.set_xlabel(style['x_label'])
    if style['y_label']:
        ax.set_ylabel(style['y_label'])
    if style['x_min'] is not None or style['x_max'] is not None:
        left, right = ax.get_xlim()
        ax.set_xlim(style['x_min'] if style['x_min'] is not None else left,
                    style['x_max'] if style['x_max'] is not None else right)
    if style['y_min'] is not None or style['y_max'] is not None:
        bottom, top = ax.get_ylim()
        ax.set_ylim(style['y_min'] if style['y_min'] is not None else bottom,
                    style['y_max'] if style['y_max'] is not None else top)
    if style['grid']:
        ax.grid(True, color=NATURE_GRID, linewidth=0.65, alpha=0.8)
    else:
        ax.grid(False)
    ax.set_axisbelow(True)
    for text in fig.findobj(lambda artist: hasattr(artist, 'set_fontfamily') and hasattr(artist, 'get_text')):
        try:
            text.set_fontfamily(_font_for_text(text.get_text(), style['font_family']))
            text.set_fontsize(style['font_size'])
        except Exception:
            continue
    ax.title.set_fontsize(style['font_size'] + 2)
    ax.xaxis.label.set_fontsize(style['font_size'])
    ax.yaxis.label.set_fontsize(style['font_size'])
    ax.xaxis.label.set_color(NATURE_TEXT)
    ax.yaxis.label.set_color(NATURE_TEXT)
    ax.tick_params(labelsize=max(7, style['font_size'] - 1), colors=NATURE_AXIS)
    for spine in ax.spines.values():
        spine.set_color('#98A2B3')
        spine.set_linewidth(0.75)
    legend = ax.get_legend()
    if legend is not None:
        if style['legend_position'] == 'none':
            legend.remove()
        else:
            try:
                legend.set_loc(style['legend_position'])
            except AttributeError:
                legend._loc = style['legend_position']
    return fig


def _annotate_selected_genes(ax, deg_df, style):
    labels = [item.strip() for item in re.split(r'[,;\n]+', style['label_genes']) if item.strip()]
    if not labels:
        return
    lookup = deg_df.set_index('gene')
    for gene in labels[:10]:
        if gene not in lookup.index:
            continue
        row = lookup.loc[gene]
        x = float(row['log2FC'])
        y = float(-np.log10(max(float(row['padj']), np.finfo(float).tiny)))
        ax.annotate(str(gene), xy=(x, y), xytext=(8, 8), textcoords='offset points',
                    fontsize=max(7, style['font_size'] - 2), color='#111827',
                    bbox={'boxstyle': 'round,pad=0.18', 'facecolor': 'white', 'edgecolor': '#98A2B3', 'linewidth': 0.55},
                    arrowprops={'arrowstyle': '-', 'color': '#667085', 'linewidth': 0.55})


def _render_bulk_volcano(data_path, style, label):
    import matplotlib.pyplot as plt
    from modules.bulk_deg import _draw_bulk_volcano

    deg_df = pd.read_csv(data_path)
    required = {'gene', 'log2FC', 'padj'}
    if not required.issubset(deg_df.columns):
        raise FigureStudioError('火山图对应的 DEG 结果缺少 gene、log2FC 或 padj 列。')
    deg_df = deg_df.copy()
    deg_df['log2FC'] = pd.to_numeric(deg_df['log2FC'], errors='coerce')
    deg_df['padj'] = pd.to_numeric(deg_df['padj'], errors='coerce').fillna(1.0)
    log_fc_cutoff = np.log2(style['fc_threshold'])
    significant = deg_df['padj'] < style['pvalue_threshold']
    deg_df['regulation'] = np.select(
        [significant & (deg_df['log2FC'] >= log_fc_cutoff),
         significant & (deg_df['log2FC'] <= -log_fc_cutoff)],
        ['Up', 'Down'], default='NS',
    )
    fig, ax = plt.subplots(figsize=(8.2, 5.8), dpi=150)
    _draw_bulk_volcano(
        ax, deg_df, pval_threshold=style['pvalue_threshold'],
        fc_threshold=style['fc_threshold'], title=style['title'] or label,
        top_n=0 if style['label_genes'] else 8, show_legend=True,
        colors={'Up': style['up_color'], 'Down': style['down_color'], 'NS': style['ns_color']},
    )
    _annotate_selected_genes(ax, deg_df, style)
    _apply_common_style(fig, ax, style)
    return fig, 'data_redraw'


def _load_sc_deg_frame(data_path, context):
    """Load one exact SC DEG comparison/cluster unit for semantic redraws."""
    deg_df = pd.read_csv(data_path)
    gene_col = 'gene' if 'gene' in deg_df.columns else 'names'
    if gene_col not in deg_df.columns:
        raise FigureStudioError('差异表达图对应的 SC DEG 结果缺少基因列。')
    fc_col = next(
        (column for column in ('log2FC', 'avg_log2FC', 'logFC') if column in deg_df.columns),
        None,
    )
    fdr_col = next(
        (column for column in ('padj', 'p.adjust', 'p_val_adj') if column in deg_df.columns),
        None,
    )
    if fc_col is None or fdr_col is None:
        raise FigureStudioError('差异表达图对应的 SC DEG 结果缺少 log2FC 或 padj 列。')
    deg_df = deg_df.copy()
    deg_df['gene'] = deg_df[gene_col].astype(str)
    deg_df['log2FC'] = pd.to_numeric(deg_df[fc_col], errors='coerce')
    deg_df['padj'] = pd.to_numeric(deg_df[fdr_col], errors='coerce').fillna(1.0)
    for key, value in context.items() or {}:
        if key in {'module', 'evidence_role', 'plot_type'}:
            continue
        if value and key in deg_df.columns:
            deg_df = deg_df[deg_df[key].astype(str) == str(value)]
    if deg_df.empty:
        raise FigureStudioError('未找到该差异表达图对应的比较单元数据。')
    return deg_df


def _render_sc_volcano(data_path, style, label, context):
    """Redraw an SC cell/pseudobulk volcano with manual gene labels."""
    import math

    from figure_engine import NatureFigureDirector

    deg_df = _load_sc_deg_frame(data_path, context)
    custom_genes = tuple(
        item.strip() for item in re.split(r'[,\n;]+', style['label_genes']) if item.strip()
    )
    director = NatureFigureDirector()
    logfc_cutoff = math.log2(style['fc_threshold'])
    spec = director.spec_from_params(
        'volcano', {},
        width='single',
        title=style['title'] or label,
        evidence_role=context.get('evidence_role', 'discovery'),
        fc_threshold=logfc_cutoff,
        fdr_threshold=style['pvalue_threshold'],
        label_n=0 if custom_genes else 8,
        label_strategy='none' if custom_genes else 'top_significant',
        label_genes=custom_genes,
        show_legend=True,
    )
    fig = director.render(spec, deg_df)
    _apply_common_style(fig, fig.axes[0], style)
    return fig, 'data_redraw'


def _ma_expression_column(frame):
    """Return the observed mean-expression field required by an MA plot."""
    for column in ('mean_expression', 'base_mean_count', 'baseMean'):
        if column in frame.columns:
            return pd.to_numeric(frame[column], errors='coerce')
    if {'mean_group1', 'mean_group2'}.issubset(frame.columns):
        return (
            pd.to_numeric(frame['mean_group1'], errors='coerce')
            + pd.to_numeric(frame['mean_group2'], errors='coerce')
        ) / 2.0
    return None


def _render_ma_frame(deg_df, style, label, *, evidence_role='comparison'):
    """Redraw a data-backed MA plot while retaining manual gene labels."""
    import math

    from figure_engine import NatureFigureDirector

    expression = _ma_expression_column(deg_df)
    if expression is None:
        raise FigureStudioError(
            'MA 图对应的 DEG 结果缺少 mean_expression、baseMean 或组均值列。'
        )
    frame = deg_df.copy()
    frame['mean_expression'] = expression
    custom_genes = tuple(
        item.strip() for item in re.split(r'[,\n;]+', style['label_genes']) if item.strip()
    )
    director = NatureFigureDirector()
    spec = director.spec_from_params(
        'ma', {}, width='single', title=style['title'] or label,
        evidence_role=evidence_role,
        fc_threshold=math.log2(style['fc_threshold']),
        fdr_threshold=style['pvalue_threshold'],
        label_n=0 if custom_genes else 6,
        label_strategy='none' if custom_genes else 'top_significant',
        label_genes=custom_genes, show_legend=True,
    )
    figure = director.render(spec, frame)
    _apply_common_style(figure, figure.axes[0], style)
    return figure, 'data_redraw'


def _render_bulk_ma(data_path, style, label):
    deg_df = pd.read_csv(data_path)
    required = {'gene', 'log2FC', 'padj'}
    if not required.issubset(deg_df.columns):
        raise FigureStudioError('MA 图对应的 DEG 结果缺少 gene、log2FC 或 padj 列。')
    deg_df = deg_df.copy()
    deg_df['log2FC'] = pd.to_numeric(deg_df['log2FC'], errors='coerce')
    deg_df['padj'] = pd.to_numeric(deg_df['padj'], errors='coerce').fillna(1.0)
    return _render_ma_frame(deg_df, style, label, evidence_role='comparison')


def _render_sc_ma(data_path, style, label, context):
    return _render_ma_frame(
        _load_sc_deg_frame(data_path, context), style, label,
        evidence_role=context.get('evidence_role', 'discovery'),
    )


def _dense_matrix(values):
    return values.toarray() if hasattr(values, 'toarray') else np.asarray(values)


def _render_bulk_heatmap(data_path, style, label, comparison_label='', default_groupby=''):
    import anndata as ad
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from modules.figure_style import NATURE_PALETTE, NATURE_TEXT
    from modules.visualization import cluster_heatmap
    from modules.native_figures import apply_sample_tick_labels

    adata = ad.read_h5ad(data_path)
    matrix = np.nan_to_num(_dense_matrix(adata.X).astype(float), nan=0.0, posinf=0.0, neginf=0.0)
    if matrix.ndim != 2 or matrix.shape[1] < 2:
        raise FigureStudioError('热图对应数据不足以重新绘制。')
    from modules.bulk_heatmap import _select_heatmap_display_samples
    try:
        sample_indices, _ = _select_heatmap_display_samples(
            adata.obs,
            sample_display_mode=style['heatmap_sample_scope'],
            groupby=style['heatmap_group_column'] or default_groupby,
            selected_groups=style['heatmap_selected_groups'],
            selected_samples=style['heatmap_selected_samples'],
            comparison_label=comparison_label,
        )
    except ValueError as exc:
        raise FigureStudioError(str(exc)) from exc
    matrix = matrix[sample_indices, :]
    display_obs = adata.obs.iloc[sample_indices].copy()
    top_n = min(style['heatmap_top_n'], matrix.shape[1])
    variances = np.nanvar(matrix, axis=0)
    gene_index = np.argsort(variances)[::-1][:top_n]
    data = matrix[:, gene_index]
    mean = data.mean(axis=0, keepdims=True)
    std = data.std(axis=0, keepdims=True)
    data = (data - mean) / np.where(std == 0, 1.0, std)
    if style['col_cluster'] and data.shape[0] > 2:
        row_order = cluster_heatmap(data, method='ward', metric='euclidean')
    else:
        row_order = list(range(data.shape[0]))
    if style['row_cluster'] and data.shape[1] > 2:
        column_order = cluster_heatmap(data.T, method='ward', metric='euclidean')
    else:
        column_order = list(range(data.shape[1]))
    data = data[np.ix_(row_order, column_order)]
    samples = [str(display_obs.index[index]) for index in row_order]
    genes = [str(adata.var_names[gene_index[index]]) for index in column_order]
    width = max(8.0, min(15.0, 5.8 + 0.16 * len(genes)))
    height = max(5.4, min(14.0, 3.8 + 0.18 * len(samples)))
    fig, ax = plt.subplots(figsize=(width, height), dpi=150)
    image = ax.imshow(data, aspect='auto', interpolation='nearest', cmap=style['heatmap_cmap'], vmin=-3, vmax=3)
    ax.set_xticks(np.arange(len(genes)), genes, rotation=45, ha='right')
    if style['heatmap_sample_labels'] == 'none':
        ax.set_yticks([])
    elif style['heatmap_sample_labels'] == 'all':
        ax.set_yticks(np.arange(len(samples)), samples)
    else:
        apply_sample_tick_labels(
            ax, samples, axis='y', max_labels=18, font_size=style['font_size'],
        )
    ax.set_title(style['title'] or label, loc='left', fontsize=12, fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(style['x_label'] or 'Gene')
    ax.set_ylabel(style['y_label'] or 'Sample')
    colorbar = fig.colorbar(image, ax=ax, fraction=0.038, pad=0.025, aspect=32)
    colorbar.outline.set_visible(False)
    colorbar.set_label('Z-score', fontsize=max(7, style['font_size'] - 1))
    annotation = style['annotation_column']
    if annotation and annotation in display_obs.columns:
        values = display_obs[annotation].astype(str).to_numpy()[row_order]
        categories = list(dict.fromkeys(values.tolist()))
        colors = {category: NATURE_PALETTE[i % len(NATURE_PALETTE)] for i, category in enumerate(categories)}
        for idx, value in enumerate(values):
            ax.add_patch(Rectangle(
                (-0.035, idx - 0.5), 0.018, 1.0, transform=ax.get_yaxis_transform(),
                clip_on=False, facecolor=colors[value], edgecolor='none', zorder=5,
            ))
        handles = [Rectangle((0, 0), 1, 1, facecolor=colors[item], edgecolor='none', label=item)
                   for item in categories]
        ax.legend(handles=handles, title=annotation, loc='upper left', bbox_to_anchor=(1.02, 1.0),
                  frameon=False, fontsize=max(7, style['font_size'] - 2))
    _apply_common_style(fig, ax, style)
    fig.tight_layout()
    return fig, 'data_redraw'


def _render_bulk_correlation(data_path, style, label):
    from modules.native_figures import correlation_heatmap_figure

    pairs = pd.read_csv(data_path)
    required = {'sample_1', 'sample_2'}
    if not required.issubset(pairs.columns):
        raise FigureStudioError('相关性结果缺少 sample_1 或 sample_2 列。')
    score_column = next(
        (column for column in pairs.columns
         if str(column).lower().endswith('_r') or str(column).lower() in {'correlation', 'r'}),
        None,
    )
    if score_column is None:
        raise FigureStudioError('相关性结果缺少 Pearson/Spearman r 数值列。')
    labels = list(dict.fromkeys(
        [str(value) for value in pairs['sample_1'].dropna().tolist()]
        + [str(value) for value in pairs['sample_2'].dropna().tolist()]
    ))
    if len(labels) < 2:
        raise FigureStudioError('相关性结果中的样本数不足。')
    index = {item: position for position, item in enumerate(labels)}
    matrix = np.eye(len(labels), dtype=float)
    for _, row in pairs.iterrows():
        first, second = str(row['sample_1']), str(row['sample_2'])
        value = _number(row.get(score_column))
        if value is None or first not in index or second not in index:
            continue
        matrix[index[first], index[second]] = value
        matrix[index[second], index[first]] = value
    group_map = {}
    for sample_column, group_column in (('sample_1', 'group_1'), ('sample_2', 'group_2')):
        if group_column not in pairs.columns:
            continue
        for sample, group in zip(pairs[sample_column], pairs[group_column]):
            if pd.notna(sample) and pd.notna(group):
                group_map.setdefault(str(sample), str(group))
    groups = [group_map.get(item, '') for item in labels]
    group_labels = groups if group_map and all(groups) else None
    method = str(score_column).replace('_r', '').capitalize() or 'Pearson'
    fig, _ = correlation_heatmap_figure(
        matrix, labels, title=style['title'] or label, method=method,
        group_labels=group_labels, colorscale=style['corr_colorscale'], cluster=True,
        mask_diagonal=True,
    )
    _apply_common_style(fig, fig.axes[0], style)
    return fig, 'data_redraw'


def _render_enrichment_frame(result_df, style, label, plot_type='enrichment_dotplot'):
    from modules.bulk_enrichment import _enrichment_figure

    if result_df.empty:
        raise FigureStudioError('富集结果表为空，无法重新绘制。')
    database = str(result_df.get('Database', pd.Series(['GO_BP'])).iloc[0])
    if plot_type != 'legacy_enrichment':
        from figure_engine import NatureFigureDirector

        director = NatureFigureDirector()
        figure_title = style['title'] or label
        # Nature templates use an English publication font.  Replace UI-only
        # Chinese labels in the title so a redraw does not produce tofu boxes
        # when a CJK font is unavailable in the export environment.
        figure_title = re.sub(r'[\u4e00-\u9fff]+', 'enrichment', str(figure_title))
        overview = plot_type == 'enrichment_overview'
        database_scope = tuple(
            value.strip() for value in str(style.get('enrichment_database_scope', '') or '')
            .replace(';', ',').replace('\n', ',').split(',') if value.strip()
        )
        spec = director.spec_from_params(
            plot_type,
            {
                'top_n': style['enrichment_overview_top_n'] if overview else style['enrichment_top_n'],
                'database_scope': database_scope,
                'pathway_selection': style['enrichment_pathway_selection'],
                'target_pathways': style['enrichment_target_pathways'],
                'gene_label_strategy': style.get('enrichment_gene_label_strategy', 'all'),
                'target_genes': style.get('enrichment_target_genes', ''),
                'max_gene_labels': style.get('enrichment_max_gene_labels', 30),
            },
            title=figure_title,
            width='double' if overview else 'single',
            formats=('svg', 'png'),
        )
        try:
            return director.render(spec, result_df), 'data_redraw'
        except (ValueError, TypeError, KeyError) as exc:
            raise FigureStudioError(str(exc)) from exc
    score_column = next((column for column in ('nes', 'NES') if column in result_df.columns), None)
    target_requested = bool(str(style.get('enrichment_target_pathways', '') or '').strip()) and style.get('enrichment_pathway_selection') in {'selected', 'selected_plus_top'}
    plot_data = result_df if target_requested else result_df.head(style['enrichment_top_n'])
    fig = _enrichment_figure(
        plot_data, title=style['title'] or label,
        database=database, score_column=score_column,
        params={
            'top_n': style['enrichment_top_n'],
            'pathway_selection': style['enrichment_pathway_selection'],
            'target_pathways': style['enrichment_target_pathways'],
        },
        ontology_colors={
            'BP': style['enrichment_bp_color'], 'MF': style['enrichment_mf_color'],
            'CC': style['enrichment_cc_color'],
            'KEGG': style['enrichment_kegg_color'], 'OTHER': style['enrichment_other_color'],
        },
    )
    if fig is None:
        raise FigureStudioError('富集结果中没有可绘制的通路。')
    # The enrichment title is a figure-level title; do not add a duplicate
    # axes title when applying the common controls.
    axes_style = {**style, 'title': ''}
    _apply_common_style(fig, fig.axes[0], axes_style)
    return fig, 'data_redraw'


def _render_bulk_enrichment(data_path, style, label, plot_type='enrichment_dotplot'):
    result_df = pd.read_csv(data_path)
    return _render_enrichment_frame(result_df, style, label, plot_type)


def _render_sc_enrichment(data_path, style, label, context, plot_type='enrichment_dotplot'):
    if plot_type == 'gsea_running':
        raise FigureStudioError(
            'GSEA 富集运行曲线无法从汇总通报表重建；请打开同一对比的 GSEA dotplot 图进行重绘。'
        )
    result_df = pd.read_csv(data_path)
    if result_df.empty:
        raise FigureStudioError('富集结果表为空，无法重新绘制。')
    for key in ('cluster', 'direction', 'method'):
        value = context.get(key)
        if value and key in result_df.columns:
            result_df = result_df[result_df[key].astype(str) == str(value)]
    if result_df.empty:
        raise FigureStudioError('未找到该富集图对应的 cluster 或方向数据。')
    return _render_enrichment_frame(result_df, style, label, plot_type)


def _font_for_pillow(size):
    from PIL import ImageFont
    try:
        from matplotlib import font_manager
        return ImageFont.truetype(font_manager.findfont('DejaVu Sans'), size)
    except Exception:
        return ImageFont.load_default()


def _render_raster_style_only(source_path, style):
    from PIL import Image, ImageDraw

    with Image.open(source_path) as image:
        image = image.convert('RGBA')
    title = style['title']
    margin = max(24, round(image.width * 0.035))
    title_height = max(0, int(style['font_size'] * 2.3)) if title else 0
    legend_height = max(0, int(style['font_size'] * 2.0)) if style['legend_text'] else 0
    canvas = Image.new('RGBA', (image.width + margin * 2, image.height + margin * 2 + title_height + legend_height), style['background'])
    draw = ImageDraw.Draw(canvas)
    if title:
        draw.text((margin, margin), title, fill='#111827', font=_font_for_pillow(int(style['font_size'] * 1.4)))
    image_y = margin + title_height
    canvas.alpha_composite(image, (margin, image_y))
    if style['legend_text']:
        legend_y = image_y + image.height + max(4, margin // 3)
        draw.text((margin, legend_y), style['legend_text'], fill='#475467', font=_font_for_pillow(int(style['font_size'])))
    return canvas.convert('RGB')


def _write_raster_svg(image, svg_path, style):
    buffer = io.BytesIO()
    image.save(buffer, format='PNG')
    encoded = base64.b64encode(buffer.getvalue()).decode('ascii')
    width, height = image.size
    Path(svg_path).write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
        f'<image width="{width}" height="{height}" href="data:image/png;base64,{encoded}"/>'
        '</svg>', encoding='utf-8',
    )


def _render_svg_style_only(source_path, svg_path, style):
    raw = Path(source_path).read_text(encoding='utf-8', errors='ignore')
    if _SAFE_SVG_RE.search(raw):
        raise FigureStudioError('SVG 含有不允许的脚本或外部资源，无法编辑。')
    try:
        root = ET.fromstring(raw)
    except ET.ParseError as exc:
        raise FigureStudioError(f'SVG 文件格式无效：{exc}') from exc
    view_box = root.attrib.get('viewBox', '0 0 1200 800').split()
    try:
        width, height = float(view_box[-2]), float(view_box[-1])
    except (ValueError, IndexError):
        width, height = 1200.0, 800.0
    title_space = 70 if style['title'] else 24
    inner = ''.join(ET.tostring(child, encoding='unicode') for child in root)
    title = f'<text x="32" y="46" font-family="{xml_escape(style["font_family"])}" font-size="{style["font_size"] + 7}" font-weight="600" fill="#111827">{xml_escape(style["title"])}</text>' if style['title'] else ''
    wrapper = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height + title_space}" '
        f'viewBox="0 0 {width} {height + title_space}"><rect width="100%" height="100%" fill="{style["background"]}"/>'
        f'{title}<g transform="translate(0 {title_space})">{inner}</g></svg>'
    )
    Path(svg_path).write_text(wrapper, encoding='utf-8')


def render_source(source, style, png_path=None, svg_path=None):
    """Render an immutable edited version and return the effective editing mode."""
    style = normalize_style(style)
    if source['edit_mode'] == 'bulk_volcano':
        fig, mode = _render_bulk_volcano(source['data_path'], style, source['label'])
    elif source['edit_mode'] == 'sc_volcano':
        fig, mode = _render_sc_volcano(
            source['data_path'], style, source['label'],
            context=source.get('sc_context') or {},
        )
    elif source['edit_mode'] == 'bulk_ma':
        fig, mode = _render_bulk_ma(source['data_path'], style, source['label'])
    elif source['edit_mode'] == 'sc_ma':
        fig, mode = _render_sc_ma(
            source['data_path'], style, source['label'],
            context=source.get('sc_context') or {},
        )
    elif source['edit_mode'] == 'bulk_heatmap':
        fig, mode = _render_bulk_heatmap(
            source['data_path'], style, source['label'],
            comparison_label=source.get('heatmap_comparison_label', ''),
            default_groupby=source.get('heatmap_groupby', ''),
        )
    elif source['edit_mode'] == 'bulk_correlation':
        fig, mode = _render_bulk_correlation(source['data_path'], style, source['label'])
    elif source['edit_mode'] in {'bulk_enrichment', 'bulk_enrichment_overview'}:
        fig, mode = _render_bulk_enrichment(
            source['data_path'], style, source['label'],
            plot_type=source.get('plot_type', 'enrichment_dotplot'),
        )
    elif source['edit_mode'] == 'sc_enrichment':
        fig, mode = _render_sc_enrichment(
            source['data_path'], style, source['label'],
            context=source.get('sc_context') or {},
            plot_type=source.get('plot_type', 'enrichment_dotplot'),
        )
    else:
        fig = None
        mode = 'style_only'

    if fig is not None:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        mpl.rcParams['svg.fonttype'] = 'none'
        mpl.rcParams['pdf.fonttype'] = 42
        if png_path:
            fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight', pad_inches=0.15, facecolor=style['background'])
        if svg_path:
            fig.savefig(svg_path, format='svg', bbox_inches='tight', pad_inches=0.15, facecolor=style['background'])
        plt.close(fig)
        return mode, style

    extension = Path(source['file_path']).suffix.lower()
    if extension in RASTER_EXTENSIONS:
        image = _render_raster_style_only(source['file_path'], style)
        if png_path:
            image.save(png_path, format='PNG', optimize=True)
        if svg_path:
            _write_raster_svg(image, svg_path, style)
        return mode, style
    if extension == '.svg':
        if svg_path:
            _render_svg_style_only(source['file_path'], svg_path, style)
        return mode, style
    raise FigureStudioError('不支持该图片格式。')


def render_preview_data_uri(source, style):
    """Render a temporary preview as a browser-safe data URI."""
    with tempfile.TemporaryDirectory(dir=Config.runtime_tmp_dir(), prefix='figure-preview-') as temp_dir:
        png_path = os.path.join(temp_dir, 'preview.png')
        svg_path = os.path.join(temp_dir, 'preview.svg')
        mode, normalized = render_source(source, style, png_path=png_path, svg_path=svg_path)
        selected = png_path if os.path.isfile(png_path) else svg_path
        if not os.path.isfile(selected):
            raise FigureStudioError('预览渲染未生成文件。')
        mime = 'image/png' if selected.endswith('.png') else 'image/svg+xml'
        encoded = base64.b64encode(Path(selected).read_bytes()).decode('ascii')
    return {'data_uri': f'data:{mime};base64,{encoded}', 'edit_mode': mode, 'style': normalized}


def save_figure_version(project_id, source, style, label=''):
    """Render a non-destructive PNG/SVG version and persist its provenance."""
    output_dir = Path(Config.project_dir(project_id)) / 'figure_studio' / 'versions'
    output_dir.mkdir(parents=True, exist_ok=True)
    version_id = uuid4().hex[:12]
    png_path = output_dir / f'{version_id}.png'
    svg_path = output_dir / f'{version_id}.svg'
    mode, normalized = render_source(source, style, png_path=str(png_path), svg_path=str(svg_path))
    # SVG-only source images intentionally have no generated PNG.  Keep that
    # absence as an empty database value rather than the current directory.
    png_output = str(png_path) if png_path.is_file() else ''
    svg_output = str(svg_path) if svg_path.is_file() else ''
    version = FigureVersion(
        id=version_id, project_id=project_id, source_kind=source['kind'], source_id=source['id'],
        # Persist the concrete capability (rather than the generic render
        # outcome) so a saved volcano/heatmap version remains editable later.
        edit_mode=(source.get('edit_mode') if source.get('edit_mode') in {
            'bulk_volcano', 'sc_volcano', 'bulk_ma', 'sc_ma', 'bulk_heatmap', 'bulk_correlation',
            'bulk_enrichment', 'bulk_enrichment_overview', 'sc_enrichment'
        } else mode),
        label=_clean_text(label) or _clean_text(normalized.get('title')) or source['label'],
        style_json=json.dumps(normalized, ensure_ascii=False),
        png_path=png_output, svg_path=svg_output,
    )
    version.save()
    return version


def validate_uploaded_image(file_storage):
    """Validate size/type before storing an uploaded image in a project folder."""
    filename = Path(file_storage.filename or '').name
    extension = Path(filename).suffix.lower()
    if extension not in {'.png', '.jpg', '.jpeg', '.svg'}:
        raise FigureStudioError('仅支持 PNG、JPG、JPEG 或 SVG 图片。')
    file_storage.stream.seek(0, os.SEEK_END)
    size = file_storage.stream.tell()
    file_storage.stream.seek(0)
    if size <= 0 or size > MAX_UPLOAD_BYTES:
        raise FigureStudioError('图片不能为空，且大小不能超过 15 MB。')
    content = file_storage.stream.read()
    file_storage.stream.seek(0)
    if extension == '.svg':
        text = content.decode('utf-8', errors='ignore')
        if _SAFE_SVG_RE.search(text) or '<svg' not in text.lower():
            raise FigureStudioError('SVG 含有不允许内容或格式无效。')
    else:
        from PIL import Image
        try:
            with Image.open(io.BytesIO(content)) as image:
                image.verify()
        except Exception as exc:
            raise FigureStudioError('上传文件不是有效图片。') from exc
    return extension
