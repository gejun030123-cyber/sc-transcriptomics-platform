import os
import json
import re
import textwrap
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


def _ensure_gsea_term_column(result_df):
    """Preserve pathway names stored in a GSEA result index."""
    result = result_df.copy()
    if 'Term' not in result.columns:
        index_name = result.index.name or 'index'
        result = result.reset_index().rename(columns={index_name: 'Term'})
    return result


ONTOLOGY_COLORS = {
    'BP': '#FDBE85',
    'MF': '#B8A9D1',
    'KEGG': '#7BC77B',
    'OTHER': '#9CB8D8',
}


def _safe_output_fragment(value, fallback='unspecified'):
    """Create a readable, filesystem-safe identifier without losing contrast identity."""
    text = re.sub(r'[^A-Za-z0-9]+', '_', str(value or '').strip()).strip('_').lower()
    return text[:96] or fallback


def _enrichment_output_key(method, database, direction=None, comparison=None):
    """Build a stable, non-overwriting stem for one enrichment result."""
    parts = ['enrichment', str(method).lower(), str(database).lower()]
    if comparison:
        parts.append(_safe_output_fragment(comparison))
    if direction:
        parts.append(str(direction).lower())
    return '_'.join(parts)


def _add_enrichment_metadata(result_df, method, database, direction=None, comparison=None):
    """Annotate a result table so different databases can be safely combined."""
    result = result_df.copy()
    metadata = {
        'Comparison': comparison or 'Unspecified',
        'Database': database,
        'Method': method,
    }
    for column, value in reversed(tuple(metadata.items())):
        if column in result.columns:
            result[column] = value
        else:
            result.insert(0, column, value)
    if direction is not None:
        result['Direction'] = direction
    elif 'Direction' not in result.columns:
        # Split-direction ORA stores the value as lower-case ``direction``
        # before the combined table is annotated; retain it rather than
        # collapsing Up and Down into an ambiguous "All" group.
        if 'direction' in result.columns:
            result['Direction'] = result['direction'].fillna('All')
        else:
            result['Direction'] = 'All'
    return result


def _enrichment_term_column(df):
    """Return the available pathway-term column, retaining index values if needed."""
    column = next((c for c in ('Term', 'Description', 'term', 'description') if c in df.columns), None)
    if column:
        return df, column
    result = df.reset_index().rename(columns={df.index.name or 'index': 'Term'})
    return result, 'Term'


def _comparison_label_from_deg_record(result_file, task):
    """Recover the selected DEG contrast without importing the web route."""
    label = str(getattr(result_file, 'label', '') or '').strip()
    match = re.search(r'\(([^()]*\bvs\b[^()]*)\)', label, flags=re.IGNORECASE)
    if match:
        return match.group(1).strip()

    filename = os.path.basename(str(getattr(result_file, 'file_path', '') or ''))
    suffix_match = re.search(r'bulk_deg_results_(\d+)\.csv$', filename, flags=re.IGNORECASE)
    if suffix_match:
        try:
            task_result = json.loads(getattr(task, 'result_json', '') or '{}')
            comparisons = list(task_result.get('comparisons') or [])
            index = int(suffix_match.group(1))
            if 0 <= index < len(comparisons):
                return str(comparisons[index]).strip()
        except (TypeError, ValueError, json.JSONDecodeError):
            pass
    return label or filename


def _backfill_legacy_enrichment_metadata(project_dir, results_dir):
    """Add provenance columns to registered enrichment tables from older runs.

    Earlier platform versions saved a generic ``enrichment_ora_results.csv``
    without the selected contrast or database.  The task record still contains
    that provenance, so recover it once before building the project-level
    integration.  Only registered, project-local tables are touched; unknown
    files remain unchanged rather than being guessed into a comparison.
    """
    try:
        from config import Config
        from models import AnalysisTask, ResultFile

        project_real = os.path.realpath(project_dir)
        project_id = Path(project_real).name
        if os.path.realpath(Config.project_dir(project_id)) != project_real:
            return []
    except (ImportError, OSError, ValueError):
        return []

    results_real = os.path.realpath(results_dir)
    tasks = AnalysisTask.get_by_project(project_id, include_branches=True)
    deg_sources = {}
    for task in tasks:
        if task.module_name != 'bulk_deg' or task.status != 'completed':
            continue
        for result_file in ResultFile.get_by_task(task.id):
            filename = os.path.basename(result_file.file_path).lower()
            if not filename.startswith('bulk_deg_results') or not filename.endswith('.csv'):
                continue
            deg_sources[os.path.realpath(result_file.file_path)] = _comparison_label_from_deg_record(
                result_file, task,
            )

    updated = []
    for task in tasks:
        if task.module_name != 'bulk_enrichment' or task.status != 'completed':
            continue
        try:
            params = json.loads(task.params_json or '{}')
        except (TypeError, ValueError, json.JSONDecodeError):
            continue
        comparison = str(params.get('input_comparison', '') or '').strip()
        if not comparison:
            input_source = str(params.get('input_source', '') or '')
            comparison = deg_sources.get(os.path.realpath(input_source), '') if input_source else ''
        if not comparison and str(params.get('custom_genes', '') or '').strip():
            comparison = 'Custom genes'
        if not comparison:
            continue

        method = str(params.get('method', 'ORA') or 'ORA')
        database = str(params.get('database', 'GO_BP') or 'GO_BP')
        for result_file in ResultFile.get_by_task(task.id):
            path = Path(result_file.file_path)
            filename = path.name
            if (
                result_file.file_type != 'csv'
                or not filename.startswith('enrichment_')
                or not filename.endswith('_results.csv')
                or filename == 'enrichment_integrated_results.csv'
                or not path.is_file()
                or os.path.realpath(path.parent) != results_real
            ):
                continue
            try:
                table = pd.read_csv(path)
            except (OSError, pd.errors.EmptyDataError):
                continue
            if table.empty or {'Comparison', 'Database', 'Method', 'Direction'}.issubset(table.columns):
                continue
            table = _add_enrichment_metadata(
                table, method=method, database=database, comparison=comparison,
            )
            table.to_csv(path, index=False)
            updated.append(str(path))
    return updated


def _integrated_enrichment_overview(group, comparison, method, direction):
    """Draw one multi-panel overview for one contrast across databases.

    Different contrasts, methods and Up/Down ORA subsets are deliberately not
    mixed: their biological question or score scale differs.  Within one group,
    every database receives its own panel so long pathway labels remain legible.
    """
    import matplotlib.pyplot as plt

    databases = [str(value) for value in group['Database'].dropna().unique().tolist()]
    if len(databases) < 2:
        return None
    is_gsea = str(method).upper() == 'GSEA'
    ncols = min(2, len(databases))
    nrows = int(np.ceil(len(databases) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.9 * ncols, 4.2 * nrows),
                             dpi=150, squeeze=False, sharex=True)
    flat_axes = axes.ravel()
    panel_colors = [
        ONTOLOGY_COLORS.get(_enrichment_ontology(database, database), ONTOLOGY_COLORS['OTHER'])
        for database in databases
    ]
    for axis, database, color in zip(flat_axes, databases, panel_colors):
        subset = group[group['Database'].astype(str) == database].copy()
        if is_gsea:
            score_column = next((column for column in ('NES', 'nes', 'ES', 'es') if column in subset.columns), None)
            if score_column is None:
                axis.set_visible(False)
                continue
            subset['_overview_score'] = pd.to_numeric(subset[score_column], errors='coerce')
            subset['_overview_rank'] = subset['_overview_score'].abs()
            x_label = 'Normalized enrichment score'
        else:
            subset['_overview_score'] = -np.log10(
                pd.to_numeric(subset['Enrichment FDR'], errors='coerce').clip(lower=1e-300)
            )
            subset['_overview_rank'] = subset['_overview_score']
            x_label = '-log10(enrichment FDR)'
        subset = subset.dropna(subset=['_overview_score']).sort_values(
            '_overview_rank', ascending=False,
        ).head(8)
        if subset.empty:
            axis.set_visible(False)
            continue
        terms = [textwrap.fill(str(term), width=34, break_long_words=False)
                 for term in subset['Term'].tolist()]
        values = subset['_overview_score'].to_numpy(dtype=float)
        positions = np.arange(len(subset))
        axis.barh(positions, values, color=color, edgecolor='none', alpha=0.92)
        axis.set_yticks(positions, terms)
        axis.invert_yaxis()
        if is_gsea and np.nanmin(values) < 0 < np.nanmax(values):
            axis.axvline(0, color='#98A2B3', linewidth=0.75, zorder=0)
        axis.set_title(database, loc='left', fontsize=10, fontweight='semibold', color='#111827')
        axis.set_xlabel(x_label, fontsize=8.5, color='#344054')
        axis.tick_params(axis='y', labelsize=7.2, length=0, pad=5, colors='#1F2937')
        axis.tick_params(axis='x', labelsize=7.5, colors='#475467')
        axis.grid(axis='x', color='#D8DEE9', linewidth=0.55, alpha=0.65, zorder=0)
        for spine in axis.spines.values():
            spine.set_visible(False)
    for axis in flat_axes[len(databases):]:
        axis.set_visible(False)
    direction_text = '' if str(direction) == 'All' else f' · {direction}'
    fig.suptitle(f'{comparison} · {method}{direction_text} pathway enrichment overview',
                 x=0.02, ha='left', fontsize=14, fontweight='semibold', color='#111827')
    fig.tight_layout(rect=(0, 0, 1, 0.93), pad=1.25)
    return fig


def _write_integrated_overview_figures(integrated, results_path):
    """Export one overview image for each comparable contrast/database group."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    plots_dir = Path(results_path).parent / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    group_columns = ['Comparison', 'Method', 'Direction']
    for (comparison, method, direction), group in integrated.groupby(group_columns, dropna=False, sort=True):
        figure = _integrated_enrichment_overview(group, str(comparison), str(method), str(direction))
        if figure is None:
            continue
        output_key = '_'.join([
            'enrichment_overview', _safe_output_fragment(method),
            _safe_output_fragment(comparison), _safe_output_fragment(direction),
        ])
        png_path = plots_dir / f'{output_key}.png'
        svg_path = plots_dir / f'{output_key}.svg'
        mpl.rcParams['svg.fonttype'] = 'none'
        figure.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.12, facecolor='white')
        figure.savefig(svg_path, format='svg', bbox_inches='tight', pad_inches=0.12, facecolor='white')
        plt.close(figure)
        label = f'{comparison} · {method}' + ('' if str(direction) == 'All' else f' · {direction}')
        outputs.extend([
            {'file_path': str(png_path), 'file_type': 'png', 'category': 'enrichment_overview',
             'label': f'{label} 多数据库富集概览'},
            {'file_path': str(svg_path), 'file_type': 'svg', 'category': 'enrichment_overview',
             'label': f'{label} 多数据库富集概览'},
        ])
    return outputs


def _write_enrichment_integration(results_dir):
    """Create project-level CSV/XLSX tables from all database-specific results.

    Each source table carries ``Comparison``, ``Database``, ``Method`` and
    ``Direction``.  The integration intentionally retains those fields instead
    of merging similarly named pathways across libraries, whose gene-set
    definitions may differ.
    """
    results_path = Path(results_dir)
    source_files = sorted(
        path for path in results_path.glob('enrichment_*_results.csv')
        if path.name not in {'enrichment_integrated_results.csv'}
    )
    tables = []
    for path in source_files:
        try:
            table = pd.read_csv(path)
        except (OSError, pd.errors.EmptyDataError):
            continue
        if table.empty or not {'Database', 'Method', 'Direction'}.issubset(table.columns):
            # Legacy generic files lack source metadata and cannot be reliably
            # assigned to a database after they have been overwritten.
            continue
        table, term_col = _enrichment_term_column(table)
        if 'Comparison' not in table.columns:
            table.insert(0, 'Comparison', 'Unspecified')
        if term_col != 'Term':
            table['Term'] = table[term_col]
        fdr_col = next((c for c in ('Adjusted P-value', 'Adjusted p-value', 'FDR', 'fdr', 'padj') if c in table.columns), None)
        pvalue_col = next((c for c in ('P-value', 'pvalue', 'pval') if c in table.columns), None)
        if 'Enrichment FDR' not in table.columns:
            table['Enrichment FDR'] = pd.to_numeric(table[fdr_col], errors='coerce') if fdr_col else np.nan
        if 'Enrichment P-value' not in table.columns:
            table['Enrichment P-value'] = pd.to_numeric(table[pvalue_col], errors='coerce') if pvalue_col else np.nan
        table['Source result'] = path.name
        tables.append(table)

    if not tables:
        return []

    integrated = pd.concat(tables, ignore_index=True, sort=False)
    preferred = [
        'Comparison', 'Database', 'Method', 'Direction', 'Term', 'Enrichment FDR',
        'Enrichment P-value', 'NES', 'Overlap', 'Genes', 'Source result',
    ]
    ordered_columns = [column for column in preferred if column in integrated.columns]
    ordered_columns.extend(column for column in integrated.columns if column not in ordered_columns)
    integrated = integrated[ordered_columns].sort_values(
        ['Enrichment FDR', 'Comparison', 'Database', 'Direction', 'Term'],
        na_position='last', kind='stable',
    )

    output_files = []
    csv_path = results_path / 'enrichment_integrated_results.csv'
    integrated.to_csv(csv_path, index=False)
    output_files.append({
        'file_path': str(csv_path), 'file_type': 'csv', 'category': 'table',
        'label': f'通路富集整合表（{len(integrated)} 条）',
    })
    output_files.extend(_write_integrated_overview_figures(integrated, results_path))

    try:
        xlsx_path = results_path / 'enrichment_integrated_results.xlsx'
        summary = (integrated.groupby(['Comparison', 'Database', 'Method', 'Direction'], dropna=False)
                   .size().reset_index(name='n_terms'))
        with pd.ExcelWriter(xlsx_path) as writer:
            integrated.to_excel(writer, sheet_name='Integrated', index=False)
            summary.to_excel(writer, sheet_name='Summary', index=False)
        output_files.append({
            'file_path': str(xlsx_path), 'file_type': 'xlsx', 'category': 'table',
            'label': '通路富集整合表（Excel）',
        })
    except (ImportError, ModuleNotFoundError, ValueError):
        # CSV is always available; Excel export remains an optional convenience.
        pass
    return output_files


def _audit_figure_text_overlap(fig):
    """Detect overlapping visible text after Matplotlib has laid out a figure."""
    from matplotlib.text import Text

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    figure_bounds = fig.bbox
    texts = []
    for artist in fig.findobj(match=Text):
        label = artist.get_text().strip()
        if not label or not artist.get_visible():
            continue
        try:
            bbox = artist.get_window_extent(renderer)
        except Exception:
            continue
        if bbox.width <= 0 or bbox.height <= 0:
            continue
        texts.append((label, bbox))

    overlap_pairs = []
    for index, (left_label, left_box) in enumerate(texts):
        for right_label, right_box in texts[index + 1:]:
            overlap_width = min(left_box.x1, right_box.x1) - max(left_box.x0, right_box.x0)
            overlap_height = min(left_box.y1, right_box.y1) - max(left_box.y0, right_box.y0)
            if overlap_width > 1 and overlap_height > 1:
                overlap_pairs.append({
                    'left': left_label[:120], 'right': right_label[:120],
                    'overlap_area_px2': round(float(overlap_width * overlap_height), 1),
                })

    out_of_bounds = [
        label[:120] for label, bbox in texts
        if bbox.x0 < figure_bounds.x0 - 1 or bbox.x1 > figure_bounds.x1 + 1
        or bbox.y0 < figure_bounds.y0 - 1 or bbox.y1 > figure_bounds.y1 + 1
    ]
    return {
        'status': 'pass' if not overlap_pairs and not out_of_bounds else 'warning',
        'n_text_elements': len(texts),
        'n_overlap_pairs': len(overlap_pairs),
        'overlap_pairs': overlap_pairs,
        'out_of_bounds_labels': out_of_bounds,
    }


def _write_figure_layout_report(results_dir, output_key, audit):
    """Persist a machine-readable layout QA report alongside enrichment tables."""
    report_path = Path(results_dir) / f'{output_key}_figure_layout_qc.json'
    report_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    return {
        'file_path': str(report_path), 'file_type': 'json', 'category': 'quality_control',
        'label': f'富集图文字重叠检测（{audit["status"]}）',
    }


def _save_enrichment_figure(analysis, fig, plots_dir, results_dir, output_key, label):
    """Export one enrichment plot and its post-layout overlap report."""
    fig._post_layout_audit = _audit_figure_text_overlap
    exported = analysis.save_matplotlib_figure(
        fig, plots_dir, f'{output_key}.png', 'enrichment', label,
        preserve_aspect=True,
    )
    audit = getattr(fig, '_post_layout_audit_result', None)
    if audit is None:
        audit = _audit_figure_text_overlap(fig)
    audit['layout_contract'] = getattr(fig, '_enrichment_layout_contract', {})
    exported.append(_write_figure_layout_report(results_dir, output_key, audit))
    return exported, audit


# OmicVerse's ``download_pathway_database`` intentionally ships only the
# databases listed below.  KEGG is not among them, yet the enrichment UI offers
# KEGG as an option.  Fetch a missing requested library directly from Enrichr,
# which is also the source format used by the bundled gene-set files.
_ENRICHR_GENESET_URL = 'https://maayanlab.cloud/Enrichr/geneSetLibrary'
_GENESET_ALIASES = {
    # OmicVerse uses the singular 2021 Human filename, while the historical UI
    # option used the plural 2019 name.
    'WikiPathways_2019_Human': 'WikiPathway_2021_Human',
}


def _download_enrichr_geneset(library_name, destination):
    """Download one Enrichr GMT-like library atomically, with basic validation."""
    query = urlencode({'mode': 'text', 'libraryName': library_name})
    request = Request(
        f'{_ENRICHR_GENESET_URL}?{query}',
        headers={'User-Agent': 'sc-transcriptomics-platform/1.0'},
    )
    try:
        with urlopen(request, timeout=30) as response:
            content = response.read()
    except Exception as exc:
        raise RuntimeError(
            f'无法下载基因集数据库 {library_name}。请检查服务器网络后重试：{exc}'
        ) from exc

    if not content or b'\t' not in content or content.lstrip().startswith(b'<'):
        raise RuntimeError(f'下载的基因集数据库 {library_name} 格式无效。')

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile('wb', dir=destination.parent, delete=False,
                            prefix=f'.{destination.name}.', suffix='.tmp') as tmp:
        tmp.write(content)
        tmp_path = Path(tmp.name)
    try:
        tmp_path.replace(destination)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _resolve_geneset_path(db_filename, organism):
    """Return an available local gene-set file, downloading only when needed."""
    candidate_names = [db_filename]
    alias = _GENESET_ALIASES.get(db_filename)
    if alias:
        candidate_names.append(alias)

    organism_lower = organism.lower()
    for name in candidate_names:
        candidate_paths = [
            Path('genesets') / f'{name}.txt',
            Path('genesets') / f'{name}_{organism}.txt',
            Path('genesets') / f'{name}_{organism_lower}.txt',
        ]
        for path in candidate_paths:
            if path.exists():
                return str(path)

    destination = Path('genesets') / f'{db_filename}.txt'
    _download_enrichr_geneset(db_filename, destination)
    return str(destination)


def _enrichment_ontology(value, database=''):
    """Normalize database labels to the compact ontology legend used in plots."""
    # Enrichment outputs sometimes set ``Gene_set`` to an internal label such
    # as ``gs_ind``.  Include the requested database as authoritative context
    # so a KEGG-only result is not shown as the generic "OTHER" category.
    text = f'{value or ""} {database or ""}'.upper()
    if 'KEGG' in text:
        return 'KEGG'
    if 'MF' in text or 'MOLECULAR FUNCTION' in text:
        return 'MF'
    if 'BP' in text or 'BIOLOGICAL PROCESS' in text:
        return 'BP'
    return 'OTHER'


def _enrichment_count(value):
    """Extract the leading overlap count from values such as ``9/132``."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return 1
    text = str(value).strip()
    try:
        return max(1, int(float(text.split('/', 1)[0])))
    except (TypeError, ValueError):
        return 1


def _enrichment_figure(result_df, title, database='GO_BP', score_column=None,
                       ontology_colors=None):
    """Create a dense but readable horizontal enrichment chart.

    The figure uses wrapped y-axis labels rather than text laid over bars.  This
    keeps long pathway names, the overlap count and the external legend from
    colliding when a user requests many terms.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    df = result_df.copy()
    if df.empty:
        return None
    term_col = next((c for c in ('Term', 'Description', 'term', 'description') if c in df.columns), None)
    if term_col is None:
        df = df.reset_index().rename(columns={df.index.name or 'index': 'Term'})
        term_col = 'Term'
    p_col = next((c for c in ('Adjusted P-value', 'Adjusted p-value', 'p.adjust', 'FDR', 'fdr', 'P-value', 'pvalue') if c in df.columns), None)
    if score_column and score_column in df.columns:
        score = pd.to_numeric(df[score_column], errors='coerce').abs()
        x_label = 'Enrichment score'
    elif p_col:
        score = -np.log10(pd.to_numeric(df[p_col], errors='coerce').clip(lower=1e-300))
        x_label = '-log10(adjusted P-value)'
    else:
        score = pd.Series(np.arange(len(df), 0, -1), index=df.index, dtype=float)
        x_label = 'Enrichment score'
    df = df.assign(_score=score).dropna(subset=['_score'])
    if df.empty:
        return None

    ontology_col = next((c for c in ('Ontology', 'ontology', 'Gene_set', 'gene_set', 'database') if c in df.columns), None)
    df['_ontology'] = [
        _enrichment_ontology(row[ontology_col] if ontology_col else database, database)
        for _, row in df.iterrows()
    ]
    # Keep ontology blocks together, then order the most significant pathways
    # first within each block as in the supplied reference figure.
    ontology_order = {'BP': 0, 'MF': 1, 'KEGG': 2, 'OTHER': 3}
    df['_ontology_rank'] = df['_ontology'].map(ontology_order).fillna(3)
    df = df.sort_values(['_ontology_rank', '_score'], ascending=[True, False])
    ontologies = df['_ontology'].tolist()
    count_col = next((c for c in ('Overlap', 'Count', 'count', 'Gene Count', 'gene_count', 'setSize', 'size') if c in df.columns), None)
    counts = [_enrichment_count(row[count_col]) if count_col else 1 for _, row in df.iterrows()]
    terms = [str(value) for value in df[term_col].tolist()]
    wrapped_terms = [textwrap.fill(term, width=42, break_long_words=False) for term in terms]
    line_counts = np.array([label.count('\n') + 1 for label in wrapped_terms], dtype=float)
    values = df['_score'].astype(float).to_numpy()
    row_steps = np.maximum(1.0, line_counts) * 1.16
    y = np.cumsum(row_steps) - row_steps / 2
    palette = {**ONTOLOGY_COLORS, **(ontology_colors or {})}
    colors = [palette.get(item, palette['OTHER']) for item in ontologies]

    height = max(4.8, min(16.0, 1.45 + 0.53 * float(line_counts.sum())))
    fig, ax = plt.subplots(figsize=(11.2, height), dpi=150)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.barh(y, values, color=colors, edgecolor='none', height=row_steps * 0.69,
            alpha=0.96, zorder=1)
    ax.set_yticks(y, wrapped_terms)
    ax.invert_yaxis()
    ax.set_xlabel(x_label, fontsize=12, color='#142a8b', labelpad=12)
    ax.set_ylabel('Pathway', fontsize=12, color='#142a8b', labelpad=18)
    fig.suptitle(title, x=0.02, y=0.985, ha='left', va='top',
                 fontsize=17, fontweight='semibold', color='#111827')
    ax.grid(axis='x', color='#d8dee9', linewidth=0.65, alpha=0.55, zorder=0)
    ax.tick_params(axis='y', labelsize=9, length=0, pad=9, colors='#111827')
    ax.tick_params(axis='x', labelsize=10, colors='#374151', width=0.6)
    for spine in ax.spines.values():
        spine.set_visible(False)

    max_score = max(float(values.max()), 1.0)
    count_x = max_score * 1.055
    ax.set_xlim(0, max_score * 1.20)
    # Do not allow the automatic locator to draw an extra terminal tick outside
    # the visible score range (for example an ``8`` tick beyond a 7.2 limit).
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6, prune='upper'))
    for yi, value, count in zip(y, values, counts):
        ax.text(min(value + max_score * 0.015, count_x), yi, f'n={count}',
                ha='left', va='center', fontsize=8.6, color='#374151', zorder=3)
    ax.text(count_x, y[0] - row_steps[0] * 0.88, 'overlap', ha='left', va='bottom',
            fontsize=8.4, color='#6b7280')

    present_ontologies = [item for item in ('BP', 'MF', 'KEGG', 'OTHER') if item in ontologies]
    if len(present_ontologies) > 1:
        ontology_handles = [Patch(facecolor=palette[item], edgecolor='none', label=item)
                            for item in present_ontologies]
        fig.legend(handles=ontology_handles, title='DATABASE', loc='upper left',
                   bbox_to_anchor=(0.84, 0.89), bbox_transform=fig.transFigure,
                   frameon=False, fontsize=9, title_fontsize=10)
        layout_right = 0.80
    else:
        layout_right = 0.98
    fig._native_layout_rect = (0.0, 0.0, layout_right, 0.94)
    fig._enrichment_layout_contract = {
        'wrapped_pathway_labels': True,
        'max_label_lines': int(line_counts.max()),
        'requested_terms': int(len(df)),
    }
    return fig


class BulkEnrichmentAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_enrichment"
    DISPLAY_NAME = "通路富集分析"
    DESCRIPTION = "GO/KEGG/WikiPathways 通路富集分析（ORA / GSEA）"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import omicverse as ov

        method = self.params.get('method', 'ORA')
        database = self.params.get('database', 'GO_BP')
        organism = self.params.get('organism', 'Human')
        pvalue_cutoff = float(self.params.get('pvalue_cutoff', 0.05))
        top_n = int(self.params.get('top_n', 20))
        input_source = self.params.get('input_source', '')
        split_direction = self.params.get('split_direction', False) in (True, 'true', 'on', '1')
        custom_genes_str = self.params.get('custom_genes', '').strip()
        comparison = str(self.params.get('input_comparison', '')).strip()
        if not comparison:
            comparison = ('Custom genes' if custom_genes_str else
                          (Path(input_source).stem if input_source else 'Unspecified'))

        self.progress(5, "加载基因列表...")

        # Validate input_source path
        if input_source:
            real_input = os.path.realpath(input_source)
            real_project = os.path.realpath(self.project_dir)
            if not real_input.startswith(real_project + os.sep) and real_input != real_project:
                raise ValueError("input_source 必须在项目目录内")

        # Load gene list from DEG CSV result
        deg_genes = []
        gene_rnk = None

        if custom_genes_str:
            deg_genes = [g.strip() for g in custom_genes_str.replace('\n', ',').split(',') if g.strip()]
            gene_rnk = None
        elif input_source and os.path.exists(input_source):
            deg_df = pd.read_csv(input_source)
            if 'regulation' in deg_df.columns:
                deg_genes = deg_df[deg_df['regulation'] != 'NS']['gene'].tolist()
                if 'log2FC' in deg_df.columns:
                    gene_rnk = deg_df[['gene', 'log2FC']].dropna().sort_values('log2FC', ascending=False)
                    gene_rnk.columns = ['gene_name', 'rank']
            elif 'gene' in deg_df.columns:
                deg_genes = deg_df['gene'].tolist()

        if not deg_genes and method == 'ORA':
            raise ValueError("未找到差异基因列表。请先运行 DEG 分析，并在 input_source 中指定结果 CSV 文件路径。")

        self.progress(15, f"加载 {database} 基因集数据库...")

        # Map database name to file
        db_map = {
            'GO_BP': 'GO_Biological_Process_2021',
            'GO_MF': 'GO_Molecular_Function_2021',
            'GO_CC': 'GO_Cellular_Component_2021',
            'KEGG': 'KEGG_2021_Human' if organism == 'Human' else 'KEGG_2019_Mouse',
            'WikiPathways': 'WikiPathways_2019_Human' if organism == 'Human' else 'WikiPathways_2019_Mouse',
            'Reactome': 'Reactome_2022',
        }

        db_filename = db_map.get(database, 'GO_Biological_Process_2021')
        organism_lower = organism.lower()

        # Resolve only the selected library.  The generic OmicVerse downloader
        # does not include KEGG and used to make a valid KEGG request fail here.
        db_path = _resolve_geneset_path(db_filename, organism)

        pathways_dict = ov.utils.geneset_prepare(db_path, organism=organism)

        self.progress(30, f"运行 {method} 富集分析...")

        result_files = []
        plots_dir = os.path.join(self.project_dir, 'plots')
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        n_sig = 0
        figure_audits = []

        if method == 'ORA':
            if split_direction and not custom_genes_str and input_source and os.path.exists(input_source):
                # Split direction: run ORA separately for Up and Down genes
                deg_df_split = pd.read_csv(input_source)
                all_split_results = []

                for direction in ['Up', 'Down']:
                    direction_genes = deg_df_split[deg_df_split['regulation'] == direction]['gene'].tolist()
                    if not direction_genes:
                        continue

                    self.progress(45 if direction == 'Up' else 60, f"运行 {direction} 基因 ORA...")
                    enr_dir = ov.bulk.geneset_enrichment(
                        gene_list=direction_genes,
                        pathways_dict=pathways_dict,
                        pvalue_type='adjust',
                        pvalue_threshold=pvalue_cutoff,
                        organism=organism_lower,
                        outdir=os.path.join(self.project_dir, f'enrichr_{direction.lower()}_tmp')
                    )

                    enr_dir['direction'] = direction
                    all_split_results.append(enr_dir)

                    output_key = _enrichment_output_key(
                        method, database, direction, comparison=comparison)
                    top_enr_dir = enr_dir.head(top_n)
                    if len(top_enr_dir) > 0:
                        fig_enrichment = _enrichment_figure(
                            top_enr_dir,
                            title=f'{comparison} · {database} ORA · {direction}-regulated genes',
                            database=database,
                        )
                        if fig_enrichment is not None:
                            exported, audit = _save_enrichment_figure(
                                self, fig_enrichment, plots_dir, results_dir, output_key,
                                f'{comparison} · {database} ORA 富集图（{direction}）',
                            )
                            result_files.extend(exported)
                            figure_audits.append(audit)

                if all_split_results:
                    combined = _add_enrichment_metadata(
                        pd.concat(all_split_results, ignore_index=True), method, database,
                        comparison=comparison,
                    )
                    output_key = _enrichment_output_key(
                        method, database, 'directional', comparison=comparison)
                    csv_path = os.path.join(results_dir, f'{output_key}_results.csv')
                    combined.to_csv(csv_path, index=False)
                    result_files.append({
                        'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                        'label': f'{comparison} · {database} 分方向 ORA 富集结果',
                    })
                    n_sig = len(combined[combined['P-value'] < pvalue_cutoff]) if 'P-value' in combined.columns else 0

            else:
                enr = ov.bulk.geneset_enrichment(
                    gene_list=deg_genes,
                    pathways_dict=pathways_dict,
                    pvalue_type='adjust',
                    pvalue_threshold=pvalue_cutoff,
                    organism=organism_lower,
                    outdir=os.path.join(self.project_dir, 'enrichr_tmp')
                )

                self.progress(60, "保存结果表...")
                enr = _add_enrichment_metadata(enr, method, database, comparison=comparison)
                output_key = _enrichment_output_key(method, database, comparison=comparison)
                csv_path = os.path.join(results_dir, f'{output_key}_results.csv')
                enr.to_csv(csv_path, index=False)
                result_files.append({
                    'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                    'label': f'{comparison} · {database} ORA 富集结果',
                })

                self.progress(70, "生成气泡图...")
                top_enr = enr.head(top_n)
                if len(top_enr) > 0:
                    fig_enrichment = _enrichment_figure(
                        top_enr,
                        title=f'{comparison} · {database} ORA enrichment',
                        database=database,
                    )
                    if fig_enrichment is not None:
                        exported, audit = _save_enrichment_figure(
                            self, fig_enrichment, plots_dir, results_dir, output_key,
                            f'{comparison} · {database} ORA 富集图',
                        )
                        result_files.extend(exported)
                        figure_audits.append(audit)

                n_sig = len(enr[enr['P-value'] < pvalue_cutoff]) if 'P-value' in enr.columns else len(enr)

        elif method == 'GSEA':
            if gene_rnk is None or len(gene_rnk) == 0:
                raise ValueError("GSEA 需要排序的基因列表。请确保 DEG 结果包含 log2FC 列。")

            pre_res = ov.bulk.geneset_enrichment_GSEA(
                gene_rnk=gene_rnk,
                pathways_dict=pathways_dict,
                processes=4,
                permutation_num=1000,
                outdir=os.path.join(self.project_dir, 'enrichr_gsea_tmp')
            )

            self.progress(71, "处理 GSEA 结果...")
            enr = _ensure_gsea_term_column(pre_res.res2d)
            # omicverse 将通路名保存在索引中；导出 index=False 前显式保留为 Term 列。
            enr_sig = enr[enr['fdr'] < pvalue_cutoff].copy()

            enr = _add_enrichment_metadata(enr, method, database, comparison=comparison)
            output_key = _enrichment_output_key(method, database, comparison=comparison)
            csv_path = os.path.join(results_dir, f'{output_key}_results.csv')
            enr.to_csv(csv_path, index=False)
            result_files.append({
                'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                'label': f'{comparison} · {database} GSEA 富集结果',
            })

            self.progress(75, "生成 GSEA 图表...")
            if len(enr_sig) > 0:
                top_gsea = enr_sig.head(top_n)
                fig_enrichment = _enrichment_figure(
                    top_gsea,
                    title=f'{comparison} · {database} GSEA enrichment',
                    database=database,
                    score_column='nes' if 'nes' in top_gsea.columns else None,
                )
                if fig_enrichment is not None:
                    exported, audit = _save_enrichment_figure(
                        self, fig_enrichment, plots_dir, results_dir, output_key,
                        f'{comparison} · {database} GSEA 富集图',
                    )
                    result_files.extend(exported)
                    figure_audits.append(audit)

            n_sig = len(enr_sig)
        else:
            raise ValueError(f"不支持的富集方法: {method}")

        self.progress(90, "保存输出...")

        # Results created before contrast-aware enrichment stored no source
        # metadata.  Repair registered legacy tables before building the shared
        # overview so an old GO result and a new KEGG result can be grouped.
        _backfill_legacy_enrichment_metadata(self.project_dir, results_dir)

        # Rebuild the project-level table after every run, retaining source
        # database/method/direction for valid cross-database comparison.
        result_files.extend(_write_enrichment_integration(results_dir))

        # Clean up temporary enrichment directories
        import shutil
        for tmp_suffix in ['_tmp', '_up_tmp', '_down_tmp', '_gsea_tmp']:
            tmp_dir = os.path.join(self.project_dir, f'enrichr{tmp_suffix}')
            if os.path.isdir(tmp_dir):
                shutil.rmtree(tmp_dir, ignore_errors=True)

        self.progress(100, "完成")

        return {
            'output_adata': None,
            'result_files': result_files,
            'summary': {
                'method': method,
                'database': database,
                'comparison': comparison,
                'organism': organism,
                'n_input_genes': len(deg_genes) if method == 'ORA' else (len(gene_rnk) if gene_rnk is not None else 0),
                'n_significant': n_sig,
                'pvalue_cutoff': pvalue_cutoff,
                'input_mode': 'custom_genes' if custom_genes_str else ('split_direction' if split_direction else 'standard'),
                'figure_layout_qc': {
                    'n_figures': len(figure_audits),
                    'n_figures_with_overlap': sum(
                        bool(audit.get('n_overlap_pairs', 0) > 0 or audit.get('out_of_bounds_labels'))
                        for audit in figure_audits
                    ),
                },
            }
        }
