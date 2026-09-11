import os
import json
import re
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from modules.base import BaseAnalysis
from modules.enrichment_statistics import (
    add_gsea_leading_edge_metrics,
    annotate_redundancy,
    gsea_leading_edge_table,
    geneset_manifest,
    load_human_genesets,
    load_human_genesets_text,
    map_gene_values,
    mapping_qc,
    prepare_deg_gene_table,
    prepare_gsea_ranking,
    run_ora_full,
    run_gsea_prerank,
    split_gene_text,
)


def _ensure_gsea_term_column(result_df):
    """Preserve pathway names stored in a GSEA result index."""
    result = result_df.copy()
    if 'Term' not in result.columns:
        index_name = result.index.name or 'index'
        result = result.reset_index().rename(columns={index_name: 'Term'})
    return result


ONTOLOGY_COLORS = {
    'BP': '#FDBE85',
    'CC': '#7A6FA8',
    'MF': '#B8A9D1',
    'KEGG': '#7BC77B',
    'OTHER': '#9CB8D8',
}


_GO_PRIORITY_DATABASES = ('GO_BP', 'GO_CC', 'GO_MF')
_PATHWAY_TRIPTYCH_DATABASES = ('KEGG', 'Reactome', 'WikiPathways')


def _as_enabled(value):
    """Interpret a persisted checkbox value without accepting arbitrary truthy text."""
    if isinstance(value, str):
        return value.strip().lower() in {'1', 'true', 'on', 'yes'}
    return value is True or value == 1


def _go_priority_config_from_params(params, *, databases, top_n):
    """Validate the optional Bulk GO-priority display contract.

    The triptych is a display layer over three independently tested GO
    ontologies.  It can therefore only be requested for a batch ORA run that
    contains BP, CC, and MF; it never widens the selected libraries or changes
    their individual FDR families.
    """
    params = params or {}
    if not _as_enabled(params.get('go_priority_enabled', False)):
        return None

    if str(params.get('method', 'ORA') or 'ORA').upper() != 'ORA':
        raise ValueError('主题优先 GO 三分区图仅适用于 ORA；GSEA 保留原有图型。')
    missing = [database for database in _GO_PRIORITY_DATABASES if database not in set(databases)]
    if missing:
        raise ValueError(
            '主题优先 GO 三分区图需要批量同时选择 GO_BP、GO_CC、GO_MF；缺少：'
            + '、'.join(missing)
        )

    # The single-cell workflow owns the deterministic priority policy.  Import
    # only after the Bulk feature is explicitly enabled, avoiding a module-level
    # workflow dependency for ordinary enrichment runs.
    from modules.sc_cell_go import _go_priority_allocation_from_params, _parse_focus_terms

    focus_terms = _parse_focus_terms(params.get('focus_terms'))
    focus_label = str(
        params.get('focus_label', '炎症与脂代谢优先') or '炎症与脂代谢优先'
    ).strip()
    if not focus_label:
        focus_label = '炎症与脂代谢优先'
    if len(focus_label) > 80:
        raise ValueError('focus_label 不能超过 80 个字符')
    try:
        display_top_n = int(float(params.get('go_priority_top_n', top_n)))
    except (TypeError, ValueError) as exc:
        raise ValueError('go_priority_top_n 必须为整数') from exc
    if display_top_n < 1 or display_top_n > 12:
        raise ValueError('go_priority_top_n 必须在 1–12 之间')
    allocation_mode, priority_allocation = _go_priority_allocation_from_params(
        params, top_n=display_top_n,
    )
    return {
        'focus_terms': focus_terms,
        'focus_label': focus_label,
        'allocation_mode': allocation_mode,
        'priority_allocation': priority_allocation,
        'top_n': display_top_n,
    }


def _pathway_triptych_config_for_batch(params, *, databases, top_n):
    """Return the Human pathway triptych display contract when applicable.

    The plot is produced automatically for an ORA batch that contains all
    three supported Human pathway databases.  It remains a display-only layer:
    every library keeps its own testing family and FDR correction.
    """
    if str((params or {}).get('method', 'ORA') or 'ORA').upper() != 'ORA':
        return None
    if not set(_PATHWAY_TRIPTYCH_DATABASES).issubset(set(databases)):
        return None
    return {'top_n': min(12, max(1, int(top_n)))}


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


def _with_gene_ratio(result_df, input_gene_count):
    """Add the ORA GeneRatio numerator/input-list denominator explicitly."""

    result = result_df.copy()
    denominator = max(1, int(input_gene_count or 0))
    count_column = next((column for column in ('Count', 'num', 'Gene Count')
                         if column in result.columns), None)
    if count_column is not None:
        counts = pd.to_numeric(result[count_column], errors='coerce')
    elif 'Overlap' in result.columns:
        counts = pd.to_numeric(
            result['Overlap'].astype(str).str.split('/', n=1).str[0], errors='coerce')
    else:
        return result
    result['GeneRatio'] = counts / float(denominator)
    result['InputGeneCount'] = denominator
    return result


def _enrichment_plot_suite(params, method):
    """Return additional Nature views requested by the enrichment UI.

    The existing dotplot remains the primary result for backwards
    compatibility, and is now also emitted through the Nature renderer so the
    complete suite shares the same SVG/PDF/PNG contract.  Incompatible views
    are never silently substituted with a generic seaborn/matplotlib chart.
    """
    if (params or {}).get('_batch_subrun'):
        # The primary dotplot is still emitted by the standard execution path.
        # Suppress extra variants so five databases do not create dozens of
        # near-duplicate images; the parent batch writes one overview.
        return []
    value = str((params or {}).get('enrichment_plot_suite', '完整 Nature 套图') or '').lower()
    core = value in {'core', 'minimal', '核心图', '核心'}
    if str(method).upper() == 'ORA':
        return ['enrichment_dotplot', 'enrichment_barplot'] if core else [
            'enrichment_dotplot', 'enrichment_barplot', 'enrichment_chord',
            'enrichment_cnetplot', 'enrichment_emapplot',
        ]
    return ['gsea', 'enrichment_barplot', 'gsea_running']


def _gsea_running_payload(pre_res, result_df, term_n=2, *, pathway_selection='top', target_pathways=()):
    """Build a lossless running-score payload from OmicVerse's GSEA object."""
    details = getattr(pre_res, 'results', None) or {}
    ranking = getattr(pre_res, 'ranking', None)
    if ranking is None or not details:
        return []
    terms = result_df.copy()
    term_col = 'Term' if 'Term' in terms.columns else terms.index.name or 'index'
    if term_col not in terms.columns:
        terms = terms.reset_index().rename(columns={terms.index.name or 'index': 'Term'})
        term_col = 'Term'
    score_col = next((column for column in ('nes', 'NES') if column in terms.columns), None)
    target_mode = bool(target_pathways) and pathway_selection in {'selected', 'selected_plus_top'}
    if target_mode:
        from figure_engine.templates.common import requested_term_rows, strip_term_id

        terms['_term'] = terms[term_col].map(str)
        terms['_display_term'] = terms['_term'].map(strip_term_id)
        targets, _, _ = requested_term_rows(terms, target_pathways)
        if pathway_selection == 'selected':
            terms = targets
        else:
            remainder = terms.loc[~terms.index.isin(targets.index)]
            terms = pd.concat([targets, remainder], ignore_index=False)
        terms = terms.head(max(1, int(term_n)))
    elif score_col:
        terms['_nes_numeric'] = pd.to_numeric(terms[score_col], errors='coerce')
        terms['_abs_nes'] = terms['_nes_numeric'].abs()
        # Keep positive and negative enrichment represented before limiting the
        # number of running-score panels.
        positive = terms[terms['_nes_numeric'] >= 0].sort_values('_abs_nes', ascending=False).head(max(1, int(np.ceil(term_n / 2))))
        negative = terms[terms['_nes_numeric'] < 0].sort_values('_abs_nes', ascending=False).head(max(0, int(np.floor(term_n / 2))))
        terms = pd.concat([positive, negative], ignore_index=True)
        if terms.empty:
            terms = result_df.copy()
    else:
        terms = terms.head(max(1, int(term_n)))
    payload = []
    for term in terms[term_col].head(max(1, int(term_n))):
        key = str(term)
        record = details.get(term) or details.get(key)
        if not record:
            continue
        hits = record.get('hit_indices', record.get('hits'))
        res = record.get('RES', record.get('res'))
        if hits is None and res is None:
            continue
        row = terms[terms[term_col].astype(str) == key].head(1)
        nes = row[score_col].iloc[0] if score_col and not row.empty else record.get('nes')
        fdr_col = next((column for column in ('fdr', 'FDR', 'FDR q-val') if column in terms.columns), None)
        fdr = row[fdr_col].iloc[0] if fdr_col and not row.empty else record.get('fdr')
        payload.append({
            'term': key,
            'ranking': ranking,
            'hit_indices': hits,
            'RES': res,
            'NES': nes,
            'FDR': fdr,
        })
    return payload


def _render_enrichment_variants(analysis, result_df, *, method, title, base_output_key,
                                label, params=None, pre_res=None,
                                semantic_warnings=()):
    """Render the fixed additional pathway views through NatureFigureDirector."""
    from figure_engine import NatureFigureDirector

    params = params or {}
    director = NatureFigureDirector()
    result_files = []
    audits = []
    suite = _enrichment_plot_suite(params, method)
    for plot_type in suite:
        data = result_df
        if plot_type == 'gsea_running':
            data = {'curves': _gsea_running_payload(
                pre_res, result_df, term_n=int(params.get('running_term_n', 2) or 2),
                pathway_selection=str(params.get('pathway_selection', 'top')),
                target_pathways=tuple(
                    str(params.get('target_pathways', '') or '').replace(';', ',').replace('\n', ',').split(',')
                ),
            )} if pre_res is not None else {'curves': []}
            if not data['curves']:
                # A summary CSV cannot reconstruct a running score.  Keep the
                # main GSEA dotplot and report the missing contract explicitly.
                continue
        try:
            spec = director.spec_from_params(
                plot_type, params,
                title=title,
                top_n=min(60, int(params.get('top_n', 12) or 12)),
                max_genes=int(params.get('enrichment_max_genes', 30) or 30),
                similarity_threshold=float(params.get('similarity_threshold', .15) or .15),
                running_term_n=int(params.get('running_term_n', 2) or 2),
                # Enrichment publication outputs always include the editable
                # vector pair and a raster fallback; UI preview formats cannot
                # silently remove PDF from the Nature contract.
                formats=('svg', 'pdf', 'png'),
            )
            figure = director.render(spec, data)
            if semantic_warnings:
                existing = getattr(figure, '_nature_semantic_warnings', None)
                if existing is not None:
                    existing.extend(str(item) for item in semantic_warnings if str(item).strip())
            suffix = {
                'enrichment_dotplot': 'dotplot',
                'gsea': 'gsea',
                'enrichment_barplot': 'barplot',
                'enrichment_chord': 'chord',
                'enrichment_cnetplot': 'cnetplot',
                'enrichment_emapplot': 'emapplot',
                'gsea_running': 'gsea_running',
            }[plot_type]
            exported, audit = _save_enrichment_figure(
                analysis, figure, os.path.join(analysis.project_dir, 'plots'),
                os.path.join(analysis.project_dir, 'results'),
                f'{base_output_key}_{suffix}', f'{label} · {suffix}',
            )
            result_files.extend(exported)
            audits.append(audit)
        except (ValueError, TypeError, KeyError) as exc:
            # A view that lacks its required scientific fields is skipped with a
            # machine-readable warning in the task summary; the valid views stay.
            audits.append({'status': 'warning', 'nature_readiness_score': 0,
                           'issues': [f'{plot_type}: {exc}'],
                           'n_overlap_pairs': 0, 'out_of_bounds_labels': []})
    return result_files, audits


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
    """Draw one Nature dotplot overview for one comparable group.

    Comparisons, methods and Up/Down subsets are deliberately grouped before
    this function is called.  The renderer then gives every database its own
    panel while keeping GeneRatio/NES, FDR and Count encodings consistent.
    """
    if 'Significant' in group.columns and group['Significant'].notna().any():
        significant = group['Significant']
        if significant.dtype != bool:
            significant = significant.astype(str).str.lower().isin({'true', '1', 'yes'})
        group = group.loc[significant].copy()
        if group.empty:
            return None
    databases = [str(value) for value in group['Database'].dropna().unique().tolist()]
    if len(databases) < 2:
        return None
    from figure_engine import NatureFigureDirector

    direction_text = '' if str(direction) == 'All' else f' · {direction}'
    plot_title = f'{comparison} · {method}{direction_text} pathway enrichment overview'
    director = NatureFigureDirector()
    spec = director.create_spec(
        'enrichment_overview', width='double', top_n=8,
        formats=('svg', 'pdf', 'png'), title=plot_title,
        database_scope=tuple(databases),
    )
    return director.render(spec, group)


def _write_integrated_overview_figures(integrated, results_path, *, plots_dir=None,
                                       output_prefix='enrichment_overview'):
    """Export one overview image for each comparable contrast/database group."""
    import matplotlib.pyplot as plt

    plots_dir = Path(plots_dir) if plots_dir else Path(results_path).parent / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    group_columns = ['Comparison', 'Method', 'Direction']
    for (comparison, method, direction), group in integrated.groupby(group_columns, dropna=False, sort=True):
        figure = _integrated_enrichment_overview(group, str(comparison), str(method), str(direction))
        if figure is None:
            continue
        output_key = '_'.join([
            output_prefix, _safe_output_fragment(method),
            _safe_output_fragment(comparison), _safe_output_fragment(direction),
        ])
        from figure_engine import export_registered_figure

        spec = getattr(figure, '_nature_spec_object', None)
        if spec is None:
            plt.close(figure)
            continue
        exported, report = export_registered_figure(
            figure, plots_dir / output_key, spec,
            category='enrichment_overview', label=(
                f'{comparison} · {method}'
                + ('' if str(direction) == 'All' else f' · {direction}')
                + ' 多数据库富集概览'
            ),
            formats=('svg', 'pdf', 'png'),
            qa_path=Path(results_path) / f'{output_key}_nature_readiness.json',
        )
        plt.close(figure)
        label = f'{comparison} · {method}' + ('' if str(direction) == 'All' else f' · {direction}')
        outputs.extend(exported)
    return outputs


def _write_go_priority_figures(integrated, results_path, *, plots_dir=None,
                               output_prefix='enrichment', config):
    """Write Bulk counterparts of the audited single-cell GO triptych.

    ``integrated`` retains complete result tables from every selected library.
    This function selects only rows already significant under their own full
    GO-ontology FDR correction.  It writes the exact selected rows and the
    per-ontology decision record before rendering, so the display priorities
    remain inspectable without rerunning enrichment.
    """
    from figure_engine import NatureFigureDirector, export_registered_figure
    from modules.sc_cell_go import _GO_PRIORITY_POLICY, _select_go_priority_rows
    import matplotlib.pyplot as plt

    plots_path = Path(plots_dir) if plots_dir else Path(results_path).parent / 'plots'
    plots_path.mkdir(parents=True, exist_ok=True)
    results_path = Path(results_path)
    reverse_database_names = {
        'GO_BP': 'GO_Biological_Process_2023',
        'GO_CC': 'GO_Cellular_Component_2023',
        'GO_MF': 'GO_Molecular_Function_2023',
    }
    source_columns = [
        'Comparison', 'Method', 'Direction', 'Database', 'Term', 'Overlap',
        'Adjusted P-value', 'Genes', 'selection_rank', 'selection_reason',
        'selection_priority', 'selection_policy',
    ]
    go_rows = integrated.loc[
        integrated.get('Database', pd.Series('', index=integrated.index)).isin(
            _GO_PRIORITY_DATABASES
        )
    ].copy()
    outputs = []
    selected_frames = []
    audit = {
        'generated_at': datetime.now(timezone.utc).isoformat(),
        'focus_label': config['focus_label'],
        'requested_focus_terms': list(config['focus_terms']),
        'statistics_note': (
            'ORA 在完整 GO 本体内分别检验，Adjusted P-value 为各完整本体内的 FDR；'
            '本审计仅记录展示排序，不重新校正。'
        ),
        'selection_policy': _GO_PRIORITY_POLICY,
        'allocation_mode': config['allocation_mode'],
        'requested_allocation': config['priority_allocation'],
        'display_cap_per_ontology': int(config['top_n']),
        'units': [],
        'warnings': [],
    }
    if go_rows.empty:
        audit['warnings'].append('整合结果中没有 GO_BP、GO_CC 或 GO_MF 行，未生成三分区图。')
    else:
        group_columns = ['Comparison', 'Method', 'Direction']
        director = NatureFigureDirector()
        for group_key, source_group in go_rows.groupby(
            group_columns, dropna=False, sort=True, observed=True,
        ):
            comparison, method, direction = (str(value) for value in group_key)
            present = set(source_group['Database'].dropna().astype(str))
            missing = [database for database in _GO_PRIORITY_DATABASES if database not in present]
            unit_audit = {
                'comparison': comparison,
                'method': method,
                'direction': direction,
                'available_databases': sorted(present),
            }
            if method.upper() != 'ORA':
                unit_audit.update({'status': 'skipped', 'reason': 'GO 三分区图只适用于 ORA'})
                audit['units'].append(unit_audit)
                continue
            if missing:
                message = (
                    f'{comparison} / {direction}: 缺少 ' + '、'.join(missing)
                    + ' 的完整结果，未生成 GO 三分区图。'
                )
                unit_audit.update({'status': 'skipped', 'reason': message})
                audit['warnings'].append(message)
                audit['units'].append(unit_audit)
                continue

            # Reuse the single-cell selector only after adapting the common
            # integrated-table metadata to its explicit completed-run contract.
            # No test family, p value, or result row is altered here.
            selection_input = source_group.copy()
            selection_input['status'] = 'completed'
            selection_input['gene_set'] = selection_input['Database'].map(reverse_database_names)
            selection_input['method'] = selection_input['Method'].astype(str)
            selection_input['direction'] = selection_input['Direction'].astype(str)
            selected, selection_audit = _select_go_priority_rows(
                selection_input, config['focus_terms'], top_n=config['top_n'],
                priority_allocation=config['priority_allocation'],
            )
            selected = selected.drop(columns=['status', 'gene_set', 'method', 'direction'], errors='ignore')
            unit_audit.update(selection_audit)
            unit_audit['status'] = 'completed'
            audit['units'].append(unit_audit)
            if not selected.empty:
                selected_frames.append(selected)

            direction_label = {
                'Up': '上调', 'Down': '下调', 'All': '全部显著 DEG',
            }.get(direction, direction)
            safe_key = '_'.join((
                output_prefix, 'go_priority', _safe_output_fragment(comparison),
                _safe_output_fragment(direction, fallback='all'),
            ))
            title_suffix = f'{comparison} · {direction_label} · {config["focus_label"]}'
            height_mm = min(225, max(132, 62 + 6 * max(len(selected), 6)))
            for plot_type, suffix, heading in (
                ('go_priority_dotplot', 'dotplot', 'GO Enrichment Plot'),
                ('go_priority_barplot', 'barplot', 'GO enrichment'),
            ):
                figure = None
                try:
                    spec = director.create_spec(
                        plot_type, width='double', height_mm=height_mm,
                        top_n=config['top_n'], formats=('svg', 'pdf', 'png'),
                        database_scope=_GO_PRIORITY_DATABASES,
                        title=f'{heading}\n{title_suffix}',
                    )
                    figure = director.render(spec, selected)
                    exported, report = export_registered_figure(
                        figure, plots_path / f'{safe_key}_{suffix}', spec,
                        category=plot_type,
                        label=(
                            f'{config["focus_label"]} GO 主题优先 · {comparison} · '
                            f'{direction_label} {"气泡图" if suffix == "dotplot" else "柱状图"}'
                            '（炎症→脂代谢→Top；全库 FDR）'
                        ),
                        qa_path=results_path / f'{safe_key}_{suffix}_nature_readiness.json',
                    )
                    outputs.extend(exported)
                    unit_audit.setdefault('figures', []).append({
                        'plot_type': plot_type,
                        'status': 'pass' if report.ready else 'warning',
                        'nature_readiness_score': report.score,
                        'issues': [issue.message for issue in report.issues],
                    })
                except Exception as exc:
                    message = f'{comparison} / {direction}: GO 主题优先 {suffix} 生成失败（{exc}）'
                    audit['warnings'].append(message)
                    unit_audit.setdefault('figures', []).append({
                        'plot_type': plot_type, 'status': 'failed', 'error': str(exc),
                    })
                finally:
                    if figure is not None:
                        plt.close(figure)

    selected_source = (
        pd.concat(selected_frames, ignore_index=True, sort=False)
        if selected_frames else pd.DataFrame(columns=source_columns)
    )
    available_columns = [column for column in source_columns if column in selected_source.columns]
    selected_source = selected_source[available_columns]
    source_path = results_path / f'{output_prefix}_go_priority_source.csv'
    selected_source.to_csv(source_path, index=False)
    outputs.append({
        'file_path': str(source_path), 'file_type': 'csv', 'category': 'table',
        'label': 'GO 主题优先三分区图来源表（炎症→脂代谢→Top；全库 FDR）',
    })
    audit_path = results_path / f'{output_prefix}_go_priority_selection_audit.json'
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    outputs.append({
        'file_path': str(audit_path), 'file_type': 'json', 'category': 'info',
        'label': 'GO 主题优先三分区图选择审计（可复现）',
    })
    return outputs


def _write_pathway_triptych_figures(integrated, results_path, *, plots_dir=None,
                                    output_prefix='enrichment', config):
    """Write KEGG/Reactome/WikiPathways Human FDR-top triptych figures."""
    from figure_engine import NatureFigureDirector, export_registered_figure
    from modules.sc_cell_go import (
        PATHWAY_TRIPTYCH_DATABASES,
        _select_pathway_triptych_rows,
    )
    import matplotlib.pyplot as plt

    results_path = Path(results_path)
    plots_path = Path(plots_dir) if plots_dir else results_path.parent / 'plots'
    plots_path.mkdir(parents=True, exist_ok=True)
    reverse_database_names = {
        'KEGG': 'KEGG_2021_Human',
        'Reactome': 'Reactome_2022',
        # Bulk uses a version-pinned 2019 cache, whereas the single-cell
        # registry uses the 2021 snapshot name.  This is an in-memory adapter
        # for the common display selector only; result provenance remains in
        # the original integrated table and run manifest.
        'WikiPathways': 'WikiPathway_2021_Human',
    }
    source_columns = [
        'Comparison', 'Method', 'Direction', 'Database', 'Term', 'Overlap',
        'Adjusted P-value', 'Genes', 'selection_rank', 'selection_reason',
        'selection_priority', 'selection_policy',
    ]
    pathway_rows = integrated.loc[
        integrated.get('Database', pd.Series('', index=integrated.index)).isin(
            _PATHWAY_TRIPTYCH_DATABASES
        )
    ].copy()
    outputs = []
    selected_frames = []
    audit = {
        'generated_at': datetime.now(timezone.utc).isoformat(),
        'statistics_note': (
            'KEGG、Reactome 与 WikiPathways Human 在完整库内分别检验，'
            'Adjusted P-value 为各库内独立 FDR；本审计仅记录 FDR Top 展示排序。'
        ),
        'selection_policy': '每个数据库仅展示已显著 term 的 FDR Top 通路。',
        'display_cap_per_database': int(config['top_n']),
        'units': [],
        'warnings': [],
    }
    if pathway_rows.empty:
        audit['warnings'].append('整合结果中没有 KEGG、Reactome 或 WikiPathways 行，未生成 Human 通路三分区图。')
    else:
        director = NatureFigureDirector()
        group_columns = ['Comparison', 'Method', 'Direction']
        for group_key, source_group in pathway_rows.groupby(
            group_columns, dropna=False, sort=True, observed=True,
        ):
            comparison, method, direction = (str(value) for value in group_key)
            present = set(source_group['Database'].dropna().astype(str))
            missing = [database for database in _PATHWAY_TRIPTYCH_DATABASES if database not in present]
            unit_audit = {
                'comparison': comparison, 'method': method, 'direction': direction,
                'available_databases': sorted(present),
            }
            if method.upper() != 'ORA':
                unit_audit.update({'status': 'skipped', 'reason': 'Human 通路三分区图只适用于 ORA'})
                audit['units'].append(unit_audit)
                continue
            if missing:
                message = (
                    f'{comparison} / {direction}: 缺少 ' + '、'.join(missing)
                    + ' 的完整结果，未生成 Human 通路三分区图。'
                )
                unit_audit.update({'status': 'skipped', 'reason': message})
                audit['warnings'].append(message)
                audit['units'].append(unit_audit)
                continue

            selection_input = source_group.copy()
            selection_input['status'] = 'completed'
            selection_input['gene_set'] = selection_input['Database'].map(reverse_database_names)
            selection_input['method'] = selection_input['Method'].astype(str)
            selection_input['direction'] = selection_input['Direction'].astype(str)
            selected, selection_audit = _select_pathway_triptych_rows(
                selection_input, top_n=config['top_n'],
            )
            selected = selected.drop(columns=['status', 'gene_set', 'method', 'direction'], errors='ignore')
            unit_audit.update(selection_audit)
            unit_audit['status'] = 'completed'
            audit['units'].append(unit_audit)
            if not selected.empty:
                selected_frames.append(selected)

            direction_label = {
                'Up': '上调', 'Down': '下调', 'All': '全部显著 DEG',
            }.get(direction, direction)
            safe_key = '_'.join((
                output_prefix, 'pathway_triptych', _safe_output_fragment(comparison),
                _safe_output_fragment(direction, fallback='all'),
            ))
            height_mm = min(225, max(132, 62 + 6 * max(len(selected), 6)))
            for plot_type, suffix, heading in (
                ('pathway_triptych_dotplot', 'dotplot', 'Human Pathway Enrichment Plot'),
                ('pathway_triptych_barplot', 'barplot', 'Human pathway enrichment'),
            ):
                figure = None
                try:
                    spec = director.create_spec(
                        plot_type, width='double', height_mm=height_mm,
                        top_n=config['top_n'], formats=('svg', 'pdf', 'png'),
                        database_scope=tuple(PATHWAY_TRIPTYCH_DATABASES.values()),
                        title=f'{heading}\n{comparison} · {direction_label}',
                    )
                    figure = director.render(spec, selected)
                    exported, report = export_registered_figure(
                        figure, plots_path / f'{safe_key}_{suffix}', spec,
                        category=plot_type,
                        label=(
                            f'KEGG/Reactome/WikiPathways Human 三分区 · {comparison} · '
                            f'{direction_label} {"气泡图" if suffix == "dotplot" else "柱状图"}'
                            '（各库独立 FDR Top）'
                        ),
                        qa_path=results_path / f'{safe_key}_{suffix}_nature_readiness.json',
                    )
                    outputs.extend(exported)
                    unit_audit.setdefault('figures', []).append({
                        'plot_type': plot_type,
                        'status': 'pass' if report.ready else 'warning',
                        'nature_readiness_score': report.score,
                        'issues': [issue.message for issue in report.issues],
                    })
                except Exception as exc:
                    message = f'{comparison} / {direction}: Human 通路三分区 {suffix} 生成失败（{exc}）'
                    audit['warnings'].append(message)
                    unit_audit.setdefault('figures', []).append({
                        'plot_type': plot_type, 'status': 'failed', 'error': str(exc),
                    })
                finally:
                    if figure is not None:
                        plt.close(figure)

    selected_source = (
        pd.concat(selected_frames, ignore_index=True, sort=False)
        if selected_frames else pd.DataFrame(columns=source_columns)
    )
    selected_source = selected_source[[
        column for column in source_columns if column in selected_source.columns
    ]]
    source_path = results_path / f'{output_prefix}_pathway_triptych_source.csv'
    selected_source.to_csv(source_path, index=False)
    outputs.append({
        'file_path': str(source_path), 'file_type': 'csv', 'category': 'table',
        'label': 'KEGG/Reactome/WikiPathways Human 三分区图来源表（各库独立 FDR Top）',
    })
    audit_path = results_path / f'{output_prefix}_pathway_triptych_selection_audit.json'
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    outputs.append({
        'file_path': str(audit_path), 'file_type': 'json', 'category': 'info',
        'label': 'KEGG/Reactome/WikiPathways Human 三分区图选择审计（可复现）',
    })
    return outputs


def _write_enrichment_integration(results_dir, *, source_paths=None, plots_dir=None,
                                  output_prefix='enrichment', go_priority_config=None,
                                  pathway_triptych_config=None):
    """Create an integrated table from explicitly selected enrichment outputs.

    Each source table carries ``Comparison``, ``Database``, ``Method`` and
    ``Direction``.  The integration intentionally retains those fields instead
    of merging similarly named pathways across libraries, whose gene-set
    definitions may differ. ``source_paths`` is required for task-scoped batch
    execution; the legacy directory scan remains only for older single runs.
    """
    results_path = Path(results_dir)
    if source_paths is None:
        source_files = sorted(
            path for path in results_path.glob('enrichment_*_results.csv')
            if path.name not in {'enrichment_integrated_results.csv'}
        )
    else:
        source_files = [Path(path) for path in source_paths]
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
        'Enrichment P-value', 'Significant', 'NES', 'GeneRatio', 'BackgroundRatio',
        'Count', 'Overlap', 'Genes', 'Source result',
    ]
    ordered_columns = [column for column in preferred if column in integrated.columns]
    ordered_columns.extend(column for column in integrated.columns if column not in ordered_columns)
    integrated = integrated[ordered_columns].sort_values(
        ['Enrichment FDR', 'Comparison', 'Database', 'Direction', 'Term'],
        na_position='last', kind='stable',
    )

    output_files = []
    csv_path = results_path / f'{output_prefix}_integrated_results.csv'
    integrated.to_csv(csv_path, index=False)
    output_files.append({
        'file_path': str(csv_path), 'file_type': 'csv', 'category': 'table',
        'label': f'通路富集整合表（{len(integrated)} 条）',
    })
    output_files.extend(_write_integrated_overview_figures(
        integrated, results_path, plots_dir=plots_dir,
        output_prefix=f'{output_prefix}_overview',
    ))
    if go_priority_config is not None:
        output_files.extend(_write_go_priority_figures(
            integrated, results_path, plots_dir=plots_dir,
            output_prefix=output_prefix, config=go_priority_config,
        ))
    if pathway_triptych_config is not None:
        output_files.extend(_write_pathway_triptych_figures(
            integrated, results_path, plots_dir=plots_dir,
            output_prefix=output_prefix, config=pathway_triptych_config,
        ))

    try:
        xlsx_path = results_path / f'{output_prefix}_integrated_results.xlsx'
        summary_source = integrated.copy()
        if 'Significant' in summary_source.columns:
            summary_source['_significant'] = (
                summary_source['Significant'].astype(str).str.lower().isin({'true', '1', 'yes'})
            )
            summary = (summary_source.groupby(
                ['Comparison', 'Database', 'Method', 'Direction'], dropna=False,
            ).agg(n_tested=('Term', 'size'), n_significant=('_significant', 'sum')).reset_index())
        else:
            summary = (summary_source.groupby(
                ['Comparison', 'Database', 'Method', 'Direction'], dropna=False,
            ).size().reset_index(name='n_tested'))
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
        'file_path': str(report_path), 'file_type': 'json', 'category': 'qc',
        'label': f'富集图文字重叠检测（{audit["status"]}）',
    }


def _save_enrichment_figure(analysis, fig, plots_dir, results_dir, output_key, label):
    """Export one enrichment plot and its post-layout overlap report."""
    nature_spec = getattr(fig, '_nature_spec_object', None)
    if nature_spec is not None:
        import matplotlib.pyplot as plt
        from figure_engine import export_registered_figure

        exported, report = export_registered_figure(
            fig, os.path.join(plots_dir, output_key), nature_spec,
            category='enrichment', label=label,
            qa_path=os.path.join(results_dir, f'{output_key}_nature_readiness.json'),
        )
        audit = {
            'status': 'pass' if report.ready else 'warning',
            'nature_readiness_score': report.score,
            'n_text_elements': report.metrics.get('text_element_count', 0),
            'n_overlap_pairs': report.metrics.get('label_overlap_count', 0),
            'overlap_pairs': [],
            'out_of_bounds_labels': [
                issue.details.get('labels', [])
                for issue in report.issues if issue.code == 'axis_clipping'
            ],
            'issues': [issue.message for issue in report.issues],
            'layout_contract': getattr(fig, '_nature_expected_size_mm', ()),
        }
        plt.close(fig)
        return exported, audit

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
    project_root = Path(__file__).resolve().parent.parent
    for name in candidate_names:
        candidate_paths = [
            Path('genesets') / f'{name}.txt',
            Path('genesets') / f'{name}_{organism}.txt',
            Path('genesets') / f'{name}_{organism_lower}.txt',
            Path('data') / 'go_gene_sets' / f'{name}.gmt',
            project_root / 'data' / 'go_gene_sets' / f'{name}.gmt',
        ]
        for path in candidate_paths:
            if path.exists():
                return str(path)

    destination = Path('genesets') / f'{db_filename}.txt'
    _download_enrichr_geneset(db_filename, destination)
    return str(destination)


def _bounded_int(value, name, minimum, maximum):
    try:
        parsed = int(float(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} 必须是整数。') from exc
    if parsed < minimum or parsed > maximum:
        raise ValueError(f'{name} 必须在 {minimum}–{maximum} 之间。')
    return parsed


def _bounded_float(value, name, minimum, maximum, *, minimum_inclusive=True):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} 必须是数字。') from exc
    lower_ok = parsed >= minimum if minimum_inclusive else parsed > minimum
    if not np.isfinite(parsed) or not lower_ok or parsed > maximum:
        left = '[' if minimum_inclusive else '('
        raise ValueError(f'{name} 必须在 {left}{minimum}, {maximum}] 范围内。')
    return parsed


def _background_from_analysis_input(input_path):
    """Recover Human symbols from the selected expression input for custom ORA."""
    if not input_path or not os.path.isfile(input_path):
        raise ValueError('没有可用于构建 ORA 背景的表达输入。')
    from modules.io_utils import read_expression_matrix

    adata = read_expression_matrix(input_path)
    try:
        table = adata.var.copy().reset_index(drop=True)
        if 'gene' in table.columns:
            replacement = 'gene_symbol' if 'gene_symbol' not in table.columns else 'source_gene'
            table = table.rename(columns={'gene': replacement})
        table.insert(0, 'gene', [str(value) for value in adata.var_names])
        prepared, identifier_map, provenance = prepare_deg_gene_table(table)
        genes = [gene for gene in prepared['_analysis_gene'].tolist() if gene]
        provenance['source'] = 'analysis_input'
        return genes, identifier_map, provenance
    finally:
        try:
            if getattr(adata, 'isbacked', False):
                adata.file.close()
        except (AttributeError, OSError):
            pass


def _plot_enrichment_rows(result, params):
    """Keep the full statistical table while limiting default plots to FDR hits."""
    target_requested = (
        bool(str((params or {}).get('target_pathways', '') or '').strip())
        and str((params or {}).get('pathway_selection', 'top')) in {'selected', 'selected_plus_top'}
    )
    if target_requested:
        return result.copy()
    if 'Significant' in result.columns:
        significant = result['Significant']
        if significant.dtype != bool:
            significant = significant.astype(str).str.lower().isin({'true', '1', 'yes'})
        return result.loc[significant].copy()
    return result.copy()


def _mapping_warnings(qc):
    warnings = []
    for label, section, threshold in (
        ('输入基因', qc.get('query', {}), 0.3),
        ('背景基因', qc.get('background', {}), 0.1),
        ('GSEA 排名', qc.get('ranking', {}), 0.3),
    ):
        rate = section.get('mapping_rate')
        if rate is not None and rate < threshold:
            warnings.append(f'{label}数据库覆盖率仅 {rate:.1%}，请检查 Human gene symbol/Ensembl 映射。')
        unresolved = int(section.get('n_unresolved_ensembl', 0) or 0)
        if unresolved:
            warnings.append(f'{label}仍有 {unresolved} 个 Ensembl ID 未映射到 symbol。')
        outside = int(section.get('n_outside_background', 0) or 0)
        if outside:
            warnings.append(f'{label}有 {outside} 个基因不在 ORA 背景中，已从检验列表排除。')
    return warnings


def _write_mapping_qc(results_dir, output_key, payload):
    path = Path(results_dir) / f'{output_key}_mapping_qc.json'
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')
    return {
        'file_path': str(path), 'file_type': 'json', 'category': 'qc',
        'label': '富集分析基因映射、背景与数据库审计',
    }


def _redundancy_audit_table(result_df):
    """Return the non-destructive significant-pathway redundancy assignment."""
    columns = [
        'Comparison', 'Database', 'Method', 'Direction', 'Term',
        'Adjusted P-value', 'fdr', 'Significant', 'Genes', 'lead_genes',
        'RedundancyCluster', 'RepresentativeTerm', 'IsRepresentative',
        'JaccardToRepresentative', 'RedundancyClusterSize',
    ]
    if result_df is None or result_df.empty or 'RedundancyCluster' not in result_df.columns:
        return pd.DataFrame(columns=columns)
    assigned = result_df['RedundancyCluster'].astype(str).str.strip() != ''
    available = [column for column in columns if column in result_df.columns]
    return result_df.loc[assigned, available].copy().reset_index(drop=True)


def _qc_workbook_frames(qc_payload):
    summaries = []
    unmapped = []
    sections = []
    if qc_payload.get('background'):
        sections.append(('Background', qc_payload['background']))
    if qc_payload.get('ranking'):
        sections.append(('GSEA ranking', qc_payload['ranking']))
    for analysis in qc_payload.get('analyses', []):
        sections.append((f"Query: {analysis.get('direction', 'All')}", analysis.get('query', {})))
    for label, section in sections:
        row = {'Section': label}
        for key, value in section.items():
            if not isinstance(value, (list, dict)):
                row[key] = value
        summaries.append(row)
        for gene in section.get('unmapped_genes', []) or []:
            unmapped.append({'Section': label, 'Gene': gene, 'Reason': 'not_in_selected_library'})
        for gene in section.get('outside_background_genes', []) or []:
            unmapped.append({'Section': label, 'Gene': gene, 'Reason': 'outside_effective_background'})
    return pd.DataFrame(summaries), pd.DataFrame(
        unmapped, columns=['Section', 'Gene', 'Reason'],
    )


def _write_enrichment_methods(results_dir, output_key, qc_payload, result_df):
    """Write a concise human-readable statistical methods and provenance record."""
    method = qc_payload['method']
    library = qc_payload['gene_set_library']
    lines = [
        'Human Bulk RNA-seq pathway enrichment — reproducibility record',
        '',
        f"Comparison: {qc_payload.get('comparison')}",
        f"Method: {method}",
        f"Database: {qc_payload.get('database')}",
        f"Gene-set snapshot: {library.get('file_name')}",
        f"Gene-set SHA256: {library.get('sha256')}",
        f"Gene-set pathways / unique genes: {library.get('n_pathways')} / {library.get('n_unique_genes')}",
        'Identifier policy: Human identifiers are upper-cased; Ensembl version suffixes are removed; '
        'ID-to-symbol conversion uses only symbol/alias columns supplied with the input.',
        '',
    ]
    if method == 'ORA':
        family = (qc_payload.get('analyses') or [{}])[0].get('testing_family', {})
        background = qc_payload.get('background', {})
        lines.extend([
            'ORA test: one-sided hypergeometric over-representation test.',
            f"Universe: {background.get('n_effective_background_genes', 0)} Human genes in both the submitted tested-gene background and selected library.",
            f"Gene-set size filter: {family.get('min_size')}–{family.get('max_size')} genes after universe intersection.",
            f"Multiple testing: Benjamini–Hochberg across all {family.get('n_pathways_tested', len(result_df))} eligible pathways, including zero-hit pathways.",
            'Effect size: Haldane–Anscombe corrected odds ratio, approximate log-OR 95% CI, and fold enrichment.',
        ])
    else:
        parameters = qc_payload.get('gsea_parameters', {})
        lines.extend([
            'GSEA test: deterministic Human pre-ranked gene-set enrichment.',
            f"Ranking metric: {parameters.get('ranking_metric')}",
            f"Permutations / seed / weight: {parameters.get('permutation_num')} / {parameters.get('seed')} / {parameters.get('weight')}",
            f"Gene-set size filter: {parameters.get('min_size')}–{parameters.get('max_size')} genes after ranking intersection.",
            'Leading-edge genes are retained per pathway with rank position and ranking score.',
        ])
    lines.extend([
        f"Significance threshold: adjusted P/FDR < {qc_payload.get('pvalue_cutoff')}",
        f"Tested / significant pathways: {len(result_df)} / {int(result_df.get('Significant', pd.Series(dtype=bool)).sum())}",
        f"Redundancy annotation: hit/leading-edge gene Jaccard ≥ {qc_payload.get('redundancy_threshold')} groups significant terms; no statistical rows are removed.",
        '',
        'Warnings:',
    ])
    warnings = qc_payload.get('warnings') or ['None']
    lines.extend(f'- {warning}' for warning in warnings)
    path = Path(results_dir) / f'{output_key}_methods.txt'
    path.write_text('\n'.join(str(line) for line in lines) + '\n', encoding='utf-8')
    return {
        'file_path': str(path), 'file_type': 'txt', 'category': 'info',
        'label': '富集分析统计方法与复现说明',
    }


def _write_enrichment_source_data(results_dir, output_key, result_df, qc_payload,
                                  redundancy_df=None, leading_edge_df=None):
    """Create a per-run Source Data workbook without replacing full CSV outputs."""
    path = Path(results_dir) / f'{output_key}_source_data.xlsx'
    mapping_summary, unmapped = _qc_workbook_frames(qc_payload)
    parameters = []
    for section_name in ('gene_set_library', 'gsea_parameters'):
        for key, value in (qc_payload.get(section_name) or {}).items():
            parameters.append({
                'Section': section_name,
                'Parameter': key,
                'Value': json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict)) else value,
            })
    parameters.extend({
        'Section': 'analysis', 'Parameter': key, 'Value': value,
    } for key, value in (
        ('method', qc_payload.get('method')),
        ('database', qc_payload.get('database')),
        ('comparison', qc_payload.get('comparison')),
        ('pvalue_cutoff', qc_payload.get('pvalue_cutoff')),
        ('redundancy_threshold', qc_payload.get('redundancy_threshold')),
    ))
    try:
        with pd.ExcelWriter(path, engine='openpyxl') as writer:
            result_df.to_excel(writer, sheet_name='Full_results', index=False)
            significant_mask = result_df.get(
                'Significant', pd.Series(False, index=result_df.index),
            )
            if significant_mask.dtype != bool:
                significant_mask = significant_mask.astype(str).str.lower().isin(
                    {'true', '1', 'yes'},
                )
            significant = result_df.loc[significant_mask]
            significant.to_excel(writer, sheet_name='Significant', index=False)
            mapping_summary.to_excel(writer, sheet_name='Mapping_QC', index=False)
            unmapped.to_excel(writer, sheet_name='Unmapped_genes', index=False)
            pd.DataFrame(parameters).to_excel(writer, sheet_name='Parameters', index=False)
            if redundancy_df is not None:
                redundancy_df.to_excel(writer, sheet_name='Redundancy_clusters', index=False)
            if leading_edge_df is not None:
                if len(leading_edge_df) <= 1_048_000:
                    leading_edge_df.to_excel(writer, sheet_name='GSEA_leading_edge', index=False)
                else:
                    pd.DataFrame([{
                        'Note': 'Leading-edge table exceeds the Excel row limit; use the registered CSV.',
                        'Rows': len(leading_edge_df),
                    }]).to_excel(writer, sheet_name='GSEA_leading_edge', index=False)
    except (ImportError, ModuleNotFoundError):
        return None
    return {
        'file_path': str(path), 'file_type': 'xlsx', 'category': 'table',
        'label': '富集分析 Source Data 工作簿',
    }


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
    if 'CC' in text or 'CELLULAR COMPONENT' in text:
        return 'CC'
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
                       ontology_colors=None, params=None):
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
    # Phase 1 publication path.  ORA uses GeneRatio/Count/FDR and GSEA uses
    # NES/size/FDR.  Legacy bar rendering below is retained only for old tables
    # that do not carry enough fields for those scientific encodings.
    target_requested = bool(str((params or {}).get('target_pathways', '') or '').strip()) and str((params or {}).get('pathway_selection', 'top')) in {'selected', 'selected_plus_top'}
    if ontology_colors is None or target_requested:
        try:
            from figure_engine import NatureFigureDirector
            from figure_engine.templates.common import compress_redundant_terms

            director = NatureFigureDirector()
            plot_type = 'gsea' if score_column and score_column in df.columns else 'enrichment'
            # Estimate the physical canvas from the number of terms that will
            # actually remain after redundancy compression.  Legacy figures
            # sized from the pre-compression table produced very tall canvases
            # with only one or two visible dots and looked unfinished at review
            # size.
            effective_terms = len(df)
            genes_column = next(
                (name for name in ('Genes', 'genes', 'geneID', 'matched_genes', 'lead_genes')
                 if name in df.columns), None,
            )
            if genes_column:
                try:
                    effective_terms = len(compress_redundant_terms(
                        df, genes_column,
                        threshold=float((params or {}).get('redundancy_threshold', 0.85)),
                        maximum=min(60, len(df)),
                    )[0])
                except Exception:
                    effective_terms = len(df)
            estimated_height = min(190.0, max(100.0, 48.0 + 12.0 * min(12, effective_terms)))
            spec = director.spec_from_params(
                plot_type, params or {},
                width='double',
                title=title,
                top_n=min(60, int((params or {}).get('top_n', len(df)) or len(df))),
                height_mm=estimated_height,
                formats=('svg', 'pdf', 'png'),
            )
            return director.render(spec, df)
        except ValueError:
            if target_requested:
                raise
            # Compatibility for historical enrichment tables with P-value only or
            # no overlap/count field.  New analysis outputs use the strict path.
            pass
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
    DESCRIPTION = "GO/KEGG/Reactome/WikiPathways Human 通路富集分析（ORA / GSEA）"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    @staticmethod
    def _batch_database_names(value):
        """Normalise the explicitly selected standard database names."""
        if isinstance(value, (list, tuple, set)):
            raw_values = value
        else:
            raw_values = re.split(r'[,;\n]+', str(value or ''))
        allowed = ('GO_BP', 'GO_MF', 'GO_CC', 'KEGG', 'Reactome', 'WikiPathways')
        selected = []
        for raw in raw_values:
            database = str(raw).strip()
            if not database:
                continue
            if database not in allowed:
                raise ValueError(f'不支持的批量 Human 基因集数据库: {database}')
            if database not in selected:
                selected.append(database)
        if not selected:
            raise ValueError('请至少选择一个批量 Human 基因集数据库。')
        return selected

    def _output_directories(self):
        """Return regular paths or a validated, task-local batch workspace."""
        batch_key = str(self.params.get('_batch_output_key', '') or '').strip()
        if not batch_key:
            return (
                os.path.join(self.project_dir, 'plots'),
                os.path.join(self.project_dir, 'results'),
            )
        safe_key = _safe_output_fragment(batch_key, fallback='batch')
        project_root = Path(self.project_dir).resolve()
        results_dir = (project_root / 'results' / 'enrichment_batches' / safe_key).resolve()
        plots_dir = (project_root / 'plots' / 'enrichment_batches' / safe_key).resolve()
        if project_root not in results_dir.parents or project_root not in plots_dir.parents:
            raise ValueError('批量富集输出目录必须位于当前项目内。')
        return str(plots_dir), str(results_dir)

    def _run_database_batch(self, input_path, databases):
        """Run one fixed DEG contrast against several independent libraries.

        A library is a statistical unit here: failures are retained in the
        summary without invalidating successful libraries, and the final
        integration receives only the result paths generated in this task.
        """
        go_priority_config = _go_priority_config_from_params(
            self.params, databases=databases,
            top_n=_bounded_int(self.params.get('top_n', 12), '展示通路数', 1, 60),
        )
        pathway_triptych_config = _pathway_triptych_config_for_batch(
            self.params, databases=databases,
            top_n=_bounded_int(self.params.get('top_n', 12), '展示通路数', 1, 60),
        )
        batch_key = _safe_output_fragment(
            self.params.get('_analysis_id', 'adhoc_batch'), fallback='adhoc_batch',
        )
        batch_params = dict(self.params)
        batch_params['_batch_output_key'] = batch_key
        workspace = BulkEnrichmentAnalysis(
            project_dir=self.project_dir, params=batch_params,
            progress_callback=self._progress,
        )
        plots_dir, results_dir = workspace._output_directories()
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        result_files = []
        database_rows = []
        source_paths = []
        total = len(databases)
        for index, database in enumerate(databases, start=1):
            self.progress(
                int(3 + (index - 1) * 88 / total),
                f'批量富集 {index}/{total}: {database}',
            )
            sub_params = dict(self.params)
            sub_params.update({
                'batch_databases': False,
                'databases': '',
                'database': database,
                '_batch_subrun': True,
                '_batch_output_key': batch_key,
            })

            def sub_progress(percent, message, *, _index=index, _database=database):
                completed = (_index - 1) + max(0, min(100, int(percent))) / 100.0
                overall = int(3 + completed * 88 / total)
                self.progress(overall, f'批量富集 {_index}/{total} · {_database}: {message}')

            try:
                sub_result = BulkEnrichmentAnalysis(
                    project_dir=self.project_dir, params=sub_params,
                    progress_callback=sub_progress,
                ).run(input_path)
            except Exception as exc:
                database_rows.append({
                    'database': database, 'status': 'failed', 'error': str(exc)[:600],
                })
                continue

            result_files.extend(sub_result.get('result_files') or [])
            summary = dict(sub_result.get('summary') or {})
            full_result_path = summary.pop('full_results_path', '')
            if full_result_path and os.path.isfile(full_result_path):
                source_paths.append(full_result_path)
            database_rows.append({
                'database': database,
                'status': 'completed',
                'n_input_genes': summary.get('n_input_genes'),
                'n_background_genes': summary.get('n_background_genes'),
                'n_ranked_genes': summary.get('n_ranked_genes'),
                'n_tested': summary.get('n_tested'),
                'n_significant': summary.get('n_significant'),
                'gene_set_sha256': summary.get('gene_set_sha256'),
                'gene_set_library': summary.get('gene_set_library'),
                'mapping_rates': summary.get('mapping_rates'),
                'mapping_warnings': summary.get('mapping_warnings', []),
            })

        summary_frame = pd.DataFrame(database_rows)
        summary_path = Path(results_dir) / f'enrichment_batch_{batch_key}_summary.csv'
        summary_frame.to_csv(summary_path, index=False)
        result_files.append({
            'file_path': str(summary_path), 'file_type': 'csv', 'category': 'table',
            'label': f'批量富集数据库汇总（{len(databases)} 个）',
        })
        manifest = {
            'comparison': str(self.params.get('input_comparison', '') or ''),
            'method': str(self.params.get('method', 'ORA') or 'ORA').upper(),
            'databases_requested': list(databases),
            'databases': database_rows,
            'fdr_scope': '每个比较 × 数据库 × 方向独立 BH/FDR；概览不重新校正。',
            'source_result_files': [Path(path).name for path in source_paths],
            'go_priority_display': ({
                'enabled': True,
                'focus_terms': list(go_priority_config['focus_terms']),
                'focus_label': go_priority_config['focus_label'],
                'allocation_mode': go_priority_config['allocation_mode'],
                'requested_allocation': go_priority_config['priority_allocation'],
                'display_cap_per_ontology': go_priority_config['top_n'],
                'statistics_note': '仅调整展示排序；各 GO 本体仍使用完整本体内独立 FDR。',
            } if go_priority_config is not None else {'enabled': False}),
            'pathway_triptych_display': ({
                'enabled': True,
                'databases': list(_PATHWAY_TRIPTYCH_DATABASES),
                'display_cap_per_database': pathway_triptych_config['top_n'],
                'statistics_note': 'KEGG、Reactome 与 WikiPathways Human 保留各自完整库内独立 FDR；三分区图不合并或重新校正 FDR。',
            } if pathway_triptych_config is not None else {'enabled': False}),
        }
        manifest_path = Path(results_dir) / f'enrichment_batch_{batch_key}_manifest.json'
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
        result_files.append({
            'file_path': str(manifest_path), 'file_type': 'json', 'category': 'qc',
            'label': '批量富集执行与统计范围说明',
        })

        if source_paths:
            integration_prefix = f'enrichment_batch_{batch_key}'
            result_files.extend(_write_enrichment_integration(
                results_dir, source_paths=source_paths, plots_dir=plots_dir,
                output_prefix=integration_prefix,
                go_priority_config=go_priority_config,
                pathway_triptych_config=pathway_triptych_config,
            ))

        completed = [row for row in database_rows if row['status'] == 'completed']
        self.progress(100, '批量富集完成')
        summary = {
            'method': manifest['method'],
            'comparison': manifest['comparison'],
            'databases_requested': list(databases),
            'database_results': database_rows,
            'n_databases_completed': len(completed),
            'n_databases_failed': len(databases) - len(completed),
            'n_tested': int(sum(row.get('n_tested') or 0 for row in completed)),
            'n_significant': int(sum(row.get('n_significant') or 0 for row in completed)),
            'fdr_scope': manifest['fdr_scope'],
            'go_priority_display': manifest['go_priority_display'],
            'pathway_triptych_display': manifest['pathway_triptych_display'],
        }
        if not completed:
            summary['error'] = '所有指定基因集数据库均未能完成；请查看批量富集审计文件。'
        return {'output_adata': None, 'result_files': result_files, 'summary': summary}

    def run(self, input_path):
        if (self.params.get('batch_databases') in (True, 'true', 'on', '1')
                and not self.params.get('_batch_subrun')):
            return self._run_database_batch(
                input_path, self._batch_database_names(self.params.get('databases', '')),
            )
        if (_as_enabled(self.params.get('go_priority_enabled', False))
                and not self.params.get('_batch_subrun')):
            raise ValueError(
                '主题优先 GO 三分区图需要启用“批量执行多个数据库”，并同时选择 GO_BP、GO_CC、GO_MF。'
            )
        method = str(self.params.get('method', 'ORA') or 'ORA').upper()
        database = str(self.params.get('database', 'GO_BP') or 'GO_BP')
        organism = str(self.params.get('organism', 'Human') or 'Human')
        if organism.lower() != 'human':
            raise ValueError('当前 Bulk RNA 富集分析仅支持 Human。')
        allowed_databases = {
            'GO_BP', 'GO_MF', 'GO_CC', 'KEGG', 'WikiPathways', 'Reactome', 'Custom_GMT',
        }
        if database not in allowed_databases:
            raise ValueError(f'不支持的 Human 基因集数据库: {database}')
        pvalue_cutoff = _bounded_float(
            self.params.get('pvalue_cutoff', 0.05), '显著性阈值', 0.0, 1.0,
            minimum_inclusive=False,
        )
        top_n = _bounded_int(self.params.get('top_n', 12), '展示通路数', 1, 60)
        input_source = str(self.params.get('input_source', '') or '').strip()
        split_direction = self.params.get('split_direction', False) in (True, 'true', 'on', '1')
        custom_genes_str = str(self.params.get('custom_genes', '') or '').strip()
        custom_background_str = str(self.params.get('custom_background', '') or '').strip()
        custom_geneset_text = str(self.params.get('custom_geneset_text', '') or '').strip()
        custom_geneset_name = str(
            self.params.get('custom_geneset_name', 'Custom Human GMT') or 'Custom Human GMT'
        ).strip()[:120]
        redundancy_threshold = _bounded_float(
            self.params.get('redundancy_threshold', 0.85), '冗余通路 Jaccard 阈值', 0.5, 1.0,
        )
        ora_min_size = _bounded_int(
            self.params.get('ora_min_size', 10), 'ORA 最小基因集', 1, 5000,
        ) if method == 'ORA' else None
        ora_max_size = _bounded_int(
            self.params.get('ora_max_size', 500), 'ORA 最大基因集', 1, 10000,
        ) if method == 'ORA' else None
        if method == 'ORA' and ora_min_size > ora_max_size:
            raise ValueError('ORA 最小基因集不能大于最大基因集。')
        if database == 'Custom_GMT' and not custom_geneset_text:
            raise ValueError('选择 Custom_GMT 时必须粘贴 Human GMT 内容。')
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
            if not os.path.isfile(real_input):
                raise ValueError('所选 DEG 结果文件不存在。')

        if method == 'GSEA' and custom_genes_str:
            raise ValueError('GSEA 需要完整排序基因表，不能使用无排序的自定义基因列表。')

        deg_df = pd.read_csv(input_source) if input_source else None
        prepared_deg = None
        identifier_map = {}
        identifier_provenance = {}
        if deg_df is not None:
            prepared_deg, identifier_map, identifier_provenance = prepare_deg_gene_table(deg_df)

        self.progress(15, f"加载 {database} 基因集数据库...")

        result_files = []
        plots_dir, results_dir = self._output_directories()
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)
        geneset_snapshot_file = None

        # Map database name to file
        db_map = {
            'GO_BP': 'GO_Biological_Process_2023',
            'GO_MF': 'GO_Molecular_Function_2023',
            'GO_CC': 'GO_Cellular_Component_2023',
            'KEGG': 'KEGG_2021_Human',
            'WikiPathways': 'WikiPathways_2019_Human',
            'Reactome': 'Reactome_2022',
        }

        if database == 'Custom_GMT':
            pathways_dict = load_human_genesets_text(
                custom_geneset_text, source_label=custom_geneset_name,
            )
            snapshot_key = _enrichment_output_key(
                method, database, comparison=comparison,
            )
            db_path = str(Path(results_dir) / f'{snapshot_key}_geneset_snapshot.txt')
            Path(db_path).write_text(custom_geneset_text.rstrip() + '\n', encoding='utf-8')
            geneset_snapshot_file = {
                'file_path': db_path, 'file_type': 'txt', 'category': 'data',
                'label': f'自定义 Human GMT 快照 · {custom_geneset_name}',
            }
            geneset_source_type = 'custom_inline_human_gmt'
        else:
            db_filename = db_map.get(database, 'GO_Biological_Process_2023')
            db_path = _resolve_geneset_path(db_filename, organism)
            pathways_dict = load_human_genesets(db_path)
            geneset_source_type = 'local_or_cached_enrichr'
        library_genes = {gene for members in pathways_dict.values() for gene in members}

        self.progress(30, f"运行 {method} 富集分析...")

        if geneset_snapshot_file is not None:
            result_files.append(geneset_snapshot_file)

        n_sig = 0
        n_tested = 0
        n_input_genes = 0
        n_background_genes = 0
        n_ranked_genes = 0
        figure_audits = []
        qc_payload = {
            'method': method,
            'database': database,
            'organism': 'Human',
            'comparison': comparison,
            'input_source': os.path.basename(input_source) if input_source else None,
            'background_comparison': self.params.get('background_comparison'),
            'identifier_mapping': identifier_provenance,
            'gene_set_library': geneset_manifest(
                db_path, pathways_dict, source_type=geneset_source_type,
            ),
            'custom_geneset_name': custom_geneset_name if database == 'Custom_GMT' else None,
            'pvalue_cutoff': pvalue_cutoff,
            'redundancy_threshold': redundancy_threshold,
            'analyses': [],
        }
        analysis_result_frame = None
        leading_edge_frame = None
        full_results_path = None

        if method == 'ORA':
            if custom_background_str:
                background_raw = split_gene_text(custom_background_str)
                background_genes = map_gene_values(background_raw, identifier_map)
                background_source = 'custom_background'
            elif prepared_deg is not None:
                primary_column = identifier_provenance['primary_column']
                background_raw = prepared_deg[primary_column].tolist()
                background_genes = prepared_deg['_analysis_gene'].tolist()
                background_source = 'all_genes_in_deg_table'
            else:
                background_genes, expression_map, expression_provenance = _background_from_analysis_input(input_path)
                identifier_map = {**expression_map, **identifier_map}
                identifier_provenance = expression_provenance
                background_raw = background_genes
                background_source = 'selected_expression_input'

            background_genes = [gene for gene in background_genes if gene]
            if not background_genes:
                raise ValueError('无法构建 ORA 背景；请选择 DEG 结果或提供自定义背景基因。')
            background_qc = mapping_qc(
                background_raw, background_genes, library_genes,
            )
            background_qc['source'] = background_source
            effective_background_genes = set(background_genes).intersection(library_genes)
            background_qc['n_submitted_background_genes'] = len(set(background_genes))
            background_qc['n_effective_background_genes'] = len(effective_background_genes)
            qc_payload['background'] = background_qc
            qc_payload['identifier_mapping'] = identifier_provenance
            n_background_genes = len(effective_background_genes)

            if prepared_deg is None and not custom_genes_str:
                raise ValueError('ORA 需要含 regulation 列的 DEG 结果或自定义基因列表。')

            primary_column = identifier_provenance.get('primary_column', 'gene')

            def ora_result(query_raw, query_genes, direction=None):
                query_genes = [gene for gene in query_genes if gene]
                query_qc = mapping_qc(
                    query_raw, query_genes, library_genes,
                    background=effective_background_genes,
                )
                full = run_ora_full(
                    query_genes, background_genes, pathways_dict,
                    pvalue_cutoff=pvalue_cutoff,
                    min_size=ora_min_size, max_size=ora_max_size,
                )
                testing_family = dict(full.attrs.get('testing_family', {}))
                full['MappedInputGeneCount'] = query_qc['n_mapped_in_background']
                full['MappedBackgroundGeneCount'] = background_qc['n_mapped_to_library']
                full['InputMappingRate'] = query_qc['mapping_rate']
                full['BackgroundMappingRate'] = background_qc['mapping_rate']
                full = annotate_redundancy(full, threshold=redundancy_threshold)
                full = _add_enrichment_metadata(
                    full, method, database, direction=direction, comparison=comparison,
                )
                assigned = full['RedundancyCluster'].astype(str).str.strip() != ''
                representatives = assigned & full['IsRepresentative'].astype(bool)
                qc_payload['analyses'].append({
                    'direction': direction or 'All',
                    'query': query_qc,
                    'testing_family': testing_family,
                    'n_tested_pathways': int(len(full)),
                    'n_significant_pathways': int(full['Significant'].sum()),
                    'n_redundancy_clusters': int(representatives.sum()),
                    'n_redundant_significant_pathways': int(assigned.sum() - representatives.sum()),
                })
                return full

            if split_direction and not custom_genes_str:
                if 'regulation' not in prepared_deg.columns:
                    raise ValueError('分方向 ORA 需要 DEG 结果包含 regulation 列。')
                # Split direction: run ORA separately for Up and Down genes
                all_split_results = []
                regulation = prepared_deg['regulation'].astype(str).str.strip().str.lower()

                for direction in ['Up', 'Down']:
                    mask = regulation == direction.lower()
                    direction_genes = prepared_deg.loc[mask, '_analysis_gene'].tolist()
                    direction_raw = prepared_deg.loc[mask, primary_column].tolist()
                    if not direction_genes:
                        continue

                    self.progress(45 if direction == 'Up' else 60, f"运行 {direction} 基因 ORA...")
                    enr_dir = ora_result(direction_raw, direction_genes, direction=direction)
                    all_split_results.append(enr_dir)

                    output_key = _enrichment_output_key(
                        method, database, direction, comparison=comparison)
                    plot_enr_dir = _plot_enrichment_rows(enr_dir, self.params)
                    if len(plot_enr_dir) > 0:
                        fig_enrichment = _enrichment_figure(
                            plot_enr_dir,
                            title=f'{comparison} · {database} ORA · {direction}-regulated genes',
                            database=database,
                            params=self.params,
                        )
                        if fig_enrichment is not None:
                            exported, audit = _save_enrichment_figure(
                                self, fig_enrichment, plots_dir, results_dir, output_key,
                                f'{comparison} · {database} ORA 富集图（{direction}）',
                            )
                            result_files.extend(exported)
                            figure_audits.append(audit)
                        variant_files, variant_audits = _render_enrichment_variants(
                            self, plot_enr_dir, method='ORA',
                            title=f'{comparison} · {database} ORA · {direction}-regulated genes',
                            base_output_key=output_key,
                            label=f'{comparison} · {database} ORA 富集图（{direction}）',
                            params=self.params,
                        )
                        result_files.extend(variant_files)
                        figure_audits.extend(variant_audits)

                if all_split_results:
                    combined = pd.concat(all_split_results, ignore_index=True)
                    output_key = _enrichment_output_key(
                        method, database, 'directional', comparison=comparison)
                    csv_path = os.path.join(results_dir, f'{output_key}_results.csv')
                    combined.to_csv(csv_path, index=False)
                    full_results_path = csv_path
                    result_files.append({
                        'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                        'label': f'{comparison} · {database} 分方向 ORA 富集结果',
                    })
                    n_sig = int(combined['Significant'].sum())
                    n_tested = int(len(combined))
                    n_input_genes = int(sum(
                        item['query']['n_in_background'] for item in qc_payload['analyses']
                    ))
                    analysis_result_frame = combined
                else:
                    raise ValueError('DEG 结果没有可用于分方向 ORA 的 Up/Down 基因。')

            else:
                if custom_genes_str:
                    query_raw = split_gene_text(custom_genes_str)
                    query_genes = map_gene_values(query_raw, identifier_map)
                else:
                    if 'regulation' not in prepared_deg.columns:
                        raise ValueError('ORA 需要 DEG 结果包含 regulation 列以识别显著基因。')
                    regulation = prepared_deg['regulation'].astype(str).str.strip().str.lower()
                    mask = regulation.isin({'up', 'down'})
                    query_raw = prepared_deg.loc[mask, primary_column].tolist()
                    query_genes = prepared_deg.loc[mask, '_analysis_gene'].tolist()
                if not query_genes:
                    raise ValueError('未找到可用于 ORA 的显著基因。')
                enr = ora_result(query_raw, query_genes)

                self.progress(60, "保存结果表...")
                output_key = _enrichment_output_key(method, database, comparison=comparison)
                csv_path = os.path.join(results_dir, f'{output_key}_results.csv')
                enr.to_csv(csv_path, index=False)
                full_results_path = csv_path
                result_files.append({
                    'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                    'label': f'{comparison} · {database} ORA 富集结果',
                })

                self.progress(70, "生成气泡图...")
                plot_enr = _plot_enrichment_rows(enr, self.params)
                if len(plot_enr) > 0:
                    fig_enrichment = _enrichment_figure(
                        plot_enr,
                        title=f'{comparison} · {database} ORA enrichment',
                        database=database,
                        params=self.params,
                    )
                    if fig_enrichment is not None:
                        exported, audit = _save_enrichment_figure(
                            self, fig_enrichment, plots_dir, results_dir, output_key,
                            f'{comparison} · {database} ORA 富集图',
                        )
                        result_files.extend(exported)
                        figure_audits.append(audit)
                    variant_files, variant_audits = _render_enrichment_variants(
                        self, plot_enr, method='ORA',
                        title=f'{comparison} · {database} ORA enrichment',
                        base_output_key=output_key,
                        label=f'{comparison} · {database} ORA 富集图',
                        params=self.params,
                    )
                    result_files.extend(variant_files)
                    figure_audits.extend(variant_audits)

                n_sig = int(enr['Significant'].sum())
                n_tested = int(len(enr))
                n_input_genes = int(qc_payload['analyses'][0]['query']['n_in_background'])
                analysis_result_frame = enr

        elif method == 'GSEA':
            if prepared_deg is None:
                raise ValueError('GSEA 需要选择一个完整 DEG 比较结果。')
            ranking_metric = str(self.params.get('ranking_metric', 'auto') or 'auto')
            gene_rnk, ranking_provenance = prepare_gsea_ranking(
                prepared_deg, ranking_metric=ranking_metric,
                identifier_map=identifier_map,
            )
            permutation_num = _bounded_int(
                self.params.get('permutation_num', 1000), 'GSEA permutation 次数', 100, 5000,
            )
            gsea_seed = _bounded_int(self.params.get('gsea_seed', 112), 'GSEA 随机种子', 0, 2147483647)
            gsea_weight = _bounded_float(self.params.get('gsea_weight', 1.0), 'GSEA weight', 0.0, 2.0)
            min_size = _bounded_int(self.params.get('gsea_min_size', 15), 'GSEA 最小基因集', 2, 5000)
            max_size = _bounded_int(self.params.get('gsea_max_size', 500), 'GSEA 最大基因集', 2, 10000)
            if min_size > max_size:
                raise ValueError('GSEA 最小基因集不能大于最大基因集。')
            primary_column = identifier_provenance['primary_column']
            rank_qc = mapping_qc(
                prepared_deg[primary_column].tolist(), gene_rnk['gene_name'].tolist(),
                library_genes,
            )
            rank_qc.update(ranking_provenance)
            ranked_gene_set = set(gene_rnk['gene_name'])
            rank_qc['n_eligible_pathways'] = sum(
                min_size <= len(set(members).intersection(ranked_gene_set)) <= max_size
                for members in pathways_dict.values()
            )
            if rank_qc['n_mapped_to_library'] == 0 or rank_qc['n_eligible_pathways'] == 0:
                raise ValueError('GSEA 排名基因与所选数据库没有满足大小阈值的通路交集。')
            qc_payload['ranking'] = rank_qc
            qc_payload['gsea_parameters'] = {
                'ranking_metric': ranking_provenance['ranking_metric_used'],
                'permutation_num': permutation_num,
                'seed': gsea_seed,
                'weight': gsea_weight,
                'min_size': min_size,
                'max_size': max_size,
                'backend': 'omicverse_numpy_prerank_direct',
            }
            try:
                from importlib.metadata import version
                qc_payload['gsea_parameters']['omicverse_version'] = version('omicverse')
            except Exception:
                qc_payload['gsea_parameters']['omicverse_version'] = 'unknown'

            pre_res = run_gsea_prerank(
                gene_rnk,
                pathways_dict,
                permutation_num=permutation_num,
                seed=gsea_seed,
                weight=gsea_weight,
                min_size=min_size,
                max_size=max_size,
            )

            self.progress(71, "处理 GSEA 结果...")
            enr = _ensure_gsea_term_column(pre_res.res2d)
            if enr.empty or 'fdr' not in enr.columns:
                raise ValueError('GSEA 没有产生可检验通路；请检查映射覆盖率和基因集大小阈值。')
            enr['Significant'] = pd.to_numeric(enr['fdr'], errors='coerce') < pvalue_cutoff
            enr = add_gsea_leading_edge_metrics(enr)
            enr = annotate_redundancy(enr, threshold=redundancy_threshold)
            enr = _add_enrichment_metadata(enr, method, database, comparison=comparison)
            enr['RankingMetric'] = ranking_provenance['ranking_metric_used']
            enr['PermutationNum'] = permutation_num
            enr['Seed'] = gsea_seed
            enr['Weight'] = gsea_weight
            enr['MinSize'] = min_size
            enr['MaxSize'] = max_size
            leading_edge_frame = gsea_leading_edge_table(enr, gene_rnk)
            output_key = _enrichment_output_key(method, database, comparison=comparison)
            csv_path = os.path.join(results_dir, f'{output_key}_results.csv')
            enr.to_csv(csv_path, index=False)
            full_results_path = csv_path
            result_files.append({
                'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                'label': f'{comparison} · {database} GSEA 富集结果',
            })
            leading_edge_path = os.path.join(
                results_dir, f'{output_key}_leading_edge.csv',
            )
            leading_edge_frame.to_csv(leading_edge_path, index=False)
            result_files.append({
                'file_path': leading_edge_path, 'file_type': 'csv', 'category': 'table',
                'label': f'{comparison} · {database} GSEA leading-edge 明细',
            })

            self.progress(75, "生成 GSEA 图表...")
            plot_enr = _plot_enrichment_rows(enr, self.params)
            if len(plot_enr) > 0:
                fig_enrichment = _enrichment_figure(
                    plot_enr,
                    title=f'{comparison} · {database} GSEA enrichment',
                    database=database,
                    score_column='nes' if 'nes' in plot_enr.columns else None,
                    params=self.params,
                )
                if fig_enrichment is not None:
                    exported, audit = _save_enrichment_figure(
                        self, fig_enrichment, plots_dir, results_dir, output_key,
                        f'{comparison} · {database} GSEA 富集图',
                    )
                    result_files.extend(exported)
                    figure_audits.append(audit)
                variant_files, variant_audits = _render_enrichment_variants(
                    self, plot_enr, method='GSEA',
                    title=f'{comparison} · {database} GSEA enrichment',
                    base_output_key=output_key,
                    label=f'{comparison} · {database} GSEA 富集图',
                    params=self.params,
                    pre_res=pre_res,
                )
                result_files.extend(variant_files)
                figure_audits.extend(variant_audits)

            n_sig = int(enr['Significant'].sum())
            n_tested = int(len(enr))
            n_ranked_genes = int(len(gene_rnk))
            n_input_genes = n_ranked_genes
            analysis_result_frame = enr
            assigned = enr['RedundancyCluster'].astype(str).str.strip() != ''
            representatives = assigned & enr['IsRepresentative'].astype(bool)
            qc_payload['gsea_parameters'].update({
                'n_redundancy_clusters': int(representatives.sum()),
                'n_redundant_significant_pathways': int(
                    assigned.sum() - representatives.sum()
                ),
                'n_leading_edge_rows': int(len(leading_edge_frame)),
            })
        else:
            raise ValueError(f"不支持的富集方法: {method}")

        self.progress(90, "保存输出...")

        warnings = _mapping_warnings(qc_payload)
        for item in qc_payload.get('analyses', []):
            warnings.extend(_mapping_warnings({'query': item.get('query', {})}))
        qc_payload['warnings'] = list(dict.fromkeys(warnings))
        qc_output_key = _enrichment_output_key(method, database, comparison=comparison)

        if analysis_result_frame is None:
            raise RuntimeError('富集分析未生成可交付的完整结果表。')
        redundancy_frame = _redundancy_audit_table(analysis_result_frame)
        redundancy_path = os.path.join(
            results_dir, f'{qc_output_key}_redundancy_clusters.csv',
        )
        redundancy_frame.to_csv(redundancy_path, index=False)
        result_files.append({
            'file_path': redundancy_path, 'file_type': 'csv', 'category': 'table',
            'label': f'{comparison} · 显著通路冗余簇审计',
        })
        methods_file = _write_enrichment_methods(
            results_dir, qc_output_key, qc_payload, analysis_result_frame,
        )
        result_files.append(methods_file)
        source_data_file = _write_enrichment_source_data(
            results_dir, qc_output_key, analysis_result_frame, qc_payload,
            redundancy_df=redundancy_frame,
            leading_edge_df=leading_edge_frame,
        )
        if source_data_file is not None:
            result_files.append(source_data_file)
        else:
            qc_payload['warnings'].append(
                '缺少 openpyxl，未生成 Source Data XLSX；完整 CSV 与 QC JSON 已保留。'
            )
        qc_payload['deliverables'] = {
            'full_results': os.path.basename(full_results_path) if full_results_path else None,
            'redundancy_audit': os.path.basename(redundancy_path),
            'leading_edge': (
                f'{qc_output_key}_leading_edge.csv' if method == 'GSEA' else None
            ),
            'source_data_workbook': (
                os.path.basename(source_data_file['file_path']) if source_data_file else None
            ),
            'methods': os.path.basename(methods_file['file_path']),
        }
        result_files.append(_write_mapping_qc(results_dir, qc_output_key, qc_payload))

        # Results created before contrast-aware enrichment stored no source
        # metadata.  Repair registered legacy tables before building the shared
        # overview so an old GO result and a new KEGG result can be grouped.
        if not self.params.get('_batch_subrun'):
            _backfill_legacy_enrichment_metadata(self.project_dir, results_dir)

            # Rebuild the legacy project-level table after a single-library
            # run. Batch runs integrate only their explicit child result paths
            # in _run_database_batch(), never the whole project directory.
            result_files.extend(_write_enrichment_integration(results_dir))

        # Clean up temporary enrichment directories
        import shutil
        for tmp_suffix in ['_tmp', '_up_tmp', '_down_tmp', '_gsea_tmp']:
            tmp_dir = os.path.join(self.project_dir, f'enrichr{tmp_suffix}')
            if os.path.isdir(tmp_dir):
                shutil.rmtree(tmp_dir, ignore_errors=True)

        self.progress(100, "完成")

        if method == 'ORA':
            input_mapping_rates = {
                item.get('direction', 'All'): item.get('query', {}).get('mapping_rate')
                for item in qc_payload.get('analyses', [])
            }
            mapping_rates = {
                'background': qc_payload.get('background', {}).get('mapping_rate'),
                'inputs': input_mapping_rates,
            }
        else:
            mapping_rates = {
                'ranking': qc_payload.get('ranking', {}).get('mapping_rate'),
            }

        return {
            'output_adata': None,
            'result_files': result_files,
            'summary': {
                'method': method,
                'database': database,
                'comparison': comparison,
                'organism': organism,
                'n_input_genes': n_input_genes,
                'n_background_genes': n_background_genes,
                'n_ranked_genes': n_ranked_genes,
                'n_tested': n_tested,
                'n_significant': n_sig,
                'pvalue_cutoff': pvalue_cutoff,
                'input_mode': 'custom_genes' if custom_genes_str else ('split_direction' if split_direction else 'standard'),
                'mapping_rates': mapping_rates,
                'mapping_warnings': qc_payload['warnings'],
                'gene_set_sha256': qc_payload['gene_set_library']['sha256'],
                'gene_set_library': qc_payload['gene_set_library']['library_name'],
                'gene_set_library_year': qc_payload['gene_set_library']['library_year'],
                'n_redundancy_clusters': int(
                    len(redundancy_frame.drop_duplicates(
                        [column for column in ('Direction', 'RedundancyCluster')
                         if column in redundancy_frame.columns]
                    )) if not redundancy_frame.empty else 0
                ),
                'n_leading_edge_rows': int(
                    len(leading_edge_frame) if leading_edge_frame is not None else 0
                ),
                # Used internally by the parent batch coordinator. The worker
                # rewrites this path when it snapshots declared result files.
                'full_results_path': full_results_path,
                'figure_layout_qc': {
                    'n_figures': len(figure_audits),
                    'n_figures_with_overlap': sum(
                        bool(audit.get('n_overlap_pairs', 0) > 0 or audit.get('out_of_bounds_labels'))
                        for audit in figure_audits
                    ),
                },
            }
        }
