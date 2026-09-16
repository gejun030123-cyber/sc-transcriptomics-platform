import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


def _save_native_figure(fig, plots_dir, stem, category, label, viz_params=None, dpi=300):
    """Save one Matplotlib canvas as both PNG and editable SVG."""
    from modules.figure_style import apply_matplotlib_style, _font_for_text, NATURE_FONT_FAMILY

    viz = viz_params or {}
    apply_matplotlib_style(fig, viz)
    for ax in getattr(fig, 'axes', []):
        ax.set_facecolor(viz.get('bg_color', 'white'))
        ax.tick_params(labelsize=max(8, float(viz.get('font_size', 12)) - 2), width=0.7)
        ax.title.set_fontfamily(_font_for_text(ax.title.get_text(), viz.get('font_family', NATURE_FONT_FAMILY)))
    result_files = []
    for fmt in ('png', 'svg'):
        path = os.path.join(plots_dir, f'{stem}.{fmt}')
        fig.savefig(path, format=fmt, dpi=dpi if fmt == 'png' else None,
                    bbox_inches='tight', pad_inches=0.15, facecolor='white')
        result_files.append({'file_path': path, 'file_type': fmt,
                             'category': category, 'label': label})
    try:
        import matplotlib.pyplot as plt
        plt.close(fig)
    except Exception:
        pass
    return result_files


VOLCANO_COLORS = {
    'Up': '#C65A5A',
    'Down': '#4C78A8',
    'NS': '#D9DEE5',
}


def _format_volcano_contrast_title(group1, group2):
    """Make common condition/stratum group names readable without changing IDs."""
    aliases = {'ctr': 'Control', 'ctrl': 'Control', 'control': 'Control'}

    def split_group(value):
        text = str(value or '').strip()
        condition, separator, stratum = text.rpartition('_')
        if separator and condition and stratum.casefold() in {'b', 'en'}:
            return condition, stratum
        return text, ''

    def display_condition(value):
        text = str(value or '').strip()
        display = aliases.get(text.casefold(), text)
        return display.replace('NH4Cl', 'NH₄Cl')

    condition1, stratum1 = split_group(group1)
    condition2, stratum2 = split_group(group2)
    if stratum1 and stratum1 == stratum2:
        return f'{display_condition(condition1)} vs {display_condition(condition2)} ({stratum1})'
    return f'{display_condition(group1)} vs {display_condition(group2)}'


def _volcano_y_limit(padj, pval_threshold=0.05, max_display=18.0):
    """Choose a linear, compact FDR range without changing stored statistics."""
    padj = np.asarray(padj, dtype=float)
    safe_padj = np.clip(np.nan_to_num(padj, nan=1.0, posinf=1.0, neginf=1.0),
                        np.finfo(float).tiny, 1.0)
    y = -np.log10(safe_padj)
    finite = y[np.isfinite(y)]
    threshold_y = -np.log10(max(float(pval_threshold), np.finfo(float).tiny))
    if finite.size == 0:
        return float(max(4.0, threshold_y + 1.5))
    # Do not apply a nonlinear tail transform: points above the fixed ceiling
    # are explicitly marked with triangles by _draw_bulk_volcano.
    if float(np.nanmax(finite)) > float(max_display):
        return float(max_display)
    return float(max(4.0, threshold_y + 1.5, np.ceil(np.nanmax(finite) + 0.6)))


def _volcano_label_positions(label_rows, y_min, y_max, x_min, x_max):
    """Seed labels with short angled leaders before the repel pass."""
    if not label_rows:
        return []

    x_span = float(x_max) - float(x_min)
    y_span = float(y_max) - float(y_min)
    x_midpoint = (float(x_min) + float(x_max)) / 2.0
    positions = []
    for direction, side in (('Up', 1), ('Down', -1)):
        side_rows = [row for row in label_rows if row['regulation'] == direction]
        for rank, row in enumerate(sorted(side_rows, key=lambda item: item['y_point'])):
            label_side = side
            # When the point is at either displayed x edge, point the label
            # inwards.  This works for a deliberately asymmetric user range.
            near_left = float(row['x_point']) <= float(x_min) + x_span * 0.22
            near_right = float(row['x_point']) >= float(x_max) - x_span * 0.22
            if near_left:
                label_side = 1
            elif near_right:
                label_side *= -1
            elif (float(row['x_point']) - x_midpoint) * label_side < 0:
                label_side *= -1
            y_side = -1 if float(row['y_point']) >= float(y_min) + y_span * 0.84 else 1
            positions.append({
                **row,
                'x_text': float(row['x_point']) + label_side * x_span * 0.08,
                'y_text': float(row['y_point']) + y_side * y_span * (0.05 + 0.02 * (rank // 2)),
                'label_side': label_side,
                'y_side': y_side,
            })

    return positions


def _draw_bulk_volcano(ax, deg_df, pval_threshold=0.05, fc_threshold=2.0,
                       title='', top_n=10, show_legend=True, y_limit=None,
                       x_limit=None, x_min=None, x_max=None, y_min=None,
                       colors=None):
    """Draw a restrained, readable Bulk RNA-seq volcano panel."""
    from matplotlib.lines import Line2D
    from modules.figure_style import NATURE_AXIS, NATURE_TEXT

    palette = {**VOLCANO_COLORS, **(colors or {})}
    log2fc = pd.to_numeric(
        pd.Series(deg_df.get('log2FC', pd.Series(dtype=float))), errors='coerce'
    ).to_numpy(dtype=float)
    padj = pd.to_numeric(
        pd.Series(deg_df.get('padj', pd.Series(dtype=float))), errors='coerce'
    ).to_numpy(dtype=float)
    regulation = deg_df.get(
        'regulation', pd.Series(index=deg_df.index, dtype=str)
    ).astype(str).to_numpy()
    safe_padj = np.clip(np.nan_to_num(padj, nan=1.0, posinf=1.0, neginf=1.0),
                        np.finfo(float).tiny, 1.0)
    y_raw = -np.log10(safe_padj)
    if y_limit is None:
        y_limit = _volcano_y_limit(padj, pval_threshold)
    y_limit = float(y_limit)
    try:
        y_min = max(0.0, float(y_min)) if y_min is not None else 0.0
    except (TypeError, ValueError):
        y_min = 0.0
    if y_min >= y_limit:
        y_min = 0.0
    y_cap_position = y_min + (y_limit - y_min) * 0.985
    y_clipped = np.isfinite(y_raw) & (y_raw > y_limit)
    y_display = np.minimum(y_raw, y_cap_position)
    log2fc_threshold = np.log2(max(float(fc_threshold), np.finfo(float).tiny))

    if x_limit is None:
        x_limit = max(2.5, float(np.ceil((log2fc_threshold + 0.25) * 2.0) / 2.0))
    x_limit = float(x_limit)
    try:
        x_min = float(x_min) if x_min is not None else -x_limit
        x_max = float(x_max) if x_max is not None else x_limit
    except (TypeError, ValueError):
        x_min, x_max = -x_limit, x_limit
    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min >= x_max:
        x_min, x_max = -x_limit, x_limit
    x_span = x_max - x_min
    # Keep the historical ±2.5 display at 98.5% of its edge while extending
    # the same relative inboard margin to asymmetric user-selected ranges.
    x_edge_inset = x_span * 0.0075
    x_clipped = np.isfinite(log2fc) & ((log2fc < x_min) | (log2fc > x_max))
    x_display = np.clip(log2fc, x_min + x_edge_inset, x_max - x_edge_inset)

    masks = {
        'NS': regulation == 'NS',
        'Down': regulation == 'Down',
        'Up': regulation == 'Up',
    }
    # Draw the dense neutral cloud first; significant points remain visible.
    for label, size, alpha, zorder in (
        ('NS', 4, 0.26, 1), ('Down', 9, 0.90, 2), ('Up', 9, 0.90, 2),
    ):
        mask = (masks[label] & np.isfinite(log2fc) & np.isfinite(y_display)
                & ~(x_clipped | y_clipped))
        if mask.any():
            ax.scatter(x_display[mask], y_display[mask], s=size,
                       color=palette[label], alpha=alpha,
                       linewidths=0, rasterized=True, zorder=zorder)

    # Triangles make both display-only truncations explicit instead of
    # widening the panel or applying a nonlinear FDR-tail transform.
    for label in ('NS', 'Down', 'Up'):
        label_mask = masks[label]
        vertical = label_mask & y_clipped
        if vertical.any():
            ax.scatter(x_display[vertical], y_display[vertical], s=17,
                       color=palette[label], marker='^', alpha=0.92,
                       linewidths=0, rasterized=True, zorder=4)
        left = label_mask & x_clipped & ~y_clipped & (log2fc < x_min)
        right = label_mask & x_clipped & ~y_clipped & (log2fc > x_max)
        if left.any():
            ax.scatter(x_display[left], y_display[left], s=17,
                       color=palette[label], marker='<', alpha=0.92,
                       linewidths=0, rasterized=True, zorder=4)
        if right.any():
            ax.scatter(x_display[right], y_display[right], s=17,
                       color=palette[label], marker='>', alpha=0.92,
                       linewidths=0, rasterized=True, zorder=4)

    ax.axhline(-np.log10(max(float(pval_threshold), np.finfo(float).tiny)),
               color='#D9DEE5', linestyle='--', linewidth=0.65, alpha=0.60, zorder=3)
    ax.axvline(log2fc_threshold, color='#D9DEE5', linestyle='--', linewidth=0.65, alpha=0.60, zorder=3)
    ax.axvline(-log2fc_threshold, color='#D9DEE5', linestyle='--', linewidth=0.65, alpha=0.60, zorder=3)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_limit)
    ax.set_title(title, loc='left', pad=8, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel(r'log$_2$(fold change)', fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(r'$-\log_{10}(\mathrm{FDR})$', fontsize=9, color=NATURE_TEXT)
    ax.grid(False)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, colors=NATURE_AXIS, width=0.7, length=3)
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name in ('left', 'bottom'))
        spine.set_color(NATURE_AXIS)
        spine.set_linewidth(0.7)

    if show_legend:
        counts = {label: int(mask.sum()) for label, mask in masks.items()}
        handles = [
            Line2D([], [], marker='o', linestyle='None', markersize=5,
                   markerfacecolor=palette['Up'], markeredgewidth=0,
                   label=f"Up  {counts['Up']:,}"),
            Line2D([], [], marker='o', linestyle='None', markersize=5,
                   markerfacecolor=palette['Down'], markeredgewidth=0,
                   label=f"Down  {counts['Down']:,}"),
            Line2D([], [], marker='o', linestyle='None', markersize=5,
                   markerfacecolor=palette['NS'], markeredgewidth=0,
                   label=f"NS  {counts['NS']:,}"),
        ]
        ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                  frameon=False, fontsize=7.5, handlelength=1.3,
                  borderaxespad=0.0)

    if top_n > 0 and not deg_df.empty:
        sig = deg_df[deg_df['regulation'].isin(['Up', 'Down'])].copy()
        readable_sig = sig[~sig['gene'].astype(str).str.match(
            r'^(?:ENS(?:G|MUSG|DARG|RNOG)\d+|[A-Za-z]+\d{8,})$', na=False,
        )]
        if readable_sig.empty:
            readable_sig = sig
        # ``top_n`` controls the result table, not visual density.  Label a
        # compact mix of FDR-leading genes and the largest effects so the
        # panel does not imply that FDR rank is the only biological priority.
        annotation_limit = min(max(int(top_n), 0), 6)
        fdr_per_side = 2 if annotation_limit >= 4 else 1
        selected = pd.concat([
            readable_sig[readable_sig['regulation'] == 'Up'].nsmallest(fdr_per_side, 'padj'),
            readable_sig[readable_sig['regulation'] == 'Down'].nsmallest(fdr_per_side, 'padj'),
            readable_sig.assign(_abs_fc=readable_sig['log2FC'].abs())
            .sort_values(['_abs_fc', 'padj'], ascending=[False, True]),
        ]).drop_duplicates(subset=['gene']).head(annotation_limit)
        label_rows = []
        for _, row in selected.iterrows():
            x_value = float(row['log2FC'])
            p_value = float(row['padj'])
            y_value = float(-np.log10(np.clip(max(p_value, np.finfo(float).tiny),
                                              np.finfo(float).tiny, 1.0)))
            label_rows.append({
                'gene': str(row['gene']),
                'regulation': str(row['regulation']),
                'x_point': float(np.clip(x_value, x_min + x_edge_inset, x_max - x_edge_inset)),
                'y_point': min(y_value, y_cap_position),
            })

        texts = []
        for label in _volcano_label_positions(label_rows, y_min, y_limit, x_min, x_max):
            texts.append(ax.annotate(
                label['gene'], (label['x_point'], label['y_point']),
                xytext=(label['x_text'], label['y_text']), textcoords='data',
                fontsize=7, ha='left' if label['label_side'] > 0 else 'right',
                va='bottom' if label['y_side'] > 0 else 'top',
                color=NATURE_TEXT, zorder=5,
                arrowprops={
                    'arrowstyle': '-', 'color': '#667085', 'linewidth': 0.5,
                    'shrinkA': 1, 'shrinkB': 1,
                },
            ))
        if texts:
            from figure_engine.templates.common import (
                adjust_labels, ensure_angled_annotation_leaders, prune_overlapping_labels,
            )
            adjust_labels(texts, ax, arrow_color='#667085')
            ensure_angled_annotation_leaders(texts, y_limit)
            prune_overlapping_labels(texts, ax)

    return {'y_limit': float(y_limit), 'y_min': float(y_min),
            'x_limit': float(max(abs(x_min), abs(x_max))),
            'x_min': float(x_min), 'x_max': float(x_max),
            'n_clipped': int(y_clipped.sum()), 'n_x_clipped': int(x_clipped.sum())}


def _parse_comparisons(comp_str):
    """解析 'A-vs-B;C-vs-D' 为 [('A','B'), ('C','D')]"""
    if not comp_str or not comp_str.strip():
        return []
    pairs = []
    for item in comp_str.replace('\n', ';').split(';'):
        item = item.strip()
        if '-vs-' in item:
            parts = item.split('-vs-', 1)
            if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                pairs.append((parts[0].strip(), parts[1].strip()))
    return pairs


def _parse_custom_groups(cg_str):
    """解析 'High=Treated_1h+Treated_3h\nLow=Ctrl' 为 {'High': ['Treated_1h','Treated_3h'], 'Low': ['Ctrl']}"""
    if not cg_str or not cg_str.strip():
        return {}
    mapping = {}
    for line in cg_str.strip().split('\n'):
        line = line.strip()
        if '=' not in line:
            continue
        name, expr = line.split('=', 1)
        name = name.strip()
        members = [m.strip() for m in expr.split('+') if m.strip()]
        if name and members:
            mapping[name] = members
    return mapping


def _json_safe(value):
    """Convert numpy scalars/containers read back from ``h5ad`` to plain Python.

    ``adata.uns`` values round-trip through HDF5 as numpy types.  Summaries are
    written to the analysis manifest with a strict ``json.dump``, so an
    un-coerced ``int64`` silently loses the whole manifest.
    """
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    return value


def _preserved_raw_counts(adata):
    """Return raw counts saved by the Bulk normalization handoff, if valid."""
    from modules.io_utils import preserved_raw_count_layer
    return preserved_raw_count_layer(adata)


def _pre_deseq2_gene_filter_provenance(adata):
    """Report whether the low-expression gene filter ran before count-based DEG.

    ``bulk_deg`` consumes the Bulk normalization handoff, which is documented as
    the step that removes low-expression genes.  Without this record an
    unchanged gene count (for example QC's ``genes_removed = 0``) can be
    mistaken for "no filtering needed", so the decision is carried into the DEG
    summary instead of living only in the upstream task's metadata.
    """
    uns = getattr(adata, 'uns', None)
    normalization = uns.get('normalization') if hasattr(uns, 'get') else None
    record = normalization.get('gene_expression_filter') if isinstance(normalization, dict) else None
    n_genes = int(adata.n_vars)
    if isinstance(record, dict):
        applied = bool(record.get('applied'))
        n_before = record.get('n_genes_before')
        n_after = record.get('n_genes_after')
        fallback = (
            f'上游 Bulk 标准化已过滤低表达基因（{n_before} → {n_after}）。'
            if n_before is not None and n_after is not None else
            '上游 Bulk 标准化已过滤低表达基因。'
        )
        return {
            'source': 'bulk_normalize',
            'applied_before_deg': applied,
            'n_genes_tested': n_genes,
            'n_genes_removed_upstream': _json_safe(record.get('n_genes_removed')),
            'record': _json_safe(record),
            'message': (
                (record.get('message') or fallback) if applied else
                '上游 Bulk 标准化未过滤低表达基因；DESeq2/edgeR/limma 将使用全部 '
                f'{n_genes} 个基因，建议先设置 min_expr_samples≥3 与 min_expr_value。'
            ),
        }
    return {
        'source': 'none',
        'applied_before_deg': False,
        'n_genes_tested': n_genes,
        'n_genes_removed_upstream': None,
        'record': None,
        'message': (
            '未检测到低表达基因过滤记录；DESeq2/edgeR/limma 将使用全部 '
            f'{n_genes} 个基因。建议先运行 Bulk 标准化（默认 CPM≥1 且至少 3 个样本表达）。'
        ),
    }


def _welch_ttest_log_expression(data, group1_samples, group2_samples, padj_method='fdr_bh'):
    """Run Welch t-test on a genes x samples matrix already on a log2 scale."""
    from scipy.stats import ttest_ind
    from statsmodels.stats.multitest import multipletests

    g1 = data[group1_samples].to_numpy(dtype=float)
    g2 = data[group2_samples].to_numpy(dtype=float)
    _, pvalues = ttest_ind(g1, g2, axis=1, equal_var=False, nan_policy='omit')
    pvalues = np.nan_to_num(np.asarray(pvalues), nan=1.0, posinf=1.0, neginf=1.0)
    padj = multipletests(pvalues, method=padj_method)[1]
    return pd.DataFrame({
        'pvalue': pvalues,
        'qvalue': padj,
        'log2FC': np.nanmean(g1, axis=1) - np.nanmean(g2, axis=1),
    }, index=data.index)


def _comparison_filter_diagnostics(deg_df, fc_threshold, padj_threshold):
    """Explain, per comparison, why genes did or did not become DEG calls."""
    log2fc_threshold = float(np.log2(max(float(fc_threshold), np.finfo(float).tiny)))
    pvalues = pd.to_numeric(deg_df.get('pvalue', pd.Series(dtype=float)), errors='coerce').to_numpy(dtype=float)
    padj = pd.to_numeric(deg_df.get('padj', pd.Series(dtype=float)), errors='coerce').to_numpy(dtype=float)
    log2fc = pd.to_numeric(deg_df.get('log2FC', pd.Series(dtype=float)), errors='coerce').to_numpy(dtype=float)

    finite_pvalue = np.isfinite(pvalues)
    finite_padj = np.isfinite(padj)
    finite_log2fc = np.isfinite(log2fc)
    # Prefer the unrounded, exact gate calls retained by _run_single_comparison.
    # Falling back to numerical reconstruction keeps this helper useful for
    # legacy result CSVs that predate these columns.
    if 'passes_padj' in deg_df.columns:
        passes_padj = deg_df['passes_padj'].fillna(False).astype(bool).to_numpy()
    else:
        passes_padj = finite_padj & (padj < float(padj_threshold))
    if 'passes_fc' in deg_df.columns:
        passes_fc = deg_df['passes_fc'].fillna(False).astype(bool).to_numpy()
    else:
        passes_fc = finite_log2fc & (np.abs(log2fc) >= log2fc_threshold)
    if 'significant' in deg_df.columns:
        both = deg_df['significant'].fillna(False).astype(bool).to_numpy()
    else:
        both = passes_padj & passes_fc

    def _minimum(values, mask):
        return float(np.min(values[mask])) if np.any(mask) else None

    abs_log2fc = np.abs(log2fc[finite_log2fc])
    return {
        'n_tested_genes': int(len(deg_df)),
        'n_raw_p_lt_0_05': int(np.sum(finite_pvalue & (pvalues < 0.05))),
        'n_raw_p_lt_threshold': int(np.sum(finite_pvalue & (pvalues < float(padj_threshold)))),
        'n_padj_lt_threshold': int(np.sum(passes_padj)),
        'n_abs_log2fc_ge_threshold': int(np.sum(passes_fc)),
        'n_both_padj_and_fc': int(np.sum(both)),
        'min_pvalue': _minimum(pvalues, finite_pvalue),
        'min_padj': _minimum(padj, finite_padj),
        'max_abs_log2fc': float(np.max(abs_log2fc)) if abs_log2fc.size else None,
        'median_abs_log2fc': float(np.median(abs_log2fc)) if abs_log2fc.size else None,
        'effective_log2fc_threshold': log2fc_threshold,
        'significance_metric': 'padj',
        'padj_threshold': float(padj_threshold),
        'multiple_testing_scope': 'per_comparison',
    }


def _build_comparison_matrix(all_deg_dfs, comparison_names):
    """Combine contrast statistics by stable feature ID, never display name.

    Gene symbols are intentionally allowed to repeat (for example aliases or
    a reference annotation with multiple features mapped to one symbol).  An
    outer merge on the user-facing ``gene`` label turns such repeats into a
    Cartesian product and can consume terabytes of memory.  ``gene_id`` is
    retained separately precisely so it remains the join key here.
    """
    if len(all_deg_dfs) != len(comparison_names):
        raise ValueError('比较结果与比较名称数量不一致，无法生成合并矩阵。')
    if not all_deg_dfs:
        return pd.DataFrame(columns=['gene_id', 'gene', 'gene_name'])

    required = {'gene_id', 'gene', 'gene_name', 'log2FC', 'padj', 'regulation'}
    matrix = None
    for comp_name, deg_df in zip(comparison_names, all_deg_dfs):
        missing = required - set(deg_df.columns)
        if missing:
            raise ValueError(
                '差异结果缺少合并比较所需列：' + '、'.join(sorted(missing))
            )
        # pyDEG deduplicates feature indices before fitting.  The defensive
        # drop here also makes a malformed third-party result unable to turn a
        # later outer merge into an unbounded Cartesian product.
        frame = deg_df.drop_duplicates(subset=['gene_id'], keep='first').copy()
        frame['gene_id'] = frame['gene_id'].astype(str)
        if matrix is None:
            matrix = frame[['gene_id', 'gene', 'gene_name']].copy()
        statistics = frame[['gene_id', 'log2FC', 'padj', 'regulation']].copy()
        statistics.columns = [
            'gene_id', f'{comp_name}_log2FC', f'{comp_name}_padj',
            f'{comp_name}_regulation',
        ]
        matrix = matrix.merge(statistics, on='gene_id', how='outer', validate='one_to_one')
    return matrix


def _run_single_comparison(adata, counts, group1_samples, group2_samples, group1, group2,
                           method, fc_threshold, pval_threshold, top_n, gene_id_to_name,
                           plots_dir, results_dir, suffix='', viz_params=None,
                           cooks_filter=True, independent_filter=True, padj_method='fdr_bh',
                           base_mean_filter=0, figure_params=None):
    """Run DEG for one comparison pair. Returns (deg_df, result_files, n_up, n_down)."""
    import omicverse as ov
    import matplotlib.pyplot as plt
    from modules.native_figures import scatter_figure

    result_files = []

    viz = viz_params or {}
    fig_width = viz.get('figure_width', 700)
    fig_height = viz.get('figure_height', 500)
    bg_color = viz.get('bg_color', 'white')
    font_family = viz.get('font_family', 'Arial')
    font_size = viz.get('font_size', 12)

    # Build count_df for OmicVerse
    count_df = pd.DataFrame(counts.T, index=adata.var_names.tolist(), columns=adata.obs.index.tolist())

    # OmicVerse pyDEG analysis
    dds = ov.bulk.pyDEG(count_df)
    n_before_dedup = count_df.shape[0]
    dds.drop_duplicates_index()
    n_after_dedup = dds.data.shape[0] if hasattr(dds, 'data') and dds.data is not None else n_before_dedup
    if n_before_dedup != n_after_dedup:
        import logging
        logging.getLogger(__name__).warning(f"[bulk_deg] 去重: {n_before_dedup} → {n_after_dedup} 基因（移除 {n_before_dedup - n_after_dedup} 个重复）")

    method_map = {
        't-test': 'ttest', 'mann-whitney': 'wilcox',
        'deseq2': 'DEseq2', 'edger': 'edgepy', 'limma': 'limma'
    }
    ov_method = method_map.get(method, 'ttest')

    if ov_method in ('edgepy', 'limma'):
        try:
            import inmoose
        except ImportError:
            raise ImportError(f"方法 {method} 需要安装 inmoose: pip install inmoose patsy")

    normalization_meta = adata.uns.get('normalization', {}) if hasattr(adata, 'uns') else {}
    is_log_transformed = bool(normalization_meta.get('is_log_transformed', False))

    # 原始计数上的简单检验沿用 pyDEG 标准化；已标准化的 log 表达矩阵不再二次标准化。
    if ov_method in ('ttest', 'wilcox') and not is_log_transformed:
        dds.normalize()

    # cooks_filter / independent_filter 仅对 DESeq2 有意义
    deg_kwargs = {'multipletests_method': padj_method}
    if ov_method == 'DEseq2':
        deg_kwargs['cooks_filter'] = cooks_filter
        deg_kwargs['independent_filter'] = independent_filter

    if ov_method == 'ttest' and is_log_transformed:
        # 连续表达值（FPKM/TPM 等）经 log2 标准化后，直接在 log 尺度上做 Welch t-test。
        # log2FC 应为两组 log2 均值之差，不能再对 log 值均值取比值。
        result = _welch_ttest_log_expression(
            dds.data, group1_samples, group2_samples, padj_method=padj_method)
    else:
        result = dds.deg_analysis(group1_samples, group2_samples, method=ov_method, **deg_kwargs)

    # Extract results
    gene_ids_list = [str(gene_id) for gene_id in result.index.tolist()]
    # ``gene`` remains the display/analysis field for backward-compatible
    # tables and enrichment.  Export the ID and supplied symbol separately so
    # a readable label never destroys the count feature identity.
    gene_name_values = [gene_id_to_name.get(gene_id, '') for gene_id in gene_ids_list]
    gene_names = [gene_name or gene_id for gene_id, gene_name in zip(
        gene_ids_list, gene_name_values)]
    log2fc = result['log2FC'].values
    pvalues = result['pvalue'].values
    # padj：优先用结果中的校正 p 值，否则自行校正
    if 'qvalue' in result.columns:
        padj = result['qvalue'].values
    elif 'padj' in result.columns:
        padj = result['padj'].values
    else:
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(np.nan_to_num(pvalues, nan=1.0), method=padj_method)
    n_genes = len(gene_names)

    # A FC entered as 2 means |log2FC| >= log2(2) = 1, not >= 2.
    log2fc_threshold = np.log2(max(float(fc_threshold), np.finfo(float).tiny))
    finite_padj = np.isfinite(padj)
    finite_log2fc = np.isfinite(log2fc)
    passes_padj = finite_padj & (padj < pval_threshold)
    passes_fc = finite_log2fc & (np.abs(log2fc) >= log2fc_threshold)
    significant = passes_padj & passes_fc
    regulation = np.full(n_genes, 'NS', dtype=object)
    regulation[significant & (log2fc >= log2fc_threshold)] = 'Up'
    regulation[significant & (log2fc <= -log2fc_threshold)] = 'Down'

    # Group means — 使用 dds 去重后的数据计算（避免重复基因均值错位）
    if not hasattr(dds, 'data') or dds.data is None:
        raise RuntimeError("pyDEG 去重后数据不可用，请检查 omicverse 版本")
    dedup_data = dds.data  # genes x samples DataFrame（去重后）
    g1_mask = dedup_data.columns.isin(group1_samples)
    g2_mask = dedup_data.columns.isin(group2_samples)
    mean1_series = dedup_data.loc[:, g1_mask].mean(axis=1)
    mean2_series = dedup_data.loc[:, g2_mask].mean(axis=1)
    mean1 = np.array([mean1_series.get(g, 0.0) for g in result.index])
    mean2 = np.array([mean2_series.get(g, 0.0) for g in result.index])

    deg_df = pd.DataFrame({
        'gene': gene_names,
        'gene_id': gene_ids_list,
        'gene_name': gene_name_values,
        'log2FC': np.round(log2fc, 4),
        'pvalue': pvalues,
        'padj': padj,
        'mean_group1': np.round(mean1, 2),
        'mean_group2': np.round(mean2, 2),
        'passes_padj': passes_padj,
        'passes_fc': passes_fc,
        'significant': significant,
        'regulation': regulation
    })
    deg_df = deg_df.sort_values('padj')

    # n_up/n_down 基于过滤前完整数据计算
    n_up = int((deg_df['regulation'] == 'Up').sum())
    n_down = int((deg_df['regulation'] == 'Down').sum())

    # 保存过滤前的完整 deg_df（用于 CSV，与火山图/MA图数据一致）
    deg_df_full = deg_df

    # 基础表达量过滤（仅影响下游箱线图，不影响 CSV 和火山图）
    if base_mean_filter > 0:
        base_mean = (deg_df['mean_group1'] + deg_df['mean_group2']) / 2
        deg_df = deg_df[base_mean >= base_mean_filter].copy()

    # Volcano plot: keep the full statistics in CSV, but use a robust display
    # ceiling so padj underflow does not flatten all points near y=0.
    file_suffix = f'_{suffix}' if suffix else ''

    # Save individual CSV — 保存过滤前的完整结果
    csv_path = os.path.join(results_dir, f'bulk_deg_results{file_suffix}.csv')
    deg_df_full.to_csv(csv_path, index=False)
    result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                         'label': f'差异表达基因列表 ({group1} vs {group2})'})

    top_genes = deg_df_full[deg_df_full['regulation'] != 'NS'].head(top_n)
    top_csv = os.path.join(results_dir, f'bulk_deg_top_genes{file_suffix}.csv')
    top_genes.to_csv(top_csv, index=False)
    result_files.append({'file_path': top_csv, 'file_type': 'csv', 'category': 'table',
                         'label': f'Top {top_n} 差异基因 ({group1} vs {group2})'})

    import matplotlib.pyplot as plt
    from figure_engine import NatureFigureDirector, export_registered_figure

    director = NatureFigureDirector()
    director_params = dict(figure_params or {'_visualization': viz})
    label_genes = director_params.get(
        'volcano_label_genes', director_params.get('plot_genes', ()))
    # Top-N governs the result table.  A compact, independent default keeps
    # the 89-mm volcano panel readable when the table exports 20+ genes.
    try:
        volcano_label_n = int(director_params.get(
            'volcano_label_n', min(max(int(top_n), 0), 4)))
    except (TypeError, ValueError):
        volcano_label_n = min(max(int(top_n), 0), 4)
    volcano_label_n = min(max(volcano_label_n, 0), 6)
    volcano_spec = director.spec_from_params(
        'volcano', director_params,
        title=_format_volcano_contrast_title(group1, group2),
        fc_threshold=float(np.log2(max(fc_threshold, np.finfo(float).tiny))),
        fdr_threshold=float(pval_threshold),
        label_n=volcano_label_n,
        label_genes=label_genes,
        show_legend=True,
    )
    # Keep the historical module contract (PNG/SVG) when no export preference
    # is supplied by a caller; an explicit static_formats/export_formats list
    # can opt into the full SVG/PDF/PNG publication bundle.
    requested_formats = (director_params.get('_visualization', {}) or {}).get(
        'static_formats', (director_params.get('_visualization', {}) or {}).get('export_formats'))
    if not requested_formats:
        volcano_spec = volcano_spec.with_updates(formats=('svg', 'png'))
    fig_vol = director.render(volcano_spec, deg_df_full)
    qa_path = os.path.join(
        results_dir, f'bulk_deg_volcano{file_suffix}_nature_readiness.json')
    exported, readiness = export_registered_figure(
        fig_vol, os.path.join(plots_dir, f'bulk_deg_volcano{file_suffix}'),
        volcano_spec, category='volcano',
        label=f'火山图 ({group1} vs {group2})', qa_path=qa_path,
    )
    result_files.extend(exported)
    if not readiness.ready:
        result_files.append({
            'file_path': '', 'file_type': 'info', 'category': 'info',
            'label': f'Volcano Nature readiness {readiness.score}/100；请查看 QA 报告。',
        })
    plt.close(fig_vol)

    # MA plot — the same scientific thresholds and palette as Volcano.
    avg_expr = (mean1 + mean2) / 2
    ma_data = deg_df_full.copy()
    ma_data['mean_expression'] = avg_expr
    ma_spec = director.spec_from_params(
        'ma', director_params,
        title=f'MA: {group1} vs {group2}',
        fc_threshold=float(np.log2(max(fc_threshold, np.finfo(float).tiny))),
        fdr_threshold=float(pval_threshold),
        label_n=volcano_label_n,
        label_genes=label_genes,
        show_legend=False,
    )
    fig_ma = director.render(ma_spec, ma_data)
    ma_exported, ma_readiness = export_registered_figure(
        fig_ma, os.path.join(plots_dir, f'bulk_deg_ma{file_suffix}'),
        ma_spec, category='ma', label=f'MA 图 ({group1} vs {group2})',
        qa_path=os.path.join(
            results_dir, f'bulk_deg_ma{file_suffix}_nature_readiness.json'),
    )
    result_files.extend(ma_exported)
    if not ma_readiness.ready:
        result_files.append({
            'file_path': '', 'file_type': 'info', 'category': 'info',
            'label': f'MA Nature readiness {ma_readiness.score}/100；请查看 QA 报告。',
        })
    plt.close(fig_ma)

    # 返回完整 deg_df 用于箱线图，过滤后 deg_df 用于下游分析
    return deg_df_full, result_files, n_up, n_down


def _run_lrt_test(adata, counts, groupby, method, pval_threshold, gene_id_to_name,
                  results_dir, padj_method='fdr_bh'):
    """Run LRT (Likelihood Ratio Test) across all groups at once. Only works with edger method."""
    import omicverse as ov
    from inmoose.edgepy import DGEList, glmLRT
    from patsy import dmatrix
    from statsmodels.stats.multitest import multipletests
    import logging
    logger = logging.getLogger(__name__)

    count_df = pd.DataFrame(counts.T, index=adata.var_names.tolist(), columns=adata.obs.index.tolist())
    dds = ov.bulk.pyDEG(count_df)
    dds.drop_duplicates_index()

    # Use the same deduplicated matrix as dds (keeps highest-sum duplicate)
    dedup_count_df = dds.data.copy()

    groups = adata.obs[groupby].astype(str)
    anno = pd.DataFrame({'group': groups.values}, index=groups.index)
    design = dmatrix("~C(group)", data=anno, return_type='dataframe')

    var = pd.DataFrame(index=dedup_count_df.index)
    var.index.name = 'gene_id'
    dge = DGEList(counts=dedup_count_df.values, samples=anno, group_col='group', genes=var)
    dge.estimateGLMCommonDisp(design=design)
    fit = dge.glmFit(design=design)
    n_coef = design.shape[1]
    lrt = glmLRT(fit, coef=list(range(1, n_coef))) if n_coef > 1 else glmLRT(fit)

    # 提取结果：兼容不同 inmoose 版本的返回格式
    try:
        if hasattr(lrt, 'table'):
            lrt_table = lrt.table
            pvalues = lrt_table['PValue'].values if 'PValue' in lrt_table.columns else lrt_table['pvalue'].values
            lrt_stat = lrt_table['F'].values if 'F' in lrt_table.columns else (lrt_table['LR'].values if 'LR' in lrt_table.columns else None)
        elif hasattr(lrt, 'PValue'):
            pvalues = np.asarray(lrt.PValue).flatten()
            lrt_stat = np.asarray(lrt.F).flatten() if hasattr(lrt, 'F') else (np.asarray(lrt.LR).flatten() if hasattr(lrt, 'LR') else None)
        elif hasattr(lrt, 'pvalue'):
            pvalues = np.asarray(lrt.pvalue).flatten()
            lrt_stat = None
        else:
            logger.warning("[bulk_deg] LRT: 无法提取 p 值")
            pvalues = np.ones(dedup_count_df.shape[0])
            lrt_stat = None
    except Exception as e:
        logger.warning(f"[bulk_deg] LRT 结果提取失败: {e}")
        pvalues = np.ones(dedup_count_df.shape[0])
        lrt_stat = None

    gene_ids = var.index.tolist()
    min_len = min(len(pvalues), len(gene_ids))
    if len(pvalues) != len(gene_ids):
        logger.warning(f"[bulk_deg] LRT 长度不匹配: pvalues={len(pvalues)}, genes={len(gene_ids)}，取交集对齐")
        pvalues = pvalues[:min_len]
        gene_ids = gene_ids[:min_len]
    if lrt_stat is not None and len(lrt_stat) != len(gene_ids):
        lrt_stat = lrt_stat[:min_len]

    _, qvalues, _, _ = multipletests(np.nan_to_num(pvalues, nan=1.0), method=padj_method)

    gene_ids = [str(gene_id) for gene_id in gene_ids]
    gene_name_values = [gene_id_to_name.get(gene_id, '') for gene_id in gene_ids]
    result_dict = {
        'gene': [gene_name or gene_id for gene_id, gene_name in zip(gene_ids, gene_name_values)],
        'gene_id': gene_ids,
        'gene_name': gene_name_values,
        'pvalue': pvalues,
        'padj': qvalues,
    }
    if lrt_stat is not None:
        result_dict['LRT_stat'] = lrt_stat
    result = pd.DataFrame(result_dict)
    result = result.sort_values('padj')
    result['significant'] = result['padj'] < pval_threshold

    lrt_csv = os.path.join(results_dir, 'bulk_deg_lrt_results.csv')
    result.to_csv(lrt_csv, index=False)

    n_sig = int(result['significant'].sum())
    return result, lrt_csv, n_sig


class BulkDEGAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_deg"
    DISPLAY_NAME = "Bulk 差异表达分析"
    DESCRIPTION = "组间差异表达基因检测：t-test / Mann-Whitney / DESeq2（基于 OmicVerse）"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        if adata.n_obs < 2:
            return "样本数不足，至少需要 2 个样本"
        if adata.n_vars == 0:
            return "基因数为 0，请检查输入数据"
        return None

    def run(self, input_path):

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)
        from modules.io_utils import infer_expression_measurement, validate_bulk_raw_counts

        requested_groupby = str(self.params.get('groupby', '') or '').strip()
        # The validated Bulk-count importer always records the design in this
        # canonical field. Keep a blank form submission usable without
        # guessing from sample-name order, while leaving any explicit groupby
        # choice authoritative.
        if not requested_groupby and 'condition' in adata.obs.columns:
            requested_groupby = 'condition'
        groupby = requested_groupby
        group1 = self.params.get('group1', '')
        group2 = self.params.get('group2', '')
        method = self.params.get('method', 't-test')
        input_measurement = infer_expression_measurement(adata, input_path)
        count_methods = {'deseq2', 'edger', 'limma'}
        preserved_counts, preserved_layer = _preserved_raw_counts(adata)
        if input_measurement != 'raw_counts' and method in count_methods and preserved_counts is None:
            raise ValueError(
                f'{method} 需要原始整数 counts；当前输入为 {input_measurement}，且未找到可用的 raw counts layer。'
                '请选择原始 counts 输入，或从保留 raw layer 的标准化输出重新分析。'
            )
        auto_log2_continuous = (
            input_measurement == 'continuous_expression'
            and method in {'t-test', 'mann-whitney'}
        )
        fc_threshold = float(self.params.get('fc_threshold', 2.0))
        pval_threshold = float(self.params.get('pval_threshold', 0.05))
        top_n = int(self.params.get('top_n', 20))
        comparisons_str = self.params.get('comparisons', '').strip()
        custom_groups_str = self.params.get('custom_groups', '').strip()
        custom_groups = _parse_custom_groups(custom_groups_str)
        comparison_pairs = _parse_comparisons(comparisons_str)

        cooks_filter = self.params.get('cooks_filter', True)
        if isinstance(cooks_filter, str):
            cooks_filter = cooks_filter.lower() in ('true', '1', 'yes', 'on')
        independent_filter = self.params.get('independent_filter', True)
        if isinstance(independent_filter, str):
            independent_filter = independent_filter.lower() in ('true', '1', 'yes', 'on')
        padj_method = self.params.get('padj_method', 'fdr_bh')
        base_mean_filter = float(self.params.get('base_mean_filter', 0))
        auto_comparisons = self.params.get('auto_comparisons', 'manual')
        reference_group = self.params.get('reference_group', '').strip()
        test_type = self.params.get('test_type', 'pairwise')
        if not np.isfinite(fc_threshold) or fc_threshold <= 0:
            raise ValueError('fc_threshold 必须是大于 0 的 Fold Change。')
        if not np.isfinite(pval_threshold) or not 0 < pval_threshold <= 1:
            raise ValueError('padj 显著性阈值必须在 (0, 1] 范围内。')
        threshold_summary = {
            # Retain the historical keys while making their meanings explicit
            # in every manifest consumed by the UI or an audit script.
            'fc_threshold': fc_threshold,
            'fc_threshold_input': fc_threshold,
            'effective_log2fc_threshold': float(np.log2(fc_threshold)),
            'pval_threshold': pval_threshold,
            'padj_threshold': pval_threshold,
            'significance_metric': 'padj',
            'padj_method': padj_method,
            'multiple_testing_scope': 'per_comparison',
        }

        self.progress(15, "构建计数矩阵...")

        # 构建基因ID→基因名映射
        gene_id_to_name = {}
        if 'gene_name' in adata.var.columns:
            for gid, gname in zip(adata.var_names, adata.var['gene_name']):
                if pd.notna(gname) and str(gname).strip():
                    gene_id_to_name[str(gid)] = str(gname).strip()

        # 构建 OmicVerse pyDEG 所需的 counts DataFrame（基因×样本）
        count_source = adata.X
        if method in count_methods and input_measurement != 'raw_counts':
            count_source = preserved_counts
            self.progress(-1, f'当前输入已标准化；{method} 使用保留的原始计数 layer “{preserved_layer}”。')
        if method in count_methods:
            validate_bulk_raw_counts(adata, matrix=count_source, context=f'{method} 差异表达')
        counts = count_source if not hasattr(count_source, 'toarray') else count_source.toarray()
        counts = counts.astype(float)
        n_inf = int(np.isinf(counts).sum())
        if n_inf > 0:
            import logging
            logging.getLogger(__name__).warning(f"[bulk_deg] 输入数据含 {n_inf} 个 inf 值（可能来自除零），已替换为 0")
        counts = np.nan_to_num(counts, nan=0.0, posinf=0.0, neginf=0.0)
        if auto_log2_continuous:
            # Continuous FPKM/TPM-like input is valid for Welch/Wilcoxon after
            # an explicit log2(x+1) transform.  Keep the transform local to
            # this DEG run and record it in the result summary instead of
            # silently mutating the uploaded matrix.
            counts = np.log2(np.clip(counts, a_min=0.0, a_max=None) + 1.0)
            normalization = dict(adata.uns.get('normalization', {}) or {})
            normalization.update({'is_log_transformed': True, 'method': 'log2_auto_for_deg'})
            adata.uns['normalization'] = normalization

        # Low-expression filtering is what keeps DESeq2/edgeR/limma from
        # spending its multiple-testing budget on genes that cannot be tested.
        # Record the upstream decision (and warn when there is none) so the DEG
        # artifact is self-contained.
        gene_filtering = _pre_deseq2_gene_filter_provenance(adata)
        if method in count_methods and not gene_filtering['applied_before_deg']:
            self.progress(-1, gene_filtering['message'])

        # 处理自动检测的分组
        if groupby == '_auto_group_':
            auto_mapping = self.params.get('_auto_group_mapping', {})
            if isinstance(auto_mapping, str):
                try:
                    auto_mapping = json.loads(auto_mapping)
                except Exception:
                    auto_mapping = {}
            if auto_mapping:
                groupby = 'auto_group'
                adata.obs[groupby] = adata.obs.index.map(lambda x: auto_mapping.get(str(x), 'unknown'))
            else:
                raise ValueError("自动分组映射数据缺失，请重新选择输入文件。")

        # 自定义合并组：在 obs 中创建临时列
        if custom_groups and groupby in adata.obs.columns:
            new_col = '_custom_group'
            adata.obs[new_col] = adata.obs[groupby].astype(str)
            for new_name, members in custom_groups.items():
                mask = adata.obs[groupby].astype(str).isin(members)
                adata.obs.loc[mask, new_col] = new_name
            groupby = new_col

        # 确定两组样本名（延迟到 single-comparison 分支处理，避免重复逻辑）。
        # Never split samples by file order when a grouping column is missing:
        # that produces plausible-looking but biologically meaningless DEG.
        from modules.io_utils import obs_grouping_info
        grouping = obs_grouping_info(
            adata, groupby, max_categories=50,
            max_numeric_categories=20, require_multiple=True,
        )
        if not grouping['valid'] and groupby not in adata.obs.columns:
            # A table upload has no obs metadata, but sample names often encode
            # a real condition (Ctrl_1_count / Treat_1_count).  Infer only when
            # every sample follows that explicit replicate pattern; never split
            # an arbitrary file ordering into synthetic groups.
            from modules.io_utils import infer_sample_group_candidates
            candidates = infer_sample_group_candidates(adata.obs.index.tolist())
            if candidates:
                groupby = 'auto_group'
                adata.obs[groupby] = adata.obs.index.map(candidates[0]['mapping'])
                grouping = obs_grouping_info(
                    adata, groupby, max_categories=50,
                    max_numeric_categories=20, require_multiple=True,
                )
                self.progress(-1, f"未找到 '{requested_groupby}'，已从样本名推断分组列 auto_group")
        if not grouping['valid']:
            raise ValueError(
                f"分组列 '{groupby}' 不是有效的分类分组列：{grouping['reason']}"
            )

        # Comparisons are labels from the currently selected grouping column;
        # they are not reusable after switching to a factorized column.  Fail
        # before launching any statistic rather than reporting the misleading
        # generic "all comparisons had too few samples" message.
        if comparison_pairs and auto_comparisons == 'manual':
            available_groups = set(adata.obs[groupby].astype(str).unique().tolist())
            requested_groups = {
                group for pair in comparison_pairs for group in pair
            }
            unknown_groups = sorted(requested_groups - available_groups)
            if unknown_groups:
                available_text = '、'.join(sorted(available_groups)[:12])
                unknown_text = '、'.join(unknown_groups[:12])
                raise ValueError(
                    f"指定比较与分组列 '{groupby}' 不一致：未找到 {unknown_text}。"
                    f"该列可用分组为 {available_text}。"
                    '请改回与比较名称对应的分组列，或更新“多组比较”中的组名。'
                )

        self.progress(25, "准备差异分析...")

        plots_dir = os.path.join(self.project_dir, 'plots')
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        # 自动生成比较列表
        if auto_comparisons in ('all_pairwise', 'vs_reference'):
            all_groups = sorted(adata.obs[groupby].astype(str).unique().tolist())
            # 注意：sorted() 对时间序列分组（如 Treated_1h, Treated_10h）可能排序不符预期
            if auto_comparisons == 'vs_reference':
                if not reference_group:
                    reference_group = all_groups[0]
                comparison_pairs = [(g, reference_group) for g in all_groups if g != reference_group]
            elif auto_comparisons == 'all_pairwise':
                comparison_pairs = []
                for i in range(len(all_groups)):
                    for j in range(i + 1, len(all_groups)):
                        g1, g2 = all_groups[i], all_groups[j]
                        # 确保参考组在第二个位置（group2）以统一 logFC 方向
                        if reference_group and g2 == reference_group:
                            comparison_pairs.append((g1, g2))
                        elif reference_group and g1 == reference_group:
                            comparison_pairs.append((g2, g1))
                        else:
                            comparison_pairs.append((g1, g2))

        # LRT 多组检验（仅 edger）
        lrt_files = []
        lrt_n_sig = 0
        if test_type == 'lrt':
            if method != 'edger':
                self.progress(27, f"LRT 仅支持 edger 方法，当前方法 {method} 将使用 pairwise 检验")
                test_type = 'pairwise'
            else:
                try:
                    import inmoose
                    self.progress(27, "运行 LRT 多组检验...")
                    lrt_result, lrt_csv, lrt_n_sig = _run_lrt_test(
                        adata, counts, groupby, method, pval_threshold,
                        gene_id_to_name, results_dir, padj_method=padj_method
                    )
                    lrt_files.append({'file_path': lrt_csv, 'file_type': 'csv', 'category': 'table', 'label': f'LRT 多组检验结果 ({lrt_n_sig} 个显著基因)'})
                except ImportError:
                    self.progress(27, "inmoose 未安装，跳过 LRT 检验")
                    test_type = 'pairwise'
                except Exception as e:
                    self.progress(27, f"LRT 检验失败: {e}，回退到 pairwise")
                    test_type = 'pairwise'

        if comparison_pairs:
            # 多组比较模式
            all_deg_dfs = []
            valid_comparisons = []
            skipped = []
            for idx, (g1, g2) in enumerate(comparison_pairs):
                pct = 25 + int(55 * idx / len(comparison_pairs))
                self.progress(pct, f"比较 {g1} vs {g2} ({idx+1}/{len(comparison_pairs)})...")

                g1_samples = list(adata.obs.index[adata.obs[groupby] == g1])
                g2_samples = list(adata.obs.index[adata.obs[groupby] == g2])
                if len(g1_samples) < 2 or len(g2_samples) < 2:
                    skipped.append(f'{g1}-vs-{g2}（样本数不足）')
                    continue

                deg_df, files, _, _ = _run_single_comparison(
                    adata, counts, g1_samples, g2_samples, g1, g2,
                    method, fc_threshold, pval_threshold, top_n,
                    gene_id_to_name, plots_dir, results_dir, suffix=str(idx),
                    viz_params=self.params.get('_visualization', {}),
                    cooks_filter=cooks_filter, independent_filter=independent_filter,
                    padj_method=padj_method, base_mean_filter=base_mean_filter,
                    figure_params=self.params)
                deg_df['comparison'] = f'{g1}-vs-{g2}'
                all_deg_dfs.append(deg_df)
                valid_comparisons.append((g1, g2))
                result_files.extend(files)

            if skipped:
                result_files.append({'file_path': '', 'file_type': 'info', 'category': 'info',
                                     'label': f'跳过的比较: {"; ".join(skipped)}'})

            if all_deg_dfs:
                merged_df = pd.concat(all_deg_dfs, ignore_index=True)
                merged_csv = os.path.join(results_dir, 'bulk_deg_all_comparisons.csv')
                merged_df.to_csv(merged_csv, index=False)
                result_files.append({'file_path': merged_csv, 'file_type': 'csv', 'category': 'table',
                                     'label': '所有比较合并结果'})

                # 多比较合并表（logFC + padj 矩阵格式）
                if len(all_deg_dfs) >= 2:
                    # 构建 logFC + padj 矩阵
                    comparison_names = [f'{g1}-vs-{g2}' for g1, g2 in valid_comparisons]
                    # Use feature IDs, not display symbols, as a matrix key.
                    # Symbols can repeat and an outer merge on them produces a
                    # Cartesian product large enough to exhaust system memory.
                    merged_matrix_df = _build_comparison_matrix(
                        all_deg_dfs, comparison_names,
                    )
                    merged_matrix_csv = os.path.join(results_dir, 'bulk_deg_merged_comparisons.csv')
                    merged_matrix_df.to_csv(merged_matrix_csv, index=False)
                    result_files.append({'file_path': merged_matrix_csv, 'file_type': 'csv', 'category': 'table', 'label': '多比较合并结果'})

                    # 共享差异基因统计
                    up_sets = [set(df[df['regulation'] == 'Up']['gene']) for df in all_deg_dfs if 'Up' in df['regulation'].values]
                    down_sets = [set(df[df['regulation'] == 'Down']['gene']) for df in all_deg_dfs if 'Down' in df['regulation'].values]
                    shared_up = len(set.intersection(*up_sets)) if len(up_sets) >= 2 else 0
                    shared_down = len(set.intersection(*down_sets)) if len(down_sets) >= 2 else 0

                    # Multi-comparison Volcano uses the Composer so each panel
                    # is the same fixed NatureVolcano template and the legend is
                    # shared instead of repeated inside every small axis.
                    import matplotlib.pyplot as plt
                    from figure_engine import NatureFigureComposer, NatureFigureDirector, export_registered_figure
                    n_comp = len(comparison_names)
                    ncols_multi = min(3, n_comp)
                    nrows_multi = int(np.ceil(n_comp / ncols_multi))
                    director_multi = NatureFigureDirector()
                    outer_spec = director_multi.spec_from_params(
                        'diagnostic', self.params, width='double', title='Multi-comparison Volcano',
                    ).with_updates(plot_type='composite', formats=('svg', 'pdf', 'png'),
                                   height_mm=82.0 if nrows_multi == 1 else 145.0)
                    panels = []
                    for m_idx, (comp_name, deg_df) in enumerate(zip(comparison_names, all_deg_dfs)):
                        panel_spec = director_multi.spec_from_params(
                            'volcano', self.params, width='single',
                            title=_format_volcano_contrast_title(g1, g2),
                            fc_threshold=np.log2(max(fc_threshold, 1e-12)),
                            fdr_threshold=pval_threshold, label_n=0, show_legend=True,
                            extra={'legend_counts': False},
                        )
                        from figure_engine.composer import FigurePanel
                        panels.append(FigurePanel(
                            plot_type='volcano', data=deg_df, spec=panel_spec,
                            label=chr(ord('a') + m_idx), row=m_idx // ncols_multi,
                            column=m_idx % ncols_multi,
                        ))
                    fig_multi = NatureFigureComposer(director_multi).compose(
                        panels, outer_spec, nrows=nrows_multi, ncols=ncols_multi,
                        shared_legend=True,
                    )
                    export_multi, report_multi = export_registered_figure(
                        fig_multi, os.path.join(plots_dir, 'bulk_deg_volcano_multi'), outer_spec,
                        category='volcano', label='多组比较 Volcano',
                        qa_path=os.path.join(results_dir, 'bulk_deg_volcano_multi_nature_readiness.json'),
                    )
                    result_files.extend(export_multi)
                    if not report_multi.ready:
                        self.progress(-1, f'多组比较 Volcano Nature readiness {report_multi.score}/100；请查看 QA 报告。')
                    plt.close(fig_multi)
                else:
                    shared_up = 0
                    shared_down = 0
                    comparison_names = [f'{g1}-vs-{g2}' for g1, g2 in comparison_pairs]

                # 箱线图（取第一个有效比较的结果生成）
                first_deg = all_deg_dfs[0]
                g1_first, g2_first = valid_comparisons[0]
                g1_samples = list(adata.obs.index[adata.obs[groupby] == g1_first])
                g2_samples = list(adata.obs.index[adata.obs[groupby] == g2_first])
                self._draw_boxplots(first_deg, g1_first, g2_first, g1_samples, g2_samples,
                                    gene_id_to_name, counts, adata, plots_dir, result_files)
            else:
                raise ValueError("所有比较均因样本数不足被跳过，请检查分组信息。")

            # LRT 文件加入结果
            result_files.extend(lrt_files)

            # 构建 per_comparison 统计
            per_comparison = {}
            for comp_name, deg_df in zip(
                [f'{g1}-vs-{g2}' for g1, g2 in valid_comparisons],
                all_deg_dfs):
                per_comparison[comp_name] = {
                    'n_up': int((deg_df['regulation'] == 'Up').sum()),
                    'n_down': int((deg_df['regulation'] == 'Down').sum()),
                    'filter_diagnostics': _comparison_filter_diagnostics(
                        deg_df, fc_threshold, pval_threshold,
                    ),
                }

            summary = {
                'method': method,
                'test_type': test_type,
                'n_comparisons': len(all_deg_dfs),
                'comparisons': [f'{g1}-vs-{g2}' for g1, g2 in valid_comparisons],
                'per_comparison': per_comparison,
                'shared_up_genes': shared_up,
                'shared_down_genes': shared_down,
                'skipped_comparisons': skipped,
                **threshold_summary,
                'lrt_n_sig': lrt_n_sig,
                'input_measurement': input_measurement,
                'auto_transform': 'log2(x+1)' if auto_log2_continuous else None,
                'requested_groupby': requested_groupby,
                'groupby': groupby,
                'gene_filtering': gene_filtering,
            }
        else:
            # 单次比较模式
            if groupby in adata.obs.columns:
                groups = adata.obs[groupby].unique().tolist()
                if not group1 or group1 not in groups:
                    group1 = groups[0]
                group2_label = group2
                if group2 == 'rest' or (not group2 or group2 not in groups):
                    group2_samples = [s for s in adata.obs.index if adata.obs.loc[s, groupby] != group1]
                    group2_label = f'rest (n={len(group2_samples)})'
                else:
                    group2_samples = list(adata.obs.index[adata.obs[groupby] == group2])
                group1_samples = list(adata.obs.index[adata.obs[groupby] == group1])
            else:
                raise ValueError(
                    f"未找到分组列 '{groupby}'；不能按样本文件顺序自动拆成 Group1/Group2。"
                )

            if len(group1_samples) < 2 or len(group2_samples) < 2:
                raise ValueError(f"样本数不足：{group1}={len(group1_samples)}个, {group2_label}={len(group2_samples)}个。每组至少需要2个样本。")

            self.progress(30, f"差异分析: {group1} vs {group2_label}...")
            deg_df, files, n_up, n_down = _run_single_comparison(
                adata, counts, group1_samples, group2_samples, group1, group2_label,
                method, fc_threshold, pval_threshold, top_n,
                gene_id_to_name, plots_dir, results_dir,
                viz_params=self.params.get('_visualization', {}),
                cooks_filter=cooks_filter, independent_filter=independent_filter,
                padj_method=padj_method, base_mean_filter=base_mean_filter,
                figure_params=self.params)
            result_files.extend(files)

            # 箱线图
            self._draw_boxplots(deg_df, group1, group2_label, group1_samples, group2_samples,
                                gene_id_to_name, counts, adata, plots_dir, result_files)

            summary = {
                'comparison': f'{group1} vs {group2_label}',
                'method': method,
                'n_genes_total': len(deg_df),
                'n_up': n_up,
                'n_down': n_down,
                'filter_diagnostics': _comparison_filter_diagnostics(
                    deg_df, fc_threshold, pval_threshold,
                ),
                **threshold_summary,
                'lrt_n_sig': lrt_n_sig,
                'input_measurement': input_measurement,
                'auto_transform': 'log2(x+1)' if auto_log2_continuous else None,
                'requested_groupby': requested_groupby,
                'groupby': groupby,
                'gene_filtering': gene_filtering,
            }

        # LRT 文件：多比较模式已在 line 551 添加，单次比较模式在此添加
        if not comparison_pairs:
            result_files.extend(lrt_files)

        self.progress(95, "保存 h5ad...")
        # 临时标记分组（不覆写已有 group 列，避免错误输入污染输出 h5ad）
        if groupby in adata.obs.columns and 'group' not in adata.obs.columns:
            adata.obs['group'] = adata.obs[groupby]
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_deg_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }

    def _draw_boxplots(self, deg_df, group1, group2, group1_samples, group2_samples,
                       gene_id_to_name, counts, adata, plots_dir, result_files):
        import matplotlib.pyplot as plt
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT, NATURE_GRID
        from figure_engine import NatureFigureDirector, export_registered_figure

        director = NatureFigureDirector()
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        nature_formats = ('svg', 'pdf', 'png')

        boxplot_n = min(int(self.params.get('boxplot_n', 5)), 20)
        # 预构建索引映射（O(1) 查找）
        obs_idx = {s: i for i, s in enumerate(adata.obs.index)}
        var_idx = {g: i for i, g in enumerate(adata.var_names)}
        # 构建双向映射：gene_name → gene_id（处理重名取首个）
        name_to_id = {}
        if gene_id_to_name:
            for gid, gname in gene_id_to_name.items():
                if gname not in name_to_id:
                    name_to_id[gname] = gid
        if 'gene_name' in adata.var.columns:
            for gid, gname in zip(adata.var_names, adata.var['gene_name']):
                gname_str = str(gname).strip()
                if gname_str and gname_str not in name_to_id:
                    name_to_id[gname_str] = str(gid)

        plot_gene_list = []
        plot_genes_str = self.params.get('plot_genes', '').strip()
        if plot_genes_str:
            plot_gene_list = [g.strip() for g in plot_genes_str.split(',') if g.strip()]

        top_de_genes = deg_df[deg_df['regulation'] != 'NS'].head(boxplot_n)
        for g in top_de_genes['gene'].tolist():
            if g not in plot_gene_list:
                plot_gene_list.append(g)

        if plot_gene_list:
            self.progress(82, f"生成 {len(plot_gene_list)} 个基因箱线图...")
        unmatched_genes = []
        for pg in plot_gene_list:
            pg_id = pg if pg in adata.var_names else name_to_id.get(pg, pg)
            if pg_id not in var_idx:
                unmatched_genes.append(pg)
                continue
            groups_values = []
            for grp_name, samples in [(group1, group1_samples), (group2, group2_samples)]:
                sample_indices = [obs_idx[s] for s in samples if s in obs_idx]
                groups_values.append(counts[sample_indices, var_idx[pg_id]])
            box_spec = director.spec_from_params(
                'diagnostic', self.params, width='single', title=f'{pg} expression',
            ).with_updates(extra={'kind': 'boxplot'}, formats=nature_formats, height_mm=68.0)
            fig_box = director.render(box_spec, {
                'kind': 'boxplot', 'groups': [str(group1), str(group2)],
                'values': groups_values, 'ylabel': 'Expression',
            })
            exported, report = export_registered_figure(
                fig_box, os.path.join(plots_dir, f'bulk_deg_box_{pg}'), box_spec,
                category='boxplot', label=f'{pg} Boxplot',
                qa_path=os.path.join(results_dir, f'bulk_deg_box_{pg}_nature_readiness.json'),
            )
            result_files.extend(exported)
            if not report.ready:
                self.progress(-1, f'{pg} 箱线图 Nature readiness {report.score}/100；请查看 QA 报告。')
            plt.close(fig_box)
        if unmatched_genes:
            import logging
            logging.getLogger(__name__).warning(f"[bulk_deg] 以下基因未找到，已跳过箱线图: {', '.join(unmatched_genes)}")
