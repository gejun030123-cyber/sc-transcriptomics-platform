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
    'Up': '#B64342',
    'Down': '#0F4D92',
    'NS': '#C7CDD6',
}


def _volcano_y_limit(padj, pval_threshold=0.05, max_display=30.0):
    """Choose a readable display ceiling without changing stored statistics."""
    padj = np.asarray(padj, dtype=float)
    safe_padj = np.clip(np.nan_to_num(padj, nan=1.0, posinf=1.0, neginf=1.0),
                        np.finfo(float).tiny, 1.0)
    y = -np.log10(safe_padj)
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return 8.0
    robust_ceiling = float(np.nanpercentile(finite, 99.5))
    threshold_y = -np.log10(max(float(pval_threshold), np.finfo(float).tiny))
    # Most publication volcano plots are more legible around 0–20.  A hard
    # ceiling prevents padj underflow (padj=0 -> ~307) from flattening the
    # entire plot; clipped points are explicitly marked by _draw_bulk_volcano.
    return float(min(max_display, max(8.0, threshold_y + 3.0,
                                    np.ceil(robust_ceiling + 1.0))))


def _volcano_label_positions(label_rows, y_limit, x_limit):
    """Place volcano labels on two non-overlapping, in-panel rails.

    Labels are allocated independently on the up- and down-regulated sides,
    then separated in data coordinates.  This keeps the result deterministic
    in PNG/SVG exports and avoids an optional layout dependency moving text
    outside the plotting area.
    """
    if not label_rows:
        return []

    min_y = float(y_limit) * 0.08
    max_y = float(y_limit) * 0.94
    min_gap = max(0.55, float(y_limit) * 0.04)
    positions = []

    for direction, side in (('Up', 1), ('Down', -1)):
        side_rows = [row for row in label_rows if row['regulation'] == direction]
        # Start with the point's y coordinate, then stack neighbouring labels
        # vertically.  Sorting low-to-high makes the final order stable.
        side_rows.sort(key=lambda row: row['y_point'])
        if not side_rows:
            continue

        y_text = []
        previous = min_y - min_gap
        for row in side_rows:
            proposed = max(min_y, min(max_y, row['y_point']))
            y_position = max(proposed, previous + min_gap)
            y_text.append(y_position)
            previous = y_position

        # Shift a full stack down together when labels near the top would
        # exceed the rail.  With the <=10 annotation cap, this preserves the
        # minimum gap in normal use; the compression fallback is defensive.
        overflow = y_text[-1] - max_y
        if overflow > 0:
            y_text = [y - overflow for y in y_text]
        if y_text[0] < min_y:
            y_text = np.linspace(min_y, max_y, len(y_text)).tolist()

        rail_x = side * float(x_limit) * 0.78
        for row, text_y in zip(side_rows, y_text):
            positions.append({**row, 'x_text': rail_x, 'y_text': float(text_y)})

    return positions


def _draw_bulk_volcano(ax, deg_df, pval_threshold=0.05, fc_threshold=2.0,
                       title='', top_n=10, show_legend=True, y_limit=None,
                       x_limit=None, colors=None):
    """Draw a restrained, readable Bulk RNA-seq volcano panel."""
    from matplotlib.lines import Line2D
    from modules.figure_style import NATURE_AXIS, NATURE_GRID, NATURE_TEXT

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
    clipped = np.isfinite(y_raw) & (y_raw > y_limit)
    y_display = np.minimum(y_raw, y_limit * 0.985)
    log2fc_threshold = np.log2(max(float(fc_threshold), np.finfo(float).tiny))

    finite_x = np.abs(log2fc[np.isfinite(log2fc)])
    if x_limit is None:
        max_abs_x = float(np.nanmax(finite_x)) if finite_x.size else log2fc_threshold * 2
        x_limit = max(2.5, log2fc_threshold * 1.8, max_abs_x * 1.08)
        x_limit = float(np.ceil(x_limit * 2.0) / 2.0)

    masks = {
        'NS': regulation == 'NS',
        'Down': regulation == 'Down',
        'Up': regulation == 'Up',
    }
    # Draw the dense neutral cloud first; significant points remain visible.
    for label, size, alpha, zorder in (
        ('NS', 8, 0.30, 1), ('Down', 11, 0.82, 2), ('Up', 11, 0.84, 2),
    ):
        mask = masks[label] & np.isfinite(log2fc) & np.isfinite(y_display)
        if mask.any():
            ax.scatter(log2fc[mask], y_display[mask], s=size,
                       color=palette[label], alpha=alpha,
                       linewidths=0, rasterized=True, zorder=zorder)

    # Triangles make the display-only y clipping explicit rather than silently
    # hiding extremely small adjusted P values.
    if clipped.any():
        ax.scatter(log2fc[clipped], y_display[clipped], s=17,
                   c=[palette.get(label, palette['NS'])
                      for label in regulation[clipped]],
                   marker='^', alpha=0.92, linewidths=0,
                   rasterized=True, zorder=4)

    ax.axhline(-np.log10(max(float(pval_threshold), np.finfo(float).tiny)),
               color='#98A2B3', linestyle='--', linewidth=0.75, zorder=3)
    ax.axvline(log2fc_threshold, color='#98A2B3', linestyle='--', linewidth=0.75, zorder=3)
    ax.axvline(-log2fc_threshold, color='#98A2B3', linestyle='--', linewidth=0.75, zorder=3)
    ax.set_xlim(-x_limit, x_limit)
    ax.set_ylim(0, y_limit)
    ax.set_title(title, loc='left', pad=8, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    ax.set_xlabel('log2(Fold Change)', fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel('-log10(adjusted P value)', fontsize=9, color=NATURE_TEXT)
    ax.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.65)
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
                   label=f"Up-regulated (n={counts['Up']})"),
            Line2D([], [], marker='o', linestyle='None', markersize=5,
                   markerfacecolor=palette['Down'], markeredgewidth=0,
                   label=f"Down-regulated (n={counts['Down']})"),
            Line2D([], [], marker='o', linestyle='None', markersize=5,
                   markerfacecolor=palette['NS'], markeredgewidth=0,
                   label=f"Not significant (n={counts['NS']})"),
            Line2D([], [], color='#98A2B3', linestyle='--', linewidth=0.75,
                   label=f"Criteria: FDR<{pval_threshold:.3g}, FC>{fc_threshold:.3g}"),
        ]
        if clipped.any():
            handles.append(Line2D([], [], marker='^', linestyle='None', markersize=5,
                                   markerfacecolor='#667085', markeredgewidth=0,
                                   label=f'{int(clipped.sum())} points clipped at y={y_limit:g}'))
        ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                  frameon=False, fontsize=7.5, handlelength=1.3,
                  borderaxespad=0.0)

    if top_n > 0 and not deg_df.empty:
        sig = deg_df[deg_df['regulation'].isin(['Up', 'Down'])].copy()
        # ``top_n`` also controls the Top-DEG table.  A table can comfortably
        # contain 20 entries, whereas annotating 20 points makes a volcano plot
        # unreadable.  Keep the plot concise without changing any output table.
        annotation_limit = min(max(int(top_n), 0), 10)
        half = max(1, int(np.ceil(annotation_limit / 2)))
        selected = pd.concat([
            sig[sig['regulation'] == 'Up'].nsmallest(half, 'padj'),
            sig[sig['regulation'] == 'Down'].nsmallest(half, 'padj'),
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
                'x_point': x_value,
                'y_point': min(y_value, y_limit * 0.985),
            })

        for label in _volcano_label_positions(label_rows, y_limit, x_limit):
            ax.annotate(
                label['gene'], (label['x_point'], label['y_point']),
                xytext=(label['x_text'], label['y_text']), textcoords='data',
                fontsize=7, ha='left' if label['regulation'] == 'Up' else 'right',
                va='center', color=NATURE_TEXT, zorder=5,
                bbox={
                    'boxstyle': 'round,pad=0.22,rounding_size=0.12',
                    'facecolor': 'white', 'edgecolor': '#D0D5DD',
                    'linewidth': 0.55, 'alpha': 0.96,
                },
                arrowprops={
                    'arrowstyle': '->', 'color': '#667085', 'linewidth': 0.5,
                    'shrinkA': 2, 'shrinkB': 2,
                },
            )

    return {'y_limit': float(y_limit), 'x_limit': float(x_limit),
            'n_clipped': int(clipped.sum())}


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


def _run_single_comparison(adata, counts, group1_samples, group2_samples, group1, group2,
                           method, fc_threshold, pval_threshold, top_n, gene_id_to_name,
                           plots_dir, results_dir, suffix='', viz_params=None,
                           cooks_filter=True, independent_filter=True, padj_method='fdr_bh',
                           base_mean_filter=0):
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
    gene_ids_list = result.index.tolist()
    if gene_id_to_name:
        gene_names = [gene_id_to_name.get(g, g) for g in gene_ids_list]
    else:
        gene_names = gene_ids_list
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

    # Regulation direction
    log2fc_threshold = np.log2(fc_threshold)
    regulation = []
    for i in range(n_genes):
        if padj[i] < pval_threshold and log2fc[i] >= log2fc_threshold:
            regulation.append('Up')
        elif padj[i] < pval_threshold and log2fc[i] <= -log2fc_threshold:
            regulation.append('Down')
        else:
            regulation.append('NS')

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
        'log2FC': np.round(log2fc, 4),
        'pvalue': pvalues,
        'padj': padj,
        'mean_group1': np.round(mean1, 2),
        'mean_group2': np.round(mean2, 2),
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
    fig_vol, ax_vol = plt.subplots(figsize=(7.5, 5.0), dpi=150)
    _draw_bulk_volcano(
        ax_vol, deg_df_full, pval_threshold=pval_threshold,
        fc_threshold=fc_threshold, title=f'火山图 ({group1} vs {group2})',
        top_n=top_n, show_legend=True,
    )
    fig_vol.tight_layout(pad=1.1)
    result_files.extend(_save_native_figure(
        fig_vol, plots_dir, f'bulk_deg_volcano{file_suffix}', 'volcano',
        f'火山图 ({group1} vs {group2})', viz,
    ))

    # MA plot
    colors = np.asarray([
        VOLCANO_COLORS.get(reg, VOLCANO_COLORS['NS'])
        for reg in regulation
    ], dtype=object)
    avg_expr = (mean1 + mean2) / 2
    fig_ma = scatter_figure(
        np.log2(avg_expr + 1), log2fc,
        title=f'MA 图 ({group1} vs {group2})',
        x_label='log2(Average Expression)', y_label='log2(Fold Change)',
        colors=colors, size=10, alpha=0.72,
    )
    fig_ma.axes[0].axhline(0, color='#98A2B3', linewidth=0.7)
    result_files.extend(_save_native_figure(
        fig_ma, plots_dir, f'bulk_deg_ma{file_suffix}', 'ma',
        f'MA 图 ({group1} vs {group2})', viz,
    ))

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

    result_dict = {
        'gene': [gene_id_to_name.get(g, g) if gene_id_to_name else g for g in gene_ids],
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
        from modules.io_utils import infer_expression_measurement

        requested_groupby = str(self.params.get('groupby', 'condition') or '').strip()
        groupby = requested_groupby
        group1 = self.params.get('group1', '')
        group2 = self.params.get('group2', '')
        method = self.params.get('method', 't-test')
        input_measurement = infer_expression_measurement(adata, input_path)
        if input_measurement != 'raw_counts' and method in {'deseq2', 'edger', 'limma'}:
            raise ValueError(f'{method} 需要原始整数 counts；当前输入为 {input_measurement}。请选择 t-test/Mann-Whitney，或从原始 counts 重新分析。')
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

        self.progress(15, "构建计数矩阵...")

        # 构建基因ID→基因名映射
        gene_id_to_name = {}
        if 'gene_name' in adata.var.columns:
            for gid, gname in zip(adata.var_names, adata.var['gene_name']):
                if pd.notna(gname) and str(gname).strip():
                    gene_id_to_name[str(gid)] = str(gname).strip()

        # 构建 OmicVerse pyDEG 所需的 counts DataFrame（基因×样本）
        counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
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
                    padj_method=padj_method, base_mean_filter=base_mean_filter)
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
                    # 构建 logFC + padj 矩阵（按 gene 列外连接，避免索引错位）
                    merged_matrix_df = all_deg_dfs[0][['gene']].copy()
                    for comp_name, deg_df in zip(comparison_names, all_deg_dfs):
                        sub = deg_df[['gene', 'log2FC', 'padj', 'regulation']].copy()
                        sub.columns = ['gene', f'{comp_name}_log2FC', f'{comp_name}_padj', f'{comp_name}_regulation']
                        merged_matrix_df = merged_matrix_df.merge(sub, on='gene', how='outer')
                    merged_matrix_csv = os.path.join(results_dir, 'bulk_deg_merged_comparisons.csv')
                    merged_matrix_df.to_csv(merged_matrix_csv, index=False)
                    result_files.append({'file_path': merged_matrix_csv, 'file_type': 'csv', 'category': 'table', 'label': '多比较合并结果'})

                    # 共享差异基因统计
                    up_sets = [set(df[df['regulation'] == 'Up']['gene']) for df in all_deg_dfs if 'Up' in df['regulation'].values]
                    down_sets = [set(df[df['regulation'] == 'Down']['gene']) for df in all_deg_dfs if 'Down' in df['regulation'].values]
                    shared_up = len(set.intersection(*up_sets)) if len(up_sets) >= 2 else 0
                    shared_down = len(set.intersection(*down_sets)) if len(down_sets) >= 2 else 0

                    # Volcano 并排展示（静态 Matplotlib 面板）
                    import matplotlib.pyplot as plt
                    from modules.figure_style import NATURE_TEXT
                    n_comp = len(comparison_names)
                    fig_multi, axes_multi = plt.subplots(
                        1, n_comp, figsize=(max(7.5, 4.1 * n_comp), 4.8),
                        dpi=150, squeeze=False,
                    )
                    all_padj = np.concatenate([
                        pd.to_numeric(df['padj'], errors='coerce').to_numpy(dtype=float)
                        for df in all_deg_dfs
                    ])
                    all_log2fc = np.concatenate([
                        pd.to_numeric(df['log2FC'], errors='coerce').to_numpy(dtype=float)
                        for df in all_deg_dfs
                    ])
                    multi_y_limit = _volcano_y_limit(all_padj, pval_threshold)
                    finite_multi_x = np.abs(all_log2fc[np.isfinite(all_log2fc)])
                    multi_x_limit = max(
                        2.5, np.log2(max(fc_threshold, np.finfo(float).tiny)) * 1.8,
                        float(np.nanmax(finite_multi_x)) * 1.08 if finite_multi_x.size else 2.5,
                    )
                    multi_x_limit = float(np.ceil(multi_x_limit * 2.0) / 2.0)
                    for m_idx, (comp_name, deg_df) in enumerate(zip(comparison_names, all_deg_dfs)):
                        ax_multi = axes_multi[0, m_idx]
                        _draw_bulk_volcano(
                            ax_multi, deg_df, pval_threshold=pval_threshold,
                            fc_threshold=fc_threshold, title=comp_name, top_n=0,
                            show_legend=(m_idx == 0), y_limit=multi_y_limit,
                            x_limit=multi_x_limit,
                        )
                        if m_idx > 0:
                            ax_multi.set_ylabel('')
                    fig_multi.suptitle('多组比较 Volcano 图', x=0.05, ha='left',
                                       fontsize=11, fontweight='semibold', color=NATURE_TEXT)
                    result_files.extend(_save_native_figure(
                        fig_multi, plots_dir, 'bulk_deg_volcano_multi', 'volcano',
                        '多组比较 Volcano', self.params.get('_visualization', {}),
                    ))
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
                'fc_threshold': fc_threshold,
                'pval_threshold': pval_threshold,
                'lrt_n_sig': lrt_n_sig,
                'input_measurement': input_measurement,
                'auto_transform': 'log2(x+1)' if auto_log2_continuous else None,
                'requested_groupby': requested_groupby,
                'groupby': groupby,
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
                padj_method=padj_method, base_mean_filter=base_mean_filter)
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
                'fc_threshold': fc_threshold,
                'pval_threshold': pval_threshold,
                'lrt_n_sig': lrt_n_sig,
                'input_measurement': input_measurement,
                'auto_transform': 'log2(x+1)' if auto_log2_continuous else None,
                'requested_groupby': requested_groupby,
                'groupby': groupby,
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
            fig_box, ax_box = plt.subplots(figsize=(5.5, 4.5), dpi=150)
            boxes = ax_box.boxplot(groups_values, tick_labels=[str(group1), str(group2)],
                                   patch_artist=True, showfliers=False)
            for patch, color in zip(boxes['boxes'], [NATURE_PALETTE[0], NATURE_PALETTE[3]]):
                patch.set_facecolor(color)
                patch.set_alpha(0.62)
                patch.set_edgecolor(color)
            for idx, vals in enumerate(groups_values, 1):
                jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) else []
                ax_box.scatter(np.full(len(vals), idx) + jitter, vals, s=14,
                               color='#667085', alpha=0.55, linewidths=0)
            ax_box.set_title(f'{pg} 表达', loc='left', fontsize=10,
                             fontweight='semibold', color=NATURE_TEXT)
            ax_box.set_ylabel('Expression', fontsize=9)
            ax_box.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
            result_files.extend(_save_native_figure(
                fig_box, plots_dir, f'bulk_deg_box_{pg}', 'boxplot',
                f'{pg} Boxplot', self.params.get('_visualization', {}),
            ))
        if unmatched_genes:
            import logging
            logging.getLogger(__name__).warning(f"[bulk_deg] 以下基因未找到，已跳过箱线图: {', '.join(unmatched_genes)}")
