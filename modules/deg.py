from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping
from modules.sc_de_utils import log1p_adata_for_cell_level_de, normalization_semantics


def _benjamini_yekutieli(pvalues):
    """Adjust p values with the dependency-safe Benjamini–Yekutieli rule."""
    import numpy as np

    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted
    indices = np.where(valid)[0]
    ordered = indices[np.argsort(values[indices])]
    n_valid = len(ordered)
    harmonic = float(np.sum(1.0 / np.arange(1, n_valid + 1)))
    running = 1.0
    for rank in range(n_valid, 0, -1):
        index = ordered[rank - 1]
        running = min(running, values[index] * n_valid * harmonic / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def _apply_by_adjustment(rank_result):
    """Replace Scanpy's BH adjusted fields with BY for every tested group."""
    pvalues = rank_result.get('pvals') if isinstance(rank_result, dict) else None
    adjusted = rank_result.get('pvals_adj') if isinstance(rank_result, dict) else None
    if pvalues is None or adjusted is None:
        return False
    groups = getattr(pvalues.dtype, 'names', None) or ()
    for group in groups:
        adjusted[group] = _benjamini_yekutieli(pvalues[group])
    return bool(groups)


VOLCANO_COLOR_MAP = {
    'Up': '#B64342',       # Nature 风格红色：上调
    'Down': '#0F4D92',     # Nature 风格深蓝：下调
    'NS': '#98A2B3',       # 中性灰：不显著
}


def classify_volcano_regulation(logfc, padj, logfc_cutoff=1.0, pval_cutoff=0.05):
    """Classify genes for the three-colour volcano plot.

    The same adjusted-p-value and absolute-logFC thresholds are used by the
    DEG count plot, CSV interpretation, and volcano colours:
    ``Up`` / ``Down`` / ``NS``.
    """
    import numpy as np

    logfc = np.asarray(logfc, dtype=float)
    padj = np.asarray(padj, dtype=float)
    labels = np.full(logfc.shape, 'NS', dtype=object)
    significant = np.isfinite(logfc) & np.isfinite(padj) & (padj < pval_cutoff)
    labels[significant & (logfc > logfc_cutoff)] = 'Up'
    labels[significant & (logfc < -logfc_cutoff)] = 'Down'
    return labels


def max_detection_fraction_by_gene(rank_result, var_names=None):
    """Map each gene to its maximum detection fraction across groups.

    Scanpy stores ``rank_genes_groups['pts']`` in original ``var_names`` order,
    while ``names`` is sorted by test rank.  Joining those structures by row
    position assigns detection fractions to the wrong genes.  This helper
    performs an index-aware lookup and retains a structured-array fallback for
    older Scanpy results.
    """
    import numpy as np
    import pandas as pd

    pts = rank_result.get('pts')
    if pts is None:
        return {}
    if isinstance(pts, pd.DataFrame):
        numeric = pts.apply(pd.to_numeric, errors='coerce')
        maxima = numeric.max(axis=1, skipna=True)
        return {
            str(gene): float(value)
            for gene, value in maxima.items() if np.isfinite(value)
        }
    values = np.asarray(pts)
    fields = values.dtype.names
    genes = [str(gene) for gene in ([] if var_names is None else var_names)]
    if fields and len(values) == len(genes):
        matrix = np.column_stack([
            np.asarray(values[field], dtype=float) for field in fields
        ])
        maxima = np.nanmax(matrix, axis=1)
        return {
            gene: float(value)
            for gene, value in zip(genes, maxima) if np.isfinite(value)
        }
    return {}


def _parse_deg_pairs(groups, raw_comparisons):
    """Parse explicit A-vs-B marker comparisons from user input."""
    import re

    group_set = {str(value) for value in groups}
    raw = str(raw_comparisons or '').strip()
    pairs = []
    if not raw:
        raise ValueError("custom 模式必须填写“指定比较（A-vs-B）”，例如 Cluster1-vs-Cluster5。")
    for item in re.split(r'[;\n]+', raw):
        item = item.strip()
        if not item:
            continue
        parts = re.split(r'\s*(?:-vs-|\s+vs\s+)\s*', item, maxsplit=1, flags=re.IGNORECASE)
        if len(parts) != 2 or not all(part.strip() for part in parts):
            raise ValueError(f"比较格式无效: '{item}'；请使用 A-vs-B。")
        experimental, control = (part.strip() for part in parts)
        if experimental not in group_set or control not in group_set:
            raise ValueError(f"比较组不存在于 '{' / '.join(sorted(group_set))}': '{item}'；"
                            "请使用与分组列取值完全一致的名称。")
        if experimental == control:
            raise ValueError(f"比较双方不能相同: '{item}'。")
        if (experimental, control) not in pairs:
            pairs.append((experimental, control))
    if not pairs:
        raise ValueError("未解析到有效比较；请使用 A-vs-B 格式，多个比较用分号或换行分隔。")
    return pairs

class DEGAnalysis(BaseAnalysis):
    MODULE_NAME = "deg"
    DISPLAY_NAME = "差异表达"
    DESCRIPTION = "差异表达基因分析（Wilcoxon 检验）"
    INPUT_REQUIRES = ['leiden']

    def run(self, input_path):
        import scanpy as sc
        import numpy as np
        import pandas as pd
        import os
        from modules.native_figures import (
            diverging_bar_figure, umap_panel_figure,
        )
        from figure_engine import NatureFigureDirector

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        adata = self.apply_scope(adata)
        requested_groupby = str(self.params.get('groupby', '') or '').strip()
        if not requested_groupby and 'celltype' in adata.obs.columns:
            requested_groupby = 'celltype'
        method = self.params.get('method', 'wilcoxon')
        n_genes = int(self.params.get('n_genes', 20))
        reference = str(self.params.get('reference', 'rest') or 'rest').strip() or 'rest'
        comparison_mode = str(self.params.get('comparison_mode', 'reference') or 'reference').strip()
        if comparison_mode not in {'reference', 'pairwise', 'custom'}:
            raise ValueError('comparison_mode 必须为 reference、pairwise 或 custom')
        comparisons_raw = str(self.params.get('comparisons', '') or '').strip()
        min_cells_pair = max(2, int(float(self.params.get('min_cells_per_group', 10) or 10)))
        pval_cutoff = float(self.params.get('pval_cutoff', 0.05))
        logfc_cutoff = float(self.params.get('logfc_cutoff', 1.0))
        min_pct = float(self.params.get('min_pct', 0.1))
        correction_method = self.params.get('correction_method', 'benjamini_hochberg')
        volcano_top_n = int(self.params.get('volcano_top_n', 10))
        volcano_genes_str = self.params.get('volcano_genes', '').strip()
        export_full_tables = self.params.get('export_full_tables', False)
        if isinstance(export_full_tables, str):
            export_full_tables = export_full_tables.strip().lower() not in {'', '0', 'false', 'no', 'off'}
        supported_methods = {'wilcoxon', 't-test', 'logreg', 't-test_overestim_var'}
        if method not in supported_methods:
            raise ValueError(f"不支持的 DEG method '{method}'")
        if n_genes < 1:
            raise ValueError('n_genes 必须至少为 1')
        if not 0 < pval_cutoff <= 1:
            raise ValueError('pval_cutoff 必须位于 (0, 1]')
        if not np.isfinite(logfc_cutoff) or logfc_cutoff < 0:
            raise ValueError('logfc_cutoff 必须为非负有限数值')
        if not 0 <= min_pct <= 1:
            raise ValueError('min_pct 必须位于 [0, 1]')

        groupby, group_info = resolve_obs_grouping(
            adata, requested_groupby, fallbacks=['celltype', 'leiden'],
            max_categories=50, max_numeric_categories=20,
            require_multiple=True,
        )
        if groupby is None:
            raise ValueError(
                f"groupby '{requested_groupby}' 不是有效的分类分组列："
                f"{group_info.get('reason', '')}"
            )
        if groupby != requested_groupby:
            self.progress(-1, f"分组列已改用 '{groupby}'：{group_info.get('requested_reason', '')}")
        if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
            adata.obs[groupby] = adata.obs[groupby].astype(str).astype('category')

        # 比较方案：reference（默认，每簇 vs 参考组/rest）、pairwise（全部两两）、
        # custom（只运行用户列出的 A-vs-B）。pairwise/custom 会先限定参与比较的组，
        # 使后续统计、图表和导出都只基于这些组，避免“选了比较但图里仍出现全部组”。
        group_labels = sorted(adata.obs[groupby].astype(str).unique().tolist())
        pairs = None
        if comparison_mode == "custom":
            pairs = _parse_deg_pairs(group_labels, comparisons_raw)
        elif comparison_mode == "pairwise":
            from itertools import combinations
            from modules.sc_batch import _looks_like_control
            raw_pairs = list(combinations(group_labels, 2))
            pairs = []
            for left, right in raw_pairs:
                if _looks_like_control(left) and not _looks_like_control(right):
                    pairs.append((right, left))
                else:
                    pairs.append((left, right))
            if len(pairs) > 200:
                raise ValueError(
                    f"pairwise 模式将生成 {len(pairs)} 个比较（超过 200 上限）；"
                    "请改用 custom 模式只保留关键比较。"
                )
        if pairs is not None:
            involved = set()
            for experimental, control in pairs:
                involved.add(experimental)
                involved.add(control)
            adata = adata[adata.obs[groupby].astype(str).isin(involved)].copy()
            if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
                adata.obs[groupby] = adata.obs[groupby].astype(str).astype('category')
            self.progress(-1, f"仅保留参与比较的 {len(involved)} 个组: {sorted(involved)}")

        # Pearson residuals are intentionally not used for fold-change
        # inference.  ``de_adata`` is a short-lived log1p reconstruction from
        # raw counts, leaving the analysis/embedding representation untouched.
        self.progress(18, "Preparing log1p expression from raw counts for interpretable log2FC...")
        de_adata = log1p_adata_for_cell_level_de(adata)
        if groupby not in de_adata.obs.columns:
            raise RuntimeError(f"DEG 临时表达对象缺少分组列 '{groupby}'")
        if not isinstance(de_adata.obs[groupby].dtype, pd.CategoricalDtype):
            de_adata.obs[groupby] = de_adata.obs[groupby].astype(str).astype('category')

        self.progress(20, f"Running DEG analysis ({method})...")
        # Scanpy 在内部已对“全部检验基因”做多重检验校正并给出 pvals_adj；
        # 避免在 Top-N 截断子集上二次校正导致 FDR 失真（子集内 m 远小于
        # 实际检验基因数，padj 会系统性偏松）。BY 不是 Scanpy 支持的
        # corr_method，因此先让 Scanpy 计算完整 p 值，再在完整检验宇宙上
        # 替换 pvals_adj。
        scanpy_corr_map = {
            'benjamini_hochberg': 'benjamini-hochberg',
            'bonferroni': 'bonferroni',
            'BY': 'benjamini-hochberg',
        }
        if correction_method not in scanpy_corr_map:
            raise ValueError(f"不支持的多重检验校正方法 '{correction_method}'")
        corr_method = scanpy_corr_map[correction_method]
        # ``sc.get.rank_genes_groups_df`` can only export the rows requested
        # here.  Keep the full tested universe in ``uns`` so that the persisted
        # ``sc_deg_full_results.csv`` is actually complete and can be used as
        # an auditable ranked input for downstream enrichment.  The top-N table
        # used for display is still limited below.
        rank_n_genes = de_adata.n_vars
        pts_kwarg = ({'pts': True}
                     if method in ('wilcoxon', 't-test', 't-test_overestim_var')
                     else {})
        reference_label = 'rest' if reference == 'rest' else str(reference)
        all_deg_frames = []
        skipped_pairs = []
        if pairs is not None:
            for index, (experimental, control) in enumerate(pairs, start=1):
                self.progress(
                    20 + int((index - 1) * 30 / max(len(pairs), 1)),
                    f"Running DEG ({method}): {experimental} vs {control}...",
                )
                counts = de_adata.obs[groupby].astype(str).value_counts()
                if (int(counts.get(experimental, 0)) < min_cells_pair
                        or int(counts.get(control, 0)) < min_cells_pair):
                    skipped_pairs.append(f"{experimental} vs {control}")
                    continue
                try:
                    sc.tl.rank_genes_groups(
                        de_adata, groupby=groupby, groups=[experimental],
                        reference=control, method=method, n_genes=rank_n_genes,
                        corr_method=corr_method, use_raw=False, **pts_kwarg)
                except KeyError as exc:
                    raise ValueError(
                        f"比较组 '{experimental} vs {control}' 无法执行: {exc}"
                    ) from exc
                pair_result = de_adata.uns['rank_genes_groups']
                if correction_method == 'BY':
                    _apply_by_adjustment(pair_result)
                pair_frame = sc.get.rank_genes_groups_df(de_adata, group=experimental)
                pair_frame['cluster'] = experimental
                pair_frame["comparison"] = f"{experimental} vs {control}"
                pair_frame['control_group'] = control
                pair_frame['max_pct_any_group'] = pair_frame['names'].astype(str).map(
                    max_detection_fraction_by_gene(pair_result, de_adata.var_names)
                )
                all_deg_frames.append(pair_frame)
            if skipped_pairs:
                self.progress(
                    -1,
                    f"跳过每组细胞数不足 {min_cells_pair} 的比较: {skipped_pairs}",
                )
            if not all_deg_frames:
                raise ValueError(
                    "所有指定比较都因每组细胞数不足而被跳过；"
                    f"请降低 min_cells_per_group（当前 {min_cells_pair}）或检查分组列。"
                )
            groups = []
            for frame in all_deg_frames:
                for value in frame['cluster'].astype(str).tolist():
                    if value not in groups:
                        groups.append(value)
        else:
            ref_kwarg = {} if reference == 'rest' else {'reference': str(reference)}
            try:
                sc.tl.rank_genes_groups(
                    de_adata, groupby=groupby, method=method,
                    n_genes=rank_n_genes, corr_method=corr_method,
                    use_raw=False, **pts_kwarg, **ref_kwarg)
            except KeyError as exc:
                raise ValueError(
                    f"参考组 '{reference}' 不存在于分组列 '{groupby}' 中: {exc}"
                ) from exc

            self.progress(50, "Extracting results...")
            result = de_adata.uns['rank_genes_groups']
            if correction_method == 'BY':
                _apply_by_adjustment(result)
            groups = list(result['names'].dtype.names)
            has_pvals = 'pvals' in result and 'pvals_adj' in result
            if not has_pvals:
                self.progress(-1, f"方法 {method} 不提供 p 值（如 logreg），"
                                  "显著性筛选与火山图按缺失 p 值处理。")

            # ``pts`` is indexed by var_names, not ranked position.  Build one
            # full rank table per group and join detection fractions by gene id.
            max_pct_expr = max_detection_fraction_by_gene(result, de_adata.var_names)
            if min_pct > 0 and not max_pct_expr:
                self.progress(-1, '当前方法未返回按基因索引的检测比例；min_pct 未应用。')
            for group in groups:
                group_frame = sc.get.rank_genes_groups_df(de_adata, group=group)
                group_frame['cluster'] = group
                group_frame["comparison"] = f"{group} vs {reference_label}"
                group_frame['control_group'] = reference_label
                group_frame['max_pct_any_group'] = group_frame['names'].astype(str).map(max_pct_expr)
                all_deg_frames.append(group_frame)

        all_deg_df = (
            pd.concat(all_deg_frames, ignore_index=True)
            if all_deg_frames else pd.DataFrame()
        )
        filtered_full = all_deg_df.copy()
        if min_pct > 0 and 'max_pct_any_group' in filtered_full.columns and not filtered_full.empty:
            filtered_full = filtered_full.loc[
                pd.to_numeric(
                    filtered_full['max_pct_any_group'], errors='coerce'
                ).ge(min_pct)
            ].copy()

        plot_deg_df = filtered_full.rename(columns={
            'names': 'gene', 'logfoldchanges': 'logfc', 'pvals': 'pval',
            'pvals_adj': 'pval_adj', 'scores': 'score',
        }).copy()
        for column in ('gene', 'logfc', 'pval', 'pval_adj', 'score', 'cluster'):
            if column not in plot_deg_df.columns:
                plot_deg_df[column] = np.nan
        plot_deg_df['regulation'] = classify_volcano_regulation(
            plot_deg_df['logfc'].to_numpy(),
            plot_deg_df['pval_adj'].to_numpy(),
            logfc_cutoff=logfc_cutoff,
            pval_cutoff=pval_cutoff,
        )

        top_frames = []
        if 'comparison' in plot_deg_df.columns:
            comparison_keys = sorted(plot_deg_df['comparison'].astype(str).unique().tolist())
            for key in comparison_keys:
                group_rows = plot_deg_df.loc[
                    plot_deg_df['comparison'].astype(str).eq(key)
                ].head(n_genes)
                top_frames.append(group_rows)
        else:
            for group in groups:
                group_rows = plot_deg_df.loc[
                    plot_deg_df['cluster'].astype(str).eq(str(group))
                ].head(n_genes)
                top_frames.append(group_rows)
        deg_df = (
            pd.concat(top_frames, ignore_index=True)
            if top_frames else plot_deg_df.head(0).copy()
        )
        deg_data = deg_df.to_dict(orient='records')

        self.progress(70, "Generating volcano plots...")
        plots_dir = self.ensure_plots_dir()
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        # Top-N 表始终落盘：这是模块的主要分析输出，不能依赖未出现在
        # schema 中的 export_full_tables 开关（GUI 默认路径永远不会写 CSV）。
        csv_path = os.path.join(results_dir, 'deg_results.csv')
        deg_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': f'Top {n_genes} DEG per Cluster'})

        marker_selection = {
            'cluster_markers': {}, 'genes': [], 'warnings': [],
            'n_data_driven': 0, 'n_classic_anchor': 0, 'n_relaxed': 0,
            'parameters': {},
        }
        custom_dotplot_str = self.params.get('custom_dotplot_genes', '').strip()
        if not custom_dotplot_str:
            try:
                from modules.annotation import select_cluster_marker_genes
                marker_selection = select_cluster_marker_genes(
                    de_adata,
                    groupby,
                    method=method,
                    n_rank_genes=max(1, min(de_adata.n_vars, max(100, int(n_genes)))),
                    min_markers=1,
                    max_markers=5,
                    padj_cutoff=pval_cutoff,
                    min_pct=min_pct,
                    min_delta_pct=0.05,
                )
                if marker_selection.get('warnings'):
                    self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
            except Exception as exc:
                self.progress(-1, f'DEG marker dotplot 选择失败，回退到排名结果：{exc}')

        volcano_units = []
        if 'comparison' in plot_deg_df.columns:
            volcano_units = sorted(plot_deg_df['comparison'].astype(str).unique().tolist())
        elif groups:
            volcano_units = [f"{groups[0]} vs {reference_label}"]
        volcano_units = volcano_units[:12]
        for unit in volcano_units:
            try:
                if 'comparison' in plot_deg_df.columns:
                    g_mask = plot_deg_df['comparison'].astype(str) == unit
                    unit_label = unit
                else:
                    g_mask = plot_deg_df['cluster'].astype(str) == str(groups[0])
                    unit_label = f"{groups[0]} vs {reference_label}"
                g_df = plot_deg_df[g_mask].copy()
                if g_df.empty:
                    self.progress(-1, f'{unit} 没有通过过滤的 DEG，跳过火山图。')
                    continue
                annotate_genes = []
                if volcano_genes_str:
                    annotate_genes = [
                        value.strip() for value in volcano_genes_str.replace('\n', ',').split(',')
                        if value.strip()
                    ]
                volcano_data = g_df.rename(columns={
                    'logfc': 'log2FC', 'pval_adj': 'padj',
                })
                director = NatureFigureDirector()
                volcano_spec = director.spec_from_params(
                    'volcano', self.params, width='single',
                    title=f'Cluster markers: {unit_label}',
                    evidence_role='discovery', fc_threshold=logfc_cutoff,
                    fdr_threshold=pval_cutoff, label_n=volcano_top_n,
                    label_genes=tuple(annotate_genes), show_legend=True,
                )
                fig = director.render(volcano_spec, volcano_data)
                fig._nature_semantic_warnings.append(
                    '该图为细胞级 cluster marker 检验；细胞不是独立生物学重复。'
                )
                import re as _re
                safe_unit = _re.sub(r'[^A-Za-z0-9_.-]+', '_', unit_label)
                result_files.extend(self.save_matplotlib_figure(
                    fig, plots_dir, f'deg_volcano_{safe_unit}.png', 'volcano',
                    f'Volcano: {unit_label}',
                ))
            except Exception as exc:
                self.progress(-1, f'DEG volcano plot failed: {exc}')

        if self.params.get('show_deg_counts_bar', True) and not plot_deg_df.empty:
            try:
                if 'comparison' in plot_deg_df.columns:
                    count_keys = sorted(plot_deg_df['comparison'].astype(str).unique().tolist())
                    up_counts = []
                    down_counts = []
                    for key in count_keys:
                        sub = plot_deg_df[plot_deg_df['comparison'].astype(str).eq(key)]
                        sig = sub['pval_adj'] < pval_cutoff
                        up_counts.append(int((sig & (sub['logfc'] > logfc_cutoff)).sum()))
                        down_counts.append(int((sig & (sub['logfc'] < -logfc_cutoff)).sum()))
                    fig_counts = diverging_bar_figure(
                        count_keys, up_counts, down_counts,
                        title=f'Significant DEG Counts (padj<{pval_cutoff}, |logFC|>{logfc_cutoff})',
                        x_label='comparison', y_label='Number of DEGs', rotation=35,
                    )
                else:
                    group_order = [str(g) for g in groups]
                    up_counts = []
                    down_counts = []
                    for g in group_order:
                        sub = plot_deg_df[plot_deg_df['cluster'].astype(str) == g]
                        sig = sub['pval_adj'] < pval_cutoff
                        up_counts.append(int((sig & (sub['logfc'] > logfc_cutoff)).sum()))
                        down_counts.append(int((sig & (sub['logfc'] < -logfc_cutoff)).sum()))
                    fig_counts = diverging_bar_figure(
                        group_order, up_counts, down_counts,
                        title=f'Significant DEG Counts (padj<{pval_cutoff}, |logFC|>{logfc_cutoff})',
                        x_label=groupby, y_label='Number of DEGs', rotation=35,
                    )
                result_files.extend(self.save_matplotlib_figure(
                    fig_counts, plots_dir, 'deg_significant_counts_bar.png',
                    'bar', 'Significant DEG Counts', formats=('png', 'svg'), dpi=300,
                ))
            except Exception as e:
                self.progress(-1, f"DEG count summary plot failed: {e}")

        # DEG Dotplot
        if self.params.get('show_dotplot', True):
            try:
                if custom_dotplot_str:
                    top_genes_list = [g.strip() for g in custom_dotplot_str.replace('\n', ',').split(',') if g.strip()]
                    not_found = [g for g in top_genes_list if g not in adata.var_names]
                    top_genes_list = [g for g in top_genes_list if g in adata.var_names]
                    dotplot_label = '自定义基因 Dotplot'
                else:
                    top_genes_list = marker_selection.get('genes', [])[:30]
                    not_found = []
                    if not top_genes_list:
                        group_key = 'comparison' if 'comparison' in plot_deg_df.columns else 'cluster'
                        for key, sub_df in plot_deg_df.groupby(group_key, observed=True):
                            top_genes_list.extend(sub_df.head(5)['gene'].astype(str).tolist())
                        top_genes_list = list(dict.fromkeys(top_genes_list))[:30]
                    dotplot_label = 'DEG Dotplot'

                if top_genes_list:
                    sc.tl.dendrogram(de_adata, groupby=groupby)
                    fig_dot = sc.pl.dotplot(de_adata, var_names=top_genes_list, groupby=groupby, return_fig=True)
                    result_files.extend(self.save_matplotlib_figure(
                        fig_dot, plots_dir, 'deg_dotplot.png', 'dotplot', dotplot_label
                    ))
                    import matplotlib.pyplot as plt
                    plt.close('all')
            except Exception as e:
                self.progress(-1, f"DEG dotplot generation failed: {e}")

        # Gene expression UMAP
        plot_genes_umap = self.params.get('plot_genes_umap', '').strip()
        if plot_genes_umap and 'X_umap' in adata.obsm:
            gene_list = [g.strip() for g in plot_genes_umap.split(',') if g.strip() and g.strip() in adata.var_names]
            for gene in gene_list[:5]:  # Limit to 5 genes
                try:
                    fig_static = self.build_publication_umap(
                        de_adata, gene, title=f'{gene} Expression (log1p counts)'
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_static, plots_dir, f'deg_gene_umap_{gene}.png', 'umap',
                        f'{gene} Expression'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close(fig_static)
                except Exception as exc:
                    self.progress(-1, f'{gene} 表达 UMAP 导出失败：{exc}')

        if self.params.get('show_top_marker_umap_panel', True) and 'X_umap' in de_adata.obsm and not plot_deg_df.empty:
            try:
                marker_limit = int(self.params.get('top_marker_umap_genes', 6))
                marker_candidates = []
                for g in groups:
                    sub = plot_deg_df[plot_deg_df['cluster'].astype(str) == str(g)]
                    sig_sub = sub[(sub['pval_adj'] < pval_cutoff) & (sub['logfc'] > logfc_cutoff)]
                    genes_for_group = sig_sub['gene'].astype(str).tolist() or sub['gene'].astype(str).tolist()
                    marker_candidates.extend(genes_for_group[:1])
                marker_genes = [g for g in dict.fromkeys(marker_candidates) if g in adata.var_names][:max(1, marker_limit)]
                if marker_genes:
                    fig_marker_umap = umap_panel_figure(
                        de_adata, marker_genes, titles=marker_genes,
                        point_size=4, opacity=0.70, ncols=min(3, len(marker_genes)),
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_marker_umap, plots_dir, 'deg_top_marker_umap_panel.png',
                        'umap', 'Top Marker Feature UMAPs', formats=('png', 'svg'), dpi=300,
                    ))
            except Exception as e:
                self.progress(-1, f"Top marker UMAP panel failed: {e}")

        # Export the complete tested table.  Top-N is only a display contract;
        # it must not truncate the ranked universe used for audit/enrichment.
        # 全量表使用与 Top-N 表一致的规范列名（gene/logfc/pval/pval_adj/…），
        # 避免同一模块的两张 CSV schema 不一致。
        if not all_deg_df.empty:
            full_export_df = all_deg_df.rename(columns={
                'names': 'gene', 'logfoldchanges': 'logfc', 'pvals': 'pval',
                'pvals_adj': 'pval_adj', 'scores': 'score',
            }).copy()
            for column in ('gene', 'logfc', 'pval', 'pval_adj', 'score', 'cluster'):
                if column not in full_export_df.columns:
                    full_export_df[column] = np.nan
            full_export_df['regulation'] = classify_volcano_regulation(
                full_export_df['logfc'].to_numpy(),
                full_export_df['pval_adj'].to_numpy(),
                logfc_cutoff=logfc_cutoff,
                pval_cutoff=pval_cutoff,
            )
            full_csv = os.path.join(results_dir, 'sc_deg_full_results.csv')
            full_export_df.to_csv(full_csv, index=False)
            result_files.append({'file_path': full_csv, 'file_type': 'csv', 'category': 'table', 'label': '完整 DEG 结果'})

        # Cluster marker heatmap
        if self.params.get('show_marker_heatmap', True) and not deg_df.empty:
            try:
                heatmap_top_n = int(self.params.get('marker_heatmap_top_n', 3))
                heatmap_genes = []
                if 'comparison' in deg_df.columns:
                    heatmap_keys = sorted(deg_df['comparison'].astype(str).unique().tolist())
                    heatmap_key_col = 'comparison'
                else:
                    heatmap_keys = list(groups)
                    heatmap_key_col = 'cluster'
                for g in heatmap_keys:
                    cluster_genes = deg_df.loc[deg_df[heatmap_key_col] == g, 'gene'].astype(str).tolist()
                    heatmap_genes.extend(cluster_genes[:max(1, heatmap_top_n)])
                heatmap_genes = [g for g in dict.fromkeys(heatmap_genes) if g in adata.var_names][:60]
                if heatmap_genes:
                    expr = de_adata[:, heatmap_genes].X
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray()
                    expr = np.asarray(expr, dtype=float)
                    group_labels = de_adata.obs[groupby].astype(str)
                    group_order = [str(g) for g in groups]
                    mean_matrix = []
                    for g in group_order:
                        mask = (group_labels == g).values
                        if mask.sum() == 0:
                            mean_matrix.append(np.zeros(len(heatmap_genes)))
                        else:
                            mean_matrix.append(expr[mask, :].mean(axis=0))
                    mean_matrix = np.asarray(mean_matrix).T
                    director = NatureFigureDirector()
                    heatmap_spec = director.spec_from_params(
                        'heatmap', self.params, width='double',
                        title=f'Top marker mean expression by {groupby}',
                        evidence_role='validation', zscore='row',
                        row_cluster=True, col_cluster=False,
                        max_row_labels=45, max_col_labels=24,
                    ).with_updates(height_mm=118.0)
                    fig_heat = director.render(heatmap_spec, {
                        'matrix': mean_matrix,
                        'gene_labels': heatmap_genes,
                        'sample_labels': group_order,
                        'colorbar_label': 'Gene-wise z-score',
                    })
                    fig_heat._nature_semantic_warnings.append(
                        '热图展示各 cluster 的平均表达，仅用于 marker/注释验证；完整细胞表达保留在 AnnData。'
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_heat, plots_dir, 'deg_marker_heatmap.png',
                        'heatmap', 'Cluster Marker Heatmap',
                    ))
            except Exception as e:
                self.progress(-1, f"Marker heatmap generation failed: {e}")

        self.progress(90, "Saving output...")
        adata.uns['marker_selection'] = {
            'cluster_key': groupby,
            'used_in_dotplot': not bool(custom_dotplot_str),
            'cluster_markers': marker_selection.get('cluster_markers', {}),
            'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
            'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
            'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
            'warnings': marker_selection.get('warnings', []),
            'parameters': marker_selection.get('parameters', {}),
        }
        output_path = self.save_output(adata, 'deg')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_groups': len(groups),
                'groups': list(groups),
                'method': method,
                'reference': reference,
                'comparison_mode': comparison_mode,
                'comparisons': (
                    sorted({f'{a} vs {b}' for a, b in pairs})
                    if pairs is not None
                    else ([f"{g} vs {reference_label}" for g in groups] if groups else [])
                ),
                'n_skipped_pairs': len(skipped_pairs),
                'min_cells_per_group': min_cells_pair,
                'scope_key': str(self.params.get('scope_key', '') or '').strip() or None,
                'scope_values': ([v.strip() for v in str(self.params.get('scope_values', '') or '').split(',') if v.strip()] or None),
                'total_deg_genes': len(deg_data),
                'n_tested_rows': int(len(all_deg_df)),
                'n_rows_after_min_pct': int(len(plot_deg_df)),
                'n_significant_rows_after_min_pct': int(
                    plot_deg_df['regulation'].isin(['Up', 'Down']).sum()
                ),
                'custom_dotplot_genes': self.params.get('custom_dotplot_genes', '').strip() or None,
                'requested_groupby': requested_groupby,
                'groupby': groupby,
                'inference_unit': 'cell',
                'statistical_status': 'exploratory_cluster_marker',
                'expression_scale': normalization_semantics(adata),
                'full_table_exported': bool(export_full_tables),
                'marker_selection': {
                    'used_in_dotplot': not bool(custom_dotplot_str),
                    'cluster_markers': marker_selection.get('cluster_markers', {}),
                    'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
                    'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
                    'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
                    'warnings': marker_selection.get('warnings', []),
                    'parameters': marker_selection.get('parameters', {}),
                },
            }
        }
