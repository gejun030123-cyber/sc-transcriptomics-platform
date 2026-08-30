from modules.base import BaseAnalysis
from modules.constants import G2M_GENES, S_GENES
from modules.io_utils import resolve_obs_grouping


VALID_HVG_FLAVORS = frozenset({'seurat_v3', 'cell_ranger', 'seurat'})
HVG_RESULT_COLUMNS = frozenset({
    'highly_variable',
    'highly_variable_rank',
    'highly_variable_rank_scanpy',
    'highly_variable_nbatches',
    'highly_variable_intersection',
    'means',
    'variances',
    'variances_norm',
    'dispersions',
    'dispersions_norm',
})


def _finite_matrix_values(matrix, max_values=200000):
    """Return a bounded, deterministic sample of finite matrix values."""
    import numpy as np
    from scipy import sparse

    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    values = np.asarray(values, dtype=float).ravel()
    if values.size > max_values:
        sample_index = np.linspace(0, values.size - 1, max_values, dtype=int)
        values = values[sample_index]
    return values[np.isfinite(values)]


def _is_count_matrix(matrix):
    """Return whether a matrix looks like non-negative integer UMI counts."""
    import numpy as np

    values = _finite_matrix_values(matrix)
    if values.size == 0:
        return False
    return bool(
        np.min(values) >= 0
        and np.all(np.isclose(values, np.rint(values), rtol=0, atol=1e-6))
    )


def _resolve_hvg_input(adata, flavor):
    """Select and validate the expression source required by a Scanpy flavor.

    ``seurat_v3`` models raw counts. The dispersion-based ``seurat`` and
    ``cell_ranger`` flavors model non-negative log-normalized expression.
    Returning ``None`` as the layer tells Scanpy to use ``adata.X``.
    """
    import numpy as np

    if flavor not in VALID_HVG_FLAVORS:
        choices = ', '.join(sorted(VALID_HVG_FLAVORS))
        raise ValueError(f"不支持的 HVG flavor: {flavor!r}；可选值为 {choices}")

    if flavor == 'seurat_v3':
        layer = 'counts' if 'counts' in adata.layers else None
        matrix = adata.layers[layer] if layer else adata.X
        if not _is_count_matrix(matrix):
            source = f"layers[{layer!r}]" if layer else 'X'
            raise ValueError(
                f"seurat_v3 要求原始非负整数 counts，但 {source} 不是计数矩阵。"
                "请保留 counts layer，或直接输入未经标准化的计数矩阵。"
            )
        return layer, layer or 'X'

    values = _finite_matrix_values(adata.X)
    if values.size == 0:
        raise ValueError(f"{flavor} 无法在空表达矩阵上选择高变基因")
    if np.min(values) < 0:
        raise ValueError(
            f"{flavor} 要求非负的 log-normalized X；当前 X 含负值，"
            "可能已经缩放或残差化。"
        )
    if _is_count_matrix(adata.X):
        raise ValueError(
            f"{flavor} 要求 log-normalized X，当前 X 看起来仍是原始整数 counts。"
            "请先运行标准化模块。"
        )
    return None, 'X'


def _clear_previous_hvg_results(var):
    """Remove stale Scanpy HVG annotations before changing flavor or batch."""
    stale_columns = [column for column in HVG_RESULT_COLUMNS if column in var.columns]
    if stale_columns:
        var.drop(columns=stale_columns, inplace=True)


def _candidate_priority(var, flavor):
    """Return gene positions ordered by Scanpy selection and variability."""
    import numpy as np
    import pandas as pd

    n_vars = len(var)
    selected = np.asarray(var['highly_variable'].fillna(False), dtype=bool)
    ordering = pd.DataFrame({'position': np.arange(n_vars, dtype=int)})
    ordering['selected_group'] = (~selected).astype(int)
    sort_columns = ['selected_group']
    ascending = [True]

    if flavor == 'seurat_v3' and 'highly_variable_rank' in var.columns:
        # 与 scanpy 的 seurat_v3 排序一致：先按中位 rank，再按支持批次数
        # 破平；此前把 n_batches 排在 rank 之前会改变选中基因集合。
        ordering['scanpy_rank'] = pd.to_numeric(
            var['highly_variable_rank'], errors='coerce'
        ).fillna(np.inf).to_numpy()
        sort_columns.append('scanpy_rank')
        ascending.append(True)

    if 'highly_variable_nbatches' in var.columns:
        ordering['n_batches'] = pd.to_numeric(
            var['highly_variable_nbatches'], errors='coerce'
        ).fillna(-1).to_numpy()
        sort_columns.append('n_batches')
        ascending.append(False)

    metric_column = (
        'variances_norm' if flavor == 'seurat_v3' else 'dispersions_norm'
    )
    if metric_column in var.columns:
        ordering['variability'] = pd.to_numeric(
            var[metric_column], errors='coerce'
        ).fillna(-np.inf).to_numpy()
        sort_columns.append('variability')
        ascending.append(False)

    sort_columns.append('position')
    ascending.append(True)
    return ordering.sort_values(
        sort_columns, ascending=ascending, kind='mergesort'
    )['position'].to_numpy(dtype=int)


def _finalize_hvg_selection(var, n_top_genes, flavor, excluded_mask, force_mask):
    """Apply exclusions/forced genes while preserving a fixed HVG count."""
    import numpy as np

    n_vars = len(var)
    target = min(int(n_top_genes), n_vars)
    if target <= 0:
        raise ValueError("n_top_genes 必须为正整数")

    excluded = np.asarray(excluded_mask, dtype=bool)
    forced = np.asarray(force_mask, dtype=bool)
    if excluded.shape != (n_vars,) or forced.shape != (n_vars,):
        raise ValueError("HVG 排除/强制包含掩码与基因数量不一致")

    # Force-include has precedence over exclusions, matching the documented UI
    # contract. An impossible request fails instead of silently exceeding n_top.
    allowed = ~excluded | forced
    forced_count = int(forced.sum())
    available_count = int(allowed.sum())
    if forced_count > target:
        raise ValueError(
            f"强制包含基因数 ({forced_count}) 超过目标 HVG 数 ({target})；"
            "请增加 n_top_genes 或减少强制基因。"
        )
    if available_count < target:
        raise ValueError(
            f"排除过滤后仅剩 {available_count} 个可用基因，无法选择 {target} 个 HVG。"
        )

    priority = _candidate_priority(var, flavor)
    forced_order = [position for position in priority if forced[position]]
    other_order = [
        position for position in priority
        if allowed[position] and not forced[position]
    ]
    selected_positions = forced_order + other_order[:target - forced_count]
    selected_set = set(selected_positions)
    final_order = [position for position in priority if position in selected_set]

    original_rank = var.get('highly_variable_rank')
    if original_rank is not None:
        var['highly_variable_rank_scanpy'] = original_rank.copy()

    var['highly_variable'] = False
    var.iloc[selected_positions, var.columns.get_loc('highly_variable')] = True
    final_rank = np.full(n_vars, np.nan, dtype=float)
    final_rank[final_order] = np.arange(len(final_order), dtype=float)
    var['highly_variable_rank'] = final_rank

    nonselected_order = [position for position in priority if position not in selected_set]
    return final_order + nonselected_order


class HVGAnalysis(BaseAnalysis):
    MODULE_NAME = "hvg"
    DISPLAY_NAME = "高变异基因选择"
    DESCRIPTION = "选择高变异基因（HVG），支持批次感知和基因过滤"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        try:
            n_hvg = int(self.params.get('n_top_genes', 2000))
        except (TypeError, ValueError):
            return "n_top_genes 必须为正整数"
        if n_hvg <= 0:
            return "n_top_genes 必须为正整数"

        flavor = str(self.params.get('hvg_flavor', 'seurat_v3') or '').strip()
        try:
            _resolve_hvg_input(adata, flavor)
        except ValueError as exc:
            return str(exc)

        requested_batch = str(self.params.get('batch_key', '') or '').strip()
        if requested_batch:
            _, batch_info = resolve_obs_grouping(
                adata, requested_batch, max_categories=50,
                max_numeric_categories=20, require_multiple=True,
            )
            if not batch_info.get('requested_valid', False):
                return (
                    f"batch_key={requested_batch!r} 无效："
                    f"{batch_info.get('requested_reason', '不是有效分类列')}"
                )

        if self.params.get('regress_cc', False):
            score_columns = {'S_score', 'G2M_score'}
            has_scores = score_columns.issubset(adata.obs.columns)
            if not has_scores and not self.params.get('cc_scoring', False):
                return (
                    "regress_cc 需要 S_score/G2M_score。请开启细胞周期评分，"
                    "或使用已包含这两列的上游结果。"
                )
            if _is_count_matrix(adata.X):
                return (
                    "regress_cc 不能作用于原始 counts。请先运行 log1p 标准化，"
                    "再进行细胞周期回归。"
                )
        return None

    def run(self, input_path):
        import os
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from modules.native_figures import scatter_figure

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)

        validation_error = self.validate_input(adata)
        if validation_error:
            raise ValueError(validation_error)

        requested_n_hvg = int(self.params.get('n_top_genes', 2000))
        n_hvg = min(requested_n_hvg, adata.n_vars)
        if n_hvg < requested_n_hvg:
            self.progress(
                -1,
                f"请求 {requested_n_hvg} 个 HVG，但数据仅有 {adata.n_vars} 个基因；"
                f"目标已调整为 {n_hvg}。",
            )

        requested_batch_key = str(self.params.get('batch_key', '') or '').strip()
        batch_key = None
        if requested_batch_key:
            batch_key, _ = resolve_obs_grouping(
                adata, requested_batch_key, max_categories=50,
                max_numeric_categories=20, require_multiple=True,
            )
            if not hasattr(adata.obs[batch_key].dtype, 'categories'):
                adata.obs[batch_key] = adata.obs[batch_key].astype(str).astype('category')

        hvg_flavor = str(self.params.get('hvg_flavor', 'seurat_v3')).strip()
        hvg_layer, hvg_input = _resolve_hvg_input(adata, hvg_flavor)
        exclude_mt = bool(self.params.get('exclude_mt_genes', False))
        exclude_cc = bool(self.params.get('exclude_cc_genes', False))
        force_genes_str = str(self.params.get('force_include_genes', '') or '').strip()
        cc_scoring = bool(self.params.get('cc_scoring', False))
        regress_cc = bool(self.params.get('regress_cc', False))

        self.progress(20, f"Selecting HVGs (flavor={hvg_flavor})...")
        hvg_kwargs = {'n_top_genes': n_hvg, 'flavor': hvg_flavor}
        if hvg_layer is not None:
            hvg_kwargs['layer'] = hvg_layer
        if batch_key:
            hvg_kwargs['batch_key'] = batch_key
        # Scanpy already performs within-batch selection and principled merging.
        # Do not repeat per-batch calls or overwrite its integrated ranking.
        _clear_previous_hvg_results(adata.var)
        sc.pp.highly_variable_genes(adata, **hvg_kwargs)

        force_genes = []
        if force_genes_str:
            parsed = [
                gene.strip()
                for gene in force_genes_str.replace('\n', ',').split(',')
                if gene.strip()
            ]
            force_genes = list(dict.fromkeys(parsed))
        var_names = list(adata.var_names.astype(str))
        var_name_to_position = {gene: index for index, gene in enumerate(var_names)}
        force_found = [gene for gene in force_genes if gene in var_name_to_position]
        force_missing = [gene for gene in force_genes if gene not in var_name_to_position]
        if force_missing:
            self.progress(
                -1,
                "以下强制包含基因未在数据中找到，已忽略："
                + ', '.join(force_missing[:20]),
            )

        mt_mask = np.asarray(
            adata.var_names.str.upper().str.startswith('MT-'), dtype=bool
        )
        cc_names = {str(gene).upper() for gene in S_GENES + G2M_GENES}
        cc_mask = np.asarray(
            [gene.upper() in cc_names for gene in var_names], dtype=bool
        )
        excluded_mask = np.zeros(adata.n_vars, dtype=bool)
        if exclude_mt:
            excluded_mask |= mt_mask
        if exclude_cc:
            excluded_mask |= cc_mask
        force_mask = np.zeros(adata.n_vars, dtype=bool)
        for gene in force_found:
            force_mask[var_name_to_position[gene]] = True

        display_order = _finalize_hvg_selection(
            adata.var, n_hvg, hvg_flavor, excluded_mask, force_mask,
        )

        cell_cycle_scored = False
        if cc_scoring:
            self.progress(50, "Scoring cell cycle...")
            var_names_set = set(var_names)
            s_in = [gene for gene in S_GENES if gene in var_names_set]
            g2m_in = [gene for gene in G2M_GENES if gene in var_names_set]
            if len(s_in) >= 5 and len(g2m_in) >= 5:
                adata_cc = adata.copy()
                if _is_count_matrix(adata_cc.X):
                    sc.pp.normalize_total(adata_cc, target_sum=1e4)
                    sc.pp.log1p(adata_cc)
                sc.tl.score_genes_cell_cycle(
                    adata_cc, s_genes=s_in, g2m_genes=g2m_in
                )
                adata.obs['S_score'] = adata_cc.obs['S_score']
                adata.obs['G2M_score'] = adata_cc.obs['G2M_score']
                adata.obs['phase'] = adata_cc.obs['phase']
                cell_cycle_scored = True
                del adata_cc
            else:
                message = (
                    "细胞周期评分所需基因不足："
                    f"S={len(s_in)}, G2M={len(g2m_in)}（各至少需要 5 个）"
                )
                if regress_cc:
                    raise ValueError(message + "；无法执行 regress_cc。")
                self.progress(-1, message + "；已跳过评分。")

        cell_cycle_regressed = False
        if regress_cc:
            missing_scores = [
                column for column in ('S_score', 'G2M_score')
                if column not in adata.obs.columns
            ]
            if missing_scores:
                raise ValueError(
                    "regress_cc 缺少评分列：" + ', '.join(missing_scores)
                )
            if _is_count_matrix(adata.X):
                raise ValueError(
                    "regress_cc 不能作用于原始 counts；请先运行 log1p 标准化。"
                )
            self.progress(65, "Regressing out cell cycle...")
            # sc.pp.regress_out 会把稀疏矩阵稠密化（约 n_cells × n_genes ×
            # 4 字节），大样本有 OOM 风险；超过安全上限时拒绝执行并给出
            # 替代方案，而不是把整机内存耗尽。
            dense_elements = int(adata.n_obs) * int(adata.n_vars)
            max_elements = int(os.environ.get('HVG_REGRESS_MAX_ELEMENTS', '200000000'))
            if dense_elements > max_elements:
                raise ValueError(
                    f"regress_cc 需要将矩阵稠密化（{adata.n_obs} 细胞 × "
                    f"{adata.n_vars} 基因 ≈ {dense_elements * 4 / 1024**3:.1f} GB），"
                    f"超过安全上限 {max_elements * 4 / 1024**3:.1f} GB。"
                    "建议：改用 pearson_residuals 标准化并在下游 regress，"
                    "或设置环境变量 HVG_REGRESS_MAX_ELEMENTS 提高上限（需自担内存风险）。"
                )
            sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
            cell_cycle_regressed = True

        hvg_metadata = dict(adata.uns.get('hvg', {}) or {})
        batch_support = {}
        if batch_key:
            if 'highly_variable_nbatches' in adata.var.columns:
                support = pd.to_numeric(
                    adata.var['highly_variable_nbatches'], errors='coerce'
                )
                distribution = {
                    str(int(key)): int(value)
                    for key, value in support.dropna().astype(int).value_counts().sort_index().items()
                }
                n_batches = int(adata.obs[batch_key].nunique())
                batch_support = {
                    'batch_key': batch_key,
                    'n_batches': n_batches,
                    'support_column': 'highly_variable_nbatches',
                    'support_distribution': distribution,
                    'n_supported_in_all_batches': int((support >= n_batches).sum()),
                    'n_supported_in_at_least_two_batches': int((support >= min(2, n_batches)).sum()),
                }
            else:
                batch_support = {
                    'batch_key': batch_key,
                    'n_batches': int(adata.obs[batch_key].nunique()),
                    'warning': 'Scanpy 未返回 highly_variable_nbatches，无法提供逐基因 batch 支持度。',
                }
        rank_column = 'highly_variable_rank' if 'highly_variable_rank' in adata.var.columns else None
        rank_available = int(pd.to_numeric(adata.var[rank_column], errors='coerce').notna().sum()) if rank_column else 0
        hvg_metadata.update({
            'flavor': hvg_flavor,
            'input': hvg_input,
            'batch_key': batch_key or '',
            'batch_selection': 'scanpy' if batch_key else 'none',
            'n_top_genes': int(n_hvg),
            'force_include_genes': force_found,
            'force_missing_genes': force_missing,
            'exclude_mt_genes': exclude_mt,
            'exclude_cc_genes': exclude_cc,
            'requested_n_top_genes': int(requested_n_hvg),
            'final_hvg_number': int(adata.var['highly_variable'].sum()),
            'highly_variable_rank_column': rank_column,
            'rank_available_genes': rank_available,
            'batch_support': batch_support,
        })
        selected_mask = adata.var['highly_variable'].fillna(False).astype(bool).to_numpy()
        selected_batch_support = {}
        if batch_key and 'highly_variable_nbatches' in adata.var.columns:
            selected_support = pd.to_numeric(
                adata.var.loc[selected_mask, 'highly_variable_nbatches'], errors='coerce'
            ).dropna().astype(int)
            selected_batch_support = {
                str(int(key)): int(value)
                for key, value in selected_support.value_counts().sort_index().items()
            }
            batch_support['selected_support_distribution'] = selected_batch_support
            batch_support['selected_supported_in_all_batches'] = int(
                (selected_support >= int(batch_support.get('n_batches', 0))).sum()
            )
            batch_support['selected_supported_in_at_least_two_batches'] = int(
                (selected_support >= min(2, int(batch_support.get('n_batches', 0)))).sum()
            )
        technical_composition = {
            'MT': int(np.asarray(
                adata.var['mt'].fillna(False) if 'mt' in adata.var.columns
                else adata.var_names.astype(str).str.upper().str.startswith('MT-'),
                dtype=bool,
            )[selected_mask].sum()),
            'Ribo': int(np.asarray(
                adata.var['ribo'].fillna(False) if 'ribo' in adata.var.columns
                else np.zeros(adata.n_vars, dtype=bool),
                dtype=bool,
            )[selected_mask].sum()),
            'Cell cycle': int(sum(
                str(gene).upper() in {str(item).upper() for item in S_GENES + G2M_GENES}
                for gene in adata.var_names[selected_mask]
            )),
        }
        technical_composition = {
            key: value for key, value in technical_composition.items() if value > 0
        }
        hvg_metadata['selected_batch_support_distribution'] = selected_batch_support
        hvg_metadata['selected_technical_composition'] = technical_composition
        adata.uns['hvg'] = hvg_metadata

        self.progress(75, "Generating HVG plot...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        var_df = adata.var.copy()
        plot_cols = [
            column for column in (
                'variances_norm', 'variances', 'dispersions_norm', 'dispersions'
            )
            if column in var_df.columns
        ]
        y_col = plot_cols[0] if plot_cols else None
        if y_col:
            x_values = (
                var_df['means'].values
                if 'means' in var_df.columns
                else np.arange(len(var_df))
            )
            colors = var_df['highly_variable'].map(
                {True: '#B64342', False: '#98A2B3'}
            ).values
            fig = scatter_figure(
                x_values, var_df[y_col].values,
                title='Highly Variable Genes', x_label='Mean', y_label=y_col,
                colors=colors, size=8, alpha=0.75,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'hvg_scatter.png', 'scatter',
                'Highly Variable Genes', formats=('png', 'svg'), dpi=300,
            ))

            if self.params.get('show_hvg_rank_plot', True):
                ranked = adata.var.iloc[display_order].copy()
                ranked['rank'] = np.arange(1, len(ranked) + 1)
                rank_colors = ranked['highly_variable'].map(
                    {True: '#B64342', False: '#98A2B3'}
                ).values
                fig_rank = scatter_figure(
                    ranked['rank'].values, ranked[y_col].values,
                    title='HVG Rank Plot', x_label='Gene rank', y_label=y_col,
                    colors=rank_colors, size=9, alpha=0.75,
                    labels=ranked.index.astype(str).tolist(), annotate_top=15,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_rank, plots_dir, 'hvg_rank_plot.png', 'scatter',
                    'HVG Rank Plot', formats=('png', 'svg'), dpi=300,
                ))

        from modules.sc_figure_diagnostics import (
            hvg_batch_support_figure,
            hvg_technical_composition_figure,
        )
        if selected_batch_support:
            fig_support = hvg_batch_support_figure(
                selected_batch_support,
                batch_support.get('n_batches', max(map(int, selected_batch_support))),
                n_hvg,
            )
            if fig_support is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_support, plots_dir, 'hvg_batch_support_selected.png',
                    'bar', 'Batch support among selected HVGs',
                    formats=('png', 'svg'), dpi=300,
                ))
        if technical_composition:
            fig_technical = hvg_technical_composition_figure(
                technical_composition, n_hvg,
            )
            if fig_technical is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_technical, plots_dir, 'hvg_technical_composition.png',
                    'bar', 'Technical composition of selected HVGs',
                    formats=('png', 'svg'), dpi=300,
                ))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'hvg')

        n_hvg_actual = int(adata.var['highly_variable'].sum())
        hvg_evidence = {
            'flavor': hvg_flavor,
            'requested_n_top_genes': int(requested_n_hvg),
            'final_hvg_number': n_hvg_actual,
            'highly_variable_rank_column': rank_column,
            'rank_available_genes': rank_available,
            'batch_support': batch_support,
        }
        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': adata.n_obs,
                'n_genes_total': adata.n_vars,
                'n_hvgs': n_hvg_actual,
                'final_hvg_number': n_hvg_actual,
                'requested_n_top_genes': requested_n_hvg,
                'hvg_flavor': hvg_flavor,
                'hvg_input': hvg_input,
                'requested_batch_key': requested_batch_key,
                'batch_key': batch_key,
                'batch_selection': 'scanpy' if batch_key else 'none',
                'hvg_evidence': hvg_evidence,
                'selected_batch_support_distribution': selected_batch_support,
                'selected_supported_in_all_batches': batch_support.get('selected_supported_in_all_batches'),
                'selected_supported_in_at_least_two_batches': batch_support.get('selected_supported_in_at_least_two_batches'),
                'selected_technical_composition': technical_composition,
                'force_include_count': len(force_found),
                'force_missing_genes': force_missing,
                'cell_cycle_scored': cell_cycle_scored,
                'cell_cycle_regressed': cell_cycle_regressed,
            },
        }
