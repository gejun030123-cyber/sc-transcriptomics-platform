import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


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


def _run_single_comparison(adata, counts, group1_samples, group2_samples, group1, group2,
                           method, fc_threshold, pval_threshold, top_n, gene_id_to_name,
                           plots_dir, results_dir, suffix='', viz_params=None,
                           cooks_filter=True, independent_filter=True, padj_method='fdr_bh',
                           base_mean_filter=0):
    """Run DEG for one comparison pair. Returns (deg_df, result_files, n_up, n_down)."""
    import plotly.graph_objects as go
    import omicverse as ov

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

    # 仅对简单统计方法调用 normalize；DESeq2/edgeR/limma 内部自行处理标准化
    if ov_method in ('ttest', 'wilcox'):
        dds.normalize()

    # cooks_filter / independent_filter 仅对 DESeq2 有意义
    deg_kwargs = {'multipletests_method': padj_method}
    if ov_method == 'DEseq2':
        deg_kwargs['cooks_filter'] = cooks_filter
        deg_kwargs['independent_filter'] = independent_filter

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

    # 保存过滤前的 top 基因用于火山图标注（与散点数据一致）
    top_genes_vol = deg_df[deg_df['regulation'] != 'NS'].head(top_n)

    # 基础表达量过滤
    if base_mean_filter > 0:
        base_mean = (deg_df['mean_group1'] + deg_df['mean_group2']) / 2
        deg_df = deg_df[base_mean >= base_mean_filter].copy()

    # Volcano plot
    file_suffix = f'_{suffix}' if suffix else ''
    neg_log_padj = -np.log10(padj + 1e-300)
    color_map = {'Up': '#e53935', 'Down': '#1a237e', 'NS': '#bdbdbd'}
    fig_vol = go.Figure()
    for reg in ['NS', 'Up', 'Down']:
        idx = [i for i in range(n_genes) if regulation[i] == reg]
        fig_vol.add_trace(go.Scattergl(
            x=log2fc[idx], y=neg_log_padj[idx], mode='markers',
            marker=dict(color=color_map[reg], size=5, opacity=0.7),
            name=f'{reg} ({len(idx)})',
            text=[gene_names[i] for i in idx],
            hovertemplate='%{text}<br>log2FC: %{x:.2f}<br>-log10(padj): %{y:.2f}'
        ))
    fig_vol.add_hline(y=-np.log10(pval_threshold), line_dash='dash', line_color='gray')
    fig_vol.add_vline(x=log2fc_threshold, line_dash='dash', line_color='gray')
    fig_vol.add_vline(x=-log2fc_threshold, line_dash='dash', line_color='gray')
    # Gene annotations（使用与散点图一致的过滤前数据）
    for _, row in top_genes_vol.iterrows():
        fig_vol.add_annotation(
            x=row['log2FC'], y=-np.log10(max(row['padj'], 1e-300)),
            text=row['gene'], showarrow=True, arrowhead=2,
            font=dict(size=max(7, font_size - 3), color='#333'), ax=20, ay=-30
        )
    fig_vol.update_layout(
        title=f'火山图 ({group1} vs {group2})',
        xaxis_title='log2(Fold Change)', yaxis_title='-log10(padj)',
        plot_bgcolor=bg_color, width=fig_width, height=fig_height,
        font=dict(family=font_family, size=font_size),
    )
    fpath = os.path.join(plots_dir, f'bulk_deg_volcano{file_suffix}.json')
    with open(fpath, 'w') as f:
        f.write(fig_vol.to_json(engine="json"))
    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'volcano',
                         'label': f'火山图 ({group1} vs {group2})'})

    # MA plot
    avg_expr = (mean1 + mean2) / 2
    fig_ma = go.Figure()
    for reg in ['NS', 'Up', 'Down']:
        idx = [i for i in range(n_genes) if regulation[i] == reg]
        fig_ma.add_trace(go.Scattergl(
            x=np.log2(avg_expr[idx] + 1), y=log2fc[idx], mode='markers',
            marker=dict(color=color_map[reg], size=5, opacity=0.7),
            name=f'{reg}', text=[gene_names[i] for i in idx],
            hovertemplate='%{text}<br>AvgExpr: %{x:.2f}<br>log2FC: %{y:.2f}'
        ))
    fig_ma.add_hline(y=0, line_dash='solid', line_color='gray', line_width=0.5)
    fig_ma.update_layout(
        title=f'MA 图 ({group1} vs {group2})', xaxis_title='log2(Average Expression)',
        yaxis_title='log2(Fold Change)',
        plot_bgcolor=bg_color, width=fig_width, height=fig_height,
        font=dict(family=font_family, size=font_size),
    )
    fpath = os.path.join(plots_dir, f'bulk_deg_ma{file_suffix}.json')
    with open(fpath, 'w') as f:
        f.write(fig_ma.to_json(engine="json"))
    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'ma',
                         'label': f'MA 图 ({group1} vs {group2})'})

    # Save individual CSV (保存过滤前的完整结果)
    csv_path = os.path.join(results_dir, f'bulk_deg_results{file_suffix}.csv')
    deg_df.to_csv(csv_path, index=False)
    result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                         'label': f'差异表达基因列表 ({group1} vs {group2})'})

    top_genes = deg_df[deg_df['regulation'] != 'NS'].head(top_n)
    top_csv = os.path.join(results_dir, f'bulk_deg_top_genes{file_suffix}.csv')
    top_genes.to_csv(top_csv, index=False)
    result_files.append({'file_path': top_csv, 'file_type': 'csv', 'category': 'table',
                         'label': f'Top {top_n} 差异基因 ({group1} vs {group2})'})

    # 返回完整 deg_df（调用方负责按需过滤）
    return deg_df, result_files, n_up, n_down


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

        groupby = self.params.get('groupby', 'condition')
        group1 = self.params.get('group1', '')
        group2 = self.params.get('group2', '')
        method = self.params.get('method', 't-test')
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

        # 确定两组样本名（延迟到 single-comparison 分支处理，避免重复逻辑）
        if groupby not in adata.obs.columns:
            n = counts.shape[0]
            half = n // 2
            if half == 0 or half == n:
                raise ValueError(f"未找到分组列 '{groupby}'，样本数不足。")
            adata.obs[groupby] = pd.Series(
                ['Group1']*half + ['Group2']*(n-half), index=adata.obs.index
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
                    from plotly.subplots import make_subplots

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

                    # Volcano 并排展示
                    import plotly.graph_objects as go
                    n_comp = len(comparison_names)
                    fig_multi = make_subplots(rows=1, cols=n_comp, subplot_titles=comparison_names,
                                               horizontal_spacing=0.05)
                    for m_idx, (comp_name, deg_df) in enumerate(zip(comparison_names, all_deg_dfs)):
                        colors = ['#e53935' if r == 'Up' else '#1565c0' if r == 'Down' else '#9e9e9e'
                                  for r in deg_df['regulation']]
                        neg_log_p = -np.log10(deg_df['padj'].values + 1e-300)
                        fig_multi.add_trace(go.Scattergl(
                            x=deg_df['log2FC'].tolist(), y=neg_log_p.tolist(),
                            mode='markers', marker=dict(color=colors, size=4),
                            text=deg_df['gene'].tolist(), showlegend=False),
                            row=1, col=m_idx + 1)
                        fig_multi.update_xaxes(title_text='log2FC', row=1, col=m_idx + 1)
                        if m_idx == 0:
                            fig_multi.update_yaxes(title_text='-log10(padj)', row=1, col=1)
                    fig_multi.update_layout(height=400, width=max(600, 300 * n_comp), title='多组比较 Volcano 图')
                    fpath = os.path.join(plots_dir, 'bulk_deg_volcano_multi.json')
                    with open(fpath, 'w') as f:
                        f.write(fig_multi.to_json(engine="json"))
                    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'volcano', 'label': '多组比较 Volcano'})
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
                'comparisons': [f'{g1}-vs-{g2}' for g1, g2 in comparison_pairs],
                'per_comparison': per_comparison,
                'shared_up_genes': shared_up,
                'shared_down_genes': shared_down,
                'skipped_comparisons': skipped,
                'fc_threshold': fc_threshold,
                'pval_threshold': pval_threshold,
                'lrt_n_sig': lrt_n_sig,
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
                n = counts.shape[0]
                half = n // 2
                if half == 0 or half == n:
                    raise ValueError(f"未找到分组列 '{groupby}'，样本数不足。")
                group1_samples = list(adata.obs.index[:half])
                group2_samples = list(adata.obs.index[half:])
                group1, group2 = "Group1", "Group2"
                group2_label = group2
                adata.obs[groupby] = pd.Series(
                    ['Group1']*half + ['Group2']*(n-half), index=adata.obs.index
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
        import plotly.graph_objects as go

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
            fig_box = go.Figure()
            for grp_name, samples in [(group1, group1_samples), (group2, group2_samples)]:
                sample_indices = [obs_idx[s] for s in samples if s in obs_idx]
                vals = counts[sample_indices, var_idx[pg_id]]
                fig_box.add_trace(go.Box(y=vals, name=str(grp_name), boxpoints='all', jitter=0.3))
            fig_box.update_layout(title=f'{pg} 表达', yaxis_title='Expression',
                                 plot_bgcolor='white', width=400, height=350)
            fpath = os.path.join(plots_dir, f'bulk_deg_box_{pg}.json')
            with open(fpath, 'w') as f: f.write(fig_box.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'boxplot', 'label': f'{pg} Boxplot'})
        if unmatched_genes:
            import logging
            logging.getLogger(__name__).warning(f"[bulk_deg] 以下基因未找到，已跳过箱线图: {', '.join(unmatched_genes)}")
