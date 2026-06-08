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
            parts = item.split('-vs-')
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
                           plots_dir, results_dir, suffix='', viz_params=None):
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
    dds.drop_duplicates_index()

    method_map = {'t-test': 'ttest', 'mann-whitney': 'wilcox', 'deseq2': 'DEseq2'}
    ov_method = method_map.get(method, 'ttest')

    dds.normalize()
    result = dds.deg_analysis(group1_samples, group2_samples, method=ov_method)

    # Extract results
    gene_ids_list = result.index.tolist()
    if gene_id_to_name:
        gene_names = [gene_id_to_name.get(g, g) for g in gene_ids_list]
    else:
        gene_names = gene_ids_list
    log2fc = result['log2FC'].values
    pvalues = result['pvalue'].values
    padj = result['qvalue'].values
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

    # Group means
    mask1_arr = np.array([s in group1_samples for s in adata.obs.index])
    mask2_arr = np.array([s in group2_samples for s in adata.obs.index])
    mean1 = counts[mask1_arr].mean(axis=0)
    mean2 = counts[mask2_arr].mean(axis=0)

    deg_df = pd.DataFrame({
        'gene': gene_names,
        'log2FC': np.round(log2fc, 4),
        'pvalue': pvalues,
        'padj': padj,
        'mean_group1': np.round(mean1[:n_genes] if len(mean1) >= n_genes else mean1, 2),
        'mean_group2': np.round(mean2[:n_genes] if len(mean2) >= n_genes else mean2, 2),
        'regulation': regulation
    })
    deg_df = deg_df.sort_values('padj')

    n_up = sum(1 for r in regulation if r == 'Up')
    n_down = sum(1 for r in regulation if r == 'Down')

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
    # Gene annotations
    top_genes_vol = deg_df[deg_df['regulation'] != 'NS'].head(top_n)
    for _, row in top_genes_vol.iterrows():
        fig_vol.add_annotation(
            x=row['log2FC'], y=-np.log10(max(row['padj'], 1e-300)),
            text=row['gene'], showarrow=True, arrowhead=2,
            font=dict(size=9, color='#333'), ax=20, ay=-30
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
    avg_expr = (mean1[:n_genes] + mean2[:n_genes]) / 2 if len(mean1) >= n_genes else (mean1 + mean2) / 2
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

    # Save individual CSV
    csv_path = os.path.join(results_dir, f'bulk_deg_results{file_suffix}.csv')
    deg_df.to_csv(csv_path, index=False)
    result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                         'label': f'差异表达基因列表 ({group1} vs {group2})'})

    top_genes = deg_df[deg_df['regulation'] != 'NS'].head(top_n)
    top_csv = os.path.join(results_dir, f'bulk_deg_top_genes{file_suffix}.csv')
    top_genes.to_csv(top_csv, index=False)
    result_files.append({'file_path': top_csv, 'file_type': 'csv', 'category': 'table',
                         'label': f'Top {top_n} 差异基因 ({group1} vs {group2})'})

    return deg_df, result_files, n_up, n_down


class BulkDEGAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_deg"
    DISPLAY_NAME = "Bulk 差异表达分析"
    DESCRIPTION = "组间差异表达基因检测：t-test / Mann-Whitney / DESeq2（基于 OmicVerse）"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
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

        self.progress(15, "构建计数矩阵...")

        # 构建基因ID→基因名映射
        gene_id_to_name = {}
        if 'gene_name' in adata.var.columns:
            for gid, gname in zip(adata.var_names, adata.var['gene_name']):
                if pd.notna(gname) and str(gname).strip():
                    gene_id_to_name[str(gid)] = str(gname).strip()

        # 构建 OmicVerse pyDEG 所需的 counts DataFrame（基因×样本）
        counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        counts = np.nan_to_num(counts.astype(float), nan=0.0, posinf=0.0, neginf=0.0)

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

        # 确定两组样本名
        if groupby in adata.obs.columns:
            groups = adata.obs[groupby].unique().tolist()
            if not group1 or group1 not in groups:
                group1 = groups[0]
            if group2 == 'rest' or (not group2 or group2 not in groups):
                group2_samples = [s for s in adata.obs.index if adata.obs.loc[s, groupby] != group1]
                group2 = f'rest (n={len(group2_samples)})'
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
            adata.obs[groupby] = pd.Series(
                ['Group1']*half + ['Group2']*(n-half), index=adata.obs.index
            )

        self.progress(25, "准备差异分析...")

        plots_dir = os.path.join(self.project_dir, 'plots')
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        if comparison_pairs:
            # 多组比较模式
            all_deg_dfs = []
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
                    viz_params=self.params.get('_visualization', {}))
                all_deg_dfs.append(deg_df)
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

                # 箱线图（取第一个比较的结果生成）
                first_deg = all_deg_dfs[0]
                g1_first, g2_first = comparison_pairs[0]
                g1_samples = list(adata.obs.index[adata.obs[groupby] == g1_first])
                g2_samples = list(adata.obs.index[adata.obs[groupby] == g2_first])
                self._draw_boxplots(first_deg, g1_first, g2_first, g1_samples, g2_samples,
                                    gene_id_to_name, counts, adata, plots_dir, result_files)
            else:
                raise ValueError("所有比较均因样本数不足被跳过，请检查分组信息。")

            summary = {
                'comparisons': [f'{g1}-vs-{g2}' for g1, g2 in comparison_pairs],
                'n_comparisons': len(all_deg_dfs),
                'skipped_comparisons': skipped,
                'method': method,
                'fc_threshold': fc_threshold,
                'pval_threshold': pval_threshold,
            }
        else:
            # 单次比较模式
            if groupby in adata.obs.columns:
                groups = adata.obs[groupby].unique().tolist()
                if not group1 or group1 not in groups:
                    group1 = groups[0]
                if group2 == 'rest' or (not group2 or group2 not in groups):
                    group2_samples = [s for s in adata.obs.index if adata.obs.loc[s, groupby] != group1]
                    group2 = f'rest (n={len(group2_samples)})'
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
                adata.obs[groupby] = pd.Series(
                    ['Group1']*half + ['Group2']*(n-half), index=adata.obs.index
                )

            self.progress(30, f"差异分析: {group1} vs {group2}...")
            deg_df, files, n_up, n_down = _run_single_comparison(
                adata, counts, group1_samples, group2_samples, group1, group2,
                method, fc_threshold, pval_threshold, top_n,
                gene_id_to_name, plots_dir, results_dir,
                viz_params=self.params.get('_visualization', {}))
            result_files.extend(files)

            # 箱线图
            self._draw_boxplots(deg_df, group1, group2, group1_samples, group2_samples,
                                gene_id_to_name, counts, adata, plots_dir, result_files)

            summary = {
                'comparison': f'{group1} vs {group2}',
                'method': method,
                'n_genes_total': len(deg_df),
                'n_up': n_up,
                'n_down': n_down,
                'fc_threshold': fc_threshold,
                'pval_threshold': pval_threshold,
            }

        self.progress(95, "保存 h5ad...")
        adata.obs['group'] = adata.obs[groupby] if groupby in adata.obs.columns else 'unknown'
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

        boxplot_n = min(5, int(self.params.get('top_n', 20)))
        name_to_id = {v: k for k, v in gene_id_to_name.items()} if gene_id_to_name else {}

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
        for pg in plot_gene_list:
            pg_id = name_to_id.get(pg, pg)
            if pg_id in adata.var_names:
                fig_box = go.Figure()
                for grp_name, samples in [(group1, group1_samples), (group2, group2_samples)]:
                    sample_idx = [list(adata.obs.index).index(s) for s in samples if s in adata.obs.index]
                    gene_idx = list(adata.var_names).index(pg_id)
                    vals = counts[sample_idx, gene_idx]
                    fig_box.add_trace(go.Box(y=vals, name=str(grp_name), boxpoints='all', jitter=0.3))
                fig_box.update_layout(title=f'{pg} 表达', yaxis_title='Expression',
                                     plot_bgcolor='white', width=400, height=350)
                fpath = os.path.join(plots_dir, f'bulk_deg_box_{pg}.json')
                with open(fpath, 'w') as f: f.write(fig_box.to_json(engine="json"))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'boxplot', 'label': f'{pg} Boxplot'})
