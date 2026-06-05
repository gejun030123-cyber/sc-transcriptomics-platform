import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkDEGAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_deg"
    DISPLAY_NAME = "Bulk 差异表达分析"
    DESCRIPTION = "组间差异表达基因检测：t-test / Mann-Whitney / DESeq2 风格"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from scipy import stats

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
        base_mean_filter = float(self.params.get('base_mean_filter', 1))

        self.progress(15, "解析分组信息...")

        # 优先使用 normalized layer（未取 log 的标准化数据），否则用 X
        is_log_transformed = False
        if 'normalized' in adata.layers:
            raw_layer = adata.layers['normalized']
            counts = raw_layer if not hasattr(raw_layer, 'toarray') else raw_layer.toarray()
            counts = counts.astype(float)
        else:
            counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
            counts = counts.astype(float)
            # 检测数据是否已经 log 变换（max < 50 通常是 log2 数据）
            if counts.max() < 50 and counts.min() >= 0:
                is_log_transformed = True
        # 清理 inf/NaN
        counts = np.nan_to_num(counts, nan=0.0, posinf=0.0, neginf=0.0)
        gene_ids = adata.var_names.tolist()
        # 优先使用基因名（symbol），如果没有则使用 gene ID
        if 'gene_name' in adata.var.columns:
            gene_names = adata.var['gene_name'].tolist()
        else:
            gene_names = gene_ids

        # BaseMean 过滤：去除低表达基因
        raw_mean = counts.mean(axis=0)
        expr_mask = raw_mean > base_mean_filter
        counts = counts[:, expr_mask]
        gene_ids = [gene_ids[i] for i in range(len(gene_ids)) if expr_mask[i]]
        gene_names = [gene_names[i] for i in range(len(gene_names)) if expr_mask[i]]

        # 处理自动检测的分组（从样本名中提取）
        if groupby == '_auto_group_':
            auto_mapping = self.params.get('_auto_group_mapping', {})
            # 处理 JSON 字符串情况
            if isinstance(auto_mapping, str):
                try:
                    auto_mapping = json.loads(auto_mapping)
                except Exception:
                    auto_mapping = {}
            if auto_mapping:
                groupby = 'auto_group'
                adata.obs[groupby] = adata.obs.index.map(
                    lambda x: auto_mapping.get(str(x), 'unknown')
                )
            else:
                raise ValueError("自动分组映射数据缺失，请重新选择输入文件。")

        if groupby in adata.obs.columns:
            groups = adata.obs[groupby].unique().tolist()
            if not group1 or group1 not in groups:
                group1 = groups[0]
            mask1 = adata.obs[groupby] == group1
            data1 = counts[mask1.values]

            if group2 == 'rest' or (not group2 or group2 not in groups):
                if group2 == 'rest':
                    # rest 模式：对照组 = 除实验组外的所有样本
                    mask2 = ~mask1
                    group2 = f'rest (n={mask2.sum()})'
                elif len(groups) < 2:
                    raise ValueError(f"分组列 '{groupby}' 中只有 {len(groups)} 个分组，无法进行差异分析。至少需要 2 个分组。")
                else:
                    group2 = groups[1]
                    mask2 = adata.obs[groupby] == group2
            else:
                mask2 = adata.obs[groupby] == group2

            data2 = counts[mask2.values]
        else:
            n = counts.shape[0]
            half = n // 2
            if half == 0 or half == n:
                raise ValueError(f"未找到分组列 '{groupby}'，且样本数 ({n}) 不足以自动分为两组。请确保数据包含分组信息或至少有 2 个样本。")
            group1, group2 = "Group1", "Group2"
            data1, data2 = counts[:half], counts[half:]
            mask1 = pd.Series([True]*half + [False]*(n-half), index=adata.obs.index)
            mask2 = pd.Series([False]*half + [True]*(n-half), index=adata.obs.index)
            adata.obs[groupby] = pd.Series(
                ['Group1']*half + ['Group2']*(n-half), index=adata.obs.index
            )

        self.progress(30, f"差异分析：{method}...")
        n_genes = counts.shape[1]

        if method == 'deseq2':
            import omicverse as ov
            # Build count matrix for pyDEG
            count_df = pd.DataFrame(counts.T, index=gene_ids, columns=adata.obs.index)
            dds_obj = ov.bulk.pyDEG(count_df)
            dds_obj.drop_duplicates_index()
            dds_obj.normalize()
            result_df = dds_obj.deg_analysis(
                list(adata.obs.index[mask1]), list(adata.obs.index[mask2]),
                method='DEseq2', alpha=pval_threshold
            )
            # Extract results
            log2fc = result_df['log2FC'].values if 'log2FC' in result_df.columns else result_df.iloc[:, 1].values
            pvalues = result_df['pvalue'].values if 'pvalue' in result_df.columns else result_df.iloc[:, 3].values
            gene_names = result_df.index.tolist()
            gene_ids = gene_names
            mean1 = np.mean(data1, axis=0)
            mean2 = np.mean(data2, axis=0)
            # Compute padj from pvalues
            from statsmodels.stats.multitest import multipletests
            try:
                _, padj, _, _ = multipletests(pvalues, method='fdr_bh')
            except Exception:
                padj = pvalues
            n_genes = len(gene_ids)
        else:
            log2fc = np.zeros(n_genes)
            pvalues = np.zeros(n_genes)

            mean1 = np.mean(data1, axis=0)
            mean2 = np.mean(data2, axis=0)

            if is_log_transformed:
                # 数据已经是 log2 变换后的，log2FC = mean1 - mean2
                log2fc = mean1 - mean2
            else:
                mean1_safe = np.where(mean1 > 0, mean1, 1e-10)
                mean2_safe = np.where(mean2 > 0, mean2, 1e-10)
                log2fc = np.log2(mean1_safe / mean2_safe)

            for i in range(n_genes):
                if method == 't-test':
                    _, p = stats.ttest_ind(data1[:, i], data2[:, i], equal_var=False)
                else:
                    _, p = stats.mannwhitneyu(data1[:, i], data2[:, i], alternative='two-sided')
                pvalues[i] = p if not np.isnan(p) else 1.0

            from statsmodels.stats.multitest import multipletests
            try:
                _, padj, _, _ = multipletests(pvalues, method='fdr_bh')
            except Exception:
                padj = pvalues

        self.progress(60, "生成结果表...")

        log2fc_threshold = np.log2(fc_threshold)
        regulation = np.where(
            (padj < pval_threshold) & (log2fc >= log2fc_threshold), 'Up',
            np.where(
                (padj < pval_threshold) & (log2fc <= -log2fc_threshold), 'Down',
                'NS'
            )
        )
        regulation = regulation.tolist()
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

        n_up = sum(1 for r in regulation if r == 'Up')
        n_down = sum(1 for r in regulation if r == 'Down')

        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        self.progress(70, "生成火山图...")
        neg_log_padj = -np.log10(padj + 1e-300)
        color_map = {'Up': '#e53935', 'Down': '#1a237e', 'NS': '#bdbdbd'}
        fig_vol = go.Figure()
        for reg in ['NS', 'Up', 'Down']:
            idx = [i for i in range(n_genes) if regulation[i] == reg]
            fig_vol.add_trace(go.Scattergl(
                x=log2fc[idx], y=neg_log_padj[idx], mode='markers',
                marker=dict(color=color_map[reg], size=5, opacity=0.7),
                name=f'{reg} ({len(idx)})',
                text=[gene_names[i] for i in idx], hovertemplate='%{text}<br>log2FC: %{x:.2f}<br>-log10(padj): %{y:.2f}'
            ))
        fig_vol.add_hline(y=-np.log10(pval_threshold), line_dash='dash', line_color='gray')
        fig_vol.add_vline(x=np.log2(fc_threshold), line_dash='dash', line_color='gray')
        fig_vol.add_vline(x=-np.log2(fc_threshold), line_dash='dash', line_color='gray')
        # 火山图基因标注：在图上标注 Top N 差异基因
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
            plot_bgcolor='white', width=700, height=500
        )
        fpath = os.path.join(plots_dir, 'bulk_deg_volcano.json')
        with open(fpath, 'w') as f: f.write(fig_vol.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'volcano', 'label': '火山图'})

        self.progress(80, "生成 MA 图...")
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
            title='MA 图', xaxis_title='log2(Average Expression)', yaxis_title='log2(Fold Change)',
            plot_bgcolor='white', width=700, height=500
        )
        fpath = os.path.join(plots_dir, 'bulk_deg_ma.json')
        with open(fpath, 'w') as f: f.write(fig_ma.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'ma', 'label': 'MA 图'})

        # 基因箱线图：为指定基因绘制分组箱线图
        plot_genes_str = self.params.get('plot_genes', '').strip()
        if plot_genes_str:
            plot_gene_list = [g.strip() for g in plot_genes_str.split(',') if g.strip()]
            for pg in plot_gene_list:
                if pg in adata.var_names:
                    fig_box = go.Figure()
                    for grp_name, grp_mask in [(group1, mask1), (group2, mask2)]:
                        vals = counts[grp_mask.values, list(adata.var_names).index(pg)]
                        fig_box.add_trace(go.Box(y=vals, name=str(grp_name), boxpoints='all', jitter=0.3))
                    fig_box.update_layout(title=f'{pg} 表达', yaxis_title='Expression',
                                         plot_bgcolor='white', width=400, height=350)
                    fpath = os.path.join(plots_dir, f'bulk_deg_box_{pg}.json')
                    with open(fpath, 'w') as f: f.write(fig_box.to_json(engine="json"))
                    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'boxplot', 'label': f'{pg} Boxplot'})

        self.progress(88, "保存差异基因 CSV...")
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        csv_path = os.path.join(results_dir, 'bulk_deg_results.csv')
        deg_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '差异表达基因列表'})

        top_genes = deg_df[deg_df['regulation'] != 'NS'].head(top_n)
        top_csv = os.path.join(results_dir, 'bulk_deg_top_genes.csv')
        top_genes.to_csv(top_csv, index=False)
        result_files.append({'file_path': top_csv, 'file_type': 'csv', 'category': 'table', 'label': f'Top {top_n} 差异基因'})

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
            'summary': {
                'comparison': f'{group1} vs {group2}',
                'method': method,
                'n_genes_total': n_genes,
                'n_up': n_up,
                'n_down': n_down,
                'fc_threshold': fc_threshold,
                'pval_threshold': pval_threshold,
            }
        }
