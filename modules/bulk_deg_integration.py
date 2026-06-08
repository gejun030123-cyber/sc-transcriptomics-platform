# modules/bulk_deg_integration.py
import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkDEGIntegrationAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_deg_integration"
    DISPLAY_NAME = "多组差异整合分析"
    DESCRIPTION = "多组比较结果整合：Upset 图、一致性评分、logFC 矩阵分析"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        self.progress(5, "加载差异分析结果...")
        results_dir = os.path.join(self.project_dir, 'results')
        if not os.path.isdir(results_dir):
            return {'output_adata': input_path, 'result_files': [],
                    'summary': {'error': 'results 目录不存在，请先运行 bulk_deg'}}

        # 扫描 bulk_deg 输出的 DEG CSV 文件
        deg_files = sorted([f for f in os.listdir(results_dir)
                     if f.startswith('bulk_deg_') and f.endswith('.csv')
                     and 'merged' not in f and 'lrt' not in f
                     and 'skipped' not in f and 'gene_filter' not in f
                     and 'all_comparisons' not in f])

        min_comparisons = int(self.params.get('min_comparisons', 2))
        consistency_n = int(self.params.get('consistency_n', 50))
        fc_threshold = float(self.params.get('fc_threshold', 2.0))
        pval_threshold = float(self.params.get('pval_threshold', 0.05))

        # 解析每个比较结果
        self.progress(15, f"解析 {len(deg_files)} 个比较结果...")
        comparisons = {}
        for f in deg_files:
            df = pd.read_csv(os.path.join(results_dir, f))
            if 'gene' not in df.columns or 'log2FC' not in df.columns:
                continue
            name = f.replace('bulk_deg_', '').replace('.csv', '')
            comparisons[name] = df

        if len(comparisons) < 2:
            return {'output_adata': input_path, 'result_files': [],
                    'summary': {'error': f'仅找到 {len(comparisons)} 个比较结果，需要至少 2 个'}}

        # 构建矩阵
        self.progress(25, "构建 logFC/padj 矩阵...")
        all_genes = set()
        for df in comparisons.values():
            all_genes.update(df['gene'].tolist())
        all_genes = sorted(all_genes)
        comp_names = sorted(comparisons.keys())

        logfc_matrix = pd.DataFrame(0.0, index=all_genes, columns=comp_names)
        padj_matrix = pd.DataFrame(1.0, index=all_genes, columns=comp_names)
        regulation_matrix = pd.DataFrame(0.0, index=all_genes, columns=comp_names)

        for name, df in comparisons.items():
            df_idx = df.set_index('gene')
            common = [g for g in all_genes if g in df_idx.index]
            if common:
                logfc_matrix.loc[common, name] = df_idx.loc[common, 'log2FC'].values
                padj_matrix.loc[common, name] = df_idx.loc[common, 'padj'].values
                reg_map = {'Up': 1, 'Down': -1}
                regulation_matrix.loc[common, name] = df_idx.loc[common, 'regulation'].map(reg_map).fillna(0).values

        # 一致性评分
        self.progress(40, "计算一致性评分...")
        log2fc_thresh = np.log2(fc_threshold)
        consistency_scores = []
        for g in all_genes:
            sig_count = int(sum(1 for c in comp_names
                               if padj_matrix.loc[g, c] < pval_threshold
                               and abs(logfc_matrix.loc[g, c]) >= log2fc_thresh))
            if sig_count >= min_comparisons:
                signs = regulation_matrix.loc[g, comp_names].values
                neg_log_p = -np.log10(padj_matrix.loc[g, comp_names].values + 1e-300)
                score = float(np.mean(signs * neg_log_p))
                consistency_scores.append({
                    'gene': g, 'consistency_score': round(score, 4),
                    'n_significant': sig_count,
                    'direction': 'Up' if score > 0 else 'Down' if score < 0 else 'Mixed',
                })
        consistency_df = pd.DataFrame(consistency_scores)
        if len(consistency_df) > 0:
            consistency_df = consistency_df.sort_values('consistency_score', ascending=False, key=abs)

        # 输出目录
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        # 1. 一致性评分表
        self.progress(50, "保存一致性评分...")
        if len(consistency_df) > 0:
            consistency_csv = os.path.join(results_dir, 'deg_integration_consistency.csv')
            consistency_df.to_csv(consistency_csv, index=False)
            result_files.append({'file_path': consistency_csv, 'file_type': 'csv', 'category': 'table', 'label': '一致性评分'})

        # 2. Upset 图（条形图展示交集模式）
        self.progress(55, "生成 Upset 图...")
        sig_sets = {}
        for c in comp_names:
            sig_genes = set(logfc_matrix.index[(padj_matrix[c] < pval_threshold) &
                                                (abs(logfc_matrix[c]) >= log2fc_thresh)])
            sig_sets[c] = sig_genes

        from itertools import combinations
        intersection_data = []
        for r in range(1, len(comp_names) + 1):
            for combo in combinations(comp_names, r):
                isect = sig_sets[combo[0]].copy()
                for c in combo[1:]:
                    isect = isect & sig_sets[c]
                if isect:
                    intersection_data.append({'sets': ' ∩ '.join(combo), 'count': len(isect), 'n_sets': len(combo)})
        intersection_data.sort(key=lambda x: x['count'], reverse=True)

        if intersection_data:
            top_intersections = intersection_data[:20]
            fig_upset = go.Figure(go.Bar(
                x=[d['sets'] for d in top_intersections],
                y=[d['count'] for d in top_intersections],
                marker_color='#1a237e'))
            fig_upset.update_layout(title='差异基因交集模式 (Upset)',
                                    xaxis_title='比较组合', yaxis_title='基因数',
                                    height=400, width=max(600, len(top_intersections)*40+200))
            fpath = os.path.join(plots_dir, 'deg_integration_upset.json')
            with open(fpath, 'w') as f: f.write(fig_upset.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'bar', 'label': 'Upset 交集图'})

        # 3. 方向一致性热图
        self.progress(65, "生成方向一致性热图...")
        if len(consistency_df) > 0:
            top_genes = consistency_df.head(consistency_n)['gene'].tolist()
            reg_subset = regulation_matrix.loc[regulation_matrix.index.isin(top_genes)]
            # 按一致性评分排序
            order = [g for g in top_genes if g in reg_subset.index]
            if order:
                reg_subset = reg_subset.loc[order]
                fig_dir = go.Figure(data=go.Heatmap(
                    z=reg_subset.values.tolist(), x=comp_names, y=reg_subset.index.tolist(),
                    colorscale=[[0, '#1565c0'], [0.5, '#f5f5f5'], [1, '#e53935']],
                    zmid=0, showscale=True, colorbar=dict(title='Direction', tickvals=[-1, 0, 1],
                                                           ticktext=['Down', 'NS', 'Up'])))
                fig_dir.update_layout(title='差异方向一致性矩阵',
                                      height=max(300, len(reg_subset)*12+100),
                                      width=max(400, len(comp_names)*80+200))
                fpath = os.path.join(plots_dir, 'deg_integration_direction_heatmap.json')
                with open(fpath, 'w') as f: f.write(fig_dir.to_json(engine="json"))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '方向一致性矩阵'})

        # 4. logFC 矩阵热图
        self.progress(75, "生成 logFC 矩阵热图...")
        if len(consistency_df) > 0:
            top_genes_fc = consistency_df.head(consistency_n)['gene'].tolist()
            logfc_subset = logfc_matrix.loc[logfc_matrix.index.isin(top_genes_fc)]
            order_fc = [g for g in top_genes_fc if g in logfc_subset.index]
            if order_fc:
                logfc_subset = logfc_subset.loc[order_fc].clip(-5, 5)
                fig_heat = go.Figure(data=go.Heatmap(
                    z=logfc_subset.values.tolist(), x=comp_names, y=logfc_subset.index.tolist(),
                    colorscale='RdBu_r', zmid=0, colorbar=dict(title='log2FC')))
                fig_heat.update_layout(title='Top 差异基因 logFC 矩阵',
                                       height=max(400, len(logfc_subset)*12+100),
                                       width=max(500, len(comp_names)*80+200))
                fpath = os.path.join(plots_dir, 'deg_integration_logfc_heatmap.json')
                with open(fpath, 'w') as f: f.write(fig_heat.to_json(engine="json"))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': 'logFC 矩阵热图'})

        # 5. 比较间 logFC 相关性热图
        self.progress(85, "生成比较间相关性热图...")
        if len(comp_names) >= 2:
            corr_mat = logfc_matrix.corr(method='pearson')
            fig_corr = go.Figure(data=go.Heatmap(
                z=corr_mat.values.tolist(), x=comp_names, y=comp_names,
                colorscale='RdBu_r', zmid=0,
                text=np.round(corr_mat.values, 2).tolist(), texttemplate='%{text}',
                colorbar=dict(title='Pearson r')))
            fig_corr.update_layout(title='比较间 logFC 相关性',
                                   width=max(400, len(comp_names)*80+200),
                                   height=max(400, len(comp_names)*80+200))
            fpath = os.path.join(plots_dir, 'deg_integration_corr.json')
            with open(fpath, 'w') as f: f.write(fig_corr.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '比较间相关性'})

        # 6. Top 一致性基因热图
        self.progress(90, "生成 Top 一致性基因热图...")
        if len(consistency_df) > 0:
            top_n = min(consistency_n, 30)
            top_consistent = consistency_df.head(top_n)
            fig_top = go.Figure(data=go.Bar(
                x=top_consistent['gene'].tolist(),
                y=top_consistent['consistency_score'].tolist(),
                marker_color=['#e53935' if s > 0 else '#1565c0' for s in top_consistent['consistency_score']]))
            fig_top.update_layout(title=f'Top {top_n} 一致性差异基因',
                                  xaxis_title='Gene', yaxis_title='Consistency Score',
                                  height=400, width=max(500, top_n*25+200))
            fpath = os.path.join(plots_dir, 'deg_integration_top_consistent.json')
            with open(fpath, 'w') as f: f.write(fig_top.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'bar', 'label': 'Top 一致性基因'})

        # 输出 h5ad（原样传递）
        self.progress(95, "保存输出...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_deg_integration_output.h5ad')
        import scanpy as sc
        if os.path.exists(input_path):
            adata = sc.read_h5ad(input_path)
            adata.write_h5ad(output_path)
        else:
            output_path = input_path

        n_consistent = len(consistency_df) if len(consistency_df) > 0 else 0
        n_up_consistent = int((consistency_df['direction'] == 'Up').sum()) if n_consistent > 0 else 0
        n_down_consistent = int((consistency_df['direction'] == 'Down').sum()) if n_consistent > 0 else 0

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_comparisons': len(comparisons),
                'comparisons': comp_names,
                'total_genes': len(all_genes),
                'n_consistent_genes': n_consistent,
                'n_up_consistent': n_up_consistent,
                'n_down_consistent': n_down_consistent,
                'min_comparisons_threshold': min_comparisons,
                'fc_threshold': fc_threshold,
                'pval_threshold': pval_threshold,
            }
        }
