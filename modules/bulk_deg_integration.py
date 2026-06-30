# modules/bulk_deg_integration.py
import os
import logging
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis

logger = logging.getLogger(__name__)


def _load_comparison_labels(project_dir):
    """Load filename→label mapping from ResultFile table (same logic as deg-comparisons API)."""
    label_map = {}
    try:
        from models import AnalysisTask, ResultFile
        results_dir = os.path.join(project_dir, 'results')
        csv_files = set(f for f in os.listdir(results_dir)
                        if f.startswith('bulk_deg_results') and f.endswith('.csv')
                        and 'merged' not in f and 'all_comparisons' not in f
                        and 'lrt' not in f and 'top_genes' not in f)
        pid = os.path.basename(project_dir)
        for t in AnalysisTask.get_by_project(pid):
            if t.module_name != 'bulk_deg':
                continue
            for rf in ResultFile.get_by_task(t.id):
                fname = os.path.basename(rf.file_path)
                if fname in csv_files and rf.label and 'vs' in rf.label:
                    lbl = rf.label
                    if '(' in lbl and ')' in lbl:
                        lbl = lbl.split('(')[-1].rstrip(')')
                    # Normalize: "moclel vs hmc3" → "moclel-vs-hmc3"
                    lbl = lbl.replace(' vs ', '-vs-').strip()
                    label_map[fname] = lbl
    except Exception as e:
        logger.warning("解析文件标签失败: %s", e)
    return label_map


class BulkDEGIntegrationAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_deg_integration"
    DISPLAY_NAME = "多组差异整合分析"
    DESCRIPTION = "多组比较结果整合：Upset 图、一致性评分、logFC 矩阵分析、表达式筛选"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import plotly.graph_objects as go

        self.progress(5, "加载差异分析结果...")
        results_dir = os.path.join(self.project_dir, 'results')
        if not os.path.isdir(results_dir):
            return {'output_adata': input_path, 'result_files': [],
                    'summary': {'error': 'results 目录不存在，请先运行 bulk_deg'}}

        # 扫描 bulk_deg 输出的 DEG CSV 文件（支持带后缀和无后缀两种格式）
        deg_files = sorted([f for f in os.listdir(results_dir)
                     if f.startswith('bulk_deg_results') and f.endswith('.csv')
                     and 'merged' not in f and 'all_comparisons' not in f
                     and 'lrt' not in f and 'top_genes' not in f])

        # 按用户选择过滤比较文件
        selected = self.params.get('selected_comparisons', '').strip()
        if selected:
            selected_names = set(s.strip() for s in selected.split(',') if s.strip())
            deg_files = [f for f in deg_files
                         if f.replace('bulk_deg_', '').replace('.csv', '') in selected_names]

        min_comparisons = int(self.params.get('min_comparisons', 2))
        # 当比较数少时自动降低阈值，避免交集过窄
        original_min = min_comparisons
        min_comparisons = min(min_comparisons, max(1, len(deg_files) - 1))
        if min_comparisons != original_min:
            self.progress(-1, f"警告: min_comparisons 从 {original_min} 自动调整为 {min_comparisons}（比较数不足）")
        consistency_n = int(self.params.get('consistency_n', 50))
        fc_threshold = float(self.params.get('fc_threshold', 2.0))
        pval_threshold = float(self.params.get('pval_threshold', 0.05))
        upset_top_n = int(self.params.get('upset_top_n', 20))
        logfc_clip = float(self.params.get('logfc_clip_range', 5.0))

        # 加载比较名映射（修复：使用真实比较名而非文件名）
        self.progress(10, "加载比较名映射...")
        label_map = _load_comparison_labels(self.project_dir)

        # 解析每个比较结果
        self.progress(15, f"解析 {len(deg_files)} 个比较结果...")
        comparisons = {}
        name_key_map = {}  # key(文件名派生) → 真实比较名
        for f in deg_files:
            df = pd.read_csv(os.path.join(results_dir, f))
            if 'gene' not in df.columns or 'log2FC' not in df.columns:
                continue
            file_key = f.replace('bulk_deg_', '').replace('.csv', '')
            # 优先使用真实比较名，否则回退到文件名
            real_name = label_map.get(f, file_key)
            comparisons[real_name] = df
            name_key_map[file_key] = real_name

        if len(comparisons) < 2:
            return {'output_adata': input_path, 'result_files': [],
                    'summary': {'error': f'仅找到 {len(comparisons)} 个比较结果，需要至少 2 个'}}

        # 构建矩阵（列名使用真实比较名）
        self.progress(25, "构建 logFC/padj 矩阵...")
        all_genes = set()
        for df in comparisons.values():
            all_genes.update(df['gene'].tolist())
        all_genes = sorted(all_genes)
        comp_names = sorted(comparisons.keys())

        logfc_matrix = pd.DataFrame(0.0, index=all_genes, columns=comp_names)
        padj_matrix = pd.DataFrame(1.0, index=all_genes, columns=comp_names)
        regulation_matrix = pd.DataFrame(0, index=all_genes, columns=comp_names, dtype=int)

        all_genes_idx = pd.Index(all_genes)
        reg_map = {'Up': 1, 'Down': -1}
        for name, df in comparisons.items():
            df_idx = df.set_index('gene')
            common = all_genes_idx.intersection(df_idx.index)
            if len(common) > 0:
                logfc_matrix.loc[common, name] = df_idx.loc[common, 'log2FC']
                padj_matrix.loc[common, name] = df_idx.loc[common, 'padj']
                regulation_matrix.loc[common, name] = df_idx.loc[common, 'regulation'].map(reg_map).fillna(0).values

        # 一致性评分
        self.progress(40, "计算一致性评分...")
        log2fc_thresh = np.log2(fc_threshold)
        consistency_scores = []
        for g in all_genes:
            pvals = padj_matrix.loc[g, comp_names].values
            abs_fcs = np.abs(logfc_matrix.loc[g, comp_names].values)
            sig_mask = (pvals < pval_threshold) & (abs_fcs >= log2fc_thresh)
            sig_count = int(sig_mask.sum())
            if sig_count >= min_comparisons:
                sig_signs = regulation_matrix.loc[g, comp_names].values[sig_mask]
                sig_neg_log_p = -np.log10(pvals[sig_mask] + 1e-300)
                if np.all(sig_signs == sig_signs[0]):
                    score = float(sig_signs[0] * np.mean(sig_neg_log_p))
                else:
                    score = 0.0
                exclude_mixed = self.params.get('exclude_mixed', False)
                if isinstance(exclude_mixed, str):
                    exclude_mixed = exclude_mixed.lower() in ('true', '1', 'yes', 'on')
                if exclude_mixed and score == 0:
                    continue
                consistency_scores.append({
                    'gene': g, 'consistency_score': round(score, 4),
                    'n_significant': sig_count,
                    'direction': 'Up' if score > 0 else 'Down' if score < 0 else 'Mixed',
                })
        consistency_df = pd.DataFrame(consistency_scores)
        if len(consistency_df) > 0:
            consistency_df = consistency_df.sort_values('consistency_score', ascending=False, key=abs)
        else:
            self.progress(-1, "警告: 无一致性基因通过阈值，请尝试放宽 fc_threshold 或 pval_threshold")

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
        self.progress(60, "生成 Upset 图...")
        sig_sets = {}
        for c in comp_names:
            sig_genes = set(logfc_matrix.index[(padj_matrix[c] < pval_threshold) &
                                                (abs(logfc_matrix[c]) >= log2fc_thresh)])
            sig_sets[c] = sig_genes

        from itertools import combinations
        intersection_data = []
        all_comp_set = set(comp_names)
        for r in range(1, len(comp_names) + 1):
            for combo in combinations(comp_names, r):
                combo_set = set(combo)
                isect = sig_sets[combo[0]].copy()
                for c in combo[1:]:
                    isect &= sig_sets[c]
                # Upset 严格模式：排除在其他比较中也显著的基因
                upset_strict = self.params.get('upset_strict', True)
                if isinstance(upset_strict, str):
                    upset_strict = upset_strict.lower() in ('true', '1', 'yes', 'on')
                if upset_strict:
                    for o in all_comp_set - combo_set:
                        isect -= sig_sets[o]
                if isect:
                    intersection_data.append({'sets': ' ∩ '.join(combo), 'count': len(isect), 'n_sets': len(combo)})
        intersection_data.sort(key=lambda x: x['count'], reverse=True)

        if intersection_data:
            top_intersections = intersection_data[:upset_top_n]
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

        # 比较差异基因数统计表
        count_rows = []
        for c in comp_names:
            n_up = int(((padj_matrix[c] < pval_threshold) & (logfc_matrix[c] >= log2fc_thresh)).sum())
            n_down = int(((padj_matrix[c] < pval_threshold) & (logfc_matrix[c] <= -log2fc_thresh)).sum())
            count_rows.append({'comparison': c, 'n_up': n_up, 'n_down': n_down, 'n_total': n_up + n_down})
        count_df = pd.DataFrame(count_rows)
        count_csv = os.path.join(results_dir, 'deg_integration_comparison_counts.csv')
        count_df.to_csv(count_csv, index=False)
        result_files.append({'file_path': count_csv, 'file_type': 'csv', 'category': 'table', 'label': '各比较差异基因数统计'})

        # Jaccard 相似度矩阵
        if len(comp_names) >= 2:
            jaccard_matrix = pd.DataFrame(1.0, index=comp_names, columns=comp_names)
            for i, c1 in enumerate(comp_names):
                for j, c2 in enumerate(comp_names):
                    if i < j:
                        intersection = len(sig_sets[c1] & sig_sets[c2])
                        union = len(sig_sets[c1] | sig_sets[c2])
                        jval = intersection / union if union > 0 else 0.0
                        jaccard_matrix.loc[c1, c2] = jval
                        jaccard_matrix.loc[c2, c1] = jval
            fig_jaccard = go.Figure(data=go.Heatmap(
                z=jaccard_matrix.values.tolist(), x=comp_names, y=comp_names,
                colorscale='Blues', zmin=0, zmax=1,
                text=np.round(jaccard_matrix.values, 3).tolist(), texttemplate='%{text}',
                colorbar=dict(title='Jaccard')))
            fig_jaccard.update_layout(
                title='差异基因 Jaccard 相似度',
                width=max(400, len(comp_names)*80+200),
                height=max(400, len(comp_names)*80+200))
            fpath = os.path.join(plots_dir, 'deg_integration_jaccard.json')
            with open(fpath, 'w') as f: f.write(fig_jaccard.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': 'Jaccard 相似度'})

        # 2b. Venn 图（2 或 3 个比较时生成）
        self.progress(62, "生成 Venn 图...")
        if len(comp_names) in (2, 3):
            only_sets = {}
            for c in comp_names:
                only = sig_sets[c].copy()
                for o in comp_names:
                    if o != c:
                        only -= sig_sets[o]
                only_sets[c] = only

            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            try:
                from matplotlib_venn import venn2, venn3
                has_venn_lib = True
            except ImportError:
                has_venn_lib = False

            if has_venn_lib:
                fig_venn, ax = plt.subplots(1, 1, figsize=(6, 5))
                if len(comp_names) == 2:
                    c0, c1 = comp_names
                    subsets = (
                        len(only_sets[c0]),
                        len(only_sets[c1]),
                        len(sig_sets[c0] & sig_sets[c1]),
                    )
                    v = venn2(subsets, set_labels=comp_names, ax=ax)
                    colors = ['#e53935', '#1565c0']
                    for i, patch_id in enumerate(['10', '01']):
                        patch = v.get_patch_by_id(patch_id)
                        if patch:
                            patch.set_color(colors[i])
                            patch.set_alpha(0.35)
                else:
                    c0, c1, c2 = comp_names
                    subsets = (
                        len(only_sets[c0]),
                        len(only_sets[c1]),
                        len(sig_sets[c0] & sig_sets[c1]),
                        len(only_sets[c2]),
                        len(sig_sets[c0] & sig_sets[c2]),
                        len(sig_sets[c1] & sig_sets[c2]),
                        len(sig_sets[c0] & sig_sets[c1] & sig_sets[c2]),
                    )
                    v = venn3(subsets, set_labels=comp_names, ax=ax)
                    colors = ['#e53935', '#1565c0', '#4caf50']
                    for i, patch_id in enumerate(['100', '010', '001']):
                        patch = v.get_patch_by_id(patch_id)
                        if patch:
                            patch.set_color(colors[i])
                            patch.set_alpha(0.35)

                ax.set_title(f'Venn: {" vs ".join(comp_names)}', fontsize=13)
                plt.tight_layout()
                fpath = os.path.join(plots_dir, 'deg_integration_venn.png')
                fig_venn.savefig(fpath, dpi=150, bbox_inches='tight', facecolor='white')
                plt.close(fig_venn)
                result_files.append({'file_path': fpath, 'file_type': 'png', 'category': 'venn', 'label': 'Venn 图'})
            else:
                self.progress(-1, "警告: matplotlib_venn 未安装，跳过 Venn 图生成")

        # 3. 方向一致性热图
        self.progress(70, "生成方向一致性热图...")
        if len(consistency_df) > 0:
            top_genes = consistency_df.head(consistency_n)['gene'].tolist()
            reg_subset = regulation_matrix.reindex(top_genes).dropna()
            if len(reg_subset) > 0:
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
        self.progress(80, "生成 logFC 矩阵热图...")
        if len(consistency_df) > 0:
            top_genes_fc = consistency_df.head(consistency_n)['gene'].tolist()
            logfc_subset = logfc_matrix.reindex(top_genes_fc).dropna()
            if len(logfc_subset) > 0:
                logfc_subset = logfc_subset.clip(-logfc_clip, logfc_clip)
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
        self.progress(87, "生成比较间相关性热图...")
        if len(comp_names) >= 2:
            sig_in_any = ((padj_matrix < pval_threshold) & (abs(logfc_matrix) >= log2fc_thresh)).any(axis=1)
            sig_logfc = logfc_matrix.loc[sig_in_any]
            corr_mat = sig_logfc.corr(method='pearson') if len(sig_logfc) > 0 else logfc_matrix.corr(method='pearson')
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

        # ===== 表达式筛选器 =====
        filter_expression = self.params.get('filter_expression', '').strip()
        if filter_expression:
            result_files = self._run_filter(
                filter_expression, comp_names, logfc_matrix, padj_matrix,
                pval_threshold, log2fc_thresh, plots_dir, results_dir, result_files)

        # 输出 h5ad（原样传递）
        self.progress(95, "保存输出...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_deg_integration_output.h5ad')
        import scanpy as sc
        if os.path.exists(input_path):
            adata = sc.read_h5ad(input_path)
        else:
            adata = sc.AnnData()
        adata.write_h5ad(output_path)

        n_consistent = len(consistency_df) if len(consistency_df) > 0 else 0
        n_up_consistent = int((consistency_df['direction'] == 'Up').sum()) if n_consistent > 0 else 0
        n_down_consistent = int((consistency_df['direction'] == 'Down').sum()) if n_consistent > 0 else 0

        summary = {
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

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }

    def _run_filter(self, filter_expression, comp_names, logfc_matrix, padj_matrix,
                    pval_threshold, log2fc_thresh, plots_dir, results_dir, result_files):
        """Execute expression filter and generate result files."""
        from modules.expression_parser import validate, evaluate, ParseError

        self.progress(91, "解析筛选表达式...")
        ast, err, warnings = validate(filter_expression, comp_names)
        if err:
            self.progress(-1, f"表达式错误: {err}")
            return result_files
        for w in warnings:
            self.progress(-1, f"警告: {w}")

        # Build gene index for evaluation
        gene_sets = {}
        for c in comp_names:
            sig_mask = (padj_matrix[c] < pval_threshold) & (abs(logfc_matrix[c]) >= log2fc_thresh)
            gene_sets[c] = set(logfc_matrix.index[sig_mask])

        try:
            filtered_genes = evaluate(ast, comp_names, gene_sets,
                                      padj_matrix, logfc_matrix,
                                      pval_threshold, log2fc_thresh)
        except ParseError as e:
            self.progress(-1, f"表达式求值错误: {e}")
            return result_files

        if not filtered_genes:
            self.progress(-1, "筛选结果为空，无基因满足条件")
            return result_files

        filtered_genes = sorted(filtered_genes)
        self.progress(92, f"筛选完成：{len(filtered_genes)} 个基因")

        # 1. Filter result CSV
        rows = []
        for g in filtered_genes:
            row = {'gene': g}
            satisfies = []
            for c in comp_names:
                p = padj_matrix.loc[g, c]
                fc = logfc_matrix.loc[g, c]
                abs_fc = abs(fc)
                sig = p < pval_threshold and abs_fc >= log2fc_thresh
                direction = 'up' if sig and fc > 0 else 'down' if sig and fc < 0 else 'ns'
                row[f'log2FC_{c}'] = round(float(fc), 4)
                row[f'padj_{c}'] = round(float(p), 6)
                row[f'direction_{c}'] = direction
                if sig:
                    satisfies.append(f"{c}:{direction}")
            row['satisfies'] = '; '.join(satisfies)
            rows.append(row)

        filter_csv = os.path.join(results_dir, 'deg_filter_results.csv')
        pd.DataFrame(rows).to_csv(filter_csv, index=False)
        result_files.append({
            'file_path': filter_csv, 'file_type': 'csv',
            'category': 'table', 'label': f'筛选结果 ({len(filtered_genes)} 基因)'
        })

        # 5. 基因列表纯文本导出
        gene_list_path = os.path.join(results_dir, 'deg_filter_gene_list.txt')
        with open(gene_list_path, 'w') as f:
            f.write('\n'.join(filtered_genes))
        result_files.append({
            'file_path': gene_list_path, 'file_type': 'txt',
            'category': 'table', 'label': '筛选基因列表 (txt)'
        })

        # 2. Filter UpSet plot (which atoms the genes satisfy)
        self.progress(93, "生成筛选 UpSet 图...")
        import plotly.graph_objects as go
        from itertools import combinations as iter_combos

        atom_sets = {}
        for c in comp_names:
            atom_sets[c] = set(g for g in filtered_genes
                               if padj_matrix.loc[g, c] < pval_threshold
                               and abs(logfc_matrix.loc[g, c]) >= log2fc_thresh)

        filter_upset_data = []
        all_c = set(comp_names)
        for r in range(1, len(comp_names) + 1):
            for combo in iter_combos(comp_names, r):
                combo_set = set(combo)
                isect = set(filtered_genes) if r == len(comp_names) else atom_sets[combo[0]].copy()
                if r < len(comp_names):
                    for c in combo[1:]:
                        isect &= atom_sets[c]
                isect = {g for g in isect if all(g in atom_sets[c] for c in combo)}
                for o in all_c - combo_set:
                    isect -= atom_sets[o]
                if isect:
                    filter_upset_data.append({
                        'sets': ' ∩ '.join(combo), 'count': len(isect), 'n_sets': len(combo)
                    })

        filter_upset_data.sort(key=lambda x: x['count'], reverse=True)
        if filter_upset_data:
            top20 = filter_upset_data[:20]
            fig_fu = go.Figure(go.Bar(
                x=[d['sets'] for d in top20],
                y=[d['count'] for d in top20],
                marker_color='#1a237e'))
            fig_fu.update_layout(title='筛选基因交集模式',
                                 xaxis_title='比较组合', yaxis_title='基因数',
                                 height=400, width=max(600, len(top20)*40+200))
            fpath = os.path.join(plots_dir, 'deg_filter_upset.json')
            with open(fpath, 'w') as f: f.write(fig_fu.to_json(engine="json"))
            result_files.append({
                'file_path': fpath, 'file_type': 'plotly_json',
                'category': 'bar', 'label': '筛选基因 Upset 图'
            })

        # 3. Filter logFC heatmap
        self.progress(94, "生成筛选基因 logFC 热图...")
        filter_show_n = int(self.params.get('filter_show_n', 80))
        show_n = min(len(filtered_genes), filter_show_n)
        show_genes = filtered_genes[:show_n]
        logfc_clip = float(self.params.get('logfc_clip_range', 5.0))
        logfc_sub = logfc_matrix.reindex(show_genes).clip(-logfc_clip, logfc_clip)
        fig_fheat = go.Figure(data=go.Heatmap(
            z=logfc_sub.values.tolist(), x=comp_names, y=show_genes,
            colorscale='RdBu_r', zmid=0, colorbar=dict(title='log2FC')))
        fig_fheat.update_layout(title='筛选基因 logFC 矩阵',
                                height=max(400, len(show_genes)*14+100),
                                width=max(500, len(comp_names)*80+200))
        fpath = os.path.join(plots_dir, 'deg_filter_logfc_heatmap.json')
        with open(fpath, 'w') as f: f.write(fig_fheat.to_json(engine="json"))
        result_files.append({
            'file_path': fpath, 'file_type': 'plotly_json',
            'category': 'heatmap', 'label': '筛选基因 logFC 热图'
        })

        # 4. Expression parse tree visualization
        self.progress(95, "生成表达式解析树...")
        tree_data = ast.to_dict()
        import json
        fpath = os.path.join(plots_dir, 'deg_filter_expression_tree.json')
        with open(fpath, 'w') as f:
            json.dump(tree_data, f, ensure_ascii=False)
        result_files.append({
            'file_path': fpath, 'file_type': 'json',
            'category': 'tree', 'label': '表达式解析树'
        })

        return result_files
