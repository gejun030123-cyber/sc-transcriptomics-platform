import os
from modules.base import BaseAnalysis

DEFAULT_TME_MARKERS = {
    'Tumor Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
    'CAF': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'FAP', 'ACTA2'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'ENG', 'CLDN5'],
    'T cells': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'TRAC'],
    'NK cells': ['NKG7', 'GNLY', 'KLRD1', 'NCAM1', 'PRF1'],
    'B cells': ['CD79A', 'CD79B', 'MS4A1', 'CD19', 'PAX5'],
    'Plasma cells': ['JCHAIN', 'MZB1', 'SDC1', 'IGHG1', 'IGKC'],
    'Monocyte/Macrophage': ['CD14', 'CD68', 'CSF1R', 'LYZ', 'S100A8', 'S100A9'],
    'Dendritic cells': ['FCER1A', 'CD1C', 'CLEC10A', 'ITGAX', 'HLA-DRA'],
    'Neutrophils': ['CSF3R', 'CXCR2', 'FCGR3B', 'S100A12'],
    'Mast cells': ['KIT', 'TPSAB1', 'TPSB2', 'HDC', 'MS4A2'],
    'Pericytes': ['RGS5', 'PDGFRB', 'NOTCH3', 'MCAM', 'ACTA2'],
    'Proliferating': ['MKI67', 'TOP2A', 'PCNA', 'STMN1', 'CDK1'],
}

DEFAULT_IMMUNE_MARKERS = {
    'T cells': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'TRAC'],
    'CD4+ T': ['CD4', 'IL7R', 'TRBC2'],
    'CD8+ T': ['CD8A', 'CD8B', 'GZMK', 'GZMA', 'CCL5'],
    'T naive': ['LEF1', 'CCR7', 'TCF7'],
    'NK cells': ['NKG7', 'GNLY', 'KLRD1', 'NCAM1', 'PRF1'],
    'B cells': ['CD79A', 'CD79B', 'MS4A1', 'CD19', 'PAX5'],
    'Plasma cells': ['JCHAIN', 'MZB1', 'SDC1', 'IGHG1', 'IGKC'],
    'Monocyte/Macrophage': ['CD14', 'CD68', 'CSF1R', 'LYZ', 'S100A8', 'S100A9'],
    'Dendritic cells': ['FCER1A', 'CD1C', 'CLEC10A', 'ITGAX', 'HLA-DRA'],
    'Neutrophils': ['CSF3R', 'CXCR2', 'FCGR3B', 'S100A12'],
    'pDC': ['GZMB', 'IL3RA', 'COBLL1', 'TCF4'],
}

DEFAULT_BLOOD_MARKERS = {
    'HSC': ['CD34', 'CD38', 'KIT', 'THY1', 'CRHBP'],
    'Erythroid': ['HBA1', 'HBA2', 'HBB', 'GYPA', 'SLC4A1'],
    'Megakaryocyte': ['PF4', 'GP9', 'ITGA2B', 'VWF', 'GP1BA'],
    'Monocyte': ['CD14', 'LYZ', 'S100A8', 'S100A9', 'VCAN'],
    'Neutrophil': ['FCGR3B', 'CSF3R', 'CXCR2', 'S100A12', 'MPO'],
    'Eosinophil': ['SIGLEC8', 'IL5RA', 'CCR3', 'EPX', 'PRG2'],
    'Basophil': ['HDC', 'MS4A2', 'KIT', 'FCER1A', 'CPA3'],
    'B cell': ['CD79A', 'MS4A1', 'CD19', 'PAX5', 'CD79B'],
    'T cell': ['CD3D', 'CD3E', 'CD2', 'TRAC', 'CD3G'],
    'NK cell': ['NKG7', 'GNLY', 'KLRD1', 'NCAM1', 'PRF1'],
    'Dendritic cell': ['FCER1A', 'CD1C', 'CLEC10A', 'ITGAX', 'HLA-DRA'],
    'pDC': ['GZMB', 'IL3RA', 'COBLL1', 'TCF4', 'IRF7'],
}

DEFAULT_PBMC_MARKERS = {
    'CD4 Naive T cells': ['IL7R', 'LTB', 'CCR7', 'TCF7', 'MALAT1', 'LEF1'],
    'CD4 Memory T cells': ['IL7R', 'LTB', 'MALAT1', 'IL32', 'AQP3', 'GPR183'],
    'CD14+ Monocytes': ['CD14', 'LYZ', 'S100A8', 'S100A9', 'LGALS3', 'FCN1'],
    'B cells': ['MS4A1', 'CD79A', 'CD79B', 'CD74', 'HLA-DRA'],
    'CD8 T cells': ['CD8A', 'CD8B', 'CCL5', 'GZMK', 'GZMA', 'TRBC2'],
    'Cytotoxic T cells': ['NKG7', 'CCL5', 'GZMB', 'PRF1', 'CTSW', 'CD3D'],
    'NK cells': ['GNLY', 'NKG7', 'KLRD1', 'PRF1', 'CTSW'],
    'FCGR3A+ Monocytes': ['FCGR3A', 'MS4A7', 'LST1', 'FCER1G', 'AIF1', 'IFITM3'],
    'Conventional DC': ['FCER1A', 'CST3', 'CD1C', 'CLEC10A', 'HLA-DPA1', 'HLA-DPB1'],
    'Plasmacytoid DC': ['GZMB', 'IRF7', 'TCF4', 'IL3RA', 'SERPINF1'],
    'Megakaryocytes': ['PPBP', 'PF4', 'SDPR', 'GNG11', 'NRGN'],
}

# A deliberately broad first-pass panel. It is suitable when tissue type is
# unknown; users should use a tissue/immune panel or custom markers afterwards
# to refine a lineage, rather than treating these labels as final subtypes.
DEFAULT_UNIVERSAL_MARKERS = {
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'KRT7', 'CDH1'],
    'Endothelial': ['PECAM1', 'VWF', 'KDR', 'EMCN', 'CLDN5'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'COL3A1'],
    'Pericyte/Smooth muscle': ['RGS5', 'PDGFRB', 'CSPG4', 'MCAM', 'ACTA2'],
    'Myeloid': ['LYZ', 'TYROBP', 'LST1', 'FCER1G', 'AIF1'],
    'T cells': ['CD3D', 'CD3E', 'TRAC', 'CD247', 'LCK'],
    'NK cells': ['NKG7', 'KLRD1', 'GNLY', 'PRF1', 'TRBC2'],
    'B cells': ['MS4A1', 'CD79A', 'CD74', 'HLA-DRA', 'CD37'],
    'Plasma cells': ['JCHAIN', 'MZB1', 'SDC1', 'DERL3', 'IGKC'],
    'Mast cells': ['TPSAB1', 'TPSB2', 'KIT', 'MS4A2', 'HDC'],
    'Cycling cells': ['MKI67', 'TOP2A', 'STMN1', 'TYMS', 'CDK1'],
}

UNIVERSAL_LABELS = {
    'Epithelial': ('Non-immune cell', 'Epithelial cell', 'Epithelial cell', 'CL:0000066'),
    'Endothelial': ('Non-immune cell', 'Endothelial cell', 'Endothelial cell', 'CL:0000115'),
    'Fibroblast': ('Stromal cell', 'Fibroblast', 'Fibroblast', 'CL:0000057'),
    'Pericyte/Smooth muscle': ('Stromal cell', 'Perivascular cell', 'Pericyte/smooth muscle cell', 'CL:0000669'),
    'Myeloid': ('Immune cell', 'Myeloid cell', 'Myeloid cell', 'CL:0000763'),
    'T cells': ('Immune cell', 'T cell', 'T cell', 'CL:0000084'),
    'NK cells': ('Immune cell', 'NK cell', 'Natural killer cell', 'CL:0000623'),
    'B cells': ('Immune cell', 'B cell', 'B cell', 'CL:0000236'),
    'Plasma cells': ('Immune cell', 'B cell', 'Plasma cell', 'CL:0000786'),
    'Mast cells': ('Immune cell', 'Mast cell', 'Mast cell', 'CL:0000097'),
    'Cycling cells': ('Cell state', 'Cycling cell', 'Cycling cell', ''),
}

MARKER_SETS = {
    'Universal': DEFAULT_UNIVERSAL_MARKERS,
    'TME': DEFAULT_TME_MARKERS,
    'Immune': DEFAULT_IMMUNE_MARKERS,
    'Blood': DEFAULT_BLOOD_MARKERS,
    'PBMC': DEFAULT_PBMC_MARKERS,
}


def resolve_marker_genes(markers, var_names, min_markers_per_type=2):
    """Match marker symbols case-insensitively and report usable coverage.

    This supports common human/mouse symbol casing differences and prevents a
    lineage with almost no measurable markers from participating in scoring.
    """
    lookup = {str(gene).upper(): str(gene) for gene in var_names}
    usable, coverage = {}, {}
    for cell_type, genes in markers.items():
        matched = list(dict.fromkeys(lookup[g.upper()] for g in genes if g.upper() in lookup))
        coverage[cell_type] = {'matched': len(matched), 'total': len(genes)}
        if len(matched) >= min_markers_per_type:
            usable[cell_type] = matched
    return usable, coverage

class AnnotationAnalysis(BaseAnalysis):
    MODULE_NAME = "annotation"
    DISPLAY_NAME = "细胞注释"
    DESCRIPTION = "基于 Marker 基因的细胞类型自动注释"
    INPUT_REQUIRES = ['leiden']

    def run(self, input_path):
        import scanpy as sc
        from modules.visualization import umap_scatter
        import json
        import numpy as np
        import plotly.graph_objects as go

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        cluster_key = self.params.get('cluster_key', 'leiden')
        resolution = self.params.get('resolution', '0.8')
        leiden_key = f'leiden_{resolution}' if f'leiden_{resolution}' in adata.obs.columns else cluster_key

        self.progress(20, "Scoring cell type markers...")
        requested_method = self.params.get('method', 'auto_marker')
        method = requested_method
        multi_evidence = method == 'multi_evidence'
        if multi_evidence:
            method = 'auto_marker'
        marker_set_name = self.params.get('marker_set', 'Universal')
        custom_markers_str = self.params.get('custom_markers', '').strip()
        confidence_method = self.params.get('confidence_method', 'none')
        mark_unknown = self.params.get('mark_unknown', True)
        merge_similar_threshold = float(self.params.get('merge_similar_threshold', 0))

        if method == 'manual' and custom_markers_str:
            self.progress(40, "Applying manual cell type mapping...")
            manual_mapping = {}
            for line in custom_markers_str.split('\n'):
                line = line.strip()
                if not line or ':' not in line:
                    continue
                cluster_id, cell_type = line.split(':', 1)
                manual_mapping[cluster_id.strip()] = cell_type.strip()

            if manual_mapping:
                adata.obs['celltype'] = adata.obs[leiden_key].astype(str).map(manual_mapping)
                adata.obs['celltype'] = adata.obs['celltype'].fillna('Unknown').astype('category')
                markers = {}
            else:
                self.progress(45, "No valid mapping found, falling back to auto_marker...")
                method = 'auto_marker'

        if method == 'celltypist':
            self.progress(46, "Running CellTypist annotation...")
            try:
                import celltypist
                from celltypist import annotate
                ct_model = self.params.get('celltypist_model', 'Immune_All_Low')
                ct_threshold = float(self.params.get('celltypist_threshold', 0.5))
                majority_voting = self.params.get('celltypist_majority_voting', True)

                predictions = annotate(
                    adata, model=ct_model,
                    majority_voting=majority_voting,
                    confidence_threshold=ct_threshold,
                )
                adata.obs['celltype'] = predictions.predicted_labels['predicted_labels'].values
                adata.obs['celltype'] = adata.obs['celltype'].astype('category')
                markers = {}
            except ImportError:
                self.progress(47, "celltypist not installed, falling back to auto_marker...")
                method = 'auto_marker'
            except Exception as e:
                self.progress(47, f"CellTypist failed: {e}, falling back to auto_marker...")
                method = 'auto_marker'

        if method == 'auto_marker':
            if custom_markers_str:
                # Parse custom markers: "CellType1:GENE1,GENE2;CellType2:GENE3,GENE4"
                # Also supports newline-separated format
                markers = {}
                for ct_genes in custom_markers_str.replace('\n', ';').split(';'):
                    ct_genes = ct_genes.strip()
                    if ':' in ct_genes:
                        ct, genes_str = ct_genes.split(':', 1)
                        markers[ct.strip()] = [g.strip() for g in genes_str.split(',') if g.strip()]
                if not markers:
                    markers = MARKER_SETS.get(marker_set_name, DEFAULT_UNIVERSAL_MARKERS)
            else:
                markers = MARKER_SETS.get(marker_set_name, DEFAULT_UNIVERSAL_MARKERS)

            min_markers = int(self.params.get('min_markers_per_type', 2))
            usable_markers, marker_coverage = resolve_marker_genes(markers, adata.var_names, min_markers)
            skipped_types = [ct for ct in markers if ct not in usable_markers]
            if skipped_types:
                self.progress(-1, f"{len(skipped_types)} 个类型的可用 marker 少于 {min_markers}，未参与自动判定。")
            for ct, available_genes in usable_markers.items():
                try:
                    sc.tl.score_genes(adata, available_genes, score_name=f'score_{ct}', use_raw=False)
                except RuntimeError as exc:
                    # Tiny targeted panels can lack Scanpy control genes. Keep
                    # the annotation usable with an explicit mean-expression
                    # fallback rather than silently omitting the lineage.
                    expr = adata[:, available_genes].X
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray()
                    adata.obs[f'score_{ct}'] = np.asarray(expr, dtype=float).mean(axis=1)
                    self.progress(-1, f'{ct} 使用平均 Marker 表达评分（{exc}）')

            self.progress(50, "Assigning cell types to clusters...")
            score_cols = [f'score_{ct}' for ct in markers if f'score_{ct}' in adata.obs.columns]
            if score_cols:
                cluster_annotations = {}
                for cluster in adata.obs[leiden_key].cat.categories:
                    mask = adata.obs[leiden_key] == cluster
                    mean_scores = {col: adata.obs.loc[mask, col].mean() for col in score_cols}
                    best_col = max(mean_scores, key=mean_scores.get)
                    best_ct = best_col.replace('score_', '')
                    cluster_annotations[cluster] = best_ct
                adata.obs['celltype'] = adata.obs[leiden_key].map(cluster_annotations).astype('category')
                min_score = float(self.params.get('min_annotation_score', 0.0))
                if min_score > 0:
                    best_scores = adata.obs[score_cols].max(axis=1)
                    if 'Unknown' not in adata.obs['celltype'].cat.categories:
                        adata.obs['celltype'] = adata.obs['celltype'].cat.add_categories(['Unknown'])
                    adata.obs.loc[best_scores < min_score, 'celltype'] = 'Unknown'
            else:
                adata.obs['celltype'] = adata.obs[leiden_key].astype(str)
            markers = usable_markers
        else:
            marker_coverage = {}

        # Multi-evidence mode keeps the rule engine interpretable, then uses
        # CellTypist as an independent human reference classifier when it is
        # available. Missing models/network never erase the marker result.
        if multi_evidence:
            adata.obs['marker_label'] = adata.obs['celltype'].astype(str)
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if score_cols:
                score_frame = adata.obs[score_cols]
                adata.obs['marker_score'] = score_frame.max(axis=1).astype(float)
                ranked = np.argsort(score_frame.values, axis=1)
                adata.obs['prediction_1'] = [score_cols[i].replace('score_', '') for i in ranked[:, -1]]
                adata.obs['prediction_1_score'] = score_frame.max(axis=1).astype(float)
                adata.obs['prediction_2'] = [score_cols[i].replace('score_', '') for i in ranked[:, -2]] if len(score_cols) > 1 else ''
                if len(score_cols) > 1:
                    adata.obs['prediction_2_score'] = score_frame.values[np.arange(adata.n_obs), ranked[:, -2]]
            try:
                from celltypist import annotate
                ct_model = self.params.get('celltypist_model', 'Immune_All_Low')
                prediction = annotate(adata, model=ct_model,
                                      majority_voting=self.params.get('celltypist_majority_voting', True))
                labels = prediction.predicted_labels
                adata.obs['celltypist_label'] = labels.iloc[:, 0].astype(str).values
                confidence = getattr(prediction, 'probability_matrix', None)
                if confidence is not None:
                    adata.obs['celltypist_score'] = confidence.max(axis=1).values
            except Exception as exc:
                self.progress(-1, f'CellTypist 不可用，保留 Marker 规则结果：{exc}')

            cluster_vote = adata.obs.groupby(leiden_key, observed=True)['marker_label'].transform(
                lambda values: values.value_counts(normalize=True).iloc[0]
            )
            adata.obs['cluster_annotation_agreement'] = cluster_vote.astype(float)
            adata.obs['annotation_status'] = np.where(
                (adata.obs.get('marker_score', 0) < float(self.params.get('min_annotation_score', 0.0))) |
                (adata.obs['cluster_annotation_agreement'] < float(self.params.get('cluster_agreement_threshold', 0.6))),
                'Unknown', 'review'
            )
            adata.obs['final_annotation'] = adata.obs['marker_label'].astype(str)
            adata.obs.loc[adata.obs['annotation_status'] == 'Unknown', 'final_annotation'] = 'Unknown'
            adata.obs['celltype'] = adata.obs['final_annotation'].astype('category')

        if marker_set_name == 'Universal' and 'celltype' in adata.obs:
            meta = adata.obs['celltype'].astype(str).map(UNIVERSAL_LABELS)
            adata.obs['cell_type_l1'] = meta.map(lambda value: value[0] if isinstance(value, tuple) else 'Unknown')
            adata.obs['cell_type_l2'] = meta.map(lambda value: value[1] if isinstance(value, tuple) else 'Unknown')
            adata.obs['cell_type_l3'] = meta.map(lambda value: value[2] if isinstance(value, tuple) else 'Unknown')
            adata.obs['cell_ontology_id'] = meta.map(lambda value: value[3] if isinstance(value, tuple) else '')
            if 'final_annotation' not in adata.obs:
                adata.obs['final_annotation'] = adata.obs['celltype'].astype(str)
            if 'annotation_status' not in adata.obs:
                adata.obs['annotation_status'] = np.where(adata.obs['celltype'].astype(str) == 'Unknown', 'Unknown', 'review')

        # 置信度评估
        confidence_col = None
        if confidence_method != 'none' and 'celltype' in adata.obs.columns:
            self.progress(60, f"Computing annotation confidence ({confidence_method})...")
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if confidence_method == 'entropy' and score_cols:
                import numpy as np
                score_matrix = adata.obs[score_cols].values
                exp_scores = np.exp(score_matrix - score_matrix.max(axis=1, keepdims=True))
                probs = exp_scores / (exp_scores.sum(axis=1, keepdims=True) + 1e-10)
                entropy = -np.sum(probs * np.log(probs + 1e-10), axis=1)
                max_entropy = np.log(len(score_cols)) if len(score_cols) > 1 else 1
                adata.obs['annotation_confidence'] = 1 - entropy / (max_entropy + 1e-10)
                confidence_col = 'annotation_confidence'
            elif confidence_method == 'score_margin' and score_cols:
                import numpy as np
                score_matrix = adata.obs[score_cols].values
                sorted_scores = np.sort(score_matrix, axis=1)
                if sorted_scores.shape[1] >= 2:
                    adata.obs['annotation_score_margin'] = sorted_scores[:, -1] - sorted_scores[:, -2]
                else:
                    adata.obs['annotation_score_margin'] = sorted_scores[:, -1]
                confidence_col = 'annotation_score_margin'

            if mark_unknown and confidence_col in adata.obs.columns:
                low_conf_mask = adata.obs[confidence_col] < 0.2
                if hasattr(adata.obs['celltype'], 'cat') and 'Unknown' not in adata.obs['celltype'].cat.categories:
                    adata.obs['celltype'] = adata.obs['celltype'].cat.add_categories(['Unknown'])
                adata.obs.loc[low_conf_mask, 'celltype'] = 'Unknown'
                n_unknown = low_conf_mask.sum()
                if n_unknown > 0:
                    self.progress(-1, f"标记 {n_unknown} 个低置信度细胞为 Unknown")

        # 相似簇合并
        if merge_similar_threshold > 0:
            import numpy as np
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if score_cols:
                cluster_profiles = {}
                for cluster in adata.obs[leiden_key].cat.categories:
                    mask = adata.obs[leiden_key] == cluster
                    cluster_profiles[cluster] = adata.obs.loc[mask, score_cols].mean().values
                clusters = list(cluster_profiles.keys())
                for i in range(len(clusters)):
                    for j in range(i + 1, len(clusters)):
                        corr = np.corrcoef(cluster_profiles[clusters[i]], cluster_profiles[clusters[j]])[0, 1]
                        if corr >= merge_similar_threshold:
                            new_ct = adata.obs.loc[adata.obs[leiden_key] == clusters[i], 'celltype'].mode()
                            new_ct = new_ct.iloc[0] if not new_ct.empty else 'Unknown'
                            adata.obs.loc[adata.obs[leiden_key] == clusters[j], 'celltype'] = new_ct
                            self.progress(-1, f"合并簇 {clusters[j]} → {clusters[i]} (r={corr:.2f})")

        self.progress(70, "Generating dotplot and UMAP...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        if marker_coverage:
            coverage_types = list(marker_coverage)
            matched = [marker_coverage[ct]['matched'] for ct in coverage_types]
            totals = [marker_coverage[ct]['total'] for ct in coverage_types]
            fig_coverage = go.Figure()
            fig_coverage.add_trace(go.Bar(name='可用 marker', x=coverage_types, y=matched, marker_color='#00897b'))
            fig_coverage.add_trace(go.Bar(name='marker 总数', x=coverage_types, y=totals, marker_color='#b0bec5'))
            fig_coverage.update_layout(title='Marker 覆盖度（当前数据）', barmode='group',
                                       xaxis=dict(tickangle=35), yaxis_title='基因数', plot_bgcolor='white',
                                       width=max(750, 80 * len(coverage_types)), height=450)
            result_files.append(self.save_plotly_json(
                fig_coverage, plots_dir, 'annotation_marker_coverage.json', 'bar', 'Marker 覆盖度'
            ))

        if self.params.get('show_celltype_composition', True) and 'celltype' in adata.obs.columns:
            ct_counts_plot = adata.obs['celltype'].astype(str).value_counts()
            ct_pct = ct_counts_plot / max(int(ct_counts_plot.sum()), 1) * 100
            fig_comp = go.Figure()
            fig_comp.add_trace(go.Bar(
                x=ct_counts_plot.index.tolist(),
                y=ct_counts_plot.values.astype(int).tolist(),
                marker_color='#00897b',
                customdata=ct_pct.round(2).values,
                text=[f'{p:.1f}%' for p in ct_pct.values],
                textposition='outside',
                hovertemplate='Cell type: %{x}<br>Cells: %{y}<br>Percent: %{customdata:.2f}%<extra></extra>',
            ))
            fig_comp.update_layout(
                title='Cell Type Composition',
                xaxis_title='Cell type',
                yaxis_title='Cell count',
                plot_bgcolor='white',
                width=max(700, 85 * max(1, len(ct_counts_plot))),
                height=460,
                xaxis=dict(tickangle=35),
            )
            result_files.append(self.save_plotly_json(
                fig_comp, plots_dir, 'annotation_celltype_composition.json',
                'bar', '细胞类型组成'
            ))

        if self.params.get('show_marker_score_heatmap', True):
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if score_cols and leiden_key in adata.obs.columns:
                cluster_labels = adata.obs[leiden_key].astype(str)
                cluster_order = sorted(cluster_labels.unique(), key=lambda x: (len(x), x))
                mean_scores = []
                for cluster in cluster_order:
                    mask = cluster_labels == cluster
                    mean_scores.append(adata.obs.loc[mask, score_cols].mean().values)
                score_matrix = np.asarray(mean_scores, dtype=float).T
                row_mean = score_matrix.mean(axis=1, keepdims=True)
                row_std = score_matrix.std(axis=1, keepdims=True) + 1e-10
                z = np.clip((score_matrix - row_mean) / row_std, -3, 3)
                fig_score_heat = go.Figure(data=go.Heatmap(
                    z=z,
                    x=cluster_order,
                    y=[c.replace('score_', '') for c in score_cols],
                    colorscale='RdBu',
                    zmid=0,
                    colorbar=dict(title='Row z-score'),
                    hovertemplate='Cluster: %{x}<br>Marker set: %{y}<br>z-score: %{z:.2f}<extra></extra>',
                ))
                fig_score_heat.update_layout(
                    title=f'Marker Score Heatmap by {leiden_key}',
                    xaxis_title='Cluster',
                    yaxis_title='Marker set',
                    plot_bgcolor='white',
                    width=max(700, 55 * max(1, len(cluster_order))),
                    height=max(450, 24 * max(1, len(score_cols))),
                )
                result_files.append(self.save_plotly_json(
                    fig_score_heat, plots_dir, 'annotation_marker_score_heatmap.json',
                    'heatmap', 'Marker Score Heatmap'
                ))

        # Dotplot for marker validation
        dotplot_genes = []
        for genes in markers.values():
            dotplot_genes.extend([g for g in genes[:2] if g in adata.var_names])
        dotplot_genes = list(dict.fromkeys(dotplot_genes))[:25]

        if dotplot_genes and 'celltype' in adata.obs.columns:
            try:
                sc.tl.dendrogram(adata, groupby='celltype')
                fig_dotplot = sc.pl.dotplot(adata, var_names=dotplot_genes, groupby='celltype', return_fig=True)
                result_files.extend(self.save_matplotlib_figure(
                    fig_dotplot, plots_dir, 'annotation_celltype_dotplot.png', 'dotplot',
                    'Cell Type Marker Dotplot'
                ))
                import matplotlib.pyplot as plt
                plt.close('all')
            except Exception as e:
                self.progress(-1, f"Annotation dotplot generation failed: {e}")

        # Marker expression box/violin plot for validating annotations
        if self.params.get('show_marker_expression_violin', True) and dotplot_genes and 'celltype' in adata.obs.columns:
            try:
                marker_genes = dotplot_genes[:8]
                idx = np.arange(adata.n_obs)
                if adata.n_obs > 5000:
                    rng = np.random.default_rng(0)
                    idx = np.sort(rng.choice(adata.n_obs, 5000, replace=False))
                expr = adata[idx, marker_genes].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray()
                expr = np.asarray(expr)
                celltypes = adata.obs['celltype'].astype(str).iloc[idx].values
                fig_marker = go.Figure()
                for gi, gene in enumerate(marker_genes):
                    fig_marker.add_trace(go.Box(
                        x=celltypes,
                        y=expr[:, gi],
                        name=gene,
                        boxpoints=False,
                    ))
                fig_marker.update_layout(
                    title='Marker Expression by Cell Type',
                    xaxis_title='Cell type',
                    yaxis_title='Expression',
                    boxmode='group',
                    plot_bgcolor='white',
                    width=max(850, 90 * max(1, adata.obs['celltype'].nunique())),
                    height=480,
                    xaxis=dict(tickangle=35),
                )
                result_files.append(self.save_plotly_json(
                    fig_marker, plots_dir, 'annotation_marker_expression_box.json',
                    'boxplot', 'Marker 表达验证图'
                ))
            except Exception as e:
                self.progress(-1, f"Marker expression plot generation failed: {e}")

        fig_json = json.dumps(umap_scatter(adata, 'celltype', title='UMAP by Cell Type'))
        fpath = os.path.join(plots_dir, 'annotation_umap_celltype.json')
        with open(fpath, 'w') as f: f.write(fig_json)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': 'UMAP by Cell Type'})
        try:
            fig_static = self.build_publication_umap(
                adata, 'celltype', title='UMAP by Cell Type'
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_static, plots_dir, 'annotation_umap_celltype.png', 'umap',
                'UMAP by Cell Type'
            ))
            import matplotlib.pyplot as plt
            plt.close(fig_static)
        except Exception as exc:
            self.progress(-1, f'静态注释 UMAP 导出失败（不影响交互图）：{exc}')

        # Annotation score margin / confidence UMAP
        if self.params.get('show_annotation_score_umap', True) and 'X_umap' in adata.obsm:
            score_col = None
            score_label = None
            if 'annotation_score_margin' in adata.obs.columns:
                score_col = 'annotation_score_margin'
                score_label = 'Annotation score margin'
            elif 'annotation_confidence' in adata.obs.columns:
                score_col = 'annotation_confidence'
                score_label = 'Annotation confidence'
            if score_col:
                coords = adata.obsm['X_umap'][:, :2]
                vals = adata.obs[score_col].astype(float).values
                fig_score = go.Figure()
                fig_score.add_trace(go.Scattergl(
                    x=coords[:, 0],
                    y=coords[:, 1],
                    mode='markers',
                    marker=dict(
                        size=4,
                        color=vals,
                        colorscale='Viridis',
                        opacity=0.75,
                        colorbar=dict(title=score_label),
                    ),
                    text=adata.obs_names.tolist(),
                    hovertemplate='%{text}<br>' + score_label + ': %{marker.color:.3f}<extra></extra>',
                ))
                fig_score.update_layout(
                    title=score_label + ' on UMAP',
                    xaxis_title='UMAP-1',
                    yaxis_title='UMAP-2',
                    plot_bgcolor='white',
                    width=700,
                    height=520,
                )
                result_files.append(self.save_plotly_json(
                    fig_score, plots_dir, 'annotation_score_umap.json',
                    'umap', score_label + ' UMAP'
                ))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'annotation')

        ct_counts = adata.obs['celltype'].value_counts().to_dict()
        self.progress(100, "Done")
        result = {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_celltypes': adata.obs['celltype'].nunique(),
                'celltype_counts': {str(k): int(v) for k, v in ct_counts.items()},
                'cluster_column': leiden_key,
                'method_used': requested_method,
                'marker_set': marker_set_name if method == 'auto_marker' or multi_evidence else None,
                'marker_coverage': marker_coverage,
            }
        }
        if 'annotation_confidence' in adata.obs.columns:
            confidence_value = round(float(adata.obs['annotation_confidence'].mean()), 3)
            result['summary']['mean_confidence'] = confidence_value
        if 'annotation_score_margin' in adata.obs.columns:
            result['summary']['mean_score_margin'] = round(float(adata.obs['annotation_score_margin'].mean()), 3)
        return result
