import os
from modules.base import BaseAnalysis

class CellCommunicationAnalysis(BaseAnalysis):
    MODULE_NAME = "cell_communication"
    DISPLAY_NAME = "细胞通讯"
    DESCRIPTION = "基于 LIANA 的细胞间通讯分析"
    INPUT_REQUIRES = ['celltype']

    def run(self, input_path):
        import scanpy as sc
        import json
        import pandas as pd

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)

        cluster_key = self.params.get('cluster_key', 'celltype')
        resource = self.params.get('resource', 'consensus')
        organism = self.params.get('organism', 'human')
        min_prop = float(self.params.get('min_prop', 0.1))
        top_n = int(self.params.get('top_n_interactions', 20))
        show_heatmap = self.params.get('show_heatmap', True)

        self.progress(20, "Running LIANA cell communication analysis...")
        try:
            import liana as li
        except ImportError:
            self.progress(100, "LIANA not installed")
            return {
                'output_adata': input_path,
                'result_files': [],
                'summary': {'error': 'LIANA is not installed. Install with: pip install liana'}
            }

        # Run LIANA rank aggregate
        self.progress(30, f"Computing interactions (resource={resource})...")
        li.mt.rank_aggregate(
            adata,
            groupby=cluster_key,
            resource_name=resource,
            organism=organism,
            min_prop=min_prop,
            verbose=False,
        )

        result_files = []
        plots_dir = self.ensure_plots_dir()
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)

        # Extract LIANA results
        liana_results = adata.uns.get('liana_res', pd.DataFrame())
        if liana_results.empty:
            self.progress(100, "No interactions found")
            return {
                'output_adata': input_path,
                'result_files': [],
                'summary': {'error': 'No significant interactions found'}
            }

        # Save full results CSV
        self.progress(60, "Saving interaction results...")
        csv_path = os.path.join(results_dir, 'cell_communication_results.csv')
        liana_results.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': 'Cell Communication Results'})

        # Top interactions
        top_interactions = liana_results.head(top_n)

        # Generate bubble plot for top interactions
        self.progress(70, "Generating interaction plots...")
        try:
            import plotly.graph_objects as go
            fig = go.Figure()
            for _, row in top_interactions.iterrows():
                source = row.get('source', '')
                target = row.get('target', '')
                ligand = row.get('ligand_complex', row.get('ligand', ''))
                receptor = row.get('receptor_complex', row.get('receptor', ''))
                lr_score = row.get('lr_means', row.get('aggregate_rank', 0))
                fig.add_trace(go.Scattergl(
                    x=[f"{ligand}-{receptor}"],
                    y=[f"{source}→{target}"],
                    mode='markers',
                    marker=dict(size=max(3, float(lr_score) * 20) if lr_score else 5,
                               color=float(lr_score) if lr_score else 0,
                               colorscale='Reds', showscale=True, colorbar=dict(title='Score')),
                    name=f"{source}→{target}",
                    showlegend=False,
                ))
            fig.update_layout(title=f'Top {top_n} Cell Communication Interactions',
                             xaxis_title='Ligand-Receptor', yaxis_title='Source→Target',
                             plot_bgcolor='white', width=800, height=500,
                             xaxis=dict(tickangle=45))
            result_files.append(self.save_plotly_json(fig, plots_dir, 'cell_communication_bubble.json', 'bubble', 'Communication Bubble Plot'))
        except Exception:
            pass

        # Heatmap of interaction counts between cell types
        if show_heatmap:
            try:
                import plotly.graph_objects as go
                if 'source' in liana_results.columns and 'target' in liana_results.columns:
                    counts = liana_results.groupby(['source', 'target']).size().reset_index(name='count')
                    pivot = counts.pivot(index='source', columns='target', values='count').fillna(0)
                    fig_hm = go.Figure(data=go.Heatmap(
                        z=pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
                        colorscale='Reds', text=pivot.values.astype(int), texttemplate='%{text}',
                    ))
                    fig_hm.update_layout(title='Communication Counts Between Cell Types',
                                        xaxis_title='Target', yaxis_title='Source',
                                        width=600, height=500)
                    result_files.append(self.save_plotly_json(fig_hm, plots_dir, 'cell_communication_heatmap.json', 'heatmap', 'Communication Heatmap'))
            except Exception:
                pass

        # Save output
        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'cell_communication')

        n_interactions = len(liana_results)
        n_cell_types = liana_results['source'].nunique() if 'source' in liana_results.columns else 0

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_interactions': n_interactions,
                'n_cell_types': n_cell_types,
                'resource': resource,
                'organism': organism,
            }
        }
