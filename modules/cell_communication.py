import os
import logging
from modules.base import BaseAnalysis
from modules.io_utils import obs_grouping_info

logger = logging.getLogger(__name__)

class CellCommunicationAnalysis(BaseAnalysis):
    MODULE_NAME = "cell_communication"
    DISPLAY_NAME = "细胞通讯"
    DESCRIPTION = "基于 LIANA 的细胞间通讯分析"
    INPUT_REQUIRES = ['celltype']

    def validate_input(self, adata):
        cluster_key = str(self.params.get('cluster_key', 'celltype') or '').strip()
        info = obs_grouping_info(
            adata, cluster_key, max_categories=50,
            max_numeric_categories=20, require_multiple=True,
        )
        if not info['valid']:
            return f"cluster_key '{cluster_key}' 不是有效的分类细胞类型列：{info['reason']}"
        return None

    def run(self, input_path):
        import scanpy as sc
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        from modules.native_figures import heatmap_figure
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT, NATURE_GRID

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)

        cluster_key = str(self.params.get('cluster_key', 'celltype') or '').strip()
        cluster_info = obs_grouping_info(
            adata, cluster_key, max_categories=50,
            max_numeric_categories=20, require_multiple=True,
        )
        if not cluster_info['valid']:
            raise ValueError(
                f"cluster_key '{cluster_key}' 不是有效的分类细胞类型列："
                f"{cluster_info['reason']}"
            )
        if not hasattr(adata.obs[cluster_key].dtype, 'categories'):
            adata.obs[cluster_key] = adata.obs[cluster_key].astype(str).astype('category')
        resource = self.params.get('resource', 'consensus')
        organism = self.params.get('organism', 'human')
        min_prop = float(self.params.get('min_prop', 0.1))
        top_n = int(self.params.get('top_n_interactions', 20))
        show_heatmap = self.params.get('show_heatmap', True)

        self.progress(20, "Running LIANA cell communication analysis...")
        try:
            import liana as li
        except ImportError:
            self.progress(-1, "LIANA 未安装，跳过细胞通讯分析")
            return {
                'output_adata': input_path,
                'result_files': [],
                'summary': {
                    'status': 'unavailable',
                    'error': 'LIANA is not installed. Install with: pip install liana',
                    'missing_dependencies': ['liana'],
                    'next_step': '安装 liana 后重新运行；当前输入数据未被覆盖。',
                }
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
            self.progress(-1, "未检测到细胞间相互作用")
            return {
                'output_adata': input_path,
                'result_files': [],
                'summary': {
                    'status': 'completed_with_warning',
                    'warning': 'No significant interactions found',
                }
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
            bubble = top_interactions.copy()
            bubble['ligand_receptor'] = bubble.apply(
                lambda row: f"{row.get('ligand_complex', row.get('ligand', ''))}-"
                f"{row.get('receptor_complex', row.get('receptor', ''))}", axis=1)
            bubble['source_target'] = bubble.apply(
                lambda row: f"{row.get('source', '')}→{row.get('target', '')}", axis=1)
            x_labels = list(dict.fromkeys(bubble['ligand_receptor'].astype(str)))
            y_labels = list(dict.fromkeys(bubble['source_target'].astype(str)))
            x_pos = {value: index for index, value in enumerate(x_labels)}
            y_pos = {value: index for index, value in enumerate(y_labels)}
            scores = []
            for _, row in bubble.iterrows():
                lr_score = row.get('lr_means', row.get('aggregate_rank', 0))
                try:
                    scores.append(float(lr_score))
                except (TypeError, ValueError):
                    scores.append(0.0)
            scores = np.nan_to_num(np.asarray(scores, dtype=float), nan=0.0)
            fig, ax = plt.subplots(figsize=(9.0, 5.8), dpi=150)
            sizes = 35 + 260 * (scores - scores.min()) / (np.ptp(scores) + 1e-9)
            artist = ax.scatter(
                [x_pos[value] for value in bubble['ligand_receptor'].astype(str)],
                [y_pos[value] for value in bubble['source_target'].astype(str)],
                s=sizes, c=scores, cmap='Reds', alpha=0.8,
                edgecolors='white', linewidths=0.45,
            )
            for _, row in bubble.iterrows():
                source = row.get('source', '')
                target = row.get('target', '')
                ligand = row.get('ligand_complex', row.get('ligand', ''))
                receptor = row.get('receptor_complex', row.get('receptor', ''))
            ax.set_xticks(np.arange(len(x_labels)), x_labels, rotation=45, ha='right', fontsize=8)
            ax.set_yticks(np.arange(len(y_labels)), y_labels, fontsize=8)
            ax.set_xlabel('Ligand-Receptor', fontsize=9)
            ax.set_ylabel('Source→Target', fontsize=9)
            ax.set_title(f'Top {top_n} Cell Communication Interactions', loc='left',
                         fontsize=10, fontweight='semibold', color=NATURE_TEXT)
            colorbar = fig.colorbar(artist, ax=ax, fraction=0.035, pad=0.025)
            colorbar.set_label('Score', fontsize=8)
            ax.grid(False)
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'cell_communication_bubble.png', 'bubble',
                'Communication Bubble Plot', formats=('png', 'svg'), dpi=300,
            ))
        except Exception as e:
            logger.warning("生成气泡图失败: %s", e)

        # Heatmap of interaction counts between cell types
        if show_heatmap:
            try:
                if 'source' in liana_results.columns and 'target' in liana_results.columns:
                    counts = liana_results.groupby(['source', 'target']).size().reset_index(name='count')
                    pivot = counts.pivot(index='source', columns='target', values='count').fillna(0)
                    fig_hm = heatmap_figure(
                        pivot.values, x_labels=pivot.columns.tolist(), y_labels=pivot.index.tolist(),
                        title='Communication Counts Between Cell Types', x_label='Target',
                        y_label='Source', colorbar_label='Interactions', vmin=0,
                        vmax=max(1.0, float(pivot.values.max())),
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_hm, plots_dir, 'cell_communication_heatmap.png',
                        'heatmap', 'Communication Heatmap', formats=('png', 'svg'), dpi=300,
                    ))
            except Exception as e:
                logger.warning("生成热图失败: %s", e)

        if self.params.get('show_network', True):
            try:
                if 'source' in liana_results.columns and 'target' in liana_results.columns:
                    edge_df = (
                        liana_results.groupby(['source', 'target'])
                        .size()
                        .reset_index(name='count')
                        .sort_values('count', ascending=False)
                    )
                    edge_df = edge_df.iloc[:max(1, top_n)].copy()
                    nodes = sorted(set(edge_df['source'].astype(str)) | set(edge_df['target'].astype(str)))
                    if nodes:
                        angles = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)
                        positions = {node: (float(np.cos(a)), float(np.sin(a))) for node, a in zip(nodes, angles)}
                        max_count = max(float(edge_df['count'].max()), 1.0)
                        fig_net, ax_net = plt.subplots(figsize=(7.5, 6.5), dpi=150)
                        for _, row in edge_df.iterrows():
                            source = str(row['source'])
                            target = str(row['target'])
                            count = float(row['count'])
                            x0, y0 = positions[source]
                            x1, y1 = positions[target]
                            ax_net.plot([x0, x1], [y0, y1], color=NATURE_PALETTE[3],
                                        alpha=0.32, linewidth=1 + 7 * count / max_count,
                                        zorder=1)
                        degrees = {node: 0 for node in nodes}
                        for _, row in edge_df.iterrows():
                            degrees[str(row['source'])] += int(row['count'])
                            degrees[str(row['target'])] += int(row['count'])
                        node_artist = ax_net.scatter(
                            [positions[n][0] for n in nodes], [positions[n][1] for n in nodes],
                            s=[max(140, min(900, 100 + 24 * degrees[n])) for n in nodes],
                            c=[degrees[n] for n in nodes], cmap='Blues',
                            edgecolors='white', linewidths=0.9, zorder=2,
                        )
                        for node in nodes:
                            x_node, y_node = positions[node]
                            ax_net.text(x_node, y_node, node, ha='center', va='center',
                                        fontsize=8, color='white', zorder=3)
                        ax_net.set_title(f'Top {len(edge_df)} Cell Communication Network',
                                         loc='left', fontsize=10, fontweight='semibold',
                                         color=NATURE_TEXT)
                        ax_net.axis('off')
                        colorbar = fig_net.colorbar(node_artist, ax=ax_net, fraction=0.035, pad=0.02)
                        colorbar.set_label('Degree', fontsize=8)
                        result_files.extend(self.save_matplotlib_figure(
                            fig_net, plots_dir, 'cell_communication_network.png',
                            'network', 'Communication Network', formats=('png', 'svg'), dpi=300,
                        ))
            except Exception as e:
                logger.warning("生成通讯网络图失败: %s", e)

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
                'cluster_key': cluster_key,
            }
        }
