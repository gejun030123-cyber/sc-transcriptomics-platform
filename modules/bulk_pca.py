import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


def _resolve_pca_color_groups(adata, requested_color_by='', auto_group_mapping=None):
    """Resolve PCA colouring without silently collapsing distinct sample groups."""
    from modules.io_utils import infer_sample_group_candidates

    requested = str(requested_color_by or '').strip()
    sample_names = [str(name) for name in adata.obs.index]
    warning = ''

    if requested and requested not in {'_auto_group_', '_auto_group', 'auto_group'}:
        if requested in adata.obs.columns:
            values = adata.obs[requested].astype(str).tolist()
            if len(set(values)) >= 2:
                return values, {
                    'requested': requested, 'used': requested, 'source': 'obs',
                    'group_counts': {str(k): int(v) for k, v in pd.Series(values).value_counts().items()},
                    'warning': '',
                }
            warning = f"obs 列 '{requested}' 只有一个分组，已尝试从样本名恢复分组。"
        else:
            warning = f"未找到 obs 列 '{requested}'，已尝试从样本名恢复分组。"

    mapping = auto_group_mapping or {}
    if isinstance(mapping, str):
        try:
            mapping = json.loads(mapping)
        except (TypeError, ValueError):
            mapping = {}
    if isinstance(mapping, dict) and mapping:
        mapped_values = [str(mapping.get(name, '')) for name in sample_names]
        if all(mapped_values) and len(set(mapped_values)) >= 2:
            adata.obs['_auto_group'] = mapped_values
            return mapped_values, {
                'requested': requested or '_auto_group_', 'used': '_auto_group',
                'source': 'sample_name_mapping',
                'group_counts': {str(k): int(v) for k, v in pd.Series(mapped_values).value_counts().items()},
                'warning': warning,
            }

    candidates = infer_sample_group_candidates(sample_names)
    if candidates:
        mapped_values = [str(candidates[0]['mapping'][name]) for name in sample_names]
        adata.obs['_auto_group'] = mapped_values
        return mapped_values, {
            'requested': requested or '_auto_group_', 'used': '_auto_group',
            'source': 'sample_name_inference',
            'group_counts': {str(k): int(v) for k, v in pd.Series(mapped_values).value_counts().items()},
            'warning': warning,
        }

    if not warning and requested in {'_auto_group_', '_auto_group', 'auto_group'}:
        warning = '样本名中未识别到可重复的分组模式，PCA 将以统一颜色显示。'
    return None, {
        'requested': requested or '_auto_group_', 'used': '', 'source': 'none',
        'group_counts': {}, 'warning': warning,
    }


def _pca_label_positions(coords):
    """Assign deterministic, separated text positions for a small PCA panel."""
    coords = np.asarray(coords, dtype=float)
    x_span = max(float(np.ptp(coords[:, 0])), 1.0)
    y_span = max(float(np.ptp(coords[:, 1])), 1.0)
    directions = [(1, 1), (1, -1), (-1, 1), (-1, -1), (0, 1), (0, -1)]
    positions = []
    for index, (x_value, y_value) in enumerate(coords):
        direction = directions[index % len(directions)]
        step = 0
        while True:
            scale = 0.025 + 0.018 * step
            candidate = (
                x_value + direction[0] * x_span * scale,
                y_value + direction[1] * y_span * scale,
            )
            if not positions or all(
                abs(candidate[0] - old_x) / x_span + abs(candidate[1] - old_y) / y_span >= 0.10
                for old_x, old_y in positions
            ):
                positions.append(candidate)
                break
            step += 1
    return positions


def _pca_embedding_figure(coords, group_values, sample_labels, title, x_label, y_label,
                          group_label='Group', show_labels=True, subtitle='',
                          max_sample_labels=8):
    """Draw a legible PCA/UMAP embedding with categorical sample identities."""
    import matplotlib.pyplot as plt
    from modules.figure_style import NATURE_AXIS, NATURE_GRID, NATURE_PALETTE, NATURE_TEXT

    coords = np.asarray(coords, dtype=float)
    labels = [str(label) for label in sample_labels]
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError('Embedding coordinates must contain at least two dimensions.')
    if len(labels) != coords.shape[0]:
        raise ValueError('Sample labels must match embedding coordinates.')

    grouped = group_values is not None and len(set(map(str, group_values))) >= 2
    fig, ax = plt.subplots(figsize=(7.8, 5.5), dpi=150)
    if grouped:
        groups = [str(value) for value in group_values]
        unique_groups = list(dict.fromkeys(groups))
        # The shared Nature palette has ten colours.  Add restrained, distinct
        # accents before cycling so 11–13 experimental groups never silently
        # reuse a legend colour.
        categorical_palette = list(NATURE_PALETTE) + ['#2E7D6E', '#A5663F', '#7561A8']
        color_map = {
            group: categorical_palette[index % len(categorical_palette)]
            for index, group in enumerate(unique_groups)
        }
        for group in unique_groups:
            mask = np.asarray([value == group for value in groups])
            ax.scatter(
                coords[mask, 0], coords[mask, 1], s=52,
                color=color_map[group], edgecolor='white', linewidth=0.75,
                alpha=0.92, label=f'{group} (n={int(mask.sum())})', zorder=3,
            )
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=50, color=NATURE_PALETTE[0],
                   edgecolor='white', linewidth=0.75, alpha=0.90, zorder=3)

    # Static PCA plots become illegible quickly when every replicate is named.
    # The legend carries the experimental grouping; reserve point labels for
    # genuinely small panels where they add information rather than collisions.
    if show_labels and len(labels) <= max_sample_labels:
        for label, (x_value, y_value), (text_x, text_y) in zip(
                labels, coords[:, :2], _pca_label_positions(coords[:, :2])):
            ax.annotate(
                label, xy=(x_value, y_value), xytext=(text_x, text_y),
                textcoords='data', ha='left' if text_x >= x_value else 'right',
                va='center', fontsize=7, color=NATURE_TEXT, zorder=5,
                bbox={'boxstyle': 'round,pad=0.18', 'facecolor': 'white',
                      'edgecolor': '#D0D5DD', 'linewidth': 0.45, 'alpha': 0.92},
                arrowprops={'arrowstyle': '-', 'color': '#98A2B3', 'linewidth': 0.45,
                            'shrinkA': 2, 'shrinkB': 2},
            )

    if np.nanmin(coords[:, 0]) < 0 < np.nanmax(coords[:, 0]):
        ax.axvline(0, color='#D0D5DD', linewidth=0.7, zorder=1)
    if np.nanmin(coords[:, 1]) < 0 < np.nanmax(coords[:, 1]):
        ax.axhline(0, color='#D0D5DD', linewidth=0.7, zorder=1)
    ax.set_title(title, loc='left', pad=20, fontsize=10,
                 fontweight='semibold', color=NATURE_TEXT)
    if subtitle:
        ax.text(0, 1.015, subtitle, transform=ax.transAxes, ha='left', va='bottom',
                fontsize=7.5, color='#667085')
    ax.set_xlabel(x_label, fontsize=9, color=NATURE_TEXT)
    ax.set_ylabel(y_label, fontsize=9, color=NATURE_TEXT)
    ax.grid(color=NATURE_GRID, linewidth=0.5, alpha=0.72)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, colors=NATURE_AXIS, width=0.7, length=3)
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name in ('left', 'bottom'))
        spine.set_color(NATURE_AXIS)
        spine.set_linewidth(0.7)
    if grouped:
        legend = ax.legend(title=group_label, loc='center left', bbox_to_anchor=(1.01, 0.5),
                           frameon=False, fontsize=7.5, title_fontsize=8,
                           handletextpad=0.5, borderaxespad=0)
        for handle in legend.legend_handles:
            handle.set_sizes([34])
        # Reserve the header band for the title and PC1+PC2 variance subtitle;
        # BaseAnalysis reuses this rect for final PNG/SVG export.
        fig._native_layout_rect = (0.0, 0.0, 0.77, 0.86)
        fig.tight_layout(rect=fig._native_layout_rect, pad=1.1)
    else:
        fig.tight_layout(pad=1.1)
    return fig


class BulkPCAAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_pca"
    DISPLAY_NAME = "Bulk PCA / UMAP 降维"
    DESCRIPTION = "对 Bulk RNA-seq 数据进行 PCA 和 UMAP 降维可视化"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        from modules.native_figures import bar_figure, line_figure
        from modules.figure_style import NATURE_PALETTE

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)
        from modules.io_utils import infer_expression_measurement

        # 清理 inf/NaN
        import numpy as _np
        adata.X = _np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

        n_comps = int(self.params.get('n_comps', 10))
        color_by = self.params.get('color_by', '')
        dimred_method = self.params.get('dimred_method', 'pca')

        self.progress(20, "标准化数据...")
        if 'normalization' not in adata.uns:
            input_measurement = infer_expression_measurement(adata, input_path)
            if input_measurement == 'raw_counts':
                sc.pp.normalize_total(adata, target_sum=1e6)
                sc.pp.log1p(adata)
                adata.uns['normalization'] = {'method': 'pca_auto_cpm_log1p', 'is_log_transformed': True}
            else:
                adata.X = np.log2(np.maximum(adata.X, 0) + 1)
                adata.uns['normalization'] = {'method': 'pca_auto_log2', 'is_log_transformed': True}
        sc.pp.scale(adata, max_value=10)

        self.progress(40, "运行 PCA...")
        actual_comps = min(n_comps, adata.n_obs - 1, adata.n_vars - 1)
        if actual_comps < 2:
            raise ValueError(f"样本数不足 ({adata.n_obs})，至少需要 3 个样本才能进行 PCA 分析。")
        sc.pp.pca(adata, n_comps=actual_comps)
        pca_variance = adata.uns['pca']['variance_ratio']

        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        color_values, color_info = _resolve_pca_color_groups(
            adata, color_by, self.params.get('_auto_group_mapping', {}))
        if color_info['warning']:
            self.progress(-1, color_info['warning'])

        # Sample names are useful only in small static panels.  Larger PCA
        # panels use the external legend so labels cannot cover nearby points.
        show_text = adata.n_obs <= 8

        self.progress(55, "生成 PCA 图...")
        pc = adata.obsm['X_pca']
        hover = adata.obs.index.tolist()
        color_descriptor = (
            'sample-name group' if color_info['source'].startswith('sample_name')
            else (color_info['used'] or 'single colour')
        )

        def _embedding_figure(coords, title, x_label, y_label):
            return _pca_embedding_figure(
                coords, color_values, hover, title, x_label, y_label,
                group_label='Group', show_labels=show_text,
                subtitle=(
                    f"Color: {color_descriptor} · "
                    f"PC1 + PC2: {(pca_variance[0] + pca_variance[1]) * 100:.1f}%"
                    if coords is pc else ''
                ),
            )

        fig_pca = _embedding_figure(
            pc, f'PCA 分析 (n={adata.n_obs})',
            f'PC1 ({pca_variance[0]*100:.1f}% variance)',
            f'PC2 ({pca_variance[1]*100:.1f}% variance)',
        )
        result_files.extend(self.save_matplotlib_figure(
            fig_pca, plots_dir, 'bulk_pca.png', 'pca', 'PCA 分析',
            formats=('png', 'svg'), dpi=300,
        ))

        # PC3 often carries a biologically meaningful secondary separation when
        # PC1+PC2 explain a moderate fraction of the total transcriptome variance.
        if actual_comps >= 3:
            for first, second, stem in [(0, 2, 'bulk_pca_pc1_pc3'), (1, 2, 'bulk_pca_pc2_pc3')]:
                fig_alt = _embedding_figure(
                    pc[:, [first, second]], f'PCA 分析：PC{first + 1} vs PC{second + 1}',
                    f'PC{first + 1} ({pca_variance[first] * 100:.1f}% variance)',
                    f'PC{second + 1} ({pca_variance[second] * 100:.1f}% variance)',
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_alt, plots_dir, stem, 'pca',
                    f'PCA：PC{first + 1} vs PC{second + 1}', formats=('png', 'svg'), dpi=300,
                ))

        self.progress(65, "生成方差解释图...")
        fig_var = bar_figure(
            [f'PC{i+1}' for i in range(actual_comps)], pca_variance * 100,
            title='PCA 方差解释比例', x_label='主成分',
            y_label='方差解释比例 (%)', rotation=45,
        )
        result_files.extend(self.save_matplotlib_figure(
            fig_var, plots_dir, 'bulk_pca_variance.png', 'pca',
            '方差解释比例', formats=('png', 'svg'), dpi=300,
        ))

        # PCA 载荷图
        if 'PCs' in adata.varm:
            loadings = adata.varm['PCs'][:, :2]
            for pc_idx, pc_name in enumerate(['PC1', 'PC2']):
                top_idx = np.argsort(np.abs(loadings[:, pc_idx]))[::-1][:10]
                values = loadings[top_idx, pc_idx]
                fig_load = bar_figure(
                    [adata.var_names[i] for i in top_idx], values,
                    title=f'{pc_name} Top 10 载荷基因', x_label='Gene',
                    y_label='Loading',
                    colors=[NATURE_PALETTE[3] if v > 0 else NATURE_PALETTE[0] for v in values],
                    rotation=45,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_load, plots_dir, f'bulk_pca_loadings_{pc_name.lower()}.png',
                    'pca', f'{pc_name} 载荷图', formats=('png', 'svg'), dpi=300,
                ))

        # 肘部图（方差累积曲线）
        cumvar = np.cumsum(pca_variance) * 100
        fig_elbow = line_figure(
            [f'PC{i+1}' for i in range(actual_comps)], [cumvar],
            title='PCA 方差累积曲线（肘部图）', x_label='主成分',
            y_label='累积方差 (%)', labels=['Cumulative'],
        )
        ax_elbow = fig_elbow.axes[0]
        ax_elbow.axhline(80, color='#98A2B3', linestyle='--', linewidth=0.8)
        ax_elbow.text(0.98, 80, '80%', transform=ax_elbow.get_yaxis_transform(),
                      ha='right', va='bottom', fontsize=7, color='#667085')
        result_files.extend(self.save_matplotlib_figure(
            fig_elbow, plots_dir, 'bulk_pca_elbow.png', 'pca', '肘部图',
            formats=('png', 'svg'), dpi=300,
        ))

        # 降维方法选择：t-SNE / UMAP / PCA-only
        if dimred_method == 'tsne' and adata.n_obs >= 3:
            from sklearn.manifold import TSNE
            self.progress(75, "运行 t-SNE...")
            perplexity = min(30, adata.n_obs - 1)
            tsne = TSNE(n_components=2, random_state=42, perplexity=max(2, perplexity))
            tsne_coords = tsne.fit_transform(adata.obsm['X_pca'])

            fig_tsne = _embedding_figure(
                tsne_coords, f't-SNE 分析 (n={adata.n_obs})', 't-SNE1', 't-SNE2',
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_tsne, plots_dir, 'bulk_tsne.png', 'tsne', 't-SNE 分析',
                formats=('png', 'svg'), dpi=300,
            ))

        elif dimred_method == 'umap' and adata.n_obs >= 10:
            self.progress(75, "运行 UMAP...")
            sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1))
            sc.tl.umap(adata)
            umap_coords = adata.obsm['X_umap']

            fig_umap = _embedding_figure(umap_coords, 'UMAP 分析', 'UMAP1', 'UMAP2')
            result_files.extend(self.save_matplotlib_figure(
                fig_umap, plots_dir, 'bulk_umap.png', 'umap', 'UMAP 分析',
                formats=('png', 'svg'), dpi=300,
            ))

        self.progress(90, "保存结果...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_pca_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_samples': adata.n_obs,
                'n_genes': adata.n_vars,
                'n_components': actual_comps,
                'pc1_variance_pct': round(float(pca_variance[0] * 100), 2),
                'pc2_variance_pct': round(float(pca_variance[1] * 100), 2),
                'pc1_pc2_variance_pct': round(float((pca_variance[0] + pca_variance[1]) * 100), 2),
                'color_by_requested': color_info['requested'],
                'color_by_used': color_info['used'],
                'grouping_source': color_info['source'],
                'group_counts': color_info['group_counts'],
                'grouping_warning': color_info['warning'],
                'dimred_method': dimred_method,
            }
        }
