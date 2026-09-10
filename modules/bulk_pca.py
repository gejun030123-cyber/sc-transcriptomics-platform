import os
import re
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


_AUTO_FACTOR_LABEL = re.compile(r'^第\s*(\d+)\s*因素(?:\s*[:：].*)?$', re.IGNORECASE)


def _auto_column_from_requested(requested):
    """Translate supported UI aliases to real, materialized ``obs`` columns."""
    normalized = str(requested or '').strip()
    if normalized in {'_auto_group_', '_auto_group', 'auto_group'}:
        return '_auto_group'
    factor_match = _AUTO_FACTOR_LABEL.match(normalized)
    if factor_match:
        return f"_auto_factor{factor_match.group(1)}"
    compact = normalized.replace('_', '').replace(' ', '').lower()
    factor_match = re.fullmatch(r'autofactor(\d+)', compact)
    if factor_match:
        return f"_auto_factor{factor_match.group(1)}"
    return normalized


def _resolve_pca_groups(adata, requested='', auto_group_mapping=None, *, role='color'):
    """Resolve color/marker values to real obs columns with an audit record."""
    from modules.io_utils import materialize_auto_sample_metadata

    requested = str(requested or '').strip()
    metadata = materialize_auto_sample_metadata(adata, auto_group_mapping)
    used = _auto_column_from_requested(requested)
    warning = ''

    if not requested:
        return None, {
            'requested': '', 'used': '', 'source': 'none', 'group_counts': {},
            'warning': '', 'display_label': '',
        }

    if used in adata.obs.columns:
        values = adata.obs[used].astype(str).tolist()
        if len(set(values)) >= 2:
            if used == '_auto_group':
                source = metadata.get('group_source') or 'sample_name_inference'
            elif used.startswith('_auto_factor'):
                source = 'sample_name_factor_inference'
            else:
                source = 'obs'
            return values, {
                'requested': requested, 'used': used, 'source': source,
                'group_counts': {str(k): int(v) for k, v in pd.Series(values).value_counts().items()},
                'warning': '', 'display_label': used,
            }
        warning = f"obs 列 '{used}' 只有一个分组，未用于 PCA {role}。"
    else:
        warning = f"未找到 obs 列 '{used}'，未用于 PCA {role}。"

    # A colour request may recover to the inferred combined group.  Marker
    # requests deliberately do not fall back: a wrong shape encoding is less
    # visible than a missing one and therefore harder to audit.
    if role == 'color' and used != '_auto_group' and '_auto_group' in adata.obs.columns:
        values = adata.obs['_auto_group'].astype(str).tolist()
        if len(set(values)) >= 2:
            return values, {
                'requested': requested, 'used': '_auto_group',
                'source': metadata.get('group_source') or 'sample_name_inference',
                'group_counts': {str(k): int(v) for k, v in pd.Series(values).value_counts().items()},
                'warning': warning, 'display_label': '_auto_group',
            }

    if role == 'color':
        if not warning:
            warning = '样本名中未识别到可重复的分组模式，PCA 将以统一颜色显示。'
        return None, {
            'requested': requested or '_auto_group_', 'used': '', 'source': 'none',
            'group_counts': {}, 'warning': warning, 'display_label': '',
        }
    return None, {
        'requested': requested, 'used': '', 'source': 'none',
        'group_counts': {}, 'warning': warning, 'display_label': '',
    }


def _resolve_pca_color_groups(adata, requested_color_by='', auto_group_mapping=None):
    """Backward-compatible color resolver used by existing callers/tests."""
    return _resolve_pca_groups(
        adata, requested_color_by, auto_group_mapping, role='color',
    )


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
        from modules.io_utils import (
            read_expression_matrix, resolve_expression_measurement, run_bulk_sample_pca,
        )
        adata = read_expression_matrix(input_path)
        input_measurement, measurement_info = resolve_expression_measurement(
            adata, input_path, self.params.get('input_measurement', 'auto'),
        )
        adata.uns['input_measurement'] = input_measurement
        adata.uns['input_measurement_provenance'] = measurement_info

        # 清理 inf/NaN
        import numpy as _np
        adata.X = _np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

        n_comps = int(self.params.get('n_comps', 10))
        color_by = self.params.get('color_by', '')
        batch_by = self.params.get('batch_by', '')
        dimred_method = self.params.get('dimred_method', 'pca')

        self.progress(20, "标准化数据...")
        if 'normalization' not in adata.uns:
            if input_measurement == 'raw_counts':
                sc.pp.normalize_total(adata, target_sum=1e6)
                sc.pp.log1p(adata)
                adata.uns['normalization'] = {'method': 'pca_auto_cpm_log1p', 'is_log_transformed': True}
            elif input_measurement == 'continuous_expression':
                if measurement_info['requires_confirmation_for_log_transform']:
                    raise ValueError(
                        '自动检测到非整数连续值，但无法仅靠数值区分线性 FPKM/TPM 与已 log 的表达矩阵。'
                        '请在“输入表达量尺度”中明确选择 continuous_expression 或 log_transformed 后重试。'
                    )
                adata.X = np.log2(np.maximum(adata.X, 0) + 1)
                adata.uns['normalization'] = {'method': 'pca_auto_log2', 'is_log_transformed': True}
            else:
                # Explicitly recognised log-expression: preserve it.  This
                # avoids direct PCA silently applying log2 a second time.
                adata.uns['normalization'] = {'method': 'pca_preserve_log_input', 'is_log_transformed': True}

        self.progress(40, "运行 PCA...")
        actual_comps = min(n_comps, adata.n_obs - 1, adata.n_vars - 1)
        if actual_comps < 2:
            raise ValueError(f"样本数不足 ({adata.n_obs})，至少需要 3 个样本才能进行 PCA 分析。")
        pca_preprocessing = run_bulk_sample_pca(adata, n_comps=actual_comps)
        pca_variance = adata.uns['pca']['variance_ratio']

        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        color_values, color_info = _resolve_pca_groups(
            adata, color_by, self.params.get('_auto_group_mapping', {}), role='color')
        batch_values, batch_info = _resolve_pca_groups(
            adata, batch_by, self.params.get('_auto_group_mapping', {}), role='marker')
        if color_info['warning']:
            self.progress(-1, color_info['warning'])
        if batch_info['warning']:
            self.progress(-1, batch_info['warning'])

        # A missing grouping is a valid exploratory PCA state.  The renderer
        # still receives a vector of real values so later legend and CSV logic
        # cannot fail on ``None``.
        if color_values is None:
            color_values = ['All'] * adata.n_obs

        # Sample names are useful only in small static panels.  Larger PCA
        # panels use the external legend so labels cannot cover nearby points.
        show_text = adata.n_obs <= 8

        self.progress(55, "生成 PCA 图...")
        pc = adata.obsm['X_pca']
        hover = adata.obs.index.tolist()
        from modules.bulk_qc import _detect_outliers_mahal
        outlier_samples, outlier_details = _detect_outliers_mahal(
            pc, hover, return_details=True,
        )

        # The publication figure is static by design.  Exporting a complete
        # sample-level score table keeps every point identifiable, including
        # PC3+ outliers that may not be visible in the PC1/PC2 panel.
        score_columns = [f'PC{i + 1}' for i in range(actual_comps)]
        pca_scores = pd.DataFrame(pc[:, :actual_comps], columns=score_columns)
        pca_scores.insert(0, 'sample_id', [str(item) for item in hover])
        pca_scores['color_group'] = [str(item) for item in color_values]
        pca_scores['marker_group'] = (
            [str(item) for item in batch_values]
            if batch_values is not None else ''
        )
        for number in (1, 2):
            column = f'_auto_factor{number}'
            pca_scores[f'factor{number}'] = (
                adata.obs[column].astype(str).tolist() if column in adata.obs.columns else ''
            )
        pca_scores['pca_outlier_distance'] = outlier_details['distances']
        pca_scores['pca_outlier_threshold'] = outlier_details['threshold']
        pca_scores['is_pca_outlier'] = pca_scores['sample_id'].isin(outlier_samples)
        pca_scores_path = os.path.join(results_dir, 'bulk_pca_scores.csv')
        pca_scores.to_csv(pca_scores_path, index=False)
        result_files.append({
            'file_path': pca_scores_path, 'file_type': 'csv', 'category': 'table',
            'label': 'PCA sample scores and metadata',
        })
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

        # The primary PCA is rendered through the deterministic publication
        # template.  Legacy helpers remain below for exploratory PC1/PC3,
        # t-SNE and UMAP panels so this Phase 1 change stays narrowly scoped.
        import matplotlib.pyplot as plt
        from figure_engine import NatureFigureDirector, export_registered_figure

        director = NatureFigureDirector()
        nature_formats = ('svg', 'pdf', 'png')
        unique_color_values = list(dict.fromkeys(str(value) for value in color_values))
        dense_grouping = len(unique_color_values) > 8
        dense_group_colors = ({group: '#4C78A8' for group in unique_color_values}
                              if dense_grouping else None)

        def _export_pca(fig, stem, label, spec, *, category='pca'):
            spec = spec.with_updates(formats=nature_formats)
            exported, report = export_registered_figure(
                fig, os.path.join(plots_dir, stem), spec,
                category=category, label=label,
                qa_path=os.path.join(results_dir, f'{stem}_nature_readiness.json'),
            )
            if not report.ready:
                self.progress(-1, f'{label} Nature readiness {report.score}/100；请查看 QA 报告。')
            plt.close(fig)
            return exported
        pca_spec = director.spec_from_params(
            'pca', self.params,
            title=f'PCA analysis (n={adata.n_obs})',
            show_legend=not dense_grouping,
        ).with_updates(
            # Small static panels can carry every sample ID; for larger panels
            # label only diagnosed outliers and use the score CSV for lookup.
            show_sample_labels=bool(self.params.get('show_sample_labels', False)) or adata.n_obs <= 12,
            outlier_labels=tuple(outlier_samples),
        )
        if dense_grouping and batch_values is None:
            batch_values = color_values
            batch_info = {
                **batch_info, 'used': color_info['used'], 'source': color_info['source'],
                'display_label': color_info['display_label'] or 'color_group',
            }
        fig_pca = director.render(pca_spec, {
            'coordinates': pc[:, :2],
            'groups': color_values,
            'batches': batch_values,
            'group_colors': dense_group_colors,
            'samples': hover,
            'explained_variance': pca_variance[:2],
            'group_label': color_info['display_label'] or 'Color group',
            'batch_label': batch_info['display_label'] or 'Marker group',
        })
        exported, readiness = export_registered_figure(
            fig_pca, os.path.join(plots_dir, 'bulk_pca'), pca_spec,
            category='pca', label='PCA 分析',
            qa_path=os.path.join(results_dir, 'bulk_pca_nature_readiness.json'),
        )
        result_files.extend(exported)
        if not readiness.ready:
            self.progress(-1, f'PCA Nature readiness {readiness.score}/100；请查看 QA 报告。')
        plt.close(fig_pca)

        # PC3 often carries a biologically meaningful secondary separation when
        # PC1+PC2 explain a moderate fraction of the total transcriptome variance.
        if actual_comps >= 3:
            for first, second, stem in [(0, 2, 'bulk_pca_pc1_pc3'), (1, 2, 'bulk_pca_pc2_pc3')]:
                alt_spec = director.spec_from_params(
                    'pca', self.params, title=f'PCA: PC{first + 1} vs PC{second + 1}',
                    show_legend=not dense_grouping,
                ).with_updates(
                    formats=nature_formats, height_mm=82.0,
                    show_sample_labels=bool(self.params.get('show_sample_labels', False)) or adata.n_obs <= 12,
                    outlier_labels=tuple(outlier_samples),
                )
                fig_alt = director.render(alt_spec, {
                    'coordinates': pc[:, [first, second]], 'groups': color_values,
                    'batches': batch_values, 'group_colors': dense_group_colors,
                    'samples': hover, 'explained_variance': pca_variance[[first, second]],
                    'x_label': f'PC{first + 1}', 'y_label': f'PC{second + 1}',
                    'group_label': color_info['display_label'] or 'Color group',
                    'batch_label': batch_info['display_label'] or 'Marker group',
                })
                result_files.extend(_export_pca(
                    fig_alt, stem, f'PCA：PC{first + 1} vs PC{second + 1}', alt_spec,
                ))

        self.progress(65, "生成方差解释图...")
        fig_var = bar_figure(
            [f'PC{i+1}' for i in range(actual_comps)], pca_variance * 100,
            title='PCA 方差解释比例', x_label='主成分',
            y_label='方差解释比例 (%)', rotation=45,
        )
        variance_spec = director.spec_from_params(
            'diagnostic', self.params, width='single', title='PCA variance explained',
        ).with_updates(extra={'kind': 'variance'}, formats=nature_formats, height_mm=70.0)
        fig_var_native = director.render(variance_spec, {
            'kind': 'variance', 'variance': pca_variance,
            'labels': [f'PC{i + 1}' for i in range(actual_comps)],
        })
        result_files.extend(_export_pca(fig_var_native, 'bulk_pca_variance', '方差解释比例', variance_spec))

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
                load_spec = director.spec_from_params(
                    'diagnostic', self.params, width='double', title=f'{pc_name} top loadings',
                ).with_updates(extra={'kind': 'bar'}, formats=nature_formats, height_mm=78.0)
                fig_load_native = director.render(load_spec, {
                    'kind': 'bar', 'labels': [adata.var_names[i] for i in top_idx],
                    'values': values, 'ylabel': 'Loading',
                })
                result_files.extend(_export_pca(
                    fig_load_native, f'bulk_pca_loadings_{pc_name.lower()}', f'{pc_name} 载荷图', load_spec,
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
        elbow_spec = director.spec_from_params(
            'diagnostic', self.params, width='single', title='PCA cumulative variance',
        ).with_updates(extra={'kind': 'variance'}, formats=nature_formats, height_mm=70.0)
        fig_elbow_native = director.render(elbow_spec, {
            'kind': 'variance', 'variance': pca_variance,
            'labels': [f'PC{i + 1}' for i in range(actual_comps)],
        })
        result_files.extend(_export_pca(fig_elbow_native, 'bulk_pca_elbow', '肘部图', elbow_spec))

        # 降维方法选择：t-SNE / UMAP / PCA-only
        if dimred_method == 'tsne' and adata.n_obs >= 3:
            from sklearn.manifold import TSNE
            self.progress(75, "运行 t-SNE...")
            perplexity = min(30, adata.n_obs - 1)
            tsne = TSNE(n_components=2, random_state=42, perplexity=max(2, perplexity))
            tsne_coords = tsne.fit_transform(adata.obsm['X_pca'])

            tsne_spec = director.spec_from_params(
                'pca', self.params, title=f't-SNE (n={adata.n_obs})', show_legend=not dense_grouping,
            ).with_updates(formats=nature_formats, height_mm=82.0)
            fig_tsne = director.render(tsne_spec, {
                'coordinates': tsne_coords, 'groups': color_values, 'batches': batch_values,
                'group_colors': dense_group_colors, 'samples': hover,
                'explained_variance': [np.nan, np.nan], 'x_label': 't-SNE1', 'y_label': 't-SNE2',
                'group_label': color_info['display_label'] or 'Color group',
                'batch_label': batch_info['display_label'] or 'Marker group',
            })
            result_files.extend(_export_pca(fig_tsne, 'bulk_tsne', 't-SNE 分析', tsne_spec, category='tsne'))

        elif dimred_method == 'umap' and adata.n_obs >= 10:
            self.progress(75, "运行 UMAP...")
            sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1))
            sc.tl.umap(adata)
            umap_coords = adata.obsm['X_umap']

            umap_spec = director.spec_from_params(
                'pca', self.params, title='UMAP', show_legend=not dense_grouping,
            ).with_updates(formats=nature_formats, height_mm=82.0)
            fig_umap = director.render(umap_spec, {
                'coordinates': umap_coords, 'groups': color_values, 'batches': batch_values,
                'group_colors': dense_group_colors, 'samples': hover,
                'explained_variance': [np.nan, np.nan], 'x_label': 'UMAP1', 'y_label': 'UMAP2',
                'group_label': color_info['display_label'] or 'Color group',
                'batch_label': batch_info['display_label'] or 'Marker group',
            })
            result_files.extend(_export_pca(fig_umap, 'bulk_umap', 'UMAP 分析', umap_spec, category='umap'))

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
                'input_measurement': input_measurement,
                'input_measurement_source': measurement_info['source'],
                'input_measurement_confidence': measurement_info['confidence'],
                'pca_preprocessing': pca_preprocessing,
                'color_by_requested': color_info['requested'],
                'color_by_used': color_info['used'],
                'color_by_source': color_info['source'],
                'color_by_warning': color_info['warning'],
                'batch_by_requested': batch_info['requested'],
                'batch_by_used': batch_info['used'],
                'batch_by_source': batch_info['source'],
                'batch_by_warning': batch_info['warning'],
                # Kept for downstream consumers of existing task summaries.
                'grouping_source': color_info['source'],
                'group_counts': color_info['group_counts'],
                'grouping_warning': color_info['warning'],
                'marker_group_counts': batch_info['group_counts'],
                'outlier_samples': outlier_samples,
                'outlier_pcs_used': outlier_details['n_components'],
                'outlier_threshold': outlier_details['threshold'],
                'dimred_method': dimred_method,
            }
        }
