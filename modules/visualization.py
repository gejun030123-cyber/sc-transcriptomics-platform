import json
import numpy as np

def umap_scatter(adata, color_key=None, basis='X_umap', max_cells=50000, title=''):
    import plotly.graph_objects as go
    coords = adata.obsm[basis]
    if coords.shape[1] > 2:
        coords = coords[:, :2]
    idx = np.arange(adata.n_obs)
    if adata.n_obs > max_cells:
        idx = np.random.choice(adata.n_obs, max_cells, replace=False)
        coords = coords[idx]
    color_vals = None
    if color_key and color_key in adata.obs.columns:
        color_vals = adata.obs[color_key].values[idx]
    fig = go.Figure()
    if color_vals is not None and hasattr(color_vals, 'dtype') and hasattr(color_vals, 'cat'):
        for cat in color_vals.cat.categories:
            mask = np.array(color_vals == cat)
            fig.add_trace(go.Scattergl(
                x=coords[mask, 0], y=coords[mask, 1],
                mode='markers', name=str(cat),
                marker=dict(size=2, opacity=0.6),
            ))
    elif color_vals is not None:
        fig.add_trace(go.Scattergl(
            x=coords[:, 0], y=coords[:, 1],
            mode='markers',
            marker=dict(size=2, opacity=0.6, color=color_vals, colorscale='Viridis',
                       colorbar=dict(title=color_key)),
        ))
    else:
        fig.add_trace(go.Scattergl(
            x=coords[:, 0], y=coords[:, 1],
            mode='markers',
            marker=dict(size=2, opacity=0.5, color='#1a237e'),
        ))
    fig.update_layout(
        title=title, xaxis_title='UMAP-1', yaxis_title='UMAP-2',
        plot_bgcolor='white', width=700, height=500,
        margin=dict(l=40, r=40, t=40, b=40)
    )
    return json.loads(fig.to_json())

def violin_plot(adata, keys, groupby=None, title=''):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    n = len(keys)
    fig = make_subplots(rows=1, cols=n, subplot_titles=keys)
    for i, key in enumerate(keys, 1):
        if key not in adata.obs.columns:
            continue
        if groupby and groupby in adata.obs.columns:
            for cat in adata.obs[groupby].cat.categories:
                mask = adata.obs[groupby] == cat
                fig.add_trace(go.Violin(
                    y=adata.obs.loc[mask, key], name=str(cat),
                    box_visible=True, meanline_visible=True
                ), row=1, col=i)
        else:
            fig.add_trace(go.Violin(
                y=adata.obs[key], name=key,
                box_visible=True, meanline_visible=True
            ), row=1, col=i)
    fig.update_layout(title=title, showlegend=bool(groupby), height=400, width=250*n)
    return json.loads(fig.to_json())

def scatter_plot(x, y, color=None, xlabel='', ylabel='', title='', hover_text=None):
    import plotly.graph_objects as go
    fig = go.Figure()
    marker = dict(size=3, opacity=0.6)
    if color is not None:
        marker['color'] = color
        marker['colorscale'] = 'Viridis'
    fig.add_trace(go.Scattergl(x=x, y=y, mode='markers', marker=marker, text=hover_text))
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel,
                     plot_bgcolor='white', width=700, height=500)
    return json.loads(fig.to_json())

def bar_plot(x, y, xlabel='', ylabel='', title=''):
    import plotly.graph_objects as go
    fig = go.Figure(go.Bar(x=x, y=y, marker_color='#1a237e'))
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel,
                     plot_bgcolor='white', width=600, height=400)
    return json.loads(fig.to_json())

def compute_gene_variability(data, metric='var'):
    """计算每个基因的变异度量。data: (samples, genes)"""
    if metric == 'var':
        return np.var(data, axis=0)
    elif metric == 'mad':
        median = np.median(data, axis=0)
        return np.median(np.abs(data - median), axis=0)
    elif metric == 'cv':
        mean = np.mean(data, axis=0)
        std = np.std(data, axis=0)
        return std / (np.abs(mean) + 1e-10)
    elif metric == 'range':
        return np.max(data, axis=0) - np.min(data, axis=0)
    else:
        raise ValueError(f"未知变异度量: {metric}，支持: var/mad/cv/range")


def transform_heatmap_data(data, row_scaling='zscore', pseudocount=1,
                           winsorize='none', clip_range=(-3, 3), missing_value='ignore'):
    """对热图数据进行标准化和变换。data: (samples, genes)"""
    # 1. 缺失值处理
    if missing_value == 'mean_fill':
        col_means = np.nanmean(data, axis=0)
        for j in range(data.shape[1]):
            mask = np.isnan(data[:, j])
            if mask.any():
                data[mask, j] = col_means[j]
    elif missing_value == 'zero_fill':
        data = np.nan_to_num(data, nan=0.0)

    # 2. Winsorize
    if winsorize != 'none':
        if winsorize == 'custom':
            pct = 0.01
        else:
            pct = float(winsorize.replace('pct', '')) / 100
        from scipy.stats.mstats import winsorize as sp_winsorize
        winsorized = sp_winsorize(data, limits=[pct, pct], axis=0)
        data = winsorized.data if hasattr(winsorized, 'data') else np.array(winsorized)

    # 3. 行标准化
    if row_scaling == 'zscore':
        mean = np.mean(data, axis=0)
        std = np.std(data, axis=0) + 1e-10
        data = (data - mean) / std
    elif row_scaling == 'center':
        data = data - np.mean(data, axis=0)

    # 4. 截断
    if clip_range is not None:
        data = np.clip(data, clip_range[0], clip_range[1])

    return data


def cluster_heatmap(data, method='ward', metric='euclidean'):
    """层次聚类，返回排序后的行索引列表。data: (n_samples, n_features)"""
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist

    if data.shape[0] <= 1:
        return list(range(data.shape[0]))

    if metric == 'pearson':
        dist = pdist(data, metric=lambda u, v: 1 - np.corrcoef(u, v)[0, 1])
    elif metric == 'spearman':
        from scipy.stats import spearmanr
        dist = pdist(data, metric=lambda u, v: 1 - spearmanr(u, v).correlation)
    elif metric == 'cosine':
        dist = pdist(data, metric='cosine')
    else:
        dist = pdist(data, metric='euclidean')

    link = linkage(dist, method=method)
    dendro = dendrogram(link, no_plot=True)
    return dendro['leaves']


DEFAULT_PALETTE = ['#1a237e', '#e53935', '#4caf50', '#ff9800', '#9c27b0',
                   '#00bcd4', '#795548', '#607d8b', '#f44336', '#3f51b5']


def build_annotation_bar(obs, columns, sample_order=None, palette=None):
    """
    为多个注释列生成颜色映射数据。

    Returns:
        fig_data: dict {col_name: {'colors': [...], 'groups': [...], 'unique': [...], 'color_map': {...}}}
        unique_groups: dict {col_name: [unique_values]}
    """
    if palette is None:
        palette = DEFAULT_PALETTE
    if sample_order is None:
        sample_order = list(range(len(obs)))

    fig_data = {}
    unique_groups = {}
    for col in columns:
        if col not in obs.columns:
            continue
        values = [str(obs[col].iloc[i]) for i in sample_order]
        uniq = sorted(set(values))
        color_map = {g: palette[i % len(palette)] for i, g in enumerate(uniq)}
        fig_data[col] = {
            'colors': [color_map[v] for v in values],
            'groups': values,
            'unique': uniq,
            'color_map': color_map,
        }
        unique_groups[col] = uniq
    return fig_data, unique_groups


def save_plotly_json(fig, plots_dir, filename, result_files,
                     file_type='plotly_json', category='heatmap', label=''):
    """保存 Plotly 图表为 JSON 并追加到 result_files 列表。"""
    import os
    fpath = os.path.join(plots_dir, filename)
    with open(fpath, 'w') as f:
        f.write(fig.to_json(engine="json"))
    result_files.append({
        'file_path': fpath, 'file_type': file_type,
        'category': category, 'label': label,
    })
