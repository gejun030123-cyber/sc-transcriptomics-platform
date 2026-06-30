"""anndata 结构检查工具 — 供 api.py 的 adata_info 和 data_info_full 共用。"""


def inspect_adata(adata, include_samples=True):
    """检查 anndata 对象结构，返回 JSON-serializable dict。

    Args:
        adata: anndata.AnnData 对象（通常以 backed='r' 模式打开）
        include_samples: 是否包含 var_names/obs_names 样本
    """
    import numpy as np

    def _dtype_str(series):
        return str(series.dtype)

    def _dtype_str_arr(arr):
        return str(arr.dtype)

    obs_info = []
    for col in adata.obs.columns:
        s = adata.obs[col]
        entry = {'name': col, 'dtype': _dtype_str(s)}
        if s.dtype == object or s.dtype.name == 'category':
            nuniq = s.nunique()
            entry['n_unique'] = int(nuniq)
            if nuniq <= 20:
                entry['values'] = [str(v) for v in sorted(s.dropna().unique().tolist())]
        obs_info.append(entry)

    var_info = []
    for col in adata.var.columns:
        s = adata.var[col]
        entry = {'name': col, 'dtype': _dtype_str(s)}
        if s.dtype == object or s.dtype.name == 'category':
            nuniq = s.nunique()
            entry['n_unique'] = int(nuniq)
            if nuniq <= 20:
                entry['values'] = [str(v) for v in sorted(s.dropna().unique().tolist())[:20]]
        var_info.append(entry)

    layers_keys = list(adata.layers.keys())

    obsm_info = []
    for k in adata.obsm.keys():
        arr = adata.obsm[k]
        obsm_info.append({'key': k, 'shape': list(arr.shape), 'dtype': _dtype_str_arr(arr)})

    varm_info = []
    for k in adata.varm.keys():
        arr = adata.varm[k]
        varm_info.append({'key': k, 'shape': list(arr.shape), 'dtype': _dtype_str_arr(arr)})

    obsp_keys = list(adata.obsp.keys())

    uns_info = []
    for k in adata.uns.keys():
        v = adata.uns[k]
        if isinstance(v, np.ndarray):
            uns_info.append({'key': k, 'type': 'ndarray', 'shape': list(v.shape), 'dtype': str(v.dtype)})
        elif isinstance(v, dict):
            uns_info.append({'key': k, 'type': 'dict', 'keys': list(v.keys())[:20]})
        else:
            uns_info.append({'key': k, 'type': type(v).__name__, 'preview': str(v)[:200]})

    info = {
        'obs_columns': obs_info,
        'var_columns': var_info,
        'layers': layers_keys,
        'obsm': obsm_info,
        'varm': varm_info,
        'obsp': obsp_keys,
        'uns': uns_info,
    }

    if include_samples:
        info['var_names_sample'] = [str(v) for v in adata.var_names[:min(20, adata.n_vars)]]
        info['obs_names_sample'] = [str(v) for v in adata.obs_names[:min(10, adata.n_obs)]]

    return info
