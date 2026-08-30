#!/usr/bin/env python
"""CellOracle 子进程执行器（虚拟敲除 / virtual KO）。

本脚本在独立的 celloracle conda/venv 环境（Python 3.9/3.10）中运行，
由平台模块 modules/virtual_ko.py 通过 subprocess 调用。平台主环境
（Python 3.12）无法 import celloracle，因此所有 CellOracle 计算都
封装在这里。

协议：
  - stdout 输出 JSON 行事件：
      {"type": "progress", "pct": <int>, "message": <str>}
      {"type": "result", "summary": <dict>, "files": [<dict>, ...]}
      {"type": "error", "message": <str>, "traceback": <str>}
  - 退出码：0 成功；非 0 失败。
  - 所有输出文件写入 config 指定的目录，路径均为绝对路径。

调用方式：
  <celloracle_env_python> celloracle_worker.py --config <job.json>
"""
import argparse
import json
import os
import sys
import traceback

# 本脚本位于平台的 modules/ 目录下，而该目录内有一个 platform/ 子包，
# 会遮蔽标准库 platform 模块（matplotlib -> uuid -> platform 会因此崩溃）。
# 把脚本所在目录从 sys.path 移除，确保标准库优先。
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if sys.path and os.path.abspath(sys.path[0]) == _SCRIPT_DIR:
    sys.path.pop(0)


def emit(event_type, **payload):
    line = json.dumps({"type": event_type, **payload},
                      ensure_ascii=False, default=str)
    sys.stdout.write(line + "\n")
    sys.stdout.flush()


def progress(pct, message):
    emit("progress", pct=max(0, min(100, int(pct))), message=str(message))


def fail(message):
    emit("error", message=str(message), traceback=traceback.format_exc())
    sys.exit(1)


def _auto_n_comps(explained_ratio, max_comps=50):
    """按 CellOracle 教程的肘部启发式自动选择 PCA 主成分数。"""
    import numpy as np
    cumsum = np.cumsum(explained_ratio)
    if len(cumsum) < 5:
        return min(len(cumsum), max_comps)
    idx = np.where(np.diff(np.diff(cumsum)) > 0.002)[0]
    if len(idx) == 0:
        n = min(20, len(cumsum))
    else:
        n = int(idx[0])
    return max(3, min(n, max_comps))


def _parse_genes(text):
    """解析基因列表（逗号/分号/换行/空格分隔）。"""
    if not text:
        return []
    import re
    parts = re.split(r"[,\s;，；]+", str(text))
    return [p.strip() for p in parts if p.strip()]


def _load_base_grn(source, path, version):
    """加载 base GRN，返回 (tfdict, meta, tf_info_df)。

    支持：
    - 内置人类 promoter base GRN（CellOracle 官方 hg19/hg38）
    - 用户上传的 TF_info matrix（parquet/csv）
    - 用户上传的 TFdict 字典 pickle
    - 用户上传的完整 Oracle 对象 pickle（取其中的 TFdict）
    """
    import pandas as pd
    import celloracle as co

    if source.startswith("builtin_human"):
        version = version or (
            "hg38_gimmemotifsv5_fpr2" if source == "builtin_human_hg38"
            else "hg19_gimmemotifsv5_fpr2")
        progress(4, "加载内置人类 promoter base GRN（%s），首次使用需要下载..." % version)
        df = co.data.load_human_promoter_base_GRN(version=version)
        return None, {"kind": "TF_info_matrix", "source": source,
                      "n_rows": int(df.shape[0]), "n_cols": int(df.shape[1])}, df

    if not path or not os.path.isfile(path):
        fail("base GRN 文件不存在: " + str(path))
    ext = os.path.splitext(path)[1].lower()
    obj = None
    pickle_warning = None
    if ext in (".parquet", ".pq"):
        obj = pd.read_parquet(path)
    elif ext in (".csv", ".tsv", ".txt", ".gz"):
        sep = "\t" if ext in (".tsv", ".txt") else ","
        obj = pd.read_csv(path, sep=sep)
    elif ext in (".pickle", ".pkl", ".gpickle", ".oracle", ".celloracle", ".links"):
        # 安全提示：pickle 反序列化可执行任意代码，仅应加载可信来源文件。
        # 平台限制这些文件只能来自项目 uploads 目录，但仍需用户知晓风险。
        pickle_warning = ("base GRN 为 pickle 格式；pickle 反序列化存在代码执行风险，"
                          "请确认文件来源可信。")
        obj = co.utility.load_pickled_object(path)
    else:
        fail("不支持的 base GRN 文件格式: " + ext)

    def _meta(**extra):
        base = {"source": path}
        if pickle_warning:
            base["pickle_security_warning"] = pickle_warning
        base.update(extra)
        return base

    # 完整 Oracle 对象：提取 TFdict
    if hasattr(obj, "TFdict") and getattr(obj, "TFdict"):
        return dict(obj.TFdict), _meta(kind="oracle_object"), None
    # 字典：{TF: [target_genes]}
    if isinstance(obj, dict):
        return dict(obj), _meta(kind="TFdict", n_tfs=len(obj)), None
    # DataFrame：TF_info matrix（含 peak_id / gene_short_name 列）
    if hasattr(obj, "columns") and "gene_short_name" in obj.columns:
        return None, _meta(kind="TF_info_matrix", n_rows=int(obj.shape[0])), obj
    fail("无法识别的 base GRN 对象: %s；请上传 CellOracle TF_info matrix（parquet/csv）、TFdict 字典 pickle 或包含 TFdict 的 Oracle 对象 pickle。" % type(obj).__name__)


def _load_links(path):
    """加载可选的 co-accessibility links 对象。"""
    import celloracle as co
    if not path or not os.path.isfile(path):
        return None
    obj = co.utility.load_pickled_object(path)
    if not hasattr(obj, "filter_links"):
        fail("links 文件不是有效的 CellOracle Links 对象")
    obj.filter_links()
    return obj


def run(config):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 100, "savefig.dpi": 300,
        "font.size": 10, "axes.titlesize": 12,
        "axes.labelsize": 10, "savefig.bbox": "tight",
    })
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import celloracle as co

    out_dir = os.path.abspath(config["output_dir"])
    plots_dir = os.path.abspath(config["plots_dir"])
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    result_files = []
    summary = {
        "module": "virtual_ko",
        "celloracle_version": getattr(co, "__version__", "unknown"),
        "genes": [], "warnings": [],
    }

    # ---------- 1. 数据加载 ----------
    progress(2, "加载输入 h5ad...")
    input_h5ad = config["input_h5ad"]
    adata = sc.read_h5ad(input_h5ad)
    if config.get("counts_layer") and config["counts_layer"] in adata.layers:
        adata.X = adata.layers[config["counts_layer"]].copy()
    else:
        candidates = [l for l in ("counts", "raw_count", "counts_raw")
                      if l in adata.layers]
        if candidates:
            adata.X = adata.layers[candidates[0]].copy()
    # 硬校验：CellOracle 的 import_anndata_as_raw_count 会把 X 当作 UMI
    # 计数建模。若没有原始 counts 层而回退到 log1p 归一化 X，GRN/模拟
    # 会在错误尺度上得到看似合理的结论。默认直接失败，除非用户在界面
    # 显式勾选 allow_non_count_fallback 承认风险。
    _allow_fallback = bool(config.get("allow_non_count_fallback", False))
    try:
        import numpy as np
        from scipy import sparse
        _x = adata.X
        if sparse.issparse(_x):
            _values = np.asarray(_x.data, dtype=float)
        else:
            _values = np.asarray(_x, dtype=float).reshape(-1)
        if _values.size > 200000:
            _values = _values[:: max(1, _values.size // 200000)][:200000]
        _values = _values[np.isfinite(_values)]
        _is_raw = (
            _values.size > 0
            and float(np.mean(_values >= 0)) >= 0.999
            and float(np.mean(np.isclose(_values, np.rint(_values)))) >= 0.995
        )
    except Exception:
        _is_raw = False
    if adata.X.min() < 0:
        fail("表达矩阵包含负值，无法作为原始 counts 输入 CellOracle。请确认输入 h5ad 保留了 counts layer。")
    if not _is_raw:
        if not _allow_fallback:
            fail(
                "未找到非负整数原始 counts 层（X 可能是 log1p 归一化或残差）。"
                "CellOracle 必须以原始 UMI counts 建模；请重新从 QC 保留 counts layer 的 h5ad 运行，"
                "或在界面显式勾选“允许非 counts 回退”（不推荐）。"
            )
        summary["warnings"].append(
            "未找到原始 counts 层，已按用户显式确认回退使用当前 X 尺度；"
            "GRN 与敲除模拟的数值解读应视为不可靠。"
        )
    # 平台 normalize 模块会在 uns 写入 'log1p' 标记；这里 X 已替换为原始
    # counts，需要清除该标记，避免 CellOracle 内部 sc.pp.log1p 产生
    # “already log-transformed”的误导性警告。
    adata.uns.pop('log1p', None)

    cluster_key = config.get("cluster_key") or ""
    if not cluster_key or cluster_key not in adata.obs.columns:
        for cand in ("celltype", "final_annotation", "cell_type",
                     "leiden", "louvain_annot"):
            if cand in adata.obs.columns:
                cluster_key = cand
                summary["warnings"].append("cluster_key 未指定或不存在，自动使用列 '%s'" % cand)
                break
    if not cluster_key or cluster_key not in adata.obs.columns:
        fail("缺少聚类/注释列。请先运行 clustering 或 annotation 模块。")
    summary["cluster_key"] = cluster_key

    embedding_name = config.get("embedding_name") or "X_umap"
    if embedding_name not in adata.obsm:
        embedding_name = next(
            (k for k in ("X_umap", "X_draw_graph_fa", "X_tsne") if k in adata.obsm),
            None)
        if not embedding_name:
            fail("缺少 2D embedding（obsm 中没有 X_umap 等键）。")
        summary["warnings"].append("自动使用 embedding '%s'" % embedding_name)
    summary["embedding_name"] = embedding_name

    max_cells = int(config.get("max_cells") or 0)
    if max_cells > 0 and adata.n_obs > max_cells:
        progress(3, "随机下采样 %d 个细胞（共 %d）..." % (max_cells, adata.n_obs))
        sc.pp.subsample(adata, n_obs=max_cells, random_state=123)

    # ---------- 2. base GRN ----------
    tfdict, meta, tf_info_df = _load_base_grn(
        config.get("base_grn_source", "upload"),
        config.get("base_grn_file", ""),
        config.get("base_grn_version", ""),
    )
    if meta and meta.get("pickle_security_warning"):
        summary["warnings"].append(meta["pickle_security_warning"])
    summary["base_grn"] = meta
    if tfdict is not None:
        tfs = set(tfdict.keys())
    else:
        tfs = set(tf_info_df.columns) - {"peak_id", "gene_short_name"}
    tfs_in_data = sorted(tfs & set(adata.var_names))
    summary["n_tfs_total"] = len(tfs)
    summary["n_tfs_in_data"] = len(tfs_in_data)
    if not tfs_in_data:
        fail("base GRN 中的 TF 与 scRNA-seq 数据没有任何交集，请检查物种/基因名是否一致。")

    # ---------- 3. 基因筛选（HVG ∪ TF ∩ data ∪ 扰动基因）----------
    if config.get("use_hvg", True) and adata.n_vars > 4000:
        hv_flag = None
        for col in ("highly_variable", "highly_variable_features"):
            if col in adata.var.columns:
                hv_flag = col
                break
        if hv_flag is not None:
            perturb_genes = _parse_genes(config.get("perturb_genes", ""))
            keep = set(adata.var_names[adata.var[hv_flag].astype(bool)])
            keep |= set(tfs_in_data)
            keep |= {g for g in perturb_genes if g in adata.var_names}
            keep = sorted(keep & set(adata.var_names))
            adata = adata[:, keep].copy()
            summary["warnings"].append(
                "基因数超过 4000，已按 HVG+TF 收敛为 %d 个基因" % adata.n_vars)

    # ---------- 4. Oracle 对象 ----------
    progress(6, "构建 Oracle 对象...")
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata, cluster_column_name=cluster_key,
        embedding_name=embedding_name)

    if tfdict is not None:
        oracle.import_TF_data(TFdict=tfdict)
    else:
        oracle.import_TF_data(TF_info_matrix=tf_info_df)

    links = _load_links(config.get("links_file", ""))
    use_cluster_specific = links is not None

    # ---------- 5. PCA + KNN imputation ----------
    progress(8, "PCA...")
    oracle.perform_PCA()
    n_comps = int(config.get("n_comps") or 0)
    if n_comps <= 0:
        n_comps = _auto_n_comps(oracle.pca.explained_variance_ratio_,
                                max_comps=int(config.get("n_comps_max") or 50))
    summary["n_comps"] = n_comps

    n_cell = oracle.adata.n_obs
    k = int(config.get("knn_k") or 0)
    if k <= 0:
        k = max(5, int(0.025 * n_cell))
    k = min(k, n_cell - 1)
    summary["knn_k"] = k

    progress(10, "KNN imputation（k=%d, n_pca_dims=%d）..." % (k, n_comps))
    n_jobs = int(config.get("n_jobs") or 4)
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True,
                          b_sight=k * 8, b_maxl=k * 4, n_jobs=n_jobs)

    # ---------- 6. GRN 推断 ----------
    grn_unit = config.get("grn_unit") or "cluster"
    alpha = float(config.get("alpha") or 10)
    progress(40, "GRN 推断（%s, alpha=%s）..." % (grn_unit, alpha))
    if use_cluster_specific:
        oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(GRN_unit=grn_unit, alpha=alpha,
                                  use_cluster_specific_TFdict=use_cluster_specific,
                                  verbose_level=0)
    summary["alpha"] = alpha
    summary["grn_unit"] = grn_unit

    # 导出 GRN 边表
    edge_cutoff = float(config.get("edge_coef_cutoff") or 0)
    top_edges = int(config.get("top_edges_per_cluster") or 500)
    edge_rows = []
    if hasattr(oracle, "coef_matrix_per_cluster"):
        coef_maps = oracle.coef_matrix_per_cluster
    elif hasattr(oracle, "coef_matrix"):
        coef_maps = {"whole": oracle.coef_matrix}
    else:
        fail("GRN 推断未产生 coef_matrix。")
    n_edges_total = 0
    for cluster, cm in coef_maps.items():
        long_df = cm.stack().reset_index()
        long_df.columns = ["target", "regulator", "coef"]
        long_df["abs_coef"] = long_df["coef"].abs()
        if edge_cutoff > 0:
            long_df = long_df[long_df["abs_coef"] >= edge_cutoff]
        long_df = long_df.sort_values("abs_coef", ascending=False).head(top_edges)
        long_df.insert(0, "cluster", str(cluster))
        n_edges_total += len(long_df)
        edge_rows.append(long_df)
    if edge_rows:
        edges_df = pd.concat(edge_rows, ignore_index=True)
        edges_path = os.path.join(out_dir, "virtual_ko_grn_edges.csv")
        edges_df.to_csv(edges_path, index=False)
        result_files.append({
            "path": edges_path, "kind": "csv", "category": "table",
            "label": "GRN 推断边表（每簇 Top 边）",
            "description": "alpha=%s, |coef|>=%s, 每簇最多 %d 条边，共 %d 条" % (alpha, edge_cutoff, top_edges, n_edges_total),
        })
        summary["n_grn_edges"] = n_edges_total
    else:
        summary["warnings"].append("没有导出任何 GRN 边（coef 全为 0 或阈值过高）")

    # ---------- 7. 扰动基因 ----------
    perturb_genes = _parse_genes(config.get("perturb_genes", ""))
    missing = [g for g in perturb_genes if g not in oracle.adata.var_names]
    if missing:
        fail("以下扰动基因不在表达矩阵中: " + ", ".join(missing))
    not_tf = [g for g in perturb_genes if g not in tfs]
    if not_tf:
        summary["warnings"].append(
            "以下基因不是 base GRN 中的 TF（仍会按目标基因方式模拟）: " + ", ".join(not_tf))
    if not perturb_genes:
        fail("请至少填写一个扰动基因。")

    n_propagation = int(config.get("n_propagation") or 3)
    n_neighbors = int(config.get("n_neighbors") or 200)
    min_mass = float(config.get("min_mass") or 0.01)
    n_grid = int(config.get("n_grid") or 40)
    top_regulated = int(config.get("top_regulated_genes") or 50)
    jobs = list(perturb_genes)
    if config.get("combine_perturbations") and len(perturb_genes) > 1:
        jobs.append("__combined__")

    per_gene_summaries = {}
    for job_i, gene in enumerate(jobs):
        label = "组合敲除" if gene == "__combined__" else gene
        base_pct = 55 + int(30 * job_i / max(len(jobs), 1))
        progress(base_pct, "模拟 %s 扰动..." % label)
        condition = ({g: 0.0 for g in perturb_genes}
                     if gene == "__combined__"
                     else {gene: 0.0})
        oracle.simulate_shift(perturb_condition=condition,
                              n_propagation=n_propagation)
        # 固定随机种子，保证随机对照（delta_embedding_random）可复现。
        np.random.seed(int(config.get("random_seed", 0)))
        oracle.estimate_transition_prob(n_neighbors=n_neighbors,
                                        knn_random=True,
                                        sampled_fraction=1)
        oracle.calculate_embedding_shift(sigma_corr=0.05)
        # p_mass 的邻居数与 transition prob 保持一致，而不是硬编码 200。
        oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=n_neighbors)
        oracle.calculate_mass_filter(min_mass=min_mass, plot=False)

        adata_o = oracle.adata
        simulated = adata_o.layers["simulated_count"]
        imputed = adata_o.layers["imputed_count"]
        delta_x = simulated - imputed
        delta_emb = np.asarray(oracle.delta_embedding)
        delta_emb_r = np.asarray(oracle.delta_embedding_random)
        delta_norm = np.linalg.norm(delta_emb, axis=1)
        delta_norm_r = np.linalg.norm(delta_emb_r, axis=1)

        # CellOracle 的 imputed/simulated 层是稀疏矩阵；np.asarray(csr)
        # 会得到 0 维 object 数组而不是数值矩阵，必须先 toarray()。
        def _dense_column(matrix, mask):
            sub = matrix[:, mask]
            return sub.toarray().ravel() if hasattr(sub, 'toarray') else np.asarray(sub).ravel()

        # 每个细胞状态表
        cell_df = pd.DataFrame({
            "cell": adata_o.obs_names,
            cluster_key: adata_o.obs[cluster_key].astype(str).values,
            "delta_embedding_x": delta_emb[:, 0],
            "delta_embedding_y": delta_emb[:, 1],
            "delta_embedding_norm": delta_norm,
            "delta_embedding_random_norm": delta_norm_r,
        })
        if gene != "__combined__":
            cell_df["imputed_%s" % gene] = _dense_column(imputed, adata_o.var_names == gene)
            cell_df["simulated_%s" % gene] = _dense_column(simulated, adata_o.var_names == gene)
        cell_csv = os.path.join(out_dir, "virtual_ko_%s_cell_state.csv" % label)
        cell_df.to_csv(cell_csv, index=False)
        result_files.append({
            "path": cell_csv, "kind": "csv", "category": "table",
            "label": "%s 每个细胞的状态偏移" % label,
        })

        # 受调控基因表（全基因 |delta| 均值排序，附方向列）
        if hasattr(delta_x, 'toarray'):
            delta_dense = delta_x.toarray()
        else:
            delta_dense = np.asarray(delta_x)
        mean_abs_delta = np.abs(delta_dense).mean(axis=0)
        mean_delta = delta_dense.mean(axis=0)
        order = np.argsort(mean_abs_delta)[::-1]
        gene_names_arr = np.asarray(adata_o.var_names).astype(str)
        top_idx = order[:top_regulated]
        perturbed_set = {str(gene)} if gene != "__combined__" else set()
        top_df = pd.DataFrame({
            "gene": gene_names_arr[top_idx],
            "mean_abs_delta_expression": mean_abs_delta[top_idx],
            # 正值为敲除后模拟表达高于 imputed 基线（上调），负值为下调。
            "mean_delta_expression": mean_delta[top_idx],
            "is_perturbed_gene": [
                "yes" if item in perturbed_set else "no"
                for item in gene_names_arr[top_idx]
            ],
        })
        reg_csv = os.path.join(out_dir, "virtual_ko_%s_regulated_genes.csv" % label)
        top_df.to_csv(reg_csv, index=False)
        result_files.append({
            "path": reg_csv, "kind": "csv", "category": "table",
            "label": "%s Top %d 受调控基因" % (label, top_regulated),
        })

        n_shifted = int((delta_norm > delta_norm_r.mean() + delta_norm_r.std()).sum())
        per_gene_summaries[gene] = {
            "label": label,
            "n_propagation": n_propagation,
            "n_cells": int(adata_o.n_obs),
            "mean_delta_embedding_norm": float(delta_norm.mean()),
            "mean_delta_embedding_random_norm": float(delta_norm_r.mean()),
            "frac_cells_shifted": float(n_shifted / max(len(delta_norm), 1)),
            "top_regulated_genes": top_df["gene"].head(10).tolist(),
        }

        # ---------- 图 ----------
        progress(base_pct + 10, "绘制 %s 结果图..." % label)
        scale = 25
        scale_sim = 0.5
        try:
            fig, ax = plt.subplots(1, 2, figsize=(13, 6))
            oracle.plot_quiver(scale=scale, ax=ax[0])
            ax[0].set_title("Simulated cell identity shift: %s KO" % label)
            oracle.plot_quiver_random(scale=scale, ax=ax[1])
            ax[1].set_title("Randomized simulation")
            _save_fig(fig, plots_dir, "virtual_ko_%s_quiver" % label,
                      result_files, "quiver 向量场（KO vs 随机对照）")

            fig, ax = plt.subplots(1, 2, figsize=(13, 6))
            oracle.plot_simulation_flow_on_grid(scale=scale_sim, ax=ax[0])
            ax[0].set_title("Simulation flow: %s KO" % label)
            oracle.plot_simulation_flow_random_on_grid(scale=scale_sim, ax=ax[1])
            ax[1].set_title("Randomized flow")
            _save_fig(fig, plots_dir, "virtual_ko_%s_flow_grid" % label,
                      result_files, "模拟流场网格（KO vs 随机对照）")

            fig, ax = plt.subplots(figsize=(8, 8))
            oracle.plot_cluster_whole(ax=ax, s=10)
            oracle.plot_simulation_flow_on_grid(scale=scale_sim, ax=ax,
                                                show_background=False)
            _save_fig(fig, plots_dir, "virtual_ko_%s_flow_clusters" % label,
                      result_files, "细胞分群 + 模拟流场")

            fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
            ax[0].hist(delta_norm_r, bins=50, alpha=0.55, color="gray",
                       label="randomized")
            ax[0].hist(delta_norm, bins=50, alpha=0.7, color="#e74c3c",
                       label=label)
            ax[0].set_xlabel("|delta embedding|")
            ax[0].set_ylabel("cell count")
            ax[0].legend(frameon=False)
            top10 = top_df.head(10).iloc[::-1]
            ax[1].barh(top10["gene"], top10["mean_abs_delta_expression"],
                       color="#4c72b0")
            ax[1].set_xlabel("mean |delta expression|")
            ax[1].set_title("Top 10 regulated genes (%s KO)" % label)
            _save_fig(fig, plots_dir, "virtual_ko_%s_delta_distribution" % label,
                      result_files, "偏移量分布 + Top 调控基因")
        except Exception as exc:
            summary["warnings"].append("%s 绘图失败: %s" % (label, exc))

        # 将模拟结果写回 adata（obs 列 + layers），供 h5ad 输出复用
        adata_o.obs["delta_embedding_norm_%s" % label] = delta_norm
        adata_o.obs["delta_embedding_random_norm_%s" % label] = delta_norm_r
        adata_o.layers["simulated_count_%s" % label] = simulated

    summary["genes"] = per_gene_summaries

    # ---------- 8. 输出 h5ad ----------
    progress(96, "保存输出 h5ad...")
    output_h5ad = os.path.join(config["intermediate_dir"],
                               "virtual_ko_output.h5ad")
    os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
    oracle.adata.write_h5ad(output_h5ad)
    summary["output_h5ad"] = output_h5ad

    if config.get("save_oracle_object"):
        oracle_path = os.path.join(out_dir, "virtual_ko.celloracle.oracle")
        oracle.to_hdf5(oracle_path)
        result_files.append({
            "path": oracle_path, "kind": "h5ad", "category": "intermediate",
            "label": "CellOracle Oracle 对象（可复用）",
        })
        # Oracle 对象实际是 CellOracle 自定义 HDF5 容器（扩展名 .oracle），
        # 平台文件类型清单里没有专用类型；用 h5ad 会误导下载/预览路由，
        # 因此声明为中间产物并注明格式。
        result_files[-1]["kind"] = "oracle"

    # ---------- 9. summary ----------
    progress(99, "汇总结果...")
    summary_path = os.path.join(out_dir, "virtual_ko_summary.json")
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2, default=str)
    result_files.append({
        "path": summary_path, "kind": "json", "category": "info",
        "label": "虚拟敲除汇总 JSON",
    })

    emit("result", summary=summary, files=result_files)
    progress(100, "虚拟敲除分析完成")


def _save_fig(fig, plots_dir, stem, result_files, label):
    import matplotlib.pyplot as plt
    png_path = os.path.join(plots_dir, "%s.png" % stem)
    svg_path = os.path.join(plots_dir, "%s.svg" % stem)
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.15)
    fig.savefig(svg_path, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    result_files.append({"path": png_path, "kind": "png", "category": "plot",
                         "label": label})
    result_files.append({"path": svg_path, "kind": "svg", "category": "plot",
                         "label": label + "（矢量）"})


def main():
    parser = argparse.ArgumentParser(description="CellOracle virtual KO worker")
    parser.add_argument("--config", required=True, help="job config JSON 路径")
    args = parser.parse_args()
    try:
        with open(args.config, encoding="utf-8") as fh:
            config = json.load(fh)
        run(config)
    except SystemExit:
        raise
    except Exception:
        emit("error", message=traceback.format_exc().splitlines()[-1],
             traceback=traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
