# SCENIC regulon 活性结果模块

## 目的与边界

`scenic` 模块补充的是 SCENIC **regulon 活性结果的本地整理与可视化**，而不是用普通 TF
表达量冒充调控活性。它读取已经在受控环境中完成的 pySCENIC/SCENIC+ AUCell 结果，因而能
稳定输出论文常用的 AUCell、RSS、regulon 网络和共活性模块图。

首期不在 Web worker 内重新跑 GRNBoost/GENIE3、cisTarget motif enrichment 或 AUCell：这些步骤
依赖版本固定的 TF 列表、ranking/motif 数据库、较大内存和可追溯的随机/并行设置。需要运行
完整 SCENIC 时，应在服务器受控环境预计算后再导入 H5AD；数据库版本、物种、TF 列表和输入
checksum 应和 H5AD 一同记录。平台不会把单细胞表达矩阵、细胞条形码、样本元数据或路径发送到
外部服务。

## 输入合同

在 H5AD 中写入：

```python
adata.obsm["X_aucell"] = auc_dataframe  # 行是 adata.obs_names，列为 regulon 名称
adata.uns["scenic_regulon_targets"] = {
    "TF1 (+)": ["TARGET1", "TARGET2"],
}
adata.uns["scenic"] = {
    "method": "pySCENIC AUCell",
    "database_version": "...",
}
```

`X_aucell` 也可改用其他 `obsm` 键；若其不是带列名的 DataFrame，必须在 `uns` 中提供与列数一一
对应的 regulon 名称列表。`scenic_regulon_targets` 可选：缺失时仍生成核心 AUCell/RSS/共活性图，
但不生成网络图。RSS 分组默认使用完成注释后的 `celltype`，因此该模块排在 `annotation` 之后。

## 输出与解释

| 输出 | 用途 | 不能说明什么 |
| --- | --- | --- |
| AUCell regulon 热图 | 展示高变异 regulon 在单细胞中的活性结构；细胞按分组排序，超大数据集确定性分层抽样展示 | TF mRNA 表达量，或条件差异显著性 |
| RSS 气泡图与条形图 | RSS = 1 − Jensen–Shannon distance，排名某 regulon 对细胞类型/cluster 的特异性 | 生物学样本重复层面的 p 值/FDR |
| regulon 靶基因网络 | 展示随输入导入的 TF–target membership，优先显示 RSS 高的 regulon | 本模块重新推断的因果/方向网络 |
| regulon 共活性相关热图 | 用 Spearman ρ 寻找共同活跃的 regulon 模块 | 单独证实 TF–target 有向边 |

所有选择后的图仍配套 CSV：完整 RSS、cell-group AUCell 汇总、热图展示选择、网络展示边和完整相关矩阵。完整逐细胞 AUCell 始终保留于输出 H5AD，而非默认导出包含 cell ID 的超大 CSV。

## 与既有转录因子模块的关系

| 模块 | 算法与输入 | 主要问题 |
| --- | --- | --- |
| `functional_state` | CollecTRI（平台冻结/项目内网络）+ `decoupler` ULM；按靶基因表达给少量指定 TF 打分 | 指定 TF 活性是否与功能通路及 sample × celltype 条件比较一致 |
| `scenic` | 已预计算的 SCENIC AUCell 与 regulon 定义 | 哪些 regulon 对细胞类型/状态最特异，哪些 regulon 共活跃 |
| `virtual_ko` | CellOracle Ridge GRN + in-silico perturbation | 某候选基因被扰动时的模拟状态转移和下游偏移 |

三者不能相互替代。尤其是 RSS 只用于细胞群特异性排序；涉及处理/疾病条件的结论仍需由
`sc_pseudobulk_deg` 与 `sc_cell_go` 在独立生物学样本层面验证。
