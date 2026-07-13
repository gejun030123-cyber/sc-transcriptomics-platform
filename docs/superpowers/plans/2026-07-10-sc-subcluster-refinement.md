# 单细胞子簇精细分析与注释补强计划

## 目标

让用户从现有聚类结果中选择一个簇，独立进行重聚类，并在同一份可追溯结果中获得子簇差异表达、marker 热图和逐子簇通路富集；同时让子簇可进入现有细胞注释模块进行更细粒度的人工或自动注释。

## 设计

新增 `subcluster` 单细胞模块。它不修改原始全样本聚类，而是复制目标簇细胞为独立 AnnData：

1. 按 `source_cluster_key` 与 `target_cluster` 精确筛选细胞，保存 `parent_cluster` 来源信息。
2. 在 PCA（优先）或已有可用表示空间重建邻居图、Leiden/Louvain 和 UMAP，结果写入 `subcluster` 列。
3. 对 `subcluster` 执行 Wilcoxon marker 检验，输出每个子簇的 Top N 结果 CSV、marker 热图、子簇 UMAP 和细胞数图。
4. 可选地为每个子簇的显著上调 marker 调用 Enrichr（GO/KEGG/Reactome 等），并输出合并表和通路气泡图。网络服务或依赖不可用时保留前述分析结果并给出明确警告。
5. 将模块登记到单细胞模块列表、注册表、流水线顺序和依赖关系；所有可调参数写入共享 schema，保证前端、表单和 API/AI 参数过滤使用同一契约。

## 注释补强

现有注释模块已支持 marker、手工映射、CellTypist、置信度和证据图。本次通过让子簇输出标准 `subcluster` 列，使其可直接作为注释模块的 `cluster_key`。同时扩充 PBMC marker 集，细分 naive/memory T、细胞毒 T、NK、单核细胞与 DC 亚群，降低粗粒度 PBMC 注释混淆。

## 验收

- `subcluster` 能在 schema/UI 中选择，且仅依赖已完成的聚类结果。
- 对合成 AnnData 的模块运行产生子簇列、DEG CSV、热图和 h5ad；若 enrichment 开启且服务可用，还产生富集 CSV/图。
- 不存在目标簇、细胞不足或分组不足时给出可读错误，不静默生成误导结果。
- 运行 `py_compile`、schema/模块针对性 pytest，并以最小 AnnData smoke run 验证产物清单。
