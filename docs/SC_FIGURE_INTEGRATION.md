# 单细胞转录组流程与 Bulk RNA 图形整合

## 1. 当前单细胞主流程

```text
数据导入
  → qc
  → normalize
  → hvg
  → dimred
  → batch_correct（可选）
  → clustering
  → subcluster（可选）
  → qc_reassess
  → annotation
  → sc_timecourse（时序项目可选）
  → deg / sc_cell_deg
  → sc_pseudobulk_deg
  → sc_cell_go（细胞级探索性 ORA 或样本级 pseudobulk ORA/GSEA）
  → trajectory / proportion / cell_communication
  → sc_csv_export
```

流程顺序和依赖以 `modules/__init__.py` 中的 `PIPELINE_ORDER`、
`PIPELINE_DEPS` 为准。条件比较优先使用 `sc_pseudobulk_deg`，因为统计重复必须是
独立生物学样本；`deg` 和 `sc_cell_deg` 的细胞级结果用于 marker 发现和探索，不能把
细胞数量当作生物学重复数。

## 2. 图形合同

**核心结论**：单细胞图形应依次证明数据质量可接受、细胞状态结构可复现、cluster
具有 marker 支持，并在真实样本重复层面验证条件差异及其通路解释。

**证据链**：

1. QC 前后分布：过滤是否改善深度、复杂度和线粒体负担。
2. PCA/UMAP/t-SNE：批次、cluster、注释和连续表达信号如何分布。
3. Marker Volcano/Heatmap：cluster 是否有一致且可解释的 marker 证据。
4. Pseudobulk Volcano/MA：条件差异是否在独立样本层面成立。
5. GO dotplot：差异基因是否支持可审查的生物过程解释。

**图形类型**：以 `quantitative grid` 为主，嵌入图是主要结构面板，QC、marker、
pseudobulk 和 GO 为从属证据。

**导出合同**：Figure Engine 图使用 89/183 mm 物理宽度、可编辑 SVG/PDF 文字、
600 DPI PNG，并输出独立 `nature_readiness.json`。完整数值始终保留在 H5AD/CSV，
图中截断标签、cluster 或 term 只影响展示。

**主要审查风险**：细胞级伪重复、UMAP 类别过多导致颜色重复、连续表达极值压缩色阶、
marker 热图标签过密、per-cluster 图数量失控、低 overlap GO term 被过度解释。

## 3. Bulk 图形迁移决策

| Bulk 图形/规范 | 单细胞落点 | 决策 | 原因 |
|---|---|---|---|
| Nature palette、字体、物理尺寸、SVG/PDF/PNG、readiness | 所有 Figure Engine 单细胞图 | 已整合 | 属于通用视觉与导出合同 |
| QC violin / diagnostics | `qc` 过滤前后指标分布 | 已整合 | 输入是 QC 指标分布，不改变统计含义 |
| PCA 风格 | `dimred` 的 cell-level PCA | 复用视觉语言，使用专用 embedding 模板 | Bulk 点代表样本，单细胞点代表细胞，不能直接复用样本 PCA 语义 |
| 新增 embedding 模板 | UMAP、t-SNE、PCA、基因表达、注释置信度 | 已整合 | 密集细胞层栅格化，文字/坐标轴保持矢量；连续和分类颜色分开处理 |
| PCA variance diagnostics | `dimred` | 已整合 | 方差解释率定义一致 |
| Volcano | `deg` cluster-versus-rest；`sc_pseudobulk_deg` 条件比较 | 已整合 | marker 图明确标记为探索性；pseudobulk 可直接作为样本级比较证据 |
| MA | `sc_pseudobulk_deg` | 已整合 | 需要样本级 mean expression，适合 pseudobulk，不用于普通 cluster marker |
| Heatmap | `deg` 的 cluster 平均 marker 表达 | 已整合 | 基因 × cluster 矩阵与模板输入一致，保留 gene-wise z-score 语义 |
| ORA/GO dotplot | `sc_cell_go` | 已整合 | 使用真实 Overlap、FDR 和 Genes；细胞级来源自动写入探索性警告 |
| 样本相关性热图 | 普通细胞矩阵 | 不迁移 | 细胞相关性会受状态连续性和细胞数支配；只适合 pseudobulk 样本矩阵 |
| Bulk normalization library-size panel | 单细胞 normalize | 暂不迁移 | 单细胞按细胞深度分布评估，不应画成样本文库柱状图 |
| UpSet | 多 cluster / 多比较 DEG | 后续候选 | 仅在用户明确选择比较集合时有意义，默认生成会迅速膨胀 |
| GSEA running curve | `sc_cell_go` 的 pseudobulk GSEA | 已整合 | 仅使用完整样本级 DEG 排名；细胞级分支不伪造 GSEA |
| WGCNA / correlation | 普通单细胞表达 | 暂不迁移 | 应先构造样本级 pseudobulk 或 metacell，再定义独立分析模块 |

## 4. 本轮实现

- 新增 `NatureEmbedding`，供 UMAP、t-SNE、cell-level PCA、基因表达和连续 QC
  指标使用；分类最多使用 20 色稳定色板，超出时写入 readiness 警告。
- `BaseAnalysis.build_publication_umap()` 保持原调用接口，但内部已经切到
  Figure Engine；因此 clustering、annotation、DEG 等现有调用点自动获得一致外观。
- `BaseAnalysis.save_matplotlib_figure()` 能识别 Figure Engine 合同，跳过会破坏物理
  尺寸的 legacy `tight bbox` 后处理，并自动注册 readiness JSON。
- `qc` 的过滤前后 violin 改用统一 diagnostics 模板，并对超大细胞集做确定性展示抽样。
- `dimred` 的 cell-level PCA、t-SNE、PCA variance 接入专用 embedding/diagnostics。
- `deg` 的首个 cluster Volcano 与 marker Heatmap 接入 Volcano/Heatmap 模板，并附加
  细胞级探索性说明。
- `sc_pseudobulk_deg` 按比较及 cluster 输出 Volcano + MA；默认每个已完成 cluster
  条件比较都有一对图。需要控制图量时，可将 `plot_max_clusters` 设为正整数，按显著 DEG
  数优先绘制前 N 个簇。无论该开关如何设置，每个可运行的 sample × cluster 单元都会
  强制输出按 condition 着色的样本 PCA、固定样本顺序的 Pearson 相关性热图、raw counts /
  cell 数 / DESeq2 size factor 三联图，以及可下载的样本 QC 表。每个完成的比较另输出
  Top DEG 的样本表达热图：blind PyDESeq2 VST（不能完成时明确标注 median-ratio log2
  回退）后按基因 z-score，列不聚类以便直接识别单 donor 驱动。
- `sc_cell_go` 保留细胞级 GO ORA 兼容分支，同时支持读取 `sc_pseudobulk_deg` 的
  样本级 DEG；pseudobulk ORA 使用每个分析单元的完整 tested-gene background，
  并支持上下调分开和预排序 GSEA。每个簇均独立计算；默认对每个有显著通路的簇出图，
  也可用 `plot_max_clusters` 限制图量。细胞级图仍明确要求用 pseudobulk 结果验证条件结论。

## 5. 建议的下一阶段

1. 为 annotation 建立“marker dotplot + confusion/置信度 + UMAP”组合图，而不是继续
   增加孤立单图。
2. 为 proportion 增加样本级组成图和适合重复测量/组成数据的统计模型，避免只按细胞总数
   做卡方检验。
3. 为 trajectory 建立“拓扑 + pseudotime + 样本覆盖 + 动态基因”组合合同，并检查每个
   条件/样本是否覆盖完整轨迹。
4. 在项目报告层按 QC → structure → annotation → sample-level comparison → pathway 的
   证据顺序自动组版，而不是按模块文件生成顺序堆图。
