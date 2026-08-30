# 细胞注释模块修订设计

## 设计取舍

本版本以类器官为主要对象，保留“自动预测 → Marker 证据 → 冲突识别 → 人工复核”的主线。
不把 CellTypist、SingleR、Azimuth 或 scANVI 作为默认最终注释器：当前数据体系缺少稳定的同协议参考时，强行映射会把新状态压到最相近的已知标签。现在仅将已下载的人类 CellTypist 模型作为可选的独立参考证据，不覆盖 Marker 最终标签。

## 输出层级

细胞类型和细胞状态分开保存：

- `cell_lineage` / `cell_type_l1`：主要谱系
- `cell_type_l2`：主要细胞类型
- `cell_type_l3`：亚型或组织内细分类型
- `dominant_cell_state`：仅用于浏览的单一主状态，如 Cycling、Stress response
- `cell_state`：兼容旧下游的主状态别名，不代表状态互斥
- `cell_state_flags` / `state_high_*`：可并存的多个状态及独立布尔标记
- `developmental_state`：类器官表达成熟度，如 `progenitor_like`、`transitional`、`mature_like`
- `annotation_source_cluster`：产生本次标签的原始 cluster；它始终逐簇保留，
  即使多个 cluster 获得相同的生物学 cell type

## 有证据时合并来源 cluster

不再根据 marker 分数或谱相关性自动合并 cluster：相近的打分轮廓可能只是共享
应激、细胞周期或代谢状态，不能等同于同一细胞类型。完成逐簇 marker、UMAP
连续性和样本分布复核后，可在 `cluster_merge_map` 中显式确认，例如
`Cycling enterocytes=8,12`。`annotation_cluster_review.csv` 始终每个 source
cluster 一行，`annotation_source_cluster` 也始终保留；因此合并标签不会丢失
来源簇。输出的 `annotation_metadata.cluster_merge_policy` 为
`manual_confirmed_only`。旧参数 `merge_similar_threshold` 会被记录并忽略。

## 初步注释尽量覆盖全部细胞

默认开启“尽量注释所有细胞”（annotate_all）：只有完全没有正向 Marker
证据（所有候选分数 ≤ 0）或显式设置最低分数门槛时才标 `Unknown`；逐细胞
一致率低但没有特异性 Marker 支持的簇，保留最佳候选标签并标记为 review
（决策原因 low_cell_agreement_label_kept_for_review），而不是直接抹成
`Unknown`。关闭该选项可恢复严格的旧规则。

### Unknown 的再回收（初始注释阶段完成）

- 零正向证据但逐细胞投票一致的簇：候选分数不高于背景且候选 Marker 几乎
  未检出时，只要逐细胞 Top1 投票一致率 ≥ 阈值且最高分与次高分有正差距，
  就保留该候选标签并标记为 review（决策原因 `unknown_recovered_by_cell_vote`），
  只有投票也混乱时才保留 `Unknown`。
- LLM 模式：LLM 判为 `Unknown` 的簇回退到 Marker 规则候选标签，写入
  `llm_unknown_fallback` 列并标记复核；`llm_annotation` 列保留模型的原始判断。
  LLM 提示词也要求：只要有任何 Marker 支持就优先给出最具体的类型，
  `Unknown` 仅保留给完全没有可支持候选的簇。
- 计数写入 summary 与 `annotation_metadata.unknown_recovery`，便于审计。

## 初步注释细化（fine_annotation）

`Universal` 面板默认开启两级“细分注释”：先判定大谱系（Epithelial、
T cells、Myeloid、B/Plasma、NK、Mast、Fibroblast、Pericyte/Smooth muscle、
Endothelial、Cycling 等），再在带细分组的谱系内部输出更细的亚型：

- T cells → CD4+ T / CD8+ T / Treg / γδ T
- Myeloid → Monocyte/Macrophage / cDC1 / cDC2 / pDC / Neutrophils
- Pericyte/Smooth muscle → Pericytes / Smooth muscle cells
- Endothelial → Endothelial / Lymphatic endothelial

亚型判定被限制在同一大谱系内，且大谱系标签始终保留在候选列表里：亚型
证据不足的簇回退到大谱系标签，不会强塞进最近的亚型，也不会变成
Unknown。这些亚型同样写入 cell_lineage / cell_type_l1..l3 /
cell_ontology_id 层级列。`Colorectal` / `Organoid` 面板本身已是两级或
亚型粒度，细分开关不会重复拆分。

## 结直肠与类器官注释

`Colorectal` 面板使用两级决策：先区分 Epithelial、T cells、Myeloid、
Fibroblast、Endothelial 和待复核的 Neural-like 谱系；仅对 broad epithelial
cluster 再区分 stem/crypt、TA/S、TA/G2M、goblet、absorptive、BEST4+、
enteroendocrine、Paneth/LYZ+、inflammatory 和 regenerative/stress 亚型。
逐簇复核表保存 top/second candidate、正向分数、负向惩罚、Marker 检出、
特异性差值和最终决策原因。低逐细胞一致率不会单独否定一个具有至少两个
cluster-specific 正向 Marker 支持的候选。

## 高分辨率多簇支持

聚类分辨率提高后簇数可能远超 50 个。注释模块把分组安全上限提高到 200
个簇，并且 LLM 辅助注释会在 cluster 数超过单批上限（默认 30、单批最大
40）时自动分批调用模型，确保每个簇都得到注释；多批结果会合并写入
llm_annotation 与 annotations_json，request_hash 在单批时为字符串、
多批时为每批哈希的列表。

## 类器官成熟度

不要求用户主观选择“早期/成熟”。系统优先读取 `culture_day`、`day`、`timepoint`、`hour` 等 `adata.obs` 元数据，并独立计算：

- 前体模块分数
- 成熟模块分数
- 增殖模块分数
- `organoid_maturity_index`

成熟度结果只作为证据和复核信息，不参与细胞类型打分。没有时间元数据时，系统明确提示仅使用表达模块。

## 证据与冲突

- 正向 Marker：当前组织 panel 和数据驱动 Marker
- 负向 Marker：降低不匹配候选分数，不单独删除标签
- Doublet：只继承 QC 实际运行 caller 的 `predicted_doublet` 与可用的 caller-specific 证据；两个 Marker 模块同时较高仅记为 `mixed_lineage_review`，不得解释成 Doublet。Scrublet 的 batch-specific 阈值仅适用于实际使用 Scrublet 的运行。
- 环境 RNA：只检查相对当前 broad lineage 的异源高普遍性 Marker；同谱系 EPCAM/KRT 等不构成污染证据，且该启发式不能替代 empty-droplet 去污染
- Cell state：从 counts 重建 library-size normalized log1p 模块表达；每个状态保存独立连续分数与布尔标记，允许重叠
- Cluster 一致性、Top1/Top2、Marker 覆盖度和人工 review 共同写入 h5ad 与复核表

## 可追溯性

每次运行记录：

- `annotation_version`
- `annotation_comment`
- `celltypist_label`、`celltypist_confidence`、`celltypist_status`、`celltypist_comparison`（启用参考时）
- 使用的 marker_set、类器官类型和 cluster 列
- Marker 覆盖度、逐簇决策证据、dominant/multi-label 状态统计、QC Doublet 来源和异源环境 RNA 统计

人工映射仍通过 `manual` 方法执行；新的输出不会覆盖旧 h5ad，后续版本可以使用不同的 `annotation_version` 区分。

## CellTypist 参考模式

界面中的“启用 CellTypist 参考交叉验证”只读取 `data/references/celltypist/` 下的本地 `.pkl`，不会在分析期间联网下载。`prob match` 是默认模式，可保留 `Unknown`、低置信度和多标签结果；CellTypist 与 Marker 的冲突写入 `celltypist_comparison`，并进入逐簇复核证据。缺少 `celltypist` 依赖或模型文件时，任务回退到 Marker 注释并记录 warning。

## 精细注释 v2：cluster marker 优先 + 证据等级（当前实现）

### 版本化面板

- 内部稳定面板 key 为 `Colorectal_refined`，版本号 `colorectal_refined_v2`
  （`REFINED_PANEL_VERSION`），别名 `Colorectal_refined_v2` 指向同一面板。
- 面板包含 17 个可独立判别的 programme（CYP3A5+、KRT20/MALRD1+、
  HSD17B2+、MTTP/RBP2+、CKB/FABP1+、NDRG1/ANKRD37+、TA/stem-like、
  S/M/G2M 周期程序、TFF1/REG4+、MUC2/FCGBP+、SPINK4/CA4+、
  DDIT3/ATF3+ ER-stress、MAML3/CHRM3+ review 等），每个保留 4–8 个 marker
  与人工确认的锚点表。
- 每次运行输出 `annotation_marker_panel_manifest.json`（面板版本、marker、
  锚点、阈值、判定规则）并写入 `adata.uns['marker_panel_manifest']`；旧 h5ad
  与旧面板文件从不改写，`annotation_version` 版本化副本机制保持不变。

### Auto 三路路由

`resolve_auto_marker_set` 不再只在 Universal / Colorectal 之间二选一：

- 缺少明确组织程序 → `Universal`
- 有肠上皮亚型程序但细分证据有限 → `Colorectal`
- 多个细分 programme 均有 marker 支持 → `Colorectal_refined`
  （panel_version = colorectal_refined_v2）

检测规则：每个 programme 至少 2 个 marker 在 ≥5% 细胞中检出，且支持程序数
≥ 3（`AUTO_PANEL_DETECTION`）。因此同类肠上皮数据首次运行即可进入精细注释；
不符合本地 refined 证据的数据仍保守使用普通面板，不会误套项目标签。

### 判定规则与证据等级

每个 cluster 独立排名差异 marker（`select_cluster_marker_genes`），再综合
programme 分数、cluster 特异 marker 数量/特异性、单强锚点 marker、cluster
内逐细胞一致率、负向 marker 冲突；UMAP/KNN 邻接只写入复核表，绝不作为自动
合并依据。

- 两个以上特异 marker 与 programme 重合 → 细分标签，`confirmed`
- 单个强锚点 marker（如 CYP3A5、HSD17B2、MALRD1）→ 细分标签，`anchor-supported`
- 有分数但证据不足 → 保留标签并标记 `provisional`（review）
- 证据不足 → `Unresolved epithelial programme: GeneA/GeneB (cluster X; review)`
  或 `Unresolved organoid programme: ...`，`unresolved-review`

证据等级写入 `annotation_evidence_tier` 列（confirmed / anchor-supported /
provisional / unresolved-review），逐簇决策写入 `annotation_decision_reason`、
`annotation_cluster_review.csv` 与 `annotation_metadata.evidence_tier_counts`。

### 永远保留来源 cluster

固定保留 `annotation_source_cluster`、`celltype_l1/l2/l3`（以及兼容别名列
`cell_type_l1/l2/l3`）、`cell_state_flags`、`annotation_evidence_tier`、
`annotation_decision_reason`。相同标签绝不自动合并来源 cluster；只有人工填写
`cluster_merge_map` 或 `manual_map_csv` 时才合并显示标签。

### 可下载、可复用的审阅文件

- `annotation_cluster_review.csv`：每个来源 cluster 一行，含候选标签、Top
  marker、支持/缺失 marker、冲突、证据等级、UMAP/KNN 邻接、决策原因与
  needs_review。
- `annotation_manual_map_template.csv`：预填自动推荐标签；`manual_label`
  默认 `KEEP`，只填写需要修改的 cluster。保存后在注释任务参数
  `manual_map_csv` 中填写路径即可应用（只更新被填写的 cluster，tier 变
  confirmed，来源 cluster 永不删除）。
- `annotation_marker_panel_manifest.json`：面板版本、marker、锚点、阈值和
  判定规则。

### 标签压缩警告

精细面板（Colorectal refined / Organoid）若明显压缩标签（解析标签数少于
来源 cluster 数的 60%）会在 summary 输出 `label_reduction_warning` 并写入
runtime_warnings，提示下载复核表逐簇检查；不会自动合并。

### 不同类器官样品的支持

`Organoid` 面板沿用同一套 cluster-marker-first 判定：没有手工锚点表的面板
使用 `derive_panel_anchor_markers` 自动派生“只属于单个 programme”的锚点，
因此 kidney/liver/lung/cerebral 等类器官样品同样获得 evidence tier、
Unresolved review 标签与审阅文件，而不依赖固定参考图谱。
