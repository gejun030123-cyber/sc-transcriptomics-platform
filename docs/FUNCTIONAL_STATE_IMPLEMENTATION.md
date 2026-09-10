# 单细胞转录调控与功能状态模块（首期范围）

## 目标与位置

`functional_state` 位于 `annotation` 后，是一个可选分析模块：

```text
annotation
  → functional_state
  → sc_pseudobulk_deg → sc_cell_go
```

它回答“某细胞类型中哪些调控/代谢/炎症程序随条件变化”，不替代现有
`sc_pseudobulk_deg` 的原始 counts 差异表达或 `sc_cell_go` 的完整排序 GSEA。

## 三层证据合同

| 层级 | 首期算法 | 保存位置 | 不可替代的含义 |
|---|---|---|---|
| Gene expression | 选定的安全 log-normalized 矩阵中的单基因表达 | `gene_expression_statistics.csv` 与 dotplot；不复制到 `obs` | TF 或效应基因本身是否表达 |
| Pathway score | `scanpy.score_genes`；仅在数值条件不满足时记录 `mean_expression_fallback` | `adata.obs[fs_pathway_*]` | 功能基因程序是否协同变化 |
| TF activity | `decoupler 2.2.0` ULM，默认使用平台冻结的人类 CollecTRI 快照 | `adata.obs[fs_tf_*]` | 下游靶基因支持的调控活性 |

TF 基因的表达不会被命名为 TF activity。TF activity 是 ULM 的斜率 t 统计量：正值表示
与正向调控靶基因上调、负值表示相反方向；它不是 TF 自身表达量。`tf_expression_vs_activity.csv`
将两者并列，保留 TF coverage/status，不能用一列替代另一列。默认相关性 X 是真实的
`PPARA activity`，Y 是 `Acute myeloid chemokine inflammation score`；若 PPARA regulon 未通过 coverage，
结果会写为不可用，绝不会以 `PPARA-associated lipid-oxidation programme score` 冒充 TF activity。后者是
`score_genes` 的脂质氧化表达模块；它不宣称全部成员均为直接 PPARA 靶基因，且与网络支持的 ULM activity 刻意不同。
默认网络是管理员在受控数据卷下载并冻结的 `collectri_human.tsv`，来源为 OmniPath 的
CollecTRI（42,990 条关系；下载时记录 `decoupler` 版本、时间、参数和 SHA-256）。每次
分析会重新校验快照，并把网络校验值、版本和 ULM 方法写入 manifest。也可显式改用当前
项目内上传的 CollecTRI/DoRothEA-like 网络；该文件同样使用 ULM，且仍拒绝符号链接和项目外路径。

## 输入、统计与安全边界

- 必须存在 `sample_id`、`condition`、`celltype`（名称可配置）；一个样本只能属于一个条件。
  `batch`/`technical_batch` 默认不能充当样本或条件，除非研究者显式确认其生物学含义；首期仅支持 Human HGNC gene symbols。
- 用户须在 `auto`、`X`、`layers[counts]` 或指定的自定义 `layers[...]` 中选择表达来源。
  自定义层不能默认为 X。原始 counts、Pearson residual 或负值矩阵会从有效 counts 构建
  `layers['_functional_log1p']`，从不覆盖 `adata.X`；没有 counts 时停止。manifest 记录该选择、
  质量检查与处理决定。
- 当前 pathway score 仅提供 `scanpy_score_genes`；网页不展示未实现的 AUCell/UCell。所选方法
  和每个实际评分方法（含 fallback）均写入 manifest。
- UMAP、细胞级 violin 与逐细胞 correlation 仅为描述性展示。逐细胞相关图显示 Spearman
  ρ、`negligible/weak/moderate/strong/very strong` 效应量等级及细胞数，刻意不显示 P；细胞
  不是独立生物学重复。
- 相关性可输出全体细胞、各 celltype 内细胞、`sample × celltype`、各 celltype 内样本级四层。
  当前对“各 celltype 内细胞”和“sample × celltype”同时保留 pooled 图及按所选 reference/
  comparison condition 分层的补充图。前者检验跨条件总体趋势，后者检查该趋势是否仅由 condition
  位移驱动；分层的逐细胞图仍纯属描述性，`sample × celltype` 分层图仍复用同一 donor 的多个
  celltype 单元，均不显示 P/FDR。各 celltype 内样本级的双变量相关为未调整的探索性结果，只在
  该层记录 exploratory BH-FDR；它不能替代调整 condition/celltype/donor 的模型。
- 正式条件统计先聚合为 `sample × celltype`，再使用 Welch 或 Mann–Whitney；默认每个
  条件至少 3 个样本，每个 sample × celltype 至少 20 个细胞，并对全部 feature×celltype
  检验做 BH-FDR。任一组少于 2 个样本不运行；每组 2 个为 `exploratory_n2`；3–4 为
  `standard_n3_4`；至少 5 为 `preferred_n5_plus`。细胞数不会填补生物学重复。
- TF coverage 是固定合同：少于 5 个检测靶基因不计算；5–9 为 `low_coverage`（热图标 `*`）；
  10–19 为 `acceptable`；至少 20 为 `good`。coverage CSV 同时给出网络靶基因数、检测靶基因数、
  比例、状态和固定门槛。
- 运行前的设计预检显示当前范围下的细胞数、每条件样本数、每 celltype 细胞数及满足门槛的
  `sample × celltype` 单元数。
- 平台冻结的 CollecTRI 不接受网页指定路径，只能由管理员维护；项目内 TF network/GMT
  只接受当前项目目录下的普通文件，拒绝符号链接与目录穿越。输入 H5AD、GMT、网络文件
  的 SHA-256 与完整评分定义写入 manifest。

## 首期内置内容与输出

内置 `Lipid metabolism & inflammation` 基因表达面板按 Regulatory TF genes、Carnitine
shuttle、Mitochondrial beta-oxidation、Peroxisomal lipid metabolism (FAO)、Intracellular epithelial FA transport、Ketogenesis、
Acute myeloid chemokine inflammation 分节记录。manifest 会记录完整章节、实际命中和缺失基因；自定义基因列表会
被单独标为 custom，不混入内置面板。

内置通路包括 PPARA-associated lipid-oxidation programme、Carnitine shuttle (mitochondrial FA entry)、
线粒体 beta-oxidation、Peroxisomal lipid metabolism (FAO)、Intracellular epithelial FA transport、Ketogenesis、
Acute myeloid chemokine inflammation、TNF-NFkB 与 Type I and II IFN response。
急性髓系趋化因子 signature 强调 cytokine/chemokine output，不代表所有炎症状态；`TNF-NFkB response`
强调 canonical NF-kB feedback/signalling targets；二者允许有限的下游基因重叠，但绝不共享完全相同的定义。小型内置或自定义 signature
仍可评分，但在 `gene_set_coverage.csv` 标示为小型/探索性证据；应以版本固定的 Hallmark、Reactome
或 GO GMT 复核。`pathway_gene_set_overlap_qc.csv` 同时给出请求基因和当前表达矩阵实际检测基因的
pairwise Jaccard overlap；任意两个名称不同但请求或实际检测基因完全相同的集会在评分前被阻止，检测
基因 Jaccard ≥0.80 会标为 review。内置面板、内联自定义集、项目内 GMT 的版本/内容 SHA-256 都会
写入 manifest。

`IBD 类器官上皮主线`预设是一个可选的研究问题导向入口。它从已冻结的 Hallmark 快照选择
inflammatory response、TNF/NF-κB、IFNα/γ、hypoxia、ROS、UPR、apoptosis、OXPHOS、fatty-acid、
cholesterol、bile-acid、WNT、E2F 与 G2/M，并额外提供 mature absorptive/enterocyte differentiation、
Stem/TA regenerative state 和 cell-cycle proliferation 三个透明的 hypothesis signature。前者按标准通路
coverage 合同评分；后者在输出中单列为 `focused_epithelial_hypothesis_signature`，不能当作通用本体或临床状态。

输出包括：sample×celltype 评分、基因/TF/通路统计、TF expression versus activity、TF 与 gene-set coverage、相关性与
重叠 QC、关键基因 dotplot、功能状态热图、配对/条件 UMAP、描述性 violin、全体及各
celltype 内的相关散点，以及 Regulator–pathway concordance matrix。相关散点和该矩阵均额外
提供 ctrl/dis 分层版本：`regulator_pathway_concordance_by_condition.csv` 与并列矩阵用于检查
condition-specific 的效应量是否与 pooled 方向一致；每个 condition 的 donor 数仍会显示，且不输出
P/FDR 或显著性星号。该矩阵默认检查 PPARA–FAO、RELA–TNF-NFkB、STAT1–IFN 等方向一致性，使用
sample×celltype 的 Spearman 效应量展示。热图把 celltype 放在分组标题、condition 放在短子列，
避免长标签相互覆盖。violin 上的 `ns/*/**/***` 严格对应 sample×celltype BH-FDR，而不是
细胞级 p 值；用户可在不改变完整表或统计的前提下指定热图行及其顺序。逐细胞评分默认仅
保存于受控 H5AD；CSV 导出需用户显式开启。violin 的 n=2 显示为探索性，不会被伪装成
显著性星号。

## 明确延期

- 在分析运行时自动联网更新 CollecTRI/DoRothEA（仅允许管理员离线下载、复核并替换冻结快照）；
- AUCell/rank-based scoring 与 mixed-effects model（确认性关联建议为
  `InflammationScore ~ TFActivity + Celltype + Condition + (1|Sample)`）；
- 自动启动或复制 pseudobulk DEG/GSEA。

这些功能需要经版本审查的本地资源或明确的重复测量设计，出现实际需求后再增加。
