# Nature Figure Engine 规范

## 1. 目标与图形合同

本规范将平台现有的“Nature 风格后处理”升级为确定性的 Nature Figure
Engine。引擎不替代统计分析，也不允许 AI 在每次绘图时临时选择字体、字号、
颜色、线宽或画布尺寸。

Phase 1 图形合同：

```text
Core conclusion:
  同一份 Bulk RNA-seq 分析结果经过固定模板渲染后，应在 89 mm 最终尺寸下
  保持科学信息清晰、视觉语言一致、文字可编辑且能够自动接受版面质量检查。
Figure archetype: quantitative grid
Target: Nature Portfolio / Nature Communications 风格的论文单图
Backend: Python (Matplotlib)
Final sizes: single column 89 mm; double column 183 mm
Primary evidence: PCA、Volcano、Heatmap、GSEA/ORA
Validation evidence: SVG/PDF/PNG 导出、文字重叠检测、最终尺寸检查
Source data: 模块 CSV/H5AD 输出
Reviewer risks:
  组合实验组使用过多颜色；热图标签拥挤；低 overlap 富集结果被过度解读；
  屏幕像素尺寸被误当作投稿尺寸；不同模块的视觉含义不一致。
```

## 2. 当前绘图架构

当前主数据流为：

```text
Analysis module
  -> 模块内部 Matplotlib / OmicVerse / legacy Plotly 绘图
  -> modules.figure_style 全局外观后处理
  -> BaseAnalysis.save_matplotlib_figure 或 legacy Plotly 静态转换
  -> PNG + SVG
  -> ResultFile / 结果页 / Figure Studio
```

主要组件：

- `modules/figure_style.py`：统一色板、字体、spine 和 Plotly layout。
- `modules/native_figures.py`：通用 heatmap、correlation、bar、scatter、UMAP 等画布。
- `modules/base.py`：保存 Matplotlib 图和兼容 Plotly 静态导出。
- `modules/reporting/static_rendering.py`：将历史 Plotly JSON 转为 Matplotlib PNG/SVG。
- `modules/figure_studio.py`：对部分平台图进行数据驱动重绘。
- `modules/bulk_pca.py`、`bulk_deg.py`、`bulk_heatmap.py`、
  `bulk_enrichment.py`：各自决定图形语义和布局。

Plotly 已不再是前端交互主输出，但仍作为历史模块的内部兼容格式存在。Phase 1
不会删除该兼容层。

## 3. 当前 nature-skills 工作方式

仓库中的 Nature 能力当前表现为共享 palette、font、spine 和静态导出规则；
外部 `nature-figure` 工作流提供“先结论、再证据、最后样式与 QA”的方法约束，
但尚未成为平台运行时的 Figure Director。

当前缺口是：分析模块仍直接创建具体 Matplotlib 图，样式层只能在图已经生成后统一
字体和轴线，无法纠正错误的图型、标签密度、实验因素编码和最终物理尺寸。

## 4. 已确认的问题

### 4.1 架构问题

1. `nature` 主题本质上仍是 palette/font 预设，不是科学图型合同。
2. 画布以 900x650 等屏幕像素表达，保存时近似除以 100 转成英寸，不是 89/183 mm
   的真实投稿尺寸。
3. 通用保存层只接受 PNG/SVG；缺少独立的 PDF 投稿验证路径。
4. 不同模块在 renderer 内混合数据选择、统计阈值、布局和样式。
5. 全局样式保存异常可能被宽泛捕获，不能形成明确的 readiness 失败报告。
6. Figure Studio 只对 Volcano、Heatmap、Correlation 和 Enrichment 提供数据驱动
   重绘；PCA、QC、MA 仍是 style-only。

### 4.2 当前真实输出问题

- PCA：组合组被映射成大量颜色，图例占据主要画布；缺少 batch marker 语义。
- Volcano：主体已经较好，但图例和长引线制造留白，标签选择只偏重 FDR。
- Heatmap：基因位于列方向导致标签碰撞；样本/表型注释条作为独立图片输出，不能
  和热图共同聚类、缩放和导出。
- ORA：真实结果中存在多个 `Count=1` 且 FDR 完全相同的通路，简单条形图呈现为
  等长柱，视觉上放大了很弱的证据。
- GSEA：当前复用横向条形图，未形成以 NES、FDR 和 gene-set size 为核心的模板。
- QC：固定 3x2 网格会保留没有信息的空面板。

## 5. 目标运行时架构

```text
Analysis result
      -> NatureFigureDirector
      -> FigureSpec
      -> deterministic Plot Template
      -> NatureStyle
      -> FigureValidator
      -> FigureComposer
      -> SVG / PDF / PNG / TIFF
```

职责边界：

- Director：只判断证据角色、模板、模式和科学参数；不直接绘图。
- FigureSpec：序列化的科学与输出合同，是 renderer 的唯一配置入口。
- Template：固定图型语法；不得回退到默认 seaborn/matplotlib publication 图。
- Style：控制真实物理尺寸、字体、颜色、线宽、边距和字体嵌入。
- Validator：在最终尺寸渲染后检查并给出分数及可执行问题。
- Composer：用 Matplotlib SubFigure 直接组合可编辑 panel，不栅格化子图，也不重新
  解释单图数据。

## 6. 当前目录

```text
figure_engine/
├── __init__.py
├── director.py
├── composer.py
├── accessibility.py
├── export.py
├── spec.py
├── validator.py
├── visual_regression.py
├── style/
│   ├── __init__.py
│   ├── nature.py
│   ├── nature_communications.py
│   └── nature_aging.py
└── templates/
    ├── __init__.py
    ├── common.py
    ├── correlation.py
    ├── pca.py
    ├── volcano.py
    ├── heatmap.py
    ├── enrichment.py
    ├── enrichment_extended.py
    ├── gene_sets.py
    ├── ma.py
    ├── upset.py
    ├── wgcna.py
    └── diagnostics.py
```

Phase 1 模板：

- `NaturePCA`
- `NatureVolcano`
- `NatureHeatmap`
- `NatureGSEA`
- `NatureEnrichmentDotplot`（GO/KEGG ORA）

Secondary Bulk RNA diagnostics now use the same engine instead of legacy static
builders: normalisation library/distribution/PCA checks, QC overview/PCA/pairs/
violin/distance, PCA variance/loadings, DEG boxplots, and time-course Q-Q,
trajectory and heatmap views.  Invariant QC metrics are omitted rather than
rendered as empty panels; sample/gene labels are density-capped at final size.

通路富集扩展模板（统一 FigureSpec + NatureStyle + Validator）：

- `NatureEnrichmentBarplot`：ORA 的 `-log10(FDR)`，GSEA 的双向 NES。
- `NatureEnrichmentChord`：真实 Genes/geneID 成员的可读性截断 chord-style membership。
- `NatureEnrichmentCnetplot`：通路方形节点、基因圆形节点和成员边。
- `NatureEnrichmentEmapplot`：基于共享基因 Jaccard 的通路相似性图。
- `NatureGSEARunning`：真实 ranking + hit_indices/RES 的归一化 running ES + hit rug。

网络图默认有严格的 pathway/gene 数量上限，并在 readiness JSON 中写明截断；不能从
只有 FDR/NES 的汇总 CSV 伪造 GSEA running 曲线。

保留现有 `modules/native_figures.py` 作为 Standard/历史图工具；Nature Portfolio
模板不再通过全局后处理改变尺寸和字号。

## 7. FigureSpec 与模式

三种模式：

- `standard`：探索与网页快速检查。
- `publication`：普通论文静态图。
- `nature_portfolio`：严格物理尺寸、固定模板和 validator 门禁。

用户可控制科学参数：FDR、log2FC、Top N、基因/通路、聚类方法、注释字段、
是否显示置信椭圆。用户不能在 Nature Portfolio 模式中任意覆盖字体、字号、轴宽、
背景或核心 palette。

`FigureSpec.fc_threshold` 统一表示 `abs(log2FC)` 阈值；现有 Bulk DEG 页面中的
普通 Fold Change 在适配层转换为 `log2(FC)`，避免隐含单位歧义。

## 8. Validator 合同

每张图输出 0-100 readiness score，并至少检查：

- canvas 与 89/183 mm 合同
- 最小字号和字体一致性
- 颜色数量与红绿冲突
- legend、label、axis、colorbar 越界或重叠
- 白色背景
- SVG/PDF 是否存在及 SVG 文字是否可编辑
- PNG DPI
- 模板产生的语义告警，例如所有 ORA term 都只有一个重叠基因

`score < 90` 或存在 error 均为失败；90 以上且无 error 才能通过投稿门禁。严格模式
下三组以上文字重叠直接记为 error。分析统计仍可成功，图形 readiness 报告必须
单独保存，不能静默吞掉布局异常；正式 release 可显式调用 `FigureValidator.require`
阻断不合格产物。

## 9. 需要修改的现有文件

Phase 1 只小范围修改：

- `modules/bulk_pca.py`：主 PCA 适配到 `NaturePCA`。
- `modules/bulk_deg.py`：主 Volcano 适配到 `NatureVolcano`。
- `modules/bulk_heatmap.py`：publication 热图适配到 `NatureHeatmap`，将注释合并。
- `modules/bulk_enrichment.py`：ORA 默认 dotplot，GSEA 使用 NES dotplot。
- `modules/schemas.py`：仅增加科学内容参数，不开放破坏视觉合同的参数。

Phase 1 不重写：

- 统计计算和结果 CSV
- 旧 Plotly 转换层
- 非 Bulk 分析模块
- Figure Studio 全部能力
- Multi-panel Composer

## 10. 验证计划

1. 单元测试 FigureSpec、尺寸换算、每种模板和 validator。
2. 用仓库已有 36 样本 Bulk RNA 参考项目生成 before/after。
3. 验证 89 mm SVG 的物理宽度和 `<text>` 可编辑文字。
4. 验证 PDF 可打开、PNG 像素/DPI、无明显文字重叠和越界。
5. 使用真实 DEG、PCA、热图表达矩阵和 ORA 表；GSEA 缺少现成真实结果时，使用
   项目 GSEA 结果 schema 构造确定性的 schema fixture，并明确标记为模板验证数据。
6. 将 readiness JSON、输入来源摘要和导出文件放入独立验证目录，不覆盖原图。

## 11. 分阶段计划

### Phase 1：单图引擎

实现 FigureSpec、Director、NatureStyle、PCA、Volcano、Heatmap、GSEA/ORA dotplot、
Validator 和真实数据验证。

### Phase 2：组合与平台交付（已完成本轮范围）

实现 FigureComposer、panel labels、共享图例、PDF/TIFF 结果注册和前端三模式控件。
Figure Studio 对 PCA/MA/QC 的数据驱动编辑保留为独立后续工作。

### Phase 3：模板扩展（核心模板与 Bulk secondary views 已完成）

本轮已增加 Correlation、MA、GSVA/ssGSEA、UpSet、WGCNA，并将 Bulk 的 QC、标准化、
PCA 辅助图、DEG 箱线图/整合矩阵和 Time-course 视图迁移到固定 diagnostics 模板。
PPI、Forest 和 Deconvolution 留待后续。

### Phase 4：投稿门禁（基础能力已完成）

本轮已加入视觉快照回归、三类色觉模拟和 90 分自动门禁；项目级视觉词汇表、
Source Data bundle、图注统计合同和目标期刊规则版本化留待后续。

## 12. Phase 1 实施与验证结果

Phase 1 已完成，并以独立 `figure_engine` 包接入 Bulk PCA、DEG、Heatmap 和
Enrichment 的主 publication 输出。原统计计算、CSV、备用 PC panel、MA plot、样本
相关性图和历史 Plotly 兼容层均未重写。

仓库参考项目 `699134b0-e1c` 的 36 个真实 Bulk RNA 样本用于最终尺寸验证；PCA、
Volcano、Heatmap 和 ORA 均直接读取现有 H5AD/CSV。仓库没有已完成的 GSEA 结果，
因此 GSEA 只使用真实 ORA term/FDR 加确定性 NES/setSize 的 schema fixture，不能
冒充真实 GSEA 生物学结论。

| 模板 | 89 mm readiness | 数据来源 | 结论 |
|---|---:|---|---|
| NaturePCA | 100/100 | 真实 PCA H5AD | pass |
| NatureVolcano | 100/100 | 真实 DEG CSV | pass |
| NatureHeatmap | 100/100 | 真实表达矩阵 + DEG 选基因 | pass |
| NatureEnrichmentDotplot | 92/100 | 真实 GO BP ORA | pass，保留 Count=1 证据警告 |
| NatureGSEA | 100/100 | schema fixture | renderer pass，非生物学验证 |

每个模板均生成 SVG、PDF、600 DPI PNG 和 readiness JSON；导出不使用
`bbox_inches='tight'`，因此不会改变 89 mm 物理画布。对比总览和逐图报告位于
`artifacts/nature_figure_phase1/`。

## 13. Phase 2 实施与验证结果

本轮新增 `NatureFigureComposer`，统一小写 panel label 和共享图例，并以 SubFigure
直接组合 PCA、MA、Correlation 和 ORA，保留 SVG/PDF 的文字与几何对象。平台前端
现在提供 Standard、Publication、Nature Portfolio 三种模式；严格模式锁定字体、
字号、颜色、线宽和背景，只开放科学参数、89/183 mm 投稿宽度与导出格式。

主保存链路与结果页现已注册 SVG、PDF、PNG、TIFF；TIFF 使用 LZW 压缩并校验 DPI，
项目图包也会包含 PDF/TIFF。新增固定模板：

- `NatureMA`
- `NatureCorrelationHeatmap`
- `NatureGSVA` / `NatureSSGSEA`
- `NatureUpSet`
- `NatureWGCNA`

本次通路富集扩展在真实项目 `bec48c50-bc9` 的 GO-BP ORA 表上完成单栏验证，输出位于
`artifacts/nature_enrichment_views/`，包括 dotplot、barplot、chord、cnetplot、emapplot
和 running-score 的 PNG/SVG/PDF/TIFF 及 QA JSON；五种 ORA 图均为 89 mm 宽、文字重叠
为 0、达到 92/100 ready，running-score 达到 100/100。该分数扣分只来自“结果被截断”
和“按命中基因 Jaccard 压缩冗余”的语义告警，不是布局失败。running-score 使用项目
真实 DEG 排名和 ORA 成员字段验证数据契约；该项目没有已完成的 GSEA 任务，因此 NES/FDR
明确显示为 n/a，不将验证曲线当作 GSEA 生物学结论。

富集结果闭环规则：ORA 结果表和绘图共同使用 `GeneRatio = Count / InputGeneCount`，
其中 `InputGeneCount` 会写入 CSV；默认 Top N 为 12，长 term 使用动态行间距，图中去掉
GO/KEGG ID 但结果表保留原始 term。命中基因 Jaccard 达到阈值（默认 0.85）的重复通路
只在图中压缩，并在 QA 中记录数量。当前输入合同没有 GO ontology 图，因此尚未启用
基于父子层级的额外压缩；需要在结果接口提供 ontology DAG 后再增加该规则。

Validator ready 阈值提升为 90，并加入色觉可分辨性、TIFF DPI、组合 panel 元数据与
显式 release gate。视觉回归使用确定性的结构指标和 256-bit dHash；色觉检查输出
protanopia、deuteranopia、tritanopia 三种预览。

真实数据验证使用项目 `699134b0-e1c`：MA、Correlation、UpSet 和组合图读取真实
DEG/H5AD/ORA 结果；GSVA/ssGSEA 和 WGCNA 因项目无已完成结果表，使用真实表达矩阵
派生的确定性 schema fixture，并明确禁止将其作为生物学结论。

| 输出 | Readiness | SVG/PDF/PNG/TIFF | 说明 |
|---|---:|:---:|---|
| NatureMA | 100/100 | yes | 真实 DEG |
| NatureCorrelationHeatmap | 100/100 | yes | 真实表达矩阵 |
| NatureGSVA | 100/100 | yes | 真实表达派生 fixture |
| NatureSSGSEA | 100/100 | yes | 真实表达派生 fixture |
| NatureUpSet | 100/100 | yes | 真实多比较 DEG |
| NatureWGCNA | 100/100 | yes | 真实表达派生 fixture |
| 2×2 Composer | 100/100 | yes | 真实 PCA/DEG/Correlation/ORA |

完整产物、before/after、三种色觉模拟和逐图 QA 位于
`artifacts/nature_figure_phase2/`。

仍未完成：Figure Studio 的 PCA/MA/QC 数据驱动编辑、PPI/Time-course/Forest 等
扩展模板、项目级 Source Data bundle 和按目标期刊版本化的规则包。

## 14. 富集图目标通路与基因标签选择

富集图不再把自动 Top N 作为唯一视图。Bulk enrichment 参数现在提供：

- `pathway_selection`: `top`、`selected` 或 `selected_plus_top`；
- `target_pathways`: 按通路名称、GO/KEGG ID 或名称片段指定，支持逗号、分号或换行；
- `gene_label_strategy`: 网络图中的 `all`、`shared`、`selected`、`none`；
- `target_genes` 与 `max_gene_labels`: 指定或限制直接显示的基因符号。

这套参数统一作用于 dotplot、barplot、GSEA、GSEA running、Chord、Cnetplot 和
Emapplot。数据库 ID 只在结果表中保留，图面使用自动换行的名称；基因节点/连线与
基因文字标签分离，因此达到标签上限时不会静默删除成员关系，而会在 readiness 报告中
明确提示。目标通路超过某种图型的单栏容量时同样会报告被截断的数量，并建议减少目标
通路或切换双栏。

真实项目 `699134b0-e1c` 已用非 Top 的 GO-BP term 验证目标通路选择，输出位于
`artifacts/enrichment_target_selection_audit/`；单栏 Chord/Cnet/Barplot 均通过
Validator，真实 ORA 的 Count=1 证据问题仍按原规则保留为语义告警。

## 15. Human Bulk enrichment P1/P2 统计与 Source Data 合同

P1 在完整 P0 统计表之上增加解释层，但不改变“完整检验表是主证据”的规则：

- ORA 默认仅检验与真实背景相交后含 10–500 个基因的通路；该过滤独立于查询命中，
  零命中但大小合格的通路仍进入同一 BH family；
- 每条 ORA 结果导出 Haldane–Anscombe 校正 odds ratio、近似 95% CI、fold enrichment、
  GeneRatio 和 BackgroundRatio；
- 显著通路按命中基因（GSEA 使用 leading-edge 基因）Jaccard 生成冗余簇、代表通路和
  相似度；该标注不删除任何统计行，图形可仅展示代表项；
- GSEA 除 pathway-level NES/FDR 外，另存每个 leading-edge 基因的排名位置和排名分数；
- GO BP/MF/CC 使用仓库固定的 Human 2023 GMT 快照，所有库记录文件名、年份、修改时间、
  pathway/gene 数和 SHA256。

P2 增加自定义与交付层：`Custom_GMT` 接受 Human GMT 文本，限制 5 MB，运行时保存原文
快照并计算 SHA256。每次 ORA/GSEA 运行均生成独立 Source Data 工作簿，包含完整结果、
显著结果、映射 QC、未映射基因、参数、冗余簇，以及 GSEA leading-edge（适用时）；同时
生成纯文本统计方法记录和机器可读 QC JSON。项目级整合表继续按 comparison、database、
method 和 direction 分层，不跨数据库混合 FDR family。

真实 Human DEG + GO-BP 2023 冒烟验证：12,987 个提交背景基因中 9,935 个进入有效注释
背景；10–500 大小过滤后检验 2,859 条通路，其中 251 条为零命中；346 条达到 FDR<0.05，
形成 326 个非破坏性冗余簇。输出效应量字段完整，Source Data 六个 ORA sheet 和方法记录
均可读取。真实排名的 50 条 GO 通路 GSEA 冒烟产生 931 条可回溯 leading-edge 基因记录。
