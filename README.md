# 单细胞与 Bulk RNA-seq AI 分析平台

这是一个基于 Flask、Scanpy、Matplotlib 和 OmicVerse 风格规范的 Web 端转录组与基因组分析平台。平台同时覆盖单细胞转录组、Bulk RNA-seq 和 WES 外显子组，支持传统按模块执行的分析流程，也支持用户通过 AI 对话触发分析、检查结果、调整参数，并围绕特定目标创建候选分支进行参数搜索。

项目当前定位不是单纯的流程封装，而是一个可交互的分析工作台：

- 用户可以上传数据，在网页中按模块执行 scRNA-seq 或 Bulk RNA-seq 分析。
- AI 助手可以读取项目状态、理解已完成步骤、调用分析工具并提出参数调整方案。
- 对“分群不满意”“必须找到最接近某种细胞类型的分群”等需求，平台提供 marker 评分、候选分支、参数 sweep 和人工采纳机制。
- 所有分析结果会落到项目目录下，包含 h5ad 中间文件、科研级静态图与必要的统计审计；单细胞 DEG 与富集的完整表只作为受控内部输入，不在界面导出。

## 核心能力

### Web 分析工作台

- 项目管理：创建项目、上传数据、查看项目状态和历史任务；主线任务完成后项目状态会自动同步为 `completed`，存在运行中任务时为 `processing`，全部失败时为 `failed`。
- AI 主工作台：项目内 `/projects/<pid>/workspace` 页面聚合单细胞、Bulk RNA、WES 分析入口、最近任务和数据文件，并预留 Bulk ATAC-seq 工作流。
- 异步任务：分析任务通过后台 worker 执行，前端可查看进度、日志和失败信息。
- 流程模板：可将已排好顺序的单细胞或 Bulk 模块及参数保存为项目模板；在分析页选择项目内输入文件后可一键提交。服务端会重新读取模板，并校验模块依赖、参数 schema、项目路径与符号链接，浏览器不能借模板注入任意模块、参数或命令。
- 参数面板：每个模块有结构化参数 schema，包含中文标签、默认值、类型和帮助说明。
- 结果管理：任务结果写入数据库，支持图表查看、CSV/XLSX 表格下载和 h5ad 中间文件下载；模块写出的每个结果文件都会登记到任务结果清单。
- 科研级图表：Scanpy/Matplotlib 原生图默认保存为 300 dpi PNG 和 SVG 矢量图，结果页优先展示 PNG，并提供 PNG/SVG 下载。
- 静态图表：统一使用 Matplotlib/OmicVerse 风格，网页默认展示 300 dpi PNG，并提供 SVG 矢量图下载。
- 图形工作台：从分析结果或用户上传图片创建非破坏性的样式版本，支持预览、PNG/SVG 导出和版本追踪。
- 可复现交付：任务 manifest、项目报告、离线图表图库和 pipeline 报告会记录参数、summary、产物和复核证据。
- 安全路径校验：API 读取和 AI 工具调用会限制在项目目录内，拒绝路径穿越和符号链接输入。

### 单细胞转录组分析

平台提供 16 个核心单细胞分析模块，并补充批量导入、细胞级 DEG、GO 富集、pseudobulk DEG 和标准 CSV 结果包等交付模块。核心顺序为：

```text
qc -> normalize -> hvg -> dimred -> batch_correct -> clustering
  -> qc_reassess / annotation -> functional_state / subcluster / sc_timecourse / deg
  -> trajectory / proportion / neighborhood_da / cell_communication / virtual_ko
```

| 模块 | 主要功能 | 典型输出 |
| --- | --- | --- |
| `qc` | 线粒体、核糖体、血红蛋白比例标记；双细胞检测（默认 scDblFinder，要求 `pyscdblfinder>=0.2.0`，缺失时停止而不回退；Scrublet 以显式参数独立运行并保留 simulated score 证据）；细胞周期评分；复杂度指标；批次自适应 QC | QC violin、counts vs genes scatter、novelty plot、cell-cycle plot、过滤前后 QC 对比、实际 caller 的 doublet score 直方图；Scrublet 额外输出 observed vs simulated doublet 分布图 |
| `normalize` | `log1p` 或 Pearson residuals 标准化，保留 counts layer | 标准化后 h5ad、library size 图、表达值分布图 |
| `hvg` | 高变异基因选择，支持批次感知、force include、排除 MT/CC 基因、细胞周期评分和回归 | HVG scatter、HVG rank plot、HVG 标记 |
| `dimred` | PCA、UMAP，可选 t-SNE/MDE，支持自动 PC 选择 | PCA variance、PCA scatter、UMAP QC 着色图 |
| `batch_correct` | Harmony、ComBat、BBKNN、Scanorama、SysVI、scVI 等批次整合入口；支持 CPU/GPU 深度模型路径 | 校正后 embedding/graph、整合前后 UMAP、batch ASW、cluster batch entropy、最大批次占比、邻居混合、图连通性、指标表 |
| `clustering` | 多分辨率 Leiden/Louvain 聚类，支持主分辨率、自动分辨率评分 | 各分辨率 UMAP、多分辨率 UMAP、cluster 标签 UMAP、cluster 细胞数图、cluster 批次组成图、分辨率 Sankey |
| `subcluster` | 对指定 cluster 进行子簇重聚类、差异表达和通路富集 | 子簇 UMAP、marker 表、热图、富集结果 |
| `qc_reassess` | 聚类后按簇评估 doublet、MT、ribo、细胞数，支持自动移除低质量簇 | 低质量簇表、按簇 QC 汇总图、doublet/MT UMAP、QC 指标 UMAP 面板、低质量簇高亮图 |
| `annotation` | 分层 cell lineage/type/subtype、独立 multi-label cell state；`Colorectal` 面板先判 broad lineage 再细分 goblet/TA/absorptive/inflammatory 等上皮亚型；逐簇保存候选、正负 Marker、决策原因；Doublet 继承 QC 实际运行的 caller，Marker 混合不冒充 Doublet；环境 RNA 仅检查异源谱系 Marker；可选本地 CellTypist 参考（不覆盖 Marker 标签） | 细胞类型 UMAP、细胞类型组成图、marker score heatmap、marker dotplot、marker 表达验证图、annotation score UMAP、成熟度 UMAP、CellTypist 参考 UMAP、逐簇复核表 |
| `functional_state` | 在选定 celltype 范围内并列评估关键基因表达、`scanpy.score_genes` 功能通路和 CollecTRI + decoupler ULM 支持的 TF 活性；正式条件比较以 sample × celltype 聚合，逐细胞图仅作描述性展示 | 基因表达 dotplot、TF/通路评分热图、condition UMAP、样本级 violin/FDR、相关性散点、TF–pathway concordance matrix、coverage/overlap QC 与 manifest |
| `sc_timecourse` | 按真实时间点进行样本级细胞组成和伪 bulk 基因动态分析，区分描述性趋势与统计推断 | 时间点 UMAP、组成曲线/热图、动态基因表和趋势图 |
| `deg` | cluster/celltype marker（探索性） | 火山图、显著 DEG 数量图、top marker UMAP 面板、dotplot、marker heatmap、基因表达 UMAP |
| `trajectory` | Diffusion Map、DPT、PAGA 拟时序 | pseudotime UMAP、pseudotime 分布图、PAGA 图、基因随拟时序变化图 |
| `proportion` | 细胞比例统计和组间比较，支持卡方、Fisher、置换检验 | 堆叠柱图、比例 heatmap、饼图、比例统计表 |
| `neighborhood_da` | 在高维整合/PCA 空间构建重叠 KNN 邻域，并以每个样本的邻域比例进行条件比较；适合连续的 stress/TA/absorptive 上皮状态 | 邻域元数据、样本级邻域比例、Mann–Whitney + 每比较 BH-FDR、DA UMAP 图 |
| `cell_communication` | 基于 LIANA 的配体-受体通讯分析 | 通讯热图、气泡图、通讯网络图、交互表 |
| `virtual_ko` | 基于 CellOracle 的 GRN 推断（Ridge 回归）与 in silico 基因敲除扰动模拟；内置人类 promoter base GRN（hg19/hg38）或上传自定义 base GRN；在独立 celloracle 环境（Python 3.9/3.10）中运行 | GRN 边表 CSV、每个基因的状态偏移 CSV、Top 受调控基因 CSV、quiver 向量场、模拟流场网格、细胞分群+流场、偏移分布图（PNG+SVG）、含模拟结果的 h5ad |
| `scenic` | 本地读取预计算的 SCENIC AUCell 及 regulon 定义；不以 TF 表达替代 regulon 活性，也不在网页运行中隐式重建网络 | AUCell regulon 热图、RSS 特异性气泡/条形图、可选 regulon 靶基因网络、regulon 共活性相关热图、完整 RSS/选择审计 CSV |

`annotation` 还提供可选的脱敏 LLM 辅助注释：仅向配置的模型发送 cluster 级 marker 摘要（不含表达矩阵、细胞条码、样本元数据或项目路径），返回结果只作为候选证据，须人工复核后采纳。

`functional_state` 的 TF activity 绝不等同于 TF 基因表达：前者仅在本地、版本冻结的网络靶基因覆盖充分时计算。它还将 ctrl/dis 分层相关性与 regulator–pathway 一致性作为探索性效应量展示，不以细胞数替代生物学重复。完整的输入安全、资源校验、统计合同和延期范围见 [功能状态模块说明](docs/FUNCTIONAL_STATE_IMPLEMENTATION.md)；标准 Hallmark/Reactome/GO 通路库的管理员同步、许可和覆盖度策略见[托管基因集注册表](docs/MANAGED_GENE_SET_REGISTRY.md)。

`scenic` 与上述 ULM TF activity 是互补模块：它展示预计算 SCENIC regulon 的 AUCell 与 RSS 细胞群特异性，而不是在网页中从表达量临时推断 GRN。输入结构、图形解释边界以及与 CellOracle/样本级 pseudobulk 的分工见 [SCENIC 分析说明](docs/SCENIC_ANALYSIS.md)。

### Bulk RNA-seq 分析

平台提供 8 个 Bulk RNA-seq 模块，典型顺序为：

```text
bulk_qc -> bulk_normalize -> bulk_pca -> bulk_deg -> bulk_heatmap -> bulk_enrichment
```

多时间点或多比较实验可追加：

```text
bulk_timecourse
bulk_deg_integration
```

| 模块 | 主要功能 | 典型输出 |
| --- | --- | --- |
| `bulk_qc` | 文库大小、检测基因数、MT/ribo 比例、Gini 复杂度、离群样本检测、自动分组 | QC 表、样本距离图、相关性热图、PCA 异常检测 |
| `bulk_normalize` | DESeq2 size factor、TMM、CPM、VST、rlog、log2 quantile | 标准化矩阵、文库大小对比 |
| `bulk_pca` | PCA、UMAP、t-SNE；按分组列着色；载荷分析 | PCA/UMAP 图、解释方差图、样本聚类 |
| `bulk_deg` | t-test、Mann-Whitney、DESeq2、edgeR、limma；pairwise 或 LRT；多比较 | DEG CSV、火山图、MA 图、基因箱线图 |
| `bulk_heatmap` | Top 变异基因、DEG、表达式筛选或手动基因热图；支持 z-score/center/winsorize | 表达热图、样本相关性热图 |
| `bulk_enrichment` | ORA/GSEA，支持 GO、KEGG、WikiPathways、Reactome 等数据库 | 富集表、barplot、dotplot、GSEA 曲线 |
| `bulk_timecourse` | 多时间点差异检测、spline F-test、轨迹聚类 | 时间趋势图、cluster profile |
| `bulk_deg_integration` | 多组 DEG 整合，支持 UpSet、Venn、一致性评分、logFC 矩阵和表达式筛选 | UpSet/Venn、logFC heatmap、共同/特异基因表 |

### WES 外显子组分析

平台在单细胞与 Bulk 表达分析之外，补充科研级 WES（全外显子组测序）能力，面向人类短读长双端测序的 germline SNV/InDel 与肿瘤 somatic SNV/InDel。WES 不进入线性 `PIPELINE_ORDER`，而是接入独立的 `WORKFLOW_REGISTRY`，由 `modules/workflows/` 提供统一的输入登记、预检、运行追踪和产物收集：

| Workflow | 说明 |
| --- | --- |
| `wes_germline` | GATK HaplotypeCaller gVCF 胚系 SNV/InDel（单样本/小家系） |
| `wes_somatic` | Mutect2 + 官方过滤链；支持 tumor-normal 与受限 tumor-only 模式 |
| `wes_vcf_normalize` | 已调用 VCF 的标准化与 VEP offline 注释（germline/somatic 模式） |

- 三类合法输入入口：FASTQ（含 SRA 本地 `fasterq-dump` 转换）、已处理 BAM/CRAM、已调用 VCF。
- 执行后端为固定版本 nf-core/sarek + Nextflow local executor；平台负责启动、查询、取消、`-resume` 和产物收集，Web 端不提供任意命令执行。
- 参考资源必须登记并通过 checksum 校验；`WES_REQUIRE_VALIDATED_REFERENCES` 默认开启，未验证的 reference bundle 无法进入生产运行。上传的 capture BED 一律登记为 `test_only`，需管理员审核后才能用于生产分析。
- 服务器 WES 数据通过 `WES_SOURCE_ROOTS` 白名单只读接入，拒绝符号链接；FASTQ、BAM、CRAM 和全量 VCF 不会发送给外部 AI。
- 项目内提供 WES 面板（`/projects/<pid>/wes`）：样本与 manifest 登记、capture kit 管理、运行状态、SRA 转换任务和产物安全下载。
- 首期边界：不做 CNV、SV、MSI、TMB、突变特征、ACMG 自动分级和大队列 joint genotyping；胚系与 somatic 使用独立 workflow contract，变异结果不混表、不共用过滤结论。详见 [WES 分析方案](docs/WES_ANALYSIS_PLAN.md) 与 [P2/P3 Runbook](docs/WES_P2_P3_RUNBOOK.md)。

### AI 对话与目标驱动 Agent

AI 助手通过 `/api/chat` 接入，支持 Anthropic 风格 API 和 OpenAI 兼容 API。工具调用分为只读工具和需要用户确认的写操作。

只读工具会自动执行：

- `get_project_status`：读取项目上传文件、已完成任务、中间文件和分支任务数量。
- `get_task_results`：读取某个任务的结果摘要和结果文件列表。
- `list_modules`：列出单细胞、Bulk 或全部模块。
- `inspect_analysis_state`：检查当前 h5ad、已完成步骤、聚类键、embedding、细胞数等。
- `inspect_adata`：查看 AnnData 结构，包括 obs/var 列、obsm、聚类列和数据维度。
- `get_cluster_summary`：按 cluster 统计细胞数、UMAP 中心、batch 分布和 QC 均值。
- `score_cell_type_signature`：用 marker 对每个 cluster 评分，寻找最接近目标细胞类型的分群。
- `list_builtin_markers`：列出内置细胞类型 marker。

写操作需要用户确认：

- `run_analysis`：运行一个分析模块。
- `propose_parameter_sweep`：为目标细胞类型生成候选参数组合。
- `run_parameter_sweep`：创建候选分支并执行参数搜索。
- `start_goal_agent`：启动目标驱动分析会话。
- `continue_goal_agent`：根据用户反馈继续调整或收敛。

适合 AI Agent 的典型需求：

```text
查看这个项目现在分析到哪一步了。
现在的分群不满意，帮我找到最接近 microglia 的 cluster。
用 PBMC marker 检查当前注释是否合理。
resolution 0.8 的单核细胞群太混，帮我设计几个候选参数。
采用评分最高的候选，后续分析从它继续。
```

#### 目标驱动分群优化流程

当用户提出“必须找到某类细胞”或“当前分群不满意”时，平台按下面的方式工作：

1. AI 读取当前项目状态和最新 h5ad。
2. AI 检查可用的 cluster key、UMAP、QC 指标和已有注释。
3. AI 使用内置或用户提供的 marker 对各 cluster 评分。
4. 如果当前分群不足以满足目标，AI 生成参数 sweep 候选，例如不同 `resolution`、`n_neighbors`、`distance_metric`。
5. 用户确认后，平台在 `data/projects/<pid>/branches/` 下创建候选分支并运行。
6. 每个候选分支被自动评分，结果显示在 Agent 实验面板。
7. 用户可以采纳某个候选；采纳后项目的当前分析基线指向该分支输出，但原主线结果不会被删除。

## 数据输入

常见输入格式：

- `.h5ad`：推荐的单细胞中间数据格式。
- 10x `.h5` / `.hdf5`：上传后可在上传页导入为标准 `.h5ad`。
- `.loom` / `.zarr`：上传或放入项目 `uploads/` 后可导入为标准 `.h5ad`。
- `.csv`、`.txt`、`.tsv`：表达矩阵或样本表。
- `.xlsx`、`.xls`：Bulk 表达矩阵。
- `.mtx.gz` / 10x 文件组合：上传 `barcodes.tsv(.gz)`、`features.tsv(.gz)` 或 `genes.tsv(.gz)`、`matrix.mtx(.gz)` 后可转换为 `.h5ad`。
- 多批次 10x ZIP：上传两个或更多分别包含 `filtered_feature_bc_matrix` 的 ZIP，上传页会为每个 ZIP 提供「样本 ID / 批次名称 / 条件」标记（条件如 对照组、疾病组，可自定义）；平台解压后写入 `obs["sample_id"]`、`obs["condition"]`、`obs["batch"]`，再合并为一个标准 `.h5ad`。未填写条件时 condition 留空，后续需补标后才能运行样本级推断。
- 服务器目录批量 10x：管理员配置 `SC_BATCH_SOURCE_ROOTS` 后，可从网页递归扫描大量 10x 目录、下载并填写样本 manifest（`sample_id,matrix_dir,condition,replicate,batch`），再合并为标准 `.h5ad`。目录名不会被自动当作生物学分组。

单细胞 DEG 按统计单位分开：`deg` 只用于 cluster marker；`sc_cell_deg` 用于同一样本不同簇、不同样本同一簇或条件间的细胞级探索性比较；具有独立生物学重复的正式条件结论使用 `sc_pseudobulk_deg`。细胞级 log2FC 始终由原始 `counts` 重建 log1p 表达，Pearson residual 不直接用于 fold change。富集页面必须选择一个已完成的 DEG 任务，系统读取其受控内部结果并核对 AnnData 版本，不按目录中“最新 CSV”猜测来源。默认交付图和 JSON 统计审计；完整 DEG/富集表只供内部下游使用，不在界面导出。详细设计与延期项见 [人类单细胞 DEG 与富集规划](docs/sc-deg-enrichment-plan.md)。

使用同一个“结果包文件夹”（默认 `sc_batch_results`）时，所有输出会整理为：

```text
results/sc_batch_results/
  01_group_proportions/      # 样本设计、细胞元数据、样本×Leiden 簇比例
  02_gene_expression/        # pseudobulk raw counts、log2(CPM+1)
  03_differential_expression/# 每个 comparison 的 DEG、比较注册表
  04_go_enrichment/          # GO 结果；未运行/不可运行时含状态表
```

项目目录结构遵循：

```text
data/projects/<project_id>/
  uploads/       # 用户上传原始文件
  intermediate/  # 各模块输出 h5ad
  plots/         # PNG/SVG 科研图
  results/       # CSV、报告、下载结果
  branches/      # AI/Agent 候选分支
  presets/       # 项目级参数预设
```

## 安装与启动

### 从 GitHub 克隆后的可运行范围

可以从 GitHub 克隆后直接启动 Web 应用、创建项目，并运行不依赖外部参考的常规
scRNA-seq / Bulk RNA-seq 流程；`data/`、SQLite 数据库、缓存和项目目录都会在首次
运行时创建。仓库**不会**也不应包含研究数据、密钥、CellTypist 模型、通路库、WES
参考或容器镜像。因此，“启动应用”和“启用所有可选分析能力”是两件不同的事：

| 能力 | 克隆 + `requirements.txt` | 额外需要 |
| --- | --- | --- |
| Web、项目管理、常规单细胞/Bulk 分析、PNG/SVG、Excel 输入/输出 | 可以 | 用户自行上传的表达数据；不附带示例人类数据 |
| Bulk `bulk_enrichment` | 首次运行可联网下载 Enrichr 基因集 | 无外网时，管理员须预置 `genesets/<library>.txt` 或 `data/go_gene_sets/<library>.gmt` |
| `sc_cell_go` 的离线 ORA/GSEA | 可以（内置 GO BP/CC/MF、Reactome 与 WikiPathways Human 快照） | KEGG 未随库分发；需由具备相应授权的管理员提供本地 GMT/TXT。自定义库可在项目内上传或设置 `SC_CELL_GO_GENE_SET_DIR` |
| `functional_state` 标准 Hallmark/Reactome/GO 与 CollecTRI TF activity | 可以（内置受版本与 SHA-256 约束的公开快照） | 可设置 `FUNCTIONAL_STATE_RESOURCE_DIR` 以使用管理员审核的替代/更新资源 |
| CellTypist 参考注释 | 包已安装，但模型未随仓库提供 | 可信 `.pkl` 模型放入 `SC_CELLTYPIST_MODEL_DIR` |
| `scenic` | 可运行 | 输入 H5AD 必须已有 AUCell `obsm` 矩阵和 regulon 定义；网页不重建 SCENIC 网络 |
| `virtual_ko` | 不可仅靠主 Python 环境运行 | 独立 Python 3.9/3.10 的 CellOracle 环境，以及 base GRN（首次使用内置 GRN 需要联网下载） |
| WES 真正执行 | 默认关闭 | Nextflow/Docker 与 Sarek 工具链、已验证的 GRCh38/capture/VEP 资源；见 [WES P2/P3 Runbook](docs/WES_P2_P3_RUNBOOK.md) |

这使克隆环境保持轻量且不分发受许可或敏感资源；缺少上述可选资源时，对应模块应给出
`unavailable`、警告或资源缺失信息，而不是把任务标为成功。

### 基础安装

建议在 Linux 上使用 Python 3.10–3.12 和独立虚拟环境。安装时需要访问 PyPI；部分
科学计算扩展（例如 `louvain`、`inmoose`）可能从源码构建，因此无预编译 wheel 的主机还
需要标准 C/C++ 构建工具和 Python 开发头文件。`requirements.txt` 是运行时依赖清单，尚未
锁定版本，生产环境应在验证后生成自己的锁定文件。

```bash
git clone https://github.com/gejun030123-cyber/sc-transcriptomics-platform.git
cd sc-transcriptomics-platform

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

# 可选：只有在需要覆盖默认路径、访问控制或 AI/WES 配置时才创建 .env
cp .env.example .env

# 初始化检查：会创建空的本地 data/ 与 instance/，不需要参考数据
python - <<'PY'
from app import create_app
create_app()
print('应用初始化成功')
PY

python app.py
```

默认访问地址为 `http://localhost:5000`。`SECRET_KEY` 未设置时会在
`instance/.secret_key` 自动生成；服务器若暴露到受控内网以外，必须在 `.env` 设置
强 `PLATFORM_ACCESS_PASSWORD`，并通过 HTTPS / VPN 等受控入口访问。

`requirements.txt` 已包括以下 Python 包；安装成功不等于其所需的本地参考也已具备：

- `openpyxl`：`.xlsx` 输入、样本 manifest 与 Excel 结果导出。
- `celltypist`：本地人类 CellTypist 参考交叉验证；模型文件需另行放置，且不会覆盖 Marker 注释。
- `liana`、`gseapy`、`decoupler==2.2.0`：分别用于细胞通讯、富集和 CollecTRI ULM TF activity。
- `harmonypy`、`bbknn`、`scanorama`、`scvi-tools`、`torch`、`kneed`：批次整合、深度模型与 Kneedle 自动 PC 选择。
- `pydeseq2`、`inmoose`、`pysam`：Bulk 统计和 WES 内容预检；WES calling 本身仍是外部工作流。
- `openai`、`anthropic`：AI 客户端；未配置 API key 时 AI 对话会明确拒绝，不影响其他模块。

下列路径有兼容回退，因而不放入基础依赖；需要相应增强结果时再安装：

```bash
# MDE 降维、模糊 c-means 时间轨迹
python -m pip install pymde fuzzy-c-means
```

### 本地资源与管理员准备

所有资源目录必须由管理员控制，使用普通文件而非符号链接；不要把 FASTQ、BAM、CRAM、
全量 VCF、项目路径或身份信息发送给外部服务。

- **CellTypist：** 将已审核的模型（例如 `Immune_All_Low.pkl`）放入
  `data/references/celltypist/`，或在 `.env` 设置 `SC_CELLTYPIST_MODEL_DIR`。
- **`sc_cell_go` 离线富集：** 克隆自带 GO BP/CC/MF、Reactome 与 WikiPathways Human 的只读快照，
  足以完成默认 ORA/GSEA。pseudobulk ORA/GSEA 强制使用本地库；KEGG 不在仓库内，必须由已获得
  对应授权的管理员放置同名 `<library>.gmt` / `<library>.txt`，自定义库也可通过该方式或
  `SC_CELL_GO_GENE_SET_DIR` 提供。
  只有兼容的细胞级 ORA 可显式选择 Enrichr，此时基因列表会发送到该在线服务。

#### KEGG Human 基因集补充（仅限已获授权的管理员）

KEGG 数据不随本仓库分发，也不能提交回公共 GitHub。请由持有适当 KEGG 许可的管理员自行从
[KEGG FTP / subscription](https://www.kegg.jp/kegg/download/) 或机构批准的内部参考资源取得数据：
学术 FTP 仅向订阅者开放，非学术使用须另行取得许可。不要把 KEGG REST API 当作公开下载渠道；
该 API 仅供学术用户学术使用，且有速率限制。

将管理员根据其许可导出的 **Human gene-symbol** GMT/TXT 保存为（UTF-8、Enrichr-style：
`term<TAB>source<TAB>GENE...`）：

```text
data/go_gene_sets/KEGG_2021_Human.gmt
```

也可在部署环境设置 `SC_CELL_GO_GENE_SET_DIR=/srv/sc-platform/references/go_gene_sets`。随后在
`sc_cell_go` 选择 **KEGG Human（需授权本地文件）** 即可。该目录中的 KEGG 文件会优先使用；未
放入的 GO、Reactome 和 WikiPathways 仍自动使用随仓库快照，因此不必复制默认资源。请在运行
manifest 中保留 KEGG release、获取日期、许可证/订阅依据和 SHA-256，且除非获得 KEGG 的明确再分发许可，
不得把该文件推送到 GitHub。
- **`functional_state`：** 克隆自带 Hallmark、Reactome、GO 与 CollecTRI v2.0 快照；标准默认
  路径是 `resources/functional_state_resources/`。管理员如需更新，设置
  `FUNCTIONAL_STATE_RESOURCE_DIR` 指向独立受控目录，再同步到该目录的 `gene_sets/`：

  ```bash
  python scripts/sync_managed_gene_sets.py \
    --resource-dir "$FUNCTIONAL_STATE_RESOURCE_DIR"
  ```

  该管理员命令从固定的官方来源下载 Hallmark、Reactome 和 GO，并记录版本、许可与
  SHA-256；细节见[托管基因集注册表](docs/MANAGED_GENE_SET_REGISTRY.md)。随库 CollecTRI
  快照同样记录来源、版本与 SHA-256。所有第三方资源、署名与许可见
  [`resources/THIRD_PARTY_DATA_NOTICES.md`](resources/THIRD_PARTY_DATA_NOTICES.md)；更新后应先
  在隔离环境验证再替换，不允许网页运行时下载。
- **`virtual_ko`：** 设置 `CELLORACLE_PYTHON` 指向独立 Python 3.9/3.10 环境；该环境
  必须安装 CellOracle 及其依赖。内置人类 promoter base GRN 首次使用会下载并缓存；无外网
  部署应预先缓存，或上传经审核的 base GRN。
- **WES：** 保持 `WES_EXECUTOR_ENABLED=false`，直到 Nextflow、Docker/Apptainer profile、
  Sarek release、GRCh38 FASTA/索引/known sites、somatic PoN、VEP cache 与 capture BED
  都已登记和 checksum 验证。它们不能通过 `pip` 或本仓库获得；完整清单与登记步骤见
  [WES P2/P3 Runbook](docs/WES_P2_P3_RUNBOOK.md)。

启动后可在网页的依赖状态页查看模块可用性；设置了 `AI_API_TOKEN` 时，访问
`/api/system/dependencies` 需要携带相应 Bearer token。每个模块会分别报告 `missing`
（无法运行）与 `optional_missing`（只影响增强能力）。

频繁开发时，可在服务器 `.env` 中临时开启自动重载：

```dotenv
PLATFORM_AUTO_RELOAD=1
```

先执行一次 `./scripts/platform-service.sh restart` 使开关生效。此后修改 Python、
HTML、CSS 或 JavaScript 文件会自动重载应用进程，浏览器刷新即可看到效果；Flask
交互式调试器仍保持关闭。安装/升级依赖或修改 `.env` 后仍应执行完整重启。
自动重载会中断进程内正在执行的分析任务，因此只应在无人运行分析的开发时段
开启；跑正式分析前将该值改为 `0` 并重启平台。

### 图像输出与下载

单细胞的 UMAP、聚类 UMAP、注释 UMAP，以及聚类/注释/差异表达 DotPlot 会同时生成两类结果：

- `PNG`：300 dpi 位图，网页默认展示，适合汇报和快速预览。
- `SVG`：矢量图，适合论文排版和后续编辑。
- 所有图表均为 Matplotlib 静态输出：PNG 用于网页预览，SVG 用于论文排版和后续编辑；不再加载 Plotly.js 或提供交互式图表。

任务结果文件类型统一由 `modules.base.VALID_RESULT_FILE_TYPES` 定义，目前包括 `csv`、`xlsx`、`png`、`svg`、`json`、`txt`、`h5ad` 和内部兼容类型 `plotly_json`。旧的 Plotly JSON 结果在登记时会尽力转换为 PNG/SVG，前端不再将其作为交互式结果展示。

新的静态图只会在任务重新运行时生成；历史任务需要重新执行对应模块才能获得 PNG/SVG。若在分析参数面板关闭某种导出格式，平台会按所选格式保存结果。

### 图形工作台

项目详情页中的“图形工作台”（`/projects/<pid>/figure-studio`）支持：

- 从已生成的火山图、热图、相关性图、富集图等结果中选择可编辑来源。
- 调整标题、字体、颜色、尺寸和版式并实时预览。
- 上传外部 PNG/JPG/SVG 图片进行样式包装；原始分析结果不会被覆盖。
- 将每次调整保存为独立版本，并下载 PNG 或 SVG。

## AI API 配置

通过环境变量或项目根目录的 `.env` 配置 AI。`.env` 已被 `.gitignore` 忽略，适合本地保存 API key：

```bash
export AI_API_KEY="your-api-key"
export AI_API_URL="https://your-openai-compatible-or-anthropic-endpoint"
export AI_MODEL="your-model-name"
export AI_API_TOKEN="optional-local-api-token"
python app.py
```

DeepSeek Anthropic-compatible 示例：

```bash
AI_API_URL="https://api.deepseek.com/anthropic"
AI_API_KEY="your-deepseek-key"
AI_MODEL="deepseek-chat"
AI_API_TOKEN=""
```

说明：

- 启动后可直接打开 `/settings/ai`（首页和 AI 主工作台均有入口）更换协议、API 地址、API Key 和模型；保存后立即生效，无需重启。API Key 只在服务端保存，页面仅显示脱敏值。
- “测试连接”会使用当前表单参数发起一次轻量服务检查，不会保存配置；确认成功后再点击“保存并应用”。“恢复环境变量默认值”会删除平台页面保存的覆盖项。
- `AI_API_KEY` 为空时，`/api/chat` 会返回未配置错误。
- `AI_API_URL` 中包含 `anthropic` 或 `claude` 时走 Anthropic Messages 格式，否则走 OpenAI compatible Chat Completions 格式。
- Anthropic-compatible 调用使用内置 HTTP 客户端，不强依赖本地安装 `anthropic` SDK。
- `AI_API_TOKEN` 为空时跳过本地 API token 认证；设置后需要请求头 `Authorization: Bearer <token>`。
- 写操作工具不会直接执行，会先返回 `proposed_tools`，前端确认后再调用 `/api/chat/approve`。
- AI 查询任务结果时，若存在 PNG/JPG/SVG 等图片，`/api/chat` 会返回 `attachments`；AI 主工作台和项目详情页会直接显示缩略图，可点击查看原图或下载。

## 使用流程

### 单细胞网页流程

1. 创建项目。
2. 上传 `.h5ad`、10x `.h5`、10x 三文件、`.loom` 或 `.zarr` 数据，并在上传页导入为标准 `.h5ad`。大量服务器 10x 样本使用“服务器目录批量 10x 导入”：先扫描、填写 manifest，再合并导入。
3. 执行 `qc`、`normalize`、`hvg`、`dimred`。
4. 执行 `clustering` 并查看多分辨率 UMAP、Sankey 和带标签 cluster UMAP。
5. 执行 `qc_reassess` 检查低质量簇。
6. 执行 `annotation`：结直肠/肠类器官选择 `Colorectal`，未知组织先用 `Universal`，其他场景选择 `Organoid`、`PBMC`、`Immune`、`Blood`、`TME` 或自定义 marker。
7. 执行 `deg`，检查 marker heatmap、火山图和统计审计。
8. 根据项目需要继续 `trajectory`、`proportion`、`neighborhood_da`、`cell_communication`。
9. 单样本/探索性项目运行 `sc_cell_deg`，有生物学重复的正式条件比较运行 `sc_pseudobulk_deg`；随后在 `sc_cell_go` 明确选择对应的 DEG 任务运行 Human ORA/GSEA。默认无需整理 CSV 数据包。

### Bulk 网页流程

1. 创建项目并上传表达矩阵。
2. 执行 `bulk_qc`，确认样本质量和自动分组。
3. 执行 `bulk_normalize`。
4. 执行 `bulk_pca` 检查分组、批次和异常样本。
5. 执行 `bulk_deg`，设置比较组、统计方法和阈值。
6. 执行 `bulk_heatmap` 或 `bulk_enrichment`。
7. 多比较场景执行 `bulk_deg_integration`。

### 流程模板一键运行

在单细胞或 Bulk 分析页中，将当前拖拽排序后的模块和参数保存为“流程模板”。之后在对应分析首页选择该模板与项目内输入文件，点击“运行流程”即可后台顺序执行。模板只保存分析配置，不保存输入路径；每次运行都会重新检查模块顺序、参数可用性和输入文件边界。长时间或重计算步骤仍由现有后台并发上限控制，失败流程可从结果页按已有机制续跑。

### WES 网页流程

1. 管理员在服务器 `.env` 中配置 `WES_SOURCE_ROOTS` 及参考资源路径，需要真正启动运行时开启 `WES_EXECUTOR_ENABLED`。
2. 打开项目详情页的 WES 面板，按入口类型登记样本 manifest（FASTQ / BAM / CRAM / VCF）。
3. 选择 workflow（胚系 / tumor-normal / 受限 tumor-only）和已验证的参考 bundle、capture kit，完成运行前预检。
4. 启动 Nextflow 运行；平台跟踪状态、PID、日志和退出码，支持取消与 `-resume`。
5. 完成后在面板查看和下载 filtered VCF、注释表、QC 与报告；所有产物经安全路径校验后提供下载。

### AI 对话流程

1. 打开项目详情页右下角 AI 对话面板。
2. 输入需求，例如“查看项目状态”或“帮我找最接近 B cell 的 cluster”。
3. AI 自动调用只读工具检查项目。
4. 如果 AI 建议运行分析或参数搜索，前端会显示待确认工具调用。
5. 用户确认后任务提交到后台执行。
6. 在 Agent 实验面板查看候选分支、评分证据和采纳按钮。

## REST API 摘要

常用接口：

| 接口 | 方法 | 说明 |
| --- | --- | --- |
| `/api/chat` | `POST` | AI 对话入口 |
| `/settings/ai` | `GET` | AI API 设置页面 |
| `/api/settings/ai` | `GET/POST` | 查看或保存 AI API 配置（Key 脱敏返回） |
| `/api/settings/ai/test` | `POST` | 测试未保存的 AI API 参数 |
| `/api/settings/ai/reset` | `POST` | 恢复环境变量默认 AI 配置 |
| `/api/chat/approve` | `POST` | 执行用户确认后的 AI 工具 |
| `/api/chat/history/<pid>` | `GET` | 获取项目聊天历史 |
| `/api/tasks/<task_id>/status` | `GET` | 查看分析任务状态 |
| `/api/projects/<pid>/tasks` | `GET` | 列出项目任务 |
| `/api/projects/<pid>/adata-info` | `GET` | 获取当前 AnnData 信息 |
| `/api/result-file/<file_id>` | `GET` | 下载 PNG、SVG、CSV 或其他结果文件 |
| `/projects/<pid>/figure-studio` | `GET` | 打开非破坏性图形工作台 |
| `/projects/<pid>/workspace` | `GET` | AI 主工作台页面 |
| `/projects/<pid>/wes` | `GET` | 项目 WES 面板（workflow、运行、参考资源、SRA 任务） |
| `/projects/<pid>/upload/import-10x-batches` | `POST` | 接收两个或更多 `batch_zip`（可带 `sample_id`/`condition` 表单），合并为带 `sample_id`/`condition`/`batch` 列的 h5ad |
| `/projects/<pid>/upload/discover-10x-directory` | `POST` | 在允许的服务器目录内扫描 10x 矩阵并生成 manifest 模板 |
| `/projects/<pid>/upload/import-10x-manifest` | `POST` | 校验已填写的 sample manifest 后异步合并任意多个 10x 样本 |
| `/api/projects/<pid>/pipeline-runs` | `POST/GET` | 创建或列出批量 pipeline run |
| `/api/projects/<pid>/pipeline-templates/<preset_id>/run` | `POST` | 按服务端保存的项目/全局流程模板启动一次 pipeline run（请求体只提供项目内输入文件） |
| `/api/pipeline-runs/<run_id>/status` | `GET` | 查看 pipeline run 状态 |
| `/api/projects/<pid>/current-context` | `GET` | 查看当前分析基线 |
| `/api/projects/<pid>/branches` | `GET/POST` | 列出或创建候选分支 |
| `/api/projects/<pid>/branches/<branch_id>/run` | `POST` | 在候选分支上运行模块 |
| `/api/projects/<pid>/branches/<branch_id>/score` | `POST` | 对候选分支评分 |
| `/api/projects/<pid>/branches/<branch_id>/accept` | `POST` | 采纳候选分支 |
| `/api/projects/<pid>/agent/jobs` | `GET` | 查看 Agent 异步任务 |

## 脚本

### PBMC3k 单细胞参考流程

```bash
python scripts/run_pbmc3k_sc_reference.py
```

该脚本用于下载或读取 PBMC3k 参考数据，运行经典单细胞分析流程，并生成结果摘要和图表 gallery。它适合检查单细胞流程是否能产生类似 Scanpy PBMC 教程的结果。

### Bulk 参考流程

```bash
python scripts/run_bulk_reference.py --help
```

该脚本用于把指定 Bulk 表达矩阵跑过核心流程，生成静态科研图、报告和结果表。常用于回归测试或离线交付分析结果。

## 测试

安装测试依赖后运行全部测试：

```bash
python -m pip install -r requirements-dev.txt
python -m pytest tests/ -q
```

常用分组：

```bash
# AI/Agent、分支、marker scoring 和 schema 相关测试
python -m pytest \
  tests/test_agent_jobs.py \
  tests/test_agent_models.py \
  tests/test_agent_orchestrator.py \
  tests/test_agent_tools.py \
  tests/test_branch_routes.py \
  tests/test_signature_scoring.py \
  tests/test_schemas.py \
  -q

# 单细胞/Bulk 模块和可视化基础测试
python -m pytest tests/test_p2_modules.py tests/test_p3_modules.py tests/test_visualization.py -q

# WES workflow、AI 设置和工作台相关测试
python -m pytest tests/test_wes_workflows.py tests/test_ai_settings.py tests/test_workspace.py -q

# 语义一致性检查
python -m pytest tests/test_semantic.py tests/test_semantic_full.py -q
```

说明：部分测试依赖可选生信包、示例数据或本地环境。新功能开发建议至少运行与改动相关的 targeted tests，再根据影响范围追加 integration/semantic tests。

## 项目结构

```text
app.py                      # Flask 应用入口和 blueprint 注册
config.py                   # 路径、AI、并发、资源配置
database.py                 # SQLite schema 初始化和迁移
models.py                   # Project、Task、ResultFile、Agent、Branch 等模型
worker.py                   # 异步分析任务执行器

routes/
  main.py                   # 首页
  projects.py               # 项目 CRUD
  upload.py                 # 文件上传和 10x 检查
  analysis.py               # 分析配置页
  results.py                # 结果展示和下载
  api.py                    # 常规 REST API
  chat.py                   # AI 对话 API
  ai_settings.py            # AI API 配置页面与设置接口
  branches.py               # Agent session、candidate branch、score、accept API
  figure_studio.py          # 图形工作台页面、预览、上传和版本下载
  workspace.py              # AI 主工作台页面
  wes.py                    # 项目 WES 面板与安全产物下载
  auth.py                   # Bearer token 鉴权装饰器

modules/
  __init__.py               # MODULE_REGISTRY、PIPELINE_ORDER、PIPELINE_DEPS
  schemas.py                # 参数 schema 和 UI 元数据
  base.py                   # 分析基类和结果保存工具
  io_utils.py               # 数据读取和 gene name remap
  visualization.py          # 统一静态图风格和颜色工具
  inspect_utils.py          # AnnData 检查工具
  expression_parser.py      # Bulk 多比较表达式解析
  ai_adapter.py             # Anthropic/OpenAI compatible AI 适配器
  ai_config.py              # 环境变量与平台页面配置的运行时合并
  ai_tools.py               # AI 工具执行后端
  agent_orchestrator.py     # 目标驱动 Agent 编排
  agent_jobs.py             # 参数 sweep 异步 job
  design_preflight.py       # 分析前实验设计和比较组检查
  figure_style.py           # 统一 Nature 风格和 Matplotlib 参数
  native_figures.py         # 原生静态科研图构建器
  figure_studio.py          # 图形来源解析、预览和版本保存
  cell_markers.py           # 内置 marker 定义
  llm_annotation.py         # 脱敏 cluster 级 LLM 辅助注释适配器
  virtual_ko.py             # CellOracle GRN 与 in silico 敲除模块
  celloracle_worker.py      # 独立 celloracle 环境的子进程 worker
  evaluators/               # 候选分支评分器
  workflows/                # WES 外部工作流：Sarek/Nextflow、预检、运行与产物
  reporting/                # manifest、项目报告、pipeline 和复核证据
  qc.py ... deg.py          # 单细胞分析模块
  bulk_qc.py ...            # Bulk 分析模块
  convert_10x.py            # 10x 转换模块

figure_engine/              # 确定性的 Nature Figure Engine
  director.py               # 科学参数 -> FigureSpec -> 固定模板
  composer.py               # SubFigure 多 panel、panel label、共享图例
  validator.py              # 90 分投稿门禁、重叠/裁切/导出 QA
  accessibility.py          # 色觉模拟与色板可分辨性检查
  style/                    # Nature/Nature Communications/Nature Aging
  templates/                # PCA、Volcano、Heatmap、GSEA/ORA、MA、Correlation、
                            # GSVA/ssGSEA、UpSet、WGCNA

templates/                  # Jinja2 页面模板
tests/                      # 单元、集成、语义和 Agent 测试
scripts/                    # 参考流程、工具脚本和平台服务管理
docs/                       # 设计、计划和审批文档
deploy/                     # 平台的 systemd 服务文件
genesets/                   # 通路基因集资源
data/                       # 本地项目数据，默认不纳入 git
```

## 仓库与本地数据边界

Git 仓库只保存源代码、模板、测试和维护文档。以下内容只属于本地运行环境，默认由 `.gitignore` 排除：

- `data/`、`instance/`、`cache/`：用户数据、SQLite 数据库、中间矩阵和运行缓存。
- `data/references/`、`data/go_gene_sets/`、`data/functional_state_resources/`、`genesets/`：
  管理员覆盖资源、参考模型与运行缓存；不随仓库提供。克隆可直接使用 `resources/` 下带许可、只读的
  基线 GO/Reactome/Hallmark/CollecTRI 快照。
- `bulk_reference_output*/`、`artifacts/`：参考流程、图形验证和导出产物，可由脚本重新生成。
- `.env`：本地密钥和部署参数，禁止提交。
- `.runtime/`、Python/pytest 缓存、覆盖率输出、安装包以及根目录下自动生成的运行报告。

需要分享分析结果时，请从项目结果页导出，或使用 `sc_csv_export` 等交付模块生成独立结果包，不要把大型 `.h5ad`、原始矩阵或含样本信息的运行目录直接提交到 Git。

### 分析产物保留与磁盘空间

当前版本**不会自动删除**已完成任务的 `results/task_artifacts/`；这保证历史结果可追溯，但重复运行大型单细胞流程会累积 `.h5ad` 快照。管理员应定期监控 `DATA_DIR` 所在数据卷，并在人工确认后处理过期产物，不能直接删除仍被分支或任务结果引用的文件。

产物去重、快照保留和带审计的清理尚未实现，后续应按[单细胞分析产物保留策略设计](docs/single-cell-artifact-retention-design.md)分阶段落地。该文档中的 `SC_ARTIFACT_*` 配置项是设计提案，**并非当前版本可直接启用的功能**。

发布到 GitHub 前，维护者应确认所有新源码、模板、脚本和文档都已加入版本控制，而不是只
存在于开发机；下面三项检查不应出现意外输出：

```bash
git diff --check
git status --short
git ls-files --others --exclude-standard
```

第二、三条中出现的功能源码必须 `git add` 后再提交；`.env`、`data/`、`instance/` 和真实
人类数据仍应保持忽略。否则 GitHub 用户得到的克隆会缺少模块，即使本机工作树能运行。

## 配置项

| 配置 | 默认值 | 说明 |
| --- | --- | --- |
| `DATA_DIR` | `<repo>/data` | 项目数据、结果、presets 存储目录 |
| `RUNTIME_TMP_DIR` | `<repo>/data/runtime_tmp` | 临时分析文件、10x 兼容转换、Python/Numba/Matplotlib 缓存；应设在非根分区的数据盘 |
| `CACHE_DIR` | `<DATA_DIR>/cache` | Scanpy 等分析读取缓存目录 |
| `DB_PATH` | `<repo>/instance/bioinfo.db` | SQLite 数据库路径 |
| `BULK_REFERENCE_OUTPUT_DIR` | `<DATA_DIR>/bulk_reference_output` | Bulk 参考流程默认输出目录 |
| `CELLMARKER_PATH` | repo 上级目录下 `CellMarker_Augmented_2021.txt` | 兼容保留项；当前默认流程不要求克隆者提供该文件 |
| `SC_CELLTYPIST_MODEL_DIR` | `<repo>/data/references/celltypist` | 可选 CellTypist `.pkl` 模型的受控目录 |
| `FUNCTIONAL_STATE_RESOURCE_DIR` | 空（使用仓库 `resources/functional_state_resources`） | 管理员冻结的 Hallmark/Reactome/GO 与 CollecTRI 替代资源根目录 |
| `SC_CELL_GO_GENE_SET_DIR` | 空（使用同一仓库快照） | `sc_cell_go` 的管理员覆盖 GMT/TXT 目录 |
| `MAX_WORKERS` | `2` | 后台分析任务并发数 |
| `MIN_FREE_RAM_GB` | `4` | 资源保护阈值 |
| `CUDA_DEVICES` | `0,1` | GPU 设备配置 |
| `SC_BATCH_SOURCE_ROOTS` | 空（关闭） | 可供网页只读批量导入的服务器数据根目录；Linux 多个根用 `:` 分隔，例如 `/home/oelab/data/GJ:/mnt/sc_data` |
| `PLATFORM_ACCESS_PASSWORD` | 空 | 受控部署的共享访问门槛；公网/跨网访问必须设置并结合 HTTPS/VPN |
| `AI_API_KEY` | 空 | AI API key |
| `AI_API_URL` | 默认兼容 Anthropic 的 URL | AI API endpoint；DeepSeek 可用 `https://api.deepseek.com/anthropic` |
| `AI_MODEL` | `mimo-v2.5-pro` | AI 模型名；DeepSeek 常用 `deepseek-chat` |
| `AI_API_TOKEN` | 空 | 本地 API Bearer token |
| `WES_SOURCE_ROOTS` | 空（关闭） | WES 服务器只读数据根目录白名单，多个根用 `:` 分隔 |
| `WES_EXECUTOR_ENABLED` | 关闭 | 是否允许启动 Nextflow WES 运行；关闭时仍可准备和审阅 run |
| `WES_NEXTFLOW_GENOME` | `GATK.GRCh38` | Sarek/iGenomes 参考键；配套资源路径见 `WES_NEXTFLOW_*` 变量 |
| `CELLORACLE_PYTHON` | 本机开发路径（必须覆盖） | virtual_ko 使用的独立 CellOracle Python 3.9/3.10 解释器路径；新部署不可依赖该默认绝对路径 |

## 功能边界

- AI 能调用平台已暴露的工具，不能绕过模块代码本身的能力边界。
- 参数 sweep 会创建候选分支并消耗计算资源，最大候选数量有限制，写操作需要用户确认。
- marker 评分用于辅助判断“最接近某细胞类型的 cluster”，不是人工注释或实验验证的替代。
- LIANA、scVI、SysVI、DESeq2/edgeR/limma 等路径依赖对应包和环境；缺失时平台会在依赖接口或任务 summary 中说明具体不可用能力和下一步安装提示。
- 大型分析输出建议保留在 `data/` 或外部结果目录，不建议直接提交到 git。
