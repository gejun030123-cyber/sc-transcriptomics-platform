# 单细胞与 Bulk RNA-seq AI 分析平台

这是一个基于 Flask、Scanpy、Matplotlib 和 OmicVerse 风格规范的 Web 端转录组分析平台。平台同时覆盖单细胞转录组和 Bulk RNA-seq，支持传统按模块执行的分析流程，也支持用户通过 AI 对话触发分析、检查结果、调整参数，并围绕特定目标创建候选分支进行参数搜索。

项目当前定位不是单纯的流程封装，而是一个可交互的分析工作台：

- 用户可以上传数据，在网页中按模块执行 scRNA-seq 或 Bulk RNA-seq 分析。
- AI 助手可以读取项目状态、理解已完成步骤、调用分析工具并提出参数调整方案。
- 对“分群不满意”“必须找到最接近某种细胞类型的分群”等需求，平台提供 marker 评分、候选分支、参数 sweep 和人工采纳机制。
- 所有分析结果会落到项目目录下，包含 h5ad 中间文件、CSV 表格和科研级静态图。

## 核心能力

### Web 分析工作台

- 项目管理：创建项目、上传数据、查看项目状态和历史任务；主线任务完成后项目状态会自动同步为 `completed`，存在运行中任务时为 `processing`，全部失败时为 `failed`。
- 异步任务：分析任务通过后台 worker 执行，前端可查看进度、日志和失败信息。
- 参数面板：每个模块有结构化参数 schema，包含中文标签、默认值、类型和帮助说明。
- 结果管理：任务结果写入数据库，支持图表查看、CSV/XLSX 表格下载和 h5ad 中间文件下载；模块写出的每个结果文件都会登记到任务结果清单。
- 科研级图表：Scanpy/Matplotlib 原生图默认保存为 300 dpi PNG 和 SVG 矢量图，结果页优先展示 PNG，并提供 PNG/SVG 下载。
- 静态图表：统一使用 Matplotlib/OmicVerse 风格，网页默认展示 300 dpi PNG，并提供 SVG 矢量图下载。
- 图形工作台：从分析结果或用户上传图片创建非破坏性的样式版本，支持预览、PNG/SVG 导出和版本追踪。
- 可复现交付：任务 manifest、项目报告、离线图表图库和 pipeline 报告会记录参数、summary、产物和复核证据。
- 安全路径校验：API 读取和 AI 工具调用会限制在项目目录内，拒绝路径穿越和符号链接输入。

### 单细胞转录组分析

平台提供 14 个单细胞模块。核心顺序为：

```text
qc -> normalize -> hvg -> dimred -> batch_correct -> clustering
  -> qc_reassess / annotation / subcluster / sc_timecourse / deg
  -> trajectory / proportion / cell_communication
```

| 模块 | 主要功能 | 典型输出 |
| --- | --- | --- |
| `qc` | 线粒体、核糖体、血红蛋白比例标记；Scrublet 双细胞检测；细胞周期评分；复杂度指标；批次自适应 QC | QC violin、counts vs genes scatter、novelty plot、cell-cycle plot、过滤前后 QC 对比、doublet score 直方图 |
| `normalize` | `log1p` 或 Pearson residuals 标准化，保留 counts layer | 标准化后 h5ad、library size 图、表达值分布图 |
| `hvg` | 高变异基因选择，支持批次感知、force include、排除 MT/CC 基因、细胞周期评分和回归 | HVG scatter、HVG rank plot、HVG 标记 |
| `dimred` | PCA、UMAP，可选 t-SNE/MDE，支持自动 PC 选择 | PCA variance、PCA scatter、UMAP QC 着色图 |
| `batch_correct` | Harmony、ComBat、BBKNN、Scanorama、SysVI、scVI 等批次整合入口；支持 CPU/GPU 深度模型路径 | 校正后 embedding/graph、整合前后 UMAP、batch ASW、cluster batch entropy、最大批次占比、邻居混合、图连通性、指标表 |
| `clustering` | 多分辨率 Leiden/Louvain 聚类，支持主分辨率、自动分辨率评分 | 各分辨率 UMAP、多分辨率 UMAP、cluster 标签 UMAP、cluster 细胞数图、cluster 批次组成图、分辨率 Sankey |
| `subcluster` | 对指定 cluster 进行子簇重聚类、差异表达和通路富集 | 子簇 UMAP、marker 表、热图、富集结果 |
| `qc_reassess` | 聚类后按簇评估 doublet、MT、ribo、细胞数，支持自动移除低质量簇 | 低质量簇表、按簇 QC 汇总图、doublet/MT UMAP、QC 指标 UMAP 面板、低质量簇高亮图 |
| `annotation` | 分层 cell lineage/type/subtype、独立 cell state；marker 自动打分、负向 marker 互斥、多证据复核；类器官自动计算前体/成熟/增殖模块和成熟度指数，并读取时间元数据；可选本地人类 CellTypist 参考交叉验证（不覆盖 Marker 标签）；注释版本/备注可追溯；doublet/环境 RNA 仅作复核证据 | 细胞类型 UMAP、细胞类型组成图、marker score heatmap、marker dotplot、marker 表达验证图、annotation score UMAP、成熟度 UMAP、CellTypist 参考 UMAP、逐簇复核表 |
| `sc_timecourse` | 按真实时间点进行样本级细胞组成和伪 bulk 基因动态分析，区分描述性趋势与统计推断 | 时间点 UMAP、组成曲线/热图、动态基因表和趋势图 |
| `deg` | Wilcoxon、t-test、logreg 等 cluster/celltype 差异表达 | DEG CSV、完整 DEG CSV、火山图、显著 DEG 数量图、top marker UMAP 面板、dotplot、marker heatmap、基因表达 UMAP |
| `trajectory` | Diffusion Map、DPT、PAGA 拟时序 | pseudotime UMAP、pseudotime 分布图、PAGA 图、基因随拟时序变化图 |
| `proportion` | 细胞比例统计和组间比较，支持卡方、Fisher、置换检验 | 堆叠柱图、比例 heatmap、饼图、比例统计表 |
| `cell_communication` | 基于 LIANA 的配体-受体通讯分析 | 通讯热图、气泡图、通讯网络图、交互表 |

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
- 两批次 10x ZIP：上传两个分别包含 `filtered_feature_bc_matrix` 的 ZIP，平台会解压、添加 `obs["batch"]`，并合并为一个标准 `.h5ad`。

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

项目提供 `requirements.txt` 作为非锁定依赖清单，但没有 `pyproject.toml` 或锁定版本文件。建议在独立 conda/venv 环境中安装依赖，并根据需要补充可选包。

```bash
git clone https://github.com/gejun030123-cyber/sc-transcriptomics-platform.git
cd sc-transcriptomics-platform

pip install flask flask-cors scanpy anndata omicverse plotly psutil scipy statsmodels patsy scikit-learn matplotlib pandas numpy
pip install gseapy pydeseq2 inmoose
```

可选依赖：

- `celltypist`：启用本地人类 CellTypist 参考交叉验证；类器官注释不把外部模型作为默认真值，冲突结果保留人工复核。
- `liana`：启用细胞通讯分析；未安装时该任务会明确返回 `unavailable` 和安装提示，不会伪装成成功。
- `gseapy`：仅用于子簇通路富集和 Bulk 富集的兼容旧路径；未安装时核心聚类、差异分析和当前 OmicVerse/本地基因集富集路径仍可使用。
- `harmonypy`、`bbknn`、`scanorama`、`scvi-tools`、`torch`：启用 Harmony、BBKNN、Scanorama、SysVI/scVI 批次整合路径；未安装时对应方法不可用。
- `kneed`：启用 Kneedle 自动 PC 选择。
- `plotly`：仅用于兼容历史分析模块的内部数据结构；网页、报告和下载结果均不再输出 Plotly 交互图。
- `openai` 或 `anthropic`：启用对应 AI API 客户端。

启动服务：

```bash
python app.py
```

默认访问地址：

```text
http://localhost:5000
```

运行后可访问 `/api/system/dependencies` 查看依赖状态。每个模块同时返回 `missing`（缺少即不可运行的依赖）和 `optional_missing`（只影响某项扩展能力的依赖），因此页面或部署检查不应仅依据整组依赖是否全部安装来判断模块是否可用。

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

- `AI_API_KEY` 为空时，`/api/chat` 会返回未配置错误。
- `AI_API_URL` 中包含 `anthropic` 或 `claude` 时走 Anthropic Messages 格式，否则走 OpenAI compatible Chat Completions 格式。
- Anthropic-compatible 调用使用内置 HTTP 客户端，不强依赖本地安装 `anthropic` SDK。
- `AI_API_TOKEN` 为空时跳过本地 API token 认证；设置后需要请求头 `Authorization: Bearer <token>`。
- 写操作工具不会直接执行，会先返回 `proposed_tools`，前端确认后再调用 `/api/chat/approve`。

## 使用流程

### 单细胞网页流程

1. 创建项目。
2. 上传 `.h5ad`、10x `.h5`、10x 三文件、`.loom` 或 `.zarr` 数据，并在上传页导入为标准 `.h5ad`。
3. 执行 `qc`、`normalize`、`hvg`、`dimred`。
4. 执行 `clustering` 并查看多分辨率 UMAP、Sankey 和带标签 cluster UMAP。
5. 执行 `qc_reassess` 检查低质量簇。
6. 执行 `annotation`，选择 `PBMC`、`Immune`、`Blood`、`TME` 或自定义 marker。
7. 执行 `deg`，检查 marker heatmap、火山图和完整 DEG CSV。
8. 根据项目需要继续 `trajectory`、`proportion`、`cell_communication`。

### Bulk 网页流程

1. 创建项目并上传表达矩阵。
2. 执行 `bulk_qc`，确认样本质量和自动分组。
3. 执行 `bulk_normalize`。
4. 执行 `bulk_pca` 检查分组、批次和异常样本。
5. 执行 `bulk_deg`，设置比较组、统计方法和阈值。
6. 执行 `bulk_heatmap` 或 `bulk_enrichment`。
7. 多比较场景执行 `bulk_deg_integration`。

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
| `/api/chat/approve` | `POST` | 执行用户确认后的 AI 工具 |
| `/api/chat/history/<pid>` | `GET` | 获取项目聊天历史 |
| `/api/tasks/<task_id>/status` | `GET` | 查看分析任务状态 |
| `/api/projects/<pid>/tasks` | `GET` | 列出项目任务 |
| `/api/projects/<pid>/adata-info` | `GET` | 获取当前 AnnData 信息 |
| `/api/result-file/<file_id>` | `GET` | 下载 PNG、SVG、CSV 或其他结果文件 |
| `/projects/<pid>/figure-studio` | `GET` | 打开非破坏性图形工作台 |
| `/projects/<pid>/upload/import-10x-batches` | `POST` | 接收两个 `batch_zip` 和两个 `batch_name`，合并为带 `batch` 列的 h5ad |
| `/api/projects/<pid>/pipeline-runs` | `POST/GET` | 创建或列出批量 pipeline run |
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

运行全部测试：

```bash
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
  branches.py               # Agent session、candidate branch、score、accept API
  figure_studio.py          # 图形工作台页面、预览、上传和版本下载
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
  ai_tools.py               # AI 工具执行后端
  agent_orchestrator.py     # 目标驱动 Agent 编排
  agent_jobs.py             # 参数 sweep 异步 job
  design_preflight.py       # 分析前实验设计和比较组检查
  figure_style.py           # 统一 Nature 风格和 Matplotlib 参数
  native_figures.py         # 原生静态科研图构建器
  figure_studio.py          # 图形来源解析、预览和版本保存
  cell_markers.py           # 内置 marker 定义
  evaluators/               # 候选分支评分器
  reporting/                # manifest、项目报告、pipeline 和复核证据
  qc.py ... deg.py          # 单细胞分析模块
  bulk_qc.py ...            # Bulk 分析模块
  convert_10x.py            # 10x 转换模块

templates/                  # Jinja2 页面模板
tests/                      # 单元、集成、语义和 Agent 测试
scripts/                    # 参考流程和工具脚本
docs/                       # 设计、计划和审批文档
genesets/                   # 通路基因集资源
data/                       # 本地项目数据，默认不纳入 git
```

## 配置项

| 配置 | 默认值 | 说明 |
| --- | --- | --- |
| `DATA_DIR` | `<repo>/data` | 项目数据、结果、presets 存储目录 |
| `RUNTIME_TMP_DIR` | `<repo>/data/runtime_tmp` | 临时分析文件、10x 兼容转换、Python/Numba/Matplotlib 缓存；应设在非根分区的数据盘 |
| `DB_PATH` | `<repo>/instance/bioinfo.db` | SQLite 数据库路径 |
| `CELLMARKER_PATH` | repo 上级目录下 `CellMarker_Augmented_2021.txt` | 可选 CellMarker 数据 |
| `MAX_WORKERS` | `2` | 后台分析任务并发数 |
| `MIN_FREE_RAM_GB` | `4` | 资源保护阈值 |
| `CUDA_DEVICES` | `0,1` | GPU 设备配置 |
| `AI_API_KEY` | 空 | AI API key |
| `AI_API_URL` | 默认兼容 Anthropic 的 URL | AI API endpoint；DeepSeek 可用 `https://api.deepseek.com/anthropic` |
| `AI_MODEL` | `mimo-v2.5-pro` | AI 模型名；DeepSeek 常用 `deepseek-chat` |
| `AI_API_TOKEN` | 空 | 本地 API Bearer token |

## 功能边界

- AI 能调用平台已暴露的工具，不能绕过模块代码本身的能力边界。
- 参数 sweep 会创建候选分支并消耗计算资源，最大候选数量有限制，写操作需要用户确认。
- marker 评分用于辅助判断“最接近某细胞类型的 cluster”，不是人工注释或实验验证的替代。
- LIANA、scVI、SysVI、DESeq2/edgeR/limma 等路径依赖对应包和环境；缺失时平台会在依赖接口或任务 summary 中说明具体不可用能力和下一步安装提示。
- 大型分析输出建议保留在 `data/` 或外部结果目录，不建议直接提交到 git。

## License

MIT
