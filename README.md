# 单细胞 & Bulk RNA-seq 生信分析平台

基于 Flask + Scanpy + OmicVerse 的 Web 端生物信息学分析平台，支持单细胞转录组和 Bulk RNA-seq 全流程分析。

## 功能特性

### 单细胞转录组分析（12 个模块）

| 模块 | 功能 | 核心工具 |
|------|------|----------|
| 质控 (QC) | MT/ribo/hb 过滤、Scrublet 双细胞检测、细胞周期评分、QC 散点图 | scanpy, omicverse |
| 标准化 | log1p / Pearson 残差标准化，文库大小分布图 | scanpy, omicverse |
| 高变异基因 | HVG 选择（批次感知）、force-include、CC 评分/回归 | scanpy |
| 降维分析 | PCA + UMAP/t-SNE，自动 PC 选择（elbow/kneedle），方差比例图 | scanpy, sklearn |
| 批次校正 | Harmony / ComBat / SysVI / Scanorama / BBKNN / scVI | omicverse, inmoose |
| 聚类分析 | 多分辨率 Leiden/Louvain 聚类、UMAP 比较、自动分辨率选择 | scanpy |
| QC 重新评估 | 聚类后 doublet/MT 检测、低质量簇标记、自动移除 | scanpy |
| 细胞注释 | CellTypist / TME marker 基因集、自定义 marker、置信度评估 | scanpy, celltypist |
| 差异表达 | Wilcoxon/t-test/logreg、火山图、dotplot、基因 UMAP | scanpy |
| 轨迹分析 | Diffusion Map + DPT / PAGA 拟时序分析 | scanpy |
| 比例分析 | 细胞比例差异（卡方/Fisher/置换检验）、堆叠柱状图/饼图 | scipy |
| 细胞通讯 | LIANA 细胞间通讯分析、气泡图、热图 | liana |

### Bulk RNA-seq 分析（8 个模块）

| 模块 | 功能 | 核心工具 |
|------|------|----------|
| 数据质控 | 文库大小、基因检测、离群值过滤、分组相关性热图、Mahalanobis 距离异常检测 | scanpy, scipy |
| 数据标准化 | DESeq2 / TMM / CPM / VST / rlog / 分位数标准化 | scipy |
| PCA / UMAP | 降维可视化、载荷图、肘部图、样本聚类 | scanpy, sklearn |
| 差异表达分析 | t-test / Mann-Whitney / DESeq2 / edgeR / limma、LRT 检验、火山图、MA 图 | omicverse, inmoose |
| 热图分析 | Top 变异/差异基因热图、多种变异度量(MAD/CV)、z-score/中心化、灵活聚类 | scipy |
| 通路富集 | ORA / GSEA（GO/KEGG/WikiPathways/Reactome） | omicverse, gseapy |
| 时序分析 | 多时间点差异基因（spline F-test）、模糊 c-means 轨迹聚类 | patsy, statsmodels |
| 多比较整合 | UpSet / Venn 图、一致性评分、logFC 矩阵、表达式筛选器（AND/OR/NOT/XOR） | plotly, matplotlib |

### AI Chat 功能

- 内置 AI 助手（支持 Anthropic / OpenAI 兼容 API）
- 4 个工具：`run_analysis`（需确认）、`get_project_status`、`get_task_results`、`list_modules`
- Token 认证（`AI_API_TOKEN` 环境变量）

### 通用功能

- 参数提示系统：每个参数附带中文说明和建议范围
- 实时进度追踪：带时间戳的运行日志
- 分析流程引导：每步完成后自动推荐下一步
- 交互式图表：Plotly.js 支持缩放、平移、导出
- 多格式输入：支持 .h5ad, .csv, .txt, .xlsx, .xls, .mtx.gz
- 自动分组检测：从样本名中自动提取实验分组

## 技术栈

- **后端**: Flask, SQLite (WAL), ThreadPoolExecutor
- **分析**: scanpy, omicverse, scipy, statsmodels, patsy, scikit-learn, inmoose, gseapy
- **前端**: Bootstrap 5, Plotly.js, Jinja2
- **数据库**: SQLite（projects, analysis_tasks, result_files）

## 安装

```bash
# 克隆仓库
git clone https://github.com/gejun030123-cyber/sc-transcriptomics-platform.git
cd sc-transcriptomics-platform

# 安装依赖
pip install flask flask-cors scanpy omicverse plotly psutil patsy statsmodels gseapy pydeseq2 scikit-learn inmoose matplotlib

# 启动服务
python app.py
```

服务将在 `http://localhost:5000` 启动。

## 使用流程

1. 创建项目 → 上传数据文件
2. 选择分析类型（单细胞 / Bulk RNA-seq）
3. 按流程顺序执行各分析模块
4. 每步完成后查看结果图表和数据表格
5. 下载 h5ad 中间文件和 CSV 结果

### 单细胞典型流程

```
质控 → 标准化 → 高变异基因 → 降维 → 聚类 → QC 重新评估 → 细胞注释 → 差异表达 → 轨迹分析 → 比例分析 → 细胞通讯
```

### Bulk RNA-seq 典型流程

```
数据质控 → 标准化 → PCA/UMAP → 差异表达 → 热图 → 通路富集 / 时序分析 / 多比较整合
```

### 多比较整合表达式筛选语法

```
# 交集：两种药物共同上调
hmc3-vs-ctrl:up AND rapa-vs-ctrl:up

# 差集：仅 hmc3 上调，rapa 无变化
hmc3-vs-ctrl:up NOT rapa-vs-ctrl:up

# 简写：所有比较共同上调
ALL:up

# 阈值覆盖
hmc3-vs-ctrl:up(padj<<0.01|logFC>2)
```

支持运算符：`AND`/`∩`、`OR`/`∪`、`NOT`/`-`、`XOR`/`△`，支持括号嵌套。

## 测试

```bash
# 运行全部测试
python -m pytest tests/ -v

# 运行集成测试
python -m pytest tests/test_integration_full.py tests/test_integration.py -v

# 运行语义断言测试
python -m pytest tests/test_semantic.py tests/test_semantic_full.py -v
```

当前测试状态：**131 passed, 9 skipped, 0 failed**（15 个测试文件，4467+ 行测试代码）

| 测试类别 | 文件数 | 测试数 | 覆盖范围 |
|---------|--------|--------|---------|
| 单元测试 | 10 | ~170 | 核心工具函数、辅助函数、模型、schemas |
| 集成测试 | 2 | ~60 | 模块端到端、深度断言、流水线、Worker |
| 语义断言 | 2 | ~65 | 静态扫描 + 运行时一致性验证 |

## 项目结构

```
├── app.py                      # Flask 应用工厂
├── config.py                   # 全局配置
├── database.py                 # SQLite 连接管理 (WAL)
├── models.py                   # ORM 模型
├── worker.py                   # 异步任务执行器
├── routes/
│   ├── main.py                 # 首页控制台
│   ├── projects.py             # 项目 CRUD
│   ├── upload.py               # 文件上传
│   ├── analysis.py             # 分析模块配置
│   ├── results.py              # 结果展示
│   ├── api.py                  # REST API
│   └── chat.py                 # AI Chat
├── modules/
│   ├── base.py                 # 分析基类
│   ├── __init__.py             # 模块注册表 + 流水线
│   ├── schemas.py              # 参数定义 + UI 元数据
│   ├── constants.py            # 基因集常量
│   ├── io_utils.py             # 数据读取
│   ├── visualization.py        # Plotly 可视化工具
│   ├── expression_parser.py    # 筛选表达式解析器
│   ├── inspect_utils.py        # 模块检查工具
│   ├── ai_adapter.py           # AI Chat 适配器
│   ├── ai_tools.py             # AI 工具执行后端
│   ├── qc.py                   # 单细胞质控
│   ├── normalize.py            # 标准化
│   ├── hvg.py                  # 高变异基因
│   ├── dimred.py               # 降维
│   ├── batch_correct.py        # 批次校正
│   ├── clustering.py           # 聚类
│   ├── qc_reassess.py          # QC 重新评估
│   ├── annotation.py           # 细胞注释
│   ├── deg.py                  # 差异表达
│   ├── trajectory.py           # 轨迹分析
│   ├── proportion.py           # 比例分析
│   ├── cell_communication.py   # 细胞通讯
│   ├── convert_10x.py          # 10x 数据转换
│   ├── bulk_qc.py              # Bulk 质控
│   ├── bulk_normalize.py       # Bulk 标准化
│   ├── bulk_pca.py             # Bulk PCA
│   ├── bulk_deg.py             # Bulk DEG
│   ├── bulk_heatmap.py         # Bulk 热图
│   ├── bulk_enrichment.py      # 通路富集
│   ├── bulk_timecourse.py      # 时序分析
│   └── bulk_deg_integration.py # 多比较整合
├── tests/                      # 测试套件（15 个文件）
├── scripts/                    # 工具脚本
├── docs/                       # 架构文档、设计文档
├── genesets/                   # 通路基因集数据库
└── templates/                  # Jinja2 HTML 模板
```

## 配置

编辑 `config.py` 修改：

- `DATA_DIR`: 数据存储目录
- `DB_PATH`: SQLite 数据库路径
- `MAX_WORKERS`: 并行任务数
- `CUDA_DEVICES`: GPU 设备号（SysVI 批次校正用）
- `AI_API_TOKEN`: AI Chat 认证 Token

## License

MIT
