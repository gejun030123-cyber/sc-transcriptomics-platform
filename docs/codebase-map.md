# 项目代码结构与关系梳理

**项目：** 单细胞 & Bulk RNA-seq 生信分析平台
**技术栈：** Flask + scanpy + omicverse + Plotly + SQLite
**入口：** `app.py` → `create_app()` → `0.0.0.0:5000`

---

## 一、目录总览

```
platform/
├── app.py                    # Flask 入口，工厂模式 create_app()
├── config.py                 # 全局配置（环境变量 + 路径工具）
├── database.py               # SQLite 连接管理 + schema 初始化
├── models.py                 # ORM 模型（Project, AnalysisTask, ResultFile）
├── worker.py                 # 后台任务执行器（ThreadPoolExecutor）
│
├── routes/                   # Flask Blueprint 路由层
│   ├── main.py               # 首页仪表盘
│   ├── projects.py           # 项目 CRUD
│   ├── upload.py             # 文件上传 + 10x 格式转换
│   ├── analysis.py           # 分析任务提交（瘦路由）
│   ├── results.py            # 结果查看 + 文件下载
│   ├── api.py                # REST API（任务状态、数据检查、预设）
│   └── chat.py               # AI 对话端点
│
├── modules/                  # 分析模块层
│   ├── __init__.py           # MODULE_REGISTRY 注册表
│   ├── base.py               # BaseAnalysis 抽象基类（模板方法）
│   ├── constants.py          # 共享常量（细胞周期基因列表）
│   ├── schemas.py            # 参数 schema + 模块元数据
│   ├── io_utils.py           # 数据读取工具（h5ad/CSV/10x）
│   ├── inspect_utils.py      # anndata 结构检查工具
│   ├── visualization.py      # Plotly 可视化工具
│   ├── expression_parser.py  # 基因集合表达式 DSL 解析器
│   ├── ai_adapter.py         # AI 对话适配器（Anthropic + OpenAI）
│   ├── ai_tools.py           # AI 工具函数调度
│   │
│   ├── qc.py                 # 单细胞：质控
│   ├── normalize.py          # 单细胞：标准化
│   ├── hvg.py                # 单细胞：高变异基因
│   ├── dimred.py             # 单细胞：降维（PCA/UMAP/t-SNE）
│   ├── batch_correct.py      # 单细胞：批次校正
│   ├── clustering.py         # 单细胞：聚类
│   ├── qc_reassess.py        # 单细胞：QC 重评估
│   ├── annotation.py         # 单细胞：细胞注释
│   ├── deg.py                # 单细胞：差异表达
│   ├── trajectory.py         # 单细胞：轨迹分析
│   ├── proportion.py         # 单细胞：比例分析
│   ├── cell_communication.py # 单细胞：细胞通讯
│   ├── convert_10x.py        # 10x Genomics 格式转换
│   │
│   ├── bulk_qc.py            # Bulk：质控
│   ├── bulk_normalize.py     # Bulk：标准化
│   ├── bulk_deg.py           # Bulk：差异表达
│   ├── bulk_pca.py           # Bulk：PCA/UMAP
│   ├── bulk_heatmap.py       # Bulk：热图
│   ├── bulk_enrichment.py    # Bulk：通路富集
│   ├── bulk_timecourse.py    # Bulk：时序分析
│   └── bulk_deg_integration.py  # Bulk：多组差异整合
│
├── templates/                # Jinja2 模板
├── tests/                    # pytest 测试（4 个文件，87 个用例）
├── genesets/                 # 基因集数据库（GO/Reactome/WikiPathways）
├── data/                     # 运行时数据（项目文件、预设、样本）
├── instance/                 # SQLite 数据库文件
└── cache/                    # h5ad 缓存
```

---

## 二、架构分层

```
Layer 4  ┌─────────────────────────────────────────────────────┐
(入口)   │  app.py — 创建 Flask app，注册 7 个 Blueprint        │
         └───────────────────────┬─────────────────────────────┘
                                 │
Layer 3  ┌───────────────────────┴─────────────────────────────┐
(路由)   │  routes/ — 7 个 Blueprint                           │
│        │  main | projects | upload | analysis | results |    │
│        │  api | chat                                         │
│        └───────┬──────────┬──────────┬──────────┬────────────┘
                │          │          │          │
Layer 2  ┌──────┴────┐ ┌───┴────┐ ┌───┴────┐ ┌───┴──────────┐
(服务)   │ models.py │ │worker  │ │schemas │ │ai_adapter    │
│        │ Project   │ │submit  │ │PARAMS  │ │chat()        │
│        │ Task      │ │_task() │ │MODULES │ │ai_tools      │
│        │ ResultFile│ │        │ │        │ │execute_tool()│
│        └─────┬─────┘ └───┬────┘ └───┬────┘ └──────┬───────┘
              │            │          │             │
Layer 1  ┌────┴────────────┴──────────┴─────────────┴──────────┐
(数据)   │  database.py — get_conn() 线程级连接                 │
│        │  config.py — Config 类 + 路径工具方法                 │
│        └─────────────────────┬───────────────────────────────┘
                              │
Layer 0  ┌────────────────────┴───────────────────────────────┐
(工具)   │  modules/io_utils.py       — 数据读取               │
         │  modules/inspect_utils.py  — anndata 检查           │
         │  modules/visualization.py  — Plotly 绑定            │
         │  modules/expression_parser.py — DSL 解析器          │
         │  modules/constants.py      — 共享常量               │
         └────────────────────────────────────────────────────┘
```

---

## 三、核心文件详解

### 3.1 入口层

#### `app.py`（37 行）

Flask 工厂函数，注册所有 Blueprint 并初始化数据库。

```
create_app()
  ├── Config 加载
  ├── CORS 启用
  ├── init_db() → database.py
  └── register_blueprint ×7
        ├── main_bp        (无前缀)
        ├── projects_bp    (/projects)
        ├── upload_bp      (/projects)
        ├── analysis_bp    (/projects)
        ├── results_bp     (/projects)
        ├── api_bp         (/api)
        └── chat_bp        (无前缀，内部路由 /api/chat)
```

---

### 3.2 配置层

#### `config.py`（~35 行）

纯配置类 + 路径工具方法，无任何外部依赖（叶节点）。

```
Config
  ├── 属性：DATA_DIR, DB_PATH, MAX_WORKERS, AI_API_KEY, ...
  └── 路径方法（@classmethod）
        ├── project_dir(pid)       → DATA_DIR/projects/{pid}
        ├── uploads_dir(pid)       → .../uploads
        ├── results_dir(pid)       → .../results
        ├── intermediate_dir(pid)  → .../intermediate
        └── plots_dir(pid)         → .../plots
```

**被引用：** 几乎所有文件（15+ 处）

---

### 3.3 数据层

#### `database.py`（~55 行）

线程级 SQLite 连接管理 + schema 定义。

```
database.py
  ├── get_conn()    → 每线程一个连接（threading.local）
  └── init_db()     → CREATE TABLE IF NOT EXISTS ×3 + INDEX ×4
```

**依赖：** `config.py`
**被引用：** `models.py`

#### `models.py`（~180 行）

自研 ORM，3 个模型类，每个都有 CRUD 方法。

```
models.py
  ├── gen_id()                → uuid4[:12]
  │
  ├── Project
  │     ├── save(), delete(), to_dict()
  │     ├── get_by_id(), get_all()
  │     ├── get_tasks()       → AnalysisTask.get_by_project()
  │     └── get_latest_adata_path()
  │
  ├── AnalysisTask
  │     ├── save(), to_dict()
  │     ├── get_by_id(), get_by_project()
  │     ├── mark_running()
  │     ├── update_progress(pct, msg, log)
  │     ├── mark_completed(output_adata, result_json)
  │     └── mark_failed(error_traceback)
  │
  └── ResultFile
        ├── save(), to_dict()
        ├── create()           → 便捷构造+保存
        ├── get_by_id(), get_by_task(), get_by_project()
```

**依赖：** `database.py`
**被引用：** 所有路由文件、`worker.py`、`ai_tools.py`、`bulk_deg_integration.py`

---

### 3.4 任务执行层

#### `worker.py`（~65 行）

后台任务执行器，使用 ThreadPoolExecutor(max_workers=2)。

```
worker.py
  ├── submit_task(task_id, project_id, module_name, params, project_dir, input_path)
  │     └── _executor.submit(_run_task, ...)
  │
  ├── _run_task(...)
  │     ├── AnalysisTask.get_by_id() → 验证任务存在
  │     ├── task.mark_running()
  │     ├── MODULE_REGISTRY[module_name] → 实例化模块
  │     ├── module.run(input_path) → 执行分析
  │     ├── ResultFile.create() × N → 保存结果文件
  │     ├── task.mark_completed() 或 task.mark_failed()
  │     └── progress_cb() → 实时更新进度
  │
  └── active_count() → 活跃任务数
```

**依赖：** `database.py`, `models.py`, `modules/__init__.py`（延迟导入）
**被引用：** `routes/analysis.py`, `routes/upload.py`, `ai_tools.py`

---

### 3.5 路由层

#### `routes/main.py`（~20 行）

首页仪表盘，显示系统状态和项目列表。

```
GET / → index()
  ├── Project.get_all()
  ├── psutil 系统指标（CPU/内存/磁盘）
  └── render_template('index.html')
```

#### `routes/projects.py`（~53 行）

项目 CRUD。

```
GET/POST /projects/new      → new()        创建项目 + 目录
GET      /projects/<pid>    → detail()     项目详情 + 文件列表
POST     /projects/<pid>/delete → delete() 删除项目 + shutil.rmtree
```

**依赖：** `models.py`, `config.py`, `modules/schemas.py`

#### `routes/upload.py`（~143 行）

文件上传 + 10x 格式检测/转换。

```
GET/POST /projects/<pid>/upload           → upload()      上传文件
GET      /projects/<pid>/upload/check-10x → check_10x()   检测 10x 文件
POST     /projects/<pid>/upload/convert-10x → convert_10x() 转换任务
```

**依赖：** `models.py`, `config.py`, `worker.py`（延迟导入）

#### `routes/analysis.py`（~100 行，优化后）

分析任务提交。**瘦路由** — 参数解析和文件扫描委托给 `schemas.py`。

```
GET /projects/<pid>/sc-analysis           → 单细胞模块列表
GET /projects/<pid>/bulk-analysis         → Bulk 模块列表
GET/POST /projects/<pid>/analyze/<module> → 参数表单 + 提交任务
  ├── parse_form_params(schema, form)     → modules/schemas.py
  ├── list_upload_files(pid, is_bulk)     → modules/schemas.py
  ├── AnalysisTask.save()
  └── submit_task() → worker.py
```

**依赖：** `models.py`, `config.py`, `modules/schemas.py`, `worker.py`（延迟导入）

#### `routes/results.py`（~93 行）

结果查看和文件下载。

```
GET /projects/<pid>/task/<task_id>        → task_detail()   结果展示
GET /projects/<pid>/results               → results_gallery()
GET /projects/<pid>/results/file/<fid>    → view_file()     文件查看
GET /projects/<pid>/download/<task_id>    → download_adata() 下载 h5ad
```

**依赖：** `models.py`, `config.py`, `modules/schemas.py`, `modules/__init__.py`

#### `routes/api.py`（~550 行，优化后）

REST API 端点集合。

```
GET  /api/system-status              → 系统状态（CPU/内存/任务数）
GET  /api/tasks/<tid>/status         → 任务进度查询
GET  /api/projects/<pid>/adata-info  → anndata 结构检查（用 inspect_utils）
GET  /api/projects/<pid>/column-values → obs 列值查询
GET  /api/projects/<pid>/deg-comparisons → DEG 比较列表
POST /api/validate-filter-expression → 表达式验证
GET  /data-info-full                 → 完整 h5ad 结构（用 inspect_utils）
GET  /api/presets                    → 预设管理（CRUD）
```

**依赖：** `models.py`, `config.py`, `worker.py`, `modules/inspect_utils.py`, `modules/io_utils.py`, `modules/expression_parser.py`

#### `routes/chat.py`（~70 行）

AI 对话端点。

```
GET  /api/chat        → 聊天页面
POST /api/chat        → 发送消息（调用 ai_adapter.chat()）
GET  /api/chat/history → 聊天历史
DELETE /api/chat/history → 清除历史
```

**依赖：** `config.py`, `modules/ai_adapter.py`（延迟导入）

---

### 3.6 分析模块层

#### `modules/__init__.py`（~84 行）

模块注册表，集中导入所有 21 个分析类。

```
MODULE_REGISTRY = {
    'qc': QCAnalysis,
    'normalize': NormalizeAnalysis,
    'hvg': HVGAnalysis,
    ...（共 21 个）
}

PIPELINE_ORDER = ['qc', 'normalize', 'hvg', 'dimred', 'batch_correct',
                  'clustering', 'qc_reassess', 'annotation', 'deg',
                  'trajectory', 'proportion', 'cell_communication']

PIPELINE_DEPS = { 'normalize': ['qc'], 'hvg': ['normalize'], ... }
```

#### `modules/base.py`（~180 行）

抽象基类，提供模板方法。

```
BaseAnalysis(ABC)
  ├── 类属性：MODULE_NAME, DISPLAY_NAME, DESCRIPTION, INPUT_REQUIRES, PARAM_SCHEMA
  ├── __init__(project_dir, params, progress_callback)
  │
  ├── 模板方法（供子类调用）
  │     ├── load_adata(input_path)        → 读取 h5ad + remap 基因名
  │     ├── save_output(adata, module_name) → 保存到 intermediate/
  │     ├── ensure_plots_dir()             → 确保 plots/ 存在
  │     └── save_plotly_json(fig, ...)     → 保存 Plotly JSON + 返回 result_file dict
  │
  ├── 工具方法
  │     ├── progress(pct, message)
  │     ├── apply_filters(adata, module_name)  → 声明式过滤
  │     ├── get_plotly_layout(title, **overrides)
  │     ├── get_viz_params()
  │     ├── export_static_fig(fig, ...)    → SVG/PNG 导出
  │     └── export_results(adata, ...)     → h5ad/CSV/loom 导出
  │
  └── 抽象方法
        ├── validate_input(adata) → 默认返回 None
        └── run(input_path) → dict   子类必须实现
```

#### `modules/constants.py`（~20 行，叶节点）

Seurat 细胞周期基因列表，被 `qc.py` 和 `hvg.py` 共享。

```
S_GENES    = ['MCM5', 'PCNA', ...]   # 43 个 S 期基因
G2M_GENES  = ['HMGB2', 'CDK1', ...]  # 53 个 G2M 期基因
```

#### `modules/schemas.py`（~415 行，叶节点）

模块元数据和参数 schema，从 `analysis.py` 提取。

```
SC_MODULE_LIST      = [{name, display, desc}, ...]   # 12 个单细胞模块
BULK_MODULE_LIST    = [{name, display, desc}, ...]   # 8 个 Bulk 模块
MODULE_DISPLAY_MAP  = {'qc': '质控', ...}
STATUS_MAP          = {'pending': '等待中', ...}
PARAM_SCHEMAS       = {'qc': [{key, label, type, default, help}, ...], ...}

parse_form_params(schema, form)  → dict    参数解析
list_upload_files(pid, is_bulk)  → list    文件列表
```

#### `modules/io_utils.py`（~253 行，叶节点）

数据读取工具，**18 个外部调用点**，项目最高扇出。

```
read_expression_matrix(file_path)  → AnnData
  自动检测格式：h5ad / Excel / TSV / CSV
  处理 count-vs-FPKM 列去重、NaN 替换、基因名 remap

remap_var_names(adata)  → AnnData
  Ensembl ID → 基因符号映射
  被 13 个单细胞模块调用

convert_10x_to_h5ad(mtx_dir, output_path, species, genome)  → AnnData
  scanpy.read_10x_mtx 包装
```

#### `modules/inspect_utils.py`（~75 行，叶节点）

anndata 结构检查，被 `api.py` 的 `adata_info` 和 `data_info_full` 共用。

```
inspect_adata(adata, include_samples=True)  → dict
  提取 obs/var 列信息、layers、obsm、varm、obsp、uns 结构
```

#### `modules/visualization.py`（~217 行，叶节点）

Plotly 可视化工具集。

```
umap_scatter(adata, color, title)       → dict   UMAP 散点图
violin_plot(adata, keys, title)         → dict   小提琴图
scatter_plot(df, x, y, color, title)    → dict   通用散点图
bar_plot(df, x, y, title)               → dict   柱状图
compute_gene_variability(adata)         → DataFrame
transform_heatmap_data(adata, genes)    → DataFrame
cluster_heatmap(data, ...)              → dict   热图
build_annotation_bar(adata, cols)       → dict   注释条
save_plotly_json(fig, path)             → None   保存 JSON
```

#### `modules/expression_parser.py`（~417 行，叶节点）

基因集合表达式 DSL 解析器，支持 AND/OR/NOT/XOR。

```
tokenize(expr)  → tokens
parse(tokens)   → AST
evaluate(ast, sets)  → set
validate(expr)  → bool

AST 节点：Atom, BinOp, UnaryNot, Shorthand
```

#### `modules/ai_adapter.py`（~333 行）

AI 对话适配器，支持 Anthropic 和 OpenAI 兼容 API。

```
chat(messages, project_id=None)  → str
  ├── 自动检测 provider（URL 匹配）
  ├── Anthropic 路径：多轮 tool-calling（最多 4 轮）
  ├── OpenAI 路径：多轮 tool-calling（最多 5 轮）
  └── 文本 tool-call 解析回退
```

**依赖：** `config.py`, `modules/ai_tools.py`（延迟导入）

#### `modules/ai_tools.py`（~178 行）

AI 工具函数调度器。

```
execute_tool(name, args, project_id)  → dict
  ├── run_analysis     → AnalysisTask + submit_task
  ├── get_project_status → 任务状态汇总
  ├── get_task_results → 结果文件列表
  └── list_modules     → MODULE_REGISTRY 键列表
```

**依赖：** `config.py`, `models.py`, `worker.py`, `modules/__init__.py`（全部延迟导入）

---

### 3.7 分析模块详情

#### 单细胞模块（12 个）

全部继承 `BaseAnalysis`，遵循统一模式：

```python
def run(self, input_path):
    adata = self.load_adata(input_path)      # 加载 + remap
    # ... 分析逻辑 ...
    plots_dir = self.ensure_plots_dir()
    result_files = []
    result_files.append(self.save_plotly_json(...))
    output_path = self.save_output(adata, MODULE_NAME)
    return {'output_adata': output_path, 'result_files': result_files, 'summary': {...}}
```

| 模块 | 文件 | 行数 | 核心功能 | 特殊依赖 |
|------|------|------|----------|----------|
| qc | qc.py | ~260 | MT/ribo/hb 过滤 + Scrublet + 细胞周期 | constants (S/G2M genes) |
| normalize | normalize.py | ~70 | log1p / Pearson 残差标准化 | omicverse |
| hvg | hvg.py | ~150 | 高变异基因选择 + 批次感知 | constants (S/G2M genes) |
| dimred | dimred.py | ~140 | PCA/UMAP/t-SNE + 自动 PC 选择 | kneed（可选） |
| batch_correct | batch_correct.py | ~170 | Harmony/ComBat/SysVI/scVI/BBKNN | scvi-tools（可选） |
| clustering | clustering.py | ~155 | Leiden/Louvain 多分辨率聚类 | annotation.DEFAULT_TME_MARKERS |
| qc_reassess | qc_reassess.py | ~160 | 聚类后 QC 重评估 | — |
| annotation | annotation.py | ~260 | 分层 cell lineage/type/subtype；Colorectal 两级谱系/上皮亚型决策；log-normalized multi-label cell state；继承 QC/Scrublet Doublet、异源 ambient 与逐簇决策证据；可选本地 CellTypist 参考 | scanpy, celltypist（可选） |
| deg | deg.py | ~200 | Wilcoxon 差异表达 + 火山图 | statsmodels |
| trajectory | trajectory.py | ~150 | 拟时序 + PAGA + 扩散图 | — |
| proportion | proportion.py | ~170 | 细胞比例分析 + 统计检验 | scipy |
| cell_communication | cell_communication.py | ~150 | LIANA 细胞通讯 | liana（可选） |

#### Bulk RNA-seq 模块（8 个）

使用 `read_expression_matrix()` 读取数据（非 `sc.read_h5ad`）。

| 模块 | 文件 | 行数 | 核心功能 |
|------|------|------|----------|
| bulk_qc | bulk_qc.py | ~280 | 文库大小/基因检测/离群值过滤 |
| bulk_normalize | bulk_normalize.py | ~300 | DESeq2/TMM/CPM/VST/rlog |
| bulk_deg | bulk_deg.py | ~650 | 差异表达（t-test/DESeq2/edgeR/limma） |
| bulk_pca | bulk_pca.py | ~180 | PCA/UMAP/t-SNE 降维 |
| bulk_heatmap | bulk_heatmap.py | ~350 | 热图 + 聚类 + 注释条 |
| bulk_enrichment | bulk_enrichment.py | ~250 | ORA/GSEA 通路富集 |
| bulk_timecourse | bulk_timecourse.py | ~300 | 时序分析 + spline + 聚类 |
| bulk_deg_integration | bulk_deg_integration.py | ~380 | UpSet 图 + 一致性评分 |

---

## 四、数据流

### 4.1 用户上传 → 分析 → 结果

```
用户浏览器
  │
  ▼
routes/upload.py::upload()
  │ POST /projects/<pid>/upload
  │ 文件保存到 data/projects/<pid>/uploads/
  ▼
routes/analysis.py::analyze()
  │ POST /projects/<pid>/analyze/<module>
  │ parse_form_params() → 构建 params dict
  │ AnalysisTask(project_id, module_name, params_json).save()
  │ submit_task() → worker.py
  ▼
worker.py::_run_task()  [ThreadPoolExecutor 后台线程]
  │ task.mark_running()
  │ MODULE_REGISTRY[module_name] → 实例化模块
  │ module.run(input_path)
  │   ├── load_adata() → 读取 h5ad
  │   ├── ... 分析逻辑 ...
  │   ├── save_plotly_json() → plots/
  │   └── save_output() → intermediate/
  │ ResultFile.create() × N → 写入数据库
  │ task.mark_completed()
  ▼
routes/results.py::task_detail()
  │ GET /projects/<pid>/task/<task_id>
  │ ResultFile.get_by_task() → 获取结果文件
  │ render_template('analysis_result.html')
  ▼
用户浏览器（Plotly 图表 + CSV 表格）
```

### 4.2 AI 对话流

```
用户浏览器
  │ POST /api/chat {message, project_id}
  ▼
routes/chat.py::chat_api()
  │ ai_adapter.chat(messages, project_id)
  ▼
modules/ai_adapter.py::chat()
  │ 自动检测 provider
  │ 构建 system prompt + tools
  │ 多轮 tool-calling 循环（最多 4-5 轮）
  │   ├── LLM API 调用
  │   ├── execute_tool(name, args, project_id) → ai_tools.py
  │   │     ├── run_analysis → submit_task() → worker.py
  │   │     ├── get_project_status → AnalysisTask 查询
  │   │     └── get_task_results → ResultFile 查询
  │   └── 将工具结果返回 LLM
  ▼
返回 AI 回复文本
```

### 4.3 数据库连接流

```
所有数据库操作 → get_conn() [database.py]
  │ threading.local() 每线程一个连接
  │ 自动设置 PRAGMA journal_mode=WAL, foreign_keys=ON
  │
  ├── Flask 请求线程 → routes/* → models.py → get_conn()
  ├── Worker 线程    → worker.py → models.py → get_conn()
  └── init_db()      → app.py → database.py → get_conn()
```

---

## 五、模块依赖图

### 5.1 核心依赖链

```
config.py (叶节点)
  └── database.py
        └── models.py
              ├── worker.py
              │     └── modules/__init__.py (MODULE_REGISTRY)
              └── routes/*.py

modules/io_utils.py (叶节点) ─────────────┐
modules/visualization.py (叶节点) ────────┤
modules/constants.py (叶节点) ────────────┤
modules/expression_parser.py (叶节点) ────┤
modules/inspect_utils.py (叶节点) ────────┤
                                          ▼
                              modules/base.py
                                    │
                    ┌───────────────┼───────────────┐
                    ▼               ▼               ▼
              12 个单细胞模块   8 个 Bulk 模块   convert_10x
```

### 5.2 跨模块耦合

```
clustering.py ──imports──→ annotation.DEFAULT_TME_MARKERS
bulk_heatmap.py ──imports──→ expression_parser (validate, evaluate)
bulk_heatmap.py ──imports──→ visualization (5 个函数)
bulk_deg_integration.py ──imports──→ expression_parser + models
ai_adapter.py ──imports──→ ai_tools.execute_tool
ai_tools.py ──imports──→ models + worker + MODULE_REGISTRY
```

---

## 六、测试覆盖

| 文件 | 测试内容 | 用例数 |
|------|----------|--------|
| test_expression_parser.py | DSL 解析器（tokenize/parse/evaluate/validate） | 47 |
| test_heatmap_helpers.py | 热图工具函数 | 24 |
| test_bulk_qc_helpers.py | Gini 系数、分组推断、离群检测 | 10 |
| test_normalize_helpers.py | TMM/VST/rlog 标准化 | 8 |

**未覆盖：** 17/21 个分析模块、所有路由、worker、models、database

---

## 七、文件行数统计

| 类别 | 文件数 | 总行数 |
|------|--------|--------|
| 入口 + 配置 + 数据 | 4 | ~290 |
| 路由 | 7 | ~750 |
| 模块工具 | 7 | ~1,300 |
| 单细胞模块 | 12 | ~2,100 |
| Bulk 模块 | 8 | ~2,700 |
| 测试 | 4 | ~780 |
| **合计** | **42** | **~7,920** |
