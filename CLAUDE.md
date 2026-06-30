# 生信分析平台 — 项目指南

## 项目概述

Flask 全栈 Web 平台，用于单细胞 RNA-seq 和 Bulk RNA-seq 转录组分析。
- **技术栈**：Python 3.10+, Flask, scanpy, omicverse, Plotly, SQLite
- **界面语言**：中文
- **入口**：`app.py` → `create_app()` → `0.0.0.0:5000`

## 目录结构

```
platform/
├── app.py              # Flask 应用工厂
├── config.py           # Config 类：路径、环境变量、AI 配置
├── database.py         # SQLite 连接管理 (WAL mode)
├── models.py           # ORM: Project, AnalysisTask, ResultFile
├── worker.py           # ThreadPoolExecutor 后台任务执行器
├── modules/            # 分析模块（21个）
│   ├── base.py         # BaseAnalysis 抽象基类
│   ├── schemas.py      # 参数定义 PARAM_SCHEMAS + 模块元数据
│   ├── __init__.py     # MODULE_REGISTRY, PIPELINE_ORDER, PIPELINE_DEPS
│   ├── ai_adapter.py   # AI Chat 适配器（Anthropic/OpenAI）
│   ├── ai_tools.py     # AI 工具执行后端
│   └── ...             # 21 个分析模块
├── routes/             # 7 个 Flask Blueprint
├── templates/          # Jinja2 HTML 模板
├── tests/              # pytest 测试
├── scripts/            # 工具脚本
└── docs/               # 架构文档、规划文档
```

## 模块系统

所有分析模块继承 `BaseAnalysis`，必须实现 `run(input_path) -> dict`。

**注册**：在 `modules/__init__.py` 的 `MODULE_REGISTRY` 中注册，key 为字符串标识。

**流水线顺序**：
- 单细胞：qc → normalize → hvg → dimred → batch_correct → clustering → qc_reassess → annotation → deg → trajectory → proportion → cell_communication
- Bulk：bulk_qc → bulk_normalize → bulk_pca → bulk_deg → bulk_heatmap → bulk_enrichment → bulk_timecourse → bulk_deg_integration

**依赖约束**：`PIPELINE_DEPS` 定义前置依赖，`validate_pipeline_order()` 校验。

## 数据流

```
用户上传 → uploads/ → [convert_10x] → uploads/*.h5ad
         → intermediate/{module}_output.h5ad（每步输出，下一步输入）
         → results/*.csv（分析结果表格）
         → plots/*.json（Plotly 可视化）
```

每个项目独立目录：`data/projects/{pid}/` 下的 `uploads/`, `intermediate/`, `results/`, `plots/`。

## 编码规范

- **语言**：代码注释用中文，变量/函数名用英文
- **模块开发**：继承 `BaseAnalysis`，参数定义在 `PARAM_SCHEMAS`
- **路由**：Blueprint 注册，单文件不超过 100 行（业务逻辑除外）
- **测试**：`tests/` 目录，`pytest tests/ -v` 运行
- **错误处理**：避免裸 `except Exception: pass`，至少记录日志
- **安全**：路径校验用 `Config._validate_pid()`，文件名用 `secure_filename()`

## 安全规则

- 所有项目 ID 必须匹配 `^[a-zA-Z0-9_-]{1,64}$`
- 文件下载前必须通过 `_validate_path()` 防止路径穿越
- AI Chat 支持 Token 认证 (`AI_API_TOKEN` 环境变量)
- `run_analysis` 工具调用需要用户确认（CONFIRM_TOOLS）

## 常用命令

```bash
# 启动服务
python app.py

# 运行测试
python -m pytest tests/ -v

# 查询任务状态
sqlite3 instance/bioinfo.db "SELECT id, module_name, status FROM analysis_tasks ORDER BY rowid DESC LIMIT 5;"

# 查询项目列表
sqlite3 instance/bioinfo.db "SELECT id, name, status FROM projects;"

# 查看系统状态
curl -s http://localhost:5000/api/system/status | python -m json.tool
```

## AI 集成

平台内置 AI Chat 功能（`routes/chat.py` + `modules/ai_adapter.py`）：
- 支持 Anthropic 和 OpenAI 兼容 API
- 4 个工具：`run_analysis`（需确认）、`get_project_status`（自动）、`get_task_results`（自动）、`list_modules`（自动）
- 工具调用循环最多 4-5 轮
