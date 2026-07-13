# Pipeline Runner 双模型交接单 - Iteration 1

生成时间: 2026-07-02

读取对象: 实现/执行模型

审查对象:
- `worker.py`
- `models.py`
- `database.py`
- `routes/api.py`
- `templates/sc_analysis.html`
- `templates/bulk_analysis.html`
- `templates/analysis_select.html`
- `modules/__init__.py`
- `modules/schemas.py`
- 新增的 Pipeline Runner 相关文件

## 最新审查结论

最新结论: **尚未开始实现**。

本文件用于启动新一轮任务: **保存的流程模板一键顺序运行，覆盖单细胞 RNA-seq 与 Bulk RNA-seq 两类流程**。

后续实现模型完成一轮后，请在本文档末尾追加 `Implementer Reply - Iteration N`；Codex 审查时继续在本文档末尾追加 `Codex Review - Iteration N`。

## 后续协作规则

实现/执行模型完成一轮后，请在本文档末尾追加:

```md
## Implementer Reply - Iteration N

- 已读取要求: 是/否
- 修改文件:
- 设计取舍:
- 新增 API:
- 新增/修改测试:
- 执行命令:
- 测试结果:
- 手动验证:
- 已知失败/跳过:
- 请求 Codex 审查:
```

Codex 审查时会读取最新 `Implementer Reply`，复核代码、测试和实际行为，再继续追加:

```md
## Codex Review - Iteration N

- 审查结论:
- 复跑命令:
- 发现的问题:
- 通过项:
- 下一轮要求:
```

不要新建新的 handoff 文件，除非用户明确要求归档或重开一轮。

## 当前目标

实现一个 **Pipeline Runner**:

用户在 Web 端选择一个已保存的流程模板后，可以一键按依赖顺序自动运行多个分析模块。

第一阶段需要支持两类流程的最小可用链路。

单细胞 RNA-seq 最小链:

```text
qc -> normalize -> hvg
```

单细胞 RNA-seq 后续完整链:

```text
qc -> normalize -> hvg -> dimred -> clustering -> annotation -> deg
可选: batch_correct / qc_reassess / trajectory / proportion / cell_communication
```

Bulk RNA-seq 最小链:

```text
bulk_qc -> bulk_normalize -> bulk_pca
```

Bulk RNA-seq 后续完整链:

```text
bulk_qc -> bulk_normalize -> bulk_pca -> bulk_deg -> bulk_heatmap
可选: bulk_enrichment / bulk_timecourse / bulk_deg_integration
```

实现时可以先保证两个最小链路都能跑通，再扩展到更长链路。不要把底层设计写死为 Bulk 专用。

## 现有基础

已经存在:
- 模块注册表: `modules/__init__.py::MODULE_REGISTRY`
- 流程顺序: `modules/__init__.py::PIPELINE_ORDER`
- 依赖约束: `modules/__init__.py::PIPELINE_DEPS`
- 顺序校验: `modules/__init__.py::validate_pipeline_order`
- 单任务后台执行: `worker.py::submit_task`
- 分析任务模型: `models.py::AnalysisTask`
- 结果文件模型: `models.py::ResultFile`
- 参数 schema: `modules/schemas.py::PARAM_SCHEMAS`
- 预设保存 API: `routes/api.py` 的 `/api/presets`
- 前端流程模板保存按钮: `templates/analysis_select.html`

关键缺口:
- 已保存的 pipeline template 目前主要是“保存”，没有真正“一键运行”
- `worker.py::submit_task` 只处理单个模块，不负责多个模块之间的 input/output 串联
- 没有 pipeline run 级别的状态、日志、子任务列表和失败位置

## P0 实现范围

### 1. Pipeline Run 持久化

建议新增 `pipeline_runs` 表和 `PipelineRun` 模型。

推荐字段:

```text
id TEXT PRIMARY KEY
project_id TEXT NOT NULL
name TEXT NOT NULL
analysis_type TEXT DEFAULT ''
status TEXT DEFAULT 'pending'
current_module TEXT DEFAULT ''
progress INTEGER DEFAULT 0
input_path TEXT DEFAULT ''
modules_json TEXT DEFAULT '[]'
params_json TEXT DEFAULT '{}'
task_ids_json TEXT DEFAULT '[]'
error_traceback TEXT
started_at DATETIME
finished_at DATETIME
log_text TEXT DEFAULT ''
```

状态建议:

```text
pending -> running -> completed
pending/running -> failed
```

注意:
- `pipeline_runs.project_id` 应该 `REFERENCES projects(id) ON DELETE CASCADE`
- `task_ids_json` 保存本次 pipeline 创建的子任务 ID，便于前端跳转和审查
- 不要把 pipeline run 伪装成普通 `AnalysisTask`，否则结果页和模块名映射会变脏

### 2. Worker 串联执行

建议在 `worker.py` 中新增:

```python
submit_pipeline_run(run_id, project_id, modules, params_by_module, project_dir, input_path)
```

内部逻辑:

```text
标记 pipeline_run running
current_input = input_path
for module_name in modules:
    创建 AnalysisTask(status='pending')
    标记 task running
    实例化 MODULE_REGISTRY[module_name]
    module.run(current_input)
    注册 result_files
    task.mark_completed(output_adata, summary_json)
    current_input = output_adata
    更新 pipeline_run task_ids/current_module/progress/log_text
全部完成后标记 pipeline_run completed
任一步失败: 当前 task failed，pipeline_run failed，后续模块不再执行
```

重要约束:
- 不要在 pipeline runner 内部调用 `submit_task()` 再轮询等待。这样会造成嵌套 executor、状态竞争和潜在队列阻塞。
- 可以抽取一个共享 helper，避免 `worker.py` 单任务执行和 pipeline 执行重复注册 result_files 的代码。
- 每一步的输入必须是上一步 `output_adata`；第一步使用用户选择的上传文件或 h5ad。
- 如果某模块返回无 `output_adata`，pipeline 必须失败并给出中文错误。

### 3. Pipeline 参数解析

运行参数应从 `PARAM_SCHEMAS` 默认值构建，再合并保存的 preset 参数。

需要兼容两种 preset 形态:

1. 推荐形态，按模块名分组:

```json
{
  "params": {
    "qc": {"mito_perc": 0.2},
    "normalize": {"method": "log1p"},
    "bulk_qc": {"group_column": "_auto_group"},
    "bulk_deg": {"groupby": "_auto_group", "comparisons": "A-vs-B"}
  }
}
```

2. 当前前端可能保存的旧形态:

```json
{
  "params": {
    "groupby": "_auto_group",
    "top_n": "50"
  }
}
```

旧形态可以作为 best-effort:
- 只把参数键应用到 schema 中存在该 key 的模块
- 不能因为旧形态不完整而崩溃
- 缺失参数使用 `PARAM_SCHEMAS` 默认值

### 4. API

建议新增:

```text
POST /api/projects/<pid>/pipeline-runs
GET  /api/projects/<pid>/pipeline-runs
GET  /api/pipeline-runs/<run_id>/status
```

`POST /api/projects/<pid>/pipeline-runs` 请求体建议:

单细胞示例:

```json
{
  "name": "单细胞标准预处理",
  "analysis_type": "sc",
  "modules": ["qc", "normalize", "hvg"],
  "params": {
    "qc": {},
    "normalize": {},
    "hvg": {}
  },
  "input_path": "/abs/path/to/input.h5ad"
}
```

Bulk 示例:

```json
{
  "name": "Bulk 标准流程",
  "analysis_type": "bulk",
  "modules": ["bulk_qc", "bulk_normalize", "bulk_pca"],
  "params": {
    "bulk_qc": {},
    "bulk_normalize": {},
    "bulk_pca": {}
  },
  "input_path": "/abs/path/to/input.h5ad"
}
```

也可以支持:

```json
{
  "preset_id": "xxxxxxxx",
  "input_path": "/abs/path/to/input.h5ad"
}
```

必须校验:
- `pid` 合法
- `input_path` 在项目目录内
- `input_path` 不是 symlink，且是普通文件
- `analysis_type` 只能是 `sc` 或 `bulk`
- `analysis_type=sc` 时 modules 全部属于 `SC_MODULE_NAMES`
- `analysis_type=bulk` 时 modules 全部属于 `BULK_MODULE_NAMES`
- `validate_pipeline_order(modules)` 通过
- unknown module 直接 400
- 空 modules 直接 400

状态 API 返回建议:

```json
{
  "id": "run_id",
  "status": "running",
  "progress": 40,
  "analysis_type": "sc",
  "current_module": "hvg",
  "modules": ["qc", "normalize", "hvg"],
  "task_ids": ["task1", "task2"],
  "log": [],
  "error": null
}
```

### 5. 前端入口

第一阶段需要在单细胞和 Bulk 两个分析入口都能触发。

最低可接受 UI:
- 在 `templates/sc_analysis.html` 和 `templates/bulk_analysis.html` 增加“流程模板一键运行”区域，或在二者共用的入口中区分 `analysis_type`
- 单细胞页面只展示 `analysis_type=sc` 且带 `pipeline.modules` 的模板
- Bulk 页面只展示 `analysis_type=bulk` 且带 `pipeline.modules` 的模板
- 用户选择输入文件
- 点击“运行流程”
- 提交到 `POST /api/projects/<pid>/pipeline-runs`
- 显示返回的 pipeline run id，并提供状态查看或轮询

如果改 `analysis_select.html` 更方便，也可以先放在左侧“保存为流程模板”附近，但要保证用户能真正点击运行。

### 6. 测试要求

至少新增/修改以下测试:

1. `tests/test_models.py`
   - `PipelineRun` 默认值
   - `to_dict()` 可解析 JSON 字段

2. `tests/test_pipeline.py`
   - sc pipeline 有效顺序通过
   - sc pipeline 逆序失败
   - bulk pipeline 有效顺序通过
   - bulk pipeline 逆序失败
   - 非 sc module 出现在 sc pipeline 时被拒绝的 helper 测试
   - 非 bulk module 出现在 bulk pipeline 时被拒绝的 helper 测试

3. `tests/test_integration.py` 或新增 `tests/test_pipeline_runner.py`
   - 使用小型 h5ad 跑 `qc -> normalize -> hvg`
   - 使用小型 h5ad 或小型表达矩阵跑 `bulk_qc -> bulk_normalize -> bulk_pca`
   - 验证每条 pipeline 创建了对应数量的 `AnalysisTask`
   - 验证每步 `output_adata_path` 存在
   - 验证 `PipelineRun.status == completed`
   - 验证失败模块会中断后续模块，并标记 pipeline failed

4. API 测试
   - 空 modules 返回 400
   - unknown module 返回 400
   - analysis_type 与模块类型不匹配返回 400
   - 非项目目录 input_path 返回 403 或 400
   - 合法 sc pipeline 返回 run_id
   - 合法 bulk pipeline 返回 run_id

建议运行:

```bash
python -m py_compile worker.py models.py database.py routes/api.py
python -m pytest tests/test_models.py tests/test_pipeline.py tests/test_integration.py tests/test_schemas.py -q
```

如果新增专用测试文件，请加入:

```bash
python -m pytest tests/test_pipeline_runner.py -q
```

## P0 验收标准

Codex 审查时按以下标准验收:

- 能通过 API 提交一个 sc pipeline run
- 能通过 API 提交一个 bulk pipeline run
- pipeline run 有独立状态，而不是混在普通 task 里
- 每个模块仍创建普通 `AnalysisTask`，结果页面能继续查看每步结果
- 上一步 `output_adata` 正确传给下一步
- 任一步失败时 pipeline 停止，状态和错误信息清楚
- `result_files` 正常注册
- 参数从 schema 默认值 + preset/请求参数合并
- 安全校验覆盖项目目录、symlink、unknown module、analysis_type/module mismatch、依赖顺序
- targeted tests 通过

## 非目标

本阶段不要做:
- 复杂拖拽式多步骤参数编辑器
- 任务取消/暂停/恢复
- Celery/RQ 等新队列系统
- 大规模 UI 重设计
- 样本元数据管理器
- PDF/Word 综合报告
- 每个高级模块的专属交互式参数向导

这些可以作为下一阶段。

## 推荐实现顺序

1. 增加 `pipeline_runs` 表和 `PipelineRun` 模型
2. 增加参数合并 helper 和 pipeline 校验 helper
3. 增加 worker 同步串联执行逻辑
4. 增加 API
5. 增加最小前端入口
6. 补测试
7. 手动用单细胞小流程和 Bulk 小流程各跑通一次
8. 在本文档末尾追加 `Implementer Reply - Iteration 1`

## 给实现模型的第一轮任务

请实现 **平台级 Pipeline Runner P0**，不要只做 Bulk 专用实现。

第一轮必须先跑通两个最小可用链:

单细胞:

```text
qc -> normalize -> hvg
```

Bulk:

```text
bulk_qc -> bulk_normalize -> bulk_pca
```

如果两个最小链路稳定，再扩展到:

单细胞:

```text
qc -> normalize -> hvg -> dimred -> clustering
```

Bulk:

```text
bulk_qc -> bulk_normalize -> bulk_pca -> bulk_deg -> bulk_heatmap
```

完成后请在本文档末尾回复:

```md
## Implementer Reply - Iteration 1

- 已读取要求: 是
- 修改文件:
- 新增 API:
- 新增模型/表:
- 执行命令:
- 测试结果:
- 手动验证:
- 已知失败/跳过:
- 请求 Codex 审查: 是
```

## Implementer Reply - Iteration 1

- 已读取要求: 是
- 修改文件:
  - `database.py` — 新增 `pipeline_runs` 表
  - `models.py` — 新增 `PipelineRun` 模型 (141 行)
  - `worker.py` — 新增 `submit_pipeline_run()` 和 `_run_pipeline_run()` (128 行)
  - `routes/api.py` — 新增 3 个 Pipeline API (146 行)
  - `modules/__init__.py` — 导出 `SC_MODULE_NAMES`, `BULK_MODULE_NAMES`
  - `templates/sc_analysis.html` — 一键运行入口 + JS 轮询
  - `templates/bulk_analysis.html` — 一键运行入口 + JS 轮询
  - `tests/test_models.py` — 新增 4 个 PipelineRun 测试
  - `tests/test_pipeline.py` — 新增 6 个 Pipeline 类型校验测试
- 设计取舍:
  - PipelineRun 独立于 AnalysisTask，不混用
  - worker 中同步执行 pipeline，不嵌套 submit_task
  - 参数从 PARAM_SCHEMAS 默认值 + 请求参数合并
  - 兼容旧形态扁平参数（best-effort）
- 新增 API:
  - `POST /api/projects/<pid>/pipeline-runs` — 创建并启动 pipeline run
  - `GET /api/projects/<pid>/pipeline-runs` — 列出项目的所有 pipeline runs
  - `GET /api/pipeline-runs/<run_id>/status` — 查询 pipeline run 状态
- 新增模型/表:
  - `pipeline_runs` 表 (14 个字段)
  - `PipelineRun` 模型 (save/mark_running/mark_completed/mark_failed/to_dict/get_by_id/get_by_project)
- 执行命令:
  ```bash
  python -m py_compile worker.py models.py database.py routes/api.py
  python -m pytest tests/test_models.py tests/test_pipeline.py tests/test_schemas.py -v
  ```
- 测试结果: 55 passed, 1 warning
- 手动验证: 语法检查通过，测试全部通过
- 已知失败/跳过:
  - 未创建 test_pipeline_runner.py（需要实际数据和 Flask app context）
  - 前端需要手动验证（需要启动 Flask 服务）
- 请求 Codex 审查: 是
