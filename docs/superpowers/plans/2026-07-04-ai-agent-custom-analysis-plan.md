# AI 目标驱动单细胞/Bulk 分析 Agent 任务计划

日期：2026-07-04

## 1. 目标

构建一个目标驱动型生信 AI agent，使用户可以通过自然语言提出分析目标，AI 能围绕目标检查当前结果、拆解问题、调整参数、创建候选分支、评估候选结果，并在用户确认后采纳结果。

第一阶段重点场景：

> 用户对单细胞分群不满意，指定必须找到最接近某种细胞类型的分群。AI 需要基于 marker、当前聚类、UMAP、batch/sample 分布和参数搜索结果，推荐最可信的候选分群。

本计划不要求 AI 直接自由编写和执行任意代码。AI 只能调用注册工具，所有高成本或写入型动作必须可审计、可确认、可回滚。

## 2. 当前项目基础

已存在能力：

- Flask Web 应用和项目管理。
- `analysis_tasks` 单模块任务模型。
- `pipeline_runs` 多模块流程模型。
- `worker.py` 后台执行单任务和 pipeline run。
- `modules/schemas.py` 中已有完整参数 schema。
- `routes/chat.py`、`modules/ai_adapter.py`、`modules/ai_tools.py` 中已有 AI chat 和工具调用基础。
- 单细胞模块：`qc`、`normalize`、`hvg`、`dimred`、`batch_correct`、`clustering`、`annotation`、`deg` 等。
- Bulk 模块：`bulk_qc`、`bulk_normalize`、`bulk_pca`、`bulk_deg`、`bulk_heatmap`、`bulk_enrichment` 等。

主要缺口：

- 当前 AI 更像“流程触发器”，不是目标驱动 agent。
- 缺少分析状态深度检查工具。
- 缺少候选分支和参数搜索机制。
- 缺少 marker/cluster 证据评分。
- 缺少 agent 步骤审计和用户采纳机制。
- 缺少前端候选实验面板。

## 3. 产品原则

1. AI 的职责是辅助实验设计和结果评估，不是无约束自动跑流程。
2. 所有写操作、高成本运行、覆盖主线结果操作必须用户确认。
3. 所有 agent 动作必须落库，能追踪“为什么这么做”。
4. 候选实验默认在 branch 中运行，不覆盖主项目最新结果。
5. 参数必须经过 `PARAM_SCHEMAS` 校验。
6. 单次参数搜索必须有候选数量上限。
7. AI 可以推荐“最接近目标的 cluster”，但不能承诺数据中一定存在该细胞类型。
8. 失败分支必须保留错误信息，不能静默跳过。

## 4. 总体架构

建议新增层次：

```text
用户自然语言
  -> AI chat / LLM tool calling
  -> agent_orchestrator
  -> agent tools / operation registry
  -> branch runner / existing worker
  -> evaluator
  -> candidate report
  -> user accept / refine / reject
```

建议新增模块：

- `modules/agent_orchestrator.py`
- `modules/agent_tools.py`
- `modules/agent_models.py` 或扩展 `models.py`
- `modules/cell_markers.py`
- `modules/evaluators/sc_cluster.py`
- `modules/evaluators/bulk_deg.py`

建议新增前端区域：

- 项目详情页 AI 聊天旁的“Agent 实验面板”。
- 候选 branch 卡片。
- 参数 diff、评分表、证据摘要、采纳按钮。

## 5. 阶段 0：修复现有 AI Chat 基础问题

### 任务 0.1：修复 `proposed_tools` 返回缺失

涉及文件：

- `routes/chat.py`
- `modules/ai_adapter.py`
- `templates/project_detail.html`

要求：

- `/api/chat` 必须返回 `proposed_tools`。
- 前端已有 `data.proposed_tools` 渲染逻辑，后端必须传出该字段。
- `run_analysis` 仍然必须用户确认后执行。
- `get_project_status`、`get_task_results`、`list_modules` 仍可自动执行。

验收：

- 用户说“帮我运行 bulk_deg”，页面出现确认卡片。
- 用户点击确认后任务提交。
- 未点击确认时不产生新的 `AnalysisTask`。
- `/api/chat` 响应至少包含：

```json
{
  "reply": "...",
  "tool_calls": [],
  "proposed_tools": [
    {
      "name": "run_analysis",
      "args": {}
    }
  ]
}
```

审批 gate：

- 确认 `run_analysis` 不会绕过确认机制。
- 确认前端不会把所有工具都当成可执行写操作。

## 6. 阶段 1：新增 Agent 数据模型

### 任务 1.1：新增 `agent_sessions`

用途：记录一次用户与 agent 的目标分析会话。

建议 schema：

```sql
CREATE TABLE IF NOT EXISTS agent_sessions (
    id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    status TEXT DEFAULT 'active',
    goal_summary TEXT DEFAULT '',
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,
    updated_at TEXT DEFAULT CURRENT_TIMESTAMP
);
```

### 任务 1.2：新增 `agent_goals`

用途：记录结构化目标，例如“寻找 microglia-like cluster”。

建议 schema：

```sql
CREATE TABLE IF NOT EXISTS agent_goals (
    id TEXT PRIMARY KEY,
    session_id TEXT NOT NULL REFERENCES agent_sessions(id) ON DELETE CASCADE,
    project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    goal_type TEXT NOT NULL,
    target_entity TEXT DEFAULT '',
    goal_json TEXT NOT NULL,
    constraints_json TEXT DEFAULT '{}',
    status TEXT DEFAULT 'draft',
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,
    updated_at TEXT DEFAULT CURRENT_TIMESTAMP
);
```

示例 `goal_json`：

```json
{
  "goal_type": "target_cluster_refinement",
  "target_cell_type": "microglia",
  "positive_markers": ["P2RY12", "TMEM119", "CX3CR1", "AIF1"],
  "negative_markers": ["LYZ", "S100A8", "FCGR3A"],
  "user_requirement": "找到最接近小胶质细胞的分群"
}
```

### 任务 1.3：新增 `agent_steps`

用途：记录 agent 每一步观察、假设、动作和结果。

建议 schema：

```sql
CREATE TABLE IF NOT EXISTS agent_steps (
    id TEXT PRIMARY KEY,
    session_id TEXT NOT NULL REFERENCES agent_sessions(id) ON DELETE CASCADE,
    goal_id TEXT REFERENCES agent_goals(id) ON DELETE SET NULL,
    step_index INTEGER DEFAULT 0,
    step_type TEXT DEFAULT '',
    thought_summary TEXT DEFAULT '',
    action_name TEXT DEFAULT '',
    action_args_json TEXT DEFAULT '{}',
    result_json TEXT DEFAULT '{}',
    status TEXT DEFAULT 'completed',
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);
```

注意：`thought_summary` 只能保存面向用户和审计的简短摘要，不保存 LLM 隐式推理链。

### 任务 1.4：新增 `analysis_branches`

用途：保存候选分支，避免候选实验覆盖主线结果。

建议 schema：

```sql
CREATE TABLE IF NOT EXISTS analysis_branches (
    id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    session_id TEXT REFERENCES agent_sessions(id) ON DELETE SET NULL,
    goal_id TEXT REFERENCES agent_goals(id) ON DELETE SET NULL,
    parent_task_id TEXT REFERENCES analysis_tasks(id) ON DELETE SET NULL,
    parent_adata_path TEXT NOT NULL,
    branch_name TEXT DEFAULT '',
    purpose TEXT DEFAULT '',
    params_json TEXT DEFAULT '{}',
    output_adata_path TEXT,
    status TEXT DEFAULT 'pending',
    accepted INTEGER DEFAULT 0,
    error_traceback TEXT,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,
    finished_at TEXT
);
```

### 任务 1.5：新增 `candidate_scores`

用途：保存候选结果评分。

建议 schema：

```sql
CREATE TABLE IF NOT EXISTS candidate_scores (
    id TEXT PRIMARY KEY,
    branch_id TEXT NOT NULL REFERENCES analysis_branches(id) ON DELETE CASCADE,
    project_id TEXT NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    evaluator_name TEXT NOT NULL,
    target_label TEXT DEFAULT '',
    score_json TEXT NOT NULL,
    total_score REAL,
    recommendation TEXT DEFAULT '',
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);
```

验收：

- 数据库初始化后自动创建所有表。
- 删除项目时相关 agent 记录级联删除。
- 不破坏现有 `projects`、`analysis_tasks`、`result_files`、`pipeline_runs`。
- 新增模型类或 helper，能完成 create/get/list/update。

审批 gate：

- 表结构必须能审计 agent 行为。
- branch 必须能追踪 parent h5ad 和参数。
- 不允许只把 agent 状态塞进聊天历史。

## 7. 阶段 2：新增只读分析状态工具

### 任务 2.1：实现 `inspect_analysis_state`

建议位置：

- `modules/agent_tools.py`
- 并在 `modules/ai_tools.py::execute_tool` 中注册。

输入：

```json
{
  "project_id": "..."
}
```

输出：

```json
{
  "project_id": "...",
  "latest_adata": "...",
  "completed_tasks": [],
  "available_cluster_keys": ["leiden", "leiden_0.8"],
  "available_embeddings": ["X_umap", "X_pca"],
  "available_annotations": ["celltype"],
  "obs_columns": [],
  "n_obs": 12345,
  "n_vars": 20000,
  "warnings": []
}
```

要求：

- 优先读取项目最新完成任务的 `output_adata_path`。
- 如果没有输出，则检查上传目录中的 h5ad。
- 返回任务历史、可用结果、可用 obs 字段。
- 不执行任何写操作。

### 任务 2.2：实现 `inspect_adata`

输入：

```json
{
  "adata_path": "..."
}
```

输出：

```json
{
  "n_obs": 10000,
  "n_vars": 20000,
  "obs_columns": [],
  "var_columns": [],
  "layers": [],
  "obsm": [],
  "cluster_keys": [],
  "embedding_keys": [],
  "sample_preview": []
}
```

要求：

- 路径必须位于项目目录内。
- 禁止符号链接。
- 大文件优先使用 backed 方式读取。
- 出错时返回结构化 error。

### 任务 2.3：实现 `get_cluster_summary`

输入：

```json
{
  "adata_path": "...",
  "cluster_key": "leiden"
}
```

输出：

```json
{
  "cluster_key": "leiden",
  "clusters": [
    {
      "cluster": "0",
      "n_cells": 1200,
      "pct_cells": 0.12,
      "qc_means": {},
      "batch_distribution": {},
      "sample_distribution": {},
      "umap_center": [1.2, -0.4]
    }
  ]
}
```

验收：

- AI 可以在不运行分析的情况下知道当前有哪些 cluster。
- 能识别 `leiden`、`louvain`、`celltype` 等常见列。
- 对缺少 cluster key 的文件返回明确提示。

审批 gate：

- 所有路径安全校验必须复用或统一实现。
- 只读工具不得创建任务、分支或文件。

## 8. 阶段 3：新增细胞类型 marker 打分能力

### 任务 3.1：新增 marker 库

建议文件：

- `modules/cell_markers.py`

首版支持：

| 细胞类型 | positive markers | negative markers |
| --- | --- | --- |
| T cell | CD3D, CD3E, TRAC | MS4A1, LST1 |
| B cell | MS4A1, CD79A, CD79B | CD3D, LST1 |
| NK | NKG7, GNLY, KLRD1 | CD3D, MS4A1 |
| Monocyte | LST1, S100A8, S100A9, FCN1 | MS4A1, CD3D |
| Macrophage | C1QA, C1QB, CD68, APOE | MS4A1, CD3D |
| Dendritic cell | FCER1A, CLEC10A, LILRA4 | CD3D, MS4A1 |
| Epithelial | EPCAM, KRT8, KRT18 | PTPRC |
| Endothelial | PECAM1, VWF, KDR | PTPRC |
| Fibroblast | COL1A1, COL1A2, DCN | PTPRC, EPCAM |
| Microglia | P2RY12, TMEM119, CX3CR1, AIF1 | S100A8, FCGR3A |
| Astrocyte | GFAP, AQP4, ALDH1L1 | PTPRC |
| Oligodendrocyte | MBP, MOG, PLP1 | PTPRC |

要求：

- 用户自定义 marker 优先。
- 内置 marker 仅作为建议。
- marker 不存在时返回 warning。

### 任务 3.2：实现 `score_cell_type_signature`

输入：

```json
{
  "adata_path": "...",
  "cluster_key": "leiden",
  "target_cell_type": "microglia",
  "positive_markers": ["P2RY12", "TMEM119", "CX3CR1", "AIF1"],
  "negative_markers": ["LYZ", "S100A8", "FCGR3A"]
}
```

输出：

```json
{
  "target_cell_type": "microglia",
  "cluster_key": "leiden",
  "best_cluster": "7",
  "confidence": "medium",
  "warnings": [],
  "cluster_scores": [
    {
      "cluster": "7",
      "n_cells": 253,
      "positive_score": 0.82,
      "negative_score": 0.12,
      "specificity_score": 0.76,
      "size_score": 0.95,
      "batch_penalty": 0.05,
      "total_score": 0.81,
      "detected_positive_markers": ["P2RY12", "TMEM119"],
      "missing_positive_markers": []
    }
  ]
}
```

评分建议：

- `positive_score`：positive marker 在 cluster 内的平均表达/检出率。
- `specificity_score`：cluster 内外表达差异。
- `negative_score`：negative marker 表达惩罚。
- `size_score`：cluster 细胞数是否合理。
- `batch_penalty`：是否被单一 batch/sample 主导。
- `total_score`：综合评分，范围 0 到 1。

验收：

- 给定 marker 后返回 cluster 排名。
- marker 全部缺失时不崩溃，返回 warning 和低 confidence。
- 能处理 gene symbol 大小写差异，至少要有明确策略。

审批 gate：

- 评分逻辑必须确定性可测试。
- 不能只让 LLM “主观判断”最佳 cluster。

## 9. 阶段 4：新增分析分支能力

### 任务 4.1：创建 branch API

新增接口：

```http
POST /api/projects/<pid>/branches
```

请求：

```json
{
  "parent_adata_path": "...",
  "branch_name": "microglia_search_res_1.2",
  "purpose": "寻找 microglia-like cluster"
}
```

响应：

```json
{
  "id": "...",
  "message": "候选分支已创建"
}
```

要求：

- `parent_adata_path` 必须在项目目录内。
- 禁止符号链接。
- branch 初始状态为 `pending`。

### 任务 4.2：在 branch 上运行模块

新增接口：

```http
POST /api/branches/<branch_id>/run
```

请求：

```json
{
  "modules": ["clustering"],
  "params": {
    "clustering": {
      "resolutions": "1.2",
      "n_neighbors": 15
    }
  }
}
```

要求：

- 复用现有 `MODULE_REGISTRY`。
- 复用 `validate_pipeline_order`。
- 参数必须通过 `PARAM_SCHEMAS` 校验。
- branch 输出写到项目 branch 专属目录。
- branch 失败不影响主线 `AnalysisTask` 最新结果。

### 任务 4.3：采纳 branch

新增接口：

```http
POST /api/branches/<branch_id>/accept
```

请求：

```json
{
  "confirm": true
}
```

要求：

- 只有用户显式确认才可采纳。
- 采纳后 `accepted=1`。
- 不删除其他 branch。
- 可选：在项目 metadata 中记录当前 accepted branch。

验收：

- 同一个 parent h5ad 能创建多个候选。
- 候选运行失败不影响项目主线。
- 采纳 branch 前项目主线结果不变。

审批 gate：

- branch 输出路径必须隔离。
- accept 操作必须是显式用户动作。

## 10. 阶段 5：实现参数搜索工具

### 任务 5.1：实现 `propose_parameter_sweep`

输入：

```json
{
  "goal_type": "target_cluster_refinement",
  "target_cell_type": "microglia",
  "current_state": {},
  "max_candidates": 6
}
```

输出：

```json
{
  "base_checkpoint": "...",
  "candidates": [
    {
      "name": "res_0.8_neighbors_15",
      "modules": ["clustering"],
      "params": {
        "clustering": {
          "resolutions": "0.8",
          "n_neighbors": 15
        }
      }
    }
  ]
}
```

首版参数空间：

- `clustering.resolutions`: `0.6`, `0.8`, `1.0`, `1.2`, `1.5`
- `clustering.n_neighbors`: `10`, `15`, `30`
- `dimred.n_comps`: `30`, `50`, `80`
- `hvg.n_top_genes`: `2000`, `3000`, `5000`
- 有 batch 且 batch 混杂明显时，可加入 `batch_correct.method=harmony`

候选生成策略：

- 如果当前已有 `X_pca` 和 `X_umap`，优先只 sweep clustering。
- 如果当前分群结构明显不稳定，再加入 dimred 参数。
- 如果目标 marker 缺失或 HVG 可能过滤目标基因，再考虑 hvg 参数。
- 默认最多 6 个候选。

### 任务 5.2：实现 `run_parameter_sweep`

输入：

```json
{
  "goal_id": "...",
  "base_checkpoint": "...",
  "candidates": [],
  "evaluator": "sc_cluster_signature"
}
```

要求：

- 逐个运行候选，避免并发打满机器。
- 每个候选创建 branch。
- 每个候选完成后自动调用 evaluator。
- 评分写入 `candidate_scores`。
- 每一步写入 `agent_steps`。

验收：

- sweep 不超过 `max_candidates`。
- 每个候选都有 params、branch、score。
- 某个候选失败不会中止全部 sweep，除非基础文件不可用。

审批 gate：

- 参数空间必须保守，不能无限组合。
- 高成本 sweep 必须用户确认。

## 11. 阶段 6：实现目标细胞分群优化 Agent

### 任务 6.1：实现 `agent_orchestrator`

建议文件：

- `modules/agent_orchestrator.py`

核心流程：

```text
1. parse_goal
2. inspect_analysis_state
3. inspect current clusters
4. score current result
5. if current result good enough:
       summarize evidence and ask user to accept / refine
   else:
       propose parameter sweep
6. ask user confirmation for sweep
7. run candidate branches
8. evaluate candidates
9. recommend best branch
10. wait for user accept / refine / stop
```

### 任务 6.2：实现 `start_goal_agent`

输入：

```json
{
  "goal_type": "target_cluster_refinement",
  "target_cell_type": "microglia",
  "positive_markers": [],
  "negative_markers": [],
  "max_candidate_runs": 6
}
```

输出：

```json
{
  "session_id": "...",
  "goal_id": "...",
  "status": "needs_confirmation",
  "summary": "...",
  "proposed_actions": []
}
```

### 任务 6.3：实现 `continue_goal_agent`

支持用户继续指令：

- “继续细分这个 cluster”
- “这个不满意，试试更高 resolution”
- “不要用 Harmony”
- “必须包含 P2RY12 和 TMEM119”
- “采用候选 B”
- “停止这个目标”

要求：

- 每次继续都加载 `agent_sessions`、`agent_goals`、`agent_steps`、候选 branch 和评分。
- 新指令必须更新 goal constraints 或创建新 step。
- 不能丢失用户约束。

验收：

- AI 能围绕同一个目标持续迭代。
- 每一步都有审计记录。
- 用户可以中途停止或采纳候选。

审批 gate：

- agent 不得绕过工具层直接执行代码。
- agent 输出必须包含下一步动作和是否需要确认。

## 12. 阶段 7：前端 Agent 实验面板

### 任务 7.1：扩展项目详情页

涉及文件：

- `templates/project_detail.html`
- 可选新增静态 JS/CSS 文件。

新增区域：

- 当前目标。
- 当前检查点。
- 候选分支列表。
- 参数 diff。
- 评分表。
- 推荐候选。
- 采纳/继续细分/丢弃按钮。

### 任务 7.2：候选结果卡片

卡片字段：

- branch 名称。
- 参数摘要。
- 最佳 cluster。
- target score。
- positive marker 命中。
- negative marker 惩罚。
- batch/sample warning。
- 查看结果。
- 采用此候选。
- 继续细分。

示例 UI 文案：

```text
候选 B：resolution=1.2, n_neighbors=15
最佳 cluster：7
Microglia score：0.81
证据：P2RY12/TMEM119/CX3CR1 表达较高，S100A8 较低
风险：cluster 7 中 batch_2 占比偏高
```

验收：

- 用户能看懂 AI 推荐理由。
- 用户能明确采纳某个候选。
- 用户不采纳时，主结果不变。

审批 gate：

- UI 不能只显示“AI 推荐”，必须显示证据和参数。
- 采纳按钮必须触发显式确认 API。

## 13. 阶段 8：安全规则

必须实现：

- 所有 agent 工具使用统一路径校验。
- 所有写操作必须区分 `auto`、`needs_confirmation`、`confirmed`。
- 所有参数都通过 schema 清洗。
- 所有 branch 默认不可覆盖主线。
- 所有高成本任务有数量上限。
- 所有 API 返回结构化 error。

禁止：

- 禁止 AI 执行任意 Python 代码。
- 禁止 AI 接收任意文件系统路径并直接读取。
- 禁止未确认自动采纳 branch。
- 禁止未确认执行参数 sweep。
- 禁止把失败候选静默隐藏。

## 14. 阶段 9：测试计划

建议新增测试文件：

- `tests/test_agent_models.py`
- `tests/test_agent_tools.py`
- `tests/test_signature_scoring.py`
- `tests/test_branch_runs.py`
- `tests/test_parameter_sweep.py`
- `tests/test_agent_orchestrator.py`

核心测试用例：

1. 非项目路径被拒绝。
2. symlink 输入被拒绝。
3. marker 不存在时返回 warning。
4. 当前 clustering 已满足目标时不启动 sweep。
5. 当前 clustering 不满足目标时生成候选参数。
6. sweep 不超过最大候选数。
7. branch 失败不影响主线。
8. accept branch 需要显式确认。
9. 参数中未知 key 被过滤。
10. 用户自定义 marker 优先于内置 marker。
11. branch 输出路径位于项目目录内。
12. agent step 能完整记录 action 和 result。

最低验收命令：

```bash
pytest tests/test_agent_models.py tests/test_agent_tools.py tests/test_signature_scoring.py
pytest tests/test_branch_runs.py tests/test_parameter_sweep.py tests/test_agent_orchestrator.py
pytest tests/test_pipeline.py tests/test_schemas.py
```

## 15. 第一轮交付范围

第一轮只交付一个闭环：

> 单细胞目标细胞分群优化：用户指定细胞类型或 marker，AI 检查当前分群，尝试少量 clustering 参数，评分候选，推荐最佳分支，用户确认后采纳。

第一轮不做：

- 不做全自动长链路科研报告。
- 不做任意自然语言转代码执行。
- 不做复杂外部数据库联网查询。
- 不做无限参数搜索。
- 不做自动覆盖主线结果。

第一轮必须包含：

- `proposed_tools` 修复。
- agent 数据表。
- 只读分析状态工具。
- marker 打分工具。
- branch 创建和运行。
- 最多 6 个候选的参数 sweep。
- 候选评分和推荐。
- 用户采纳 branch。
- 基础前端候选卡片。
- 关键测试。

## 16. 后续扩展

第二轮可扩展：

- 指定 cluster 的二次细分。
- annotation 修正和细胞类型重命名建议。
- batch correction 策略自动比较。
- QC 阈值调整建议。
- marker dotplot/feature plot 自动生成。

第三轮可扩展到 Bulk RNA：

- 用户指定比较目标，AI 自动组织 DEG。
- 多组比较一致性筛选。
- 交集/并集/方向性基因集合。
- 富集分析和通路解释。
- 热图基因集自动选择。

## 17. 审批清单

每个实现 PR 或代码批次必须回答：

1. 本次新增了哪些 agent 工具？
2. 哪些工具是只读，哪些工具会写入或运行任务？
3. 写入工具是否需要用户确认？
4. 参数是否通过 schema 校验？
5. 路径是否限制在项目目录内？
6. branch 是否隔离主线结果？
7. 失败结果是否可见？
8. 是否新增测试覆盖？
9. 是否影响现有 pipeline run？
10. 用户是否能理解 AI 推荐依据？

审批结论只能是：

- `approved`：可以合并。
- `changes_requested`：必须修改后复审。
- `blocked`：设计方向不符合安全或产品目标，需要重新设计。

---

## 18. Codex 审批回复

审批日期：2026-07-04

审批对象：本任务计划文档，即目标驱动型单细胞/Bulk 分析 Agent 的设计与分阶段实施方案。

审批结论：**`approved`**

本计划可以作为后续实现工作的基准文档。方案方向符合项目目标：不是把 AI 做成固定流程触发器，而是围绕用户目标进行状态检查、参数试验、候选分支评估、证据展示和用户确认采纳。

### 通过理由

1. **产品目标清晰**

   计划明确第一轮只聚焦“单细胞目标细胞分群优化”闭环，避免一开始扩展到过大的科研自动化范围。这个收敛是合理的。

2. **Agent 行为边界明确**

   文档明确禁止任意代码执行、禁止未确认高成本任务、禁止未确认覆盖主线结果。这些规则能有效降低 AI agent 失控风险。

3. **Branch 隔离原则正确**

   候选实验默认在 branch 中运行，用户采纳前不覆盖主线结果。这是该系统能否可信使用的核心约束，必须在实现中严格保持。

4. **审计链路充分**

   `agent_sessions`、`agent_goals`、`agent_steps`、`analysis_branches`、`candidate_scores` 的建模方向合理，能追踪目标、动作、参数、结果和推荐依据。

5. **工具分层合理**

   计划区分只读工具、写入工具、高成本工具和采纳工具，并要求写操作走确认机制。这符合当前平台已有 AI 工具调用模型。

6. **测试要求具体**

   测试计划覆盖路径安全、symlink 拒绝、marker 缺失、参数上限、branch 失败隔离、显式采纳和 schema 清洗等关键风险点。

### 实现时必须保持的硬性要求

以下要求属于后续代码审批的硬门槛，不满足则应给出 `changes_requested`：

1. 未采纳 branch 不能成为项目 latest h5ad。
2. branch 内任务不能污染主线 `AnalysisTask` 最新结果查询。
3. 所有 agent 写操作必须需要用户确认。
4. 所有路径必须限制在项目目录内，并拒绝 symlink。
5. 所有参数必须通过 `PARAM_SCHEMAS` 清洗后才能执行。
6. `run_parameter_sweep` 必须有执行期硬上限，不能只在候选生成阶段限制。
7. `goal_id`、`session_id`、`branch_id` 必须校验 project 归属。
8. 用户采纳 branch 必须是显式动作，不能由 AI 自动完成。
9. 测试命令必须在干净环境中可复现。
10. 前端必须展示推荐依据，而不是只展示“AI 推荐”。

### 非阻塞建议

- 第一轮实现中可以先只允许 `clustering` sweep，等闭环稳定后再加入 `dimred`、`hvg`、`batch_correct`。
- marker scoring 首版可以保持确定性算法，但建议尽早对表达值做归一化或 rank/percentile 转换，避免 raw mean expression 导致评分饱和。
- Branch API 最好设计为 project-scoped，例如 `/api/projects/<pid>/branches/<branch_id>/...`，便于做归属校验。
- 如果没有完整用户系统，至少应让 agent 写接口复用当前 AI token 或同等保护机制。

### 审批说明

本次 `approved` 仅表示设计计划可以作为实现依据，不代表任何具体代码实现自动通过。后续实现仍必须按第 17 节审批清单逐项复审。
