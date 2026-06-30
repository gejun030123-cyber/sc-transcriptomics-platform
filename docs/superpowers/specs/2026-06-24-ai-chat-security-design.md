# AI 对话模块安全加固与确认执行机制

## 背景

AI 对话模块（`routes/chat.py` + `modules/ai_adapter.py` + `modules/ai_tools.py`）存在以下问题：
- Prompt injection 风险（文本解析 `tool` 块可被注入）
- `_run_analysis` 参数无校验，可注入任意参数
- 聊天历史内存无限增长、多进程不共享
- AI API 调用无超时
- 工具调用自动执行无用户确认
- 接口无认证

## 设计决策

| 决策 | 选择 | 理由 |
|---|---|---|
| 分析任务执行模式 | 前端确认后执行 | 安全 + 用户可控 |
| 聊天历史存储 | 内存 + LRU 大小限制 | 简单，无外部依赖 |
| 认证方式 | API Token（环境变量） | 轻量，无用户系统 |
| 工具调用解析 | 仅原生 tool_use | 彻底消除 prompt injection |

## 架构变更

### 1. 移除 prompt injection 路径

**文件**: `modules/ai_adapter.py`

- 删除 `_parse_tool_calls_from_text()` 函数
- 删除 `_strip_tool_blocks()` 函数
- 删除 `SYSTEM_PROMPT` 中 `## 工具调用方式` 段落（````tool` 块格式说明）
- `_chat_anthropic()` 中移除文本解析降级逻辑（lines 251-267）
- Anthropic 路径仅处理 `response.content` 中的 `tool_use` 块

### 2. 确认执行机制

**当前流程**:
```
用户消息 → /api/chat → AI 返回 tool_use → 后端立即执行 → 返回结果
```

**新流程**:
```
用户消息 → /api/chat → AI 返回 tool_use → 后端不执行 → 返回 proposed_tools
                                                          ↓
                                               前端显示确认卡片
                                                          ↓
                                               用户点确认 → POST /api/chat/approve
                                                          ↓
                                               后端执行 → 返回结果
```

#### 2.1 `/api/chat` 返回格式变更

```json
{
  "reply": "我将用 DESeq2 方法对 hmc3 vs ctrl 进行差异分析...",
  "proposed_tools": [
    {
      "name": "run_analysis",
      "args": {"module_name": "bulk_deg", "params": {"method": "deseq2"}},
      "description": "执行 bulk_deg 差异表达分析"
    }
  ],
  "tool_calls": []
}
```

- `proposed_tools`: AI 提议但未执行的工具调用列表（仅 `run_analysis` 需确认）
- `tool_calls`: 已自动执行的工具调用（`get_project_status`、`get_task_results`、`list_modules` 仍自动执行）

#### 2.2 工具分类

| 工具 | 行为 | 理由 |
|---|---|---|
| `run_analysis` | 需确认 | 触发计算资源，影响数据 |
| `get_project_status` | 自动执行 | 只读查询，无副作用 |
| `get_task_results` | 自动执行 | 只读查询，无副作用 |
| `list_modules` | 自动执行 | 只读查询，无副作用 |

#### 2.3 新增 `POST /api/chat/approve` 端点

```python
@chat_bp.route('/api/chat/approve', methods=['POST'])
@require_ai_token
def approve_tool():
    data = request.get_json()
    tool_name = data.get('tool_name')
    args = data.get('args', {})
    project_id = data.get('project_id')
    
    if tool_name != 'run_analysis':
        return jsonify({"error": "只有 run_analysis 需要确认"}), 400
    
    # 参数校验
    validation_error = _validate_analysis_params(args)
    if validation_error:
        return jsonify({"error": validation_error}), 400
    
    result = execute_tool(tool_name, args, project_id)
    return jsonify(result)
```

#### 2.4 参数校验

`_run_analysis` 新增 `_validate_analysis_params(args)` 校验：
- `module_name` 必须在 `MODULE_REGISTRY` 中
- `params` 的键必须在对应模块的 `PARAM_SCHEMAS` 中
- `params` 的值类型必须匹配 schema 定义
- 未知键直接移除（不报错，避免 AI 幻觉参数干扰）

### 3. API Token 认证

**文件**: `config.py`, `routes/chat.py`

```python
# config.py 新增
AI_API_TOKEN = os.environ.get('AI_API_TOKEN', '')

# routes/chat.py 装饰器
def require_ai_token(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if Config.AI_API_TOKEN:
            auth = request.headers.get('Authorization', '')
            if not auth.startswith('Bearer ') or auth[7:] != Config.AI_API_TOKEN:
                return jsonify({"error": "认证失败"}), 401
        return f(*args, **kwargs)
    return decorated
```

- `AI_API_TOKEN` 为空时跳过认证（向后兼容）
- 所有 `/api/chat*` 路由加 `@require_ai_token`

### 4. 聊天历史 LRU 限制

**文件**: `routes/chat.py`

```python
MAX_HISTORY_PER_PROJECT = 50

class ChatHistoryStore:
    def __init__(self, max_per_project=MAX_HISTORY_PER_PROJECT):
        self._store = {}
        self._max = max_per_project

    def get(self, project_id):
        return self._store.get(project_id, [])

    def set(self, project_id, messages):
        if len(messages) > self._max:
            messages = messages[-self._max:]
        self._store[project_id] = messages

    def clear(self, project_id):
        self._store.pop(project_id, None)

_chat_histories = ChatHistoryStore()
```

- 每个项目最多保留 50 条消息（user + assistant）
- 超出时丢弃最早的对话轮次
- `tool_result` 类型的消息不单独计数（作为 assistant 消息的一部分）

### 5. API 超时

**文件**: `modules/ai_adapter.py`

```python
# Anthropic
client = anthropic.Anthropic(
    base_url=Config.AI_API_URL,
    api_key=Config.AI_API_KEY,
    timeout=60.0,  # 新增
)

# OpenAI
client = OpenAI(
    base_url=Config.AI_API_URL,
    api_key=Config.AI_API_KEY,
    timeout=60.0,  # 新增
)
```

## 不变的部分

- `_chat_openai` 的工具循环逻辑保持不变（已原生支持 tool_calls）
- `TOOLS_ANTHROPIC` / `TOOLS_OPENAI` 工具定义不变
- `SYSTEM_PROMPT` 的能力描述和回复规则不变
- 前端聊天 UI 组件结构不变（新增确认卡片样式）

## 文件变更清单

| 文件 | 变更类型 | 内容 |
|---|---|---|
| `config.py` | 修改 | 新增 `AI_API_TOKEN` |
| `routes/chat.py` | 修改 | Token 认证装饰器、聊天历史 LRU、approve 端点、proposed_tools 返回 |
| `modules/ai_adapter.py` | 修改 | 移除文本解析、超时、工具分类逻辑 |
| `modules/ai_tools.py` | 修改 | 新增 `_validate_analysis_params` |
| `templates/base.html` 或相关 JS | 修改 | 前端确认卡片（AI 消息中的 proposed_tools 渲染） |
