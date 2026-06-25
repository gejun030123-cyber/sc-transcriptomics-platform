# modules/ai_adapter.py
"""AI 对话适配器 — 支持 Anthropic (Claude) 和 OpenAI 兼容 API"""
import json
from config import Config


# 工具定义：AI 可调用的平台功能
TOOLS_ANTHROPIC = [
    {
        "name": "run_analysis",
        "description": "执行一个分析模块。返回任务 ID 和状态。",
        "input_schema": {
            "type": "object",
            "properties": {
                "module_name": {
                    "type": "string",
                    "description": "分析模块名称，如 bulk_deg, bulk_heatmap, bulk_normalize 等"
                },
                "input_path": {
                    "type": "string",
                    "description": "输入数据文件路径（h5ad/csv）。留空则自动查找最新可用文件。"
                },
                "params": {
                    "type": "object",
                    "description": "分析参数，如 {\"method\": \"deseq2\", \"fc_threshold\": 1.5}。留空使用默认值。"
                }
            },
            "required": ["module_name"]
        }
    },
    {
        "name": "get_project_status",
        "description": "获取当前项目的完整状态：已上传文件、已完成任务、可用的中间文件。",
        "input_schema": {
            "type": "object",
            "properties": {},
            "required": []
        }
    },
    {
        "name": "get_task_results",
        "description": "获取某个分析任务的结果摘要，包括统计数据和生成的图表列表。",
        "input_schema": {
            "type": "object",
            "properties": {
                "task_id": {
                    "type": "string",
                    "description": "任务 ID"
                }
            },
            "required": ["task_id"]
        }
    },
    {
        "name": "list_modules",
        "description": "列出所有可用的分析模块及其说明和参数。",
        "input_schema": {
            "type": "object",
            "properties": {
                "pipeline_type": {
                    "type": "string",
                    "enum": ["bulk", "sc", "all"],
                    "description": "模块类型：bulk=Bulk RNA-seq, sc=单细胞, all=全部"
                }
            },
            "required": []
        }
    }
]

TOOLS_OPENAI = [
    {"type": "function", "function": {**t, "parameters": t["input_schema"]}}
    for t in TOOLS_ANTHROPIC
]

# 自动执行的只读工具（无需用户确认）
AUTO_EXEC_TOOLS = {'get_project_status', 'get_task_results', 'list_modules'}
# 需要用户确认的工具
CONFIRM_TOOLS = {'run_analysis'}


# System prompt
SYSTEM_PROMPT = """你是一个生信分析助手，帮助用户进行 RNA-seq 数据分析。你的平台是基于 Flask 的 Web 应用，支持单细胞和 Bulk RNA-seq 全流程分析。

## 你的能力
1. **自然语言触发分析**：用户说"对 hmc3 和 ctrl 做差异分析"，你调用 run_analysis(module_name="bulk_deg", params={{...}})
2. **结果解读**：用户问"哪些基因在所有药物中共同上调？"，你调用 get_task_results 查看结果并解读
3. **参数建议**：用户问"用哪个方法好？"，你根据数据情况给出建议
4. **项目概览**：用户问"现在分析到哪一步了？"，你调用 get_project_status 了解情况

## 可用模块
- bulk_qc: 数据质控
- bulk_normalize: 数据标准化 (DESeq2/TMM/CPM/VST/rlog)
- bulk_deg: 差异表达分析 (t-test/mann-whitney/deseq2/edger/limma)
- bulk_pca: PCA/UMAP 降维
- bulk_heatmap: 热图可视化
- bulk_enrichment: 通路富集 (GO/KEGG/WikiPathways)
- bulk_timecourse: 时序分析
- bulk_deg_integration: 多组差异整合

## 回复规则
- 用中文回复
- 简洁明了，不啰嗦
- 执行分析时，先确认参数再执行（如"我将用 DESeq2 方法对 hmc3 vs ctrl 进行差异分析，FC 阈值 2.0，padj 阈值 0.05，确认执行吗？"）
- 分析完成后，简要总结关键结果
- 如果用户要求不明确，主动询问关键参数"""


def _is_anthropic():
    """判断是否使用 Anthropic API"""
    url = Config.AI_API_URL.lower()
    return 'anthropic' in url or 'claude' in url


def chat(messages, project_id=None):
    """
    与 LLM 对话，支持工具调用循环。自动检测 Anthropic / OpenAI 格式。
    """
    tool_calls_log = []

    if _is_anthropic():
        return _chat_anthropic(messages, project_id, tool_calls_log)
    else:
        return _chat_openai(messages, project_id, tool_calls_log)


def _chat_anthropic(messages, project_id, tool_calls_log):
    """Anthropic API 格式对话"""
    import anthropic
    from modules.ai_tools import execute_tool

    client = anthropic.Anthropic(
        base_url=Config.AI_API_URL,
        api_key=Config.AI_API_KEY,
        timeout=60.0,
    )

    proposed_tools = []

    api_messages = []
    for msg in messages:
        if msg.get("role") in ("user", "assistant"):
            content = msg.get("content", "")
            if isinstance(content, list):
                # Anthropic 格式的 content blocks
                api_messages.append({"role": msg["role"], "content": content})
            else:
                api_messages.append({"role": msg["role"], "content": str(content)})

    if not api_messages:
        api_messages.append({"role": "user", "content": "你好"})

    # 第一次尝试：带 tools 调用
    try:
        response = client.messages.create(
            model=Config.AI_MODEL,
            system=SYSTEM_PROMPT,
            messages=api_messages,
            tools=TOOLS_ANTHROPIC,
            max_tokens=2048,
        )

        has_tool_use = False
        tool_results = []
        text_content = ""

        for block in response.content:
            if block.type == "text":
                text_content += block.text
            elif block.type == "tool_use":
                has_tool_use = True
                func_name = block.name
                args = block.input if isinstance(block.input, dict) else {}
                tool_calls_log.append({"name": func_name, "args": args})
                if func_name in CONFIRM_TOOLS:
                    proposed_tools.append({"name": func_name, "args": args})
                    result_str = json.dumps({"status": "pending_confirmation", "message": "等待用户确认"})
                else:
                    result = execute_tool(func_name, args, project_id)
                    result_str = json.dumps(result, ensure_ascii=False, default=str)
                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": result_str,
                })

        if has_tool_use:
            # 工具调用成功，继续对话循环
            api_messages.append({"role": "assistant", "content": response.content})
            api_messages.append({"role": "user", "content": tool_results})
            for _ in range(4):
                response = client.messages.create(
                    model=Config.AI_MODEL,
                    system=SYSTEM_PROMPT,
                    messages=api_messages,
                    tools=TOOLS_ANTHROPIC,
                    max_tokens=2048,
                )
                new_text = ""
                new_tool_results = []
                has_more_tools = False
                for block in response.content:
                    if block.type == "text":
                        new_text += block.text
                    elif block.type == "tool_use":
                        has_more_tools = True
                        fn = block.name
                        ar = block.input if isinstance(block.input, dict) else {}
                        tool_calls_log.append({"name": fn, "args": ar})
                        if fn in CONFIRM_TOOLS:
                            proposed_tools.append({"name": fn, "args": ar})
                            rs = json.dumps({"status": "pending_confirmation", "message": "等待用户确认"})
                        else:
                            r = execute_tool(fn, ar, project_id)
                            rs = json.dumps(r, ensure_ascii=False, default=str)
                        new_tool_results.append({"type": "tool_result", "tool_use_id": block.id, "content": rs})
                if not has_more_tools:
                    all_msgs = messages + [{"role": "assistant", "content": new_text}]
                    return {"reply": new_text, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools}
                api_messages.append({"role": "assistant", "content": response.content})
                api_messages.append({"role": "user", "content": new_tool_results})

            all_msgs = messages + [{"role": "assistant", "content": text_content}]
            return {"reply": text_content, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools}

        # 纯文本回复
        all_msgs = messages + [{"role": "assistant", "content": text_content}]
        return {"reply": text_content, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools}

    except Exception as e:
        return {"reply": f"AI 调用失败: {str(e)}", "tool_calls": tool_calls_log, "messages": messages, "proposed_tools": proposed_tools}


def _chat_openai(messages, project_id, tool_calls_log):
    """OpenAI 兼容 API 格式对话"""
    from openai import OpenAI
    from modules.ai_tools import execute_tool

    client = OpenAI(
        base_url=Config.AI_API_URL,
        api_key=Config.AI_API_KEY,
        timeout=60.0,
    )

    proposed_tools = []
    all_messages = [{"role": "system", "content": SYSTEM_PROMPT}] + messages

    for _ in range(5):
        response = client.chat.completions.create(
            model=Config.AI_MODEL,
            messages=all_messages,
            tools=TOOLS_OPENAI,
            tool_choice="auto",
            max_tokens=2048,
        )

        msg = response.choices[0].message
        all_messages.append(msg.model_dump())

        if not msg.tool_calls:
            return {
                "reply": msg.content or "",
                "tool_calls": tool_calls_log,
                "messages": all_messages[1:],
                "proposed_tools": proposed_tools,
            }

        for tc in msg.tool_calls:
            func_name = tc.function.name
            try:
                args = json.loads(tc.function.arguments)
            except json.JSONDecodeError:
                args = {}

            tool_calls_log.append({"name": func_name, "args": args})

            if func_name in CONFIRM_TOOLS:
                proposed_tools.append({"name": func_name, "args": args})
                result_str = json.dumps({"status": "pending_confirmation", "message": "等待用户确认"})
            else:
                result = execute_tool(func_name, args, project_id)
                result_str = json.dumps(result, ensure_ascii=False, default=str)

            all_messages.append({
                "role": "tool",
                "tool_call_id": tc.id,
                "content": result_str,
            })

    return {
        "reply": "抱歉，处理过程过于复杂，请简化您的请求。",
        "tool_calls": tool_calls_log,
        "messages": all_messages[1:],
        "proposed_tools": proposed_tools,
    }
