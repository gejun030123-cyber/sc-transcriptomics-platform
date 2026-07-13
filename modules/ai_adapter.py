# modules/ai_adapter.py
"""AI 对话适配器 — 支持 Anthropic Messages 和 OpenAI 兼容 API"""
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
                    "description": "输入数据文件路径（h5ad/csv/tsv/txt/xlsx/xls）。留空则自动查找最新可用文件。"
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
    },
    {
        "name": "recommend_analysis_config",
        "description": "[只读] 根据实际输入数据的 count/连续表达类型、样本数、分组结构、已完成步骤和用户目标，推荐分析方法与参数。返回数据证据、推荐值、理由、替代方案、前置步骤和风险；不执行分析。",
        "input_schema": {
            "type": "object",
            "properties": {
                "module_name": {
                    "type": "string",
                    "description": "要决策的分析模块，如 bulk_normalize、bulk_deg、clustering、batch_correct"
                },
                "input_path": {
                    "type": "string",
                    "description": "可选输入路径。留空时按模块自动选择原始上传或最新中间文件"
                },
                "objective": {
                    "type": "string",
                    "description": "用户目标或约束，如保留 B/En 差异、优先发现稀有细胞、避免过度校正"
                }
            },
            "required": ["module_name"]
        }
    },
    {
        "name": "inspect_analysis_state",
        "description": "[只读] 检查项目当前分析状态：最新 h5ad、已完成任务、可用聚类键、可用嵌入、细胞数等。用于了解分析进展。",
        "input_schema": {
            "type": "object",
            "properties": {},
            "required": []
        }
    },
    {
        "name": "inspect_adata",
        "description": "[只读] 检查 AnnData 或 Bulk 表达矩阵的结构：样本数、基因数、样本名及可用元数据。",
        "input_schema": {
            "type": "object",
            "properties": {
                "adata_path": {
                    "type": "string",
                    "description": "h5ad 或 Bulk 表格文件路径。留空则自动查找最新文件。"
                }
            },
            "required": []
        }
    },
    {
        "name": "get_cluster_summary",
        "description": "[只读] 获取指定聚类键的每个 cluster 概况：细胞数、UMAP中心、batch分布、QC均值。",
        "input_schema": {
            "type": "object",
            "properties": {
                "adata_path": {
                    "type": "string",
                    "description": "h5ad 文件路径"
                },
                "cluster_key": {
                    "type": "string",
                    "description": "聚类键，如 leiden, louvain。默认 leiden。"
                }
            },
            "required": ["adata_path"]
        }
    },
    {
        "name": "score_cell_type_signature",
        "description": "[只读] 基于 marker 基因对每个 cluster 进行细胞类型签名评分。返回最佳 cluster 及置信度。用于寻找最接近特定细胞类型的分群。",
        "input_schema": {
            "type": "object",
            "properties": {
                "adata_path": {
                    "type": "string",
                    "description": "h5ad 文件路径。留空则自动查找最新。"
                },
                "cluster_key": {
                    "type": "string",
                    "description": "聚类键，默认 leiden"
                },
                "target_cell_type": {
                    "type": "string",
                    "description": "目标细胞类型：microglia, t_cell, b_cell, nk, monocyte, macrophage, dendritic_cell, epithelial, endothelial, fibroblast, astrocyte, oligodendrocyte"
                },
                "positive_markers": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "用户自定义 positive marker 基因（优先于内置库）"
                },
                "negative_markers": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "用户自定义 negative marker 基因（优先于内置库）"
                }
            },
            "required": ["target_cell_type"]
        }
    },
    {
        "name": "list_builtin_markers",
        "description": "[只读] 列出所有内置细胞类型及其 marker 基因定义。",
        "input_schema": {
            "type": "object",
            "properties": {},
            "required": []
        }
    },
    {
        "name": "propose_parameter_sweep",
        "description": "[需确认] 为优化特定细胞类型分群生成参数搜索候选列表。基于当前分析状态智能选择参数空间。",
        "input_schema": {
            "type": "object",
            "properties": {
                "goal_type": {
                    "type": "string",
                    "description": "目标类型，如 target_cluster_refinement",
                    "default": "target_cluster_refinement"
                },
                "target_cell_type": {
                    "type": "string",
                    "description": "目标细胞类型，如 microglia"
                },
                "max_candidates": {
                    "type": "integer",
                    "description": "最大候选数，默认6，上限12"
                }
            },
            "required": ["target_cell_type"]
        }
    },
    {
        "name": "run_parameter_sweep",
        "description": "[需确认] 执行参数搜索：创建候选分支、运行聚类、自动评分。需要 goal_id 和候选列表。",
        "input_schema": {
            "type": "object",
            "properties": {
                "goal_id": {
                    "type": "string",
                    "description": "目标 ID（来自 start_goal_agent 的返回值）"
                },
                "base_checkpoint": {
                    "type": "string",
                    "description": "基础 h5ad 路径（所有候选的起点）"
                },
                "candidates": {
                    "type": "array",
                    "items": {"type": "object"},
                    "description": "候选参数列表（来自 propose_parameter_sweep）"
                }
            },
            "required": ["goal_id", "base_checkpoint", "candidates"]
        }
    },
    {
        "name": "start_goal_agent",
        "description": "[需确认] 启动目标驱动的分析 Agent 会话。自动检查当前分群状态、评分、决定是否需要参数搜索。",
        "input_schema": {
            "type": "object",
            "properties": {
                "target_cell_type": {
                    "type": "string",
                    "description": "目标细胞类型：microglia, t_cell, b_cell, nk, monocyte, macrophage, dendritic_cell, epithelial, endothelial, fibroblast, astrocyte, oligodendrocyte"
                },
                "goal_type": {
                    "type": "string",
                    "description": "目标类型，默认 target_cluster_refinement",
                    "default": "target_cluster_refinement"
                },
                "positive_markers": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "用户自定义 positive marker 基因"
                },
                "negative_markers": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "用户自定义 negative marker 基因"
                },
                "user_requirement": {
                    "type": "string",
                    "description": "用户自然语言需求（如'找到最接近小胶质细胞的分群'）"
                },
                "max_candidate_runs": {
                    "type": "integer",
                    "description": "最大候选运行数，默认6"
                }
            },
            "required": ["target_cell_type"]
        }
    },
    {
        "name": "continue_goal_agent",
        "description": "[需确认] 继续目标分析会话。支持指令：采用候选X、继续细分、不满意/更高resolution、停止。",
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {
                    "type": "string",
                    "description": "会话 ID"
                },
                "instruction": {
                    "type": "string",
                    "description": "用户指令：采用候选X / 继续细分 / 不满意，试试更高resolution / 停止 / 必须包含某基因 / 不用某方法"
                }
            },
            "required": ["session_id", "instruction"]
        }
    }
]

TOOLS_OPENAI = [
    {"type": "function", "function": {**t, "parameters": t["input_schema"]}}
    for t in TOOLS_ANTHROPIC
]

# 自动执行的只读工具（无需用户确认）
AUTO_EXEC_TOOLS = {
    'get_project_status', 'get_task_results', 'list_modules',
    'inspect_analysis_state', 'inspect_adata', 'get_cluster_summary',
    'score_cell_type_signature', 'list_builtin_markers',
    'recommend_analysis_config',
}
# 需要用户确认的工具
CONFIRM_TOOLS = {'run_analysis', 'propose_parameter_sweep', 'run_parameter_sweep',
                 'start_goal_agent', 'continue_goal_agent'}
# 注：accept_branch 仅通过前端 Branch API 调用（POST /api/branches/<id>/accept），
# 不作为 AI 工具暴露，确保用户在前端显式操作采纳。


# System prompt
SYSTEM_PROMPT = """你是一个生信分析助手，帮助用户进行 RNA-seq 数据分析。你的平台是基于 Flask 的 Web 应用，支持单细胞和 Bulk RNA-seq 全流程分析。

## 你的能力
1. **自然语言触发分析**：用户说"对 hmc3 和 ctrl 做差异分析"，你调用 run_analysis(module_name="bulk_deg", params={{...}})
2. **结果解读**：用户问"哪些基因在所有药物中共同上调？"，你调用 get_task_results 查看结果并解读
3. **方法与参数决策**：调用 recommend_analysis_config，用实际数据画像区分 raw counts、连续表达值和已 log 数据，再推荐方法、参数和分组
4. **项目概览**：用户问"现在分析到哪一步了？"，你调用 get_project_status 了解情况
5. **分析状态检查**：使用 inspect_analysis_state 查看当前聚类/嵌入/注释信息
6. **细胞类型打分**：使用 score_cell_type_signature 对已有分群进行 marker 评分
7. **目标优化 Agent**：用户指定细胞类型，你使用 start_goal_agent 自动检查、生成候选参数、评分并推荐最佳分群

## 目标优化 Agent 使用流程
当用户表示对分群不满意或想找特定细胞类型时：
1. 调用 start_goal_agent(target_cell_type="microglia", user_requirement="用户原话")
2. 如果返回 needs_confirmation，解释当前评分和建议的候选参数，请用户确认
3. 用户确认后，调用 run_parameter_sweep 执行搜索
4. 展示结果：最佳候选、评分、推荐 cluster
5. 用户选择采纳 branch 或继续调整

## 方法与参数推进流程
当用户询问“怎么分析”“用什么方法”“帮我设置参数”或要求运行新模块时：
1. 先调用 recommend_analysis_config(module_name, objective)，不得只凭文件名或通用经验猜测。
2. 向用户展示：数据类型证据、推荐方法、完整参数、分组/比较、替代方法、前置步骤与风险。
3. 如 should_run=false，停止提交并说明缺少的元数据或不适用原因。
4. 只有用户确认后，才使用 recommend_analysis_config 返回的 input_path 和 recommended_params 调用 run_analysis。
5. 不得把 DESeq2/edgeR 用于 FPKM/TPM 连续值；不得对已 log 数据重复标准化；不得在无时间列时推荐时序分析。

## 可用模块（完整列表）
单细胞：qc, normalize, hvg, dimred, batch_correct, clustering, qc_reassess, annotation, deg, trajectory, proportion, cell_communication
Bulk：bulk_qc, bulk_normalize, bulk_deg, bulk_pca, bulk_heatmap, bulk_enrichment, bulk_timecourse, bulk_deg_integration

## 回复规则
- 用中文回复
- 简洁明了，不啰嗦
- 执行分析/参数搜索等写操作时，先确认参数再执行
- 分析完成后，简要总结关键结果
- 如果用户要求不明确，主动询问关键参数
- 对需要确认的工具调用，明确告知用户需要批准才能执行"""


def _is_anthropic():
    """判断是否使用 Anthropic API"""
    url = Config.AI_API_URL.lower()
    return 'anthropic' in url or 'claude' in url


def _provider_name():
    url = Config.AI_API_URL.lower()
    if 'deepseek' in url:
        return 'deepseek-anthropic' if _is_anthropic() else 'deepseek-openai'
    return 'anthropic' if _is_anthropic() else 'openai-compatible'


def _safe_api_url():
    """Return API URL without query strings for display/logging."""
    return Config.AI_API_URL.split('?', 1)[0]


def get_ai_config_status():
    """Expose non-secret AI configuration status for UI diagnostics."""
    return {
        "configured": bool(Config.AI_API_KEY),
        "provider": _provider_name(),
        "mode": "anthropic_messages" if _is_anthropic() else "openai_chat_completions",
        "api_url": _safe_api_url(),
        "model": Config.AI_MODEL,
        "requires_local_token": bool(Config.AI_API_TOKEN),
    }


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
    """Anthropic Messages API 格式对话.

    这里直接用 HTTP 调用，避免 DeepSeek Anthropic-compatible endpoint
    依赖本地安装 anthropic SDK。
    """
    from modules.ai_tools import execute_tool

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
        response = _anthropic_messages_create(api_messages)
        content_blocks = response.get("content", [])

        has_tool_use = False
        tool_results = []
        text_content = ""

        for block in content_blocks:
            if block.get("type") == "text":
                text_content += block.get("text", "")
            elif block.get("type") == "tool_use":
                has_tool_use = True
                func_name = block.get("name", "")
                args = block.get("input") if isinstance(block.get("input"), dict) else {}
                tool_calls_log.append({"name": func_name, "args": args})
                if func_name in CONFIRM_TOOLS:
                    proposed_tools.append({"name": func_name, "args": args})
                    result_str = json.dumps({"status": "pending_confirmation", "message": "等待用户确认"})
                else:
                    result = execute_tool(func_name, args, project_id)
                    result_str = json.dumps(result, ensure_ascii=False, default=str)
                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": block.get("id"),
                    "content": result_str,
                })

        if has_tool_use:
            # 工具调用成功，继续对话循环
            api_messages.append({"role": "assistant", "content": content_blocks})
            api_messages.append({"role": "user", "content": tool_results})
            for _ in range(4):
                response = _anthropic_messages_create(api_messages)
                content_blocks = response.get("content", [])
                new_text = ""
                new_tool_results = []
                has_more_tools = False
                for block in content_blocks:
                    if block.get("type") == "text":
                        new_text += block.get("text", "")
                    elif block.get("type") == "tool_use":
                        has_more_tools = True
                        fn = block.get("name", "")
                        ar = block.get("input") if isinstance(block.get("input"), dict) else {}
                        tool_calls_log.append({"name": fn, "args": ar})
                        if fn in CONFIRM_TOOLS:
                            proposed_tools.append({"name": fn, "args": ar})
                            rs = json.dumps({"status": "pending_confirmation", "message": "等待用户确认"})
                        else:
                            r = execute_tool(fn, ar, project_id)
                            rs = json.dumps(r, ensure_ascii=False, default=str)
                        new_tool_results.append({"type": "tool_result", "tool_use_id": block.get("id"), "content": rs})
                if not has_more_tools:
                    all_msgs = messages + [{"role": "assistant", "content": new_text}]
                    return {"reply": new_text, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools}
                api_messages.append({"role": "assistant", "content": content_blocks})
                api_messages.append({"role": "user", "content": new_tool_results})

            reply_text = new_text or "工具调用已执行，但回复生成超出轮次限制。请查看任务状态了解结果。"
            all_msgs = messages + [{"role": "assistant", "content": reply_text}]
            return {"reply": reply_text, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools}

        # 纯文本回复
        all_msgs = messages + [{"role": "assistant", "content": text_content}]
        return {"reply": text_content, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools}

    except Exception as e:
        return {"reply": f"AI 调用失败: {str(e)}", "tool_calls": tool_calls_log, "messages": messages, "proposed_tools": proposed_tools}


def _anthropic_messages_endpoint():
    base_url = Config.AI_API_URL.rstrip("/")
    if base_url.endswith("/v1/messages") or base_url.endswith("/messages"):
        return base_url
    return base_url + "/v1/messages"


def _anthropic_messages_create(api_messages):
    import requests

    headers = {
        "Content-Type": "application/json",
        "x-api-key": Config.AI_API_KEY,
        "Authorization": f"Bearer {Config.AI_API_KEY}",
        "anthropic-version": "2023-06-01",
    }
    payload = {
        "model": Config.AI_MODEL,
        "system": SYSTEM_PROMPT,
        "messages": api_messages,
        "tools": TOOLS_ANTHROPIC,
        "max_tokens": 2048,
    }
    response = requests.post(
        _anthropic_messages_endpoint(),
        headers=headers,
        json=payload,
        timeout=60.0,
    )
    if response.status_code >= 400:
        detail = response.text[:500]
        try:
            body = response.json()
            detail = body.get("error", body)
        except ValueError:
            pass
        raise RuntimeError(f"Anthropic-compatible API 返回 {response.status_code}: {detail}")
    return response.json()


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
