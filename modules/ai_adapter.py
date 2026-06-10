# modules/ai_adapter.py
"""AI 对话适配器 — 支持任意 OpenAI 兼容 API"""
import json
from config import Config


def get_client():
    """获取 OpenAI 兼容客户端"""
    from openai import OpenAI
    return OpenAI(
        base_url=Config.AI_API_URL,
        api_key=Config.AI_API_KEY,
    )


# 工具定义：AI 可调用的平台功能
TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "run_analysis",
            "description": "执行一个分析模块。返回任务 ID 和状态。",
            "parameters": {
                "type": "object",
                "properties": {
                    "module_name": {
                        "type": "string",
                        "description": "分析模块名称，如 bulk_deg, bulk_heatmap, bulk_normalize 等"
                    },
                    "input_path": {
                        "type": "string",
                        "description": "输入数据文件路径（h5ad/csv）。留空则使用项目中最新可用的中间文件。"
                    },
                    "params": {
                        "type": "object",
                        "description": "分析参数，如 {\"method\": \"deseq2\", \"fc_threshold\": 1.5}。留空使用默认值。"
                    }
                },
                "required": ["module_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "get_project_status",
            "description": "获取当前项目的完整状态：已上传文件、已完成任务、可用的中间文件。",
            "parameters": {
                "type": "object",
                "properties": {},
                "required": []
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "get_task_results",
            "description": "获取某个分析任务的结果摘要，包括统计数据和生成的图表列表。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "任务 ID"
                    }
                },
                "required": ["task_id"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "list_modules",
            "description": "列出所有可用的分析模块及其说明和参数。",
            "parameters": {
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
    }
]


# System prompt
SYSTEM_PROMPT = """你是一个生信分析助手，帮助用户进行 RNA-seq 数据分析。你的平台是基于 Flask 的 Web 应用，支持单细胞和 Bulk RNA-seq 全流程分析。

## 你的能力
1. **自然语言触发分析**：用户说"对 hmc3 和 ctrl 做差异分析"，你调用 run_analysis(module_name="bulk_deg", params={...})
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


def chat(messages, project_id=None):
    """
    与 LLM 对话，支持工具调用循环。

    Args:
        messages: 对话历史 [{"role": "user", "content": "..."}]
        project_id: 当前项目 ID（工具执行时需要）

    Returns:
        {"reply": str, "tool_calls": list, "messages": list}
    """
    from modules.ai_tools import execute_tool

    client = get_client()
    all_messages = [{"role": "system", "content": SYSTEM_PROMPT}] + messages
    tool_calls_log = []

    for _ in range(5):  # 最多 5 轮工具调用循环
        response = client.chat.completions.create(
            model=Config.AI_MODEL,
            messages=all_messages,
            tools=TOOLS,
            tool_choice="auto",
            max_tokens=2048,
        )

        msg = response.choices[0].message
        all_messages.append(msg.model_dump())

        if not msg.tool_calls:
            # 无工具调用，返回最终回复
            return {
                "reply": msg.content or "",
                "tool_calls": tool_calls_log,
                "messages": all_messages[1:],  # 去掉 system prompt
            }

        # 执行工具调用
        for tc in msg.tool_calls:
            func_name = tc.function.name
            try:
                args = json.loads(tc.function.arguments)
            except json.JSONDecodeError:
                args = {}

            tool_calls_log.append({"name": func_name, "args": args})

            try:
                result = execute_tool(func_name, args, project_id)
                result_str = json.dumps(result, ensure_ascii=False, default=str)
            except Exception as e:
                result_str = json.dumps({"error": str(e)}, ensure_ascii=False)

            all_messages.append({
                "role": "tool",
                "tool_call_id": tc.id,
                "content": result_str,
            })

    # 超过循环次数
    return {
        "reply": "抱歉，处理过程过于复杂，请简化您的请求。",
        "tool_calls": tool_calls_log,
        "messages": all_messages[1:],
    }
