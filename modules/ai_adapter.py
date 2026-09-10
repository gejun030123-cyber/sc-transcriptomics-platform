# modules/ai_adapter.py
"""AI 对话适配器 — 支持 Anthropic Messages 和 OpenAI 兼容 API"""
import json
from config import Config
from modules.ai_config import get_effective_ai_config


# 工具定义：AI 可调用的平台功能
TOOLS_ANTHROPIC = [
    {
        "name": "run_analysis",
        "description": "执行一个分析模块。返回任务 ID 和状态。params 可使用网页表单中的全部用户可调参数；不确定键名、选项或联动条件时，先调用 get_module_parameters 获取当前模块的实际参数契约。",
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
        "name": "get_module_parameters",
        "description": "[只读] 获取某个分析模块与网页表单完全一致的用户可调参数、默认值、可选项、数值范围和联动条件。用户要求修改任何参数时先调用此工具，随后把用户确认的原值传给 run_analysis。",
        "input_schema": {
            "type": "object",
            "properties": {
                "module_name": {
                    "type": "string",
                    "description": "模块名称，如 sc_cell_go、bulk_deg、clustering"
                }
            },
            "required": ["module_name"]
        }
    },
    {
        "name": "run_pipeline",
        "description": "[需确认] 一次性提交多个有依赖关系的分析模块。后台会按 modules 的顺序串联前一步输出；不要用多个 run_analysis 代替全流程。返回 pipeline run ID，可查询进度和失败模块。",
        "input_schema": {
            "type": "object",
            "properties": {
                "analysis_type": {
                    "type": "string",
                    "enum": ["sc", "bulk"],
                    "description": "sc=单细胞，bulk=Bulk RNA-seq"
                },
                "modules": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "按依赖顺序排列的模块；例如 sc: qc, normalize, hvg, dimred, clustering, annotation"
                },
                "params": {
                    "type": "object",
                    "description": "按模块名分组的参数，例如 {\"clustering\": {\"resolutions\": \"1.0\"}}"
                },
                "input_path": {
                    "type": "string",
                    "description": "可选输入文件；留空时自动使用项目当前数据"
                },
                "name": {
                    "type": "string",
                    "description": "可选流程名称"
                }
            },
            "required": ["analysis_type", "modules"]
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
        "description": "获取某个分析任务的结果摘要，包括统计数据、生成的图表列表；如果存在 PNG/JPG/SVG 等图片，还会返回可在当前对话中展示的图片附件。",
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
                    "description": "要决策的分析模块，如 sc_pseudobulk_deg、functional_state、proportion、neighborhood_da、trajectory 或 bulk_deg"
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
        "name": "search_pathway_terms",
        "description": "[只读] 把自然语言生物学主题（如'脂代谢和炎症'）映射到本地基因集的具体通路 term。先向用户展示命中的 term 列表并确认，再用于 sc_cell_go 的 focus_terms 参数。不接触表达数据。",
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "自然语言主题或关键词，如 '脂代谢和炎症'、'lipid metabolism'"
                },
                "libraries": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "可选，限定检索的本地基因集库，如 ['GO_Biological_Process_2023','KEGG_2021_Human']；留空检索全部库"
                },
                "limit": {
                    "type": "integer",
                    "description": "最多返回的 term 数，默认 60"
                },
                "list_themes": {
                    "type": "boolean",
                    "description": "设为 true 时只列出平台支持的主题词典，不做检索"
                }
            },
            "required": []
        }
    },
    {
        "name": "read_task_table",
        "description": "[只读] 读取任务结果表（CSV/TSV/TXT/JSON/XLSX）的筛选摘要：支持列裁剪、关键词包含过滤、FDR 阈值和行数上限。用于在对话中解读富集、DEG、QC 等结果表内容。",
        "input_schema": {
            "type": "object",
            "properties": {
                "file_id": {
                    "type": "string",
                    "description": "结果文件 ID（来自 get_task_results 的 result_files）"
                },
                "columns": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "可选，只返回这些列"
                },
                "contains": {
                    "type": "string",
                    "description": "可选，行级关键词包含过滤（不区分大小写）"
                },
                "fdr_max": {
                    "type": "number",
                    "description": "可选，按 FDR/Adjusted P-value 列过滤小于该值的行"
                },
                "limit": {
                    "type": "integer",
                    "description": "最多返回行数，默认 50，上限 200"
                }
            },
            "required": ["file_id"]
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
    },
    {
        "name": "list_wes_workflows",
        "description": "[只读] 列出平台登记的 WES 工作流、支持的输入类型、模式和执行器状态。不会启动分析。",
        "input_schema": {
            "type": "object",
            "properties": {},
            "required": []
        }
    },
    {
        "name": "list_wes_references",
        "description": "[只读] 列出管理员登记的 WES reference asset，返回 assembly、bundle 版本和资源类型，不返回参考序列内容。",
        "input_schema": {
            "type": "object",
            "properties": {
                "assembly": {"type": "string", "description": "可选，如 GRCh38"},
                "bundle_version": {"type": "string", "description": "可选 reference bundle 版本"},
                "status": {"type": "string", "description": "可选，默认 registered"}
            },
            "required": []
        }
    },
    {
        "name": "inspect_wes_manifest",
        "description": "[只读] 对 WES 样本 manifest 做输入、样本角色、参考资源和文件路径预检查；不会启动外部工作流。优先使用 manifest_id，也可传项目内 manifest_path 或内联 manifest。",
        "input_schema": {
            "type": "object",
            "properties": {
                "workflow_key": {
                    "type": "string",
                    "description": "WES 工作流键，如 wes_germline、wes_somatic、wes_annotate_only"
                },
                "manifest_id": {
                    "type": "string",
                    "description": "已登记的项目 manifest ID"
                },
                "manifest_path": {
                    "type": "string",
                    "description": "项目目录内 JSON/CSV/TSV manifest 路径"
                },
                "manifest": {
                    "type": "object",
                    "description": "内联 manifest 对象；仅用于预检查，不会持久化"
                },
                "check_files": {
                    "type": "boolean",
                    "description": "是否检查输入文件存在，默认 true"
                },
                "check_content": {
                    "type": "boolean",
                    "description": "是否使用已安装的 pysam 做 BAM/CRAM/VCF/FASTQ 内容级检查，默认 false"
                }
            },
            "required": ["workflow_key"]
        }
    },
    {
        "name": "inspect_wes_preflight",
        "description": "[只读] WES manifest 预检查别名，用于确认样本配对、文件类型和参考资源；不执行 Nextflow。",
        "input_schema": {
            "type": "object",
            "properties": {
                "workflow_key": {"type": "string"},
                "manifest_id": {"type": "string"},
                "manifest_path": {"type": "string"},
                "manifest": {"type": "object"},
                "check_files": {"type": "boolean"},
                "check_content": {"type": "boolean"}
            },
            "required": ["workflow_key"]
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
    'get_module_parameters',
    'inspect_analysis_state', 'inspect_adata', 'get_cluster_summary',
    'score_cell_type_signature', 'list_builtin_markers',
    'recommend_analysis_config',
    'search_pathway_terms', 'read_task_table',
    'list_wes_workflows', 'list_wes_references', 'inspect_wes_manifest', 'inspect_wes_preflight',
}
# 需要用户确认的工具
CONFIRM_TOOLS = {'run_analysis', 'run_pipeline', 'propose_parameter_sweep', 'run_parameter_sweep',
                 'start_goal_agent', 'continue_goal_agent'}
# 注：accept_branch 仅通过前端 Branch API 调用（POST /api/branches/<id>/accept），
# 不作为 AI 工具暴露，确保用户在前端显式操作采纳。


# System prompt
SYSTEM_PROMPT = """你是一个生信分析助手，帮助用户进行 RNA-seq 和 WES 数据分析。你的平台是基于 Flask 的 Web 应用，支持单细胞、Bulk RNA-seq 全流程，以及 WES 输入预检查与工作流登记。

## 你的能力
1. **自然语言触发分析**：用户说"对 hmc3 和 ctrl 做差异分析"，你调用 run_analysis(module_name="bulk_deg", params={{...}})
2. **结果解读**：用户问"哪些基因在所有药物中共同上调？"，你调用 get_task_results 查看结果并解读
3. **方法与参数决策**：调用 recommend_analysis_config，用实际数据画像区分 raw counts、连续表达值和已 log 数据，再推荐方法、参数和分组
4. **项目概览**：用户问"现在分析到哪一步了？"，你调用 get_project_status 了解情况
5. **分析状态检查**：使用 inspect_analysis_state 查看当前聚类/嵌入/注释信息
6. **细胞类型打分**：使用 score_cell_type_signature 对已有分群进行 marker 评分
7. **目标优化 Agent**：用户指定细胞类型，你使用 start_goal_agent 自动检查、生成候选参数、评分并推荐最佳分群
8. **WES 输入审阅**：使用 list_wes_workflows、list_wes_references 和 inspect_wes_manifest 检查清单、配对、参考资源和文件能力；这些工具只读，不启动 WES。
9. **主题驱动富集**：用户表达生物学主题（如"脂代谢和炎症"）时，使用 search_pathway_terms 映射到具体通路，再用 sc_cell_go 的 focus_terms 运行主题优先 GO 富集
10. **参数契约查询**：用户要求修改网页中的任意参数时，使用 get_module_parameters 获取与网页一致的实际键名、范围和联动条件；不得凭记忆编造参数名或丢弃用户已确认的值。

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
5. 不得把 DESeq2/edgeR 用于 FPKM/TPM 连续值；不得对已 log 数据重复标准化；不得在无时间列时推荐时序分析；不得把单个细胞当作多时间点的独立生物学重复。

## 用户手动参数调整
当用户说“把某参数改为…/按网页参数运行/调整阈值、图形、聚类、模型或导出设置”时：
1. 先调用 get_module_parameters(module_name=...)，以返回的 key、options、min/max、show_if/depends_on 为唯一参数契约；这覆盖所有模块的所有网页可调字段。
2. 将用户指定值与必要的控制字段一起回显；不得把用户值替换成推荐默认值，也不得把受 show_if 控制的有效字段静默丢弃。
3. 对会启动任务的改动，获得用户确认后再调用 run_analysis，并原样传入已确认的用户可调参数；服务端仍会执行范围、联动和方法兼容性校验。
4. 若参数不在契约中或违反范围/联动，解释原因并让用户选择有效值；不得猜测近似参数。

## 全流程执行
当用户要求“全流程/一键完成/从头跑到结果”时：先检查数据和必要元数据，并给出模块顺序与关键参数；获得一次确认后，必须调用 **run_pipeline** 一次性提交整个流程，不能逐个调用 run_analysis。流程在后台按顺序等待每一步完成后再执行下一步；回复中说明可通过 pipeline run 状态查看进度和失败位置。若设计检查表明后续模块缺少分组、比较或时间元数据，只提交可安全执行的核心流程，并明确说明未提交的模块和原因。

## 可用模块（完整列表）
单细胞：qc, normalize, hvg, dimred, batch_correct, clustering, qc_reassess, annotation, functional_state, sc_timecourse, deg, sc_pseudobulk_deg, sc_cell_go, trajectory, proportion, neighborhood_da, cell_communication
Bulk：bulk_qc, bulk_normalize, bulk_deg, bulk_pca, bulk_heatmap, bulk_enrichment, bulk_timecourse, bulk_deg_integration

## IBD / 对照类器官单细胞下游流程
当用户的目标是“炎症/应激增强、TA/增殖改变、成熟吸收/代谢 enterocyte 下降”或相近 IBD 类器官问题时，按以下证据层次推荐，不要把不同统计单位混为一谈：
1. 同一细胞类型内部 disease-vs-control 表达变化：先推荐 sc_pseudobulk_deg；它要求原始 counts、sample_id、condition 和 celltype/annotation，按 sample × celltype 聚合。不得以 sc_cell_deg 的细胞级 p 值替代正式结论。
2. 炎症、TNF/NF-kB、IFN、hypoxia、ROS、UPR、apoptosis、FAO/peroxisome/OXPHOS、胆固醇/脂代谢、WNT、细胞周期等状态：推荐 functional_state；IBD/肠炎目标可选择 analysis_focus=ibd_organoid_epithelial。先检查本地冻结基因集是否可用，并将 cell-level 图解释为描述性，样本 × celltype 比较才是正式统计。
3. 细胞组成改变：推荐 proportion，并显式设置 analysis_unit=sample、sample_key、condition_key 和 celltype groupby；不得把细胞数卡方检验作为最终生物学重复证据。
4. 连续上皮状态中 disease 富集的局部区域：推荐 neighborhood_da。只在存在 X_pca/Harmony/scVI 等高维表示、且各条件有独立样本时使用；不可用二维 UMAP 距离代替，也不得称为完整 R/Milo 分析。
5. Stem/progenitor → TA → absorptive/metabolic 或 secretory 分化：推荐 trajectory 的 PAGA/DPT，并要求用户确认起始 Stem/progenitor 状态。RNA velocity 只有在明确存在 spliced/unspliced 信息且平台具备专用工作流时才能声称可做；当前 trajectory 不能替代 velocity。
6. 每个主要 celltype 的 DEG 完成后，使用 sc_cell_go 绑定该明确的 pseudobulk 任务 ID 做 Hallmark/Reactome/GO/KEGG 富集；不得扫描“最新 CSV”猜测 DEG 来源。优先关注 lipid absorption、cholesterol/FAO/peroxisome 与 inflammatory/TNF-NF-kB/hypoxia。
7. 上皮为主的类器官中，cell_communication 是低优先级；除非用户明确要做或存在 immune/stromal 细胞，不能把通讯当作主结论。
8. 每次要执行上述模块前，先调用 recommend_analysis_config；它会依据当前元数据和表示返回 should_run。若缺少原始 counts、sample/condition、注释或高维表示，说明缺口并停止提交。用户确认后再调用 run_analysis。

## 主题驱动富集流程
当用户希望富集结果聚焦某个生物学主题（如"跑完差异表达了，我要关于脂代谢和炎症的富集"）时：
1. 调用 search_pathway_terms(query="用户主题原话")；不确定平台支持哪些主题时，先传 list_themes=true 查看。
2. 把命中的通路 term 按库分组展示给用户，说明每个库命中数量，请用户确认或增删；不要未经确认直接提交分析。
3. 用户确认后，调用 run_analysis(module_name="sc_cell_go", params={{"deg_source_task_id": "<已完成的DEG任务ID>", "method": "ORA", "focus_terms": ["term1", "term2"]}})。focus_terms 不会改变全量统计：富集照常在全部通路上运行，FDR 在全库上校正；GO 图在每个本体内优先炎症、再脂代谢/确认主题，并以普通 Top 通路补足，输出可审计的三分区气泡图与柱状图。
4. 用户要求调整展示名额时，先调用 get_module_parameters(module_name="sc_cell_go")；提交时必须保留 go_priority_allocation_mode="custom" 及相应的 go_inflammation_slots、go_lipid_slots、go_confirmed_theme_slots、go_top_pathway_slots。四者之和不得超过 plot_top_n；这些参数只改变展示名额，不改变 ORA 或 FDR。
5. 任务完成后，用 get_task_results 找到富集结果表 file_id，再调用 read_task_table(file_id=..., contains=..., fdr_max=...) 摘要主题结果并解读。
6. 解读时明确：主题优先图不是独立的校正检验，仍要结合完整结果表；某主题未命中时仍展示 FDR Top 通路，某本体没有任何显著通路时才如实标为空，不得夸大。

## 回复规则
- 用中文回复
- 简洁明了，不啰嗦
- 执行分析/参数搜索等写操作时，先确认参数再执行
- 分析完成后，简要总结关键结果
- 如果用户要求不明确，主动询问关键参数
- 对需要确认的工具调用，明确告知用户需要批准才能执行"""


def _is_anthropic():
    """判断是否使用 Anthropic API"""
    settings = get_effective_ai_config()
    provider = settings.get("provider", "auto")
    if provider == "anthropic":
        return True
    if provider == "openai":
        return False
    url = settings["api_url"].lower()
    return 'anthropic' in url or 'claude' in url


def _provider_name():
    settings = get_effective_ai_config()
    url = settings["api_url"].lower()
    if 'deepseek' in url:
        return 'deepseek-anthropic' if _is_anthropic() else 'deepseek-openai'
    return 'anthropic' if _is_anthropic() else 'openai-compatible'


def _safe_api_url():
    """Return API URL without query strings for display/logging."""
    return get_effective_ai_config()["api_url"].split('?', 1)[0]


def get_ai_config_status():
    """Expose non-secret AI configuration status for UI diagnostics."""
    settings = get_effective_ai_config()
    return {
        "configured": bool(settings["api_key"]),
        "provider": _provider_name(),
        "mode": "anthropic_messages" if _is_anthropic() else "openai_chat_completions",
        "api_url": _safe_api_url(),
        "model": settings["model"],
        "requires_local_token": bool(Config.AI_API_TOKEN),
    }


def _merge_attachments(target, result):
    """Collect safe result-image metadata from tool responses for the UI."""
    if not isinstance(result, dict):
        return
    for attachment in result.get("attachments", []) or []:
        if not isinstance(attachment, dict) or not attachment.get("url"):
            continue
        key = attachment.get("id") or attachment["url"]
        if not any((item.get("id") or item.get("url")) == key for item in target):
            target.append({
                "id": attachment.get("id", ""),
                "task_id": attachment.get("task_id", ""),
                "label": attachment.get("label", "分析图"),
                "category": attachment.get("category", "plot"),
                "type": attachment.get("type", "png"),
                "url": attachment["url"],
            })


def chat(messages, project_id=None):
    """
    与 LLM 对话，支持工具调用循环。自动检测 Anthropic / OpenAI 格式。
    """
    tool_calls_log = []
    # Attachment metadata is for the browser only; do not send the extra
    # fields back to OpenAI-compatible clients as message properties.
    model_messages = [
        {key: value for key, value in message.items() if key != "attachments"}
        for message in messages
    ]

    if _is_anthropic():
        return _chat_anthropic(model_messages, project_id, tool_calls_log)
    else:
        return _chat_openai(model_messages, project_id, tool_calls_log)


def _chat_anthropic(messages, project_id, tool_calls_log):
    """Anthropic Messages API 格式对话.

    这里直接用 HTTP 调用，避免 DeepSeek Anthropic-compatible endpoint
    依赖本地安装 anthropic SDK。
    """
    from modules.ai_tools import execute_tool

    proposed_tools = []
    attachments = []

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
                    _merge_attachments(attachments, result)
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
                            _merge_attachments(attachments, r)
                            rs = json.dumps(r, ensure_ascii=False, default=str)
                        new_tool_results.append({"type": "tool_result", "tool_use_id": block.get("id"), "content": rs})
                if not has_more_tools:
                    all_msgs = messages + [{"role": "assistant", "content": new_text}]
                    return {"reply": new_text, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools, "attachments": attachments}
                api_messages.append({"role": "assistant", "content": content_blocks})
                api_messages.append({"role": "user", "content": new_tool_results})

            reply_text = new_text or "工具调用已执行，但回复生成超出轮次限制。请查看任务状态了解结果。"
            all_msgs = messages + [{"role": "assistant", "content": reply_text}]
            return {"reply": reply_text, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools, "attachments": attachments}

        # 纯文本回复
        all_msgs = messages + [{"role": "assistant", "content": text_content}]
        return {"reply": text_content, "tool_calls": tool_calls_log, "messages": all_msgs, "proposed_tools": proposed_tools, "attachments": attachments}

    except Exception as e:
        return {"reply": f"AI 调用失败: {str(e)}", "tool_calls": tool_calls_log, "messages": messages, "proposed_tools": proposed_tools, "attachments": attachments}


def _anthropic_messages_endpoint():
    base_url = get_effective_ai_config()["api_url"].rstrip("/")
    if base_url.endswith("/v1/messages") or base_url.endswith("/messages"):
        return base_url
    if base_url.endswith("/v1"):
        return base_url + "/messages"
    return base_url + "/v1/messages"


def _anthropic_messages_create(api_messages):
    import requests

    settings = get_effective_ai_config()

    headers = {
        "Content-Type": "application/json",
        "x-api-key": settings["api_key"],
        "Authorization": f"Bearer {settings['api_key']}",
        "anthropic-version": "2023-06-01",
    }
    payload = {
        "model": settings["model"],
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

    settings = get_effective_ai_config()

    client = OpenAI(
        base_url=settings["api_url"],
        api_key=settings["api_key"],
        timeout=60.0,
    )

    proposed_tools = []
    attachments = []
    all_messages = [{"role": "system", "content": SYSTEM_PROMPT}] + messages

    for _ in range(5):
        response = client.chat.completions.create(
            model=settings["model"],
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
                "attachments": attachments,
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
                _merge_attachments(attachments, result)
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
        "attachments": attachments,
    }
