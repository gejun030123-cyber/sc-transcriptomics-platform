"""主题 → 通路 term 的确定性映射（Theme lexicon）。

把中文/英文的自然语言生物学主题（如"脂代谢""炎症"）映射到平台本地
GMT/TXT 基因集快照中的具体通路 term，供 AI 的 ``search_pathway_terms``
工具和 ``sc_cell_go`` 的 ``focus_terms`` 参数使用。

约束：

* 匹配完全在服务器本地进行，对 term 名称做确定性关键词匹配；不接触
  表达数据，不向外部服务发送任何内容。
* 基因集文件与 ``sc_cell_go`` 实际富集运行使用同一来源
  （``data/go_gene_sets`` 平台快照），保证返回的 term 名称可以直接作为
  ``focus_terms`` 精确命中富集结果表。
* 找不到主题或 term 时如实返回空结果，绝不猜测或伪造通路名称。
"""

import re
from pathlib import Path

from config import Config

# 与 modules/schemas.py 中 sc_cell_go 的本地基因集选项保持一致。
DEFAULT_LIBRARIES = (
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
    "WikiPathway_2021_Human",
)

_LIBRARY_NAME_RE = re.compile(r"[A-Za-z0-9_.-]+")

# 每个主题包含中文标签、中文关键词（用于把用户输入归到主题）和英文关键词
# （用于匹配本地基因集的英文 term 名称）。英文关键词刻意保持宽松的子串
# 语义，例如 "inflamm" 可同时命中 inflammation / inflammatory。
THEMES = {
    "lipid_metabolism": {
        "label": "脂代谢",
        "description": "脂肪酸、胆固醇、磷脂、脂蛋白等脂质代谢与转运",
        "keywords_en": [
            "lipid", "fatty acid", "acyl", "sterol", "cholesterol",
            "triglyceride", "phospholipid", "sphingolipid", "ceramide",
            "lipoprotein", "lipase", "glycerolipid", "eicosanoid",
            "prostaglandin", "leukotriene", "ketone", "ppar",
        ],
        "keywords_zh": ["脂", "脂质", "脂肪酸", "胆固醇", "脂代谢", "脂质代谢", "甘油"],
    },
    "inflammation": {
        "label": "炎症",
        "description": "炎症反应、细胞因子、干扰素与炎症小体相关通路",
        "keywords_en": [
            "inflamm", "cytokine", "chemokine", "interleukin", "interferon",
            "tumor necrosis factor", "nf-kappa", "toll-like", "nod-like",
            "nod-like receptor", "il-1", "il-6", "il-17", "acute phase",
            "inflammasome",
        ],
        "keywords_zh": ["炎症", "发炎", "炎性", "细胞因子", "干扰素", "白介素"],
    },
    "immune_response": {
        "label": "免疫应答",
        "description": "免疫细胞活化、抗原呈递、抗体与补体通路",
        "keywords_en": [
            "immune", "immunoglobulin", "antigen", "leukocyte", "lymphocyte",
            "t cell", "b cell", "nk cell", "mhc", "complement",
            "fc receptor", "t cell receptor", "b cell receptor", "major histocompatibility",
        ],
        "keywords_zh": ["免疫", "抗原", "补体", "抗体", "淋巴细胞"],
    },
    "hypoxia": {
        "label": "缺氧",
        "description": "缺氧应答与 HIF 相关通路",
        "keywords_en": ["hypoxia", "hypoxic", "hif-1", "oxygen level", "anaerobic"],
        "keywords_zh": ["缺氧", "低氧"],
    },
    "oxidative_stress": {
        "label": "氧化应激",
        "description": "活性氧、抗氧化与过氧化相关通路",
        "keywords_en": [
            "oxidative stress", "reactive oxygen", "antioxidant", "peroxide",
            "glutathione", "cellular response to oxidative",
        ],
        "keywords_zh": ["氧化应激", "活性氧", "抗氧化"],
    },
    "apoptosis": {
        "label": "凋亡",
        "description": "程序性细胞死亡与 caspase 相关通路",
        "keywords_en": [
            "apopt", "cell death", "caspase", "bcl-2", "necroptosis", "ferroptosis",
        ],
        "keywords_zh": ["凋亡", "程序性死亡", "细胞死亡", "坏死"],
    },
    "emt": {
        "label": "EMT 与转移",
        "description": "上皮-间质转化、细胞黏附、细胞外基质与侵袭",
        "keywords_en": [
            "epithelial to mesenchymal", "epithelial mesenchymal", "cadherin",
            "extracellular matrix", "collagen", "focal adhesion",
            "cell adhesion", "migration", "invasion", "metastasis",
        ],
        "keywords_zh": ["上皮间质", "emt", "转移", "侵袭", "细胞黏附", "细胞外基质"],
    },
    "cell_cycle": {
        "label": "细胞周期",
        "description": "细胞周期、有丝分裂与 DNA 复制通路",
        "keywords_en": [
            "cell cycle", "mitotic", "mitosis", "dna replication",
            "chromosome segregation", "g1/s", "g2/m", "e2f", "cyclin",
        ],
        "keywords_zh": ["细胞周期", "有丝分裂", "增殖"],
    },
}

# 中文关键词只在用户输入中做主题归类，不直接匹配英文 term 名称；
# 真正用于 term 匹配的是主题的英文关键词加查询中的英文词。
_ASCII_TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9_-]{2,}")


def list_themes():
    """返回全部支持的主题（供 AI 工具与前端展示）。"""
    return [
        {
            "key": key,
            "label": spec["label"],
            "description": spec.get("description", ""),
            "keywords": list(dict.fromkeys(spec["keywords_zh"] + spec["keywords_en"])),
            "keywords_en": list(spec["keywords_en"]),
            "keywords_zh": list(spec["keywords_zh"]),
        }
        for key, spec in THEMES.items()
    ]


def resolve_themes(text):
    """把自然语言文本解析为命中的主题 key 列表（可命中多个）。"""
    lowered = str(text or "").lower()
    if not lowered.strip():
        return []
    matched = []
    for key, spec in THEMES.items():
        triggers = [item.lower() for item in spec["keywords_zh"] + spec["keywords_en"]]
        if any(trigger in lowered for trigger in triggers):
            matched.append(key)
    return matched


def _search_keywords(query):
    """返回用于匹配 term 名称的关键词集合。"""
    matched = resolve_themes(query)
    keywords = set()
    for key in matched:
        keywords.update(item.lower() for item in THEMES[key]["keywords_en"])
    # 保留查询中的英文词，让未登记主题的英文直查也能工作。
    keywords.update(token.lower() for token in _ASCII_TOKEN_RE.findall(str(query or "")))
    return matched, {item for item in keywords if len(item) >= 3}


def _managed_gene_set_dir():
    """与 sc_cell_go 相同的平台本地基因集目录。"""
    return (Path(Config.DATA_DIR) / "go_gene_sets").resolve()


def _resolve_libraries(libraries):
    if not libraries:
        return list(DEFAULT_LIBRARIES)
    if isinstance(libraries, str):
        libraries = [item.strip() for item in re.split(r"[,;\n]+", libraries) if item.strip()]
    resolved = []
    for library in libraries:
        library = str(library).strip()
        if library in {".", ".."} or not _LIBRARY_NAME_RE.fullmatch(library):
            raise ValueError(f"基因集名称 '{library}' 无效；只允许字母、数字、点、下划线和连字符。")
        if library not in resolved:
            resolved.append(library)
    return resolved


def load_library_terms(library):
    """读取一个平台管理的本地基因集快照。

    路径及符号链接校验复用 ``sc_cell_go`` 的实际富集实现，确保检索到的
    term 可直接用于同一运行的 ``focus_terms``。读取失败时由调用方决定
    是向用户报错还是记录为单库错误。
    """
    from modules.sc_cell_go import _local_gene_set_path, _read_local_gene_sets

    resolved = _resolve_libraries([library])
    path = _local_gene_set_path(resolved[0], str(_managed_gene_set_dir()))
    return _read_local_gene_sets(path)


def _match_score(term_lower, keywords):
    return sum(1 for keyword in keywords if keyword in term_lower)


def search_terms(query, libraries=None, limit=60):
    """在本地基因集快照中做确定性主题/关键词检索。

    返回 dict：themes（命中主题）、keywords（实际用于匹配的关键词）、
    terms（按关键词命中数排序的 term 列表，每项含 term、library、n_genes、
    example_genes、matched_keywords）与 library_errors（单个库缺失时不中断，
    逐库记录原因）。
    """
    query = str(query or "").strip()
    matched_themes, keywords = _search_keywords(query)
    requested = _resolve_libraries(libraries)
    capped_limit = max(1, min(int(limit if limit is not None else 60), 200))

    result_terms = []
    library_errors = {}
    if keywords:
        for library in requested:
            try:
                gene_sets = load_library_terms(library)
            except (FileNotFoundError, ValueError) as exc:
                library_errors[library] = str(exc)
                continue
            for term, genes in gene_sets.items():
                term_lower = term.lower()
                hits = [keyword for keyword in sorted(keywords) if keyword in term_lower]
                if hits:
                    result_terms.append({
                        "term": term,
                        "library": library,
                        "n_genes": len(genes),
                        "example_genes": genes[:6],
                        "matched_keywords": hits,
                        "_score": len(hits),
                    })
    result_terms.sort(key=lambda item: (-item["_score"], item["library"], item["term"]))
    for item in result_terms:
        item.pop("_score", None)
    truncated = len(result_terms) > capped_limit
    return {
        "query": query,
        "themes": matched_themes,
        "keywords": sorted(keywords),
        "libraries_searched": requested,
        "library_errors": library_errors,
        "n_terms": min(len(result_terms), capped_limit),
        "n_terms_total": len(result_terms),
        "truncated": truncated,
        "terms": result_terms[:capped_limit],
    }
