# modules/evaluators/__init__.py
"""评估器模块 — 对分析候选结果进行确定性评分"""

from .sc_cluster import (
    score_cluster_signature,
    score_cell_type_signature,
)
