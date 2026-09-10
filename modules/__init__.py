from .qc import QCAnalysis
from .normalize import NormalizeAnalysis
from .hvg import HVGAnalysis
from .dimred import DimredAnalysis
from .batch_correct import BatchCorrectAnalysis
from .clustering import ClusteringAnalysis
from .subcluster import SubclusterAnalysis
from .qc_reassess import QCReassessAnalysis
from .annotation import AnnotationAnalysis
from .functional_state import FunctionalStateAnalysis
from .scenic import ScenicAnalysis
from .deg import DEGAnalysis
from .trajectory import TrajectoryAnalysis
from .sc_timecourse import SCTimecourseAnalysis
from .proportion import ProportionAnalysis
from .neighborhood_da import NeighborhoodDAAnalysis
from .cell_communication import CellCommunicationAnalysis
from .virtual_ko import VirtualKOAnalysis
from .sc_batch import SCBatchImport
from .sc_batch_export import SCBatchCSVExport, SCPseudobulkDEG
from .sc_cell_deg import SCCellLevelDEG
from .sc_cell_go import SCCellGOEnrichment
from .bulk_qc import BulkQCAnalysis
from .bulk_normalize import BulkNormalizeAnalysis
from .bulk_deg import BulkDEGAnalysis
from .bulk_pca import BulkPCAAnalysis
from .bulk_heatmap import BulkHeatmapAnalysis
from .bulk_enrichment import BulkEnrichmentAnalysis
from .bulk_timecourse import BulkTimecourseAnalysis
from .bulk_deg_integration import BulkDEGIntegrationAnalysis
from .convert_10x import Convert10x
from .schemas import SC_MODULE_NAMES, BULK_MODULE_NAMES

MODULE_REGISTRY = {
    # 单细胞分析模块
    'qc': QCAnalysis,
    'normalize': NormalizeAnalysis,
    'hvg': HVGAnalysis,
    'dimred': DimredAnalysis,
    'batch_correct': BatchCorrectAnalysis,
    'clustering': ClusteringAnalysis,
    'subcluster': SubclusterAnalysis,
    'qc_reassess': QCReassessAnalysis,
    'annotation': AnnotationAnalysis,
    'functional_state': FunctionalStateAnalysis,
    'scenic': ScenicAnalysis,
    'deg': DEGAnalysis,
    'trajectory': TrajectoryAnalysis,
    'sc_timecourse': SCTimecourseAnalysis,
    'proportion': ProportionAnalysis,
    'neighborhood_da': NeighborhoodDAAnalysis,
    'cell_communication': CellCommunicationAnalysis,
    'virtual_ko': VirtualKOAnalysis,
    'sc_batch_import': SCBatchImport,
    'sc_cell_deg': SCCellLevelDEG,
    'sc_cell_go': SCCellGOEnrichment,
    'sc_pseudobulk_deg': SCPseudobulkDEG,
    'sc_csv_export': SCBatchCSVExport,
    # Bulk RNA-seq 分析模块
    'bulk_qc': BulkQCAnalysis,
    'bulk_normalize': BulkNormalizeAnalysis,
    'bulk_deg': BulkDEGAnalysis,
    'bulk_pca': BulkPCAAnalysis,
    'bulk_heatmap': BulkHeatmapAnalysis,
    'bulk_enrichment': BulkEnrichmentAnalysis,
    'bulk_timecourse': BulkTimecourseAnalysis,
    'bulk_deg_integration': BulkDEGIntegrationAnalysis,
    # 数据导入模块
    'convert_10x': Convert10x,
}

PIPELINE_ORDER = [
    'qc', 'normalize', 'hvg', 'dimred', 'batch_correct', 'clustering', 'subcluster',
    'qc_reassess', 'annotation', 'functional_state', 'scenic', 'sc_timecourse', 'deg', 'sc_cell_deg', 'sc_pseudobulk_deg', 'sc_cell_go', 'trajectory', 'proportion', 'neighborhood_da', 'cell_communication', 'virtual_ko', 'sc_csv_export',
    'bulk_qc', 'bulk_normalize', 'bulk_pca', 'bulk_deg', 'bulk_heatmap', 'bulk_enrichment', 'bulk_timecourse',
    'bulk_deg_integration',
]

# 模块依赖约束：value 中的模块必须在 key 之前执行
PIPELINE_DEPS = {
    # 单细胞
    'normalize': ['qc'],
    'hvg': ['normalize'],
    'dimred': ['hvg'],
    'batch_correct': ['hvg'],
    'clustering': ['dimred'],
    'subcluster': ['clustering'],
    'qc_reassess': ['clustering'],
    'annotation': ['clustering'],
    'functional_state': ['annotation'],
    'scenic': ['annotation'],
    'sc_timecourse': ['clustering'],
    'deg': ['clustering'],
    'trajectory': ['clustering'],
    'proportion': ['clustering'],
    'neighborhood_da': ['clustering'],
    'cell_communication': ['annotation'],
    'virtual_ko': ['clustering'],
    'sc_pseudobulk_deg': ['qc'],
    'sc_cell_deg': ['clustering'],
    'sc_csv_export': ['qc'],
    # Bulk RNA-seq
    'bulk_normalize': ['bulk_qc'],
    'bulk_pca': ['bulk_normalize'],
    'bulk_deg': ['bulk_normalize'],
    'bulk_heatmap': ['bulk_deg'],
    'bulk_enrichment': ['bulk_deg'],
    'bulk_timecourse': ['bulk_normalize'],
    'bulk_deg_integration': ['bulk_deg'],
}

def validate_pipeline_order(modules):
    """校验模块执行顺序是否满足依赖约束。
    返回 (is_valid: bool, errors: list[str])。
    """
    errors = []
    seen = set()
    for mod in modules:
        deps = PIPELINE_DEPS.get(mod, [])
        for dep in deps:
            if dep not in seen and dep in modules:
                errors.append(f"'{mod}' 依赖 '{dep}'，但 '{dep}' 未在其之前执行")
        seen.add(mod)
    return len(errors) == 0, errors
