from .qc import QCAnalysis
from .preprocess import PreprocessAnalysis
from .dimred import DimredAnalysis
from .batch_correct import BatchCorrectAnalysis
from .clustering import ClusteringAnalysis
from .annotation import AnnotationAnalysis
from .deg import DEGAnalysis
from .trajectory import TrajectoryAnalysis
from .proportion import ProportionAnalysis
from .bulk_qc import BulkQCAnalysis
from .bulk_normalize import BulkNormalizeAnalysis
from .bulk_deg import BulkDEGAnalysis
from .bulk_pca import BulkPCAAnalysis
from .bulk_heatmap import BulkHeatmapAnalysis

MODULE_REGISTRY = {
    # 单细胞分析模块
    'qc': QCAnalysis,
    'preprocess': PreprocessAnalysis,
    'dimred': DimredAnalysis,
    'batch_correct': BatchCorrectAnalysis,
    'clustering': ClusteringAnalysis,
    'annotation': AnnotationAnalysis,
    'deg': DEGAnalysis,
    'trajectory': TrajectoryAnalysis,
    'proportion': ProportionAnalysis,
    # Bulk RNA-seq 分析模块
    'bulk_qc': BulkQCAnalysis,
    'bulk_normalize': BulkNormalizeAnalysis,
    'bulk_deg': BulkDEGAnalysis,
    'bulk_pca': BulkPCAAnalysis,
    'bulk_heatmap': BulkHeatmapAnalysis,
}

PIPELINE_ORDER = [
    'qc', 'preprocess', 'dimred', 'batch_correct', 'clustering',
    'annotation', 'deg', 'trajectory', 'proportion',
    'bulk_qc', 'bulk_normalize', 'bulk_deg', 'bulk_pca', 'bulk_heatmap',
]
