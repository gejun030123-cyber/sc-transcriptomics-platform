from .qc import QCAnalysis
from .preprocess import PreprocessAnalysis
from .dimred import DimredAnalysis
from .batch_correct import BatchCorrectAnalysis
from .clustering import ClusteringAnalysis
from .annotation import AnnotationAnalysis
from .deg import DEGAnalysis
from .trajectory import TrajectoryAnalysis
from .proportion import ProportionAnalysis

MODULE_REGISTRY = {
    'qc': QCAnalysis,
    'preprocess': PreprocessAnalysis,
    'dimred': DimredAnalysis,
    'batch_correct': BatchCorrectAnalysis,
    'clustering': ClusteringAnalysis,
    'annotation': AnnotationAnalysis,
    'deg': DEGAnalysis,
    'trajectory': TrajectoryAnalysis,
    'proportion': ProportionAnalysis,
}

PIPELINE_ORDER = [
    'qc', 'preprocess', 'dimred', 'batch_correct', 'clustering',
    'annotation', 'deg', 'trajectory', 'proportion'
]
