"""Phase 1 固定图型模板。"""

from .enrichment import NatureEnrichmentDotplot, NatureGSEA
from .enrichment_overview import NatureEnrichmentOverview
from .enrichment_extended import (
    NatureEnrichmentBarplot,
    NatureEnrichmentChord,
    NatureEnrichmentCnetplot,
    NatureEnrichmentEmapplot,
    NatureGSEARunning,
)
from .correlation import NatureCorrelationHeatmap
from .gene_sets import NatureGSVA, NatureSSGSEA
from .heatmap import NatureHeatmap
from .ma import NatureMA
from .pca import NaturePCA
from .upset import NatureUpSet
from .volcano import NatureVolcano
from .wgcna import NatureWGCNA
from .diagnostics import NatureDiagnostic

__all__ = [
    'NaturePCA', 'NatureVolcano', 'NatureHeatmap',
    'NatureGSEA', 'NatureEnrichmentDotplot', 'NatureEnrichmentOverview',
    'NatureEnrichmentBarplot', 'NatureEnrichmentChord',
    'NatureEnrichmentCnetplot', 'NatureEnrichmentEmapplot',
    'NatureGSEARunning',
    'NatureMA', 'NatureCorrelationHeatmap', 'NatureGSVA', 'NatureSSGSEA',
    'NatureUpSet', 'NatureWGCNA',
    'NatureDiagnostic',
]
