"""FigureSpec：科学参数和投稿输出的唯一绘图合同。"""

from dataclasses import asdict, dataclass, field, replace
from typing import Any, Mapping, Optional, Sequence, Tuple


VALID_MODES = frozenset({'standard', 'publication', 'nature_portfolio'})
VALID_WIDTHS = frozenset({'single', 'double'})
VALID_STYLES = frozenset({'nature', 'nature_communications', 'nature_aging'})
WIDTH_MM = {'single': 89.0, 'double': 183.0}


def _string_tuple(values: Optional[Sequence[Any]]) -> Tuple[str, ...]:
    if values is None:
        return ()
    if isinstance(values, str):
        values = values.replace('\n', ',').replace(';', ',').split(',')
    return tuple(dict.fromkeys(str(value).strip() for value in values if str(value).strip()))


@dataclass
class FigureSpec:
    """可序列化且与 renderer 解耦的图形说明。

    `fc_threshold` 始终表示绝对 log2 fold-change 阈值。适配旧页面时应先把普通
    Fold Change 转换为 log2，避免同一个字段在不同图中有不同单位。
    """

    plot_type: str
    mode: str = 'nature_portfolio'
    width: str = 'single'
    style: str = 'nature'
    title: str = ''
    evidence_role: str = 'discovery'
    height_mm: Optional[float] = None
    dpi: int = 600
    formats: Tuple[str, ...] = ('svg', 'pdf', 'png')

    # DEG / Volcano
    fc_threshold: float = 1.0
    fdr_threshold: float = 0.05
    label_strategy: str = 'top_significant'
    label_n: int = 8
    label_genes: Tuple[str, ...] = ()
    show_legend: bool = False

    # PCA
    confidence_ellipse: bool = False
    ellipse_level: float = 0.95
    show_sample_labels: bool = False
    outlier_labels: Tuple[str, ...] = ()

    # Heatmap
    zscore: str = 'row'
    row_cluster: bool = True
    col_cluster: bool = True
    cluster_method: str = 'average'
    distance_metric: str = 'euclidean'
    max_row_labels: int = 30
    max_col_labels: int = 14
    color_limit: Optional[float] = None

    # GSEA / ORA
    top_n: int = 12
    # Multi-database enrichment overview.  ``top_n`` is interpreted per
    # database panel for this plot type; the complete integrated table remains
    # available as source data.
    database_scope: Tuple[str, ...] = ()
    term_wrap: int = 34
    term_max_chars: int = 82
    max_genes: int = 30
    similarity_threshold: float = 0.15
    redundancy_threshold: float = 0.85
    running_term_n: int = 2
    # Enrichment target selection and membership-label controls.  These are
    # scientific display controls; the complete enrichment table is never
    # discarded when a figure is narrowed to user-selected terms.
    pathway_selection: str = 'top'
    pathway_terms: Tuple[str, ...] = ()
    gene_label_strategy: str = 'all'
    gene_labels: Tuple[str, ...] = ()
    max_gene_labels: int = 30

    # MA / correlation / gene-set scores / network summaries
    correlation_method: str = 'pearson'
    mask_diagonal: bool = True
    annotate_cells: bool = False
    top_intersections: int = 20
    score_scale: str = 'row'

    # 模板可以读取的非视觉科学选项；不允许借此覆盖 style rcParams。
    extra: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.plot_type = str(self.plot_type).strip().lower()
        self.mode = str(self.mode).strip().lower()
        self.width = str(self.width).strip().lower()
        self.style = str(self.style).strip().lower()
        if not self.plot_type:
            raise ValueError('FigureSpec.plot_type 不能为空')
        if self.mode not in VALID_MODES:
            raise ValueError(f'不支持的绘图模式: {self.mode}')
        if self.width not in VALID_WIDTHS:
            raise ValueError(f'不支持的投稿宽度: {self.width}')
        if self.style not in VALID_STYLES:
            raise ValueError(f'不支持的 Nature style: {self.style}')
        if self.height_mm is not None and not 35 <= float(self.height_mm) <= 230:
            raise ValueError('height_mm 必须在 35-230 mm 之间')
        if not 72 <= int(self.dpi) <= 1200:
            raise ValueError('dpi 必须在 72-1200 之间')
        if not 0 < float(self.fdr_threshold) <= 1:
            raise ValueError('fdr_threshold 必须在 (0, 1]')
        if float(self.fc_threshold) < 0:
            raise ValueError('fc_threshold 必须大于等于 0')
        if not 0 < float(self.ellipse_level) < 1:
            raise ValueError('ellipse_level 必须在 (0, 1)')
        if self.cluster_method not in {'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'}:
            raise ValueError(f'不支持的聚类方法: {self.cluster_method}')
        self.label_n = max(0, min(30, int(self.label_n)))
        self.top_n = max(1, min(60, int(self.top_n)))
        self.max_genes = max(4, min(80, int(self.max_genes)))
        self.running_term_n = max(1, min(4, int(self.running_term_n)))
        self.pathway_selection = str(self.pathway_selection).strip().lower()
        if self.pathway_selection not in {'top', 'selected', 'selected_plus_top'}:
            raise ValueError('pathway_selection 必须为 top、selected 或 selected_plus_top')
        self.gene_label_strategy = str(self.gene_label_strategy).strip().lower()
        if self.gene_label_strategy not in {'all', 'shared', 'selected', 'none'}:
            raise ValueError('gene_label_strategy 必须为 all、shared、selected 或 none')
        self.pathway_terms = _string_tuple(self.pathway_terms)
        self.database_scope = _string_tuple(self.database_scope)
        self.gene_labels = _string_tuple(self.gene_labels)
        self.max_gene_labels = max(0, min(80, int(self.max_gene_labels)))
        self.similarity_threshold = float(self.similarity_threshold)
        if not 0 < self.similarity_threshold <= 1:
            raise ValueError('similarity_threshold 必须在 (0, 1]')
        self.redundancy_threshold = float(self.redundancy_threshold)
        if not 0 < self.redundancy_threshold <= 1:
            raise ValueError('redundancy_threshold 必须在 (0, 1]')
        self.top_intersections = max(1, min(40, int(self.top_intersections)))
        self.correlation_method = str(self.correlation_method).strip().lower()
        if self.correlation_method not in {'pearson', 'spearman'}:
            raise ValueError(f'不支持的相关性方法: {self.correlation_method}')
        self.score_scale = str(self.score_scale).strip().lower()
        if self.score_scale not in {'row', 'none', 'precomputed'}:
            raise ValueError(f'不支持的 gene-set score scale: {self.score_scale}')
        self.max_row_labels = max(0, int(self.max_row_labels))
        self.max_col_labels = max(0, int(self.max_col_labels))
        self.formats = tuple(dict.fromkeys(
            str(value).lower() for value in self.formats
            if str(value).lower() in {'svg', 'pdf', 'png', 'tiff'}
        ))
        if not self.formats:
            raise ValueError('至少需要一个有效导出格式')
        self.label_genes = _string_tuple(self.label_genes)
        self.outlier_labels = _string_tuple(self.outlier_labels)
        self.extra = dict(self.extra or {})

    @property
    def width_mm(self) -> float:
        return WIDTH_MM[self.width]

    @property
    def width_inches(self) -> float:
        return self.width_mm / 25.4

    def with_updates(self, **updates: Any) -> 'FigureSpec':
        return replace(self, **updates)

    def to_dict(self) -> dict:
        payload = asdict(self)
        payload['width_mm'] = self.width_mm
        payload['label_genes'] = list(self.label_genes)
        payload['outlier_labels'] = list(self.outlier_labels)
        payload['pathway_terms'] = list(self.pathway_terms)
        payload['database_scope'] = list(self.database_scope)
        payload['gene_labels'] = list(self.gene_labels)
        payload['formats'] = list(self.formats)
        return payload

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> 'FigureSpec':
        values = dict(payload)
        values.pop('width_mm', None)
        return cls(**values)
