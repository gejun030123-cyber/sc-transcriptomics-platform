import json
import os

from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping
from modules.llm_annotation import LLMAnnotationError, run_llm_cluster_annotation
from modules.sc_de_utils import _matrix_is_raw_counts, log1p_expression_from_counts


# CellTypist is intentionally optional.  The downloaded reference models live
# outside project outputs so they can be shared by all projects without
# consuming the system temporary directory.  The environment override keeps
# deployments with a separate data volume reproducible.
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CELLTYPIST_MODEL_DIR = os.environ.get(
    'SC_CELLTYPIST_MODEL_DIR',
    os.path.join(_PROJECT_ROOT, 'data', 'references', 'celltypist'),
)
CELLTYPIST_MODEL_OPTIONS = [
    'Immune_All_Low.pkl',
    'Immune_All_High.pkl',
    'Cells_Intestinal_Tract.pkl',
    'Developing_Human_Organs.pkl',
    'Developing_Human_Brain.pkl',
    'Cells_Fetal_Lung.pkl',
    'Cells_Lung_Airway.pkl',
    'Human_Lung_Atlas.pkl',
    'Nuclei_Lung_Airway.pkl',
    'Healthy_Human_Liver.pkl',
    'Adult_Human_PancreaticIslet.pkl',
    'Fetal_Human_Pancreas.pkl',
    'Healthy_Adult_Heart.pkl',
    'Adult_Human_Skin.pkl',
    'Fetal_Human_Skin.pkl',
    'Adult_Human_Vascular.pkl',
    'Cells_Adult_Breast.pkl',
    'Pan_Fetal_Human.pkl',
]


def resolve_celltypist_model(model_name=None, model_dir=None):
    """Resolve a downloaded CellTypist model inside the managed model folder.

    CellTypist models are pickle files.  A request-controlled absolute path,
    traversal component, or symlink would therefore turn annotation into an
    arbitrary pickle loader.  Custom models remain supported by placing a
    basename-only ``.pkl`` file in the administrator-configured directory.
    """
    requested = str(model_name or 'Immune_All_Low.pkl').strip()
    if (
        not requested
        or requested != os.path.basename(requested)
        or '/' in requested
        or '\\' in requested
        or not requested.lower().endswith('.pkl')
    ):
        raise ValueError('CellTypist 模型只允许填写托管模型目录内的 .pkl 文件名。')
    model_root = os.path.realpath(str(model_dir or CELLTYPIST_MODEL_DIR))
    candidate = os.path.join(model_root, requested)
    if os.path.islink(candidate):
        raise ValueError('CellTypist 模型不能是符号链接。')
    if not os.path.isfile(candidate):
        raise FileNotFoundError(
            f"CellTypist 模型不存在: {candidate}。请先将选定的 .pkl 放入 {model_root}。"
        )
    resolved = os.path.realpath(candidate)
    try:
        contained = os.path.commonpath([model_root, resolved]) == model_root
    except ValueError:
        contained = False
    if not contained:
        raise ValueError('CellTypist 模型必须位于托管模型目录内。')
    return resolved


def _annotation_expression_matrix(adata, target_sum=10_000.0):
    """Return a non-negative annotation matrix without mutating ``adata``.

    Marker scoring and marker ranking require an expression scale, not signed
    Pearson residuals or z-scores.  Prefer protected counts, then a compatible
    ``raw`` matrix, and only finally the current ``X``.  Integer-valued sources
    are treated as raw counts and normalized exactly once.
    """
    import numpy as np

    if 'counts' in adata.layers:
        matrix = adata.layers['counts']
        source = 'layers["counts"]'
    elif adata.raw is not None:
        raw_names = {str(gene) for gene in adata.raw.var_names}
        if all(str(gene) in raw_names for gene in adata.var_names):
            matrix = adata.raw[:, list(adata.var_names)].X
            source = 'adata.raw.X'
        else:
            matrix = adata.X
            source = 'adata.X'
    else:
        matrix = adata.X
        source = 'adata.X'

    values = matrix.data if hasattr(matrix, 'tocsr') else np.asarray(matrix)
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size and not np.isfinite(values).all():
        raise ValueError(f'{source} 含 NaN/Inf，不能用于细胞注释。')
    if values.size and (values < 0).any():
        raise ValueError(
            f'{source} 含负值，且没有可恢复的非负 counts/raw 表达矩阵；'
            '不能在 Pearson residuals 或 scaled X 上进行 Marker 注释。'
        )
    if _matrix_is_raw_counts(matrix):
        return (
            log1p_expression_from_counts(matrix, target_sum=target_sum),
            f'{source} → library-size normalized log1p',
        )
    copied = matrix.copy() if hasattr(matrix, 'copy') else np.asarray(matrix, dtype=float).copy()
    return copied, f'{source}（assumed normalized non-negative expression）'


def _celltypist_label_family(label):
    """Map labels to conservative broad families for cross-method review."""
    value = str(label or '').strip().lower()
    if value in {'', 'unknown', 'unassigned', 'heterogeneous', 'nan'}:
        return set()
    families = {
        'immune': ('immune', 't cell', 't-cell', 'b cell', 'b-cell', 'nk',
                   'monocyte', 'macrophage', 'dendritic', 'neutrophil',
                   'mast', 'plasma', 'lymph', 'myeloid', 'thymocyte'),
        'epithelial': ('epithelial', 'enterocyte', 'goblet', 'paneth', 'tuft',
                       'ductal', 'acinar', 'hepatocyte', 'cholangiocyte',
                       'alveolar', 'club', 'ciliated', 'basal', 'beta cell',
                       'alpha cell', 'cardiomyocyte'),
        'endothelial': ('endothelial', 'vascular', 'capillary'),
        'mesenchymal': ('fibroblast', 'stromal', 'pericyte', 'smooth muscle',
                        'mesenchymal', 'stellate'),
        'neural': ('neuron', 'neural', 'astro', 'oligo', 'microglia',
                   'radial glia', 'progenitor', 'ependymal', 'choroid'),
        'cycling': ('cycling', 'proliferating', 'progenitor'),
    }
    return {
        family for family, tokens in families.items()
        if any(token in value for token in tokens)
    }


def _celltypist_comparison_status(marker_label, reference_label):
    """Return an auditable comparison status without changing final labels."""
    marker = str(marker_label or '').strip()
    reference = str(reference_label or '').strip()
    if reference.lower() in {'', 'unknown', 'unassigned', 'heterogeneous', 'nan'}:
        return 'unassigned'
    if '|' in reference:
        return 'multi_label'
    if marker.lower() in {'', 'unknown', 'nan'}:
        return 'marker_unknown'
    if marker.lower() == reference.lower():
        return 'agree'
    if _celltypist_label_family(marker) & _celltypist_label_family(reference):
        return 'lineage_agree'
    return 'conflict'


def run_celltypist_reference(
    adata,
    model_name=None,
    model_dir=None,
    mode='prob match',
    p_thres=0.5,
    majority_voting=False,
    over_clustering=None,
):
    """Run a local CellTypist model as reference evidence.

    The function is deliberately side-effect free with respect to the input
    AnnData.  It creates a temporary query matrix, preferring raw counts from
    ``layers['counts']`` and normalising it for CellTypist.  It returns labels
    and probabilities; the caller decides how to combine them with Marker
    evidence.  Missing optional dependencies or models are surfaced as clear
    exceptions so the main Marker path can gracefully fall back.
    """
    import importlib
    import numpy as np
    import pandas as pd
    import sys

    requested_model_dir = str(model_dir or CELLTYPIST_MODEL_DIR)
    # CellTypist creates its own cache on import.  Point that cache to the
    # project data volume before importing it; otherwise read-only or quota-
    # limited home directories can make an otherwise valid local model fail.
    os.environ.setdefault(
        'CELLTYPIST_FOLDER',
        os.path.join(requested_model_dir, '.celltypist_runtime'),
    )
    os.makedirs(os.environ['CELLTYPIST_FOLDER'], exist_ok=True)
    # A direct helper call can happen before the application Config module has
    # initialised Scanpy/Numba caches.  Keep those optional compilation files
    # on the project data volume as well, avoiding read-only home failures.
    os.environ.setdefault(
        'NUMBA_CACHE_DIR',
        os.path.join(os.environ['CELLTYPIST_FOLDER'], 'numba_cache'),
    )
    os.makedirs(os.environ['NUMBA_CACHE_DIR'], exist_ok=True)
    if 'scanpy' not in sys.modules:
        os.environ.setdefault('NUMBA_DISABLE_JIT', '1')
    model_path = resolve_celltypist_model(model_name, requested_model_dir)
    celltypist = importlib.import_module('celltypist')

    # Avoid copying the full AnnData object (including large obsm/uns fields).
    # CellTypist only needs a cell-by-gene matrix and gene symbols.
    from anndata import AnnData
    matrix, expression_source = _annotation_expression_matrix(adata)
    query = AnnData(
        X=matrix,
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        var=pd.DataFrame(index=adata.var_names.copy()),
    )

    requested_mode = str(mode or 'prob match').strip().lower()
    if requested_mode in {'prob match', 'prob_match', 'probability'}:
        mode = 'prob match'
    elif requested_mode in {'best match', 'best_match'}:
        mode = 'best match'
    else:
        raise ValueError("CellTypist mode 必须为 'prob match' 或 'best match'。")
    p_thres = float(p_thres)
    if not np.isfinite(p_thres) or not 0 <= p_thres <= 1:
        raise ValueError('CellTypist p_thres 必须位于 [0, 1]。')
    kwargs = {
        'model': model_path,
        'mode': mode,
        'p_thres': p_thres,
        'majority_voting': bool(majority_voting),
    }
    if majority_voting and over_clustering is not None:
        kwargs['over_clustering'] = over_clustering
    prediction = celltypist.annotate(query, **kwargs)

    labels_table = getattr(prediction, 'predicted_labels', None)
    if labels_table is None:
        raise RuntimeError('CellTypist 返回结果缺少 predicted_labels。')
    if hasattr(labels_table, 'columns'):
        preferred_label = 'majority_voting' if majority_voting and 'majority_voting' in labels_table.columns else 'predicted_labels'
        if preferred_label in labels_table.columns:
            labels = labels_table[preferred_label].astype(str).to_numpy()
        else:
            labels = labels_table.iloc[:, 0].astype(str).to_numpy()
        if 'conf_score' in labels_table.columns:
            confidence = pd.to_numeric(labels_table['conf_score'], errors='coerce').to_numpy(dtype=float)
        else:
            confidence = None
    else:
        labels = np.asarray(labels_table).astype(str).reshape(-1)
        confidence = None

    if len(labels) != adata.n_obs:
        raise RuntimeError(
            f'CellTypist 返回 {len(labels)} 个标签，但输入有 {adata.n_obs} 个细胞。'
        )

    probability_matrix = getattr(prediction, 'probability_matrix', None)
    probability = None
    top2 = np.repeat('', adata.n_obs).astype(object)
    if probability_matrix is not None:
        probability_values = np.asarray(probability_matrix, dtype=float)
        if probability_values.ndim == 2 and probability_values.shape[0] == adata.n_obs:
            probability = np.nanmax(probability_values, axis=1)
            if confidence is None:
                confidence = probability
            if probability_values.shape[1] > 1 and hasattr(probability_matrix, 'columns'):
                order = np.argsort(np.nan_to_num(probability_values, nan=-np.inf), axis=1)
                top2 = np.asarray(probability_matrix.columns)[order[:, -2]].astype(object)
    if confidence is None:
        confidence = np.full(adata.n_obs, np.nan, dtype=float)
    confidence = np.asarray(confidence, dtype=float)
    labels = np.where(
        np.isin(np.char.lower(labels.astype(str)), ['unassigned', 'heterogeneous', 'nan']),
        'Unknown', labels,
    ).astype(object)
    status = np.where(
        labels == 'Unknown', 'unassigned',
        np.where(np.char.find(labels.astype(str), '|') >= 0, 'multi_label',
                 np.where(np.isfinite(confidence) & (confidence < p_thres), 'low_confidence', 'reference')),
    ).astype(object)
    return {
        'labels': labels,
        'top2': top2,
        'confidence': confidence,
        'status': status,
        'probability_max': probability if probability is not None else confidence,
        'model_path': model_path,
        'model_name': os.path.basename(model_path),
        'mode': mode,
        'p_thres': p_thres,
        'majority_voting': bool(majority_voting),
        'expression_source': expression_source,
        'package_version': str(getattr(celltypist, '__version__', 'unknown')),
    }

DEFAULT_TME_MARKERS = {
    'Tumor Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
    'CAF': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'FAP', 'ACTA2'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'ENG', 'CLDN5'],
    'T cells': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'TRAC'],
    'NK cells': ['NKG7', 'GNLY', 'KLRD1', 'NCAM1', 'PRF1'],
    'B cells': ['CD79A', 'CD79B', 'MS4A1', 'CD19', 'PAX5'],
    'Plasma cells': ['JCHAIN', 'MZB1', 'SDC1', 'IGHG1', 'IGKC'],
    'Monocyte/Macrophage': ['CD14', 'CD68', 'CSF1R', 'LYZ', 'S100A8', 'S100A9'],
    'Dendritic cells': ['FCER1A', 'CD1C', 'CLEC10A', 'ITGAX', 'HLA-DRA'],
    'Neutrophils': ['CSF3R', 'CXCR2', 'FCGR3B', 'S100A12'],
    'Mast cells': ['KIT', 'TPSAB1', 'TPSB2', 'HDC', 'MS4A2'],
    'Pericytes': ['RGS5', 'PDGFRB', 'NOTCH3', 'MCAM', 'ACTA2'],
    'Proliferating': ['MKI67', 'TOP2A', 'PCNA', 'STMN1', 'CDK1'],
}

DEFAULT_IMMUNE_MARKERS = {
    'T cells': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'TRAC'],
    'CD4+ T': ['CD4', 'IL7R', 'TRBC2'],
    'CD8+ T': ['CD8A', 'CD8B', 'GZMK', 'GZMA', 'CCL5'],
    'T naive': ['LEF1', 'CCR7', 'TCF7'],
    'NK cells': ['NKG7', 'GNLY', 'KLRD1', 'NCAM1', 'PRF1'],
    'B cells': ['CD79A', 'CD79B', 'MS4A1', 'CD19', 'PAX5'],
    'Plasma cells': ['JCHAIN', 'MZB1', 'SDC1', 'IGHG1', 'IGKC'],
    'Monocyte/Macrophage': ['CD14', 'CD68', 'CSF1R', 'LYZ', 'S100A8', 'S100A9'],
    'Dendritic cells': ['FCER1A', 'CD1C', 'CLEC10A', 'ITGAX', 'HLA-DRA'],
    'Neutrophils': ['CSF3R', 'CXCR2', 'FCGR3B', 'S100A12'],
    'pDC': ['GZMB', 'IL3RA', 'COBLL1', 'TCF4'],
}

DEFAULT_BLOOD_MARKERS = {
    'HSC': ['CD34', 'CD38', 'KIT', 'THY1', 'CRHBP'],
    'Erythroid': ['HBA1', 'HBA2', 'HBB', 'GYPA', 'SLC4A1'],
    'Megakaryocyte': ['PF4', 'GP9', 'ITGA2B', 'VWF', 'GP1BA'],
    'Monocyte': ['CD14', 'LYZ', 'S100A8', 'S100A9', 'VCAN'],
    'Neutrophil': ['FCGR3B', 'CSF3R', 'CXCR2', 'S100A12', 'MPO'],
    'Eosinophil': ['SIGLEC8', 'IL5RA', 'CCR3', 'EPX', 'PRG2'],
    'Basophil': ['HDC', 'MS4A2', 'KIT', 'FCER1A', 'CPA3'],
    'B cell': ['CD79A', 'MS4A1', 'CD19', 'PAX5', 'CD79B'],
    'T cell': ['CD3D', 'CD3E', 'CD2', 'TRAC', 'CD3G'],
    'NK cell': ['NKG7', 'GNLY', 'KLRD1', 'NCAM1', 'PRF1'],
    'Dendritic cell': ['FCER1A', 'CD1C', 'CLEC10A', 'ITGAX', 'HLA-DRA'],
    'pDC': ['GZMB', 'IL3RA', 'COBLL1', 'TCF4', 'IRF7'],
}

DEFAULT_PBMC_MARKERS = {
    'CD4 Naive T cells': ['IL7R', 'LTB', 'CCR7', 'TCF7', 'MALAT1', 'LEF1'],
    'CD4 Memory T cells': ['IL7R', 'LTB', 'MALAT1', 'IL32', 'AQP3', 'GPR183'],
    'CD14+ Monocytes': ['CD14', 'LYZ', 'S100A8', 'S100A9', 'LGALS3', 'FCN1'],
    'B cells': ['MS4A1', 'CD79A', 'CD79B', 'CD74', 'HLA-DRA'],
    'CD8 T cells': ['CD8A', 'CD8B', 'CCL5', 'GZMK', 'GZMA', 'TRBC2'],
    'Cytotoxic T cells': ['NKG7', 'CCL5', 'GZMB', 'PRF1', 'CTSW', 'CD3D'],
    'NK cells': ['GNLY', 'NKG7', 'KLRD1', 'PRF1', 'CTSW'],
    'FCGR3A+ Monocytes': ['FCGR3A', 'MS4A7', 'LST1', 'FCER1G', 'AIF1', 'IFITM3'],
    'Conventional DC': ['FCER1A', 'CST3', 'CD1C', 'CLEC10A', 'HLA-DPA1', 'HLA-DPB1'],
    'Plasmacytoid DC': ['GZMB', 'IRF7', 'TCF4', 'IL3RA', 'SERPINF1'],
    'Megakaryocytes': ['PPBP', 'PF4', 'SDPR', 'GNG11', 'NRGN'],
}

# A deliberately broad first-pass panel. It is suitable when tissue type is
# unknown; users should use a tissue/immune panel or custom markers afterwards
# to refine a lineage, rather than treating these labels as final subtypes.
DEFAULT_UNIVERSAL_MARKERS = {
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'KRT7', 'CDH1'],
    'Endothelial': ['PECAM1', 'VWF', 'KDR', 'EMCN', 'CLDN5'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'COL3A1'],
    'Pericyte/Smooth muscle': ['RGS5', 'PDGFRB', 'CSPG4', 'MCAM', 'ACTA2'],
    'Myeloid': ['LYZ', 'TYROBP', 'LST1', 'FCER1G', 'AIF1'],
    'T cells': ['CD3D', 'CD3E', 'TRAC', 'CD247', 'LCK'],
    'NK cells': ['NKG7', 'KLRD1', 'GNLY', 'PRF1', 'TRBC2'],
    'B cells': ['MS4A1', 'CD79A', 'CD74', 'HLA-DRA', 'CD37'],
    'Plasma cells': ['JCHAIN', 'MZB1', 'SDC1', 'DERL3', 'IGKC'],
    'Mast cells': ['TPSAB1', 'TPSB2', 'KIT', 'MS4A2', 'HDC'],
    'Cycling cells': ['MKI67', 'TOP2A', 'STMN1', 'TYMS', 'CDK1'],
}

# Colorectal organoids need a hierarchical panel: broad non-epithelial
# lineages must remain visible, while epithelial clusters need a second pass
# that can resolve biologically meaningful intestinal subtypes.  Keeping the
# two levels explicit prevents a goblet signature from competing directly
# with a broad ``Epithelial`` signature and prevents immune/stromal cells from
# being forced into the nearest epithelial subtype.
COLORECTAL_LINEAGE_MARKERS = {
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'KRT7', 'CDH1'],
    'T cells': ['CD3D', 'CD3E', 'TRAC', 'CD247', 'LCK'],
    'Myeloid': ['LST1', 'TYROBP', 'FCER1G', 'CSF1R', 'AIF1', 'STAB1', 'MS4A6A'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'COL6A2'],
    'Endothelial': ['PECAM1', 'VWF', 'EMCN', 'KDR', 'CLDN5', 'MMRN1'],
    'Neural-like (review)': ['TUBB3', 'ELAVL3', 'ELAVL4', 'SNAP25', 'RIMS2',
                             'RIMBP2', 'CACNA1A', 'RBFOX3'],
}

COLORECTAL_EPITHELIAL_SUBTYPE_MARKERS = {
    'Stem/crypt-like epithelial': ['LGR5', 'ASCL2', 'OLFM4', 'SMOC2', 'SOX9'],
    'TA/S-phase epithelial': ['TYMS', 'PCLAF', 'PCNA', 'MCM5', 'MCM6', 'FEN1', 'GINS2'],
    'TA/G2M epithelial': ['MKI67', 'TOP2A', 'CENPF', 'UBE2C', 'BIRC5', 'ASPM', 'PLK1'],
    'Goblet-like epithelial': ['FCGBP', 'CLCA1', 'MUC2', 'SPINK4', 'TFF3', 'AGR2', 'ITLN1'],
    'Absorptive/enterocyte-like epithelial': ['ALPI', 'VIL1', 'FABP1', 'SI', 'KRT20',
                                               'CA1', 'CA2', 'CEACAM7'],
    'BEST4+ absorptive epithelial': ['BEST4', 'OTOP2', 'CA7', 'GUCA2A', 'SPIB'],
    'Enteroendocrine-like epithelial': ['CHGA', 'CHGB', 'NEUROD1', 'PAX6', 'PCSK1'],
    'Paneth/LYZ+ secretory epithelial': ['LYZ', 'MMP7', 'DEFA5', 'DEFA6', 'REG1A'],
    'Inflammatory epithelial': ['LCN2', 'DUOX2', 'NOS2', 'REG1A', 'REG3A', 'CXCL8'],
    'Regenerative/stress epithelial': ['KRT19', 'TFF2', 'CLDN4', 'TACSTD2', 'ANXA1',
                                        'KRT17', 'GDF15', 'CA9'],
}

DEFAULT_COLORECTAL_MARKERS = {
    **COLORECTAL_LINEAGE_MARKERS,
    **COLORECTAL_EPITHELIAL_SUBTYPE_MARKERS,
}

# 结直肠上皮的本地细分面板。这里的标签刻意对应可重复识别的 marker
# programme，而不是把多个不同 cluster 收进 “High metabolic” 之类的宽标签。
# 它适用于同一模型或相近的肠上皮数据；低证据 programme 会显式进入 review，
# 不能被解释为跨模型的参考图谱真值。
COLORECTAL_REFINED_MARKERS = {
    'CYP3A5+ enterocytes': ['CYP3A5', 'PTPRR', 'FAM13A', 'CADPS', 'ATXN1'],
    'KRT20/MALRD1+ absorptive enterocytes': ['KRT20', 'MALRD1', 'NR1H4', 'CA9', 'TM4SF4'],
    'HSD17B2+ absorptive enterocytes': ['HSD17B2', 'C11orf86', 'SDCBP2', 'PCK1', 'LGALS1'],
    'MTTP/RBP2+ absorptive enterocytes': ['MTTP', 'RBP2', 'GSTA1', 'MAF', 'TM4SF4'],
    'CKB/FABP1+ metabolic enterocytes': ['CKB', 'FABP1', 'FXYD3', 'PHGR1', 'S100A6'],
    'NDRG1/ANKRD37+ hypoxia-response enterocytes': ['NDRG1', 'ANKRD37', 'EFNA1', 'MXD1', 'PNRC1'],
    'MAML3/CHRM3+ epithelial cells (review)': ['MAML3', 'CHRM3', 'PGAP1', 'LINC01811', 'PVT1'],
    'TFF1/REG4+ secretory cells': ['TFF1', 'REG4', 'TFF3', 'SLPI', 'AGR2'],
    'MUC2/FCGBP+ Goblet cells': ['MUC2', 'FCGBP', 'CLCA1', 'GALNT8', 'TFF3'],
    'SPINK4/CA4+ secretory cells': ['SPINK4', 'CA4', 'CCL25', 'GSTA1', 'CLCA1'],
    'DDIT3/ATF3+ ER-stress enterocytes': ['DDIT3', 'ATF3', 'PPP1R15A', 'GADD45B', 'XBP1', 'HSPA5'],
    'TA/stem-like cells': ['CLDN2', 'ASCL2', 'LGR5', 'OLFM4', 'SOX9', 'RAMP1'],
    'S-phase enterocytes': ['FEN1', 'DTL', 'MCM10', 'GINS2', 'BRIP1', 'TYMS', 'PCNA', 'MCM5'],
    'MIR924HG/DIAPH3+ cycling cells': ['MIR924HG', 'DIAPH3', 'POLQ', 'CENPP', 'LINC01572'],
    'Histone-rich cycling enterocytes': ['HIST1H1B', 'HIST1H2AJ', 'HIST2H2AC', 'HIST1H2AH', 'HIST1H3B'],
    'BIRC5/CDKN3+ cycling enterocytes': ['BIRC5', 'CDKN3', 'UBE2C', 'PCLAF', 'HIST1H1D'],
    'CDC20/CCNB1+ cycling enterocytes': ['CDC20', 'CCNB1', 'DLGAP5', 'CENPF', 'ASPM', 'TOP2A'],
}

# One highly specific, independently ranked marker can be more informative
# than a shared score-gene module.  The anchor rule is intentionally limited
# to local, cluster-specific signatures and only applies to the refined panel.
COLORECTAL_REFINED_ANCHOR_MARKERS = {
    'CYP3A5+ enterocytes': {'CYP3A5'},
    'KRT20/MALRD1+ absorptive enterocytes': {'MALRD1'},
    'HSD17B2+ absorptive enterocytes': {'HSD17B2'},
    'MTTP/RBP2+ absorptive enterocytes': {'MTTP', 'RBP2'},
    'CKB/FABP1+ metabolic enterocytes': {'CKB', 'FABP1'},
    'NDRG1/ANKRD37+ hypoxia-response enterocytes': {'NDRG1', 'ANKRD37'},
    'MAML3/CHRM3+ epithelial cells (review)': {'MAML3', 'CHRM3'},
    'TFF1/REG4+ secretory cells': {'TFF1', 'REG4'},
    'MUC2/FCGBP+ Goblet cells': {'MUC2', 'FCGBP'},
    'SPINK4/CA4+ secretory cells': {'SPINK4', 'CA4'},
    'DDIT3/ATF3+ ER-stress enterocytes': {'DDIT3', 'ATF3', 'PPP1R15A'},
    'TA/stem-like cells': {'ASCL2', 'LGR5', 'OLFM4'},
    'S-phase enterocytes': {'FEN1', 'MCM10', 'GINS2'},
    'MIR924HG/DIAPH3+ cycling cells': {'MIR924HG', 'DIAPH3'},
    'Histone-rich cycling enterocytes': {'HIST1H1B', 'HIST1H2AJ'},
    'BIRC5/CDKN3+ cycling enterocytes': {'BIRC5', 'CDKN3'},
    'CDC20/CCNB1+ cycling enterocytes': {'CDC20', 'CCNB1'},
}

# 本次结直肠数据的人工复核显示名称。细分 marker programme 仍保存在
# ``annotation_marker_programme`` 中，避免把 marker 证据或原始 programme
# 名称改写掉；这个表只统一面向科研图、结果表和 UMAP 图例的最终显示标签。
# ``CLDN2/NME1`` 的映射是当前项目中用户确认的 cluster 3 标签，不适用于
# 其他 unresolved programme。
COLORECTAL_REFINED_DISPLAY_LABELS = {
    'SPINK4/CA4+ secretory cells': 'Goblet cells',
    'TFF1/REG4+ secretory cells': 'REG4+ secretory cells',
    'NDRG1/ANKRD37+ hypoxia-response enterocytes': 'ER-Stress Enterocytes',
    'CDC20/CCNB1+ cycling enterocytes': 'G2/M TA cells',
    'Histone-rich cycling enterocytes': 'Cycling Enterocytes',
    'MAML3/CHRM3+ epithelial cells (review)': 'MAML3/CHRM3+ Epithelial cells',
    'CYP3A5+ enterocytes': 'CYP3A5+ Enterocytes',
    'S-phase enterocytes': 'S-phase TA cells',
    'HSD17B2+ absorptive enterocytes': 'Absorptive Enterocytes',
    'Unresolved epithelial programme: CLDN2/NME1 (cluster 3; review)': (
        'Inflammatory stress cells'
    ),
    'CKB/FABP1+ metabolic enterocytes': 'High metabolic cells',
}

# ---------------------------------------------------------------------------
# 版本化精细 Marker 面板契约
# ---------------------------------------------------------------------------
# 面板名是稳定的内部 key（Colorectal_refined），版本号记录在 manifest 与
# annotation_metadata 中。旧 h5ad / 旧面板文件从不改写；每次运行只新增版本化
# 副本。Colorectal_refined_v2 是同一个面板的显式别名，方便用户在下游配置里
# 指名“细标签优先”的判定规则。
REFINED_PANEL_VERSION = 'colorectal_refined_v2'
PANEL_ALIASES = {'Colorectal_refined_v2': 'Colorectal_refined'}

# 自动注释的证据等级。``provisional`` 表示“有分数/Marker 支持但不足以确认”，
# ``anchor-supported`` 表示“单个强特异锚点 Marker 支持”，``unresolved-review``
# 表示“证据不足，仅保留独立 cluster 程序描述，等待人工复核”。
EVIDENCE_TIERS = ('confirmed', 'anchor-supported', 'provisional', 'unresolved-review')

# Auto 模式三路路由的检测阈值。细分面板判定要求每个 programme 至少有
# min_detected_genes_per_programme 个 marker 在 detection_fraction_threshold
# 比例的细胞中检出，且支持的程序数达到 min_supported_programmes。
AUTO_PANEL_DETECTION = {
    'detection_fraction_threshold': 0.05,
    'min_detected_genes_per_programme': 2,
    'min_supported_programmes': 3,
}


def normalize_marker_set_name(name):
    """把用户可见的别名映射回内部稳定面板 key。"""
    return PANEL_ALIASES.get(str(name or '').strip(), str(name or '').strip())


def curated_annotation_display_label(label, marker_set_name):
    """Return the reviewed display label without discarding marker provenance."""
    marker_set_name = normalize_marker_set_name(marker_set_name)
    value = str(label or 'Unknown')
    if marker_set_name == 'Colorectal_refined':
        return COLORECTAL_REFINED_DISPLAY_LABELS.get(value, value)
    return value


def apply_curated_annotation_display_labels(adata, marker_set_name):
    """Apply reviewed display aliases and retain the original programme label.

    This is intentionally a label-only, post-decision step.  It neither
    changes marker scores nor evidence tiers; researchers can always inspect
    ``annotation_marker_programme`` to see the programme that produced a
    displayed cell-type name.
    """
    import numpy as np

    if 'celltype' not in adata.obs.columns:
        return 0
    marker_set_name = normalize_marker_set_name(marker_set_name)
    if marker_set_name != 'Colorectal_refined':
        return 0
    programme_labels = adata.obs['celltype'].astype(str)
    display_labels = programme_labels.map(
        lambda value: curated_annotation_display_label(value, marker_set_name)
    )
    changed = display_labels.ne(programme_labels)
    adata.obs['annotation_marker_programme'] = programme_labels.astype('category')
    adata.obs['annotation_display_label_source'] = np.where(
        changed, 'curated_current_run_map', 'marker_programme',
    )
    adata.obs['celltype'] = display_labels.astype('category')
    return int(changed.sum())


def derive_panel_anchor_markers(markers):
    """按面板结构自动派生“单程序专属”锚点 marker。

    一个 marker 只出现在一个 programme 中时，它天然是比共享分数模块更强的
    判别信号。该派生规则与数据无关、可审计，用于 Organoid 等没有手工锚点表
    的面板；Colorectal refined 仍优先使用人工确认的锚点表。
    """
    counts = {}
    for cell_type, genes in (markers or {}).items():
        for gene in genes:
            counts.setdefault(str(gene).strip().upper(), set()).add(str(cell_type))
    derived = {}
    for gene, owners in counts.items():
        if len(owners) == 1:
            derived.setdefault(next(iter(owners)), set()).add(gene)
    return {str(key): set(value) for key, value in derived.items()}


def build_marker_panel_manifest(
    marker_set_name,
    organoid_type='intestinal',
    markers=None,
    anchors=None,
    negative_markers=None,
    params=None,
):
    """构造一次运行可下载、可复用的 marker panel manifest。

    清单记录面板版本、每个 programme 的 marker/锚点/负向 marker、判定规则与
    阈值。保存为 annotation_marker_panel_manifest.json，并写入 h5ad uns，使
    旧输出也能回查当时用的面板。
    """
    markers = dict(markers or {})
    marker_set_name = normalize_marker_set_name(marker_set_name)
    params = dict(params or {})
    if marker_set_name == 'Colorectal_refined':
        panel_version = REFINED_PANEL_VERSION
        display_name = 'Colorectal refined v2（cluster marker 优先的肠上皮程序面板）'
    elif marker_set_name == 'Colorectal':
        panel_version = 'colorectal_v1'
        display_name = 'Colorectal（大谱系 + 肠上皮亚型两级面板）'
    elif marker_set_name == 'Organoid':
        panel_version = f'organoid_{organoid_type}_v1'
        display_name = f'Organoid {organoid_type}（本地谱系面板）'
    else:
        panel_version = f'{marker_set_name.lower()}_v1'
        display_name = f'{marker_set_name} marker panel'
    return {
        'panel_version': panel_version,
        'marker_set': marker_set_name,
        'display_name': display_name,
        'organoid_type': organoid_type if marker_set_name == 'Organoid' else None,
        'n_programmes': len(markers),
        'decision_rules': {
            'rule_1': '两个以上 cluster 特异 marker 与 programme 重合：输出细分标签（confirmed）。',
            'rule_2': '单个强特异锚点 marker（如 CYP3A5、HSD17B2、MALRD1）：输出细分标签，标记 anchor-supported。',
            'rule_3': '细分证据不足：输出 Unresolved epithelial/organoid programme: GeneA/GeneB (cluster X; review)，不伪装成确定细胞类型。',
            'rule_4': '来源 cluster 永不自动合并；只有人工填写的 cluster_merge_map / manual_map_csv 才合并显示标签。',
            'rule_5': 'UMAP/KNN 邻接仅作为复核证据写入审阅文件，不作为自动合并依据。',
            'rule_6': '共享的高分代谢/应激模块不能覆盖有独立 cluster marker 的程序。',
        },
        'thresholds': {
            'signature_overlap_min': 2,
            'single_anchor_rule': True,
            'anchor_min_pct': float(params.get('marker_min_pct', 0.10)),
            'marker_min_delta_pct': float(params.get('marker_min_delta_pct', 0.05)),
            'cluster_agreement_threshold': float(params.get('cluster_agreement_threshold', 0.6)),
            'min_annotation_score': float(params.get('min_annotation_score', 0.0)),
            'auto_panel_detection': dict(AUTO_PANEL_DETECTION),
        },
        'programmes': {
            str(cell_type): {
                'markers': list(dict.fromkeys(genes)),
                'anchors': sorted(
                    str(gene) for gene in (anchors or {}).get(cell_type, set())
                ),
                'negative_markers': list(dict.fromkeys(
                    (negative_markers or {}).get(cell_type, [])
                )),
            }
            for cell_type, genes in markers.items()
        },
    }


def compute_cluster_neighborhood_evidence(adata, cluster_key, k=5):
    """计算 cluster 间 UMAP/KNN 邻接，仅作为复核证据。

    该结果绝不参与自动合并：它只回答“哪个来源 cluster 与当前 cluster 在
    UMAP 或 KNN 图上相邻”，供人工确认过渡态/可能合并对象。缺失嵌入或邻接
    图时返回空列表并保持可审计。
    """
    import numpy as np

    if cluster_key not in adata.obs.columns:
        return {}
    labels = adata.obs[cluster_key].astype(str)
    clusters = sorted(labels.unique().tolist(), key=lambda value: (len(value), value))

    umap_adjacency = {}
    if 'X_umap' in adata.obsm:
        try:
            coords = np.asarray(adata.obsm['X_umap'], dtype=float)
            if coords.shape[0] == adata.n_obs and coords.shape[1] >= 2:
                centroids = {}
                for cluster in clusters:
                    mask = labels.to_numpy() == cluster
                    centroids[cluster] = coords[mask].mean(axis=0)
                for cluster in clusters:
                    distances = [
                        (other, float(np.linalg.norm(centroids[cluster] - centroids[other])))
                        for other in clusters if other != cluster
                    ]
                    distances.sort(key=lambda item: (item[1], item[0]))
                    umap_adjacency[cluster] = [other for other, _ in distances[:k]]
        except Exception:
            umap_adjacency = {}

    knn_adjacency = {}
    try:
        if 'connectivities' in adata.obsp:
            connectivity = adata.obsp['connectivities'].tocsr()
            if connectivity.shape == (adata.n_obs, adata.n_obs):
                for cluster in clusters:
                    mask = labels.to_numpy() == cluster
                    rows = connectivity[mask]
                    total = float(rows.sum()) if rows.nnz else 0.0
                    if total <= 0:
                        knn_adjacency[cluster] = []
                        continue
                    column_sums = np.asarray(rows.sum(axis=0)).reshape(-1)
                    neighbours = []
                    for other in clusters:
                        if other == cluster:
                            continue
                        other_mask = labels.to_numpy() == other
                        weight = float(column_sums[other_mask].sum()) / total
                        if weight > 0:
                            neighbours.append((other, weight))
                    neighbours.sort(key=lambda item: (-item[1], item[0]))
                    knn_adjacency[cluster] = [other for other, _ in neighbours[:k]]
    except Exception:
        knn_adjacency = {}

    return {
        cluster: {
            'umap_adjacent_clusters': list(umap_adjacency.get(cluster, [])),
            'knn_adjacent_clusters': list(knn_adjacency.get(cluster, [])),
        }
        for cluster in clusters
    }


def parse_annotation_manual_map_csv(path_or_text):
    """解析人工审阅后的 annotation_manual_map_template.csv。

    要求 CSV 至少包含 ``cluster`` 与 ``manual_label`` 两列；其余列为说明性
    内容。返回 ``{cluster_id: final_label}``，未填写或等于 ``KEEP``/原标签的
    条目保持不变。路径优先，其次是文本内容，便于测试与 API 直接传字符串。
    """
    import pandas as pd

    if isinstance(path_or_text, (str, os.PathLike)) and os.path.isfile(str(path_or_text)):
        frame = pd.read_csv(str(path_or_text), dtype=str)
    else:
        import io
        frame = pd.read_csv(io.StringIO(str(path_or_text)), dtype=str)
    if frame is None or frame.empty:
        return {}
    frame.columns = [str(column).strip() for column in frame.columns]
    if 'cluster' not in frame.columns or 'manual_label' not in frame.columns:
        raise ValueError(
            'annotation_manual_map_template.csv 必须包含 cluster 与 manual_label 两列。'
        )
    mapping = {}
    for _, row in frame.iterrows():
        cluster = str(row['cluster']).strip()
        label = str(row.get('manual_label', '') or '').strip()
        if not cluster or not label or label.lower() in {'keep', 'nan'} or label == 'nan':
            continue
        mapping[cluster] = label
    return mapping

# Tissue-specific organoid panels.  These are intentionally compact, lineage-
# oriented first-pass signatures rather than a reference atlas.  Organoids
# vary across biological conditions, so the result stays
# in the review path and exposes matched-gene coverage rather than treating a
# fixed panel as a universal truth.
ORGANOID_MARKER_SETS = {
    'intestinal': {
        'Intestinal stem cells': ['LGR5', 'ASCL2', 'OLFM4', 'SMOC2'],
        'Proliferating cells': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
        'Enterocytes': ['ALPI', 'VIL1', 'FABP1', 'SI'],
        'Goblet cells': ['MUC2', 'SPINK4', 'TFF3', 'AGR2'],
        'Enteroendocrine cells': ['CHGA', 'CHGB', 'NEUROD1', 'PAX6'],
        'Paneth cells': ['LYZ', 'MMP7', 'DEFA5', 'DEFA6'],
        'Tuft cells': ['POU2F3', 'TRPM5', 'SOX9', 'AVIL'],
        'BEST4+ enterocytes': ['BEST4', 'OTOP2', 'CA7', 'GUCA2A'],
    },
    'cerebral': {
        'Neural progenitors/radial glia': ['SOX2', 'PAX6', 'VIM', 'HES1'],
        'Intermediate neural progenitors': ['EOMES', 'NEUROD1', 'TBR1', 'HES6'],
        'Proliferating neural progenitors': ['MKI67', 'TOP2A', 'CENPF', 'PCNA'],
        'Neurons': ['DCX', 'TUBB3', 'MAP2', 'RBFOX3'],
        'Excitatory neurons': ['SATB2', 'SLC17A7', 'CAMK2A', 'TBR1'],
        'Interneurons': ['DLX1', 'DLX2', 'GAD1', 'GAD2'],
        'Astroglia/astrocytes': ['GFAP', 'S100B', 'AQP4', 'ALDH1L1'],
        'OPC/oligodendrocytes': ['OLIG1', 'OLIG2', 'SOX10', 'PDGFRA'],
        'Choroid plexus': ['TTR', 'OTX2', 'KCNJ13', 'AQP1'],
    },
    'kidney': {
        'Nephron progenitors': ['SIX2', 'CITED1', 'EYA1', 'LHX1'],
        'Podocytes': ['NPHS1', 'NPHS2', 'PODXL', 'WT1'],
        'Proximal tubule': ['LRP2', 'SLC3A1', 'SLC34A1', 'HNF4A'],
        'Loop of Henle': ['SLC12A1', 'UMOD', 'CLDN10', 'POU3F3'],
        'Distal tubule': ['GATA3', 'KCNJ10', 'SLC12A3', 'TRPM6'],
        'Collecting duct/ureteric epithelium': ['AQP2', 'KRT8', 'KRT18', 'EPCAM'],
        'Stromal cells': ['COL1A1', 'COL3A1', 'DCN', 'COL1A2'],
        'Endothelial cells': ['EMCN', 'KDR', 'PECAM1', 'VWF'],
    },
    'liver': {
        'Hepatocyte-like cells': ['HNF4A', 'ALB', 'APOA1', 'TTR', 'SERPINA1'],
        'Hepatoblast/progenitor cells': ['AFP', 'KRT19', 'EPCAM', 'SOX9'],
        'Cholangiocyte-like cells': ['KRT7', 'KRT19', 'EPCAM', 'KRT18', 'SOX9'],
        'Stellate/mesenchymal cells': ['COL1A1', 'COL3A1', 'COL1A2', 'DCN'],
        'Endothelial cells': ['EMCN', 'KDR', 'PECAM1', 'VWF'],
        'Proliferating cells': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
    'lung': {
        'Alveolar type 1 cells': ['PDPN', 'AGER', 'CAV1', 'EMP2'],
        'Alveolar type 2 cells': ['SFTPC', 'SFTPA1', 'SFTPB', 'ABCA3'],
        'Club cells': ['SCGB1A1', 'KRT19', 'KRT8', 'KRT18'],
        'Ciliated cells': ['FOXJ1', 'PIFO', 'TPPP3', 'CAPS'],
        'Basal cells': ['TP63', 'KRT5', 'KRT14', 'KRT17'],
        'Goblet cells': ['MUC5AC', 'BPIFB1', 'AGR2', 'SPDEF'],
        'Proliferating cells': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
        'Mesenchymal cells': ['COL1A1', 'COL3A1', 'DCN', 'LUM'],
    },
    'pancreatic': {
        'Ductal cells': ['KRT19', 'KRT8', 'KRT18', 'KRT7'],
        'Ductal progenitors': ['SOX9', 'SPP1', 'KRT19', 'KRT8'],
        'Endocrine/islet cells': ['CHGA', 'CHGB', 'ISL1', 'NEUROD1'],
        'Beta cells': ['INS', 'IAPP', 'PCSK1', 'MAFA'],
        'Alpha cells': ['GCG', 'LOXL4', 'IRX2', 'TTR'],
        'Acinar cells': ['PRSS1', 'PRSS2', 'REG1A', 'CPA1'],
        'Stellate cells': ['COL1A1', 'COL3A1', 'COL1A2', 'COL6A1'],
        'Proliferating cells': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
    'cardiac': {
        'Cardiomyocytes': ['TNNT2', 'MYH6', 'NKX2-5', 'ACTN2', 'TNNI3'],
        'Cardiac progenitors': ['ISL1', 'TBX5', 'HAND1', 'GATA4'],
        'Endothelial cells': ['EMCN', 'PECAM1', 'VWF', 'KDR'],
        'Smooth muscle/pericytes': ['ACTA2', 'TAGLN', 'RGS5', 'PDGFRB'],
        'Fibroblasts': ['COL1A1', 'COL3A1', 'DCN', 'LUM'],
        'Epicardial cells': ['WT1', 'TBX18', 'TCF21', 'KRT19'],
        'Proliferating cells': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
}

ORGANOID_MARKER_SET_LABELS = {
    'intestinal': 'Intestinal',
    'cerebral': 'Cerebral/brain',
    'kidney': 'Kidney',
    'liver': 'Liver',
    'lung': 'Lung',
    'pancreatic': 'Pancreatic',
    'cardiac': 'Cardiac',
}

# Organoid maturity is measured from the observed expression and, when
# present, the sample's own time metadata.  These compact modules are not
# used to rename a cell type; they produce a continuous review signal.
ORGANOID_MATURITY_MARKERS = {
    'intestinal': {
        'progenitor': ['LGR5', 'ASCL2', 'OLFM4', 'SMOC2'],
        'mature': ['ALPI', 'VIL1', 'FABP1', 'SI', 'MUC2', 'CHGA', 'LYZ'],
        'cycling': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
    'cerebral': {
        'progenitor': ['SOX2', 'PAX6', 'VIM', 'HES1', 'EOMES'],
        'mature': ['DCX', 'TUBB3', 'MAP2', 'RBFOX3', 'GFAP', 'S100B'],
        'cycling': ['MKI67', 'TOP2A', 'CENPF', 'PCNA'],
    },
    'kidney': {
        'progenitor': ['SIX2', 'CITED1', 'EYA1', 'LHX1'],
        'mature': ['NPHS1', 'NPHS2', 'SLC34A1', 'UMOD', 'SLC12A3', 'AQP2'],
        'cycling': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
    'liver': {
        'progenitor': ['AFP', 'KRT19', 'EPCAM', 'SOX9'],
        'mature': ['ALB', 'HNF4A', 'APOA1', 'TTR', 'SERPINA1', 'KRT7'],
        'cycling': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
    'lung': {
        'progenitor': ['TP63', 'KRT5', 'KRT14', 'KRT17', 'SCGB1A1'],
        'mature': ['AGER', 'PDPN', 'SFTPC', 'SFTPA1', 'FOXJ1', 'PIFO'],
        'cycling': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
    'pancreatic': {
        'progenitor': ['SOX9', 'SPP1', 'KRT19', 'KRT8'],
        'mature': ['INS', 'IAPP', 'GCG', 'PRSS1', 'CPA1', 'CHGA'],
        'cycling': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
    'cardiac': {
        'progenitor': ['ISL1', 'TBX5', 'HAND1', 'GATA4'],
        'mature': ['TNNT2', 'MYH6', 'ACTN2', 'TNNI3', 'PECAM1', 'VWF'],
        'cycling': ['MKI67', 'TOP2A', 'PCNA', 'STMN1'],
    },
}

# Negative evidence is intentionally conservative.  A listed gene lowers a
# candidate score only when it is detectably enriched; it never deletes a
# label by itself.  This is especially important for transitional organoid
# states where two lineages can be biologically real rather than technical
# contamination.
DEFAULT_NEGATIVE_MARKERS = {
    'Epithelial': ['PTPRC', 'LYZ', 'COL1A1', 'PECAM1'],
    'Endothelial': ['EPCAM', 'KRT8', 'KRT18', 'COL1A1', 'PTPRC'],
    'Fibroblast': ['EPCAM', 'KRT8', 'KRT18', 'PECAM1', 'PTPRC'],
    'Pericyte/Smooth muscle': ['EPCAM', 'KRT8', 'KRT18', 'PTPRC'],
    'Myeloid': ['EPCAM', 'KRT8', 'KRT18', 'COL1A1', 'PECAM1'],
    'T cells': ['EPCAM', 'KRT8', 'KRT18', 'LYZ', 'COL1A1'],
    'NK cells': ['EPCAM', 'KRT8', 'KRT18', 'COL1A1', 'CD3D'],
    'B cells': ['EPCAM', 'KRT8', 'KRT18', 'LYZ', 'COL1A1'],
    'Plasma cells': ['EPCAM', 'KRT8', 'KRT18', 'LYZ', 'COL1A1'],
    'Mast cells': ['EPCAM', 'KRT8', 'KRT18', 'LYZ', 'COL1A1'],
    'Cycling cells': [],
}

COLORECTAL_NEGATIVE_MARKERS = {
    'Epithelial': ['PTPRC', 'LST1', 'COL1A1', 'PECAM1'],
    'T cells': ['EPCAM', 'KRT8', 'KRT18', 'LST1', 'CSF1R'],
    'Myeloid': ['EPCAM', 'KRT8', 'KRT18', 'CD3D', 'TRAC'],
    'Fibroblast': ['EPCAM', 'KRT8', 'KRT18', 'PTPRC', 'PECAM1'],
    'Endothelial': ['EPCAM', 'KRT8', 'KRT18', 'PTPRC', 'COL1A1'],
    'Neural-like (review)': ['PTPRC', 'EPCAM', 'COL1A1', 'PECAM1'],
}
for _epithelial_subtype in COLORECTAL_EPITHELIAL_SUBTYPE_MARKERS:
    COLORECTAL_NEGATIVE_MARKERS[_epithelial_subtype] = [
        'PTPRC', 'LST1', 'CSF1R', 'COL1A1', 'PECAM1',
    ]
for _refined_label in COLORECTAL_REFINED_MARKERS:
    COLORECTAL_NEGATIVE_MARKERS[_refined_label] = [
        'PTPRC', 'LST1', 'CSF1R', 'COL1A1', 'PECAM1',
    ]

ORGANOID_NEGATIVE_MARKERS = {
    'intestinal': {
        'Intestinal stem cells': ['MUC2', 'ALPI', 'CHGA', 'LYZ'],
        'Enterocytes': ['MUC2', 'CHGA', 'LYZ', 'MKI67'],
        'Goblet cells': ['ALPI', 'CHGA', 'LYZ', 'MKI67'],
        'Enteroendocrine cells': ['MUC2', 'ALPI', 'LYZ', 'MKI67'],
        'Paneth cells': ['MUC2', 'ALPI', 'CHGA', 'MKI67'],
        'Tuft cells': ['MUC2', 'ALPI', 'CHGA', 'LYZ'],
        'BEST4+ enterocytes': ['MUC2', 'CHGA', 'LYZ', 'MKI67'],
        'Proliferating cells': ['MUC2', 'ALPI', 'CHGA', 'LYZ'],
    },
    'cerebral': {
        'Neural progenitors/radial glia': ['DCX', 'MAP2', 'RBFOX3', 'GFAP'],
        'Intermediate neural progenitors': ['MAP2', 'RBFOX3', 'GFAP'],
        'Proliferating neural progenitors': ['DCX', 'MAP2', 'RBFOX3', 'GFAP'],
        'Neurons': ['SOX2', 'PAX6', 'MKI67', 'GFAP'],
        'Excitatory neurons': ['GAD1', 'GAD2', 'SOX2', 'MKI67'],
        'Interneurons': ['SLC17A7', 'CAMK2A', 'SOX2', 'MKI67'],
        'Astroglia/astrocytes': ['DCX', 'MAP2', 'RBFOX3', 'MKI67'],
        'OPC/oligodendrocytes': ['DCX', 'MAP2', 'GFAP', 'MKI67'],
    },
    'kidney': {
        'Nephron progenitors': ['NPHS1', 'NPHS2', 'UMOD', 'SLC12A3'],
        'Podocytes': ['SIX2', 'CITED1', 'SLC34A1', 'AQP2'],
        'Proximal tubule': ['NPHS1', 'NPHS2', 'UMOD', 'AQP2'],
        'Loop of Henle': ['SLC34A1', 'NPHS1', 'NPHS2', 'AQP2'],
        'Distal tubule': ['SLC34A1', 'UMOD', 'NPHS1', 'AQP2'],
        'Collecting duct/ureteric epithelium': ['SLC34A1', 'UMOD', 'NPHS1'],
    },
    'liver': {
        'Hepatocyte-like cells': ['KRT7', 'KRT19', 'COL1A1', 'MKI67'],
        'Hepatoblast/progenitor cells': ['ALB', 'APOA1', 'COL1A1', 'MKI67'],
        'Cholangiocyte-like cells': ['ALB', 'APOA1', 'COL1A1', 'MKI67'],
        'Stellate/mesenchymal cells': ['ALB', 'EPCAM', 'KRT19', 'PECAM1'],
        'Endothelial cells': ['ALB', 'KRT19', 'COL1A1', 'EPCAM'],
        'Proliferating cells': ['ALB', 'KRT19', 'COL1A1', 'PECAM1'],
    },
    'lung': {
        'Alveolar type 1 cells': ['SFTPC', 'SCGB1A1', 'FOXJ1', 'MKI67'],
        'Alveolar type 2 cells': ['AGER', 'PDPN', 'FOXJ1', 'MKI67'],
        'Club cells': ['SFTPC', 'AGER', 'FOXJ1', 'MKI67'],
        'Ciliated cells': ['SFTPC', 'AGER', 'SCGB1A1', 'MKI67'],
        'Basal cells': ['SFTPC', 'AGER', 'SCGB1A1', 'FOXJ1'],
        'Goblet cells': ['SFTPC', 'AGER', 'FOXJ1', 'MKI67'],
        'Proliferating cells': ['SFTPC', 'AGER', 'SCGB1A1', 'COL1A1'],
        'Mesenchymal cells': ['SFTPC', 'AGER', 'SCGB1A1', 'EPCAM'],
    },
    'pancreatic': {
        'Ductal cells': ['INS', 'GCG', 'PRSS1', 'COL1A1'],
        'Ductal progenitors': ['INS', 'GCG', 'PRSS1', 'COL1A1'],
        'Endocrine/islet cells': ['KRT19', 'KRT8', 'PRSS1', 'COL1A1'],
        'Beta cells': ['GCG', 'KRT19', 'PRSS1', 'COL1A1'],
        'Alpha cells': ['INS', 'KRT19', 'PRSS1', 'COL1A1'],
        'Acinar cells': ['INS', 'GCG', 'KRT19', 'COL1A1'],
        'Stellate cells': ['INS', 'GCG', 'KRT19', 'EPCAM'],
    },
    'cardiac': {
        'Cardiomyocytes': ['PECAM1', 'COL1A1', 'RGS5', 'MKI67'],
        'Cardiac progenitors': ['TNNT2', 'MYH6', 'COL1A1', 'MKI67'],
        'Endothelial cells': ['TNNT2', 'MYH6', 'COL1A1', 'EPCAM'],
        'Smooth muscle/pericytes': ['TNNT2', 'MYH6', 'PECAM1', 'EPCAM'],
        'Fibroblasts': ['TNNT2', 'MYH6', 'PECAM1', 'EPCAM'],
        'Epicardial cells': ['TNNT2', 'MYH6', 'PECAM1', 'MKI67'],
        'Proliferating cells': ['TNNT2', 'MYH6', 'COL1A1', 'PECAM1'],
    },
}

UNIVERSAL_LABELS = {
    'Epithelial': ('Non-immune cell', 'Epithelial cell', 'Epithelial cell', 'CL:0000066'),
    'Endothelial': ('Non-immune cell', 'Endothelial cell', 'Endothelial cell', 'CL:0000115'),
    'Fibroblast': ('Stromal cell', 'Fibroblast', 'Fibroblast', 'CL:0000057'),
    'Pericyte/Smooth muscle': ('Stromal cell', 'Perivascular cell', 'Pericyte/smooth muscle cell', 'CL:0000669'),
    'Myeloid': ('Immune cell', 'Myeloid cell', 'Myeloid cell', 'CL:0000763'),
    'T cells': ('Immune cell', 'T cell', 'T cell', 'CL:0000084'),
    'NK cells': ('Immune cell', 'NK cell', 'Natural killer cell', 'CL:0000623'),
    'B cells': ('Immune cell', 'B cell', 'B cell', 'CL:0000236'),
    'Plasma cells': ('Immune cell', 'B cell', 'Plasma cell', 'CL:0000786'),
    'Mast cells': ('Immune cell', 'Mast cell', 'Mast cell', 'CL:0000097'),
    'Cycling cells': ('Cell state', 'Cycling cell', 'Cycling cell', ''),
}

# 初步注释的细化（fine-grained）第二级面板。fine_annotation 模式在
# DEFAULT_UNIVERSAL_MARKERS 大谱系判定之后，只对带细分组的谱系再做一次
# 亚型判定；大谱系标签本身始终保留在候选列表里，因此亚型证据不足的簇会
# 回退到大谱系标签，而不是被强行塞进最近的亚型或被标成 Unknown。
UNIVERSAL_FINE_SUBTYPE_MARKERS = {
    # T / NK 轴
    'CD4+ T cells': ['CD4', 'IL7R', 'LEF1', 'CCR7', 'TCF7'],
    'CD8+ T cells': ['CD8A', 'CD8B', 'GZMK', 'GZMA', 'CCL5'],
    'Treg cells': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'BATF'],
    'Gamma-delta T cells': ['TRDC', 'TRGC1', 'TRGC2', 'TRDV2'],
    # Myeloid 轴
    'Monocyte/Macrophage': ['CD14', 'LYZ', 'CSF1R', 'CD68', 'FCN1'],
    'cDC1': ['CLEC9A', 'XCR1', 'CADM1', 'BATF3'],
    'cDC2': ['CD1C', 'FCER1A', 'CLEC10A', 'HLA-DPA1'],
    'pDC': ['IL3RA', 'TCF4', 'IRF7', 'LILRA4', 'GZMB'],
    'Neutrophils': ['FCGR3B', 'CSF3R', 'CXCR2', 'S100A12'],
    # 间质/血管轴（拆分原来的合并标签）
    'Pericytes': ['RGS5', 'PDGFRB', 'NOTCH3', 'MCAM', 'KCNJ8'],
    'Smooth muscle cells': ['ACTA2', 'TAGLN', 'MYH11', 'CNN1'],
    'Lymphatic endothelial': ['LYVE1', 'PROX1', 'FLT4', 'PDPN'],
}

# 亚型 → 所属大谱系；用于把第二级判定限制在同一谱系内，并防止同谱系
# 亚型被当成“谱系混合/环境 RNA”证据。
UNIVERSAL_FINE_LINEAGE_MAP = {
    'CD4+ T cells': 'T cells',
    'CD8+ T cells': 'T cells',
    'Treg cells': 'T cells',
    'Gamma-delta T cells': 'T cells',
    'Monocyte/Macrophage': 'Myeloid',
    'cDC1': 'Myeloid',
    'cDC2': 'Myeloid',
    'pDC': 'Myeloid',
    'Neutrophils': 'Myeloid',
    'Pericytes': 'Pericyte/Smooth muscle',
    'Smooth muscle cells': 'Pericyte/Smooth muscle',
    'Lymphatic endothelial': 'Endothelial',
}

# 亚型标签在 cell_lineage / cell_type_l1..l3 / cell_ontology_id
# 层级列中的稳定映射，避免细化标签落入 "Unresolved lineage" 兜底。
UNIVERSAL_FINE_LABELS = {
    'CD4+ T cells': ('Immune cell', 'T cell', 'CD4+ T cell', 'CL:0000624'),
    'CD8+ T cells': ('Immune cell', 'T cell', 'CD8+ T cell', 'CL:0000625'),
    'Treg cells': ('Immune cell', 'T cell', 'Regulatory T cell', 'CL:0000815'),
    'Gamma-delta T cells': ('Immune cell', 'T cell', 'Gamma-delta T cell', 'CL:0000798'),
    'Monocyte/Macrophage': ('Immune cell', 'Myeloid cell', 'Monocyte/macrophage', 'CL:0000576'),
    'cDC1': ('Immune cell', 'Dendritic cell', 'Conventional type 1 dendritic cell', 'CL:0000990'),
    'cDC2': ('Immune cell', 'Dendritic cell', 'Conventional type 2 dendritic cell', 'CL:0000991'),
    'pDC': ('Immune cell', 'Dendritic cell', 'Plasmacytoid dendritic cell', 'CL:0000784'),
    'Neutrophils': ('Immune cell', 'Myeloid cell', 'Neutrophil', 'CL:0000775'),
    'Pericytes': ('Stromal cell', 'Perivascular cell', 'Pericyte', 'CL:0000669'),
    'Smooth muscle cells': ('Stromal cell', 'Perivascular cell', 'Smooth muscle cell', 'CL:0000192'),
    'Lymphatic endothelial': ('Non-immune cell', 'Endothelial cell', 'Lymphatic endothelial cell', 'CL:0002138'),
}

COLORECTAL_HIERARCHY = {
    'Epithelial': ('Epithelial lineage', 'Colorectal epithelium',
                   'Epithelial subtype unresolved'),
    'T cells': ('Immune lineage', 'T cells', 'T cells'),
    'Myeloid': ('Immune lineage', 'Myeloid cells', 'Myeloid cells'),
    'Fibroblast': ('Mesenchymal lineage', 'Fibroblasts', 'Fibroblasts'),
    'Endothelial': ('Endothelial lineage', 'Endothelial cells', 'Endothelial cells'),
    'Neural-like (review)': ('Neural lineage', 'Neural-like cells',
                             'Neural-like cells (review)'),
    'Stem/crypt-like epithelial': ('Epithelial lineage', 'Colorectal epithelium',
                                    'Stem/crypt-like epithelial'),
    'TA/S-phase epithelial': ('Epithelial lineage', 'Cycling/TA epithelium',
                              'TA/S-phase epithelial'),
    'TA/G2M epithelial': ('Epithelial lineage', 'Cycling/TA epithelium',
                          'TA/G2M epithelial'),
    'Goblet-like epithelial': ('Epithelial lineage', 'Secretory epithelium',
                               'Goblet-like epithelial'),
    'Absorptive/enterocyte-like epithelial': ('Epithelial lineage', 'Absorptive epithelium',
                                               'Absorptive/enterocyte-like epithelial'),
    'BEST4+ absorptive epithelial': ('Epithelial lineage', 'Absorptive epithelium',
                                     'BEST4+ absorptive epithelial'),
    'Enteroendocrine-like epithelial': ('Epithelial lineage', 'Secretory epithelium',
                                        'Enteroendocrine-like epithelial'),
    'Paneth/LYZ+ secretory epithelial': ('Epithelial lineage', 'Secretory epithelium',
                                         'Paneth/LYZ+ secretory epithelial'),
    'Inflammatory epithelial': ('Epithelial lineage', 'Inflammatory epithelium',
                                'Inflammatory epithelial'),
    'Regenerative/stress epithelial': ('Epithelial lineage', 'Regenerative epithelium',
                                        'Regenerative/stress epithelial'),
}

COLORECTAL_REFINED_HIERARCHY = {
    'High metabolic enterocytes': ('Epithelial lineage', 'Absorptive epithelium',
                                   'High metabolic enterocytes'),
    'CYP3A5+ enterocytes': ('Epithelial lineage', 'Absorptive epithelium',
                            'CYP3A5+ enterocytes'),
    'KRT20/MALRD1+ absorptive enterocytes': (
        'Epithelial lineage', 'Absorptive epithelium',
        'KRT20/MALRD1+ absorptive enterocytes',
    ),
    'HSD17B2+ absorptive enterocytes': (
        'Epithelial lineage', 'Absorptive epithelium',
        'HSD17B2+ absorptive enterocytes',
    ),
    'MTTP/RBP2+ absorptive enterocytes': (
        'Epithelial lineage', 'Absorptive epithelium',
        'MTTP/RBP2+ absorptive enterocytes',
    ),
    'CKB/FABP1+ metabolic enterocytes': (
        'Epithelial lineage', 'Absorptive epithelium',
        'CKB/FABP1+ metabolic enterocytes',
    ),
    'NDRG1/ANKRD37+ hypoxia-response enterocytes': (
        'Epithelial lineage', 'Stress-response epithelium',
        'NDRG1/ANKRD37+ hypoxia-response enterocytes',
    ),
    'MAML3/CHRM3+ epithelial cells (review)': (
        'Epithelial lineage', 'Unresolved epithelial programme',
        'MAML3/CHRM3+ epithelial cells (review)',
    ),
    'TFF1/REG4+ secretory cells': (
        'Epithelial lineage', 'Secretory epithelium',
        'TFF1/REG4+ secretory cells',
    ),
    'MUC2/FCGBP+ Goblet cells': (
        'Epithelial lineage', 'Secretory epithelium',
        'MUC2/FCGBP+ Goblet cells',
    ),
    'SPINK4/CA4+ secretory cells': (
        'Epithelial lineage', 'Secretory epithelium',
        'SPINK4/CA4+ secretory cells',
    ),
    'DDIT3/ATF3+ ER-stress enterocytes': (
        'Epithelial lineage', 'Stress-response epithelium',
        'DDIT3/ATF3+ ER-stress enterocytes',
    ),
    'TA/stem-like cells': (
        'Epithelial lineage', 'Transit-amplifying epithelium',
        'TA/stem-like cells',
    ),
    'S-phase enterocytes': (
        'Epithelial lineage', 'Cycling epithelium', 'S-phase enterocytes',
    ),
    'MIR924HG/DIAPH3+ cycling cells': (
        'Epithelial lineage', 'Cycling epithelium',
        'MIR924HG/DIAPH3+ cycling cells',
    ),
    'Histone-rich cycling enterocytes': (
        'Epithelial lineage', 'Cycling epithelium',
        'Histone-rich cycling enterocytes',
    ),
    'BIRC5/CDKN3+ cycling enterocytes': (
        'Epithelial lineage', 'Cycling epithelium',
        'BIRC5/CDKN3+ cycling enterocytes',
    ),
    'CDC20/CCNB1+ cycling enterocytes': (
        'Epithelial lineage', 'Cycling epithelium',
        'CDC20/CCNB1+ cycling enterocytes',
    ),
    'MTTP/MALRD1+ absorptive enterocytes': (
        'Epithelial lineage', 'Absorptive epithelium',
        'MTTP/MALRD1+ absorptive enterocytes',
    ),
    'KRT20/CA9+ absorptive enterocytes': (
        'Epithelial lineage', 'Absorptive epithelium',
        'KRT20/CA9+ absorptive enterocytes',
    ),
    'Absorptive enterocytes': ('Epithelial lineage', 'Absorptive epithelium',
                               'Absorptive enterocytes'),
    'Stress/secretory enterocytes': ('Epithelial lineage', 'Secretory/stress epithelium',
                                     'Stress/secretory enterocytes'),
    'ER stress enterocytes': ('Epithelial lineage', 'Stress epithelium',
                              'ER stress enterocytes'),
    'Metabolic enterocytes': ('Epithelial lineage', 'Absorptive epithelium',
                              'Metabolic enterocytes'),
    'RGS10/MT1G+ metabolic enterocytes': (
        'Epithelial lineage', 'Absorptive epithelium',
        'RGS10/MT1G+ metabolic enterocytes',
    ),
    'CKB/FABP1+ high-metabolic enterocytes': (
        'Epithelial lineage', 'Absorptive epithelium',
        'CKB/FABP1+ high-metabolic enterocytes',
    ),
    'TA cells': ('Epithelial lineage', 'Transit-amplifying epithelium', 'TA cells'),
    'Inflammatory stress cells': ('Epithelial lineage', 'Inflammatory epithelium',
                                  'Inflammatory stress cells'),
    'TA/S-phase cells': ('Epithelial lineage', 'Transit-amplifying epithelium',
                         'TA/S-phase cells'),
    'Goblet cells': ('Epithelial lineage', 'Secretory epithelium', 'Goblet cells'),
    'Goblet-like cells (low-n review)': (
        'Epithelial lineage', 'Secretory epithelium',
        'Goblet-like cells (low-n review)',
    ),
    'Cycling enterocytes': ('Epithelial lineage', 'Cycling enterocytes',
                            'Cycling enterocytes'),
    # Curated display labels for the current colorectal refined run.  Detailed
    # programme labels above remain supported so historical outputs retain
    # their original hierarchy.
    'REG4+ secretory cells': ('Epithelial lineage', 'Secretory epithelium',
                              'REG4+ secretory cells'),
    'ER-Stress Enterocytes': ('Epithelial lineage', 'Stress-response epithelium',
                              'ER-Stress Enterocytes'),
    'G2/M TA cells': ('Epithelial lineage', 'Cycling/TA epithelium',
                      'G2/M TA cells'),
    'Cycling Enterocytes': ('Epithelial lineage', 'Cycling epithelium',
                            'Cycling Enterocytes'),
    'MAML3/CHRM3+ Epithelial cells': (
        'Epithelial lineage', 'Unresolved epithelial programme',
        'MAML3/CHRM3+ Epithelial cells',
    ),
    'CYP3A5+ Enterocytes': ('Epithelial lineage', 'Absorptive epithelium',
                            'CYP3A5+ Enterocytes'),
    'S-phase TA cells': ('Epithelial lineage', 'Cycling/TA epithelium',
                         'S-phase TA cells'),
    'Absorptive Enterocytes': ('Epithelial lineage', 'Absorptive epithelium',
                               'Absorptive Enterocytes'),
    'High metabolic cells': ('Epithelial lineage', 'Absorptive epithelium',
                             'High metabolic cells'),
}

# Organoid labels are kept separate from transient states.  The hierarchy is
# intentionally broad: it provides a stable parent label without pretending
# that a compact organoid panel can resolve every fine subtype.
ORGANOID_HIERARCHY = {
    'intestinal': {
        'Intestinal stem cells': ('Epithelial lineage', 'Intestinal epithelium', 'Intestinal stem cells'),
        'Proliferating cells': ('Cell state', 'Cycling cells', 'Proliferating cells'),
        'Enterocytes': ('Epithelial lineage', 'Intestinal epithelium', 'Enterocytes'),
        'Goblet cells': ('Epithelial lineage', 'Intestinal epithelium', 'Goblet cells'),
        'Enteroendocrine cells': ('Epithelial lineage', 'Intestinal epithelium', 'Enteroendocrine cells'),
        'Paneth cells': ('Epithelial lineage', 'Intestinal epithelium', 'Paneth cells'),
        'Tuft cells': ('Epithelial lineage', 'Intestinal epithelium', 'Tuft cells'),
        'BEST4+ enterocytes': ('Epithelial lineage', 'Intestinal epithelium', 'BEST4+ enterocytes'),
    },
    'cerebral': {
        'Neural progenitors/radial glia': ('Neural lineage', 'Neural progenitors', 'Radial glia'),
        'Intermediate neural progenitors': ('Neural lineage', 'Neural progenitors', 'Intermediate neural progenitors'),
        'Proliferating neural progenitors': ('Cell state', 'Cycling neural cells', 'Proliferating neural progenitors'),
        'Neurons': ('Neural lineage', 'Neurons', 'Neurons'),
        'Excitatory neurons': ('Neural lineage', 'Neurons', 'Excitatory neurons'),
        'Interneurons': ('Neural lineage', 'Neurons', 'Interneurons'),
        'Astroglia/astrocytes': ('Neural lineage', 'Glia', 'Astroglia/astrocytes'),
        'OPC/oligodendrocytes': ('Neural lineage', 'Glia', 'OPC/oligodendrocytes'),
        'Choroid plexus': ('Epithelial lineage', 'Choroid plexus', 'Choroid plexus'),
    },
    'kidney': {
        'Nephron progenitors': ('Epithelial lineage', 'Nephron', 'Nephron progenitors'),
        'Podocytes': ('Epithelial lineage', 'Nephron', 'Podocytes'),
        'Proximal tubule': ('Epithelial lineage', 'Nephron', 'Proximal tubule'),
        'Loop of Henle': ('Epithelial lineage', 'Nephron', 'Loop of Henle'),
        'Distal tubule': ('Epithelial lineage', 'Nephron', 'Distal tubule'),
        'Collecting duct/ureteric epithelium': ('Epithelial lineage', 'Collecting duct', 'Collecting duct/ureteric epithelium'),
        'Stromal cells': ('Mesenchymal lineage', 'Stromal cells', 'Stromal cells'),
        'Endothelial cells': ('Endothelial lineage', 'Endothelial cells', 'Endothelial cells'),
    },
    'liver': {
        'Hepatocyte-like cells': ('Epithelial lineage', 'Hepatic epithelium', 'Hepatocyte-like cells'),
        'Hepatoblast/progenitor cells': ('Epithelial lineage', 'Hepatic epithelium', 'Hepatoblast/progenitor cells'),
        'Cholangiocyte-like cells': ('Epithelial lineage', 'Biliary epithelium', 'Cholangiocyte-like cells'),
        'Stellate/mesenchymal cells': ('Mesenchymal lineage', 'Stellate cells', 'Stellate/mesenchymal cells'),
        'Endothelial cells': ('Endothelial lineage', 'Endothelial cells', 'Endothelial cells'),
        'Proliferating cells': ('Cell state', 'Cycling cells', 'Proliferating cells'),
    },
    'lung': {
        'Alveolar type 1 cells': ('Epithelial lineage', 'Alveolar epithelium', 'Alveolar type 1 cells'),
        'Alveolar type 2 cells': ('Epithelial lineage', 'Alveolar epithelium', 'Alveolar type 2 cells'),
        'Club cells': ('Epithelial lineage', 'Airway epithelium', 'Club cells'),
        'Ciliated cells': ('Epithelial lineage', 'Airway epithelium', 'Ciliated cells'),
        'Basal cells': ('Epithelial lineage', 'Airway epithelium', 'Basal cells'),
        'Goblet cells': ('Epithelial lineage', 'Airway epithelium', 'Goblet cells'),
        'Proliferating cells': ('Cell state', 'Cycling cells', 'Proliferating cells'),
        'Mesenchymal cells': ('Mesenchymal lineage', 'Mesenchymal cells', 'Mesenchymal cells'),
    },
    'pancreatic': {
        'Ductal cells': ('Epithelial lineage', 'Pancreatic epithelium', 'Ductal cells'),
        'Ductal progenitors': ('Epithelial lineage', 'Pancreatic epithelium', 'Ductal progenitors'),
        'Endocrine/islet cells': ('Epithelial lineage', 'Endocrine cells', 'Endocrine/islet cells'),
        'Beta cells': ('Epithelial lineage', 'Endocrine cells', 'Beta cells'),
        'Alpha cells': ('Epithelial lineage', 'Endocrine cells', 'Alpha cells'),
        'Acinar cells': ('Epithelial lineage', 'Acinar cells', 'Acinar cells'),
        'Stellate cells': ('Mesenchymal lineage', 'Stellate cells', 'Stellate cells'),
        'Proliferating cells': ('Cell state', 'Cycling cells', 'Proliferating cells'),
    },
    'cardiac': {
        'Cardiomyocytes': ('Cardiac lineage', 'Cardiomyocytes', 'Cardiomyocytes'),
        'Cardiac progenitors': ('Cardiac lineage', 'Cardiac progenitors', 'Cardiac progenitors'),
        'Endothelial cells': ('Endothelial lineage', 'Endothelial cells', 'Endothelial cells'),
        'Smooth muscle/pericytes': ('Mesenchymal lineage', 'Perivascular cells', 'Smooth muscle/pericytes'),
        'Fibroblasts': ('Mesenchymal lineage', 'Fibroblasts', 'Fibroblasts'),
        'Epicardial cells': ('Cardiac lineage', 'Epicardial cells', 'Epicardial cells'),
        'Proliferating cells': ('Cell state', 'Cycling cells', 'Proliferating cells'),
    },
}

CELL_STATE_MARKERS = {
    'Cycling': ['MKI67', 'TOP2A', 'STMN1', 'PCNA', 'TYMS'],
    'Interferon response': ['ISG15', 'IFIT1', 'IFIT3', 'MX1', 'IFI6'],
    'Stress response': ['FOS', 'JUN', 'DUSP1', 'HSPA1A', 'DNAJB1'],
    'Hypoxia response': ['CA9', 'VEGFA', 'LDHA', 'HIF1A', 'NDRG1'],
    'Apoptosis-like': ['BAX', 'BBC3', 'PMAIP1', 'DDIT3', 'JUN'],
}

MARKER_SETS = {
    'Universal': DEFAULT_UNIVERSAL_MARKERS,
    'Colorectal': DEFAULT_COLORECTAL_MARKERS,
    'Colorectal_refined': COLORECTAL_REFINED_MARKERS,
    'TME': DEFAULT_TME_MARKERS,
    'Immune': DEFAULT_IMMUNE_MARKERS,
    'Blood': DEFAULT_BLOOD_MARKERS,
    'PBMC': DEFAULT_PBMC_MARKERS,
}


def get_marker_set(marker_set_name='Universal', organoid_type='intestinal'):
    """Return a flat marker panel for the requested tissue context.

    ``Organoid`` is resolved separately from the legacy flat ``MARKER_SETS``
    mapping so existing callers remain backwards compatible.  Invalid
    organoid types deliberately fall back to the intestinal panel; callers can
    detect the invalid request by checking ``organoid_type`` against
    ``ORGANOID_MARKER_SETS`` and emit a user-facing warning.
    """
    marker_set_name = normalize_marker_set_name(marker_set_name)
    if marker_set_name == 'Organoid':
        return ORGANOID_MARKER_SETS.get(str(organoid_type).strip().lower(),
                                        ORGANOID_MARKER_SETS['intestinal'])
    return MARKER_SETS.get(marker_set_name, DEFAULT_UNIVERSAL_MARKERS)


def _auto_panel_detection(adata, programmes, selection):
    """统计一个面板里有多少 programme 在当前数据中达到检出阈值。"""
    import numpy as np

    n_cells = max(int(getattr(adata, 'n_obs', 0)), 1)
    matrix = adata.layers['counts'] if 'counts' in adata.layers else adata.X
    var_lookup = {str(gene).upper(): index for index, gene in enumerate(adata.var_names)}
    supported = []
    evidence = {}
    for programme, genes in programmes.items():
        positions = [var_lookup[str(gene).upper()] for gene in genes
                     if str(gene).upper() in var_lookup]
        if positions:
            values = matrix[:, positions]
            if hasattr(values, 'getnnz'):
                detected = np.asarray(values.getnnz(axis=0)).reshape(-1) / n_cells
            else:
                detected = (np.asarray(values) > 0).mean(axis=0)
            expressed = int(
                (detected >= selection['detection_fraction_threshold']).sum()
            )
        else:
            expressed = 0
        evidence[programme] = {
            'matched_genes': len(positions),
            'detected_genes': expressed,
        }
        if expressed >= selection['min_detected_genes_per_programme']:
            supported.append(programme)
    return supported, evidence


def resolve_auto_marker_set(adata, requested_marker_set='Auto', organoid_type='intestinal'):
    """用户选择 ``Auto`` 时，以证据强度决定三路组织注释面板。

    - 缺少明确组织程序 → Universal
    - 有肠上皮程序但细分证据有限 → Colorectal
    - 多个细分 programme 均有 marker 支持 → Colorectal refined v2

    常规 10x 的全基因清单本身不能作为组织证据；只有表达层面检测到足够多
    的独立程序才会升级面板。升级只改变面板选择，从不合并来源 cluster。
    """
    requested = normalize_marker_set_name(requested_marker_set or 'Auto') or 'Auto'
    selection = {
        'requested': requested,
        'applied': requested,
        'mode': 'explicit',
        'detection_fraction_threshold': AUTO_PANEL_DETECTION['detection_fraction_threshold'],
        'min_detected_genes_per_programme': AUTO_PANEL_DETECTION['min_detected_genes_per_programme'],
        'min_supported_programmes': AUTO_PANEL_DETECTION['min_supported_programmes'],
        'programme_evidence': {},
        'refined_programme_evidence': {},
    }
    if requested != 'Auto':
        return requested, selection

    selection['mode'] = 'auto_expression_evidence'
    broad_supported, broad_evidence = _auto_panel_detection(
        adata, COLORECTAL_EPITHELIAL_SUBTYPE_MARKERS, selection,
    )
    refined_supported, refined_evidence = _auto_panel_detection(
        adata, COLORECTAL_REFINED_MARKERS, selection,
    )
    selection['programme_evidence'] = broad_evidence
    selection['supported_programmes'] = broad_supported
    selection['n_supported_programmes'] = len(broad_supported)
    selection['refined_programme_evidence'] = refined_evidence
    selection['supported_refined_programmes'] = refined_supported
    selection['n_supported_refined_programmes'] = len(refined_supported)
    selection['panel_version'] = None

    min_supported = selection['min_supported_programmes']
    if (
        len(refined_supported) >= min_supported
        and len(broad_supported) >= min_supported
    ):
        selection['applied'] = 'Colorectal_refined'
        selection['panel_version'] = REFINED_PANEL_VERSION
        selection['reason'] = (
            '检测到多个可独立判别的肠上皮细分 programme，'
            '自动使用 colorectal_refined_v2 面板进行 cluster marker 优先的精细注释。'
        )
    elif len(broad_supported) >= min_supported:
        selection['applied'] = 'Colorectal'
        selection['reason'] = (
            '检测到多个肠上皮亚型 Marker 程序，但细分 programme 证据有限，'
            '使用 Colorectal 面板进行保守二级注释。'
        )
    else:
        selection['applied'] = 'Universal'
        selection['reason'] = (
            '未检测到足够的肠上皮亚型程序，使用 Universal 大谱系面板；'
            '不同类器官/组织样品请显式选择 Organoid 或对应面板。'
        )
    return selection['applied'], selection


def get_negative_marker_set(marker_set_name='Universal', organoid_type='intestinal'):
    """Return conservative anti-marker modules for the selected tissue panel."""
    marker_set_name = normalize_marker_set_name(marker_set_name)
    if marker_set_name in {'Colorectal', 'Colorectal_refined'}:
        return {
            cell_type: list(genes)
            for cell_type, genes in COLORECTAL_NEGATIVE_MARKERS.items()
        }
    if str(marker_set_name) == 'Organoid':
        organoid_type = str(organoid_type).strip().lower()
        panel = ORGANOID_NEGATIVE_MARKERS.get(organoid_type, {})
        return {cell_type: list(genes) for cell_type, genes in panel.items()}
    return {
        cell_type: list(genes)
        for cell_type, genes in DEFAULT_NEGATIVE_MARKERS.items()
        if cell_type in get_marker_set(marker_set_name)
    }


def detect_maturity_time_key(adata, requested=None):
    """Find an observed culture/time column without inventing a stage label."""
    requested = str(requested or '').strip()
    if requested and requested.lower() not in {'auto', 'none'}:
        return requested if requested in adata.obs.columns else None
    preferred = ('timepoint', 'time_point', 'culture_day', 'day', 'days', 'hour', 'week', 'stage')
    technical = {'batch', 'technical_batch', 'sequencing_batch', 'library_batch'}
    lower = {str(column).lower(): str(column) for column in adata.obs.columns}
    for name in preferred:
        column = lower.get(name)
        if not column or name in technical:
            continue
        values = adata.obs[column].dropna().astype(str).unique()
        if 2 <= len(values) <= 50:
            return column
    return None


def _coerce_maturity_time(values):
    """Convert common D7/7d/48h labels to a numeric descriptive time axis."""
    import re
    import numpy as np
    import pandas as pd

    numeric = pd.to_numeric(values, errors='coerce')
    if numeric.notna().all():
        return numeric.to_numpy(dtype=float)
    converted = []
    for value in values.astype(str):
        match = re.search(r'(-?\d+(?:\.\d+)?)\s*(day|d|hour|h|week|w)?', value.lower())
        if not match:
            converted.append(np.nan)
            continue
        number = float(match.group(1))
        unit = match.group(2) or 'day'
        if unit in {'hour', 'h'}:
            number /= 24.0
        elif unit in {'week', 'w'}:
            number *= 7.0
        converted.append(number)
    return np.asarray(converted, dtype=float)


def compute_organoid_maturity_evidence(adata, organoid_type='intestinal', time_key=None):
    """Estimate progenitor-to-mature state from expression and observed time metadata.

    The result is descriptive and independent of the cell-type assignment.  It
    never uses a user-selected ``early``/``mature`` prior to change scores.
    """
    import numpy as np
    import pandas as pd

    organoid_type = str(organoid_type or 'intestinal').strip().lower()
    modules = ORGANOID_MATURITY_MARKERS.get(
        organoid_type, ORGANOID_MATURITY_MARKERS['intestinal']
    )
    var_lookup = {str(gene).upper(): str(gene) for gene in adata.var_names}
    matched = {}
    coverage = {}
    for module, genes in modules.items():
        usable = list(dict.fromkeys(
            var_lookup[str(gene).upper()] for gene in genes
            if str(gene).upper() in var_lookup
        ))
        matched[module] = usable
        coverage[module] = {'matched': len(usable), 'total': len(genes)}

    expression_source = ''
    raw_scores = {}
    for module, genes in matched.items():
        raw_scores[module], module_source = _mean_log_normalized_expression_for_genes(
            adata, genes,
        )
        if genes and not expression_source:
            expression_source = module_source

    def _robust_scale(values):
        values = np.asarray(values, dtype=float)
        positive = values[np.isfinite(values) & (values > 0)]
        if positive.size == 0:
            return np.zeros(adata.n_obs, dtype=float)
        scale = float(np.nanpercentile(positive, 95))
        return np.clip(values / (scale + 1e-12), 0.0, 1.0)

    progenitor = _robust_scale(raw_scores.get('progenitor', []))
    mature = _robust_scale(raw_scores.get('mature', []))
    cycling = _robust_scale(raw_scores.get('cycling', []))
    maturity_index = np.clip(mature - progenitor - 0.20 * cycling, -1.0, 1.0)
    has_signal = (progenitor + mature + cycling) > 0
    state = np.full(adata.n_obs, 'maturity_unresolved', dtype=object)
    state[has_signal & (maturity_index <= -0.25)] = 'progenitor_like'
    state[has_signal & (maturity_index >= 0.25)] = 'mature_like'
    state[has_signal & (maturity_index > -0.25) & (maturity_index < 0.25)] = 'transitional'

    detected_time_key = detect_maturity_time_key(adata, time_key)
    time_numeric = np.full(adata.n_obs, np.nan, dtype=float)
    time_labels = np.repeat('', adata.n_obs).astype(object)
    trend_rho = np.nan
    if detected_time_key:
        time_labels = adata.obs[detected_time_key].astype(str).to_numpy()
        time_numeric = _coerce_maturity_time(adata.obs[detected_time_key])
        valid = np.isfinite(time_numeric) & np.isfinite(maturity_index)
        if valid.sum() >= 3 and np.unique(time_numeric[valid]).size >= 2:
            trend_rho = float(pd.Series(time_numeric[valid]).corr(
                pd.Series(maturity_index[valid]), method='spearman'
            ))

    warnings = []
    if not matched.get('progenitor') or not matched.get('mature'):
        warnings.append('成熟度模块覆盖不足；仅输出可用模块，不对成熟状态作强解释。')
    if not detected_time_key:
        warnings.append('未检测到 timepoint/day/hour/stage 等元数据列，成熟度仅基于表达模块。')
    elif np.isfinite(trend_rho):
        warnings.append('成熟度与时间点的相关性为描述性结果，不能替代样本级纵向统计。')
    return {
        'progenitor_score': progenitor,
        'mature_score': mature,
        'cycling_score': cycling,
        'maturity_index': maturity_index,
        'maturity_state': state,
        'time_key': detected_time_key,
        'expression_source': expression_source or 'no matched maturity genes',
        'time_labels': time_labels,
        'time_numeric': time_numeric,
        'time_trend_spearman': trend_rho,
        'marker_coverage': coverage,
        'warnings': warnings,
    }


def hierarchy_for_label(label, marker_set_name='Universal', organoid_type='intestinal'):
    """Return stable lineage/type/subtype fields without merging cell states."""
    label = str(label or 'Unknown')
    marker_set_name = normalize_marker_set_name(marker_set_name)
    if label == 'Unknown':
        return ('Unknown', 'Unknown', 'Unknown', '')
    if marker_set_name == 'Universal':
        values = UNIVERSAL_LABELS.get(label)
        if values:
            return values
        # 细化注释模式下的亚型标签（CD4+/CD8+ T、cDC/pDC、周细胞/平滑肌、
        # 淋巴管内皮等）同样映射到稳定层级列。
        values = UNIVERSAL_FINE_LABELS.get(label)
        if values:
            return values
    if marker_set_name == 'Colorectal':
        values = COLORECTAL_HIERARCHY.get(label)
        if values:
            return (*values, '')
    if marker_set_name == 'Colorectal_refined':
        values = COLORECTAL_REFINED_HIERARCHY.get(label)
        if values:
            return (*values, '')
        if label.startswith('Unresolved epithelial programme:'):
            return (
                'Epithelial lineage', 'Unresolved epithelial programme', label, '',
            )
    if marker_set_name == 'Organoid':
        values = ORGANOID_HIERARCHY.get(str(organoid_type).lower(), {}).get(label)
        if values:
            return (*values, '')
    # 任何 “Unresolved ... programme:” 标签都是证据不足的独立 cluster 程序，
    # 必须保留为 review 而不是伪装成确定细胞类型。
    if label.startswith('Unresolved ') and ' programme:' in label:
        return ('Unresolved lineage', 'Unresolved programme', label, '')
    # Scenario panels do not yet have a curated ontology hierarchy; preserve
    # the label rather than fabricating a lineage or ontology ID.
    return ('Unresolved lineage', label, label, '')


def compute_cell_state_evidence(adata, threshold=0.35):
    """Compute independent, multi-label cell-state modules.

    ``dominant_state`` is a convenience view only.  ``state_flags`` and the
    per-module boolean masks are authoritative for overlapping biology such as
    simultaneous interferon, stress and hypoxia responses.
    """
    import numpy as np

    var_lookup = {str(gene).upper(): str(gene) for gene in adata.var_names}
    scores = {}
    coverage = {}
    expression_source = ''
    for state_name, genes in CELL_STATE_MARKERS.items():
        matched = list(dict.fromkeys(
            var_lookup[str(gene).upper()] for gene in genes
            if str(gene).upper() in var_lookup
        ))
        coverage[state_name] = {'matched': len(matched), 'total': len(genes)}
        values, module_source = _mean_log_normalized_expression_for_genes(
            adata, matched,
        )
        if matched and not expression_source:
            expression_source = module_source
        positive = values[np.isfinite(values) & (values > 0)]
        scale = float(np.nanpercentile(positive, 95)) if positive.size else 0.0
        scores[state_name] = np.clip(values / (scale + 1e-12), 0.0, 1.0) if scale else np.zeros(adata.n_obs)

    score_names = list(scores)
    score_matrix = np.column_stack([scores[name] for name in score_names]) if score_names else np.zeros((adata.n_obs, 0))
    top_index = score_matrix.argmax(axis=1) if score_matrix.shape[1] else np.zeros(adata.n_obs, dtype=int)
    top_score = score_matrix.max(axis=1) if score_matrix.shape[1] else np.zeros(adata.n_obs)
    dominant = np.asarray([
        score_names[index] if value >= float(threshold) else 'baseline'
        for index, value in zip(top_index, top_score)
    ], dtype=object)
    flags = np.asarray([
        ';'.join(name for name in score_names if scores[name][index] >= float(threshold)) or 'baseline'
        for index in range(adata.n_obs)
    ], dtype=object)
    high_flags = {
        name: np.asarray(scores[name] >= float(threshold), dtype=bool)
        for name in score_names
    }
    return {
        'scores': scores,
        'dominant_state': dominant,
        # Backward-compatible alias for older callers.  It is explicitly the
        # dominant state, not the multi-label representation.
        'primary_state': dominant,
        'state_flags': flags,
        'high_flags': high_flags,
        'coverage': coverage,
        'threshold': float(threshold),
        'expression_source': expression_source or 'no matched state genes',
    }


def parse_marker_mapping(text):
    """Parse ``CellType:GENE1,GENE2`` mappings from form/API text."""
    parsed = {}
    for item in str(text or '').replace('\n', ';').split(';'):
        item = item.strip()
        if ':' not in item:
            continue
        cell_type, genes = item.split(':', 1)
        genes = [gene.strip() for gene in genes.split(',') if gene.strip()]
        if cell_type.strip() and genes:
            parsed[cell_type.strip()] = genes
    return parsed


def parse_cluster_merge_mapping(text):
    """Parse explicitly confirmed source-cluster merges.

    Automatic score-profile correlation is not sufficient evidence that two
    clusters are biologically identical.  This parser deliberately accepts
    only a small, reviewable mapping of final label to source-cluster IDs:
    ``Mature absorptive enterocytes=4,10``.  The caller validates that every
    ID exists in the selected cluster column and retains that column intact.
    """
    mapping = {}
    assigned_clusters = {}
    for raw_line in str(text or '').splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if '=' not in line:
            raise ValueError(
                'cluster_merge_map 每行必须使用“合并后标签=ClusterID1,ClusterID2”格式。'
            )
        label, raw_clusters = line.split('=', 1)
        label = label.strip()
        cluster_ids = [value.strip() for value in raw_clusters.split(',') if value.strip()]
        if not label or len(label) > 120:
            raise ValueError('cluster_merge_map 的合并后标签不能为空且最多 120 个字符。')
        if label.lower() == 'unknown':
            raise ValueError('cluster_merge_map 不能把来源 cluster 合并为 Unknown。')
        if len(cluster_ids) < 2:
            raise ValueError('每个 cluster_merge_map 条目至少需要两个来源 cluster。')
        if label in mapping:
            raise ValueError(f'cluster_merge_map 中的标签 “{label}” 重复。')
        duplicates = [
            cluster for cluster in cluster_ids if cluster in assigned_clusters
        ]
        if duplicates:
            raise ValueError(
                '同一来源 cluster 只能出现在一个合并条目中：'
                + ', '.join(sorted(set(duplicates)))
            )
        mapping[label] = cluster_ids
        assigned_clusters.update({cluster: label for cluster in cluster_ids})
    return mapping


def _mean_expression_for_genes(adata, genes, layer=None):
    """Return per-cell mean expression for an already matched gene list."""
    import numpy as np

    if not genes:
        return np.zeros(adata.n_obs, dtype=float)
    view = adata[:, list(dict.fromkeys(genes))]
    expr = view.layers[layer] if layer and layer in view.layers else view.X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray()
    return np.asarray(expr, dtype=float).mean(axis=1).ravel()


def _mean_log_normalized_expression_for_genes(adata, genes, target_sum=10000.0):
    """Return module expression on a library-normalized log1p scale.

    Raw UMI counts are useful for detection evidence but not for comparing
    state-module magnitudes across cells with different library sizes.  When a
    counts layer exists, normalize only the requested genes using the full
    per-cell library size, so this remains sparse and does not copy the whole
    expression matrix.  The function never mutates ``adata``.
    """
    import numpy as np

    if not genes:
        return np.zeros(adata.n_obs, dtype=float), 'no matched genes'

    genes = list(dict.fromkeys(genes))
    # 只有确认是原始整数 counts 的层才能按 counts 归一化；第三方 h5ad
    # 的 counts 层可能是已归一化浮点，再 normalize+log1p 会造成双重变换。
    if 'counts' in adata.layers and _matrix_is_raw_counts(adata.layers['counts']):
        matrix = adata[:, genes].layers['counts'].copy()
        library_size = np.asarray(
            adata.layers['counts'].sum(axis=1), dtype=float
        ).reshape(-1)
        scale = np.divide(
            float(target_sum), library_size,
            out=np.zeros_like(library_size, dtype=float), where=library_size > 0,
        )
        if hasattr(matrix, 'tocsr'):
            matrix = matrix.astype(float).multiply(scale[:, None]).tocsr()
            matrix.data = np.log1p(matrix.data)
            values = np.asarray(matrix.mean(axis=1), dtype=float).reshape(-1)
        else:
            normalized = np.asarray(matrix, dtype=float) * scale[:, None]
            values = np.log1p(normalized).mean(axis=1)
        return values, 'layers["counts"] → library-size normalized log1p'

    source = 'adata.X'
    matrix = adata[:, genes].X
    if adata.raw is not None:
        raw_names = {str(gene) for gene in adata.raw.var_names}
        if all(str(gene) in raw_names for gene in genes):
            raw_matrix = adata.raw[:, genes].X
            raw_values = raw_matrix.data if hasattr(raw_matrix, 'data') else np.asarray(raw_matrix)
            if np.asarray(raw_values).size == 0 or np.nanmin(raw_values) >= 0:
                matrix = raw_matrix
                source = 'adata.raw.X'
    stored = matrix.data if hasattr(matrix, 'data') else np.asarray(matrix)
    stored = np.asarray(stored, dtype=float)
    is_nonnegative = stored.size == 0 or np.nanmin(stored) >= 0
    is_integer = stored.size == 0 or np.allclose(stored, np.rint(stored), atol=1e-6)
    if is_nonnegative and is_integer:
        full = adata.raw.X if source == 'adata.raw.X' else adata.X
        library_size = np.asarray(full.sum(axis=1), dtype=float).reshape(-1)
        scale = np.divide(
            float(target_sum), library_size,
            out=np.zeros_like(library_size, dtype=float), where=library_size > 0,
        )
        if hasattr(matrix, 'tocsr'):
            matrix = matrix.astype(float).multiply(scale[:, None]).tocsr()
            matrix.data = np.log1p(matrix.data)
            values = np.asarray(matrix.mean(axis=1), dtype=float).reshape(-1)
        else:
            values = np.log1p(np.asarray(matrix, dtype=float) * scale[:, None]).mean(axis=1)
        return values, source + ' → library-size normalized log1p'

    if hasattr(matrix, 'toarray'):
        matrix = matrix.toarray()
    values = np.asarray(matrix, dtype=float).mean(axis=1).reshape(-1)
    if source == 'adata.X':
        prepared_source = str(
            adata.uns.get('_annotation_temporary_expression_source', '') or ''
        )
        if prepared_source:
            return values, prepared_source
    return values, source + '（assumed normalized expression）'


def compute_annotation_quality_evidence(
    adata,
    score_cols,
    marker_genes_by_type=None,
    *,
    lineage_mixture_threshold=0.30,
    doublet_threshold=None,
    assigned_labels=None,
    label_families=None,
    ambient_threshold=0.35,
    ambient_prevalence=0.50,
):
    """Collect QC doublet provenance, lineage mixture and ambient evidence.

    Marker-score ambiguity is *not* a doublet caller.  Doublet status is
    inherited from QC doublet detection when available; otherwise it is explicitly
    ``not_evaluated``.  The former annotation heuristic is retained under the
    separate ``lineage_mixture_*`` fields as a review signal.

    Ambient evidence is only raised for markers foreign to the assigned broad
    lineage.  A widespread epithelial marker is therefore expected in an
    epithelial organoid cell and cannot, by itself, be called ambient RNA.
    """
    import numpy as np
    import pandas as pd

    # Backward compatibility for saved API payloads from before marker
    # ambiguity was correctly separated from Doublet evidence.
    if doublet_threshold is not None:
        lineage_mixture_threshold = doublet_threshold

    qc_labels = None
    for column in ('predicted_doublet', 'predicted_doublets'):
        if column in adata.obs.columns:
            qc_labels = adata.obs[column].fillna(False).astype(bool).to_numpy()
            break
    qc_summary = adata.uns.get('qc_doublets') or adata.uns.get('qc_scrublet', {})
    if not isinstance(qc_summary, dict):
        qc_summary = {}
    if qc_labels is not None:
        doublet_status = np.where(
            qc_labels, 'qc_predicted_doublet', 'qc_not_predicted'
        ).astype(object)
        doublet_source = 'qc_predicted_doublet'
    elif qc_summary.get('available') and qc_summary.get('filter_doublets'):
        qc_labels = np.zeros(adata.n_obs, dtype=bool)
        doublet_status = np.repeat('qc_not_predicted', adata.n_obs).astype(object)
        doublet_source = f"qc_{qc_summary.get('method', 'doublet')}_filtered_output"
    else:
        doublet_status = np.repeat('not_evaluated', adata.n_obs).astype(object)
        doublet_source = 'not_available'
    if 'doublet_score' in adata.obs.columns:
        doublet_score = pd.to_numeric(
            adata.obs['doublet_score'], errors='coerce'
        ).to_numpy(dtype=float)
    else:
        doublet_score = np.full(adata.n_obs, np.nan, dtype=float)

    score_cols = [str(col) for col in score_cols if col in adata.obs.columns]
    if not score_cols:
        return {
            'doublet_score': doublet_score,
            'doublet_status': doublet_status,
            'doublet_source': doublet_source,
            'lineage_mixture_score': np.zeros(adata.n_obs),
            'lineage_mixture_status': np.repeat('not_evaluated', adata.n_obs),
            'ambient_score': np.zeros(adata.n_obs),
            'ambient_status': np.repeat('not_evaluated', adata.n_obs),
            'top1': np.repeat('', adata.n_obs), 'top2': np.repeat('', adata.n_obs),
            'top1_score': np.zeros(adata.n_obs), 'top2_score': np.zeros(adata.n_obs),
            'ambient_genes': [], 'ambient_genes_by_label': {},
            'ambient_expression_source': 'not_available',
            'warnings': ['没有可用 marker score，未计算谱系混合/环境 RNA 证据。'],
        }

    score_matrix = np.asarray(adata.obs[score_cols].values, dtype=float)
    score_matrix = np.nan_to_num(score_matrix, nan=0.0, posinf=0.0, neginf=0.0)
    order = np.argsort(score_matrix, axis=1)
    top1_idx = order[:, -1]
    top2_idx = order[:, -2] if len(score_cols) > 1 else top1_idx
    top1_score = score_matrix[np.arange(adata.n_obs), top1_idx]
    top2_score = score_matrix[np.arange(adata.n_obs), top2_idx]
    # Softmax makes the score comparable across panels with different numbers
    # of cell types.  2 * min(p1, p2) is 1 for an exactly balanced two-way call.
    shifted = score_matrix - score_matrix.max(axis=1, keepdims=True)
    exp_scores = np.exp(np.clip(shifted, -30, 30))
    probabilities = exp_scores / (exp_scores.sum(axis=1, keepdims=True) + 1e-12)
    lineage_mixture_score = 2.0 * np.minimum(
        probabilities[np.arange(adata.n_obs), top1_idx],
        probabilities[np.arange(adata.n_obs), top2_idx],
    )
    # Two negative/near-background scores are ambiguity, not evidence of a
    # doublet.  Require both modules to have a positive score before raising
    # the doublet signal.
    positive_pair = (top1_score > 0) & (top2_score > 0)
    lineage_mixture_score = np.where(positive_pair, lineage_mixture_score, 0.0)
    if len(score_cols) < 2:
        lineage_mixture_score = np.zeros(adata.n_obs, dtype=float)
        lineage_mixture_status = np.repeat('not_evaluated', adata.n_obs)
    else:
        lineage_mixture_status = np.where(
            lineage_mixture_score >= float(lineage_mixture_threshold),
            'mixed_lineage_review', 'no_strong_lineage_mixture',
        )

    marker_genes_by_type = marker_genes_by_type or {}
    label_families = {
        str(key): str(value) for key, value in (label_families or {}).items()
    }
    var_lookup = {str(gene).upper(): str(gene) for gene in adata.var_names}
    marker_to_types = {}
    for cell_type, genes in marker_genes_by_type.items():
        for gene in genes or []:
            var_gene = var_lookup.get(str(gene).upper())
            if var_gene:
                marker_to_types.setdefault(var_gene, set()).add(str(cell_type))
    all_markers = sorted({
        var_lookup[str(gene).upper()]
        for genes in marker_genes_by_type.values() for gene in (genes or [])
        if str(gene).upper() in var_lookup
    })
    if all_markers:
        detection_matrix = (
            adata[:, all_markers].layers['counts']
            if 'counts' in adata.layers else adata[:, all_markers].X
        )
        global_fraction = _nonzero_fraction(detection_matrix)
        ambient_expression_source = (
            'layers["counts"] detection' if 'counts' in adata.layers
            else 'adata.X detection'
        )
    else:
        detection_matrix = None
        global_fraction = np.array([])
        ambient_expression_source = 'not_available'
    # A marker shared by several candidate types is not useful as ambient
    # evidence.  Restrict the heuristic to a lineage-specific marker that is
    # detected in many cells but is not part of the winning module.
    ambient_genes = [
        gene for gene, fraction in zip(all_markers, global_fraction)
        if float(fraction) >= float(ambient_prevalence)
        and len(marker_to_types.get(gene, set())) == 1
    ]
    top1_labels = np.asarray(
        [score_cols[i].replace('score_', '') for i in top1_idx], dtype=object,
    )
    if assigned_labels is None:
        assigned_labels = top1_labels
    else:
        assigned_labels = np.asarray(pd.Series(assigned_labels).astype(str), dtype=object)
        if len(assigned_labels) != adata.n_obs:
            raise ValueError('assigned_labels 长度必须与 adata.n_obs 一致。')

    ambient_genes_by_label = {}
    if ambient_genes:
        ambient_matrix = (
            adata[:, ambient_genes].layers['counts']
            if 'counts' in adata.layers else adata[:, ambient_genes].X
        )
        if hasattr(ambient_matrix, 'toarray'):
            ambient_matrix = ambient_matrix.toarray()
        ambient_detected = np.asarray(ambient_matrix, dtype=float) > 0
        ambient_score = np.zeros(adata.n_obs, dtype=float)
        evaluated = np.ones(adata.n_obs, dtype=bool)
        for index, assigned_label in enumerate(assigned_labels):
            if str(assigned_label).lower() in {'', 'unknown', 'nan'}:
                evaluated[index] = False
                continue
            assigned_family = label_families.get(
                str(assigned_label), str(assigned_label)
            )
            eligible = [
                gene_index for gene_index, gene in enumerate(ambient_genes)
                if assigned_family not in {
                    label_families.get(marker_type, marker_type)
                    for marker_type in marker_to_types.get(gene, set())
                }
            ]
            if eligible:
                ambient_score[index] = float(np.mean(ambient_detected[index, eligible]))
                ambient_genes_by_label.setdefault(
                    str(assigned_label), [ambient_genes[i] for i in eligible]
                )
    else:
        ambient_score = np.zeros(adata.n_obs, dtype=float)
        evaluated = np.zeros(adata.n_obs, dtype=bool)
    ambient_status = np.where(
        ~evaluated, 'not_evaluated',
        np.where(
            ambient_score >= float(ambient_threshold),
            'ambient_signal', 'no_strong_ambient_signal',
        ),
    ).astype(object)
    warnings = []
    if ambient_genes:
        warnings.append(
            '环境 RNA 仅按“异源谱系 marker 检出”做启发式筛查；同谱系 marker 不计入。'
            '若需定量去污染，请提供 empty droplets 或运行专门去污染步骤。'
        )
    return {
        'doublet_score': doublet_score,
        'doublet_status': doublet_status,
        'doublet_source': doublet_source,
        'lineage_mixture_score': lineage_mixture_score,
        'lineage_mixture_status': lineage_mixture_status,
        'ambient_score': ambient_score,
        'ambient_status': ambient_status,
        'top1': top1_labels,
        'top2': np.asarray([score_cols[i].replace('score_', '') for i in top2_idx], dtype=object),
        'top1_score': top1_score,
        'top2_score': top2_score,
        'ambient_genes': ambient_genes,
        'ambient_genes_by_label': ambient_genes_by_label,
        'ambient_expression_source': ambient_expression_source,
        'warnings': warnings,
    }


def resolve_marker_genes(markers, var_names, min_markers_per_type=2):
    """Match marker symbols case-insensitively and report usable coverage.

    This supports common human/mouse symbol casing differences and prevents a
    lineage with almost no measurable markers from participating in scoring.
    """
    lookup = {str(gene).upper(): str(gene) for gene in var_names}
    usable, coverage = {}, {}
    for cell_type, genes in markers.items():
        matched = list(dict.fromkeys(lookup[g.upper()] for g in genes if g.upper() in lookup))
        coverage[cell_type] = {'matched': len(matched), 'total': len(genes)}
        if len(matched) >= min_markers_per_type:
            usable[cell_type] = matched
    return usable, coverage


DEFAULT_MARKER_EXCLUDE = {
    'MALAT1', 'ACTB', 'GAPDH', 'B2M', 'UBC', 'PPIA', 'EEF1A1',
    'HSP90AA1', 'HSP90AB1', 'TMSB10', 'TMSB4X', 'RPLP0', 'RPS18',
}


def _nonzero_fraction(matrix):
    """Return per-gene detection fractions for dense or sparse matrices."""
    import numpy as np

    n_rows = max(int(matrix.shape[0]), 1)
    if hasattr(matrix, 'getnnz'):
        return np.asarray(matrix.getnnz(axis=0), dtype=float).ravel() / n_rows
    array = np.asarray(matrix)
    return np.asarray((array > 0).sum(axis=0), dtype=float).ravel() / n_rows


def select_cluster_marker_genes(
    adata,
    cluster_key,
    classic_markers=None,
    *,
    method='wilcoxon',
    n_rank_genes=200,
    min_markers=2,
    max_markers=5,
    padj_cutoff=0.05,
    min_pct=0.10,
    min_delta_pct=0.05,
    anchor_min_pct=0.10,
    exclude_genes=None,
    allow_relaxed_fallback=True,
):
    """Select evidence-backed, cluster-specific marker genes for a dotplot.

    The selector deliberately separates two roles:

    * ``data_driven`` genes pass adjusted-P, positive logFC, detection-rate and
      specificity filters from an independent per-cluster ranking.
    * ``classic_anchor`` genes come from a supplied marker panel and are kept
      only when detected in the cluster, so the figure remains biologically
      interpretable without pretending that a fixed panel is a DEG result.

    Genes are assigned globally so the same marker is not repeatedly displayed
    for several clusters.  If a small dataset has no gene passing the strict
    thresholds, a ranked positive-effect fallback is recorded explicitly in the
    returned metadata instead of being silently presented as significant.
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc

    def _safe_float(value, default):
        try:
            value = float(value)
            return default if not np.isfinite(value) else value
        except (TypeError, ValueError):
            return default

    def _safe_int(value, default):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return default

    method = str(method or 'wilcoxon')
    n_rank_genes = _safe_int(n_rank_genes, 200)
    min_markers = _safe_int(min_markers, 2)
    max_markers = _safe_int(max_markers, 5)
    padj_cutoff = min(1.0, max(0.0, _safe_float(padj_cutoff, 0.05)))
    min_pct = min(1.0, max(0.0, _safe_float(min_pct, 0.10)))
    min_delta_pct = min(1.0, max(0.0, _safe_float(min_delta_pct, 0.05)))
    anchor_min_pct = min(1.0, max(0.0, _safe_float(anchor_min_pct, 0.10)))

    empty = {
        'cluster_markers': {},
        'marker_details': {},
        'genes': [],
        'warnings': [],
        'parameters': {
            'method': method, 'padj_cutoff': float(padj_cutoff),
            'min_pct': float(min_pct), 'min_delta_pct': float(min_delta_pct),
            'min_markers': int(min_markers), 'max_markers': int(max_markers),
            'anchor_min_pct': float(anchor_min_pct),
            'n_rank_genes': int(n_rank_genes),
        },
    }
    if cluster_key not in adata.obs.columns:
        empty['warnings'].append(f"聚类列 '{cluster_key}' 不存在")
        return empty

    labels = adata.obs[cluster_key].astype(str)
    groups = sorted(labels.dropna().unique().tolist(), key=lambda x: (len(x), x))
    if len(groups) < 2:
        empty['warnings'].append('至少需要两个 cluster 才能进行差异 marker 排名')
        return empty

    min_markers = max(1, int(min_markers))
    max_markers = max(min_markers, min(10, int(max_markers)))
    n_rank_genes = max(1, min(int(n_rank_genes), int(adata.n_vars)))
    empty['parameters'].update({
        'min_markers': int(min_markers),
        'max_markers': int(max_markers),
        'n_rank_genes': int(n_rank_genes),
    })
    work = adata.copy()
    work.obs[cluster_key] = labels.astype('category')
    rank_key = '_nature_marker_rank'
    ranked_by_group = {}
    if method == 'logreg':
        # Scanpy's logistic-regression backend stores a single binary class
        # per call.  Run one-vs-rest explicitly so every cluster receives a
        # result instead of silently dropping all but the first category.
        empty['warnings'].append('logreg 不提供 adjusted P 值，候选将标记为 ranked_relaxed')
        for group in groups:
            binary = work.copy()
            binary.obs['_nature_marker_binary'] = pd.Categorical(
                np.where(labels.values == str(group), str(group), '__rest__'),
                # Scanpy's logreg result keeps the first categorical field;
                # placing rest first makes its coefficient correspond to the
                # target-vs-rest effect rather than the inverse sign.
                categories=['__rest__', str(group)],
            )
            try:
                sc.tl.rank_genes_groups(
                    binary, groupby='_nature_marker_binary', method='logreg',
                    n_genes=n_rank_genes, key_added=rank_key, use_raw=False,
                )
                result_field = binary.uns[rank_key]['names'].dtype.names[0]
                ranked_by_group[str(group)] = sc.get.rank_genes_groups_df(
                    binary, group=result_field, key=rank_key,
                )
            except Exception as exc:
                empty['warnings'].append(f'cluster {group} logreg 排名失败：{exc}')
                ranked_by_group[str(group)] = pd.DataFrame()
    else:
        try:
            sc.tl.rank_genes_groups(
                work, groupby=cluster_key, method=str(method),
                n_genes=n_rank_genes, key_added=rank_key, use_raw=False,
            )
        except Exception as exc:
            empty['warnings'].append(f'cluster marker 排名失败：{exc}')
            return empty

    gene_names = [str(gene) for gene in work.var_names]
    gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
    exclude = {str(gene).upper() for gene in (exclude_genes or [])}
    exclude.update(DEFAULT_MARKER_EXCLUDE)

    def _is_excluded(gene):
        upper = str(gene).upper()
        return (
            upper in exclude
            or upper.startswith(('MT-', 'MT.'))
            or upper.startswith('RPS')
            or upper.startswith('RPL')
        )

    anchor_set = set()
    if isinstance(classic_markers, dict):
        for genes in classic_markers.values():
            anchor_set.update(str(gene).upper() for gene in genes)
    elif classic_markers:
        anchor_set.update(str(gene).upper() for gene in classic_markers)
    anchor_set -= exclude

    group_candidates = {}
    result_groups = list(work.obs[cluster_key].cat.categories.astype(str))
    for group in result_groups:
        group_mask = labels.values == str(group)
        rest_mask = ~group_mask
        pct_in = _nonzero_fraction(work[group_mask].X)
        pct_out = _nonzero_fraction(work[rest_mask].X) if rest_mask.any() else np.zeros(work.n_vars)
        try:
            ranked = (
                ranked_by_group.get(str(group), pd.DataFrame())
                if method == 'logreg'
                else sc.get.rank_genes_groups_df(work, group=group, key=rank_key)
            )
        except Exception as exc:
            empty['warnings'].append(f"cluster {group} marker 提取失败：{exc}")
            group_candidates[str(group)] = []
            continue

        candidates = {}
        for _, row in ranked.iterrows():
            gene = str(row.get('names', ''))
            if gene not in gene_to_idx or _is_excluded(gene):
                continue
            idx = gene_to_idx[gene]
            pct_group = float(pct_in[idx])
            pct_rest = float(pct_out[idx])
            delta_pct = pct_group - pct_rest
            raw_logfc = float(pd.to_numeric(row.get('logfoldchanges', np.nan), errors='coerce'))
            padj = float(pd.to_numeric(row.get('pvals_adj', np.nan), errors='coerce'))
            score = float(pd.to_numeric(row.get('scores', np.nan), errors='coerce'))
            logfc = raw_logfc
            if not np.isfinite(logfc) and method == 'logreg':
                # For logreg, score is the available signed effect ranking.
                logfc = score
            if not np.isfinite(logfc):
                logfc = 0.0
            if not np.isfinite(padj):
                padj = 1.0
            if not np.isfinite(score):
                score = 0.0
            effect_score = (
                max(logfc, 0.0) * max(delta_pct, 0.0)
                * max(-np.log10(max(padj, 1e-300)), 0.0)
            )
            is_anchor = (
                gene.upper() in anchor_set
                and pct_group >= float(anchor_min_pct)
                and (logfc > 0 or delta_pct > 0)
            )
            strict = (
                padj <= float(padj_cutoff)
                and logfc > 0
                and pct_group >= float(min_pct)
                and delta_pct >= float(min_delta_pct)
            )
            if strict or is_anchor:
                candidates[gene] = {
                    'gene': gene, 'score': effect_score,
                    'logfc': logfc, 'padj': padj,
                    'pct_in': pct_group, 'pct_out': pct_rest,
                    'delta_pct': delta_pct,
                    'source': 'classic_anchor' if is_anchor else 'data_driven',
                    'strict': bool(strict),
                }

        strict_data = [item for item in candidates.values() if item['source'] == 'data_driven' and item['strict']]
        if not strict_data and allow_relaxed_fallback:
            for _, row in ranked.iterrows():
                gene = str(row.get('names', ''))
                if gene not in gene_to_idx or _is_excluded(gene) or gene in candidates:
                    continue
                idx = gene_to_idx[gene]
                pct_group = float(pct_in[idx])
                pct_rest = float(pct_out[idx])
                raw_logfc = float(pd.to_numeric(row.get('logfoldchanges', np.nan), errors='coerce'))
                logfc = raw_logfc
                if not np.isfinite(logfc) and method == 'logreg':
                    logfc = float(pd.to_numeric(row.get('scores', np.nan), errors='coerce'))
                if not np.isfinite(logfc) or logfc <= 0 or pct_group < float(min_pct):
                    continue
                delta_pct = pct_group - pct_rest
                candidates[gene] = {
                    'gene': gene, 'score': max(logfc, 0.0) * max(delta_pct, 0.0),
                    'logfc': logfc, 'padj': float(pd.to_numeric(row.get('pvals_adj', 1.0), errors='coerce')),
                    'pct_in': pct_group, 'pct_out': pct_rest,
                    'delta_pct': delta_pct, 'source': 'ranked_relaxed', 'strict': False,
                }
            if any(item['source'] == 'ranked_relaxed' for item in candidates.values()):
                empty['warnings'].append(
                    f'cluster {group} 没有基因同时满足严格阈值，已标记 ranked_relaxed 作为探索性兜底'
                )
        group_candidates[str(group)] = list(candidates.values())

    selected = {str(group): [] for group in result_groups}
    used_genes = set()
    # Reserve at most one classic marker per cluster, globally unique.
    anchor_pairs = sorted(
        ((group, item) for group, items in group_candidates.items() for item in items if item['source'] == 'classic_anchor'),
        key=lambda pair: pair[1]['score'], reverse=True,
    )
    for group, item in anchor_pairs:
        if selected[group] or item['gene'] in used_genes:
            continue
        selected[group].append(item)
        used_genes.add(item['gene'])

    def _candidate_priority(item):
        source_bonus = {
            'data_driven': 1.0,
            'classic_anchor': 0.8,
            'ranked_relaxed': 0.3,
        }.get(item['source'], 0)
        return float(item['score']) + source_bonus

    # First give every cluster a fair chance to reach the requested minimum.
    # A single cluster with many very strong genes must not consume the whole
    # panel before a small cluster receives any evidence.
    for group in selected:
        while len(selected[group]) < min_markers:
            available = [
                item for item in group_candidates.get(group, [])
                if item['gene'] not in used_genes
            ]
            if not available:
                break
            item = max(available, key=_candidate_priority)
            selected[group].append(item)
            used_genes.add(item['gene'])

    # Fill remaining slots globally by evidence strength, while respecting the
    # per-cluster maximum and preserving the global no-duplicate guarantee.
    while True:
        pairs = []
        for group, items in group_candidates.items():
            if len(selected[group]) >= max_markers:
                continue
            for item in items:
                if item['gene'] not in used_genes:
                    pairs.append((_candidate_priority(item), group, item))
        if not pairs:
            break
        _, group, item = max(pairs, key=lambda pair: pair[0])
        selected[group].append(item)
        used_genes.add(item['gene'])

    for group in selected:
        selected[group].sort(key=lambda item: (item['source'] != 'classic_anchor', -item['score'], item['gene']))
        if len(selected[group]) < min_markers:
            empty['warnings'].append(
                f'cluster {group} 最终仅保留 {len(selected[group])} 个 marker，低于目标 {min_markers} 个'
            )
    empty['cluster_markers'] = {group: [item['gene'] for item in items] for group, items in selected.items()}
    empty['marker_details'] = selected
    empty['genes'] = list(dict.fromkeys(gene for items in empty['cluster_markers'].values() for gene in items))
    empty['n_data_driven'] = sum(item['source'] == 'data_driven' for items in selected.values() for item in items)
    empty['n_classic_anchor'] = sum(item['source'] == 'classic_anchor' for items in selected.values() for item in items)
    empty['n_relaxed'] = sum(item['source'] == 'ranked_relaxed' for items in selected.values() for item in items)
    return empty


def build_marker_validation_matrix(adata, genes, groupby='celltype', max_cells=5000):
    """Prepare mean expression and detection fractions for annotation review.

    Annotation often receives an AnnData whose ``X`` has already been scaled for
    PCA.  Using that matrix for an expression plot produces negative values and
    makes the plot look like a collection of arbitrary score boxes.  Prefer the
    raw UMI layer and create a display-only library-size-normalized log1p matrix;
    this never mutates ``adata`` or changes the annotation scores.
    """
    import numpy as np
    import pandas as pd

    if groupby not in adata.obs.columns:
        return None
    available = {str(gene): str(gene) for gene in adata.var_names}
    selected_genes = []
    for gene in genes or []:
        actual = available.get(str(gene))
        if actual and actual not in selected_genes:
            selected_genes.append(actual)
    if not selected_genes:
        return None

    # Plotting all cells is exact for ordinary projects but needlessly inflates
    # a static figure for very large atlases.  A deterministic subset keeps the
    # preview reproducible while preserving the original cell labels.
    n_obs = int(adata.n_obs)
    indices = np.arange(n_obs)
    if n_obs > int(max_cells):
        rng = np.random.default_rng(0)
        indices = np.sort(rng.choice(n_obs, int(max_cells), replace=False))

    def _dense(matrix):
        if hasattr(matrix, 'toarray'):
            matrix = matrix.toarray()
        return np.asarray(matrix, dtype=float)

    source = ''
    display_warning = ''
    gene_indices = [available[gene] for gene in selected_genes]
    if 'counts' in adata.layers:
        count_matrix = _dense(adata[indices, gene_indices].layers['counts'])
        library_size = _dense(adata[indices].layers['counts'].sum(axis=1)).ravel()
        normalized = np.divide(
            count_matrix, library_size[:, None],
            out=np.zeros_like(count_matrix), where=library_size[:, None] > 0,
        ) * 10000.0
        display_matrix = np.log1p(normalized)
        source = 'layers["counts"] → library-size normalized log1p'
    else:
        # Fall back to raw/normalized X only when a counts layer is unavailable.
        # The fallback is labelled so reviewers know that the display is not
        # based on a recoverable UMI layer.
        matrix = None
        if adata.raw is not None:
            raw_lookup = {str(gene): index for index, gene in enumerate(adata.raw.var_names)}
            if all(gene in raw_lookup for gene in selected_genes):
                matrix = _dense(adata.raw[indices, selected_genes].X)
                source = 'adata.raw.X'
        if matrix is None:
            matrix = _dense(adata[indices, gene_indices].X)
            source = 'adata.X'
        if np.nanmin(matrix) < 0:
            display_warning = '未找到 counts 层；当前表达矩阵含负值，图中按非负部分显示。'
            matrix = np.clip(matrix, 0.0, None)
            source += '（scaled，已截断负值）'
        elif np.allclose(matrix, np.rint(matrix)):
            library_size = matrix.sum(axis=1)
            normalized = np.divide(
                matrix, library_size[:, None],
                out=np.zeros_like(matrix), where=library_size[:, None] > 0,
            ) * 10000.0
            matrix = np.log1p(normalized)
            source += ' → library-size normalized log1p'
        display_matrix = matrix

    labels = pd.Series(adata.obs[groupby].astype(str).iloc[indices].values)
    category_order = labels.value_counts().index.tolist()
    if 'Unknown' in category_order:
        category_order = [item for item in category_order if item != 'Unknown'] + ['Unknown']
    means = []
    detections = []
    for category in category_order:
        mask = labels.to_numpy() == category
        group_values = display_matrix[mask]
        means.append(np.nanmean(group_values, axis=0))
        detections.append(np.mean(group_values > 0, axis=0))
    return {
        'mean_expression': np.asarray(means, dtype=float),
        'detection_fraction': np.asarray(detections, dtype=float),
        'categories': category_order,
        'genes': selected_genes,
        'expression_source': source,
        'display_warning': display_warning,
        'n_cells_plotted': int(len(indices)),
    }


def build_marker_decisions(
    adata,
    cluster_key,
    markers,
    *,
    marker_set_name='Universal',
    multi_evidence=False,
    agreement_threshold=0.60,
    min_annotation_score=0.0,
    marker_min_pct=0.10,
    marker_min_delta_pct=0.05,
    broad_markers=None,
    subtype_markers=None,
    subtype_lineage_map=None,
    cluster_marker_selection=None,
    annotate_all=True,
):
    """Build hierarchical per-cell predictions and auditable cluster calls.

    Broad lineage is selected first and, when a subtype panel is supplied
    (Colorectal epithelium, or Universal fine mode), the subtype is decided in
    a second pass.  The broad label itself stays in the subtype candidate
    list, so a cluster without decisive subtype evidence falls back to the
    lineage label instead of being forced into the nearest subtype.

    annotate_all implements the "尽量注释所有细胞" rule: a low per-cell
    vote alone no longer erases a cluster label.  Without specific marker
    support the label is kept for review (reason
    low_cell_agreement_label_kept_for_review); only a total absence of
    positive evidence (all candidate scores <= 0), an explicit
    min_annotation_score gate, or the legacy strict mode still produce
    Unknown.
    """
    import numpy as np
    import pandas as pd

    marker_set_name = normalize_marker_set_name(marker_set_name)
    # 细标签优先的判定规则用于 Colorectal refined 与 Organoid 面板：cluster 的
    # 独立排名 marker 优先于共享分数模块。锚点来源：精细面板用人工确认的锚点表，
    # Organoid 等面板用“只属于单个 programme”的自动派生锚点。
    refined_like = marker_set_name in {'Colorectal_refined', 'Organoid'}
    if refined_like:
        panel_anchor_markers = (
            COLORECTAL_REFINED_ANCHOR_MARKERS
            if marker_set_name == 'Colorectal_refined'
            else derive_panel_anchor_markers(markers)
        )
    else:
        panel_anchor_markers = {}
    # UMAP/KNN 邻接只进复核表，绝不参与自动合并。
    neighborhood = compute_cluster_neighborhood_evidence(adata, cluster_key)

    score_types = [
        cell_type for cell_type in markers
        if f'score_{cell_type}' in adata.obs.columns
    ]
    if not score_types:
        return {
            'cluster_annotations': {}, 'cluster_evidence': {},
            'prediction_1': np.repeat('', adata.n_obs),
            'prediction_2': np.repeat('', adata.n_obs),
            'prediction_1_score': np.zeros(adata.n_obs),
            'prediction_2_score': np.zeros(adata.n_obs),
            'broad_prediction': np.repeat('', adata.n_obs),
        }

    def rank_frame(candidate_types):
        candidate_types = [item for item in candidate_types if item in score_types]
        if not candidate_types:
            return None
        columns = [f'score_{item}' for item in candidate_types]
        values = np.nan_to_num(
            adata.obs[columns].to_numpy(dtype=float), nan=-np.inf,
        )
        order = np.argsort(values, axis=1)
        top = order[:, -1]
        second = order[:, -2] if len(columns) > 1 else top
        row = np.arange(adata.n_obs)
        return {
            'types': candidate_types,
            'columns': columns,
            'values': values,
            'top_label': np.asarray([candidate_types[i] for i in top], dtype=object),
            'second_label': np.asarray([candidate_types[i] for i in second], dtype=object),
            'top_score': values[row, top],
            'second_score': values[row, second],
        }

    all_rank = rank_frame(score_types)
    broad_rank = None
    subtype_rank = None
    subtype_groups = {}
    if marker_set_name == 'Colorectal':
        broad_rank = rank_frame(COLORECTAL_LINEAGE_MARKERS)
        subtype_rank = rank_frame(COLORECTAL_EPITHELIAL_SUBTYPE_MARKERS)
        subtype_lineage_map = subtype_lineage_map or {
            subtype: 'Epithelial' for subtype in COLORECTAL_EPITHELIAL_SUBTYPE_MARKERS
        }
    elif broad_markers is not None and subtype_markers is not None:
        broad_rank = rank_frame(list(broad_markers))
        subtype_rank = rank_frame(list(subtype_markers))
    if subtype_rank is not None:
        lineage_map = subtype_lineage_map or {}
        for subtype in subtype_rank['types']:
            parent = lineage_map.get(subtype, '')
            if parent:
                subtype_groups.setdefault(parent, []).append(subtype)

    if broad_rank is not None:
        prediction_1 = broad_rank['top_label'].copy()
        prediction_2 = broad_rank['second_label'].copy()
        prediction_1_score = broad_rank['top_score'].copy()
        prediction_2_score = broad_rank['second_score'].copy()
        broad_prediction = broad_rank['top_label'].copy()
        if subtype_rank is not None:
            # Per-cell second pass, restricted to the broad prediction's own
            # subtype group.  The broad label is appended as an explicit
            # fallback so weak subtype evidence keeps the lineage label.
            for parent, subtypes in subtype_groups.items():
                if parent not in broad_rank['types']:
                    continue
                mask = broad_prediction == parent
                if not bool(mask.any()):
                    continue
                candidates = subtypes + [parent]
                columns = [f'score_{item}' for item in candidates]
                values = np.nan_to_num(
                    adata.obs.loc[mask, columns].to_numpy(dtype=float), nan=-np.inf,
                )
                order = np.argsort(values, axis=1)
                row = np.arange(int(mask.sum()))
                top = order[:, -1]
                second = order[:, -2] if len(candidates) > 1 else top
                prediction_1[mask] = np.asarray(
                    [candidates[i] for i in top], dtype=object,
                )
                prediction_2[mask] = np.asarray(
                    [candidates[i] for i in second], dtype=object,
                )
                prediction_1_score[mask] = values[row, top]
                prediction_2_score[mask] = values[row, second]
    else:
        prediction_1 = all_rank['top_label']
        prediction_2 = all_rank['second_label']
        prediction_1_score = all_rank['top_score']
        prediction_2_score = all_rank['second_score']
        broad_prediction = prediction_1.copy()

    labels = adata.obs[cluster_key].astype(str)
    # ``score_genes`` is useful for broad lineages but its control-gene
    # subtraction makes raw scores from compact, heterogeneous state modules
    # non-comparable.  A local refined panel therefore gives priority to the
    # cluster's independently ranked, data-driven markers.  This is not a
    # hard-coded cluster-ID map: the current cluster has to reproduce at least
    # two signature genes before the refined label can override score ranking.
    refined_cluster_markers = {}
    if refined_like:
        refined_cluster_markers = (
            (cluster_marker_selection or {}).get('cluster_markers', {}) or {}
        )

    def _signature_overlap(cluster, cell_type):
        selected = refined_cluster_markers.get(str(cluster), [])
        selected_by_upper = {
            str(gene).strip().upper(): str(gene).strip()
            for gene in selected if str(gene).strip()
        }
        matched = []
        for gene in markers.get(cell_type, []):
            matched_gene = selected_by_upper.get(str(gene).strip().upper())
            if matched_gene and matched_gene not in matched:
                matched.append(matched_gene)
        return matched

    cluster_annotations = {}
    cluster_evidence = {}
    expression = adata.layers['counts'] if 'counts' in adata.layers else adata.X
    var_positions = {str(gene): index for index, gene in enumerate(adata.var_names)}

    for cluster in sorted(labels.unique().tolist(), key=lambda value: (len(value), value)):
        mask = labels.to_numpy() == str(cluster)
        if broad_rank is not None:
            broad_means = {
                cell_type: float(np.nanmean(adata.obs.loc[mask, f'score_{cell_type}']))
                for cell_type in broad_rank['types']
            }
            broad_candidate = max(broad_means, key=broad_means.get)
            group_subtypes = subtype_groups.get(broad_candidate, [])
            if group_subtypes and subtype_rank is not None:
                # 第二级亚型判定只能在同一层级的亚型之间比较。大谱系
                # marker（EPCAM/KRT8 等）天然比亚型模块更广，若直接把
                # ``Epithelial`` 加回排序候选，会系统性压过任何真实亚型，
                # 让 Colorectal 结果退化为全部 Epithelial。仅当全部亚型
                # 分数都不为正（没有相对背景的亚型证据）才回退到大谱系。
                subtype_means = {
                    cell_type: float(np.nanmean(adata.obs.loc[mask, f'score_{cell_type}']))
                    for cell_type in group_subtypes
                }
                if subtype_means and max(subtype_means.values()) > 0:
                    candidate_means = subtype_means
                else:
                    candidate_means = {broad_candidate: broad_means[broad_candidate]}
            else:
                candidate_means = broad_means
        else:
            broad_candidate = ''
            candidate_means = {
                cell_type: float(np.nanmean(adata.obs.loc[mask, f'score_{cell_type}']))
                for cell_type in score_types
            }

        signature_overlaps = {
            cell_type: _signature_overlap(cluster, cell_type)
            for cell_type in candidate_means
        }
        signature_anchor_overlaps = {
            cell_type: [
                gene for gene in signature_overlaps[cell_type]
                if str(gene).strip().upper() in {
                    str(anchor).strip().upper()
                    for anchor in panel_anchor_markers.get(cell_type, set())
                }
            ]
            for cell_type in candidate_means
        }
        max_signature_overlap = max(
            (len(values) for values in signature_overlaps.values()), default=0,
        )
        has_refined_anchor = any(signature_anchor_overlaps.values())
        use_refined_signature = (
            refined_like
            and (max_signature_overlap >= 2 or has_refined_anchor)
        )
        if use_refined_signature:
            # First compare independently selected signature markers.  A
            # single predefined anchor (for example CYP3A5 or HSD17B2) is
            # allowed to tie a two-marker signature, but cannot outrank a
            # three-marker programme.  This prevents shared score modules
            # from flattening real local states while remaining conservative
            # when only a weak single-marker hint is present.
            ordered = sorted(
                candidate_means,
                key=lambda cell_type: (
                    max(
                        len(signature_overlaps[cell_type]),
                        2 if signature_anchor_overlaps[cell_type] else 0,
                    ),
                    bool(signature_anchor_overlaps[cell_type]),
                    len(signature_overlaps[cell_type]),
                    len(signature_overlaps[cell_type]) / max(
                        len(markers.get(cell_type, [])), 1,
                    ),
                    candidate_means[cell_type],
                ),
                reverse=True,
            )
        else:
            ordered = sorted(candidate_means, key=candidate_means.get, reverse=True)
        candidate = ordered[0]
        second_candidate = ordered[1] if len(ordered) > 1 else ''
        candidate_score = float(candidate_means[candidate])
        second_score = float(candidate_means.get(second_candidate, np.nan))
        score_margin = (
            candidate_score - second_score if np.isfinite(second_score) else np.nan
        )
        agreement = float(np.mean(prediction_1[mask] == candidate))

        candidate_genes = [
            gene for gene in markers.get(candidate, []) if gene in var_positions
        ]
        supported_markers = []
        missing_markers = []
        mean_detection = 0.0
        mean_specificity = 0.0
        if candidate_genes:
            positions = [var_positions[gene] for gene in candidate_genes]
            target = expression[mask][:, positions]
            rest = expression[~mask][:, positions]
            target_pct = _nonzero_fraction(target)
            rest_pct = _nonzero_fraction(rest) if int((~mask).sum()) else np.zeros(len(positions))
            delta = target_pct - rest_pct
            for gene, target_value, delta_value in zip(candidate_genes, target_pct, delta):
                if (
                    float(target_value) >= float(marker_min_pct)
                    and float(delta_value) >= float(marker_min_delta_pct)
                ):
                    supported_markers.append(gene)
                else:
                    missing_markers.append(gene)
            mean_detection = float(np.mean(target_pct))
            mean_specificity = float(np.mean(delta))
        strong_marker_support = len(supported_markers) >= 2

        raw_score_column = f'raw_marker_{candidate}'
        negative_score_column = f'negative_score_{candidate}'
        positive_score = (
            float(pd.to_numeric(adata.obs.loc[mask, raw_score_column], errors='coerce').mean())
            if raw_score_column in adata.obs.columns else np.nan
        )
        negative_marker_score = (
            float(pd.to_numeric(adata.obs.loc[mask, negative_score_column], errors='coerce').mean())
            if negative_score_column in adata.obs.columns else 0.0
        )
        negative_penalty = (
            max(float(positive_score - candidate_score), 0.0)
            if np.isfinite(positive_score) else 0.0
        )

        provisional_refined_label = ''
        if (
            refined_like
            and not use_refined_signature
            and not strong_marker_support
        ):
            # Do not silently assign multiple unsupported clusters to the
            # same broad programme.  Preserve the source cluster and show the
            # independently ranked genes as an honest, reviewable provisional
            # identity.  A later explicit merge can still use cluster_merge_map.
            source_signature = [
                str(gene).strip() for gene in refined_cluster_markers.get(str(cluster), [])
                if str(gene).strip()
            ][:2]
            if source_signature:
                prefix = (
                    'Unresolved epithelial programme: '
                    if marker_set_name == 'Colorectal_refined'
                    else 'Unresolved organoid programme: '
                )
                provisional_refined_label = (
                    prefix
                    + '/'.join(source_signature)
                    + f' (cluster {cluster}; review)'
                )

        final_label = candidate
        if min_annotation_score > 0 and candidate_score < min_annotation_score:
            final_label = 'Unknown'
            decision_reason = 'below_min_annotation_score'
        elif provisional_refined_label:
            final_label = provisional_refined_label
            decision_reason = 'refined_unresolved_cluster_programme_review'
        elif (not np.isfinite(candidate_score) or candidate_score <= 0) and mean_detection <= 0:
            # 完全没有正向 Marker 证据（候选分数不高于背景且候选 Marker
            # 几乎未检出）。为了“一开始就尽量注释完成”，只要逐细胞投票
            # 高度一致且不是平分（margin > 0），就保留该标签并标记为
            # review；只有投票也混乱时才保留 Unknown。负向 Marker 惩罚把
            # 分数压到 0 附近但 Marker 真实检出的簇不会进入此分支。
            if (
                annotate_all
                and np.isfinite(score_margin) and float(score_margin) > 0
                and agreement >= float(agreement_threshold)
            ):
                decision_reason = 'unknown_recovered_by_cell_vote'
            else:
                final_label = 'Unknown'
                decision_reason = 'no_positive_marker_evidence'
        elif multi_evidence and agreement < float(agreement_threshold):
            if strong_marker_support:
                decision_reason = 'low_cell_agreement_but_cluster_markers_support_candidate'
            elif annotate_all:
                # 尽量注释所有细胞：一致率低但没有特异性 Marker 支持时，
                # 保留最佳候选标签并标记为 review，而不是直接抹成 Unknown。
                decision_reason = 'low_cell_agreement_label_kept_for_review'
            else:
                final_label = 'Unknown'
                decision_reason = 'low_cell_agreement_without_specific_marker_support'
        else:
            decision_reason = (
                'hierarchical_colorectal_marker_decision'
                if marker_set_name == 'Colorectal'
                else 'hierarchical_marker_decision'
                if broad_rank is not None and subtype_rank is not None
                else 'cluster_mean_marker_decision'
            )

        resolved_signature_overlap = signature_overlaps.get(candidate, [])
        resolved_anchor_overlap = signature_anchor_overlaps.get(candidate, [])
        negative_conflict = bool(negative_penalty > 0)
        if final_label == 'Unknown' or provisional_refined_label:
            evidence_tier = 'unresolved-review'
        elif refined_like and use_refined_signature:
            if len(resolved_signature_overlap) >= 2:
                evidence_tier = 'confirmed'
            elif resolved_anchor_overlap:
                evidence_tier = 'anchor-supported'
            else:
                evidence_tier = 'provisional'
        elif strong_marker_support:
            evidence_tier = 'confirmed'
        elif decision_reason in (
            'low_cell_agreement_label_kept_for_review',
            'unknown_recovered_by_cell_vote',
        ):
            evidence_tier = 'provisional'
        else:
            evidence_tier = 'unresolved-review' if final_label == 'Unknown' else 'provisional'
        cluster_annotations[str(cluster)] = final_label
        cluster_neighborhood = neighborhood.get(str(cluster), {}) or {}
        cluster_evidence[str(cluster)] = {
            'broad_candidate': broad_candidate or candidate,
            'top_candidate': candidate,
            'second_candidate': second_candidate,
            'positive_marker_score': positive_score,
            'negative_marker_score': negative_marker_score,
            'negative_marker_penalty': negative_penalty,
            'negative_marker_conflict': negative_conflict,
            'adjusted_candidate_score': candidate_score,
            'second_candidate_score': second_score,
            'score_margin': score_margin,
            'cell_prediction_agreement': agreement,
            'candidate_marker_detection_fraction': mean_detection,
            'candidate_marker_specificity_delta': mean_specificity,
            'detected_positive_markers': supported_markers,
            'missing_positive_markers': missing_markers,
            'strong_marker_support': bool(strong_marker_support),
            'refined_signature_overlap': resolved_signature_overlap,
            'refined_signature_overlap_count': len(resolved_signature_overlap),
            'refined_signature_anchor_overlap': resolved_anchor_overlap,
            'refined_signature_anchor_mode_applied': bool(resolved_anchor_overlap),
            'refined_signature_mode_applied': bool(use_refined_signature),
            'provisional_cluster_label': provisional_refined_label,
            'evidence_tier': evidence_tier,
            'umap_adjacent_clusters': list(
                cluster_neighborhood.get('umap_adjacent_clusters', [])
            ),
            'knn_adjacent_clusters': list(
                cluster_neighborhood.get('knn_adjacent_clusters', [])
            ),
            'final_label': final_label,
            'decision_reason': (
                'refined_cluster_marker_signature_decision'
                if use_refined_signature and decision_reason == 'cluster_mean_marker_decision'
                else decision_reason
            ),
        }

    return {
        'cluster_annotations': cluster_annotations,
        'cluster_evidence': cluster_evidence,
        'prediction_1': prediction_1,
        'prediction_2': prediction_2,
        'prediction_1_score': prediction_1_score,
        'prediction_2_score': prediction_2_score,
        'broad_prediction': broad_prediction,
    }


def build_annotation_cluster_review(
    adata, cluster_key, marker_selection=None, neighborhood_evidence=None,
):
    """Create a compact, auditable cluster-level annotation review table.

    ``neighborhood_evidence`` 来自 compute_cluster_neighborhood_evidence，
    只作为人工复核参考（哪个来源 cluster 在 UMAP/KNN 上相邻），绝不参与自动
    合并。
    """
    import numpy as np
    import pandas as pd

    if cluster_key not in adata.obs.columns:
        return pd.DataFrame()

    marker_map = (marker_selection or {}).get('cluster_markers', {}) or {}

    def mode_or_unknown(values):
        values = values.dropna().astype(str)
        if values.empty:
            return 'Unknown'
        return str(values.value_counts().index[0])

    rows = []
    labels = adata.obs[cluster_key].astype(str)
    for cluster in sorted(labels.unique().tolist(), key=lambda value: (len(value), value)):
        mask = labels == cluster
        frame = adata.obs.loc[mask]
        celltype = mode_or_unknown(frame['celltype']) if 'celltype' in frame else 'Unknown'
        final_label = mode_or_unknown(frame['final_annotation']) if 'final_annotation' in frame else celltype
        annotation_status = mode_or_unknown(frame['annotation_status']) if 'annotation_status' in frame else (
            'Unknown' if celltype == 'Unknown' else 'review'
        )
        celltype_values = frame['celltype'].astype(str) if 'celltype' in frame else pd.Series([], dtype=str)
        if 'cluster_annotation_agreement' in frame.columns:
            label_agreement = float(pd.to_numeric(
                frame['cluster_annotation_agreement'], errors='coerce'
            ).mean())
        else:
            label_agreement = float((celltype_values == celltype).mean()) if len(celltype_values) else np.nan
        row = {
            'cluster': str(cluster),
            'n_cells': int(mask.sum()),
            'assigned_celltype': celltype,
            'final_annotation': final_label,
            'marker_programme': (
                mode_or_unknown(frame['annotation_marker_programme'])
                if 'annotation_marker_programme' in frame else celltype
            ),
            'display_label_source': (
                mode_or_unknown(frame['annotation_display_label_source'])
                if 'annotation_display_label_source' in frame else 'marker_programme'
            ),
            'annotation_status': annotation_status,
            'evidence_tier': (
                mode_or_unknown(frame['annotation_evidence_tier'])
                if 'annotation_evidence_tier' in frame else ''
            ),
            'cluster_label_agreement': round(label_agreement, 4) if np.isfinite(label_agreement) else np.nan,
            'top_markers': ', '.join(marker_map.get(str(cluster), [])),
        }
        for source, target in (
            ('annotation_broad_candidate', 'broad_candidate'),
            ('annotation_top_candidate', 'top_candidate'),
            ('annotation_second_candidate', 'second_candidate'),
            ('annotation_detected_positive_markers', 'detected_positive_markers'),
            ('annotation_missing_positive_markers', 'missing_positive_markers'),
            ('annotation_decision_reason', 'final_decision_reason'),
            ('annotation_refined_signature_overlap', 'refined_signature_overlap'),
            ('annotation_refined_signature_anchor_overlap', 'refined_anchor_overlap'),
            ('annotation_negative_marker_conflict', 'negative_marker_conflict'),
        ):
            row[target] = mode_or_unknown(frame[source]) if source in frame.columns else ''
        neighborhood = (neighborhood_evidence or {}).get(str(cluster), {}) or {}
        row['umap_adjacent_clusters'] = ', '.join(
            str(gene) for gene in neighborhood.get('umap_adjacent_clusters', [])
        )
        row['knn_adjacent_clusters'] = ', '.join(
            str(gene) for gene in neighborhood.get('knn_adjacent_clusters', [])
        )
        for source, target in (
            ('annotation_candidate_score', 'top_candidate_score'),
            ('annotation_second_candidate_score', 'second_candidate_score'),
            ('annotation_cluster_score_margin', 'candidate_score_margin'),
            ('annotation_positive_marker_score', 'positive_marker_score'),
            ('annotation_negative_marker_penalty', 'negative_marker_penalty'),
            ('annotation_candidate_marker_detection_fraction', 'candidate_marker_detection_fraction'),
            ('annotation_candidate_marker_specificity_delta', 'candidate_marker_specificity_delta'),
        ):
            row[target] = (
                round(float(pd.to_numeric(frame[source], errors='coerce').mean()), 4)
                if source in frame.columns else np.nan
            )
        if 'celltypist_label' in frame.columns:
            celltypist_labels = frame['celltypist_label'].astype(str)
            celltypist_status = frame.get(
                'celltypist_status', pd.Series('', index=frame.index)
            ).astype(str)
            celltypist_comparison = frame.get(
                'celltypist_comparison', pd.Series('', index=frame.index)
            ).astype(str)
            row['celltypist_reference_label'] = mode_or_unknown(celltypist_labels)
            row['celltypist_conflict_fraction'] = round(
                float((celltypist_comparison == 'conflict').mean()), 4
            )
            row['celltypist_unassigned_fraction'] = round(
                float((celltypist_status == 'unassigned').mean()), 4
            )
            row['celltypist_multi_label_fraction'] = round(
                float((celltypist_status == 'multi_label').mean()), 4
            )
        if 'llm_annotation' in frame.columns:
            row['llm_annotation'] = mode_or_unknown(frame['llm_annotation'])
            row['llm_confidence'] = mode_or_unknown(
                frame.get('llm_confidence', pd.Series('unknown', index=frame.index))
            )
            row['llm_rationale'] = mode_or_unknown(
                frame.get('llm_rationale', pd.Series('', index=frame.index))
            )
            row['llm_review_note'] = mode_or_unknown(
                frame.get('llm_review_note', pd.Series('', index=frame.index))
            )
            if 'marker_label' in frame.columns:
                row['llm_marker_comparison'] = (
                    'agree'
                    if row['llm_annotation'].lower() == mode_or_unknown(frame['marker_label']).lower()
                    else 'different'
                )
        for source, target in (
            ('marker_score', 'mean_marker_score'),
            ('annotation_confidence', 'mean_confidence'),
            ('annotation_score_margin', 'mean_score_margin'),
            ('annotation_doublet_score', 'mean_doublet_score'),
            ('annotation_lineage_mixture_score', 'mean_lineage_mixture_score'),
            ('annotation_ambient_score', 'mean_ambient_score'),
        ):
            if source in frame.columns:
                row[target] = round(float(pd.to_numeric(frame[source], errors='coerce').mean()), 4)
            else:
                row[target] = np.nan
        row['needs_review'] = bool(
            row.get('evidence_tier') == 'unresolved-review'
            or celltype == 'Unknown'
            or final_label == 'Unknown'
            or annotation_status.lower() == 'unknown'
            or (np.isfinite(label_agreement) and label_agreement < 0.8)
            or (np.isfinite(row['mean_confidence']) and row['mean_confidence'] < 0.2)
            or (np.isfinite(row['mean_score_margin']) and row['mean_score_margin'] < 0.05)
            or annotation_status.lower().startswith('review_')
        )
        if 'celltypist_conflict_fraction' in row:
            row['needs_review'] = bool(
                row['needs_review']
                or row['celltypist_conflict_fraction'] > 0
                or row['celltypist_unassigned_fraction'] > 0
                or row['celltypist_multi_label_fraction'] > 0
            )
        if row.get('llm_confidence') in {'low', 'unknown'}:
            row['needs_review'] = True
        rows.append(row)
    return pd.DataFrame(rows)


def build_llm_cluster_summaries(adata, cluster_key, marker_selection, marker_decisions):
    """Create the cluster-only evidence allowed to leave the platform.

    The helper intentionally has no access to AnnData matrices, obs metadata,
    project paths, annotation comments, or cell identifiers.  It produces the
    small marker summary consumed by ``run_llm_cluster_annotation``.
    """
    import pandas as pd

    if cluster_key not in adata.obs.columns:
        raise ValueError(f"聚类列 '{cluster_key}' 不存在，无法准备 LLM 注释。")
    marker_map = (marker_selection or {}).get('cluster_markers', {}) or {}
    decision_map = (marker_decisions or {}).get('cluster_evidence', {}) or {}
    labels = adata.obs[cluster_key].astype(str)
    summaries = []
    for cluster in sorted(labels.unique().tolist(), key=lambda value: (len(value), value)):
        evidence = decision_map.get(str(cluster), {}) or {}
        support = evidence.get('detected_positive_markers', []) or []
        summaries.append({
            'cluster': str(cluster),
            'n_cells': int((labels == str(cluster)).sum()),
            'top_markers': list(marker_map.get(str(cluster), [])),
            'marker_panel_candidate': evidence.get('final_label', 'Unknown'),
            'marker_panel_support': ', '.join(str(gene) for gene in support),
        })
    if not summaries:
        raise ValueError('没有可用于 LLM 注释的 cluster。')
    return summaries

class AnnotationAnalysis(BaseAnalysis):
    MODULE_NAME = "annotation"
    DISPLAY_NAME = "细胞注释"
    DESCRIPTION = "基于 Marker 基因的细胞类型自动注释"
    INPUT_REQUIRES = ['leiden']

    def run(self, input_path):
        import scanpy as sc
        import numpy as np
        import pandas as pd
        from modules.native_figures import (
            bar_figure,
            grouped_bar_figure,
            heatmap_figure,
            marker_dotplot_figure,
            umap_figure,
        )

        runtime_warnings = []
        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        # 输入可能已带人工/上游 celltype：先保留旧标签供审计，避免本模块
        # 静默覆盖用户标注后无从追溯。
        if 'celltype' in adata.obs.columns and 'previous_celltype' not in adata.obs.columns:
            adata.obs['previous_celltype'] = adata.obs['celltype'].astype(str)
            runtime_warnings.append(
                '输入已包含 celltype 列，已保存为 previous_celltype；'
                '本次运行将写入新的 celltype。'
            )
        requested_cluster_key = str(self.params.get('cluster_key', 'leiden') or '').strip()
        resolution = self.params.get('resolution', '0.8')
        # ``cluster_key`` is an explicit user choice and must never be silently
        # replaced by the separate resolution convenience field.  The latter
        # only selects ``leiden_<resolution>`` when the user left the cluster
        # field on its generic ``leiden`` default.  Otherwise a task can show
        # annotation evidence for a different partition from the one selected
        # in the UI (for example leiden_0.6 requested but leiden_0.8 used).
        resolution_key = f'leiden_{resolution}'
        if requested_cluster_key in {'', 'leiden'} and resolution_key in adata.obs.columns:
            preferred_cluster_key = resolution_key
        else:
            preferred_cluster_key = requested_cluster_key
        # 高分辨率聚类可能产生远多于 50 个簇；注释按簇独立判定，因此放宽
        # 分组安全上限，让所有簇都能进入初步注释而不是被拦在门外。
        leiden_key, cluster_info = resolve_obs_grouping(
            adata, preferred_cluster_key,
            fallbacks=[requested_cluster_key, 'leiden', 'leiden_0.8', 'leiden_0.6', 'leiden_1.0'],
            max_categories=200, max_numeric_categories=200,
            require_multiple=False,
        )
        if leiden_key is None:
            raise ValueError(
                f"cluster_key '{requested_cluster_key}' 不是有效的分类聚类列："
                f"{cluster_info.get('reason', '')}"
            )
        if leiden_key != preferred_cluster_key:
            self.progress(-1, f"聚类列已改用 '{leiden_key}'：{cluster_info.get('requested_reason', '')}")
        if not isinstance(adata.obs[leiden_key].dtype, pd.CategoricalDtype):
            adata.obs[leiden_key] = adata.obs[leiden_key].astype(str).astype('category')

        self.progress(20, "Scoring cell type markers...")
        requested_method = str(self.params.get('method', 'auto_marker') or 'auto_marker').strip()
        supported_methods = {'multi_evidence', 'auto_marker', 'llm_assisted', 'manual'}
        if requested_method not in supported_methods:
            raise ValueError(
                f"不支持的注释方法 '{requested_method}'；当前仅支持 "
                "multi_evidence、auto_marker、llm_assisted 和 manual。"
            )
        use_llm_annotation = requested_method == 'llm_assisted'
        # LLM labels are based on an independent, cluster-level marker summary.
        # Keep the existing marker scoring path as the traceable evidence layer.
        method = 'auto_marker' if use_llm_annotation else requested_method
        annotation_version = str(self.params.get('annotation_version', 'v1') or 'v1').strip() or 'v1'
        annotation_comment = str(self.params.get('annotation_comment', '') or '').strip()
        multi_evidence = method == 'multi_evidence'
        if multi_evidence:
            method = 'auto_marker'
        use_celltypist_reference = bool(self.params.get('use_celltypist_reference', False))
        if use_celltypist_reference and requested_method == 'manual':
            runtime_warnings.append('CellTypist 参考交叉验证仅在 auto_marker/multi_evidence 中启用；manual 模式已跳过。')
            use_celltypist_reference = False
        # Keep the PCA/neighbour representation in the persisted object, but
        # run every marker score/rank/expression plot on a biologically
        # interpretable non-negative expression scale.  In particular, signed
        # Pearson residuals must never drive cell-type labels.
        original_expression = adata.X
        annotation_expression, annotation_expression_source = _annotation_expression_matrix(adata)
        adata.X = annotation_expression
        had_temporary_source = '_annotation_temporary_expression_source' in adata.uns
        previous_temporary_source = adata.uns.get('_annotation_temporary_expression_source')
        adata.uns['_annotation_temporary_expression_source'] = annotation_expression_source
        requested_marker_set_name = normalize_marker_set_name(
            self.params.get('marker_set', 'Auto') or 'Auto'
        )
        if requested_marker_set_name not in set(MARKER_SETS).union({'Organoid', 'Auto'}):
            raise ValueError(f"不支持的 Marker 基因集 '{requested_marker_set_name}'。")
        organoid_type = str(self.params.get('organoid_type', 'intestinal') or 'intestinal').strip().lower()
        custom_markers_str = self.params.get('custom_markers', '').strip()
        if requested_marker_set_name == 'Auto' and custom_markers_str:
            # 自定义 Marker 是用户明确指定的完整面板，不能与自动选择的组织
            # 层级混用。
            marker_set_name = 'Universal'
            marker_panel_selection = {
                'requested': 'Auto', 'applied': marker_set_name,
                'mode': 'custom_markers',
                'reason': '已提供自定义 Marker，未自动选择组织面板。',
            }
        else:
            marker_set_name, marker_panel_selection = resolve_auto_marker_set(
                adata, requested_marker_set_name, organoid_type=organoid_type,
            )
        if requested_marker_set_name == 'Auto':
            self.progress(-1, marker_panel_selection['reason'])
        confidence_method = self.params.get('confidence_method', 'none')
        if confidence_method not in {'none', 'entropy', 'score_margin'}:
            raise ValueError(f"不支持的 confidence_method '{confidence_method}'。")
        marker_selection_method = str(
            self.params.get('marker_selection_method', 'wilcoxon') or 'wilcoxon'
        )
        if marker_selection_method not in {'wilcoxon', 'logreg'}:
            raise ValueError(
                f"不支持的 marker_selection_method '{marker_selection_method}'。"
            )
        if requested_method == 'manual' and not custom_markers_str:
            raise ValueError('manual 注释必须填写 ClusterID:CellType 映射。')
        # Confidence is a review signal by default.  Automatic relabelling is
        # opt-in and, in multi-evidence mode, additionally requires an
        # internally inconsistent cluster so a coherent cluster is not erased
        # merely because entropy is low across a broad universal panel.
        mark_unknown = self.params.get('mark_unknown', False)
        confidence_cutoff = float(self.params.get('confidence_cutoff', 0.2))
        confidence_cutoff = min(1.0, max(0.0, confidence_cutoff))
        # Score-profile correlation is not a biologically sufficient merging
        # rule. Retire the old automatic threshold but retain explicit,
        # reviewable cluster maps for cases with corroborating evidence.
        legacy_merge_similar_threshold = float(
            self.params.get('merge_similar_threshold', 0)
        )
        if (
            not np.isfinite(legacy_merge_similar_threshold)
            or not 0 <= legacy_merge_similar_threshold <= 1
        ):
            raise ValueError('merge_similar_threshold 必须位于 [0, 1]。')
        if legacy_merge_similar_threshold > 0:
            runtime_warnings.append(
                '已忽略旧版“相似簇合并阈值”：相关性不足以自动合并细胞类型；'
                '如已完成逐簇复核，请使用 cluster_merge_map 明确指定。'
            )
        confirmed_cluster_merges = parse_cluster_merge_mapping(
            self.params.get('cluster_merge_map', '')
        )
        available_clusters = set(adata.obs[leiden_key].astype(str).unique())
        unknown_merge_clusters = sorted({
            cluster
            for cluster_ids in confirmed_cluster_merges.values()
            for cluster in cluster_ids
            if cluster not in available_clusters
        })
        if unknown_merge_clusters:
            raise ValueError(
                f'cluster_merge_map 包含不存在于 {leiden_key} 的 cluster：'
                + ', '.join(unknown_merge_clusters)
            )
        # 初步注释细化：Universal 面板先判大谱系，再在带细分组的谱系内
        # 输出亚型标签；亚型证据不足时回退到大谱系。Colorectal/Organoid
        # 面板本身已是两级/亚型粒度，不再重复拆分。
        fine_annotation = bool(self.params.get('fine_annotation', True))
        annotate_all = bool(self.params.get('annotate_all', True))
        # 自定义 Marker 面板被视为用户显式指定的完整方案，不叠加 Universal
        # 亚型拆分；custom_markers_str 在此处已确定（manual 校验之前）。
        fine_mode = bool(
            fine_annotation
            and marker_set_name not in {'Colorectal', 'Organoid'}
            and marker_set_name == 'Universal'
            and not custom_markers_str
        )
        if fine_annotation and marker_set_name in {'TME', 'Immune', 'Blood', 'PBMC'}:
            runtime_warnings.append(
                '当前 Marker 面板已自带亚型粒度标签，细分注释开关不重复拆分；'
                '如需要 CD4+/CD8+ T、cDC/pDC 等细化，请改用 Universal 面板。'
            )
        if marker_set_name == 'Organoid' and organoid_type not in ORGANOID_MARKER_SETS:
            runtime_warnings.append(
                f"未知类器官类型 '{organoid_type}'，已回退到 intestinal marker panel。"
            )
            organoid_type = 'intestinal'

        negative_marker_text = str(self.params.get('negative_markers', '') or '').strip()
        negative_marker_weight = min(
            2.0, max(0.0, float(self.params.get('negative_marker_weight', 0.5)))
        )
        markers = {}
        negative_markers = {}
        marker_coverage = {}
        negative_marker_coverage = {}
        marker_decisions = {
            'cluster_annotations': {}, 'cluster_evidence': {},
        }
        marker_selection = None

        if method == 'manual' and custom_markers_str:
            self.progress(40, "Applying manual cell type mapping...")
            manual_mapping = {}
            for cluster_id, cell_type in parse_marker_mapping(custom_markers_str).items():
                # parse_marker_mapping 的 RHS 是按逗号拆分的 list；
                # manual 格式为 ClusterID:CellType，取第一个元素即完整类型名。
                # 若存 list 进 mapping，map() 后 Series 元素为 list，
                # .astype('category') 会抛 TypeError: unhashable type: 'list'。
                manual_mapping[cluster_id] = (cell_type[0] if cell_type else 'Unknown')

            if manual_mapping:
                adata.obs['celltype'] = adata.obs[leiden_key].astype(str).map(manual_mapping)
                adata.obs['celltype'] = adata.obs['celltype'].fillna('Unknown').astype('category')
                markers = {}
            else:
                self.progress(45, "No valid mapping found, falling back to auto_marker...")
                method = 'auto_marker'

        if method == 'auto_marker':
            if custom_markers_str:
                # Parse custom markers: "CellType1:GENE1,GENE2;CellType2:GENE3,GENE4"
                markers = parse_marker_mapping(custom_markers_str)
                if not markers:
                    markers = dict(get_marker_set(marker_set_name, organoid_type))
            else:
                markers = dict(get_marker_set(marker_set_name, organoid_type))

            if fine_mode:
                # 把第二级亚型面板并入打分（不修改模块级常量），大谱系
                # 面板仍作为第一级判定保留。
                for subtype, genes in UNIVERSAL_FINE_SUBTYPE_MARKERS.items():
                    markers.setdefault(subtype, list(genes))

            negative_markers = get_negative_marker_set(marker_set_name, organoid_type)
            if negative_marker_text:
                custom_negative = parse_marker_mapping(negative_marker_text)
                if custom_negative:
                    negative_markers.update(custom_negative)
            if fine_mode:
                # 亚型继承所属大谱系的负向 Marker，防止跨谱系污染信号
                # 把亚型分数抬高。
                negative_markers = dict(negative_markers)
                for subtype, parent in UNIVERSAL_FINE_LINEAGE_MAP.items():
                    negative_markers.setdefault(
                        subtype, list(negative_markers.get(parent, [])),
                    )

            min_markers = int(self.params.get('min_markers_per_type', 2))
            usable_markers, marker_coverage = resolve_marker_genes(markers, adata.var_names, min_markers)
            usable_negative_markers, negative_marker_coverage = resolve_marker_genes(
                negative_markers, adata.var_names, 1,
            )
            skipped_types = [ct for ct in markers if ct not in usable_markers]
            if skipped_types:
                message = f"{len(skipped_types)} 个类型的可用 marker 少于 {min_markers}，未参与自动判定。"
                runtime_warnings.append(message)
                self.progress(-1, message)
            for ct, available_genes in usable_markers.items():
                try:
                    sc.tl.score_genes(adata, available_genes, score_name=f'score_{ct}', use_raw=False)
                except RuntimeError as exc:
                    # Tiny targeted panels can lack Scanpy control genes. Keep
                    # the annotation usable with an explicit mean-expression
                    # fallback rather than silently omitting the lineage.
                    expr = adata[:, available_genes].X
                    if hasattr(expr, 'toarray'):
                        expr = expr.toarray()
                    adata.obs[f'score_{ct}'] = np.asarray(expr, dtype=float).mean(axis=1)
                    self.progress(-1, f'{ct} 使用平均 Marker 表达评分（{exc}）')
                # Keep the observed score for auditability; maturity evidence
                # is reported separately and never changes cell-type scoring.
                adata.obs[f'raw_marker_{ct}'] = pd.to_numeric(
                    adata.obs[f'score_{ct}'], errors='coerce'
                ).fillna(0.0)
                adata.obs[f'score_{ct}'] = adata.obs[f'raw_marker_{ct}']

            for ct, available_genes in usable_negative_markers.items():
                negative_name = f'negative_{ct}'
                try:
                    sc.tl.score_genes(adata, available_genes, score_name=negative_name, use_raw=False)
                    negative_values = pd.to_numeric(
                        adata.obs[negative_name], errors='coerce'
                    ).fillna(0.0).to_numpy(dtype=float)
                except RuntimeError as exc:
                    negative_values = _mean_expression_for_genes(adata, available_genes)
                    self.progress(-1, f'{ct} 负向 marker 使用平均表达评分（{exc}）')
                negative_values = np.maximum(negative_values, 0.0)
                adata.obs[f'negative_score_{ct}'] = negative_values
                if f'score_{ct}' in adata.obs.columns:
                    adata.obs[f'score_{ct}'] = (
                        pd.to_numeric(adata.obs[f'score_{ct}'], errors='coerce').fillna(0.0)
                        - negative_marker_weight * negative_values
                    )

            self.progress(50, "Assigning cell types to clusters...")
            score_cols = [f'score_{ct}' for ct in markers if f'score_{ct}' in adata.obs.columns]
            if score_cols:
                if marker_set_name in {'Colorectal_refined', 'Organoid'}:
                    # 细标签优先契约（Colorectal refined 与 Organoid 面板）
                    # 同时依赖独立排名的 cluster marker 与模块分数。在最终
                    # 判定前先计算，避免丰度高的共享上皮/代谢程序覆盖有自身
                    # 独立 signature 的簇。
                    self.progress(48, 'Ranking cluster-specific markers for refined annotation...')
                    marker_selection = select_cluster_marker_genes(
                        adata,
                        leiden_key,
                        classic_markers=usable_markers,
                        method=marker_selection_method,
                        n_rank_genes=int(self.params.get('marker_rank_genes', 200)),
                        min_markers=int(self.params.get('marker_min_per_cluster', 2)),
                        max_markers=int(self.params.get('marker_max_per_cluster', 5)),
                        padj_cutoff=float(self.params.get('marker_padj_cutoff', 0.05)),
                        min_pct=float(self.params.get('marker_min_pct', 0.10)),
                        min_delta_pct=float(self.params.get('marker_min_delta_pct', 0.05)),
                    )
                    if marker_selection.get('warnings'):
                        runtime_warnings.extend(marker_selection['warnings'])
                        self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
                marker_decisions = build_marker_decisions(
                    adata,
                    leiden_key,
                    usable_markers,
                    marker_set_name=marker_set_name,
                    multi_evidence=multi_evidence,
                    agreement_threshold=float(
                        self.params.get('cluster_agreement_threshold', 0.6)
                    ),
                    min_annotation_score=float(
                        self.params.get('min_annotation_score', 0.0)
                    ),
                    marker_min_pct=float(self.params.get('marker_min_pct', 0.10)),
                    marker_min_delta_pct=float(
                        self.params.get('marker_min_delta_pct', 0.05)
                    ),
                    broad_markers=(
                        {ct: genes for ct, genes in usable_markers.items()
                         if ct not in UNIVERSAL_FINE_SUBTYPE_MARKERS}
                        if fine_mode else None
                    ),
                    subtype_markers=(
                        {ct: genes for ct, genes in usable_markers.items()
                         if ct in UNIVERSAL_FINE_SUBTYPE_MARKERS}
                        if fine_mode else None
                    ),
                    subtype_lineage_map=(
                        UNIVERSAL_FINE_LINEAGE_MAP if fine_mode else None
                    ),
                    cluster_marker_selection=marker_selection,
                    annotate_all=annotate_all,
                )
                cluster_annotations = marker_decisions['cluster_annotations']
                cluster_values = adata.obs[leiden_key].astype(str)
                adata.obs['celltype'] = cluster_values.map(
                    cluster_annotations
                ).fillna('Unknown').astype('category')
                evidence_map = marker_decisions.get('cluster_evidence', {})
                for column, evidence_key, default in (
                    ('annotation_broad_candidate', 'broad_candidate', ''),
                    ('annotation_top_candidate', 'top_candidate', ''),
                    ('annotation_second_candidate', 'second_candidate', ''),
                    ('annotation_candidate_score', 'adjusted_candidate_score', np.nan),
                    ('annotation_second_candidate_score', 'second_candidate_score', np.nan),
                    ('annotation_cluster_score_margin', 'score_margin', np.nan),
                    ('cluster_annotation_agreement', 'cell_prediction_agreement', np.nan),
                    ('annotation_positive_marker_score', 'positive_marker_score', np.nan),
                    ('annotation_negative_marker_penalty', 'negative_marker_penalty', np.nan),
                    ('annotation_candidate_marker_detection_fraction', 'candidate_marker_detection_fraction', np.nan),
                    ('annotation_candidate_marker_specificity_delta', 'candidate_marker_specificity_delta', np.nan),
                    ('annotation_decision_reason', 'decision_reason', ''),
                    ('annotation_evidence_tier', 'evidence_tier', ''),
                ):
                    mapping = {
                        str(cluster): values.get(evidence_key, default)
                        for cluster, values in evidence_map.items()
                    }
                    adata.obs[column] = cluster_values.map(mapping).fillna(default)
                for column, evidence_key in (
                    ('annotation_detected_positive_markers', 'detected_positive_markers'),
                    ('annotation_missing_positive_markers', 'missing_positive_markers'),
                    ('annotation_refined_signature_overlap', 'refined_signature_overlap'),
                    ('annotation_refined_signature_anchor_overlap', 'refined_signature_anchor_overlap'),
                    ('annotation_umap_adjacent_clusters', 'umap_adjacent_clusters'),
                    ('annotation_knn_adjacent_clusters', 'knn_adjacent_clusters'),
                ):
                    mapping = {
                        str(cluster): ', '.join(values.get(evidence_key, []))
                        for cluster, values in evidence_map.items()
                    }
                    adata.obs[column] = cluster_values.map(mapping).fillna('')
                for column, evidence_key, default in (
                    ('annotation_negative_marker_conflict', 'negative_marker_conflict', False),
                ):
                    mapping = {
                        str(cluster): values.get(evidence_key, default)
                        for cluster, values in evidence_map.items()
                    }
                    adata.obs[column] = cluster_values.map(mapping).fillna(default)
            else:
                # 没有任何可用 marker 时绝不能把 Leiden 簇 ID 当作细胞类型；
                # 一律标记 Unknown 并在警告中说明原因。
                adata.obs['celltype'] = pd.Series(
                    'Unknown', index=adata.obs.index, dtype='category',
                )
                adata.obs['annotation_evidence_tier'] = 'unresolved-review'
                adata.obs['annotation_decision_reason'] = 'no_usable_marker_scores'
                runtime_warnings.append(
                    '没有可用的 marker 评分，全部细胞标记为 Unknown；'
                    '请检查 marker 基因集、物种或基因名（区分大小写）。'
                )
                self.progress(
                    -1, '没有可用的 marker 评分，全部细胞标记为 Unknown',
                )
            markers = usable_markers
        else:
            marker_coverage = {}
            negative_marker_coverage = {}
            negative_markers = {}

        # Quality evidence is recorded for both auto_marker and
        # multi_evidence.  Scrublet provenance is inherited from QC; marker
        # ambiguity is kept separately as lineage-mixture evidence.
        score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
        quality_markers = markers if isinstance(markers, dict) else {}
        quality_score_cols = score_cols
        quality_assigned_labels = adata.obs.get('celltype')
        quality_label_families = {}
        if marker_set_name == 'Colorectal':
            quality_markers = {
                label: markers[label]
                for label in COLORECTAL_LINEAGE_MARKERS
                if label in markers
            }
            quality_score_cols = [
                f'score_{label}' for label in quality_markers
                if f'score_{label}' in adata.obs.columns
            ]
            quality_assigned_labels = adata.obs.get(
                'annotation_broad_candidate', quality_assigned_labels,
            )
            quality_label_families = {
                label: 'Epithelial'
                for label in COLORECTAL_EPITHELIAL_SUBTYPE_MARKERS
            }
            quality_label_families['Epithelial'] = 'Epithelial'
        elif marker_set_name == 'Colorectal_refined':
            # 精细参考面板的所有标签都属于肠上皮；代谢、应激和细胞周期
            # 状态不能被当成跨谱系混合或环境 RNA 证据。
            quality_markers = markers
            quality_score_cols = [
                f'score_{label}' for label in quality_markers
                if f'score_{label}' in adata.obs.columns
            ]
            quality_label_families = {
                label: 'Epithelial' for label in quality_markers
            }
        elif fine_mode:
            # 细化模式在 broad 谱系层面评估谱系混合与环境 RNA：同谱系的
            # CD4+/CD8+ T 等亚型不算混合证据，亚型 Marker 也不被误判为
            # 异源环境 RNA。
            broad_labels = set(UNIVERSAL_FINE_LINEAGE_MAP.values())
            quality_markers = {
                label: markers[label]
                for label in markers
                if label in broad_labels or label in UNIVERSAL_FINE_SUBTYPE_MARKERS
            }
            quality_score_cols = [
                f'score_{label}' for label in broad_labels
                if f'score_{label}' in adata.obs.columns
            ]
            quality_assigned_labels = adata.obs.get(
                'annotation_broad_candidate', quality_assigned_labels,
            )
            quality_label_families = dict(UNIVERSAL_FINE_LINEAGE_MAP)
            for broad_label in broad_labels:
                quality_label_families.setdefault(broad_label, broad_label)
        quality_evidence = compute_annotation_quality_evidence(
            adata,
            quality_score_cols,
            quality_markers,
            lineage_mixture_threshold=float(self.params.get(
                'lineage_mixture_threshold',
                self.params.get('doublet_score_threshold', 0.30),
            )),
            assigned_labels=quality_assigned_labels,
            label_families=quality_label_families,
            ambient_threshold=float(self.params.get('ambient_score_threshold', 0.35)),
            ambient_prevalence=float(self.params.get('ambient_prevalence', 0.50)),
        )
        adata.obs['annotation_doublet_score'] = quality_evidence['doublet_score']
        adata.obs['annotation_doublet_status'] = quality_evidence['doublet_status']
        adata.obs['annotation_lineage_mixture_score'] = quality_evidence['lineage_mixture_score']
        adata.obs['annotation_lineage_mixture_status'] = quality_evidence['lineage_mixture_status']
        adata.obs['annotation_ambient_score'] = quality_evidence['ambient_score']
        adata.obs['annotation_ambient_status'] = quality_evidence['ambient_status']
        adata.obs['annotation_top1'] = quality_evidence['top1']
        adata.obs['annotation_top2'] = quality_evidence['top2']
        adata.obs['annotation_top1_score'] = quality_evidence['top1_score']
        adata.obs['annotation_top2_score'] = quality_evidence['top2_score']
        if quality_evidence.get('warnings'):
            runtime_warnings.extend(quality_evidence['warnings'])

        maturity_evidence = None
        if marker_set_name in {'Organoid', 'Colorectal', 'Colorectal_refined'}:
            maturity_evidence = compute_organoid_maturity_evidence(
                adata,
                organoid_type=(
                    'intestinal' if marker_set_name in {'Colorectal', 'Colorectal_refined'} else organoid_type
                ),
                time_key=self.params.get('maturity_time_key', ''),
            )
            for column, values in (
                ('organoid_progenitor_score', maturity_evidence['progenitor_score']),
                ('organoid_mature_score', maturity_evidence['mature_score']),
                ('organoid_cycling_score', maturity_evidence['cycling_score']),
                ('organoid_maturity_index', maturity_evidence['maturity_index']),
                ('organoid_maturity_state', maturity_evidence['maturity_state']),
                ('organoid_time_label', maturity_evidence['time_labels']),
                ('organoid_time_numeric', maturity_evidence['time_numeric']),
            ):
                adata.obs[column] = values
            adata.uns['organoid_maturity_evidence'] = {
                'organoid_type': organoid_type,
                'time_key': maturity_evidence['time_key'],
                'expression_source': maturity_evidence['expression_source'],
                'time_trend_spearman': (
                    float(maturity_evidence['time_trend_spearman'])
                    if np.isfinite(maturity_evidence['time_trend_spearman']) else None
                ),
                'marker_coverage': maturity_evidence['marker_coverage'],
                'warnings': maturity_evidence['warnings'],
            }
            runtime_warnings.extend(maturity_evidence['warnings'])

        state_evidence = compute_cell_state_evidence(
            adata, threshold=float(self.params.get('state_score_threshold', 0.35))
        )
        for state_name, values in state_evidence['scores'].items():
            column_name = 'state_score_' + state_name.lower().replace(' ', '_').replace('-', '_')
            adata.obs[column_name] = values
            high_column = 'state_high_' + state_name.lower().replace(' ', '_').replace('-', '_')
            adata.obs[high_column] = state_evidence['high_flags'][state_name]
        adata.obs['dominant_cell_state'] = state_evidence['dominant_state']
        # Backward-compatible alias; the name is retained for downstream
        # readers, while metadata now states that this is the dominant view.
        adata.obs['cell_state'] = state_evidence['dominant_state']
        adata.obs['cell_state_flags'] = state_evidence['state_flags']
        adata.uns['cell_state_evidence'] = {
            'coverage': state_evidence['coverage'],
            'threshold': state_evidence['threshold'],
            'expression_source': state_evidence['expression_source'],
            'dominant_state_column': 'dominant_cell_state',
            'multi_label_column': 'cell_state_flags',
            'multi_label_boolean_prefix': 'state_high_',
        }

        # Multi-evidence mode keeps the rule engine interpretable and records
        # independent per-cell marker predictions for cluster-level review.
        if multi_evidence:
            cluster_values = adata.obs[leiden_key].astype(str)
            cluster_evidence = marker_decisions.get('cluster_evidence', {})
            top_candidate_map = {
                str(cluster): values.get('top_candidate', 'Unknown')
                for cluster, values in cluster_evidence.items()
            }
            adata.obs['marker_label'] = cluster_values.map(
                top_candidate_map
            ).fillna('Unknown')
            for key in (
                'prediction_1', 'prediction_2',
                'prediction_1_score', 'prediction_2_score',
            ):
                values = marker_decisions.get(key)
                if values is not None and len(values) == adata.n_obs:
                    adata.obs[key] = values
            if 'prediction_1_score' in adata.obs.columns:
                adata.obs['marker_score'] = pd.to_numeric(
                    adata.obs['prediction_1_score'], errors='coerce'
                )
            adata.uns['annotation_agreement_source'] = (
                'hierarchical_prediction_1'
                if marker_set_name == 'Colorectal' or fine_mode else 'prediction_1'
            )
            adata.obs['final_annotation'] = adata.obs['celltype'].astype(str)
            adata.obs['annotation_status'] = np.where(
                adata.obs['final_annotation'].eq('Unknown'), 'Unknown', 'review'
            )

        # In LLM mode, rank cluster markers before plotting and send only this
        # compact evidence to the configured provider.  No raw expression,
        # cell identifiers, donor/sample metadata, project paths, or free-form
        # annotation comments can reach the request builder.
        llm_evidence = None
        if use_llm_annotation:
            self.progress(56, 'Preparing de-identified cluster marker summaries for LLM annotation...')
            marker_selection = select_cluster_marker_genes(
                adata,
                leiden_key,
                classic_markers=markers if isinstance(markers, dict) else {},
                method=marker_selection_method,
                n_rank_genes=int(self.params.get('marker_rank_genes', 200)),
                min_markers=int(self.params.get('marker_min_per_cluster', 2)),
                max_markers=int(self.params.get('marker_max_per_cluster', 5)),
                padj_cutoff=float(self.params.get('marker_padj_cutoff', 0.05)),
                min_pct=float(self.params.get('marker_min_pct', 0.10)),
                min_delta_pct=float(self.params.get('marker_min_delta_pct', 0.05)),
            )
            if marker_selection.get('warnings'):
                runtime_warnings.extend(marker_selection['warnings'])
                self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
            cluster_summaries = build_llm_cluster_summaries(
                adata, leiden_key, marker_selection, marker_decisions,
            )
            self.progress(58, 'Running LLM cluster annotation...')
            try:
                llm_evidence = run_llm_cluster_annotation(
                    cluster_summaries,
                    species=self.params.get('llm_species', 'human'),
                    tissue_context=self.params.get('llm_tissue_context', ''),
                    max_clusters=self.params.get('llm_max_clusters', 30),
                )
            except LLMAnnotationError as exc:
                # This is the selected primary annotation method.  Do not
                # silently relabel the output as an LLM result when the model
                # service is unavailable or returned incomplete JSON.
                raise RuntimeError(f'LLM 注释失败：{exc}') from exc

            marker_labels = adata.obs['celltype'].astype(str).copy()
            llm_by_cluster = {
                str(item['cluster']): item for item in llm_evidence['annotations']
            }
            cluster_values = adata.obs[leiden_key].astype(str)
            adata.obs['marker_label'] = marker_labels
            adata.obs['llm_annotation'] = cluster_values.map(
                {cluster: item['cell_type'] for cluster, item in llm_by_cluster.items()}
            ).fillna('Unknown')
            adata.obs['llm_confidence'] = cluster_values.map(
                {cluster: item['confidence'] for cluster, item in llm_by_cluster.items()}
            ).fillna('unknown')
            adata.obs['llm_rationale'] = cluster_values.map(
                {cluster: item['rationale'] for cluster, item in llm_by_cluster.items()}
            ).fillna('')
            adata.obs['llm_review_note'] = cluster_values.map(
                {cluster: item['review_note'] for cluster, item in llm_by_cluster.items()}
            ).fillna('')
            adata.obs['llm_marker_comparison'] = np.where(
                adata.obs['llm_annotation'].astype(str).str.lower()
                == marker_labels.astype(str).str.lower(),
                'agree', 'different',
            )
            # LLM 对某些簇判为 Unknown 时，回退到 Marker 规则候选标签，
            # 保证初步注释尽量不留 Unknown；llm_annotation 仍保留模型的
            # 原始判断，llm_unknown_fallback 标记被回退的细胞供复核。
            llm_values = adata.obs['llm_annotation'].astype(str)
            llm_unknown_mask = llm_values.str.lower().eq('unknown')
            marker_fallback_mask = llm_unknown_mask & (marker_labels != 'Unknown')
            adata.obs['llm_unknown_fallback'] = marker_fallback_mask
            fallback_labels = np.where(
                marker_fallback_mask, marker_labels, llm_values,
            )
            adata.obs['celltype'] = pd.Series(
                fallback_labels, index=adata.obs.index,
            ).astype('category')
            n_llm_fallback = int(marker_fallback_mask.sum())
            if n_llm_fallback:
                runtime_warnings.append(
                    f'LLM 将 {n_llm_fallback} 个细胞判为 Unknown，'
                    '已回退到 Marker 规则候选标签并标记复核（llm_unknown_fallback）。'
                )
                self.progress(-1, f'{n_llm_fallback} 个 LLM-Unknown 细胞回退到 Marker 候选标签')
            adata.uns['llm_annotation'] = {
                key: value for key, value in llm_evidence.items()
                if key not in {'annotations', 'cluster_input'}
            }
            # AnnData cannot persist ``list[dict]`` values in ``uns``.  Keep
            # the complete cluster-level audit losslessly as JSON instead of
            # letting h5py coerce (or reject) the nested Python objects.
            adata.uns['llm_annotation']['annotations_json'] = json.dumps(
                llm_evidence['annotations'], ensure_ascii=False, default=str,
            )
            adata.uns['llm_annotation']['annotation_count'] = int(
                len(llm_evidence['annotations'])
            )
            if llm_evidence.get('cluster_input') is not None:
                adata.uns['llm_annotation']['cluster_input_json'] = json.dumps(
                    llm_evidence['cluster_input'], ensure_ascii=False, default=str,
                )

        # CellTypist is an optional external reference, never the final label
        # source.  It runs after Marker evidence so conflicts can be displayed
        # without changing the existing Unknown/review logic.
        celltypist_evidence = None
        if use_celltypist_reference:
            self.progress(58, 'Running local CellTypist reference...')
            try:
                celltypist_evidence = run_celltypist_reference(
                    adata,
                    model_name=self.params.get('celltypist_model', 'Immune_All_Low.pkl'),
                    mode=self.params.get('celltypist_mode', 'prob match'),
                    p_thres=float(self.params.get('celltypist_p_threshold', 0.5)),
                    majority_voting=bool(self.params.get('celltypist_majority_voting', False)),
                    over_clustering=adata.obs[leiden_key].astype(str).to_numpy(),
                )
                reference_labels = pd.Series(
                    celltypist_evidence['labels'].astype(str), index=adata.obs.index,
                )
                marker_reference = (
                    adata.obs['marker_label'].astype(str)
                    if 'marker_label' in adata.obs.columns
                    else adata.obs['celltype'].astype(str)
                )
                comparison = [
                    _celltypist_comparison_status(marker, reference)
                    for marker, reference in zip(marker_reference, reference_labels)
                ]
                adata.obs['celltypist_label'] = reference_labels.astype(str)
                adata.obs['celltypist_top2'] = pd.Series(
                    celltypist_evidence['top2'].astype(str), index=adata.obs.index,
                )
                adata.obs['celltypist_confidence'] = celltypist_evidence['confidence']
                adata.obs['celltypist_status'] = pd.Series(
                    celltypist_evidence['status'].astype(str), index=adata.obs.index,
                )
                adata.obs['celltypist_comparison'] = pd.Series(comparison, index=adata.obs.index)
                adata.obs['celltypist_agreement'] = pd.Series(
                    np.isin(comparison, ['agree', 'lineage_agree']), index=adata.obs.index,
                )
                reference_counts = {
                    str(key): int(value)
                    for key, value in reference_labels.value_counts().items()
                }
                status_counts = {
                    str(key): int(value)
                    for key, value in pd.Series(celltypist_evidence['status']).value_counts().items()
                }
                comparison_counts = {
                    str(key): int(value)
                    for key, value in pd.Series(comparison).value_counts().items()
                }
                celltypist_evidence['label_counts'] = reference_counts
                celltypist_evidence['status_counts'] = status_counts
                celltypist_evidence['comparison_counts'] = comparison_counts
                adata.uns['celltypist_reference'] = {
                    key: value for key, value in celltypist_evidence.items()
                    if key not in {'labels', 'top2', 'confidence', 'status', 'probability_max'}
                }
                if comparison_counts.get('conflict', 0):
                    runtime_warnings.append(
                        f"CellTypist 与 Marker 存在 {comparison_counts['conflict']} 个细胞的标签冲突，已保留人工复核。"
                    )
            except ModuleNotFoundError:
                message = '未安装 celltypist 依赖，已跳过 CellTypist 参考；Marker 注释仍正常完成。'
                runtime_warnings.append(message)
                self.progress(-1, message)
            except Exception as exc:
                message = f'CellTypist 参考执行失败，已回退到 Marker 注释：{exc}'
                runtime_warnings.append(message)
                self.progress(-1, message)

        # 置信度评估
        confidence_col = None
        if confidence_method != 'none' and 'celltype' in adata.obs.columns:
            self.progress(60, f"Computing annotation confidence ({confidence_method})...")
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if confidence_method == 'entropy' and score_cols:
                import numpy as np
                score_matrix = adata.obs[score_cols].values
                exp_scores = np.exp(score_matrix - score_matrix.max(axis=1, keepdims=True))
                probs = exp_scores / (exp_scores.sum(axis=1, keepdims=True) + 1e-10)
                entropy = -np.sum(probs * np.log(probs + 1e-10), axis=1)
                max_entropy = np.log(len(score_cols)) if len(score_cols) > 1 else 1
                adata.obs['annotation_confidence'] = 1 - entropy / (max_entropy + 1e-10)
                confidence_col = 'annotation_confidence'
            elif confidence_method == 'score_margin' and score_cols:
                import numpy as np
                score_matrix = adata.obs[score_cols].values
                sorted_scores = np.sort(score_matrix, axis=1)
                if sorted_scores.shape[1] >= 2:
                    adata.obs['annotation_score_margin'] = sorted_scores[:, -1] - sorted_scores[:, -2]
                else:
                    adata.obs['annotation_score_margin'] = sorted_scores[:, -1]
                confidence_col = 'annotation_score_margin'

            if mark_unknown and confidence_col in adata.obs.columns:
                raw_low_conf_mask = adata.obs[confidence_col] < confidence_cutoff
                low_conf_mask = raw_low_conf_mask.copy()
                if multi_evidence and 'cluster_annotation_agreement' in adata.obs.columns:
                    agreement_cutoff = float(self.params.get('cluster_agreement_threshold', 0.6))
                    coherent_cluster = adata.obs['cluster_annotation_agreement'] >= agreement_cutoff
                    low_conf_mask &= ~coherent_cluster
                    retained_for_review = int((raw_low_conf_mask & coherent_cluster).sum())
                    if retained_for_review:
                        runtime_warnings.append(
                            f'{retained_for_review} 个低置信度但 cluster 标签一致的细胞保留为 review，未自动改为 Unknown。'
                        )
                if hasattr(adata.obs['celltype'], 'cat') and 'Unknown' not in adata.obs['celltype'].cat.categories:
                    adata.obs['celltype'] = adata.obs['celltype'].cat.add_categories(['Unknown'])
                adata.obs.loc[low_conf_mask, 'celltype'] = 'Unknown'
                n_unknown = low_conf_mask.sum()
                if n_unknown > 0:
                    runtime_warnings.append(
                        f'有 {int(n_unknown)} 个细胞因置信度低于 {confidence_cutoff:g} 且 cluster 标签不一致而被标记为 Unknown。'
                    )
                    self.progress(-1, f"标记 {n_unknown} 个低置信度细胞为 Unknown")

        # Only a human-reviewed map may combine final annotation labels. The
        # input cluster column remains unmodified and is also written to
        # ``annotation_source_cluster`` below, preserving all evidence for a
        # later split/review.
        if confirmed_cluster_merges:
            source_clusters = adata.obs[leiden_key].astype(str)
            labels = adata.obs['celltype'].astype(str).copy()
            merge_group = pd.Series('', index=adata.obs.index, dtype=object)
            for merged_label, cluster_ids in confirmed_cluster_merges.items():
                mask = source_clusters.isin(cluster_ids)
                labels.loc[mask] = merged_label
                merge_group.loc[mask] = merged_label
            adata.obs['celltype'] = labels.astype('category')
            adata.obs['annotation_merge_group'] = merge_group.astype('category')
            self.progress(
                -1,
                f'按已确认映射合并 {sum(len(ids) for ids in confirmed_cluster_merges.values())} 个来源 cluster；原 cluster 已保留。',
            )

        # The reviewed display names apply to the current colorectal refined
        # run only.  They make every downstream cell-type plot use the same
        # terminology while retaining the pre-display marker programme in a
        # separate, auditable obs column.
        curated_display_labels_applied = apply_curated_annotation_display_labels(
            adata, marker_set_name,
        )
        if curated_display_labels_applied:
            self.progress(
                -1,
                f'已应用 {curated_display_labels_applied} 个细胞的复核显示标签；'
                '原始 marker programme 已保存在 annotation_marker_programme。',
            )

        # 审阅 CSV 应用：annotation_manual_map_template.csv 只改人工确认的
        # cluster 的显示标签，来源 cluster 列与逐簇证据从不删除。应用后
        # 覆盖 evidence tier 为 confirmed，并保留原始 decision_reason 供审计。
        manual_map_csv_path = str(self.params.get('manual_map_csv', '') or '').strip()
        manual_map_applied = {}
        if manual_map_csv_path and os.path.isfile(manual_map_csv_path):
            manual_map_applied = parse_annotation_manual_map_csv(manual_map_csv_path)
            available_clusters = set(adata.obs[leiden_key].astype(str).unique())
            unknown_clusters = sorted(set(manual_map_applied) - available_clusters)
            if unknown_clusters:
                raise ValueError(
                    'annotation_manual_map_template.csv 包含不存在的 cluster：'
                    + ', '.join(unknown_clusters[:8]),
                )
            if manual_map_applied:
                source_clusters = adata.obs[leiden_key].astype(str)
                labels = adata.obs['celltype'].astype(str).copy()
                previous_reasons = adata.obs['annotation_decision_reason'].astype(str).copy()
                for cluster_id, final_label in manual_map_applied.items():
                    mask = source_clusters == str(cluster_id)
                    labels.loc[mask] = final_label
                    adata.obs.loc[mask, 'annotation_decision_reason'] = (
                        'manual_map_applied:' + previous_reasons[mask].astype(str)
                    )
                    adata.obs.loc[mask, 'annotation_evidence_tier'] = 'confirmed'
                adata.obs['celltype'] = labels.astype('category')
                adata.obs['annotation_manual_map_applied'] = source_clusters.isin(
                    list(manual_map_applied)
                )
                prior_display_source = adata.obs.get(
                    'annotation_display_label_source',
                    pd.Series('marker_programme', index=adata.obs.index),
                )
                adata.obs['annotation_display_label_source'] = np.where(
                    adata.obs['annotation_manual_map_applied'],
                    'manual_cluster_map', prior_display_source,
                )
                self.progress(
                    -1,
                    f'已应用人工审阅 CSV：{len(manual_map_applied)} 个来源 cluster 的显示标签已更新；来源 cluster 已保留。',
                )

        self.progress(70, "Generating dotplot and UMAP...")
        plots_dir = self.ensure_plots_dir()
        result_files = []
        from modules.sc_figure_diagnostics import (
            composition_by_group,
            composition_by_group_figure,
            composition_heatmap_figure,
            condition_composition_figure,
            resolve_donor_key,
            state_score_by_group_figure,
            state_score_heatmap_figure,
        )
        donor_key, _ = resolve_donor_key(
            adata,
            self.params.get('donor_key', self.params.get('sample_key', '')),
            require_multiple=True,
        )
        annotation_donor_composition = {}
        annotation_condition_key = str(self.params.get('condition_key', '') or '').strip()

        # Marker dotplots are evidence displays, not fixed marker-panel
        # illustrations.  Rank genes independently for every cluster, then
        # keep a small number of data-driven markers plus at most one classic
        # anchor per cluster.  The returned metadata is written to summary so
        # the exact selection can be reviewed or reproduced.
        if marker_selection is None:
            marker_selection = select_cluster_marker_genes(
                adata,
                leiden_key,
                classic_markers=markers if isinstance(markers, dict) else {},
                method=marker_selection_method,
                n_rank_genes=int(self.params.get('marker_rank_genes', 200)),
                min_markers=int(self.params.get('marker_min_per_cluster', 2)),
                max_markers=int(self.params.get('marker_max_per_cluster', 5)),
                padj_cutoff=float(self.params.get('marker_padj_cutoff', 0.05)),
                min_pct=float(self.params.get('marker_min_pct', 0.10)),
                min_delta_pct=float(self.params.get('marker_min_delta_pct', 0.05)),
            )
            if marker_selection.get('warnings'):
                runtime_warnings.extend(marker_selection.get('warnings', []))
                self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
        dotplot_genes = marker_selection.get('genes', [])[:40]
        adata.uns['marker_selection'] = {
            'cluster_key': leiden_key,
            'cluster_markers': marker_selection.get('cluster_markers', {}),
            'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
            'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
            'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
            'warnings': marker_selection.get('warnings', []),
            'parameters': marker_selection.get('parameters', {}),
        }

        if marker_coverage:
            coverage_programmes = list(marker_coverage)
            coverage_types = [
                curated_annotation_display_label(programme, marker_set_name)
                for programme in coverage_programmes
            ]
            matched = [marker_coverage[ct]['matched'] for ct in coverage_programmes]
            totals = [marker_coverage[ct]['total'] for ct in coverage_programmes]
            fig_coverage = grouped_bar_figure(
                coverage_types,
                [('可用 marker', matched), ('marker 总数', totals)],
                title='Marker coverage by displayed cell type',
                x_label='Displayed cell type', y_label='Gene count',
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_coverage, plots_dir, 'annotation_marker_coverage.png',
                'bar', 'Marker coverage by displayed cell type', formats=('png', 'svg'), dpi=300,
            ))

        if self.params.get('show_celltype_composition', True) and 'celltype' in adata.obs.columns:
            ct_counts_plot = adata.obs['celltype'].astype(str).value_counts()
            ct_pct = ct_counts_plot / max(int(ct_counts_plot.sum()), 1) * 100
            fig_comp = bar_figure(
                ct_counts_plot.index.tolist(), ct_counts_plot.values,
                title='Cell Type Composition', x_label='Cell type',
                y_label='Cell count',
                annotations=[f'{int(n):,} ({p:.1f}%)'
                             for n, p in zip(ct_counts_plot.values, ct_pct.values)],
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_comp, plots_dir, 'annotation_celltype_composition.png',
                'bar', '细胞类型组成', formats=('png', 'svg'), dpi=300,
            ))

            if donor_key:
                count_table, fraction_table = composition_by_group(
                    adata.obs, donor_key, 'celltype', max_labels=30,
                )
                if fraction_table is not None and not fraction_table.empty:
                    annotation_donor_composition = {
                        'group_key': donor_key,
                        'counts': {
                            str(group): {str(label): int(value) for label, value in row.items()}
                            for group, row in count_table.to_dict(orient='index').items()
                        },
                        'fractions': {
                            str(group): {str(label): float(value) for label, value in row.items()}
                            for group, row in fraction_table.to_dict(orient='index').items()
                        },
                    }
                    fig_donor = composition_by_group_figure(
                        fraction_table,
                        group_label=donor_key,
                        label_label='Cell type',
                        title='Cell-type composition by donor / sample',
                    )
                    if fig_donor is not None:
                        result_files.extend(self.save_matplotlib_figure(
                            fig_donor, plots_dir, 'annotation_donor_celltype_composition.png',
                            'bar', 'Donor × cell type composition',
                            formats=('png', 'svg'), dpi=300,
                        ))
                    fig_donor_heatmap = composition_heatmap_figure(
                        fraction_table,
                        group_label=donor_key,
                        label_label='Cell type',
                        title='Cell-type fraction by donor / sample',
                    )
                    if fig_donor_heatmap is not None:
                        result_files.extend(self.save_matplotlib_figure(
                            fig_donor_heatmap, plots_dir, 'annotation_donor_celltype_heatmap.png',
                            'heatmap', 'Donor × cell type fraction heatmap',
                            formats=('png', 'svg'), dpi=300,
                        ))

            # The condition-level composition is useful even when a donor
            # column is unavailable: it is a descriptive, within-condition
            # 100% stacked plot and does not treat cells as replicates.
            if annotation_condition_key and annotation_condition_key in adata.obs.columns:
                fig_condition = condition_composition_figure(
                    adata.obs, donor_key, annotation_condition_key, 'celltype',
                )
                if fig_condition is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        fig_condition, plots_dir, 'annotation_condition_celltype_composition.png',
                        'bar', 'Condition-level cell type composition (100% stacked)',
                        formats=('png', 'svg'), dpi=300,
                    ))

        if self.params.get('show_marker_score_heatmap', True):
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if score_cols and leiden_key in adata.obs.columns:
                cluster_labels = adata.obs[leiden_key].astype(str)
                cluster_order = sorted(cluster_labels.unique(), key=lambda x: (len(x), x))
                mean_scores = []
                for cluster in cluster_order:
                    mask = cluster_labels == cluster
                    mean_scores.append(adata.obs.loc[mask, score_cols].mean().values)
                score_matrix = np.asarray(mean_scores, dtype=float).T
                row_mean = score_matrix.mean(axis=1, keepdims=True)
                row_std = score_matrix.std(axis=1, keepdims=True) + 1e-10
                z = np.clip((score_matrix - row_mean) / row_std, -3, 3)
                fig_score_heat = heatmap_figure(
                    z,
                    x_labels=cluster_order,
                    y_labels=[curated_annotation_display_label(
                        c.replace('score_', ''), marker_set_name,
                    ) for c in score_cols],
                    title=f'Candidate marker-programme scores by {leiden_key}',
                    x_label='Source cluster', y_label='Candidate marker programme',
                    colorbar_label='Row z-score',
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_score_heat, plots_dir, 'annotation_marker_score_heatmap.png',
                    'heatmap', 'Candidate marker-programme score heatmap', formats=('png', 'svg'), dpi=300,
                ))

        # Preserve the raw candidate score view and add a decision-oriented
        # heatmap with final label and top-vs-second margin under each cluster.
        decision_score_cols = [column for column in adata.obs.columns if column.startswith('score_')]
        if decision_score_cols and leiden_key in adata.obs.columns:
            decision_clusters = sorted(adata.obs[leiden_key].astype(str).unique(), key=lambda value: (len(value), value))
            decision_labels = adata.obs[leiden_key].astype(str)
            decision_matrix = pd.DataFrame(index=decision_score_cols, columns=decision_clusters, dtype=float)
            for cluster in decision_clusters:
                mask = decision_labels == cluster
                decision_matrix[cluster] = adata.obs.loc[mask, decision_score_cols].mean().to_numpy(dtype=float)
            final_labels_by_cluster = {}
            if 'celltype' in adata.obs.columns:
                for cluster in decision_clusters:
                    values = adata.obs.loc[
                        decision_labels == cluster, 'celltype'
                    ].dropna().astype(str)
                    if not values.empty:
                        final_labels_by_cluster[str(cluster)] = str(
                            values.value_counts().index[0]
                        )
            decision_display_label_map = {
                str(column).replace('score_', ''): curated_annotation_display_label(
                    str(column).replace('score_', ''), marker_set_name,
                )
                for column in decision_score_cols
            }
            from modules.sc_figure_diagnostics import annotation_decision_heatmap
            fig_decision = annotation_decision_heatmap(
                decision_matrix, decision_clusters,
                decisions=marker_decisions.get('cluster_evidence', {}),
                final_labels=final_labels_by_cluster,
                display_label_map=decision_display_label_map,
            )
            if fig_decision is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_decision, plots_dir, 'annotation_cluster_decision_heatmap.png',
                    'heatmap', 'Source cluster × final cell-type decision scores',
                    formats=('png', 'svg'), dpi=300,
                ))

        # A fixed lineage panel is intentionally grouped by source cluster, not
        # final annotation, so the validation does not become circular.
        if marker_set_name == 'Colorectal' and leiden_key in adata.obs.columns:
            canonical_genes = list(dict.fromkeys(
                gene for genes in DEFAULT_COLORECTAL_MARKERS.values() for gene in genes
            ))
            try:
                canonical_validation = build_marker_validation_matrix(
                    adata, canonical_genes, groupby=leiden_key, max_cells=5000,
                )
                if canonical_validation:
                    fig_canonical = marker_dotplot_figure(
                        canonical_validation['mean_expression'],
                        canonical_validation['detection_fraction'],
                        canonical_validation['categories'], canonical_validation['genes'],
                        title='Canonical colorectal marker panel by cluster',
                        x_label='Source cluster', y_label='Canonical marker gene',
                        value_label='Mean log1p expression',
                        detection_label='% cells expressing',
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_canonical, plots_dir, 'annotation_canonical_marker_dotplot.png',
                        'dotplot', 'Canonical colorectal marker panel',
                        formats=('png', 'svg'), dpi=300,
                    ))
            except Exception as exc:
                self.progress(-1, f'结直肠 canonical marker 图生成失败：{exc}')

            fixed_lineage_genes = list(dict.fromkeys(
                gene for lineage, genes in COLORECTAL_LINEAGE_MARKERS.items()
                if lineage != 'Epithelial' for gene in genes
            ))
            try:
                lineage_validation = build_marker_validation_matrix(
                    adata, fixed_lineage_genes, groupby=leiden_key, max_cells=5000,
                )
                if lineage_validation:
                    fig_lineage = marker_dotplot_figure(
                        lineage_validation['mean_expression'],
                        lineage_validation['detection_fraction'],
                        lineage_validation['categories'], lineage_validation['genes'],
                        title='Non-epithelial lineage validation by cluster',
                        x_label='Source cluster', y_label='Lineage marker',
                        value_label='Mean log1p expression',
                        detection_label='% cells expressing',
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_lineage, plots_dir, 'annotation_non_epithelial_validation_dotplot.png',
                        'dotplot', 'Non-epithelial lineage validation',
                        formats=('png', 'svg'), dpi=300,
                    ))
            except Exception as exc:
                self.progress(-1, f'非上皮 marker 验证图生成失败：{exc}')

            epithelial_labels = adata.obs['celltype'].astype(str).where(
                adata.obs['celltype'].astype(str).isin(COLORECTAL_EPITHELIAL_SUBTYPE_MARKERS),
                'Non-epithelial',
            )
            adata.obs['epithelial_subtype'] = epithelial_labels.astype('category')
            try:
                fig_epithelial = self.build_publication_umap(
                    adata, 'epithelial_subtype', title='Colorectal epithelial subtype UMAP',
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_epithelial, plots_dir, 'annotation_epithelial_subtype_umap.png',
                    'umap', 'Epithelial subtype UMAP', formats=('png', 'svg'), dpi=300,
                ))
            except Exception as exc:
                self.progress(-1, f'上皮亚型 UMAP 导出失败：{exc}')

        state_score_columns = [
            column for column in adata.obs.columns if column.startswith('state_score_')
        ]
        if state_score_columns:
            state_group_key = donor_key or ('epithelial_subtype' if 'epithelial_subtype' in adata.obs.columns else 'celltype')
            fig_state_heatmap = state_score_heatmap_figure(
                adata.obs, state_score_columns, state_group_key,
            )
            if fig_state_heatmap is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_state_heatmap, plots_dir, 'annotation_cell_state_score_heatmap.png',
                    'heatmap', 'Cell-state score heatmap', formats=('png', 'svg'), dpi=300,
                ))
            fig_state_group = state_score_by_group_figure(
                adata.obs, state_score_columns, state_group_key,
            )
            if fig_state_group is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_state_group, plots_dir, 'annotation_cell_state_by_group.png',
                    'violin', 'Cell-state scores by donor / group', formats=('png', 'svg'), dpi=300,
                ))
            if 'epithelial_subtype' in adata.obs.columns and state_group_key != 'epithelial_subtype':
                fig_state_subtype = state_score_heatmap_figure(
                    adata.obs, state_score_columns, 'epithelial_subtype',
                    title='Cell-state score by epithelial subtype',
                )
                if fig_state_subtype is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        fig_state_subtype, plots_dir, 'annotation_cell_state_by_epithelial_subtype.png',
                        'heatmap', 'Cell-state score × epithelial subtype',
                        formats=('png', 'svg'), dpi=300,
                    ))
            if 'X_umap' in adata.obsm:
                for column in state_score_columns:
                    try:
                        fig_state_umap = self.build_publication_umap(
                            adata, column, title=f'{column.replace("state_score_", "")} on UMAP',
                        )
                        result_files.extend(self.save_matplotlib_figure(
                            fig_state_umap, plots_dir,
                            f'annotation_{column}_umap.png', 'umap',
                            f'{column} UMAP', formats=('png', 'svg'), dpi=300,
                        ))
                    except Exception as exc:
                        self.progress(-1, f'{column} UMAP 导出失败：{exc}')

        # Cluster-specific marker dotplot for marker validation.  Annotation
        # labels remain available in the companion expression plot; grouping
        # this evidence plot by the source cluster prevents a cell-type label
        # from hiding mixed or uncertain clusters.
        if dotplot_genes and leiden_key in adata.obs.columns:
            try:
                sc.tl.dendrogram(adata, groupby=leiden_key)
                fig_dotplot = sc.pl.dotplot(
                    adata, var_names=dotplot_genes, groupby=leiden_key,
                    return_fig=True,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_dotplot, plots_dir, 'annotation_celltype_dotplot.png', 'dotplot',
                    'Cluster-specific Marker Dotplot'
                ))
                import matplotlib.pyplot as plt
                plt.close('all')
            except Exception as e:
                self.progress(-1, f"Annotation dotplot generation failed: {e}")

        # Marker expression validation.  This is intentionally a dot plot:
        # colour represents mean display expression and dot size represents the
        # fraction of cells expressing the marker.  The old grouped boxplot
        # used the already-scaled ``adata.X`` and produced a duplicated legend
        # (one entry per box), which made the evidence unreadable.
        marker_validation = None
        if self.params.get('show_marker_expression_violin', True) and dotplot_genes and 'celltype' in adata.obs.columns:
            try:
                marker_validation = build_marker_validation_matrix(
                    adata, dotplot_genes[:16], groupby='celltype', max_cells=5000,
                )
                if marker_validation:
                    fig_marker = marker_dotplot_figure(
                        marker_validation['mean_expression'],
                        marker_validation['detection_fraction'],
                        marker_validation['categories'],
                        marker_validation['genes'],
                        title='Marker Expression Validation',
                        x_label='Annotated cell type', y_label='Marker gene',
                        value_label='Mean log1p expression',
                        detection_label='% cells expressing',
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_marker, plots_dir, 'annotation_marker_expression_dotplot.png',
                        'dotplot', 'Marker 表达验证图', formats=('png', 'svg'), dpi=300,
                    ))
            except Exception as e:
                self.progress(-1, f"Marker expression plot generation failed: {e}")

        try:
            fig_static = self.build_publication_umap(
                adata, 'celltype', title='UMAP by Cell Type'
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_static, plots_dir, 'annotation_umap_celltype.png', 'umap',
                'UMAP by Cell Type'
            ))
            import matplotlib.pyplot as plt
            plt.close(fig_static)
        except Exception as exc:
            self.progress(-1, f'注释 UMAP 导出失败：{exc}')

        if celltypist_evidence is not None and 'celltypist_label' in adata.obs.columns and 'X_umap' in adata.obsm:
            try:
                fig_celltypist = self.build_publication_umap(
                    adata, 'celltypist_label', title='CellTypist Reference Labels'
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_celltypist, plots_dir, 'annotation_celltypist_reference_umap.png',
                    'umap', 'CellTypist Reference Labels UMAP', formats=('png', 'svg'),
                ))
                import matplotlib.pyplot as plt
                plt.close(fig_celltypist)
            except Exception as exc:
                self.progress(-1, f'CellTypist UMAP 导出失败：{exc}')

        # Annotation score margin / confidence UMAP
        if self.params.get('show_annotation_score_umap', True) and 'X_umap' in adata.obsm:
            score_col = None
            score_label = None
            if 'annotation_score_margin' in adata.obs.columns:
                score_col = 'annotation_score_margin'
                score_label = 'Annotation score margin'
            elif 'annotation_confidence' in adata.obs.columns:
                score_col = 'annotation_confidence'
                score_label = 'Annotation confidence'
            if score_col:
                fig_score = self.build_publication_umap(
                    adata, score_col, title=score_label + ' on UMAP'
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_score, plots_dir, 'annotation_score_umap.png',
                    'umap', score_label + ' UMAP', formats=('png', 'svg'), dpi=300,
                ))

        if maturity_evidence is not None and 'X_umap' in adata.obsm:
            try:
                fig_maturity = self.build_publication_umap(
                    adata, 'organoid_maturity_index',
                    title='Organoid Maturity Index on UMAP',
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_maturity, plots_dir, 'annotation_organoid_maturity_umap.png',
                    'umap', 'Organoid Maturity Index UMAP', formats=('png', 'svg'), dpi=300,
                ))
            except Exception as exc:
                self.progress(-1, f'成熟度 UMAP 导出失败：{exc}')

        # Confidence filtering and explicitly confirmed merges can modify
        # ``celltype`` after multi-evidence fields were created. Keep the
        # review-facing labels synchronized while preserving each source
        # cluster in a separate column and review-table row.
        if 'celltype' in adata.obs.columns:
            final_values = adata.obs['celltype'].astype(str)
            adata.obs['final_annotation'] = final_values
            doublet_flag = adata.obs.get(
                'annotation_doublet_status', pd.Series('', index=adata.obs.index)
            ).astype(str).eq('qc_predicted_doublet')
            mixture_flag = adata.obs.get(
                'annotation_lineage_mixture_status',
                pd.Series('', index=adata.obs.index),
            ).astype(str).eq('mixed_lineage_review')
            ambient_flag = adata.obs.get(
                'annotation_ambient_status', pd.Series('', index=adata.obs.index)
            ).astype(str).eq('ambient_signal')
            llm_low_confidence = adata.obs.get(
                'llm_confidence', pd.Series('', index=adata.obs.index)
            ).astype(str).str.lower().isin({'low', 'unknown'})
            review_status = []
            for unknown, qc_doublet, mixture, ambient, llm_low in zip(
                final_values.eq('Unknown'), doublet_flag, mixture_flag, ambient_flag,
                llm_low_confidence,
            ):
                if unknown:
                    review_status.append('Unknown')
                    continue
                reasons = []
                if qc_doublet:
                    reasons.append('qc_doublet')
                if mixture:
                    reasons.append('lineage_mixture')
                if ambient:
                    reasons.append('ambient')
                if llm_low:
                    reasons.append('llm_low_confidence')
                review_status.append('review_' + '_'.join(reasons) if reasons else 'review')
            adata.obs['annotation_status'] = review_status
            hierarchy = [
                hierarchy_for_label(label, marker_set_name, organoid_type)
                for label in final_values
            ]
            adata.obs['cell_lineage'] = [item[0] for item in hierarchy]
            adata.obs['cell_type_l1'] = [item[0] for item in hierarchy]
            adata.obs['cell_type_l2'] = [item[1] for item in hierarchy]
            adata.obs['cell_type_l3'] = [item[2] for item in hierarchy]
            adata.obs['cell_ontology_id'] = [item[3] for item in hierarchy]
            # 与需求文档中 celltype_l1/l2/l3 命名一致的别名列；cell_type_* 是
            # 平台既有列名，两者保持同步，兼容不同读取方。
            adata.obs['celltype_l1'] = [item[0] for item in hierarchy]
            adata.obs['celltype_l2'] = [item[1] for item in hierarchy]
            adata.obs['celltype_l3'] = [item[2] for item in hierarchy]
            adata.obs['annotation_source'] = (
                'manual' if str(method) == 'manual'
                else 'llm_cluster_evidence' if llm_evidence is not None
                else 'marker_evidence'
            )
            # Keep the source cluster explicit even when several clusters
            # legitimately receive the same biological cell-type label.  This
            # makes a repeated label auditable rather than silently treating
            # it as one combined annotation group.
            adata.obs['annotation_source_cluster'] = (
                str(leiden_key) + '=' + adata.obs[leiden_key].astype(str)
            ).astype('category')
            if maturity_evidence is not None:
                adata.obs['developmental_state'] = adata.obs['organoid_maturity_state'].astype(str)
            else:
                adata.obs['developmental_state'] = 'not_assessed'

        # Annotation computations above intentionally used the temporary
        # expression matrix.  Restore the input analysis representation before
        # writing h5ad so downstream embeddings and normalization provenance do
        # not silently change.
        adata.X = original_expression
        if had_temporary_source:
            adata.uns['_annotation_temporary_expression_source'] = previous_temporary_source
        else:
            adata.uns.pop('_annotation_temporary_expression_source', None)
        panel_anchors = (
            COLORECTAL_REFINED_ANCHOR_MARKERS
            if marker_set_name == 'Colorectal_refined'
            else derive_panel_anchor_markers(markers if isinstance(markers, dict) else {})
        )
        marker_panel_manifest = build_marker_panel_manifest(
            marker_set_name, organoid_type, markers, panel_anchors,
            negative_markers, self.params,
        )
        adata.uns['marker_panel_manifest'] = marker_panel_manifest
        evidence_tier_counts = {}
        if 'annotation_evidence_tier' in adata.obs.columns:
            evidence_tier_counts = {
                str(key): int(value)
                for key, value in adata.obs['annotation_evidence_tier'].astype(str).value_counts().items()
            }
        adata.uns['annotation_metadata'] = {
            'annotation_version': annotation_version,
            'panel_version': marker_panel_manifest['panel_version'],
            'annotation_comment': annotation_comment,
            'method': requested_method,
            'requested_marker_set': requested_marker_set_name,
            'marker_set': marker_set_name,
            'marker_panel_selection': marker_panel_selection,
            'fine_annotation': bool(fine_annotation),
            'fine_mode_applied': bool(fine_mode),
            'annotate_all': bool(annotate_all),
            'curated_display_labels_applied_cells': int(curated_display_labels_applied),
            'curated_display_label_map': (
                COLORECTAL_REFINED_DISPLAY_LABELS
                if marker_set_name == 'Colorectal_refined' else {}
            ),
            'unknown_recovery': {
                'cell_vote_recovered_clusters': sum(
                    1 for values in marker_decisions.get('cluster_evidence', {}).values()
                    if values.get('decision_reason') == 'unknown_recovered_by_cell_vote'
                ),
                'llm_unknown_marker_fallback_cells': int(
                    adata.obs['llm_unknown_fallback'].sum()
                ) if 'llm_unknown_fallback' in adata.obs else 0,
            },
            'organoid_type': (
                'intestinal' if marker_set_name in {'Colorectal', 'Colorectal_refined'}
                else organoid_type if marker_set_name == 'Organoid' else None
            ),
            'cluster_key': leiden_key,
            'cluster_merge_policy': 'manual_confirmed_only',
            'legacy_merge_similar_threshold_requested': legacy_merge_similar_threshold,
            'confirmed_cluster_merges_json': json.dumps(
                confirmed_cluster_merges, ensure_ascii=False, sort_keys=True,
            ),
            'annotation_expression_source': annotation_expression_source,
            'source': (
                'manual' if method == 'manual'
                else 'llm_cluster_evidence' if llm_evidence is not None
                else 'marker_evidence'
            ),
            'doublet_evidence_source': quality_evidence.get('doublet_source'),
            'ambient_evidence_rule': 'foreign_lineage_marker_detection',
            'celltypist_reference_enabled': bool(celltypist_evidence is not None),
            'celltypist_model': (
                celltypist_evidence.get('model_name') if celltypist_evidence else None
            ),
            'llm_annotation_enabled': bool(llm_evidence is not None),
            'llm_annotation_model': (
                llm_evidence.get('model') if llm_evidence else None
            ),
            'llm_annotation_prompt_version': (
                llm_evidence.get('prompt_version') if llm_evidence else None
            ),
            'llm_annotation_request_hash': (
                llm_evidence.get('request_hash') if llm_evidence else None
            ),
            'evidence_tier_counts': evidence_tier_counts,
            'evidence_tiers': list(EVIDENCE_TIERS),
            'manual_map_csv_applied': {
                'clusters': sorted(str(cluster) for cluster in manual_map_applied),
                'n_clusters': len(manual_map_applied),
            },
        }
        adata.uns['annotation_cluster_decisions'] = marker_decisions.get(
            'cluster_evidence', {}
        )

        # 邻接证据只用于人工复核，不参与自动合并。
        neighborhood_evidence = compute_cluster_neighborhood_evidence(adata, leiden_key)
        cluster_review = build_annotation_cluster_review(
            adata, leiden_key, marker_selection,
            neighborhood_evidence=neighborhood_evidence,
        )
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        if not cluster_review.empty:
            review_path = os.path.join(results_dir, 'annotation_cluster_review.csv')
            cluster_review.to_csv(review_path, index=False)
            result_files.append({
                'file_path': review_path, 'file_type': 'csv', 'category': 'table',
                'label': 'Annotation Cluster Review',
            })
            # 可下载、可复用的人工审阅模板：预填自动推荐标签，用户只需修改
            # 少数争议 cluster 的 manual_label，再通过 manual_map_csv 参数应用。
            template_columns = [
                column for column in (
                    'cluster', 'n_cells', 'final_annotation', 'marker_programme',
                    'display_label_source', 'evidence_tier',
                    'top_markers', 'broad_candidate', 'top_candidate',
                    'candidate_marker_detection_fraction',
                    'candidate_marker_specificity_delta',
                    'negative_marker_conflict', 'final_decision_reason',
                    'umap_adjacent_clusters', 'knn_adjacent_clusters',
                ) if column in cluster_review.columns
            ]
            template = cluster_review[template_columns].copy()
            template = template.rename(columns={'final_annotation': 'recommended_label'})
            if 'manual_label' not in template.columns:
                template.insert(
                    template.columns.get_loc('recommended_label') + 1,
                    'manual_label', 'KEEP',
                )
            template['notes'] = (
                'manual_label 默认 KEEP = 保持推荐标签；把要改的 cluster 的 '
                'manual_label 改成目标标签后保存 CSV，并在注释任务参数 '
                'manual_map_csv 中填写该文件路径即可应用（只更新被填写的 cluster）。'
            )
            template_path = os.path.join(
                results_dir, 'annotation_manual_map_template.csv',
            )
            template.to_csv(template_path, index=False)
            result_files.append({
                'file_path': template_path, 'file_type': 'csv', 'category': 'table',
                'label': 'Annotation Manual Map Template',
            })
        if markers or str(method) != 'manual':
            manifest_path = os.path.join(
                results_dir, 'annotation_marker_panel_manifest.json',
            )
            with open(manifest_path, 'w', encoding='utf-8') as handle:
                json.dump(
                    marker_panel_manifest, handle, ensure_ascii=False, indent=2,
                )
            result_files.append({
                'file_path': manifest_path, 'file_type': 'json', 'category': 'table',
                'label': 'Annotation Marker Panel Manifest',
            })

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'annotation')
        # 保留版本化副本：schema 承诺“不会覆盖输出中的旧 h5ad”。
        # 标准名 annotation_output.h5ad 继续被下游流水线读取，同时把本次
        # 注释结果按 annotation_version 另存一份，v2 不会覆盖 v1。
        try:
            import re as _re
            safe_version = _re.sub(r'[^A-Za-z0-9_.-]', '_', str(annotation_version))
            if safe_version and safe_version.lower() not in {'annotation_output'}:
                versioned_name = f'annotation_{safe_version}_output.h5ad'
                intermediate_dir = os.path.dirname(output_path)
                adata.write_h5ad(os.path.join(intermediate_dir, versioned_name))
        except Exception as exc:
            self.progress(-1, f'版本化注释副本保存失败（不影响主输出）: {exc}')

        ct_counts = adata.obs['celltype'].value_counts().to_dict()
        final_label_values = adata.obs['celltype'].astype(str)
        unresolved_label_mask = final_label_values.str.startswith('Unresolved ') | final_label_values.eq('Unknown')
        n_resolved_labels = int(final_label_values.loc[~unresolved_label_mask].nunique())
        cluster_label_modes = adata.obs.groupby(
            leiden_key, observed=True, sort=False,
        )['celltype'].agg(lambda series: series.astype(str).mode().iloc[0])
        n_unresolved_clusters = int((
            cluster_label_modes.str.startswith('Unresolved ')
            | cluster_label_modes.eq('Unknown')
        ).sum())
        label_reduction_warning = None
        if marker_set_name in {'Colorectal_refined', 'Organoid'} and len(cluster_label_modes) >= 4:
            # 精细面板预期标签数接近来源 cluster 数；只有明显压缩（少于
            # 来源簇数的 60%）才主动警告，避免对健康的 6/7 配置误报。
            collapse_threshold = max(2, int(np.ceil(len(cluster_label_modes) * 0.6)))
            if n_resolved_labels < collapse_threshold:
                message = (
                    f'精细注释出现标签压缩：{len(cluster_label_modes)} 个来源 cluster 只得到 '
                    f'{n_resolved_labels} 个细分标签（{n_unresolved_clusters} 个未解析）。',
                    '系统不会自动合并来源 cluster；请下载 annotation_cluster_review.csv 复核。',
                )
                label_reduction_warning = {
                    'triggered': True,
                    'n_source_clusters': int(len(cluster_label_modes)),
                    'n_resolved_labels': n_resolved_labels,
                    'n_unresolved_clusters': n_unresolved_clusters,
                    'message': ' '.join(message),
                }
                runtime_warnings.append(label_reduction_warning['message'])
        self.progress(100, "Done")
        result = {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_celltypes': adata.obs['celltype'].nunique(),
                'celltype_counts': {str(k): int(v) for k, v in ct_counts.items()},
                'cluster_column': leiden_key,
                'requested_cluster_key': requested_cluster_key,
                'n_source_clusters': int(adata.obs[leiden_key].nunique()),
                'n_final_labels': int(adata.obs['celltype'].nunique()),
                'n_resolved_labels': n_resolved_labels,
                'n_unresolved_clusters': n_unresolved_clusters,
                'label_reduction_warning': label_reduction_warning,
                'n_celltype_l1': int(adata.obs['cell_type_l1'].nunique()) if 'cell_type_l1' in adata.obs else 0,
                'n_celltype_l2': int(adata.obs['cell_type_l2'].nunique()) if 'cell_type_l2' in adata.obs else 0,
                'n_celltype_l3': int(adata.obs['cell_type_l3'].nunique()) if 'cell_type_l3' in adata.obs else 0,
                'panel_version': marker_panel_manifest['panel_version'],
                'evidence_tier_counts': evidence_tier_counts,
                'manual_map_csv_applied': {
                    'clusters': sorted(str(cluster) for cluster in manual_map_applied),
                    'n_clusters': len(manual_map_applied),
                },
                'cluster_merge_policy': 'manual_confirmed_only',
                'confirmed_cluster_merges': confirmed_cluster_merges,
                'method_used': requested_method,
                'annotation_version': annotation_version,
                'annotation_comment': annotation_comment,
                'annotation_expression_source': annotation_expression_source,
                'requested_marker_set': requested_marker_set_name,
                'marker_set': marker_set_name if method == 'auto_marker' or multi_evidence else None,
                'marker_panel_selection': marker_panel_selection,
                'fine_annotation': bool(fine_annotation),
                'fine_mode_applied': bool(fine_mode),
                'annotate_all': bool(annotate_all),
                'n_cell_vote_recovered': sum(
                    1 for values in marker_decisions.get('cluster_evidence', {}).values()
                    if values.get('decision_reason') == 'unknown_recovered_by_cell_vote'
                ),
                'n_llm_unknown_fallback': int(
                    adata.obs['llm_unknown_fallback'].sum()
                ) if 'llm_unknown_fallback' in adata.obs else 0,
                'organoid_type': (
                    'intestinal' if marker_set_name in {'Colorectal', 'Colorectal_refined'}
                    else organoid_type if marker_set_name == 'Organoid' else None
                ),
                'organoid_type_label': (
                    ORGANOID_MARKER_SET_LABELS.get(organoid_type, organoid_type)
                    if marker_set_name == 'Organoid'
                    else 'Colorectal/intestinal organoid'
                    if marker_set_name in {'Colorectal', 'Colorectal_refined'} else None
                ),
                'celltypist_reference': (
                    {
                        'enabled': True,
                        'model_name': celltypist_evidence.get('model_name'),
                        'model_path': celltypist_evidence.get('model_path'),
                        'package_version': celltypist_evidence.get('package_version'),
                        'mode': celltypist_evidence.get('mode'),
                        'p_threshold': celltypist_evidence.get('p_thres'),
                        'majority_voting': celltypist_evidence.get('majority_voting'),
                        'expression_source': celltypist_evidence.get('expression_source'),
                        'label_counts': celltypist_evidence.get('label_counts', {}),
                        'status_counts': celltypist_evidence.get('status_counts', {}),
                        'comparison_counts': celltypist_evidence.get('comparison_counts', {}),
                    }
                    if celltypist_evidence is not None else None
                ),
                'llm_annotation': (
                    {
                        'enabled': True,
                        'provider': llm_evidence.get('provider'),
                        'model': llm_evidence.get('model'),
                        'prompt_version': llm_evidence.get('prompt_version'),
                        'request_hash': llm_evidence.get('request_hash'),
                        'tissue_context_provided': bool(
                            llm_evidence.get('tissue_context_provided')
                        ),
                        'sent_data_policy': llm_evidence.get('sent_data_policy'),
                        'annotations': llm_evidence.get('annotations', []),
                        'global_note': llm_evidence.get('global_note', ''),
                    }
                    if llm_evidence is not None else None
                ),
                'maturity_evidence': (
                    {
                        'time_key': maturity_evidence.get('time_key'),
                        'expression_source': maturity_evidence.get('expression_source'),
                        'time_trend_spearman': (
                            float(maturity_evidence.get('time_trend_spearman'))
                            if np.isfinite(maturity_evidence.get('time_trend_spearman', np.nan)) else None
                        ),
                        'marker_coverage': maturity_evidence.get('marker_coverage', {}),
                        'state_counts': {
                            str(key): int(value)
                            for key, value in pd.Series(
                                maturity_evidence.get('maturity_state', [])
                            ).value_counts().items()
                        },
                        'warnings': maturity_evidence.get('warnings', []),
                    }
                    if maturity_evidence is not None else None
                ),
                'cell_state_evidence': {
                    'threshold': state_evidence.get('threshold', 0.35),
                    'coverage': state_evidence.get('coverage', {}),
                    'expression_source': state_evidence.get('expression_source', 'adata.X'),
                    'dominant_state_column': 'dominant_cell_state',
                    'multi_label_column': 'cell_state_flags',
                    'dominant_state_counts': {
                        str(key): int(value)
                        for key, value in adata.obs['dominant_cell_state'].astype(str).value_counts().items()
                    },
                    'multi_label_state_counts': {
                        str(state_name): int(np.asarray(flags, dtype=bool).sum())
                        for state_name, flags in state_evidence.get('high_flags', {}).items()
                    },
                    'n_cells_with_multiple_states': int(np.sum(
                        np.column_stack(list(state_evidence.get('high_flags', {}).values())).sum(axis=1) > 1
                    )) if state_evidence.get('high_flags') else 0,
                    # Backward-compatible alias: explicitly mirrors the
                    # exclusive dominant-state view.
                    'state_counts': {
                        str(key): int(value)
                        for key, value in adata.obs['cell_state'].astype(str).value_counts().items()
                    },
                },
                'marker_coverage': marker_coverage,
                'donor_key': donor_key,
                'condition_key': annotation_condition_key,
                'donor_celltype_composition': annotation_donor_composition,
                'negative_marker_coverage': negative_marker_coverage,
                'negative_marker_weight': negative_marker_weight,
                'cluster_decisions': marker_decisions.get('cluster_evidence', {}),
                'marker_selection': {
                    'cluster_markers': marker_selection.get('cluster_markers', {}),
                    'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
                    'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
                    'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
                    'warnings': marker_selection.get('warnings', []),
                    'parameters': marker_selection.get('parameters', {}),
                },
                'marker_validation': {
                    'genes': marker_validation.get('genes', []) if marker_validation else [],
                    'categories': marker_validation.get('categories', []) if marker_validation else [],
                    'expression_source': marker_validation.get('expression_source', '') if marker_validation else '',
                    'display_warning': marker_validation.get('display_warning', '') if marker_validation else '',
                    'n_cells_plotted': int(marker_validation.get('n_cells_plotted', 0)) if marker_validation else 0,
                },
                'runtime_warnings': list(dict.fromkeys(runtime_warnings)),
                'agreement_source': adata.uns.get('annotation_agreement_source', 'marker_label'),
                'quality_evidence': {
                    'doublet_source': quality_evidence.get('doublet_source', 'not_available'),
                    'qc_doublets': adata.uns.get('qc_doublets') or adata.uns.get('qc_scrublet', {}),
                    'lineage_mixture_threshold': float(self.params.get(
                        'lineage_mixture_threshold',
                        self.params.get('doublet_score_threshold', 0.30),
                    )),
                    'ambient_threshold': float(self.params.get('ambient_score_threshold', 0.35)),
                    'ambient_prevalence': float(self.params.get('ambient_prevalence', 0.50)),
                    'ambient_genes': quality_evidence.get('ambient_genes', []),
                    'ambient_genes_by_label': quality_evidence.get('ambient_genes_by_label', {}),
                    'ambient_expression_source': quality_evidence.get('ambient_expression_source', ''),
                    'n_qc_predicted_doublet': int(
                        (adata.obs['annotation_doublet_status'].astype(str) == 'qc_predicted_doublet').sum()
                    ) if 'annotation_doublet_status' in adata.obs else 0,
                    # Backward-compatible count, now strictly QC-derived.
                    'n_suspect_doublet': int(
                        (adata.obs['annotation_doublet_status'].astype(str) == 'qc_predicted_doublet').sum()
                    ) if 'annotation_doublet_status' in adata.obs else 0,
                    'n_lineage_mixture_review': int(
                        (adata.obs['annotation_lineage_mixture_status'].astype(str) == 'mixed_lineage_review').sum()
                    ) if 'annotation_lineage_mixture_status' in adata.obs else 0,
                    'n_ambient_signal': int(
                        (adata.obs['annotation_ambient_status'].astype(str) == 'ambient_signal').sum()
                    ) if 'annotation_ambient_status' in adata.obs else 0,
                },
                'n_review_clusters': int(cluster_review['needs_review'].sum()) if not cluster_review.empty else 0,
                'n_clusters_reviewed': int(len(cluster_review)),
            }
        }
        if 'annotation_confidence' in adata.obs.columns:
            confidence_value = round(float(adata.obs['annotation_confidence'].mean()), 3)
            result['summary']['mean_confidence'] = confidence_value
        if 'annotation_score_margin' in adata.obs.columns:
            result['summary']['mean_score_margin'] = round(float(adata.obs['annotation_score_margin'].mean()), 3)
        return result
