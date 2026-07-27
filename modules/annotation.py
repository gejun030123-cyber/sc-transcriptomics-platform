import os

from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping


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
    """Resolve a downloaded CellTypist model without triggering a download."""
    requested = str(model_name or 'Immune_All_Low.pkl').strip()
    model_dir = str(model_dir or CELLTYPIST_MODEL_DIR)
    if os.path.isabs(requested):
        candidate = requested
    else:
        # Only a basename is accepted for the configured model directory.  A
        # user may still provide an explicit absolute path for a custom model.
        candidate = os.path.join(model_dir, os.path.basename(requested))
    if not os.path.isfile(candidate):
        raise FileNotFoundError(
            f"CellTypist 模型不存在: {candidate}。请先将选定的 .pkl 放入 {model_dir}。"
        )
    return os.path.realpath(candidate)


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
    celltypist = importlib.import_module('celltypist')
    model_path = resolve_celltypist_model(model_name, requested_model_dir)

    # Avoid copying the full AnnData object (including large obsm/uns fields).
    # CellTypist only needs a cell-by-gene matrix and gene symbols.
    from anndata import AnnData
    if 'counts' in adata.layers:
        matrix = adata.layers['counts'].copy()
        expression_source = 'layers["counts"] → library-size normalized log1p'
    else:
        matrix = adata.X.copy() if hasattr(adata.X, 'copy') else adata.X
        expression_source = 'adata.X'
    query = AnnData(
        X=matrix,
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        var=pd.DataFrame(index=adata.var_names.copy()),
    )
    if 'counts' in adata.layers:
        import scanpy as sc
        sc.pp.normalize_total(query, target_sum=10000.0, inplace=True)
        sc.pp.log1p(query)
    else:
        values = query.X.data if hasattr(query.X, 'data') else np.asarray(query.X)
        if np.nanmin(values) < 0:
            raise ValueError(
                'CellTypist 需要非负 counts 或 log1p 表达矩阵；当前 adata.X 含负值，'
                '请从保存 counts 层的结果重新运行。'
            )

    mode = str(mode or 'prob match').strip().lower()
    mode = 'prob match' if mode in {'prob match', 'prob_match', 'probability'} else 'best match'
    p_thres = min(1.0, max(0.0, float(p_thres)))
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
    if str(marker_set_name) == 'Organoid':
        return ORGANOID_MARKER_SETS.get(str(organoid_type).strip().lower(),
                                        ORGANOID_MARKER_SETS['intestinal'])
    return MARKER_SETS.get(marker_set_name, DEFAULT_UNIVERSAL_MARKERS)


def get_negative_marker_set(marker_set_name='Universal', organoid_type='intestinal'):
    """Return conservative anti-marker modules for the selected tissue panel."""
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

    expression_layer = 'counts' if 'counts' in adata.layers else None
    raw_scores = {
        module: _mean_expression_for_genes(adata, genes, layer=expression_layer)
        for module, genes in matched.items()
    }

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
        'expression_source': (
            f'layers["{expression_layer}"]' if expression_layer else 'adata.X'
        ),
        'time_labels': time_labels,
        'time_numeric': time_numeric,
        'time_trend_spearman': trend_rho,
        'marker_coverage': coverage,
        'warnings': warnings,
    }


def hierarchy_for_label(label, marker_set_name='Universal', organoid_type='intestinal'):
    """Return stable lineage/type/subtype fields without merging cell states."""
    label = str(label or 'Unknown')
    if label == 'Unknown':
        return ('Unknown', 'Unknown', 'Unknown', '')
    if marker_set_name == 'Universal':
        values = UNIVERSAL_LABELS.get(label)
        if values:
            return values
    if marker_set_name == 'Organoid':
        values = ORGANOID_HIERARCHY.get(str(organoid_type).lower(), {}).get(label)
        if values:
            return (*values, '')
    # Scenario panels do not yet have a curated ontology hierarchy; preserve
    # the label rather than fabricating a lineage or ontology ID.
    return ('Unresolved lineage', label, label, '')


def compute_cell_state_evidence(adata, threshold=0.35):
    """Compute orthogonal cell-state modules (cycling, stress, IFN, etc.)."""
    import numpy as np

    var_lookup = {str(gene).upper(): str(gene) for gene in adata.var_names}
    expression_layer = 'counts' if 'counts' in adata.layers else None
    scores = {}
    coverage = {}
    for state_name, genes in CELL_STATE_MARKERS.items():
        matched = list(dict.fromkeys(
            var_lookup[str(gene).upper()] for gene in genes
            if str(gene).upper() in var_lookup
        ))
        coverage[state_name] = {'matched': len(matched), 'total': len(genes)}
        values = _mean_expression_for_genes(adata, matched, layer=expression_layer)
        positive = values[np.isfinite(values) & (values > 0)]
        scale = float(np.nanpercentile(positive, 95)) if positive.size else 0.0
        scores[state_name] = np.clip(values / (scale + 1e-12), 0.0, 1.0) if scale else np.zeros(adata.n_obs)

    score_names = list(scores)
    score_matrix = np.column_stack([scores[name] for name in score_names]) if score_names else np.zeros((adata.n_obs, 0))
    top_index = score_matrix.argmax(axis=1) if score_matrix.shape[1] else np.zeros(adata.n_obs, dtype=int)
    top_score = score_matrix.max(axis=1) if score_matrix.shape[1] else np.zeros(adata.n_obs)
    primary = np.asarray([
        score_names[index] if value >= float(threshold) else 'baseline'
        for index, value in zip(top_index, top_score)
    ], dtype=object)
    flags = np.asarray([
        ';'.join(name for name in score_names if scores[name][index] >= float(threshold)) or 'baseline'
        for index in range(adata.n_obs)
    ], dtype=object)
    return {
        'scores': scores,
        'primary_state': primary,
        'state_flags': flags,
        'coverage': coverage,
        'threshold': float(threshold),
        'expression_source': f'layers["{expression_layer}"]' if expression_layer else 'adata.X',
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


def compute_annotation_quality_evidence(
    adata,
    score_cols,
    marker_genes_by_type=None,
    *,
    doublet_threshold=0.30,
    ambient_threshold=0.35,
    ambient_prevalence=0.50,
):
    """Compute doublet and ambient-RNA evidence without forcing labels.

    The doublet score is high when two independent marker modules receive
    similar positive scores.  Ambient evidence is a heuristic based on marker
    genes detected in a large fraction of cells and should be interpreted as a
    review signal unless empty-droplet data are available.
    """
    import numpy as np

    score_cols = [str(col) for col in score_cols if col in adata.obs.columns]
    if not score_cols:
        return {
            'doublet_score': np.zeros(adata.n_obs),
            'doublet_status': np.repeat('not_evaluated', adata.n_obs),
            'ambient_score': np.zeros(adata.n_obs),
            'ambient_status': np.repeat('not_evaluated', adata.n_obs),
            'top1': np.repeat('', adata.n_obs), 'top2': np.repeat('', adata.n_obs),
            'top1_score': np.zeros(adata.n_obs), 'top2_score': np.zeros(adata.n_obs),
            'ambient_genes': [], 'warnings': ['没有可用 marker score，未计算 doublet/环境 RNA 证据。'],
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
    doublet_score = 2.0 * np.minimum(
        probabilities[np.arange(adata.n_obs), top1_idx],
        probabilities[np.arange(adata.n_obs), top2_idx],
    )
    # Two negative/near-background scores are ambiguity, not evidence of a
    # doublet.  Require both modules to have a positive score before raising
    # the doublet signal.
    positive_pair = (top1_score > 0) & (top2_score > 0)
    doublet_score = np.where(positive_pair, doublet_score, 0.0)
    if len(score_cols) < 2:
        doublet_score = np.zeros(adata.n_obs, dtype=float)
        doublet_status = np.repeat('not_evaluated', adata.n_obs)
    else:
        doublet_status = np.where(
            doublet_score >= float(doublet_threshold), 'suspect_doublet', 'single_or_uncertain'
        )

    marker_genes_by_type = marker_genes_by_type or {}
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
    global_fraction = _nonzero_fraction(adata[:, all_markers].X) if all_markers else np.array([])
    # A marker shared by several candidate types is not useful as ambient
    # evidence.  Restrict the heuristic to a lineage-specific marker that is
    # detected in many cells but is not part of the winning module.
    ambient_genes = [
        gene for gene, fraction in zip(all_markers, global_fraction)
        if float(fraction) >= float(ambient_prevalence)
        and len(marker_to_types.get(gene, set())) == 1
    ]
    if ambient_genes:
        ambient_matrix = adata[:, ambient_genes].X
        if hasattr(ambient_matrix, 'toarray'):
            ambient_matrix = ambient_matrix.toarray()
        ambient_matrix = np.asarray(ambient_matrix, dtype=float)
        top1_labels = np.asarray([score_cols[i].replace('score_', '') for i in top1_idx])
        ambient_values = np.zeros(adata.n_obs, dtype=float)
        for index, top1_label in enumerate(top1_labels):
            eligible = [
                gene_index for gene_index, gene in enumerate(ambient_genes)
                if top1_label not in marker_to_types.get(gene, set())
            ]
            if eligible:
                ambient_values[index] = float(np.mean(ambient_matrix[index, eligible]))
        max_value = float(np.nanpercentile(ambient_values, 95)) if np.any(ambient_values > 0) else 0.0
        ambient_score = np.clip(ambient_values / (max_value + 1e-12), 0.0, 1.0)
    else:
        ambient_score = np.zeros(adata.n_obs, dtype=float)
    ambient_status = np.where(
        ambient_score >= float(ambient_threshold), 'ambient_signal', 'no_strong_ambient_signal'
    )
    warnings = []
    if ambient_genes:
        warnings.append(
            '环境 RNA 仅按高普遍性 marker 做启发式筛查；若需定量去污染，请提供 empty droplets 或运行专门去污染步骤。'
        )
    return {
        'doublet_score': doublet_score,
        'doublet_status': doublet_status,
        'ambient_score': ambient_score,
        'ambient_status': ambient_status,
        'top1': np.asarray([score_cols[i].replace('score_', '') for i in top1_idx], dtype=object),
        'top2': np.asarray([score_cols[i].replace('score_', '') for i in top2_idx], dtype=object),
        'top1_score': top1_score,
        'top2_score': top2_score,
        'ambient_genes': ambient_genes,
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


def build_annotation_cluster_review(adata, cluster_key, marker_selection=None):
    """Create a compact, auditable cluster-level annotation review table."""
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
        label_agreement = float((celltype_values == celltype).mean()) if len(celltype_values) else np.nan
        row = {
            'cluster': str(cluster),
            'n_cells': int(mask.sum()),
            'assigned_celltype': celltype,
            'final_annotation': final_label,
            'annotation_status': annotation_status,
            'cluster_label_agreement': round(label_agreement, 4) if np.isfinite(label_agreement) else np.nan,
            'top_markers': ', '.join(marker_map.get(str(cluster), [])),
        }
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
        for source, target in (
            ('marker_score', 'mean_marker_score'),
            ('annotation_confidence', 'mean_confidence'),
            ('annotation_score_margin', 'mean_score_margin'),
            ('annotation_doublet_score', 'mean_doublet_score'),
            ('annotation_ambient_score', 'mean_ambient_score'),
        ):
            if source in frame.columns:
                row[target] = round(float(pd.to_numeric(frame[source], errors='coerce').mean()), 4)
            else:
                row[target] = np.nan
        row['needs_review'] = bool(
            celltype == 'Unknown'
            or final_label == 'Unknown'
            or annotation_status.lower() == 'unknown'
            or (np.isfinite(label_agreement) and label_agreement < 0.8)
            or (np.isfinite(row['mean_confidence']) and row['mean_confidence'] < 0.2)
            or (np.isfinite(row['mean_score_margin']) and row['mean_score_margin'] < 0.05)
            or annotation_status.lower() in {
                'review_doublet', 'review_ambient', 'review_doublet_ambient',
            }
        )
        if 'celltypist_conflict_fraction' in row:
            row['needs_review'] = bool(
                row['needs_review']
                or row['celltypist_conflict_fraction'] > 0
                or row['celltypist_unassigned_fraction'] > 0
                or row['celltypist_multi_label_fraction'] > 0
            )
        rows.append(row)
    return pd.DataFrame(rows)

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
        requested_cluster_key = str(self.params.get('cluster_key', 'leiden') or '').strip()
        resolution = self.params.get('resolution', '0.8')
        preferred_cluster_key = f'leiden_{resolution}' if f'leiden_{resolution}' in adata.obs.columns else requested_cluster_key
        leiden_key, cluster_info = resolve_obs_grouping(
            adata, preferred_cluster_key,
            fallbacks=[requested_cluster_key, 'leiden', 'leiden_0.8', 'leiden_0.6', 'leiden_1.0'],
            max_categories=50, max_numeric_categories=20,
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
        supported_methods = {'multi_evidence', 'auto_marker', 'manual'}
        if requested_method not in supported_methods:
            raise ValueError(
                f"不支持的注释方法 '{requested_method}'；当前仅支持 "
                "multi_evidence、auto_marker 和 manual。"
            )
        method = requested_method
        annotation_version = str(self.params.get('annotation_version', 'v1') or 'v1').strip() or 'v1'
        annotation_comment = str(self.params.get('annotation_comment', '') or '').strip()
        multi_evidence = method == 'multi_evidence'
        if multi_evidence:
            method = 'auto_marker'
        use_celltypist_reference = bool(self.params.get('use_celltypist_reference', False))
        if use_celltypist_reference and requested_method == 'manual':
            runtime_warnings.append('CellTypist 参考交叉验证仅在 auto_marker/multi_evidence 中启用；manual 模式已跳过。')
            use_celltypist_reference = False
        marker_set_name = str(self.params.get('marker_set', 'Universal') or 'Universal')
        organoid_type = str(self.params.get('organoid_type', 'intestinal') or 'intestinal').strip().lower()
        custom_markers_str = self.params.get('custom_markers', '').strip()
        confidence_method = self.params.get('confidence_method', 'none')
        # Confidence is a review signal by default.  Automatic relabelling is
        # opt-in and, in multi-evidence mode, additionally requires an
        # internally inconsistent cluster so a coherent cluster is not erased
        # merely because entropy is low across a broad universal panel.
        mark_unknown = self.params.get('mark_unknown', False)
        confidence_cutoff = float(self.params.get('confidence_cutoff', 0.2))
        confidence_cutoff = min(1.0, max(0.0, confidence_cutoff))
        merge_similar_threshold = float(self.params.get('merge_similar_threshold', 0))
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

        if method == 'manual' and custom_markers_str:
            self.progress(40, "Applying manual cell type mapping...")
            manual_mapping = {}
            for cluster_id, cell_type in parse_marker_mapping(custom_markers_str).items():
                manual_mapping[cluster_id] = cell_type

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
                    markers = get_marker_set(marker_set_name, organoid_type)
            else:
                markers = get_marker_set(marker_set_name, organoid_type)

            negative_markers = get_negative_marker_set(marker_set_name, organoid_type)
            if negative_marker_text:
                custom_negative = parse_marker_mapping(negative_marker_text)
                if custom_negative:
                    negative_markers.update(custom_negative)

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
                cluster_annotations = {}
                for cluster in adata.obs[leiden_key].cat.categories:
                    mask = adata.obs[leiden_key] == cluster
                    mean_scores = {col: adata.obs.loc[mask, col].mean() for col in score_cols}
                    best_col = max(mean_scores, key=mean_scores.get)
                    best_ct = best_col.replace('score_', '')
                    cluster_annotations[cluster] = best_ct
                adata.obs['celltype'] = adata.obs[leiden_key].map(cluster_annotations).astype('category')
                min_score = float(self.params.get('min_annotation_score', 0.0))
                if min_score > 0:
                    best_scores = adata.obs[score_cols].max(axis=1)
                    if 'Unknown' not in adata.obs['celltype'].cat.categories:
                        adata.obs['celltype'] = adata.obs['celltype'].cat.add_categories(['Unknown'])
                    adata.obs.loc[best_scores < min_score, 'celltype'] = 'Unknown'
            else:
                adata.obs['celltype'] = adata.obs[leiden_key].astype(str)
            markers = usable_markers
        else:
            marker_coverage = {}
            negative_marker_coverage = {}
            negative_markers = {}

        # Quality evidence is recorded for both auto_marker and
        # multi_evidence.  It is deliberately a review signal: a suspected
        # doublet or ambient signal is not automatically discarded.
        score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
        quality_evidence = compute_annotation_quality_evidence(
            adata,
            score_cols,
            markers if isinstance(markers, dict) else {},
            doublet_threshold=float(self.params.get('doublet_score_threshold', 0.30)),
            ambient_threshold=float(self.params.get('ambient_score_threshold', 0.35)),
            ambient_prevalence=float(self.params.get('ambient_prevalence', 0.50)),
        )
        adata.obs['annotation_doublet_score'] = quality_evidence['doublet_score']
        adata.obs['annotation_doublet_status'] = quality_evidence['doublet_status']
        adata.obs['annotation_ambient_score'] = quality_evidence['ambient_score']
        adata.obs['annotation_ambient_status'] = quality_evidence['ambient_status']
        adata.obs['annotation_top1'] = quality_evidence['top1']
        adata.obs['annotation_top2'] = quality_evidence['top2']
        adata.obs['annotation_top1_score'] = quality_evidence['top1_score']
        adata.obs['annotation_top2_score'] = quality_evidence['top2_score']
        if quality_evidence.get('warnings'):
            runtime_warnings.extend(quality_evidence['warnings'])

        maturity_evidence = None
        if marker_set_name == 'Organoid':
            maturity_evidence = compute_organoid_maturity_evidence(
                adata,
                organoid_type=organoid_type,
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
        adata.obs['cell_state'] = state_evidence['primary_state']
        adata.obs['cell_state_flags'] = state_evidence['state_flags']
        adata.uns['cell_state_evidence'] = {
            'coverage': state_evidence['coverage'],
            'threshold': state_evidence['threshold'],
            'expression_source': state_evidence['expression_source'],
        }

        # Multi-evidence mode keeps the rule engine interpretable and records
        # independent per-cell marker predictions for cluster-level review.
        if multi_evidence:
            adata.obs['marker_label'] = adata.obs['celltype'].astype(str)
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if score_cols:
                score_frame = adata.obs[score_cols]
                adata.obs['marker_score'] = score_frame.max(axis=1).astype(float)
                ranked = np.argsort(score_frame.values, axis=1)
                adata.obs['prediction_1'] = [score_cols[i].replace('score_', '') for i in ranked[:, -1]]
                adata.obs['prediction_1_score'] = score_frame.max(axis=1).astype(float)
                adata.obs['prediction_2'] = [score_cols[i].replace('score_', '') for i in ranked[:, -2]] if len(score_cols) > 1 else ''
                if len(score_cols) > 1:
                    adata.obs['prediction_2_score'] = score_frame.values[np.arange(adata.n_obs), ranked[:, -2]]
            # ``marker_label`` is assigned once per cluster by the marker
            # scorer, so it is always perfectly coherent by construction.  In
            # multi-evidence mode use the independent per-cell top prediction
            # when available; this makes the agreement score actually detect a
            # mixed/ambiguous cluster instead of reporting 1.0 for every one.
            agreement_key = 'prediction_1' if 'prediction_1' in adata.obs.columns else 'marker_label'
            cluster_vote = adata.obs.groupby(leiden_key, observed=True)[agreement_key].transform(
                lambda values: values.value_counts(normalize=True).iloc[0]
            )
            adata.obs['cluster_annotation_agreement'] = cluster_vote.astype(float)
            adata.uns['annotation_agreement_source'] = agreement_key
            adata.obs['annotation_status'] = np.where(
                (adata.obs.get('marker_score', 0) < float(self.params.get('min_annotation_score', 0.0))) |
                (adata.obs['cluster_annotation_agreement'] < float(self.params.get('cluster_agreement_threshold', 0.6))),
                'Unknown', 'review'
            )
            adata.obs['final_annotation'] = adata.obs['marker_label'].astype(str)
            adata.obs.loc[adata.obs['annotation_status'] == 'Unknown', 'final_annotation'] = 'Unknown'
            adata.obs['celltype'] = adata.obs['final_annotation'].astype('category')

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

        # 相似簇合并
        if merge_similar_threshold > 0:
            import numpy as np
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if score_cols:
                cluster_profiles = {}
                for cluster in adata.obs[leiden_key].cat.categories:
                    mask = adata.obs[leiden_key] == cluster
                    cluster_profiles[cluster] = adata.obs.loc[mask, score_cols].mean().values
                clusters = list(cluster_profiles.keys())
                for i in range(len(clusters)):
                    for j in range(i + 1, len(clusters)):
                        corr = np.corrcoef(cluster_profiles[clusters[i]], cluster_profiles[clusters[j]])[0, 1]
                        if corr >= merge_similar_threshold:
                            new_ct = adata.obs.loc[adata.obs[leiden_key] == clusters[i], 'celltype'].mode()
                            new_ct = new_ct.iloc[0] if not new_ct.empty else 'Unknown'
                            adata.obs.loc[adata.obs[leiden_key] == clusters[j], 'celltype'] = new_ct
                            self.progress(-1, f"合并簇 {clusters[j]} → {clusters[i]} (r={corr:.2f})")

        self.progress(70, "Generating dotplot and UMAP...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        # Marker dotplots are evidence displays, not fixed marker-panel
        # illustrations.  Rank genes independently for every cluster, then
        # keep a small number of data-driven markers plus at most one classic
        # anchor per cluster.  The returned metadata is written to summary so
        # the exact selection can be reviewed or reproduced.
        marker_selection = select_cluster_marker_genes(
            adata,
            leiden_key,
            classic_markers=markers if isinstance(markers, dict) else {},
            method=self.params.get('marker_selection_method', 'wilcoxon'),
            n_rank_genes=int(self.params.get('marker_rank_genes', 200)),
            min_markers=int(self.params.get('marker_min_per_cluster', 2)),
            max_markers=int(self.params.get('marker_max_per_cluster', 5)),
            padj_cutoff=float(self.params.get('marker_padj_cutoff', 0.05)),
            min_pct=float(self.params.get('marker_min_pct', 0.10)),
            min_delta_pct=float(self.params.get('marker_min_delta_pct', 0.05)),
        )
        dotplot_genes = marker_selection.get('genes', [])[:40]
        if marker_selection.get('warnings'):
            runtime_warnings.extend(marker_selection.get('warnings', []))
            self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
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
            coverage_types = list(marker_coverage)
            matched = [marker_coverage[ct]['matched'] for ct in coverage_types]
            totals = [marker_coverage[ct]['total'] for ct in coverage_types]
            fig_coverage = grouped_bar_figure(
                coverage_types,
                [('可用 marker', matched), ('marker 总数', totals)],
                title='Marker 覆盖度（当前数据）',
                x_label='Cell type', y_label='基因数',
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_coverage, plots_dir, 'annotation_marker_coverage.png',
                'bar', 'Marker 覆盖度', formats=('png', 'svg'), dpi=300,
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
                    y_labels=[c.replace('score_', '') for c in score_cols],
                    title=f'Marker Score Heatmap by {leiden_key}',
                    x_label='Cluster', y_label='Marker set',
                    colorbar_label='Row z-score',
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_score_heat, plots_dir, 'annotation_marker_score_heatmap.png',
                    'heatmap', 'Marker Score Heatmap', formats=('png', 'svg'), dpi=300,
                ))

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

        # Confidence filtering and optional cluster merging can modify
        # ``celltype`` after multi-evidence fields were created.  Keep the
        # review-facing labels synchronized so a CSV never reports a stale
        # final label as confirmed.
        if 'celltype' in adata.obs.columns:
            final_values = adata.obs['celltype'].astype(str)
            adata.obs['final_annotation'] = final_values
            base_status = np.where(
                final_values.eq('Unknown'), 'Unknown', 'review'
            )
            doublet_flag = adata.obs.get(
                'annotation_doublet_status', pd.Series('', index=adata.obs.index)
            ).astype(str).eq('suspect_doublet')
            ambient_flag = adata.obs.get(
                'annotation_ambient_status', pd.Series('', index=adata.obs.index)
            ).astype(str).eq('ambient_signal')
            adata.obs['annotation_status'] = np.where(
                final_values.eq('Unknown'), 'Unknown',
                np.where(doublet_flag & ambient_flag, 'review_doublet_ambient',
                         np.where(doublet_flag, 'review_doublet',
                                  np.where(ambient_flag, 'review_ambient', base_status)))
            )
            hierarchy = [
                hierarchy_for_label(label, marker_set_name, organoid_type)
                for label in final_values
            ]
            adata.obs['cell_lineage'] = [item[0] for item in hierarchy]
            adata.obs['cell_type_l1'] = [item[0] for item in hierarchy]
            adata.obs['cell_type_l2'] = [item[1] for item in hierarchy]
            adata.obs['cell_type_l3'] = [item[2] for item in hierarchy]
            adata.obs['cell_ontology_id'] = [item[3] for item in hierarchy]
            adata.obs['annotation_source'] = (
                'manual' if str(method) == 'manual' else 'marker_evidence'
            )
            if maturity_evidence is not None:
                adata.obs['developmental_state'] = adata.obs['organoid_maturity_state'].astype(str)
            else:
                adata.obs['developmental_state'] = 'not_assessed'

        adata.uns['annotation_metadata'] = {
            'annotation_version': annotation_version,
            'annotation_comment': annotation_comment,
            'method': requested_method,
            'marker_set': marker_set_name,
            'organoid_type': organoid_type if marker_set_name == 'Organoid' else None,
            'cluster_key': leiden_key,
            'source': 'manual' if method == 'manual' else 'marker_evidence',
            'celltypist_reference_enabled': bool(celltypist_evidence is not None),
            'celltypist_model': (
                celltypist_evidence.get('model_name') if celltypist_evidence else None
            ),
        }

        cluster_review = build_annotation_cluster_review(
            adata, leiden_key, marker_selection,
        )
        if not cluster_review.empty:
            results_dir = os.path.join(self.project_dir, 'results')
            os.makedirs(results_dir, exist_ok=True)
            review_path = os.path.join(results_dir, 'annotation_cluster_review.csv')
            cluster_review.to_csv(review_path, index=False)
            result_files.append({
                'file_path': review_path, 'file_type': 'csv', 'category': 'table',
                'label': 'Annotation Cluster Review',
            })

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'annotation')

        ct_counts = adata.obs['celltype'].value_counts().to_dict()
        self.progress(100, "Done")
        result = {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_celltypes': adata.obs['celltype'].nunique(),
                'celltype_counts': {str(k): int(v) for k, v in ct_counts.items()},
                'cluster_column': leiden_key,
                'requested_cluster_key': requested_cluster_key,
                'method_used': requested_method,
                'annotation_version': annotation_version,
                'annotation_comment': annotation_comment,
                'marker_set': marker_set_name if method == 'auto_marker' or multi_evidence else None,
                'organoid_type': organoid_type if marker_set_name == 'Organoid' else None,
                'organoid_type_label': (
                    ORGANOID_MARKER_SET_LABELS.get(organoid_type, organoid_type)
                    if marker_set_name == 'Organoid' else None
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
                    'state_counts': {
                        str(key): int(value)
                        for key, value in adata.obs['cell_state'].astype(str).value_counts().items()
                    },
                },
                'marker_coverage': marker_coverage,
                'negative_marker_coverage': negative_marker_coverage,
                'negative_marker_weight': negative_marker_weight,
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
                    'doublet_threshold': float(self.params.get('doublet_score_threshold', 0.30)),
                    'ambient_threshold': float(self.params.get('ambient_score_threshold', 0.35)),
                    'ambient_prevalence': float(self.params.get('ambient_prevalence', 0.50)),
                    'ambient_genes': quality_evidence.get('ambient_genes', []),
                    'n_suspect_doublet': int(
                        (adata.obs['annotation_doublet_status'].astype(str) == 'suspect_doublet').sum()
                    ) if 'annotation_doublet_status' in adata.obs else 0,
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
