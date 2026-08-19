"""Human-only statistical helpers for Bulk RNA-seq pathway enrichment.

The web module deliberately keeps these routines independent from OmicVerse so
ORA can be tested offline and so the exact background, multiple-testing family,
and identifier losses are explicit in exported results.
"""

from __future__ import annotations

import hashlib
import importlib.util
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


GENE_COLUMNS = (
    'gene', 'Gene', 'gene_symbol', 'GeneSymbol', 'gene_name', 'GeneName',
    'symbol', 'Symbol',
)
SYMBOL_COLUMNS = (
    'gene_symbol', 'GeneSymbol', 'gene_name', 'GeneName', 'symbol', 'Symbol',
)
ID_COLUMNS = (
    'gene_id', 'GeneID', 'gene_ids', 'ensembl_gene_id', 'ensembl_id',
)
ALIAS_COLUMNS = ('alias', 'aliases', 'gene_alias', 'gene_aliases', 'synonyms')
STATISTIC_COLUMNS = (
    'stat', 'statistic', 'wald_stat', 'waldStatistic', 't_stat', 't_statistic',
    'tvalue', 't_value', 'score',
)
LOG2FC_COLUMNS = ('log2FC', 'log2fc', 'logFC', 'logfc')
PVALUE_COLUMNS = ('pvalue', 'P-value', 'pval', 'PValue', 'padj')


def normalise_human_gene(value) -> str:
    """Return a canonical matching token for a human gene identifier."""
    if value is None:
        return ''
    try:
        if bool(pd.isna(value)):
            return ''
    except (TypeError, ValueError):
        pass
    text = str(value).strip().strip('"').strip("'")
    if not text or text.lower() in {'nan', 'none', 'null'}:
        return ''
    # Ensembl version suffixes are not part of the stable gene identifier.
    if re.fullmatch(r'ENSG\d+(?:\.\d+)?', text, flags=re.IGNORECASE):
        text = text.split('.', 1)[0]
    return text.upper()


def split_gene_text(value) -> list[str]:
    """Split comma/newline/semicolon separated identifiers without guessing prose."""
    if value is None:
        return []
    return [part.strip() for part in re.split(r'[,;\s]+', str(value)) if part.strip()]


def _first_column(frame: pd.DataFrame, candidates, *, exclude=()) -> str | None:
    excluded = set(exclude)
    return next((column for column in candidates
                 if column in frame.columns and column not in excluded), None)


def build_identifier_map(frame: pd.DataFrame, primary_column: str | None = None):
    """Build an offline ID/alias -> HGNC-symbol map from columns already supplied.

    No online annotation service is consulted.  This makes the mapping
    reproducible and ensures identifiers that cannot be resolved from the
    submitted data are reported rather than silently guessed.
    """
    primary_column = primary_column or _first_column(frame, GENE_COLUMNS)
    if primary_column is None:
        raise ValueError('结果表缺少 gene/gene_symbol 列。')
    symbol_column = _first_column(frame, SYMBOL_COLUMNS, exclude=(primary_column,))
    id_columns = [column for column in ID_COLUMNS if column in frame.columns]
    alias_columns = [column for column in ALIAS_COLUMNS if column in frame.columns]
    mapping: dict[str, str] = {}
    conflicts = set()

    for _, row in frame.iterrows():
        primary = normalise_human_gene(row.get(primary_column))
        symbol = normalise_human_gene(row.get(symbol_column)) if symbol_column else primary
        if not symbol:
            symbol = primary
        if not symbol:
            continue
        tokens = [primary, symbol]
        tokens.extend(normalise_human_gene(row.get(column)) for column in id_columns)
        for column in alias_columns:
            tokens.extend(normalise_human_gene(item)
                          for item in split_gene_text(row.get(column)))
        for token in (item for item in tokens if item):
            previous = mapping.setdefault(token, symbol)
            if previous != symbol:
                conflicts.add(token)
    return mapping, {
        'primary_column': primary_column,
        'symbol_column': symbol_column,
        'id_columns': id_columns,
        'alias_columns': alias_columns,
        'n_identifier_map_entries': len(mapping),
        'n_identifier_conflicts': len(conflicts),
        'identifier_conflicts': sorted(conflicts),
    }


def map_gene_values(values, identifier_map=None) -> list[str]:
    """Map values through a supplied offline map, then Human-normalise them."""
    identifier_map = identifier_map or {}
    mapped = []
    for value in values:
        token = normalise_human_gene(value)
        if token:
            mapped.append(identifier_map.get(token, token))
    return mapped


def prepare_deg_gene_table(frame: pd.DataFrame):
    """Attach ``_analysis_gene`` and return mapping provenance for a DEG table."""
    primary = _first_column(frame, GENE_COLUMNS)
    if primary is None:
        raise ValueError('DEG 结果缺少 gene/gene_symbol 列。')
    identifier_map, provenance = build_identifier_map(frame, primary)
    result = frame.copy()
    result['_analysis_gene'] = [
        identifier_map.get(normalise_human_gene(value), normalise_human_gene(value))
        for value in result[primary]
    ]
    provenance['n_rows'] = int(len(result))
    provenance['n_rows_without_gene'] = int((result['_analysis_gene'] == '').sum())
    return result, identifier_map, provenance


def _parse_human_geneset_lines(lines, source_label) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for raw_line in lines:
        line = raw_line.rstrip('\r\n')
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) >= 3:
            term, genes = parts[0].strip(), parts[2:]
        elif len(parts) == 2:
            term, genes = parts[0].strip(), re.split(r'[,;\s]+', parts[1].strip())
        else:
            continue
        members = list(dict.fromkeys(
            gene for gene in (normalise_human_gene(item) for item in genes) if gene
        ))
        if term and members:
            result[term] = list(dict.fromkeys(result.get(term, []) + members))
    if not result:
        raise ValueError(f'基因集没有可用通路：{source_label}')
    return result


def load_human_genesets(path) -> dict[str, list[str]]:
    """Load a GMT/Enrichr TXT file and canonicalise members as Human symbols."""
    with Path(path).open(encoding='utf-8') as handle:
        return _parse_human_geneset_lines(handle, path)


def load_human_genesets_text(text, source_label='custom Human GMT') -> dict[str, list[str]]:
    """Parse an inline Human GMT snapshot without network or ID guessing."""
    if len(str(text or '').encode('utf-8')) > 5 * 1024 * 1024:
        raise ValueError('自定义 GMT 不能超过 5 MB。')
    return _parse_human_geneset_lines(str(text or '').splitlines(), source_label)


def geneset_manifest(path, pathways: dict[str, list[str]], *, source_type='local_file') -> dict:
    """Return version/checksum evidence for the exact local gene-set resource."""
    source = Path(path)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    union = {gene for members in pathways.values() for gene in members}
    year = re.search(r'(?<!\d)(20\d{2})(?!\d)', source.stem)
    return {
        'file_name': source.name,
        'library_name': source.stem,
        'library_year': int(year.group(1)) if year else None,
        'source_type': source_type,
        'modified_utc': datetime.fromtimestamp(
            source.stat().st_mtime, tz=timezone.utc,
        ).isoformat(),
        'sha256': digest,
        'n_pathways': len(pathways),
        'n_unique_genes': len(union),
    }


def mapping_qc(raw_values, mapped_values, library_genes, *, background=None) -> dict:
    """Describe duplicate, annotation and background losses for one gene list."""
    raw_nonempty = [str(value).strip() for value in raw_values
                    if normalise_human_gene(value)]
    mapped_nonempty = [normalise_human_gene(value) for value in mapped_values
                       if normalise_human_gene(value)]
    unique = set(mapped_nonempty)
    library_genes = set(library_genes)
    mapped_to_library = unique.intersection(library_genes)
    unresolved = sorted(unique - library_genes)
    background_set = set(background or [])
    outside_background = sorted(unique - background_set) if background is not None else []
    mapped_in_background = (mapped_to_library.intersection(background_set)
                            if background is not None else mapped_to_library)
    denominator = len(unique)
    return {
        'n_input_rows': len(raw_nonempty),
        'n_unique_genes': denominator,
        'n_duplicates_collapsed': max(0, len(mapped_nonempty) - denominator),
        'n_mapped_to_library': len(mapped_to_library),
        'n_mapped_in_background': len(mapped_in_background),
        'mapping_rate': (len(mapped_to_library) / denominator) if denominator else 0.0,
        'n_unmapped': len(unresolved),
        'unmapped_genes': unresolved,
        'n_unresolved_ensembl': sum(gene.startswith('ENSG') for gene in unresolved),
        'n_outside_background': len(outside_background),
        'n_in_background': denominator - len(outside_background),
        'outside_background_genes': outside_background,
    }


def _bh_adjust(pvalues) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    if values.size == 0:
        return values
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return output


def run_ora_full(query_genes, background_genes, pathways, pvalue_cutoff=0.05,
                 min_size=1, max_size=None):
    """Run Human ORA against a real tested-gene universe and retain all terms.

    Every pathway with at least one member in the supplied background belongs to
    the BH family, including pathways with zero query hits (P=1).  This is the
    auditable full-table contract; plots may select the significant subset.
    """
    background_submitted = set(map_gene_values(background_genes))
    library_genes = {
        gene for members in pathways.values() for gene in map_gene_values(members)
    }
    # Only annotated members of the submitted tested-gene universe are eligible
    # for the hypergeometric table. Unmapped genes remain visible in QC but must
    # not inflate N or k.
    background = background_submitted.intersection(library_genes)
    query_submitted = set(map_gene_values(query_genes))
    query = query_submitted.intersection(background)
    if not background:
        raise ValueError('ORA 背景基因集为空。')
    if not query:
        raise ValueError('ORA 输入基因均不在背景基因集中。')
    if not 0 < float(pvalue_cutoff) <= 1:
        raise ValueError('显著性阈值必须在 (0, 1]。')
    min_size = int(min_size)
    max_size = len(background) if max_size is None else int(max_size)
    if min_size < 1 or max_size < min_size:
        raise ValueError('ORA 通路大小范围无效。')

    bg_size = len(background)
    query_size = len(query)
    rows = []
    n_background_supported = 0
    n_excluded_by_size = 0
    for term in sorted(pathways):
        category = set(map_gene_values(pathways[term])).intersection(background)
        category_size = len(category)
        if category_size == 0:
            continue
        n_background_supported += 1
        if category_size < min_size or category_size > max_size:
            n_excluded_by_size += 1
            continue
        hits = sorted(query.intersection(category))
        overlap = len(hits)
        pvalue = float(hypergeom.sf(overlap - 1, bg_size, category_size, query_size))
        # Exact 2x2 table with Haldane-Anscombe correction for zero cells.
        a = overlap
        b = query_size - overlap
        c = category_size - overlap
        d = bg_size - a - b - c
        odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        log_or_se = np.sqrt(
            1 / (a + 0.5) + 1 / (b + 0.5)
            + 1 / (c + 0.5) + 1 / (d + 0.5)
        )
        ci_low = np.exp(np.log(odds_ratio) - 1.96 * log_or_se)
        ci_high = np.exp(np.log(odds_ratio) + 1.96 * log_or_se)
        combined = -np.log(max(pvalue, 1e-300)) * np.log(max(odds_ratio, 1e-300))
        rows.append({
            'Gene_set': 'gs_ind',
            'Term': str(term),
            'Overlap': f'{overlap}/{category_size}',
            'P-value': pvalue,
            'Odds Ratio': float(odds_ratio),
            'Odds Ratio CI95 Low': float(ci_low),
            'Odds Ratio CI95 High': float(ci_high),
            'Combined Score': float(combined),
            'Genes': ';'.join(hits),
            'Count': overlap,
            'GeneRatio': overlap / float(query_size),
            'BackgroundRatio': category_size / float(bg_size),
            'FoldEnrichment': ((overlap / float(query_size))
                               / (category_size / float(bg_size))),
            'GeneSetSize': category_size,
            'GeneSetSizeMin': min_size,
            'GeneSetSizeMax': max_size,
            'SubmittedInputGeneCount': len(query_submitted),
            'InputGeneCount': query_size,
            'SubmittedBackgroundGeneCount': len(background_submitted),
            'BackgroundGeneCount': bg_size,
        })
    if not rows:
        if n_background_supported and n_excluded_by_size == n_background_supported:
            raise ValueError(
                f'当前背景支持的 {n_background_supported} 条通路均被大小范围 '
                f'{min_size}–{max_size} 排除。'
            )
        raise ValueError('所选数据库在当前背景基因集中没有可检验通路。')
    result = pd.DataFrame(rows)
    result['Adjusted P-value'] = _bh_adjust(result['P-value'].to_numpy())
    result['Significant'] = result['Adjusted P-value'] < float(pvalue_cutoff)
    result['logp'] = -np.log10(result['Adjusted P-value'].clip(lower=1e-300))
    result['logc'] = np.log(result['Odds Ratio'].clip(lower=1e-300))
    result['num'] = result['Count']
    result['fraction'] = result['BackgroundRatio']
    result.attrs['testing_family'] = {
        'n_pathways_in_library': len(pathways),
        'n_pathways_supported_by_background': n_background_supported,
        'n_pathways_excluded_by_size': n_excluded_by_size,
        'n_pathways_tested': len(result),
        'min_size': min_size,
        'max_size': max_size,
    }
    preferred = [
        'Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Significant',
        'Odds Ratio', 'Odds Ratio CI95 Low', 'Odds Ratio CI95 High',
        'FoldEnrichment', 'Combined Score', 'Genes', 'Count', 'GeneRatio',
        'BackgroundRatio', 'GeneSetSize', 'GeneSetSizeMin', 'GeneSetSizeMax',
        'SubmittedInputGeneCount', 'InputGeneCount',
        'SubmittedBackgroundGeneCount', 'BackgroundGeneCount', 'logp', 'logc',
        'num', 'fraction',
    ]
    output = result[preferred].sort_values(
        ['Adjusted P-value', 'P-value', 'Term'], kind='stable',
    ).reset_index(drop=True)
    output.attrs.update(result.attrs)
    return output


def _genes_from_value(value) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (list, tuple, set, np.ndarray)):
        values = value
    else:
        values = re.split(r'[;,/|]', str(value))
    return list(dict.fromkeys(
        normalise_human_gene(item) for item in values if normalise_human_gene(item)
    ))


def annotate_redundancy(result: pd.DataFrame, threshold=0.85) -> pd.DataFrame:
    """Annotate significant pathways with auditable hit-gene Jaccard clusters.

    This never removes statistical rows. Representatives are chosen in the
    existing table order, which should already be FDR-first.
    """
    threshold = float(threshold)
    if not 0 <= threshold <= 1:
        raise ValueError('冗余阈值必须在 [0, 1]。')
    output = result.copy()
    for column, default in (
        ('RedundancyCluster', ''), ('RepresentativeTerm', ''),
        ('IsRepresentative', False), ('JaccardToRepresentative', np.nan),
        ('RedundancyClusterSize', 0),
    ):
        output[column] = default
    term_column = _metric_column(output, ('Term', 'Description', 'pathway', 'term'))
    gene_column = _metric_column(
        output, ('Genes', 'genes', 'geneID', 'matched_genes', 'lead_genes'),
    )
    if term_column is None or gene_column is None or output.empty:
        return output
    if 'Significant' in output.columns:
        significant = output['Significant']
        if significant.dtype != bool:
            significant = significant.astype(str).str.lower().isin({'true', '1', 'yes'})
    else:
        significant = pd.Series(True, index=output.index)

    representatives = []
    assignments = {}
    for index, row in output.loc[significant].iterrows():
        genes = set(_genes_from_value(row[gene_column]))
        if not genes:
            continue
        best = None
        for cluster_index, representative in enumerate(representatives, start=1):
            union = genes | representative['genes']
            similarity = len(genes & representative['genes']) / len(union) if union else 0.0
            if similarity >= threshold and (best is None or similarity > best[1]):
                best = (cluster_index, similarity, representative)
        if best is None:
            representatives.append({
                'term': str(row[term_column]), 'genes': genes, 'index': index,
            })
            cluster_index = len(representatives)
            assignments[index] = (cluster_index, 1.0, str(row[term_column]), True)
        else:
            cluster_index, similarity, representative = best
            assignments[index] = (
                cluster_index, similarity, representative['term'], False,
            )

    sizes = {}
    for cluster_index, _, _, _ in assignments.values():
        sizes[cluster_index] = sizes.get(cluster_index, 0) + 1
    for index, (cluster_index, similarity, representative, is_representative) in assignments.items():
        output.at[index, 'RedundancyCluster'] = f'R{cluster_index:03d}'
        output.at[index, 'RepresentativeTerm'] = representative
        output.at[index, 'IsRepresentative'] = is_representative
        output.at[index, 'JaccardToRepresentative'] = float(similarity)
        output.at[index, 'RedundancyClusterSize'] = sizes[cluster_index]
    return output


def add_gsea_leading_edge_metrics(result: pd.DataFrame) -> pd.DataFrame:
    """Add explicit leading-edge counts/ratios to a full GSEA result table."""
    output = result.copy()
    lead_column = _metric_column(
        output, ('lead_genes', 'Lead_genes', 'leading_edge', 'core_enrichment'),
    )
    size_column = _metric_column(
        output, ('matched_size', 'geneset_size', 'setSize', 'size'),
    )
    if lead_column is None:
        output['LeadingEdgeGeneCount'] = 0
        output['LeadingEdgeRatio'] = np.nan
        return output
    counts = output[lead_column].map(lambda value: len(_genes_from_value(value)))
    output['LeadingEdgeGeneCount'] = counts.astype(int)
    sizes = (pd.to_numeric(output[size_column], errors='coerce')
             if size_column else pd.Series(np.nan, index=output.index))
    output['LeadingEdgeRatio'] = counts / sizes.replace(0, np.nan)
    return output


def gsea_leading_edge_table(result: pd.DataFrame, ranking: pd.DataFrame) -> pd.DataFrame:
    """Expand pathway-level leading edges to one traceable pathway-gene row."""
    columns = [
        'Term', 'Gene', 'RankPosition', 'RankingScore', 'NES', 'FDR', 'Significant',
    ]
    term_column = _metric_column(result, ('Term', 'Description', 'pathway', 'term'))
    lead_column = _metric_column(
        result, ('lead_genes', 'Lead_genes', 'leading_edge', 'core_enrichment'),
    )
    if term_column is None or lead_column is None or result.empty:
        return pd.DataFrame(columns=columns)
    rank_frame = ranking.reset_index(drop=True).copy()
    gene_column = _metric_column(rank_frame, ('gene_name', 'gene', 'Gene'))
    rank_column = _metric_column(rank_frame, ('rank', 'score', 'stat'))
    if gene_column is None or rank_column is None:
        raise ValueError('GSEA leading-edge 导出需要 gene_name/rank 排名表。')
    rank_lookup = {
        normalise_human_gene(gene): (position + 1, float(score))
        for position, (gene, score) in enumerate(
            zip(rank_frame[gene_column], rank_frame[rank_column])
        ) if normalise_human_gene(gene) and np.isfinite(float(score))
    }
    nes_column = _metric_column(result, ('nes', 'NES'))
    fdr_column = _metric_column(result, ('fdr', 'FDR', 'Adjusted P-value'))
    rows = []
    for _, pathway in result.iterrows():
        significant = pathway.get('Significant', False)
        if not isinstance(significant, (bool, np.bool_)):
            significant = str(significant).lower() in {'true', '1', 'yes'}
        for gene in _genes_from_value(pathway[lead_column]):
            position, score = rank_lookup.get(gene, (np.nan, np.nan))
            rows.append({
                'Term': str(pathway[term_column]),
                'Gene': gene,
                'RankPosition': position,
                'RankingScore': score,
                'NES': pathway[nes_column] if nes_column else np.nan,
                'FDR': pathway[fdr_column] if fdr_column else np.nan,
                'Significant': bool(significant),
            })
    return pd.DataFrame(rows, columns=columns)


def _metric_column(frame: pd.DataFrame, candidates) -> str | None:
    lookup = {str(column).lower(): column for column in frame.columns}
    return next((lookup[candidate.lower()] for candidate in candidates
                 if candidate.lower() in lookup), None)


def prepare_gsea_ranking(frame: pd.DataFrame, ranking_metric='auto', identifier_map=None):
    """Create a finite, unique Human pre-ranked table with explicit provenance."""
    primary = _first_column(frame, GENE_COLUMNS)
    if primary is None:
        raise ValueError('GSEA 需要 gene/gene_symbol 列。')
    identifier_map = identifier_map or build_identifier_map(frame, primary)[0]
    metric = str(ranking_metric or 'auto').strip()
    statistic_col = _metric_column(frame, STATISTIC_COLUMNS)
    log2fc_col = _metric_column(frame, LOG2FC_COLUMNS)
    pvalue_col = _metric_column(frame, PVALUE_COLUMNS)

    if metric == 'auto':
        source = statistic_col or log2fc_col
        metric_used = 'statistic' if statistic_col else 'log2FC'
    elif metric == 'statistic':
        source = statistic_col
        metric_used = 'statistic'
    elif metric == 'log2FC':
        source = log2fc_col
        metric_used = 'log2FC'
    elif metric == 'signed_log10_pvalue':
        source = None
        metric_used = metric
    else:
        raise ValueError(f'不支持的 GSEA 排序指标: {metric}')

    if metric_used == 'statistic' and not source:
        raise ValueError('DEG 结果没有可用的 Wald/t statistic 列。')
    if metric_used == 'log2FC' and not source:
        raise ValueError('DEG 结果没有 log2FC 列。')
    if metric_used == 'signed_log10_pvalue' and (not log2fc_col or not pvalue_col):
        raise ValueError('signed_log10_pvalue 需要 log2FC 和 pvalue/padj 列。')

    genes = [identifier_map.get(normalise_human_gene(value), normalise_human_gene(value))
             for value in frame[primary]]
    if metric_used == 'signed_log10_pvalue':
        fold_change = pd.to_numeric(frame[log2fc_col], errors='coerce').to_numpy(dtype=float)
        pvalues = pd.to_numeric(frame[pvalue_col], errors='coerce').to_numpy(dtype=float)
        scores = np.sign(fold_change) * -np.log10(np.clip(pvalues, 1e-300, 1.0))
    else:
        scores = pd.to_numeric(frame[source], errors='coerce').to_numpy(dtype=float)

    ranking = pd.DataFrame({'gene_name': genes, 'rank': scores})
    valid = (ranking['gene_name'] != '') & np.isfinite(ranking['rank'])
    n_invalid = int((~valid).sum())
    ranking = ranking.loc[valid].copy()
    n_before_dedup = len(ranking)
    ranking['_abs_rank'] = ranking['rank'].abs()
    ranking = (ranking.sort_values(['_abs_rank', 'rank'], ascending=[False, False], kind='stable')
               .drop_duplicates('gene_name', keep='first'))
    n_duplicates = n_before_dedup - len(ranking)
    ranking = ranking.drop(columns='_abs_rank').sort_values('rank', ascending=False, kind='stable')
    n_tied = int(ranking['rank'].duplicated(keep=False).sum())
    if ranking.empty:
        raise ValueError('GSEA 排名清洗后没有有效基因。')
    return ranking.reset_index(drop=True), {
        'ranking_metric_requested': metric,
        'ranking_metric_used': metric_used,
        'ranking_source_column': source or f'{log2fc_col} + {pvalue_col}',
        'n_ranked_input_rows': int(len(frame)),
        'n_ranked_genes': int(len(ranking)),
        'n_invalid_ranking_rows_removed': n_invalid,
        'n_duplicate_genes_collapsed': int(n_duplicates),
        'n_genes_with_tied_scores': n_tied,
    }


def run_gsea_prerank(ranking, pathways, *, permutation_num=1000, seed=112,
                     weight=1.0, min_size=15, max_size=500):
    """Run OmicVerse's deterministic NumPy GSEA without importing ``ov.bulk``.

    OmicVerse's public lazy ``bulk`` import also imports Scanpy modules that are
    unrelated to preranked GSEA and can fail in read-only/packaged runtimes when
    Numba attempts to create a cache locator.  The installed pure-NumPy backend
    is a self-contained module, so load that exact implementation directly and
    record the package version in the caller's QC manifest.
    """
    from importlib.metadata import PackageNotFoundError, distribution

    try:
        package_root = Path(distribution('omicverse').locate_file('omicverse'))
    except PackageNotFoundError as exc:
        raise RuntimeError('GSEA 需要已安装的 omicverse NumPy backend。') from exc
    module_path = package_root / 'bulk' / '_gsea_numpy.py'
    if not module_path.is_file():
        raise RuntimeError(f'未找到 OmicVerse NumPy GSEA backend：{module_path}')
    spec = importlib.util.spec_from_file_location('_platform_omicverse_gsea_numpy', module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError('无法加载 OmicVerse NumPy GSEA backend。')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.prerank(
        ranking, pathways,
        permutation_num=int(permutation_num), seed=int(seed), weight=float(weight),
        min_size=int(min_size), max_size=int(max_size), progress=False,
    )
