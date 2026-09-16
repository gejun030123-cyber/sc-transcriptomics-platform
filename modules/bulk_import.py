"""Validated import contract for Bulk RNA-seq raw-count tables.

The generic project uploader deliberately accepts many scientific file types.
That is useful for exploratory work, but a DESeq2-ready Bulk input needs a
stricter contract: a genes-by-samples *integer* count matrix plus explicit
sample metadata. This module keeps that contract in one place and writes a
portable ``.h5ad`` whose ``obs`` table supplies the downstream design.
"""

from __future__ import annotations

import hashlib
import gzip
import os
import re
from collections import Counter

import anndata as ad
import numpy as np
import pandas as pd


TABULAR_EXTENSIONS = frozenset({'.csv', '.tsv', '.txt', '.xlsx', '.xls'})
GENE_ANNOTATION_EXTENSIONS = frozenset({
    *TABULAR_EXTENSIONS, '.gtf', '.gff', '.gff3',
    '.gtf.gz', '.gff.gz', '.gff3.gz',
})

# featureCounts, HTSeq-adjacent exports, and vendor count tables commonly
# retain these feature annotations before the actual sample columns. They are
# metadata, never sample counts. Restrict this to explicit header names: an
# arbitrary non-numeric sample column must still fail loudly rather than be
# silently dropped.
_ANNOTATION_COLUMN_NAMES = frozenset({
    'gene_name', 'genename', 'gene_symbol', 'genesymbol', 'symbol',
    'description', 'gene_description', 'chromosome', 'chrom', 'chr', 'seqname', 'seq_name', 'sequence_name',
    'start', 'end', 'strand', 'length', 'gene_length', 'gene_size',
    'exonic_gene_sizes', 'gc', 'gc_content', 'biotype', 'feature',
    'feature_type', 'annotation',
})
_GENE_NAME_COLUMN_NAMES = frozenset({
    'gene_name', 'genename', 'gene_symbol', 'genesymbol', 'symbol',
})
_GENE_ID_COLUMN_NAMES = frozenset({
    'gene_id', 'geneid', 'ensembl_gene_id', 'ensembl_id', 'feature_id',
})
_CHROMOSOME_COLUMN_NAMES = frozenset({
    'chromosome', 'chrom', 'chr', 'seqname', 'seq_name', 'sequence_name',
})


def _normalise_column_name(column):
    """Normalise a header only for matching known annotation aliases."""
    import re
    return re.sub(r'[^a-z0-9]+', '_', str(column).strip().lower()).strip('_')


def _featurecounts_sample_alias(header):
    """Return the sample ID encoded by a common featureCounts BAM header.

    STAR output is commonly passed to featureCounts as
    ``mapping/S1.Aligned.sortedByCoord.out.bam`` while the experimental design
    correctly calls the sample ``S1``. This is an import-only alias candidate;
    it is used only when every resulting ID is unique and exactly matches the
    supplied sample sheet.
    """
    import re

    name = re.split(r'[\\/]', str(header).strip())[-1]
    name = re.sub(r'\.(?:bam|sam|cram)$', '', name, flags=re.IGNORECASE)
    name = re.sub(
        r'(?:[._-]aligned)?[._-]sortedbycoord[._-]out$', '', name,
        flags=re.IGNORECASE,
    )
    name = re.sub(r'[._-]aligned[._-]out$', '', name, flags=re.IGNORECASE)
    return name.strip()


def tabular_extension_allowed(filename):
    """Return whether a browser-uploaded Bulk table has a supported suffix."""
    return os.path.splitext(str(filename or ''))[1].lower() in TABULAR_EXTENSIONS


def gene_annotation_extension_allowed(filename):
    """Return whether an optional gene-ID mapping file has a safe format.

    The mapping must come from the exact reference used for counting.  This
    function intentionally accepts only local tabular/GTF/GFF files and never
    falls back to a web lookup, which could silently use the wrong species or
    genome build.
    """
    lower = str(filename or '').strip().lower()
    return any(lower.endswith(extension) for extension in GENE_ANNOTATION_EXTENSIONS)


def file_sha256(path, *, chunk_size=1024 * 1024):
    """Calculate a streaming checksum without loading a sequencing table at once."""
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _read_table(path, label):
    """Read one supported tabular file while preserving its header row."""
    suffix = os.path.splitext(path)[1].lower()
    try:
        if suffix in {'.xlsx', '.xls'}:
            return pd.read_excel(path)
        # The Python sniffer copes with CSV, TSV, and common plain-text exports.
        # featureCounts prepends provenance such as ``# Program:featureCounts``;
        # it is not a header or a sample column and must be skipped before
        # identifying the actual Geneid/Chr/... table header.
        return pd.read_csv(path, sep=None, engine='python', comment='#')
    except (OSError, ValueError, pd.errors.ParserError, UnicodeDecodeError) as exc:
        raise ValueError(f'无法读取{label}；请提供 CSV、TSV、TXT 或 Excel 文件。') from exc


def _valid_gene_name(value):
    """Return a cleaned gene symbol/name, or an empty string for missing data."""
    if value is None:
        return ''
    try:
        if bool(pd.isna(value)):
            return ''
    except (TypeError, ValueError):
        pass
    value = str(value).strip()
    return '' if value.lower() in {'', 'nan', 'none', 'null'} else value


def _gene_identifier_without_version(value):
    """Remove only an unambiguous trailing Ensembl-style version component."""
    return re.sub(r'\.(\d+)$', '', str(value).strip())


def _build_gene_name_lookup(rows, *, source_label, id_column, name_column):
    """Build exact and version-tolerant local ID -> symbol lookup tables."""
    exact = {}
    conflicts = {}
    for raw_id, raw_name in rows:
        gene_id = str(raw_id).strip()
        gene_name = _valid_gene_name(raw_name)
        if not gene_id or not gene_name:
            continue
        existing = exact.get(gene_id)
        if existing and existing != gene_name:
            conflicts.setdefault(gene_id, {existing}).add(gene_name)
            continue
        exact[gene_id] = gene_name
    if conflicts:
        examples = '、'.join(sorted(conflicts)[:5])
        raise ValueError(
            f'{source_label}中同一 gene_id 对应多个 gene_name：{examples}。'
            '请使用与 counts 相同参考版本的注释文件。'
        )
    if not exact:
        raise ValueError(f'{source_label}未读取到有效的 gene_id 与 gene_name 映射。')

    versionless_candidates = {}
    for gene_id, gene_name in exact.items():
        versionless_candidates.setdefault(_gene_identifier_without_version(gene_id), set()).add(gene_name)
    versionless = {
        gene_id: next(iter(names)) for gene_id, names in versionless_candidates.items()
        if len(names) == 1
    }
    return exact, versionless, {
        'gene_id_column': str(id_column),
        'gene_name_column': str(name_column),
        'n_annotation_rows': int(len(exact)),
        'n_ambiguous_versionless_ids': int(
            sum(len(names) > 1 for names in versionless_candidates.values())
        ),
    }


def _parse_gtf_attributes(raw_attributes):
    """Read standard GTF ``key \"value\";`` and GFF3 ``key=value`` fields."""
    values = {}
    for piece in str(raw_attributes).strip().strip(';').split(';'):
        piece = piece.strip()
        if not piece:
            continue
        if '=' in piece:
            key, value = piece.split('=', 1)
        else:
            match = re.match(r'([^\s]+)\s+(.+)$', piece)
            if not match:
                continue
            key, value = match.groups()
        values[key.strip().lower()] = value.strip().strip('"')
    return values


def _read_gtf_gene_mapping(path):
    """Parse gene_id/gene_name from a local GTF or GFF annotation file."""
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    rows = []
    try:
        with opener(path, 'rt', encoding='utf-8', errors='replace') as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) < 9:
                    continue
                # Prefer gene features so repeated exons/transcripts cannot
                # introduce a conflicting parent mapping.
                if fields[2].strip().lower() != 'gene':
                    continue
                attrs = _parse_gtf_attributes(fields[8])
                gene_id = attrs.get('gene_id') or attrs.get('id', '')
                if gene_id.lower().startswith('gene:'):
                    gene_id = gene_id.split(':', 1)[1]
                gene_name = (
                    attrs.get('gene_name') or attrs.get('gene_symbol')
                    or attrs.get('name') or attrs.get('gene') or ''
                )
                rows.append((gene_id, gene_name))
    except OSError as exc:
        raise ValueError('无法读取基因注释 GTF/GFF 文件。') from exc
    return _build_gene_name_lookup(
        rows, source_label='基因注释 GTF/GFF 文件',
        id_column='gene_id/ID', name_column='gene_name/Name',
    )


def _read_tabular_gene_mapping(path):
    """Parse an explicit two-column (or richer) local annotation table."""
    frame = _clean_columns(_read_table(path, '基因注释表'), '基因注释表')
    normalised = {_normalise_column_name(column): column for column in frame.columns}
    id_column = next(
        (normalised[name] for name in _GENE_ID_COLUMN_NAMES if name in normalised),
        None,
    )
    name_column = next(
        (normalised[name] for name in _GENE_NAME_COLUMN_NAMES if name in normalised),
        None,
    )
    if not id_column or not name_column:
        raise ValueError(
            '基因注释表必须包含 gene_id（或 GeneID/ensembl_gene_id）和 '
            'gene_name（或 gene_symbol/GeneSymbol）两列。'
        )
    return _build_gene_name_lookup(
        zip(frame[id_column], frame[name_column]), source_label='基因注释表',
        id_column=id_column, name_column=name_column,
    )


def _read_gene_annotation_mapping(path):
    """Read a supplied local mapping and record its reproducible provenance."""
    lower = str(path).lower()
    is_gtf_or_gff = lower.endswith(('.gtf', '.gff', '.gff3', '.gtf.gz', '.gff.gz', '.gff3.gz'))
    if is_gtf_or_gff:
        exact, versionless, details = _read_gtf_gene_mapping(path)
        source_format = 'gtf_gff'
    else:
        exact, versionless, details = _read_tabular_gene_mapping(path)
        source_format = 'tabular'
    return exact, versionless, {
        **details,
        'source_format': source_format,
        'filename': os.path.basename(path),
        'sha256': file_sha256(path),
    }


def _clean_columns(frame, label):
    columns = [str(column).strip() for column in frame.columns]
    if not columns or any(not column for column in columns):
        raise ValueError(f'{label}包含空列名。')
    duplicates = sorted({column for column, count in Counter(columns).items() if count > 1})
    if duplicates:
        raise ValueError(f'{label}列名重复：' + '、'.join(duplicates[:5]))
    result = frame.copy()
    result.columns = columns
    return result


def _clean_required_text(values, label):
    result = values.astype('string').str.strip()
    if result.isna().any() or result.eq('').any():
        raise ValueError(f'{label}不能包含空值。')
    return result.astype(str)


def _as_raw_integer_counts(frame):
    """Validate and return a finite, non-negative integer-valued count matrix."""
    numeric = frame.apply(pd.to_numeric, errors='coerce')
    invalid = numeric.isna().to_numpy()
    if invalid.any():
        row, column = np.argwhere(invalid)[0]
        raise ValueError(
            '原始 counts 矩阵必须全部为非负整数；'
            f'第 {int(row) + 2} 行、样本列 “{frame.columns[int(column)]}” 不是有效整数。'
        )
    values = numeric.to_numpy(dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError('原始 counts 矩阵不能包含 NaN、inf 或 -inf。')
    if (values < 0).any():
        raise ValueError('原始 counts 矩阵不能包含负数。')
    if not np.isclose(values, np.round(values), rtol=0.0, atol=1e-8).all():
        raise ValueError('检测到非整数表达值；TPM、FPKM、CPM 或 log 表达值不能用于 DESeq2。')
    if values.size and float(values.max()) >= float(np.iinfo(np.int64).max):
        raise ValueError('counts 数值超出可安全导入的整数范围。')
    return np.rint(values).astype(np.int64)


def build_bulk_counts_adata(count_matrix_path, metadata_path, gene_annotation_path=None):
    """Build a DESeq2-ready ``AnnData`` from raw counts and a sample sheet.

    Required count-table shape is ``gene_id, sample_1, sample_2, ...``.
    The sample sheet must contain exactly one row per count-table sample and
    the required columns ``sample_id`` and ``condition``. Extra metadata
    columns, such as ``batch`` or ``sex``, are retained in ``obs``. An
    optional local mapping table/GTF can supply gene symbols when a standard
    featureCounts output contains only ``Geneid``.
    """
    count_table = _clean_columns(_read_table(count_matrix_path, 'counts 矩阵'), 'counts 矩阵')
    if count_table.shape[0] == 0 or count_table.shape[1] < 3:
        raise ValueError('counts 矩阵至少需要 1 列 gene_id 和 2 个样本列。')

    gene_column = count_table.columns[0]
    gene_ids = _clean_required_text(count_table.iloc[:, 0], f'基因 ID 列 “{gene_column}”')
    duplicate_genes = gene_ids[gene_ids.duplicated()].unique().tolist()
    if duplicate_genes:
        raise ValueError('基因 ID 不能重复：' + '、'.join(duplicate_genes[:5]))

    annotation_columns = [
        column for column in count_table.columns[1:]
        if _normalise_column_name(column) in _ANNOTATION_COLUMN_NAMES
    ]
    count_columns = [
        column for column in count_table.columns[1:]
        if column not in annotation_columns
    ]
    if len(count_columns) < 2:
        raise ValueError('counts 矩阵至少需要两个样本 count 列。')
    source_sample_headers = [str(column).strip() for column in count_columns]
    if any(not sample_id for sample_id in source_sample_headers):
        raise ValueError('counts 矩阵包含空样本列名。')
    duplicate_samples = sorted({item for item, count in Counter(source_sample_headers).items() if count > 1})
    if duplicate_samples:
        raise ValueError('counts 矩阵样本列名重复：' + '、'.join(duplicate_samples[:5]))
    counts = _as_raw_integer_counts(count_table[count_columns])

    metadata = _clean_columns(_read_table(metadata_path, '样本信息表'), '样本信息表')
    required = {'sample_id', 'condition'}
    missing_columns = sorted(required - set(metadata.columns))
    if missing_columns:
        raise ValueError('样本信息表缺少必填列：' + '、'.join(missing_columns))
    if metadata.empty:
        raise ValueError('样本信息表不能为空。')
    metadata = metadata.copy()
    metadata['sample_id'] = _clean_required_text(metadata['sample_id'], 'sample_id')
    metadata['condition'] = _clean_required_text(metadata['condition'], 'condition')
    duplicate_metadata = metadata.loc[metadata['sample_id'].duplicated(), 'sample_id'].unique().tolist()
    if duplicate_metadata:
        raise ValueError('样本信息表的 sample_id 不能重复：' + '、'.join(duplicate_metadata[:5]))

    sample_ids = list(source_sample_headers)
    sample_id_source = 'count_matrix_header'
    matrix_samples = set(sample_ids)
    metadata_samples = set(metadata['sample_id'])
    missing_metadata = [sample_id for sample_id in sample_ids if sample_id not in metadata_samples]
    extra_metadata = [sample_id for sample_id in metadata['sample_id'] if sample_id not in matrix_samples]
    if missing_metadata or extra_metadata:
        aliases = [_featurecounts_sample_alias(header) for header in source_sample_headers]
        alias_is_unambiguous = (
            all(aliases)
            and len(set(aliases)) == len(aliases)
            and set(aliases) == metadata_samples
        )
        if alias_is_unambiguous:
            sample_ids = aliases
            sample_id_source = 'featurecounts_alignment_filename_alias'
            matrix_samples = set(sample_ids)
            missing_metadata = []
            extra_metadata = []
    if missing_metadata or extra_metadata:
        parts = []
        if missing_metadata:
            parts.append('矩阵中缺少样本信息：' + '、'.join(missing_metadata[:5]))
        if extra_metadata:
            parts.append('样本信息表中不存在于矩阵的样本：' + '、'.join(extra_metadata[:5]))
        raise ValueError('；'.join(parts))

    # A direct DESeq2 contrast needs replicate-based dispersion estimation.
    # The platform's DEG runner has the same lower bound, so reject a bad
    # design during upload instead of allowing a deceptively successful import.
    condition_counts = metadata['condition'].value_counts()
    singleton_conditions = condition_counts[condition_counts < 2].index.tolist()
    if len(condition_counts) < 2:
        raise ValueError('样本信息表至少需要两个 condition 组。')
    if singleton_conditions:
        raise ValueError('每个 condition 至少需要 2 个生物学重复；样本不足的组：' + '、'.join(singleton_conditions[:5]))

    obs = metadata.set_index('sample_id', drop=False).loc[sample_ids].copy()
    # AnnData/HDF5 writes homogeneous string fields reliably. Numeric
    # covariates stay numeric; arbitrary object columns become explicit text.
    for column in obs.columns:
        if pd.api.types.is_object_dtype(obs[column]) or pd.api.types.is_string_dtype(obs[column]):
            obs[column] = obs[column].fillna('').astype(str)
    obs.index = pd.Index(sample_ids, name='sample_id')
    var = pd.DataFrame(index=pd.Index(gene_ids.tolist(), name=gene_column))
    # Keep this explicit even when the index is Geneid.  It makes every DEG
    # result auditable after a human-readable symbol is added as ``gene``.
    var['gene_id'] = gene_ids.tolist()
    gene_name_column = next(
        (column for column in annotation_columns
         if _normalise_column_name(column) in _GENE_NAME_COLUMN_NAMES),
        None,
    )
    chromosome_column = next(
        (column for column in annotation_columns
         if _normalise_column_name(column) in _CHROMOSOME_COLUMN_NAMES),
        None,
    )
    gene_names = [''] * len(gene_ids)
    gene_name_source = 'none'
    if gene_name_column:
        gene_names = [_valid_gene_name(value) for value in count_table[gene_name_column]]
        if any(gene_names):
            gene_name_source = 'count_matrix'

    annotation_info = None
    if gene_annotation_path:
        exact_lookup, versionless_lookup, annotation_info = _read_gene_annotation_mapping(
            gene_annotation_path)
        mapped_exact = 0
        mapped_versionless = 0
        for index, gene_id in enumerate(gene_ids):
            # A symbol embedded in the exact count table remains authoritative.
            # The optional annotation fills only blank values, avoiding a
            # silent mismatch between counting and display annotations.
            if gene_names[index]:
                continue
            gene_id = str(gene_id)
            gene_name = exact_lookup.get(gene_id)
            if gene_name:
                mapped_exact += 1
            else:
                gene_name = versionless_lookup.get(_gene_identifier_without_version(gene_id), '')
                if gene_name:
                    mapped_versionless += 1
            if gene_name:
                gene_names[index] = gene_name
        external_mapped = mapped_exact + mapped_versionless
        if external_mapped == 0:
            raise ValueError(
                '基因注释文件与 counts 的 gene_id 没有匹配项。'
                '请确认它来自相同物种和相同参考基因组/注释版本。'
            )
        annotation_info.update({
            'n_matched_exact': int(mapped_exact),
            'n_matched_versionless': int(mapped_versionless),
            'n_matched_total': int(external_mapped),
        })
        if external_mapped:
            gene_name_source = (
                'count_matrix_plus_annotation_file'
                if gene_name_column else 'annotation_file'
            )

    n_gene_names = int(sum(bool(value) for value in gene_names))
    if n_gene_names:
        var['gene_name'] = gene_names
    if chromosome_column:
        chromosome_values = [
            _valid_gene_name(value) for value in count_table[chromosome_column]
        ]
        if any(chromosome_values):
            # Use a canonical field name even if the source calls it ``Chr``;
            # this survives h5ad handoff and lets Bulk QC detect chrM/MT.
            var['chromosome'] = chromosome_values
    adata = ad.AnnData(X=counts.T, obs=obs, var=var)
    adata.uns['input_measurement'] = 'raw_counts'
    adata.uns['bulk_import'] = {
        'contract': 'bulk_raw_counts_with_sample_metadata_v1',
        'measurement': 'raw_counts',
        'gene_id_column': str(gene_column),
        'ignored_annotation_columns': [str(column) for column in annotation_columns],
        'sample_id_source': sample_id_source,
        'sample_metadata_required_columns': ['sample_id', 'condition'],
        'n_samples': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'gene_name_source': gene_name_source,
        'n_gene_names': n_gene_names,
        'chromosome_source': 'count_matrix' if chromosome_column else 'none',
        'condition_sizes': {str(key): int(value) for key, value in condition_counts.sort_index().items()},
        'count_matrix_sha256': file_sha256(count_matrix_path),
        'sample_metadata_sha256': file_sha256(metadata_path),
    }
    if annotation_info:
        adata.uns['bulk_import']['gene_annotation'] = annotation_info
    return adata
