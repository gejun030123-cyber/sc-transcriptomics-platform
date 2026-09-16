"""Bulk raw-count upload contract and DESeq2 handoff regression tests."""

import io
import os

import anndata as ad
import numpy as np
import pandas as pd


COUNTS = b"""gene_id,Ctrl_1,Ctrl_2,Treat_1,Treat_2
GAPDH,100,110,220,230
TP53,20,22,44,46
EGFR,0,1,7,8
"""
COUNTS_WITH_VERSIONED_ENSEMBL_IDS = b"""gene_id,Ctrl_1,Ctrl_2,Treat_1,Treat_2
ENSG000001.5,100,110,220,230
ENSG000002,20,22,44,46
"""
METADATA = b"""sample_id,condition,batch
Treat_2,Treat,B2
Ctrl_1,Ctrl,B1
Treat_1,Treat,B2
Ctrl_2,Ctrl,B1
"""
FEATURECOUNTS = b"""Geneid,Chr,Start,End,Strand,Length,gene_name,Ctrl_1,Ctrl_2,Treat_1,Treat_2
ENSG000001,1,100,200,+,101,GAPDH,100,110,220,230
ENSG000002,1,300,400,-,101,TP53,20,22,44,46
"""
FEATURECOUNTS_WITH_PREAMBLE = b"""# Program:featureCounts v2.0.6; Command:\"featureCounts\"
""" + FEATURECOUNTS
FEATURECOUNTS_WITH_BAM_HEADERS = b"""Geneid,Chr,Start,End,Strand,Length,mapping/Ctr_B_1.Aligned.sortedByCoord.out.bam,mapping/Ctr_B_2.Aligned.sortedByCoord.out.bam,mapping/Ctr_B_3.Aligned.sortedByCoord.out.bam,mapping/Ctr_En_1.Aligned.sortedByCoord.out.bam,mapping/Ctr_En_2.Aligned.sortedByCoord.out.bam,mapping/Ctr_En_3.Aligned.sortedByCoord.out.bam
ENSG000001,1,100,200,+,101,100,110,120,220,230,240
ENSG000002,1,300,400,-,101,20,22,23,44,46,48
"""
FEATURECOUNTS_BAM_METADATA = b"""sample_id,condition
Ctr_B_1,Ctrl_B
Ctr_B_2,Ctrl_B
Ctr_B_3,Ctrl_B
Ctr_En_1,Ctrl_En
Ctr_En_2,Ctrl_En
Ctr_En_3,Ctrl_En
"""
GENE_ANNOTATION_GTF = b"""##description: test annotation
chr1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id \"ENSG000001\"; gene_name \"GAPDH\";
chr1\ttest\tgene\t300\t400\t.\t-\t.\tgene_id \"ENSG000002\"; gene_name \"TP53\";
"""
GENE_ANNOTATION_TABLE = b"""gene_id,gene_name
ENSG000001,GAPDH
ENSG000002,TP53
"""


def test_bulk_count_import_creates_deseq2_ready_h5ad(test_project):
    from app import create_app
    from config import Config
    from modules.design_preflight import preflight_blockers
    from modules.io_utils import read_expression_matrix, resolve_expression_measurement

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(COUNTS), 'raw_counts.csv'),
                'sample_metadata': (io.BytesIO(METADATA), 'sample_sheet.csv'),
            },
            content_type='multipart/form-data',
        )

    assert response.status_code == 201
    payload = response.get_json()
    assert payload['ok'] is True
    assert payload['n_samples'] == 4
    assert payload['n_genes'] == 3
    assert payload['conditions'] == {'Ctrl': 2, 'Treat': 2}
    assert payload['next_url'].endswith('/analyze/bulk_qc')

    output_path = os.path.join(Config.uploads_dir(test_project), payload['output_file'])
    imported = ad.read_h5ad(output_path)
    assert imported.obs_names.tolist() == ['Ctrl_1', 'Ctrl_2', 'Treat_1', 'Treat_2']
    assert imported.obs['condition'].tolist() == ['Ctrl', 'Ctrl', 'Treat', 'Treat']
    assert imported.obs['batch'].tolist() == ['B1', 'B1', 'B2', 'B2']
    assert imported.var_names.tolist() == ['GAPDH', 'TP53', 'EGFR']
    assert np.array_equal(imported.X, np.asarray([
        [100, 20, 0], [110, 22, 1], [220, 44, 7], [230, 46, 8],
    ]))
    assert imported.uns['bulk_import']['contract'] == 'bulk_raw_counts_with_sample_metadata_v1'
    assert len(imported.uns['bulk_import']['count_matrix_sha256']) == 64

    loaded = read_expression_matrix(output_path)
    measurement, details = resolve_expression_measurement(loaded, output_path)
    assert measurement == 'raw_counts'
    assert details['source'] == 'input_metadata'
    assert preflight_blockers(
        loaded, 'bulk_deg', {'method': 'deseq2'}, output_path,
    ) == []


def test_bulk_count_import_rejects_noninteger_expression(test_project):
    from app import create_app

    app = create_app()
    app.config['TESTING'] = True
    bad_counts = COUNTS.replace(b'220', b'22.5', 1)
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(bad_counts), 'tpm_like.csv'),
                'sample_metadata': (io.BytesIO(METADATA), 'sample_sheet.csv'),
            },
            content_type='multipart/form-data',
        )

    assert response.status_code == 400
    assert '非整数表达值' in response.get_json()['error']


def test_bulk_count_import_accepts_featurecounts_annotation_columns(test_project):
    from app import create_app
    from config import Config

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(FEATURECOUNTS), 'featurecounts.csv'),
                'sample_metadata': (io.BytesIO(METADATA), 'sample_sheet.csv'),
            },
            content_type='multipart/form-data',
        )

    assert response.status_code == 201
    path = os.path.join(Config.uploads_dir(test_project), response.get_json()['output_file'])
    imported = ad.read_h5ad(path)
    assert imported.shape == (4, 2)
    assert imported.var['gene_name'].tolist() == ['GAPDH', 'TP53']
    assert imported.var['chromosome'].tolist() == ['1', '1']
    assert list(imported.uns['bulk_import']['ignored_annotation_columns']) == [
        'Chr', 'Start', 'End', 'Strand', 'Length', 'gene_name',
    ]
    assert imported.uns['bulk_import']['chromosome_source'] == 'count_matrix'


def test_bulk_count_import_skips_featurecounts_provenance_preamble(test_project):
    from app import create_app

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(FEATURECOUNTS_WITH_PREAMBLE), 'featurecounts.txt'),
                'sample_metadata': (io.BytesIO(METADATA), 'sample_sheet.csv'),
            },
            content_type='multipart/form-data',
        )

    assert response.status_code == 201
    assert response.get_json()['n_samples'] == 4


def test_bulk_count_import_matches_featurecounts_bam_headers_to_sample_sheet(test_project):
    from app import create_app
    from config import Config

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(FEATURECOUNTS_WITH_BAM_HEADERS), 'featurecounts.csv'),
                'sample_metadata': (io.BytesIO(FEATURECOUNTS_BAM_METADATA), 'sample_sheet.csv'),
            },
            content_type='multipart/form-data',
        )

    assert response.status_code == 201
    path = os.path.join(Config.uploads_dir(test_project), response.get_json()['output_file'])
    imported = ad.read_h5ad(path)
    assert imported.obs_names.tolist() == [
        'Ctr_B_1', 'Ctr_B_2', 'Ctr_B_3', 'Ctr_En_1', 'Ctr_En_2', 'Ctr_En_3',
    ]
    assert imported.uns['bulk_import']['sample_id_source'] == 'featurecounts_alignment_filename_alias'


def test_bulk_count_import_maps_featurecounts_gene_ids_from_gtf(test_project):
    """A standard featureCounts Geneid-only table can retain readable symbols."""
    from app import create_app
    from config import Config

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(FEATURECOUNTS_WITH_BAM_HEADERS), 'featurecounts.csv'),
                'sample_metadata': (io.BytesIO(FEATURECOUNTS_BAM_METADATA), 'sample_sheet.csv'),
                'gene_annotation': (io.BytesIO(GENE_ANNOTATION_GTF), 'reference.gtf'),
            },
            content_type='multipart/form-data',
        )

    assert response.status_code == 201
    payload = response.get_json()
    assert payload['n_gene_names'] == 2
    assert payload['gene_name_source'] == 'annotation_file'
    path = os.path.join(Config.uploads_dir(test_project), payload['output_file'])
    imported = ad.read_h5ad(path)
    assert imported.var['gene_id'].tolist() == ['ENSG000001', 'ENSG000002']
    assert imported.var['gene_name'].tolist() == ['GAPDH', 'TP53']
    assert imported.uns['bulk_import']['gene_annotation']['source_format'] == 'gtf_gff'
    assert imported.uns['bulk_import']['gene_annotation']['n_matched_total'] == 2


def test_bulk_count_import_maps_versioned_gene_ids_from_annotation_table(tmp_path):
    """Only a trailing Ensembl version is tolerated when matching local maps."""
    from modules.bulk_import import build_bulk_counts_adata

    counts_path = tmp_path / 'counts.csv'
    metadata_path = tmp_path / 'samples.csv'
    annotation_path = tmp_path / 'genes.csv'
    counts_path.write_bytes(COUNTS_WITH_VERSIONED_ENSEMBL_IDS)
    metadata_path.write_bytes(METADATA)
    annotation_path.write_bytes(GENE_ANNOTATION_TABLE)

    imported = build_bulk_counts_adata(
        str(counts_path), str(metadata_path), str(annotation_path),
    )

    assert imported.var_names.tolist()[0] == 'ENSG000001.5'
    assert imported.var['gene_name'].tolist() == ['GAPDH', 'TP53']
    annotation = imported.uns['bulk_import']['gene_annotation']
    assert annotation['source_format'] == 'tabular'
    assert annotation['n_matched_exact'] == 1
    assert annotation['n_matched_versionless'] == 1


def test_imported_counts_run_deseq2_without_manual_groupby(test_project):
    """The h5ad emitted by the upload endpoint reaches the actual DESeq2 runner."""
    from app import create_app
    from config import Config
    from modules.bulk_deg import BulkDEGAnalysis

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(COUNTS), 'raw_counts.csv'),
                'sample_metadata': (io.BytesIO(METADATA), 'sample_sheet.csv'),
            },
            content_type='multipart/form-data',
        )
    assert response.status_code == 201
    input_path = os.path.join(Config.uploads_dir(test_project), response.get_json()['output_file'])

    result = BulkDEGAnalysis(
        Config.project_dir(test_project),
        {
            'method': 'deseq2',
            # groupby is intentionally omitted: imported condition is canonical.
            'group1': 'Treat',
            'group2': 'Ctrl',
            'fc_threshold': 1.5,
            'pval_threshold': 0.05,
        },
        lambda *_: None,
    ).run(input_path)

    assert result['summary']['method'] == 'deseq2'
    assert result['summary']['groupby'] == 'condition'
    result_table = pd.read_csv(os.path.join(
        Config.results_dir(test_project), 'bulk_deg_results.csv'))
    assert {'gene', 'gene_id', 'gene_name'}.issubset(result_table.columns)
    assert set(result_table['gene_id']) == {'GAPDH', 'TP53', 'EGFR'}


def test_bulk_count_import_requires_replicates(test_project):
    from app import create_app

    app = create_app()
    app.config['TESTING'] = True
    bad_metadata = b"""sample_id,condition
Ctrl_1,Ctrl
Ctrl_2,Ctrl
Treat_1,Treat
Treat_2,Other
"""
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(COUNTS), 'raw_counts.csv'),
                'sample_metadata': (io.BytesIO(bad_metadata), 'sample_sheet.csv'),
            },
            content_type='multipart/form-data',
        )

    assert response.status_code == 400
    assert '至少需要 2 个生物学重复' in response.get_json()['error']


def test_bulk_deg_uses_preserved_raw_counts_after_normalization(tmp_path):
    from modules.bulk_deg import _preserved_raw_counts

    raw = np.asarray([[10, 2], [12, 3], [20, 5], [24, 6]], dtype=int)
    normalized = ad.AnnData(
        X=np.log2(raw + 1),
        obs=pd.DataFrame({'condition': ['Ctrl', 'Ctrl', 'Treat', 'Treat']}, index=['C1', 'C2', 'T1', 'T2']),
        var=pd.DataFrame(index=['G1', 'G2']),
        layers={'raw': raw},
    )
    layer, name = _preserved_raw_counts(normalized)
    assert name == 'raw'
    assert np.array_equal(layer, raw)


def test_pre_deseq2_gene_filter_provenance_tracks_normalization_filter():
    """The DEG summary must state whether low-expression genes were removed."""
    from modules.bulk_deg import _pre_deseq2_gene_filter_provenance

    adata = ad.AnnData(X=np.ones((3, 2)), var=pd.DataFrame(index=['G1', 'G2']))

    missing = _pre_deseq2_gene_filter_provenance(adata)
    assert missing['source'] == 'none'
    assert missing['applied_before_deg'] is False
    assert missing['n_genes_tested'] == 2
    assert '未检测到低表达基因过滤记录' in missing['message']

    adata.uns['normalization'] = {'gene_expression_filter': {
        'applied': True,
        'n_genes_before': 100,
        'n_genes_after': 40,
        'n_genes_removed': 60,
        'message': '标准化阶段已过滤低表达基因。',
    }}
    filtered = _pre_deseq2_gene_filter_provenance(adata)
    assert filtered['source'] == 'bulk_normalize'
    assert filtered['applied_before_deg'] is True
    assert filtered['n_genes_removed_upstream'] == 60
    assert filtered['message'] == '标准化阶段已过滤低表达基因。'


def test_pre_deseq2_gene_filter_provenance_is_json_safe_after_h5ad_roundtrip(tmp_path):
    """uns values return as numpy scalars; the manifest writer uses strict json."""
    import json
    from modules.bulk_deg import _pre_deseq2_gene_filter_provenance

    adata = ad.AnnData(X=np.ones((3, 2)), var=pd.DataFrame(index=['G1', 'G2']))
    adata.uns['normalization'] = {'gene_expression_filter': {
        'applied': True,
        'n_genes_before': 100,
        'n_genes_after': 40,
        'n_genes_removed': 60,
        'minimum_expression_value': 1.0,
    }}
    path = tmp_path / 'normalized.h5ad'
    adata.write_h5ad(path)
    restored = ad.read_h5ad(path)

    record = _pre_deseq2_gene_filter_provenance(restored)

    assert record['applied_before_deg'] is True
    assert record['n_genes_removed_upstream'] == 60
    assert isinstance(record['n_genes_removed_upstream'], int)
    # Strict serialisation is what the analysis manifest performs.
    json.dumps({'gene_filtering': record}, ensure_ascii=False)


def test_deseq2_uses_raw_layer_after_bulk_normalization(test_project):
    """QC/normalization output remains a valid DESeq2 input via its raw layer."""
    from app import create_app
    from config import Config
    from modules.bulk_deg import BulkDEGAnalysis
    from modules.bulk_normalize import BulkNormalizeAnalysis

    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-bulk-counts',
            data={
                'count_matrix': (io.BytesIO(COUNTS), 'raw_counts.csv'),
                'sample_metadata': (io.BytesIO(METADATA), 'sample_sheet.csv'),
            },
            content_type='multipart/form-data',
        )
    input_path = os.path.join(Config.uploads_dir(test_project), response.get_json()['output_file'])

    normalized = BulkNormalizeAnalysis(
        Config.project_dir(test_project),
        {'method': 'deseq2', 'min_expr_samples': 0},
        lambda *_: None,
    ).run(input_path)
    normalized_path = normalized['output_adata']
    assert 'raw' in ad.read_h5ad(normalized_path).layers

    messages = []
    result = BulkDEGAnalysis(
        Config.project_dir(test_project),
        {'method': 'deseq2', 'group1': 'Treat', 'group2': 'Ctrl'},
        lambda _, message: messages.append(message),
    ).run(normalized_path)
    assert result['summary']['method'] == 'deseq2'
    assert any('保留的原始计数 layer' in message for message in messages)
    # min_expr_samples=0 disabled the upstream filter: the DEG summary must say
    # so rather than letting an unfiltered gene set pass silently.
    gene_filtering = result['summary']['gene_filtering']
    assert gene_filtering['source'] == 'bulk_normalize'
    assert gene_filtering['applied_before_deg'] is False
    assert gene_filtering['n_genes_tested'] > 0
    assert any('未过滤低表达基因' in message for message in messages)
