"""P0 regression tests for auditable Human Bulk enrichment statistics."""

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


def test_human_identifier_mapping_uses_supplied_symbol_and_alias_columns():
    from modules.enrichment_statistics import prepare_deg_gene_table

    frame = pd.DataFrame({
        'gene': ['ENSG00000141510.18', 'MYC'],
        'gene_symbol': ['TP53', 'MYC'],
        'aliases': ['P53;BCC7', 'C-MYC'],
        'log2FC': [2.0, -1.5],
    })
    prepared, identifier_map, provenance = prepare_deg_gene_table(frame)

    assert prepared['_analysis_gene'].tolist() == ['TP53', 'MYC']
    assert identifier_map['ENSG00000141510'] == 'TP53'
    assert identifier_map['P53'] == 'TP53'
    assert identifier_map['C-MYC'] == 'MYC'
    assert provenance['symbol_column'] == 'gene_symbol'


def test_ora_uses_real_background_and_retains_zero_hit_terms():
    from modules.enrichment_statistics import run_ora_full

    result = run_ora_full(
        query_genes=['A', 'B'],
        background_genes=['A', 'B', 'C', 'D', 'E'],
        pathways={
            'query pathway': ['A', 'B'],
            'zero hit pathway': ['C', 'D'],
            'outside background': ['Z'],
        },
        pvalue_cutoff=0.4,
    )

    assert result['Term'].tolist() == ['query pathway', 'zero hit pathway']
    assert result['SubmittedBackgroundGeneCount'].eq(5).all()
    assert result['BackgroundGeneCount'].eq(4).all()
    assert result['InputGeneCount'].eq(2).all()
    assert result.loc[result['Term'] == 'zero hit pathway', 'P-value'].item() == 1.0
    assert result.loc[result['Term'] == 'zero hit pathway', 'Overlap'].item() == '0/2'
    # Both background-supported pathways belong to the BH family.
    assert result.loc[
        result['Term'] == 'query pathway', 'Adjusted P-value'
    ].item() == pytest.approx(1 / 3)
    assert result['Significant'].tolist() == [True, False]
    enriched = result.loc[result['Term'] == 'query pathway'].iloc[0]
    assert enriched['FoldEnrichment'] == pytest.approx(2.0)
    assert enriched['Odds Ratio CI95 Low'] < enriched['Odds Ratio']
    assert enriched['Odds Ratio CI95 High'] > enriched['Odds Ratio']


def test_ora_size_filter_defines_bh_family_before_query_hits():
    from modules.enrichment_statistics import run_ora_full

    result = run_ora_full(
        query_genes=['A', 'B'],
        background_genes=['A', 'B', 'C', 'D', 'E', 'F'],
        pathways={
            'too small but hit': ['A'],
            'eligible hit': ['A', 'B'],
            'eligible zero': ['C', 'D'],
            'too large': ['A', 'B', 'C', 'D', 'E'],
        },
        pvalue_cutoff=0.5,
        min_size=2,
        max_size=4,
    )

    assert set(result['Term']) == {'eligible hit', 'eligible zero'}
    assert result.attrs['testing_family'] == {
        'n_pathways_in_library': 4,
        'n_pathways_supported_by_background': 4,
        'n_pathways_excluded_by_size': 2,
        'n_pathways_tested': 2,
        'min_size': 2,
        'max_size': 4,
    }


def test_redundancy_annotation_keeps_all_rows_and_marks_representatives():
    from modules.enrichment_statistics import annotate_redundancy

    frame = pd.DataFrame({
        'Term': ['first', 'near duplicate', 'distinct', 'not significant'],
        'Genes': ['A;B;C', 'A;B;C;D', 'X;Y', 'A;B;C'],
        'Significant': [True, True, True, False],
    })
    annotated = annotate_redundancy(frame, threshold=0.7)

    assert len(annotated) == len(frame)
    assert annotated['RedundancyCluster'].tolist() == ['R001', 'R001', 'R002', '']
    assert annotated['IsRepresentative'].tolist() == [True, False, True, False]
    assert annotated.loc[1, 'RepresentativeTerm'] == 'first'
    assert annotated.loc[0, 'RedundancyClusterSize'] == 2


def test_inline_human_gmt_parser_and_manifest_are_reproducible(tmp_path):
    from modules.enrichment_statistics import (
        geneset_manifest, load_human_genesets_text,
    )

    text = '# Human custom v1\nPath A\tdescription\ttp53\tMYC\nPath B\tna\tEGFR\n'
    pathways = load_human_genesets_text(text)
    snapshot = tmp_path / 'Lab_Human_2026.txt'
    snapshot.write_text(text, encoding='utf-8')
    manifest = geneset_manifest(snapshot, pathways, source_type='custom_inline_human_gmt')

    assert pathways == {'Path A': ['TP53', 'MYC'], 'Path B': ['EGFR']}
    assert manifest['library_year'] == 2026
    assert manifest['source_type'] == 'custom_inline_human_gmt'
    assert len(manifest['sha256']) == 64


def test_gsea_ranking_prefers_statistic_and_collapses_duplicate_symbols():
    from modules.enrichment_statistics import prepare_deg_gene_table, prepare_gsea_ranking

    frame = pd.DataFrame({
        'gene': ['ENSG1', 'ENSG2', 'B', 'C'],
        'gene_symbol': ['A', 'A', 'B', 'C'],
        'stat': [2.0, -4.0, np.nan, 1.0],
        'log2FC': [8.0, 7.0, 6.0, 5.0],
    })
    _, identifier_map, _ = prepare_deg_gene_table(frame)
    ranking, qc = prepare_gsea_ranking(frame, 'auto', identifier_map)

    assert qc['ranking_metric_used'] == 'statistic'
    assert qc['n_duplicate_genes_collapsed'] == 1
    assert qc['n_invalid_ranking_rows_removed'] == 1
    assert dict(ranking.values.tolist()) == {'C': 1.0, 'A': -4.0}


def test_direct_omicverse_numpy_gsea_backend_runs_without_bulk_lazy_import():
    from modules.enrichment_statistics import run_gsea_prerank

    ranking = pd.DataFrame({
        'gene_name': ['A', 'B', 'C', 'D', 'E'],
        'rank': [3.0, 2.0, 0.1, -1.5, -2.5],
    })
    result = run_gsea_prerank(
        ranking,
        {'positive': ['A', 'B'], 'negative': ['D', 'E']},
        permutation_num=100, seed=7, min_size=2, max_size=4,
    )

    assert set(result.res2d.index) == {'positive', 'negative'}
    assert {'es', 'nes', 'pval', 'fdr', 'lead_genes'}.issubset(result.res2d.columns)


def _write_tiny_inputs(tmp_path):
    geneset = tmp_path / 'human_tiny.txt'
    geneset.write_text(
        'query pathway\t\tA\tB\n'
        'zero hit pathway\t\tC\tD\n',
        encoding='utf-8',
    )
    deg = tmp_path / 'bulk_deg_results.csv'
    pd.DataFrame({
        'gene': ['A', 'B', 'C', 'D', 'E'],
        'log2FC': [2.0, 1.5, -0.2, 0.1, 0.0],
        'stat': [4.0, 3.0, -1.0, 0.5, 0.0],
        'pvalue': [.001, .005, .4, .8, 1.0],
        'padj': [.01, .02, .5, .9, 1.0],
        'regulation': ['Up', 'Up', 'NS', 'NS', 'NS'],
    }).to_csv(deg, index=False)
    return geneset, deg


def _disable_enrichment_figures(monkeypatch):
    import modules.bulk_enrichment as enrichment

    monkeypatch.setattr(enrichment, '_enrichment_figure', lambda *args, **kwargs: None)
    monkeypatch.setattr(
        enrichment, '_render_enrichment_variants',
        lambda *args, **kwargs: ([], []),
    )


def test_bulk_ora_end_to_end_exports_full_table_and_mapping_qc(tmp_path, monkeypatch):
    import modules.bulk_enrichment as enrichment

    geneset, deg = _write_tiny_inputs(tmp_path)
    monkeypatch.setattr(enrichment, '_resolve_geneset_path', lambda *args: str(geneset))
    _disable_enrichment_figures(monkeypatch)
    analysis = enrichment.BulkEnrichmentAnalysis(
        project_dir=str(tmp_path),
        params={
            'method': 'ORA', 'database': 'GO_BP', 'organism': 'Human',
            'pvalue_cutoff': 0.4, 'input_source': str(deg),
            'input_comparison': 'Treat vs Ctrl', 'enrichment_plot_suite': '核心图',
            'ora_min_size': 1, 'ora_max_size': 4,
        },
        progress_callback=lambda *_: None,
    )

    output = analysis.run(str(deg))
    table = pd.read_csv(
        tmp_path / 'results' / 'enrichment_ora_go_bp_treat_vs_ctrl_results.csv'
    )
    qc = json.loads((
        tmp_path / 'results' / 'enrichment_ora_go_bp_treat_vs_ctrl_mapping_qc.json'
    ).read_text(encoding='utf-8'))

    assert len(table) == 2
    assert table['SubmittedBackgroundGeneCount'].eq(5).all()
    assert table['BackgroundGeneCount'].eq(4).all()
    assert set(table['Significant'].astype(str)) == {'True', 'False'}
    assert qc['background']['source'] == 'all_genes_in_deg_table'
    assert qc['analyses'][0]['n_tested_pathways'] == 2
    assert output['summary']['n_tested'] == 2
    assert output['summary']['n_significant'] == 1
    assert output['summary']['mapping_rates'] == {
        'background': pytest.approx(0.8),
        'inputs': {'All': pytest.approx(1.0)},
    }
    mapping_output = next(
        item for item in output['result_files']
        if item['file_path'].endswith('_mapping_qc.json')
    )
    assert mapping_output['category'] == 'qc'
    assert (tmp_path / 'results' / 'enrichment_ora_go_bp_treat_vs_ctrl_methods.txt').is_file()
    assert (tmp_path / 'results' / 'enrichment_ora_go_bp_treat_vs_ctrl_source_data.xlsx').is_file()
    redundancy = pd.read_csv(
        tmp_path / 'results' / 'enrichment_ora_go_bp_treat_vs_ctrl_redundancy_clusters.csv'
    )
    assert redundancy['RepresentativeTerm'].tolist() == ['query pathway']


def test_bulk_gsea_forwards_reproducible_parameters_and_exports_all_terms(
        tmp_path, monkeypatch):
    import modules.bulk_enrichment as enrichment

    geneset, deg = _write_tiny_inputs(tmp_path)
    monkeypatch.setattr(enrichment, '_resolve_geneset_path', lambda *args: str(geneset))
    _disable_enrichment_figures(monkeypatch)
    captured = {}

    def fake_gsea(*args, **kwargs):
        captured.update(kwargs)
        result = pd.DataFrame({
            'es': [0.7, -0.4], 'nes': [2.1, -1.2],
            'pval': [.001, .2], 'fdr': [.01, .3],
            'geneset_size': [2, 2], 'matched_size': [2, 2],
            'lead_genes': ['A;B', 'C;D'],
        }, index=pd.Index(['query pathway', 'zero hit pathway'], name='Term'))
        return SimpleNamespace(res2d=result, results={}, ranking=args[0])

    monkeypatch.setattr(enrichment, 'run_gsea_prerank', fake_gsea)
    analysis = enrichment.BulkEnrichmentAnalysis(
        project_dir=str(tmp_path),
        params={
            'method': 'GSEA', 'database': 'GO_BP', 'organism': 'Human',
            'pvalue_cutoff': 0.05, 'input_source': str(deg),
            'input_comparison': 'Treat vs Ctrl', 'ranking_metric': 'statistic',
            'permutation_num': 2000, 'gsea_seed': 77, 'gsea_weight': 1.5,
            'gsea_min_size': 2, 'gsea_max_size': 4,
        },
        progress_callback=lambda *_: None,
    )

    output = analysis.run(str(deg))
    table = pd.read_csv(
        tmp_path / 'results' / 'enrichment_gsea_go_bp_treat_vs_ctrl_results.csv'
    )
    leading_edge = pd.read_csv(
        tmp_path / 'results' / 'enrichment_gsea_go_bp_treat_vs_ctrl_leading_edge.csv'
    )

    assert captured['permutation_num'] == 2000
    assert captured['seed'] == 77
    assert captured['weight'] == 1.5
    assert captured['min_size'] == 2
    assert captured['max_size'] == 4
    assert len(table) == 2
    assert table['RankingMetric'].eq('statistic').all()
    assert table['Seed'].eq(77).all()
    assert table['Significant'].tolist() == [True, False]
    assert table['LeadingEdgeGeneCount'].tolist() == [2, 2]
    assert set(leading_edge['Gene']) == {'A', 'B', 'C', 'D'}
    assert leading_edge.loc[leading_edge['Gene'] == 'A', 'RankPosition'].item() == 1
    assert output['summary']['n_tested'] == 2
    assert output['summary']['n_ranked_genes'] == 5
    assert output['summary']['n_leading_edge_rows'] == 4


def test_bulk_ora_accepts_custom_human_gmt_and_snapshots_source(tmp_path, monkeypatch):
    import modules.bulk_enrichment as enrichment

    _, deg = _write_tiny_inputs(tmp_path)
    _disable_enrichment_figures(monkeypatch)
    custom = 'Custom A\tdescription\tA\tB\nCustom B\tdescription\tC\tD\n'
    analysis = enrichment.BulkEnrichmentAnalysis(
        project_dir=str(tmp_path),
        params={
            'method': 'ORA', 'database': 'Custom_GMT', 'organism': 'Human',
            'custom_geneset_name': 'Lab panel v1', 'custom_geneset_text': custom,
            'pvalue_cutoff': 0.4, 'input_source': str(deg),
            'input_comparison': 'Treat vs Ctrl', 'ora_min_size': 1,
            'ora_max_size': 4,
        },
        progress_callback=lambda *_: None,
    )

    output = analysis.run(str(deg))
    snapshot = (
        tmp_path / 'results'
        / 'enrichment_ora_custom_gmt_treat_vs_ctrl_geneset_snapshot.txt'
    )
    qc = json.loads((
        tmp_path / 'results'
        / 'enrichment_ora_custom_gmt_treat_vs_ctrl_mapping_qc.json'
    ).read_text(encoding='utf-8'))

    assert snapshot.read_text(encoding='utf-8') == custom
    assert qc['gene_set_library']['source_type'] == 'custom_inline_human_gmt'
    assert qc['custom_geneset_name'] == 'Lab panel v1'
    assert any(item['file_path'] == str(snapshot) for item in output['result_files'])


def test_bulk_batch_enrichment_is_task_scoped_and_keeps_database_fdr_separate(
        tmp_path, monkeypatch):
    """A batch must integrate only its explicit child tables, never history."""
    import modules.bulk_enrichment as enrichment

    geneset, deg = _write_tiny_inputs(tmp_path)
    _disable_enrichment_figures(monkeypatch)
    monkeypatch.setattr(enrichment, '_resolve_geneset_path', lambda *args: str(geneset))
    monkeypatch.setattr(enrichment, '_write_integrated_overview_figures', lambda *args, **kwargs: [])

    historical = tmp_path / 'results' / 'enrichment_ora_reactome_old_results.csv'
    historical.parent.mkdir(exist_ok=True)
    pd.DataFrame({
        'Comparison': ['Old'], 'Database': ['Reactome'], 'Method': ['ORA'],
        'Direction': ['All'], 'Term': ['must not be included'],
        'Adjusted P-value': [0.001],
    }).to_csv(historical, index=False)

    output = enrichment.BulkEnrichmentAnalysis(
        project_dir=str(tmp_path),
        params={
            'method': 'ORA', 'batch_databases': True,
            'databases': 'GO_BP,KEGG', '_analysis_id': 'batch-001',
            'organism': 'Human', 'pvalue_cutoff': 0.4,
            'input_source': str(deg), 'input_comparison': 'Treat vs Ctrl',
            'ora_min_size': 1, 'ora_max_size': 4,
        },
        progress_callback=lambda *_: None,
    ).run(str(deg))

    results_dir = tmp_path / 'results' / 'enrichment_batches' / 'batch_001'
    integrated = pd.read_csv(
        results_dir / 'enrichment_batch_batch_001_integrated_results.csv'
    )
    go_table = pd.read_csv(results_dir / 'enrichment_ora_go_bp_treat_vs_ctrl_results.csv')
    kegg_table = pd.read_csv(results_dir / 'enrichment_ora_kegg_treat_vs_ctrl_results.csv')

    assert output['summary']['n_databases_completed'] == 2
    assert output['summary']['n_databases_failed'] == 0
    assert set(integrated['Database']) == {'GO_BP', 'KEGG'}
    assert 'must not be included' not in set(integrated['Term'])
    assert not (tmp_path / 'results' / 'enrichment_integrated_results.csv').exists()
    assert integrated.loc[integrated['Database'] == 'GO_BP', 'Enrichment FDR'].tolist() == go_table['Adjusted P-value'].tolist()
    assert integrated.loc[integrated['Database'] == 'KEGG', 'Enrichment FDR'].tolist() == kegg_table['Adjusted P-value'].tolist()


def test_bulk_batch_enrichment_retains_partial_library_failures(tmp_path, monkeypatch):
    import modules.bulk_enrichment as enrichment

    geneset, deg = _write_tiny_inputs(tmp_path)
    _disable_enrichment_figures(monkeypatch)
    monkeypatch.setattr(enrichment, '_write_integrated_overview_figures', lambda *args, **kwargs: [])

    def resolve_geneset(name, organism):
        if name == 'KEGG_2021_Human':
            raise FileNotFoundError('KEGG snapshot is not installed')
        return str(geneset)

    monkeypatch.setattr(enrichment, '_resolve_geneset_path', resolve_geneset)
    output = enrichment.BulkEnrichmentAnalysis(
        project_dir=str(tmp_path),
        params={
            'method': 'ORA', 'batch_databases': True,
            'databases': 'GO_BP,KEGG', '_analysis_id': 'partial-001',
            'organism': 'Human', 'pvalue_cutoff': 0.4,
            'input_source': str(deg), 'input_comparison': 'Treat vs Ctrl',
            'ora_min_size': 1, 'ora_max_size': 4,
        },
        progress_callback=lambda *_: None,
    ).run(str(deg))

    rows = {row['database']: row for row in output['summary']['database_results']}
    assert output['summary']['n_databases_completed'] == 1
    assert output['summary']['n_databases_failed'] == 1
    assert rows['GO_BP']['status'] == 'completed'
    assert rows['KEGG']['status'] == 'failed'
    assert 'snapshot is not installed' in rows['KEGG']['error']


def test_batch_database_form_selection_preserves_all_checked_values():
    from werkzeug.datastructures import MultiDict
    from modules.schemas import PARAM_SCHEMAS, parse_form_params

    params = parse_form_params(PARAM_SCHEMAS['bulk_enrichment'], MultiDict([
        ('batch_databases', 'on'), ('databases', 'GO_BP'), ('databases', 'KEGG'),
    ]))

    assert params['batch_databases'] is True
    assert params['databases'] == 'GO_BP,KEGG'
    assert 'database' not in params


def test_enrichment_schema_is_human_only_and_gsea_controls_are_conditional():
    from modules.schemas import PARAM_SCHEMAS

    schema = {item['key']: item for item in PARAM_SCHEMAS['bulk_enrichment']}
    assert schema['organism']['options'] == ['Human']
    assert 'Custom_GMT' in schema['database']['options']
    assert schema['batch_databases']['default'] is True
    assert schema['databases']['type'] == 'static_multiselect'
    assert schema['custom_geneset_text']['show_if'] == {'database': 'Custom_GMT'}
    assert schema['ora_min_size']['show_if'] == {'method': 'ORA'}
    assert schema['ranking_metric']['show_if'] == {'method': 'GSEA'}
    assert schema['custom_genes']['show_if'] == {'method': 'ORA'}
    assert schema['split_direction']['show_if'] == {'method': 'ORA'}
    assert schema['go_priority_enabled']['show_if'] == {
        'batch_databases': True, 'method': 'ORA',
    }
    assert schema['focus_terms']['show_if'] == {'go_priority_enabled': True}
    assert schema['go_priority_top_n']['max'] == 12
    assert schema['go_priority_allocation_mode']['show_if'] == {
        'go_priority_enabled': True,
    }


def test_bulk_go_priority_config_requires_all_three_go_ontologies():
    import pytest
    from modules.bulk_enrichment import BulkEnrichmentAnalysis, _go_priority_config_from_params

    config = _go_priority_config_from_params(
        {
            'go_priority_enabled': True,
            'method': 'ORA',
            'focus_terms': 'lipid catabolic process',
            'go_priority_top_n': 6,
        },
        databases=['GO_BP', 'GO_CC', 'GO_MF'], top_n=8,
    )

    assert config['focus_terms'] == ['lipid catabolic process']
    assert config['allocation_mode'] == 'balanced'
    assert config['top_n'] == 6
    with pytest.raises(ValueError, match='GO_BP、GO_CC、GO_MF'):
        _go_priority_config_from_params(
            {'go_priority_enabled': True, 'method': 'ORA'},
            databases=['GO_BP', 'GO_CC'], top_n=8,
        )
    with pytest.raises(ValueError, match='批量执行多个数据库'):
        BulkEnrichmentAnalysis(
            project_dir='.', params={'go_priority_enabled': True, 'batch_databases': False},
            progress_callback=lambda *_args: None,
        ).run('')
