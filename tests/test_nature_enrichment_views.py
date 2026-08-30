"""Regression tests for the extended Nature pathway-enrichment views."""

from pathlib import Path

import numpy as np
import pandas as pd


def _ora_fixture():
    return pd.DataFrame({
        'Term': [
            'Extracellular structure organization',
            'Cell migration',
            'Inflammatory response',
            'Mitochondrial organization',
            'Lipid transport',
        ],
        'Adjusted P-value': [1e-8, 2e-6, .002, .01, .03],
        'Overlap': ['8/120', '6/120', '5/120', '4/120', '3/120'],
        'Genes': [
            'COL1A1;COL1A2;SPP1;DCN;LUM;MMP2;MMP9;VCAN',
            'SPP1;VCAN;ITGB1;CXCR4;MMP2;ACTB',
            'STAT1;IRF1;CXCR4;SPP1;NFKB1',
            'COX5A;ATP5F1;NDUFA1;NDUFS2',
            'APOA1;APOB;LPL',
        ],
    })


def test_extended_ora_views_validate_at_single_column(tmp_path):
    from figure_engine import NatureFigureDirector, export_figure

    director = NatureFigureDirector()
    for plot_type in ('barplot', 'chord', 'cnetplot', 'emapplot'):
        spec = director.create_spec(plot_type, width='single', formats=('svg', 'pdf', 'png'))
        figure = director.render(spec, _ora_fixture())
        paths, report = export_figure(figure, tmp_path / plot_type, spec)
        assert report.ready, (plot_type, report.to_dict())
        assert set(paths) == {'svg', 'pdf', 'png'}
        assert '<text' in Path(paths['svg']).read_text(encoding='utf-8')
        import matplotlib.pyplot as plt
        plt.close(figure)


def test_multidatabase_overview_uses_facets_and_shared_encodings(tmp_path):
    from figure_engine import NatureFigureDirector, export_figure

    rows = []
    databases = {
        'GO_BP': ['matrix organization (GO:0001)', 'cell migration (GO:0002)'],
        'GO_CC': ['extracellular region (GO:0003)', 'membrane (GO:0004)'],
        'GO_MF': ['collagen binding (GO:0005)', 'receptor activity (GO:0006)'],
        'KEGG': ['ECM interaction (KEGG:hsa0007)', 'focal adhesion (KEGG:hsa0008)'],
    }
    for database, terms in databases.items():
        for index, term in enumerate(terms):
            rows.append({
                'Database': database, 'Method': 'ORA', 'Direction': 'All',
                'Comparison': 'demo', 'Term': term,
                'Adjusted P-value': 10 ** (-4 + index),
                'GeneRatio': .02 + .01 * index, 'num': 10 + index * 8,
                'Genes': f'{database}_{index}_A;{database}_{index}_B',
            })
    director = NatureFigureDirector()
    spec = director.create_spec(
        'enrichment_overview', width='double', top_n=4,
        formats=('svg', 'pdf', 'png'), title='demo · ORA overview',
    )
    figure = director.render(spec, pd.DataFrame(rows))
    assert getattr(figure, '_nature_panel_grid')['panels'] == 4
    assert figure.legends
    assert {text.get_text() for text in figure.legends[0].get_texts()} == {'10', '14', '18'}
    paths, report = export_figure(figure, tmp_path / 'overview', spec)
    assert report.ready, report.to_dict()
    assert report.metrics['panel_grid']['ncols'] == 2
    assert '<text' in Path(paths['svg']).read_text(encoding='utf-8')
    import matplotlib.pyplot as plt
    plt.close(figure)


def test_go_focus_triptych_keeps_an_explicit_empty_ontology_panel():
    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame({
        'Database': ['GO_BP', 'GO_MF'], 'Method': ['ORA', 'ORA'],
        'Direction': ['Up', 'Up'], 'Comparison': ['demo', 'demo'],
        'Term': ['lipid transport', 'cholesterol binding'],
        'Adjusted P-value': [.003, .02], 'Overlap': ['5/100', '3/100'],
        'Genes': ['APOE;APOB', 'APOA1;APOC3'],
    })
    spec = NatureFigureDirector().create_spec(
        'enrichment_overview', width='double', height_mm=180, top_n=6,
        formats=('png',), database_scope=('GO_BP', 'GO_CC', 'GO_MF'),
        extra={'facet_layout': 'one_column', 'include_empty_databases': True},
    )
    figure = NatureFigureDirector().render(spec, frame)
    assert getattr(figure, '_nature_panel_grid') == {'nrows': 3, 'ncols': 1, 'panels': 3}
    assert any('No FDR-significant' in text.get_text() for axis in figure.axes for text in axis.texts)
    import matplotlib.pyplot as plt
    plt.close(figure)


def test_go_focus_triptych_keeps_overlapping_selected_terms():
    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame({
        'Database': ['GO_BP', 'GO_BP'], 'Method': ['ORA', 'ORA'],
        'Direction': ['Up', 'Up'], 'Comparison': ['demo', 'demo'],
        'Term': ['inflammatory response', 'cytokine-mediated signalling'],
        'Adjusted P-value': [.003, .02], 'Overlap': ['5/100', '4/100'],
        # A high Jaccard similarity would normally reduce this to one term.
        'Genes': ['A;B;C;D;E', 'A;B;C;D;F'],
    })
    spec = NatureFigureDirector().create_spec(
        'enrichment_overview', width='double', height_mm=140, top_n=6,
        formats=('png',), database_scope=('GO_BP', 'GO_CC', 'GO_MF'),
        extra={
            'facet_layout': 'one_column', 'include_empty_databases': True,
            'disable_redundancy_compression': True,
        },
    )
    figure = NatureFigureDirector().render(spec, frame)
    labels = [label.get_text() for label in figure.axes[0].get_yticklabels()]
    assert labels == ['inflammatory response', 'cytokine-mediated signalling']
    import matplotlib.pyplot as plt
    plt.close(figure)


def test_go_focus_triptych_renders_truthful_empty_panels_when_no_term_is_significant():
    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame(columns=[
        'Database', 'Method', 'Direction', 'Comparison', 'Term',
        'Adjusted P-value', 'Overlap', 'Genes',
    ])
    spec = NatureFigureDirector().create_spec(
        'enrichment_overview', width='double', height_mm=180, top_n=6,
        formats=('png',), database_scope=('GO_BP', 'GO_CC', 'GO_MF'),
        extra={'facet_layout': 'one_column', 'include_empty_databases': True},
    )
    figure = NatureFigureDirector().render(spec, frame)
    assert getattr(figure, '_nature_panel_grid') == {'nrows': 3, 'ncols': 1, 'panels': 3}
    assert sum(
        'No FDR-significant' in text.get_text()
        for axis in figure.axes for text in axis.texts
    ) == 3
    import matplotlib.pyplot as plt
    plt.close(figure)


def test_gsea_running_uses_real_rank_and_hit_contract(tmp_path):
    from figure_engine import NatureFigureDirector, export_figure

    ranking = np.linspace(2.4, -2.1, 360)
    hits = np.array([4, 8, 15, 23, 31, 46, 59, 75, 91])
    spec = NatureFigureDirector().create_spec(
        'gsea_running', width='single', formats=('svg', 'pdf', 'png'), running_term_n=1,
    )
    figure = NatureFigureDirector().render(spec, {
        'curves': [{
            'term': 'Interferon gamma response',
            'ranking': ranking,
            'hit_indices': hits,
            'NES': 2.31,
            'FDR': .004,
        }],
    })
    paths, report = export_figure(figure, tmp_path / 'gsea_running', spec)
    assert report.ready, report.to_dict()
    line = figure.axes[0].lines[0]
    assert np.nanmax(np.abs(line.get_ydata())) <= 1.05
    assert len(figure.axes[1].collections) == 1
    assert len(figure.axes[1].collections[0].get_segments()) == len(hits)
    import matplotlib.pyplot as plt
    plt.close(figure)


def test_enrichment_suite_keeps_dotplot_and_adds_fixed_views():
    from modules.bulk_enrichment import _enrichment_plot_suite

    assert _enrichment_plot_suite({'enrichment_plot_suite': '完整 Nature 套图'}, 'ORA') == [
        'enrichment_dotplot', 'enrichment_barplot', 'enrichment_chord',
        'enrichment_cnetplot', 'enrichment_emapplot',
    ]
    assert _enrichment_plot_suite({'enrichment_plot_suite': '核心图'}, 'ORA') == [
        'enrichment_dotplot', 'enrichment_barplot',
    ]
    assert _enrichment_plot_suite({'enrichment_plot_suite': '核心图'}, 'GSEA') == [
        'gsea', 'enrichment_barplot', 'gsea_running',
    ]


def test_gene_ratio_and_term_redundancy_are_explicit(tmp_path):
    from figure_engine.templates.common import compress_redundant_terms, strip_term_id
    from modules.bulk_enrichment import _with_gene_ratio
    from figure_engine import NatureFigureDirector, export_figure

    frame = pd.DataFrame({
        'Term': ['A (GO:0000001)', 'B (GO:0000002)'],
        'Count': [8, 7],
        'Genes': ['A;B;C;D;E;F;G;H', 'A;B;C;D;E;F;G'],
    })
    corrected = _with_gene_ratio(frame, 40)
    assert corrected['GeneRatio'].tolist() == [0.2, 0.175]
    assert corrected['InputGeneCount'].tolist() == [40, 40]
    assert strip_term_id(frame.loc[0, 'Term']) == 'A'
    selected, removed = compress_redundant_terms(frame, 'Genes', threshold=.85, maximum=2)
    assert len(selected) == 1
    assert removed == 1

    invalid = pd.DataFrame({'Term': ['invalid'], 'FDR': [.01], 'GeneRatio': [1.2], 'Count': [12]})
    director = NatureFigureDirector()
    figure = director.render(director.create_spec('enrichment', top_n=1), invalid)
    _, report = export_figure(figure, tmp_path / 'invalid_ratio', director.create_spec('enrichment', top_n=1))
    assert report.ready is False
    assert any(issue.code == 'semantic_error' for issue in report.issues)
    import matplotlib.pyplot as plt
    plt.close(figure)


def test_gsea_top_n_balances_both_nes_directions_and_strips_ids():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    frame = pd.DataFrame({
        'Term': [
            'Positive A (GO:0000001)', 'Positive B (GO:0000002)',
            'Positive C (GO:0000003)', 'Negative A (GO:0000004)',
            'Negative B (GO:0000005)', 'Negative C (GO:0000006)',
        ],
        'NES': [3.2, 2.8, 2.5, -3.0, -2.6, -2.2],
        'FDR': [.001, .002, .003, .001, .002, .003],
        'setSize': [30, 25, 20, 28, 23, 18],
    })
    spec = NatureFigureDirector().create_spec('gsea', top_n=4)
    figure = NatureFigureDirector().render(spec, frame)
    labels = [text.get_text() for text in figure.axes[0].get_yticklabels()]
    assert len(labels) == 4
    assert any('Positive' in label for label in labels)
    assert any('Negative' in label for label in labels)
    assert not any('GO:' in label for label in labels)
    plt.close(figure)


def test_enrichment_barplot_exposes_score_and_gene_count_at_each_bar_end():
    import matplotlib.pyplot as plt

    from figure_engine import NatureFigureDirector

    director = NatureFigureDirector()
    figure = director.render(director.create_spec('barplot', top_n=3), _ora_fixture())
    labels = [text.get_text() for text in figure.axes[0].texts]
    assert any('8.00' in label and 'n=8' in label for label in labels)
    assert figure.axes[0].get_xlim()[1] > 8.0 * 1.4
    plt.close(figure)


def test_target_pathways_and_gene_labels_apply_to_network_views(tmp_path):
    import matplotlib.pyplot as plt
    from figure_engine import NatureFigureDirector, export_figure

    frame = _ora_fixture().copy()
    frame.loc[4, 'Term'] = 'Lipid transport (KEGG:00001)'
    director = NatureFigureDirector()
    spec = director.create_spec(
        'cnetplot', width='single', pathway_selection='selected',
        pathway_terms=('KEGG:00001',), gene_label_strategy='all',
        max_gene_labels=20, formats=('svg', 'pdf', 'png'),
    )
    figure = director.render(spec, frame)
    labels = [text.get_text() for text in figure.axes[0].texts]
    assert any('Lipid transport' in label for label in labels)
    assert any('APOA1' in label for label in labels)
    assert not any('P1  Extracellular' in label for label in labels)
    assert any('用户选择' in warning for warning in figure._nature_semantic_warnings)
    paths, report = export_figure(figure, tmp_path / 'target_cnet', spec)
    assert report.ready, report.to_dict()
    assert '<text' in Path(paths['svg']).read_text(encoding='utf-8')
    plt.close(figure)


def test_director_maps_pathway_and_gene_selection_controls():
    from figure_engine import NatureFigureDirector

    spec = NatureFigureDirector().spec_from_params('gsea', {
        'pathway_selection': 'selected_plus_top',
        'target_pathways': 'term A\nGO:0000002',
        'gene_label_strategy': 'selected',
        'target_genes': 'COL1A1, SPP1',
        'max_gene_labels': 11,
    })
    assert spec.pathway_selection == 'selected_plus_top'
    assert spec.pathway_terms == ('term A', 'GO:0000002')
    assert spec.gene_label_strategy == 'selected'
    assert spec.gene_labels == ('COL1A1', 'SPP1')
    assert spec.max_gene_labels == 11

    overview = NatureFigureDirector().spec_from_params('enrichment_overview', {
        'database_scope': 'GO_BP,KEGG', 'top_n': 6,
    })
    assert overview.width == 'double'
    assert overview.database_scope == ('GO_BP', 'KEGG')


def test_selected_gene_is_prioritised_before_network_node_cap():
    import matplotlib.pyplot as plt
    from figure_engine import NatureFigureDirector

    spec = NatureFigureDirector().create_spec(
        'cnetplot', width='single', max_genes=4, max_gene_labels=4,
        gene_label_strategy='selected', gene_labels=('MMP9',),
    )
    figure = NatureFigureDirector().render(spec, _ora_fixture())
    labels = [text.get_text() for text in figure.axes[0].texts]
    assert 'MMP9' in labels
    plt.close(figure)
