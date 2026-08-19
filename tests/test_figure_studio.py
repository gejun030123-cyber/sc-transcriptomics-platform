"""Coverage for non-destructive Figure Studio rendering and project scoping."""

import base64
import io
import os

import pytest

from app import create_app


def _write_png(path, color=(245, 247, 250)):
    from PIL import Image
    Image.new('RGB', (160, 100), color).save(path, format='PNG')


def _make_volcano_source(project_id):
    from config import Config
    from models import AnalysisTask, ResultFile

    results_dir = Config.results_dir(project_id)
    os.makedirs(results_dir, exist_ok=True)
    png_path = os.path.join(results_dir, 'bulk_deg_volcano.png')
    _write_png(png_path)
    csv_path = os.path.join(results_dir, 'bulk_deg_results.csv')
    with open(csv_path, 'w', encoding='utf-8') as handle:
        handle.write(
            'gene,log2FC,padj\n'
            'TP53,2.6,0.001\nVEGFA,1.9,0.004\nIL6,0.3,0.8\nMYC,-2.4,0.002\n'
        )
    task = AnalysisTask(project_id=project_id, module_name='bulk_deg', status='completed')
    task.save()
    return ResultFile.create(
        task_id=task.id, project_id=project_id, file_type='png', category='plot',
        label='Volcano plot', file_path=png_path,
    )


def test_volcano_version_reclassifies_without_recomputing_deg(test_project):
    from modules.figure_studio import render_preview_data_uri, resolve_source, save_figure_version

    source_file = _make_volcano_source(test_project)
    source = resolve_source(test_project, 'result_file', source_file.id)
    assert source['edit_mode'] == 'bulk_volcano'

    style = {
        'pvalue_threshold': 0.01, 'fc_threshold': 1.5,
        'up_color': '#b64342', 'down_color': '#0f4d92', 'ns_color': '#c7cdd6',
        'label_genes': 'TP53, MYC', 'title': 'Adjusted volcano',
    }
    preview = render_preview_data_uri(source, style)
    assert preview['edit_mode'] == 'data_redraw'
    assert preview['data_uri'].startswith('data:image/png;base64,')

    version = save_figure_version(test_project, source, style, 'threshold-adjusted volcano')
    assert version.edit_mode == 'bulk_volcano'
    assert os.path.isfile(version.png_path)
    assert os.path.isfile(version.svg_path)
    assert 'TP53' in open(version.svg_path, encoding='utf-8').read()


def test_saved_volcano_version_keeps_data_backed_editing(test_project):
    from modules.figure_studio import list_sources, resolve_source, save_figure_version

    source_file = _make_volcano_source(test_project)
    original = resolve_source(test_project, 'result_file', source_file.id)
    version = save_figure_version(test_project, original, {'label_genes': 'TP53'}, 'editable volcano')
    reopened = resolve_source(test_project, 'version', version.id)

    assert reopened['edit_mode'] == 'bulk_volcano'
    source_item = next(item for item in list_sources(test_project) if item['kind'] == 'version')
    assert source_item['edit_mode'] == 'bulk_volcano'


def test_heatmap_version_is_redrawn_from_output_adata(test_project):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from config import Config
    from models import AnalysisTask, ResultFile
    from modules.figure_studio import render_preview_data_uri, resolve_source, save_figure_version

    intermediate_dir = Config.intermediate_dir(test_project)
    results_dir = Config.results_dir(test_project)
    os.makedirs(intermediate_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    output_path = os.path.join(intermediate_dir, 'bulk_heatmap_output.h5ad')
    ad.AnnData(
        X=np.asarray([[1, 3, 8, 2], [2, 4, 7, 3], [9, 2, 1, 8], [8, 1, 2, 9]], dtype=float),
        obs=pd.DataFrame({'group': ['ctrl', 'ctrl', 'treat', 'treat']}, index=['s1', 's2', 's3', 's4']),
        var=pd.DataFrame(index=['G1', 'G2', 'G3', 'G4']),
    ).write_h5ad(output_path)
    plot_path = os.path.join(results_dir, 'bulk_heatmap.png')
    _write_png(plot_path)
    task = AnalysisTask(
        project_id=test_project, module_name='bulk_heatmap', status='completed', output_adata_path=output_path,
        params_json='{"groupby":"group","deg_comparison_label":"treat vs ctrl"}',
    )
    task.save()
    result = ResultFile.create(
        task_id=task.id, project_id=test_project, file_type='png', category='heatmap',
        label='Bulk heatmap', file_path=plot_path,
    )

    source = resolve_source(test_project, 'result_file', result.id)
    assert source['edit_mode'] == 'bulk_heatmap'
    assert source['heatmap_options']['group_columns']['group'] == ['ctrl', 'treat']
    preview = render_preview_data_uri(source, {'heatmap_cmap': 'viridis', 'annotation_column': 'group'})
    assert preview['edit_mode'] == 'data_redraw'
    version = save_figure_version(test_project, source, {'heatmap_cmap': 'viridis', 'annotation_column': 'group'})
    assert version.edit_mode == 'bulk_heatmap'
    assert os.path.isfile(version.png_path)
    assert os.path.isfile(version.svg_path)

    subset = save_figure_version(test_project, source, {
        'heatmap_cmap': 'viridis', 'heatmap_sample_labels': 'all',
        'heatmap_sample_scope': 'selected_groups', 'heatmap_group_column': 'group',
        'heatmap_selected_groups': 'ctrl',
    }, 'ctrl-only heatmap')
    svg_text = open(subset.svg_path, encoding='utf-8').read()
    assert 's1' in svg_text and 's2' in svg_text
    assert 's3' not in svg_text and 's4' not in svg_text


def test_enrichment_version_uses_selected_ontology_colours(test_project):
    from config import Config
    from models import AnalysisTask, ResultFile
    from modules.figure_studio import render_preview_data_uri, resolve_source, save_figure_version

    results_dir = Config.results_dir(test_project)
    os.makedirs(results_dir, exist_ok=True)
    plot_path = os.path.join(results_dir, 'enrichment_ora_kegg.png')
    _write_png(plot_path)
    table_path = os.path.join(results_dir, 'enrichment_ora_kegg_results.csv')
    with open(table_path, 'w', encoding='utf-8') as handle:
        handle.write(
            'Database,Method,Direction,Term,Adjusted P-value,Overlap\n'
            'KEGG,ORA,All,Fatty acid metabolism,0.001,4/120\n'
            'KEGG,ORA,All,PPAR signaling pathway,0.008,3/95\n'
        )
    task = AnalysisTask(project_id=test_project, module_name='bulk_enrichment', status='completed')
    task.save()
    result = ResultFile.create(task_id=task.id, project_id=test_project, file_type='png',
                               category='enrichment', label='KEGG ORA 富集图', file_path=plot_path)

    source = resolve_source(test_project, 'result_file', result.id)
    assert source['edit_mode'] == 'bulk_enrichment'
    preview = render_preview_data_uri(source, {'enrichment_kegg_color': '#123456', 'enrichment_top_n': 10})
    assert preview['edit_mode'] == 'data_redraw'
    version = save_figure_version(test_project, source, {'enrichment_kegg_color': '#123456'})
    assert version.edit_mode == 'bulk_enrichment'
    assert '#123456' in open(version.svg_path, encoding='utf-8').read().lower()


def test_enrichment_chord_source_keeps_semantic_renderer_and_target_selection(test_project):
    from config import Config
    from models import AnalysisTask, ResultFile
    from modules.figure_studio import render_preview_data_uri, resolve_source, save_figure_version

    results_dir = Config.results_dir(test_project)
    os.makedirs(results_dir, exist_ok=True)
    plot_path = os.path.join(results_dir, 'enrichment_ora_go_bp_nh4cl_chord.png')
    _write_png(plot_path)
    table_path = os.path.join(results_dir, 'enrichment_ora_go_bp_nh4cl_results.csv')
    with open(table_path, 'w', encoding='utf-8') as handle:
        handle.write(
            'Database,Method,Direction,Term,Adjusted P-value,Overlap,Genes\n'
            'GO_BP,ORA,All,Extracellular matrix organization (GO:0005581),0.001,4/120,COL1A1;MMP2\n'
            'GO_BP,ORA,All,Cell migration (GO:0016477),0.008,3/120,SPP1;MMP2\n'
        )
    task = AnalysisTask(project_id=test_project, module_name='bulk_enrichment', status='completed')
    task.save()
    result = ResultFile.create(task_id=task.id, project_id=test_project, file_type='png',
                               category='enrichment', label='GO BP ORA · chord', file_path=plot_path)

    source = resolve_source(test_project, 'result_file', result.id)
    assert source['edit_mode'] == 'bulk_enrichment'
    assert source['plot_type'] == 'enrichment_chord'
    preview = render_preview_data_uri(source, {
        'enrichment_pathway_selection': 'selected',
        'enrichment_target_pathways': 'GO:0005581',
        'enrichment_gene_label_strategy': 'all',
    })
    assert preview['edit_mode'] == 'data_redraw'
    version = save_figure_version(test_project, source, {
        'enrichment_pathway_selection': 'selected',
        'enrichment_target_pathways': 'GO:0005581',
    }, 'selected chord pathway')
    svg_text = open(version.svg_path, encoding='utf-8').read()
    assert 'Extracellular' in svg_text
    assert 'Cell migration' not in svg_text


def test_enrichment_dotplot_source_stays_dotplot_in_figure_studio(test_project):
    """The dotplot preview must not silently fall back to the legacy barplot."""
    from config import Config
    from models import AnalysisTask, ResultFile
    from modules.figure_studio import render_preview_data_uri, resolve_source, save_figure_version

    results_dir = Config.results_dir(test_project)
    os.makedirs(results_dir, exist_ok=True)
    # Current analysis runs use the unsuffixed primary image stem; its
    # GeneRatio/Count table is what lets Figure Studio distinguish it from
    # legacy Overlap-only barplots.
    plot_path = os.path.join(results_dir, 'enrichment_ora_go_bp.png')
    _write_png(plot_path)
    table_path = os.path.join(results_dir, 'enrichment_ora_go_bp_results.csv')
    with open(table_path, 'w', encoding='utf-8') as handle:
        handle.write(
            'Database,Method,Direction,Term,Adjusted P-value,GeneRatio,Count,Genes\n'
            'GO_BP,ORA,All,Extracellular matrix organization (GO:0005581),0.001,0.05,12,COL1A1;MMP2\n'
            'GO_BP,ORA,All,Cell migration (GO:0016477),0.008,0.03,8,SPP1;MMP2\n'
        )
    task = AnalysisTask(project_id=test_project, module_name='bulk_enrichment', status='completed')
    task.save()
    result = ResultFile.create(task_id=task.id, project_id=test_project, file_type='png',
                               category='enrichment', label='GO BP ORA · dotplot', file_path=plot_path)

    source = resolve_source(test_project, 'result_file', result.id)
    assert source['edit_mode'] == 'bulk_enrichment'
    assert source['plot_type'] == 'enrichment_dotplot'
    preview = render_preview_data_uri(source, {
        'enrichment_pathway_selection': 'selected',
        'enrichment_target_pathways': 'GO:0005581',
    })
    assert preview['edit_mode'] == 'data_redraw'
    version = save_figure_version(test_project, source, {
        'enrichment_pathway_selection': 'selected',
        'enrichment_target_pathways': 'GO:0005581',
    }, 'selected dotplot pathway')
    svg_text = open(version.svg_path, encoding='utf-8').read()
    assert 'Gene ratio' in svg_text
    assert 'Cell migration' not in svg_text
    assert '-log10(adjusted P-value)' not in svg_text


def test_enrichment_overview_source_uses_integrated_table_and_database_scope(test_project):
    from config import Config
    from models import AnalysisTask, ResultFile
    from modules.figure_studio import render_preview_data_uri, resolve_source

    results_dir = Config.results_dir(test_project)
    plots_dir = Config.plots_dir(test_project)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    plot_path = os.path.join(plots_dir, 'enrichment_overview_ora_demo_all.png')
    _write_png(plot_path)
    table_path = os.path.join(results_dir, 'enrichment_integrated_results.csv')
    with open(table_path, 'w', encoding='utf-8') as handle:
        handle.write(
            'Comparison,Database,Method,Direction,Term,Enrichment FDR,GeneRatio,num,Genes\n'
            'demo,GO_BP,ORA,All,Matrix organization (GO:0001),0.001,0.05,12,A;B\n'
            'demo,KEGG,ORA,All,ECM interaction (KEGG:hsa0007),0.004,0.03,8,C;D\n'
        )
    task = AnalysisTask(project_id=test_project, module_name='bulk_enrichment', status='completed')
    task.save()
    result = ResultFile.create(task_id=task.id, project_id=test_project, file_type='png',
                               category='enrichment_overview', label='多数据库富集概览', file_path=plot_path)

    source = resolve_source(test_project, 'result_file', result.id)
    assert source['edit_mode'] == 'bulk_enrichment_overview'
    assert source['plot_type'] == 'enrichment_overview'
    preview = render_preview_data_uri(source, {
        'enrichment_database_scope': 'KEGG',
        'enrichment_pathway_selection': 'selected',
        'enrichment_target_pathways': 'KEGG:hsa0007',
    })
    assert preview['edit_mode'] == 'data_redraw'


def test_correlation_version_supports_colormap_editing(test_project):
    from config import Config
    from models import AnalysisTask, ResultFile
    from modules.figure_studio import render_preview_data_uri, resolve_source, save_figure_version

    results_dir = Config.results_dir(test_project)
    os.makedirs(results_dir, exist_ok=True)
    plot_path = os.path.join(results_dir, 'bulk_corr_heatmap.png')
    _write_png(plot_path)
    pairs_path = os.path.join(results_dir, 'bulk_correlation_pairs.csv')
    with open(pairs_path, 'w', encoding='utf-8') as handle:
        handle.write(
            'sample_1,sample_2,group_1,group_2,pearson_r\n'
            's1,s2,ctrl,ctrl,0.98\ns1,s3,ctrl,treat,0.72\ns2,s3,ctrl,treat,0.75\n'
        )
    task = AnalysisTask(project_id=test_project, module_name='bulk_heatmap', status='completed')
    task.save()
    result = ResultFile.create(task_id=task.id, project_id=test_project, file_type='png',
                               category='heatmap', label='样本相关性热图', file_path=plot_path)

    source = resolve_source(test_project, 'result_file', result.id)
    assert source['edit_mode'] == 'bulk_correlation'
    preview = render_preview_data_uri(source, {'corr_colorscale': 'viridis'})
    assert preview['edit_mode'] == 'data_redraw'
    version = save_figure_version(test_project, source, {'corr_colorscale': 'viridis'})
    assert version.edit_mode == 'bulk_correlation'
    assert os.path.isfile(version.png_path)


def test_style_only_image_exports_png_and_svg_without_semantic_rewrite(test_project):
    from config import Config
    from models import FigureAsset
    from modules.figure_studio import resolve_source, save_figure_version

    upload_dir = os.path.join(Config.project_dir(test_project), 'figure_studio', 'uploads')
    os.makedirs(upload_dir, exist_ok=True)
    original_path = os.path.join(upload_dir, 'uploaded.png')
    _write_png(original_path, color=(220, 230, 240))
    asset = FigureAsset.create(test_project, 'uploaded image', 'png', original_path)

    version = save_figure_version(
        test_project, resolve_source(test_project, 'asset', asset.id),
        {'title': 'Appearance-only title', 'background': '#FFFFFF', 'legend_text': 'Source: uploaded image'},
    )
    assert version.edit_mode == 'style_only'
    assert os.path.isfile(version.png_path)
    assert os.path.isfile(version.svg_path)
    assert version.png_path != '.'
    assert 'data:image/png;base64,' in open(version.svg_path, encoding='utf-8').read()


@pytest.fixture
def figure_client(tmp_path, monkeypatch):
    from config import Config

    data_dir = tmp_path / 'data'
    data_dir.mkdir()
    monkeypatch.setattr(Config, 'DATA_DIR', str(data_dir))
    monkeypatch.setattr(Config, 'DB_PATH', str(tmp_path / 'figure_studio.db'))
    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client


def _make_project(project_id):
    from config import Config
    from database import get_conn

    os.makedirs(Config.project_dir(project_id), exist_ok=True)
    conn = get_conn()
    try:
        conn.execute('INSERT INTO projects (id, name, status) VALUES (?, ?, ?)', (project_id, project_id, 'active'))
        conn.commit()
    finally:
        conn.close()


def test_upload_preview_and_save_are_project_scoped(figure_client):
    _make_project('figure_project')
    _make_project('other_project')
    from PIL import Image

    page = figure_client.get('/projects/figure_project/figure-studio')
    assert page.status_code == 200
    assert '图形美化工作台'.encode() in page.data

    buffer = io.BytesIO()
    Image.new('RGB', (24, 24), (200, 210, 220)).save(buffer, format='PNG')
    buffer.seek(0)
    upload = figure_client.post(
        '/projects/figure_project/figure-studio/upload',
        data={'image': (buffer, 'figure.png')}, content_type='multipart/form-data',
    )
    assert upload.status_code == 201
    source = upload.get_json()['source']

    preview = figure_client.post('/projects/figure_project/figure-studio/preview', json={
        'source_kind': source['kind'], 'source_id': source['id'], 'style': {'title': 'Uploaded'},
    })
    assert preview.status_code == 200
    assert preview.get_json()['edit_mode'] == 'style_only'

    saved = figure_client.post('/projects/figure_project/figure-studio/save', json={
        'source_kind': source['kind'], 'source_id': source['id'], 'style': {'title': 'Uploaded'},
    })
    assert saved.status_code == 201
    version = saved.get_json()['version']
    assert figure_client.get(version['png_url']).status_code == 200
    assert figure_client.get(f'/projects/other_project/figure-studio/versions/{version["id"]}/png').status_code == 302
