"""Regression tests for selecting one contrast before pathway enrichment."""

import json
import os

import pytest

from app import create_app


@pytest.fixture
def client(tmp_path, monkeypatch):
    from config import Config

    data_dir = tmp_path / 'data'
    data_dir.mkdir()
    monkeypatch.setattr(Config, 'DATA_DIR', str(data_dir))
    monkeypatch.setattr(Config, 'DB_PATH', str(tmp_path / 'test.db'))
    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as test_client:
        yield test_client


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


def _make_deg_table(project_id, suffix, label):
    from config import Config
    from models import AnalysisTask, ResultFile

    results_dir = Config.results_dir(project_id)
    os.makedirs(results_dir, exist_ok=True)
    path = os.path.join(results_dir, f'bulk_deg_results{suffix}.csv')
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('gene,log2FC,padj,regulation\nTP53,2.2,0.01,Up\nMYC,-1.8,0.02,Down\n')
    task = AnalysisTask(project_id=project_id, module_name='bulk_deg', status='completed')
    task.save()
    return ResultFile.create(task_id=task.id, project_id=project_id, file_type='csv',
                             category='table', label=label, file_path=path)


def test_enrichment_lists_each_deg_contrast_and_excludes_merged_table(client):
    from config import Config
    from models import AnalysisTask, ResultFile
    from routes.analysis import _bulk_deg_enrichment_sources

    _make_project('enrich_project')
    first = _make_deg_table('enrich_project', '_0', '差异表达基因列表 (Ctrl vs Treat)')
    second = _make_deg_table('enrich_project', '_1', '差异表达基因列表 (Dose vs Treat)')

    # A merged table is not a valid single biological contrast for ORA/GSEA.
    merged_path = os.path.join(Config.results_dir('enrich_project'), 'bulk_deg_merged_comparisons.csv')
    with open(merged_path, 'w', encoding='utf-8') as handle:
        handle.write('gene,log2FC\nTP53,2.2\n')
    merged_task = AnalysisTask(project_id='enrich_project', module_name='bulk_deg', status='completed')
    merged_task.save()
    ResultFile.create(task_id=merged_task.id, project_id='enrich_project', file_type='csv',
                      category='table', label='多比较合并结果', file_path=merged_path)

    sources = _bulk_deg_enrichment_sources('enrich_project')
    assert {item['path'] for item in sources} == {first.file_path, second.file_path}

    page = client.get('/projects/enrich_project/analyze/bulk_enrichment')
    assert page.status_code == 200
    assert 'Ctrl vs Treat'.encode() in page.data
    assert 'Dose vs Treat'.encode() in page.data
    assert '多比较合并结果'.encode() not in page.data


def test_enrichment_rejects_blank_deg_choice_before_submitting_task(client):
    from models import AnalysisTask

    _make_project('enrich_project')
    _make_deg_table('enrich_project', '_0', '差异表达基因列表 (Ctrl vs Treat)')

    response = client.post('/projects/enrich_project/analyze/bulk_enrichment', data={
        'input_source': '', 'custom_genes': '',
    }, follow_redirects=True)

    assert response.status_code == 200
    assert '请选择一个 DEG 比较结果'.encode() in response.data
    assert len(AnalysisTask.get_by_project('enrich_project')) == 1


def test_enrichment_backfills_legacy_result_provenance(client):
    from config import Config
    from models import AnalysisTask, ResultFile
    from modules.bulk_enrichment import _backfill_legacy_enrichment_metadata

    project_id = 'enrich_project'
    _make_project(project_id)
    deg_result = _make_deg_table(
        project_id, '_0', '差异表达基因列表 (Ctrl vs Treat)',
    )
    legacy_task = AnalysisTask(
        project_id=project_id,
        module_name='bulk_enrichment',
        status='completed',
        params_json=json.dumps({
            'input_source': deg_result.file_path,
            'method': 'ORA',
            'database': 'GO_BP',
        }),
    )
    legacy_task.save()
    legacy_path = os.path.join(Config.results_dir(project_id), 'enrichment_ora_results.csv')
    with open(legacy_path, 'w', encoding='utf-8') as handle:
        handle.write('Term,Adjusted P-value\nresponse to virus,0.01\n')
    ResultFile.create(
        task_id=legacy_task.id, project_id=project_id, file_type='csv', category='table',
        label='ORA 富集结果', file_path=legacy_path,
    )

    updated = _backfill_legacy_enrichment_metadata(
        Config.project_dir(project_id), Config.results_dir(project_id),
    )

    assert updated == [legacy_path]
    with open(legacy_path, encoding='utf-8') as handle:
        header = handle.readline().strip().split(',')
        row = handle.readline().strip().split(',')
    assert header[:4] == ['Comparison', 'Database', 'Method', 'Term']
    assert row[:3] == ['Ctrl vs Treat', 'GO_BP', 'ORA']
    assert 'Direction' in header


def test_heatmap_form_exposes_display_scope_and_group_picker(client):
    _make_project('enrich_project')

    response = client.get('/projects/enrich_project/analyze/bulk_heatmap')

    assert response.status_code == 200
    assert '热图展示样本范围'.encode() in response.data
    assert '仅 DEG 两组'.encode() in response.data
    assert b'dynamic-multiselect' in response.data
