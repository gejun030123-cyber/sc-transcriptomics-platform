"""Regression tests for expression-matrix to h5ad imports."""

import json
import os

import pandas as pd


def _write_expression_matrix(path):
    pd.DataFrame(
        {'cell_a': [1, 0, 3], 'cell_b': [0, 2, 1]},
        index=['GeneA', 'GeneB', 'GeneC'],
    ).to_csv(path)


def test_expression_matrix_is_offered_for_single_cell_import(test_project):
    from config import Config
    from routes.upload import check_importable_sc_files

    uploads_dir = Config.uploads_dir(test_project)
    _write_expression_matrix(os.path.join(uploads_dir, 'counts.csv'))
    # A 10x metadata component must not be misidentified as a standalone
    # expression matrix.
    with open(os.path.join(uploads_dir, 'barcodes.tsv'), 'w', encoding='utf-8') as handle:
        handle.write('cell_a\ncell_b\n')

    result = check_importable_sc_files(uploads_dir)

    assert result == {
        'has_importable': True,
        'files': [{
            'name': 'counts.csv',
            'format': 'expression_matrix',
            'size_mb': 0.0,
        }],
    }


def test_expression_matrix_import_route_submits_conversion_task(test_project, monkeypatch):
    from app import create_app
    from config import Config
    from models import AnalysisTask
    import worker

    _write_expression_matrix(os.path.join(Config.uploads_dir(test_project), 'counts.csv'))
    submitted = {}

    def capture_submit(*args):
        submitted['args'] = args
        return True

    monkeypatch.setattr(Config, 'PLATFORM_ACCESS_PASSWORD', '')
    monkeypatch.setattr(worker, 'submit_task', capture_submit)
    app = create_app()
    app.config['TESTING'] = True

    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-sc',
            data={'source_file': 'counts.csv', 'input_format': 'expression_matrix'},
        )

    assert response.status_code == 200
    task = AnalysisTask.get_by_id(response.get_json()['task_id'])
    assert task.module_name == 'convert_10x'
    assert json.loads(task.params_json)['input_format'] == 'expression_matrix'
    assert submitted['args'][2] == 'convert_10x'
    assert submitted['args'][3]['input_format'] == 'expression_matrix'


def test_upload_page_autostarts_conversion_for_complete_10x_upload(test_project, monkeypatch):
    from app import create_app
    from config import Config

    monkeypatch.setattr(Config, 'PLATFORM_ACCESS_PASSWORD', '')
    app = create_app()
    app.config['TESTING'] = True

    with app.test_client() as client:
        response = client.get(f'/projects/{test_project}/upload')

    html = response.get_data(as_text=True)
    assert 'check10x({autoConvert: true})' in html
    assert 'function check10x({autoConvert = false} = {})' in html
    assert 'startConvert10x();' in html


def test_expression_matrix_import_writes_standard_h5ad_with_sample_id(tmp_path):
    import anndata
    from modules.convert_10x import Convert10x

    project_dir = tmp_path / 'project'
    uploads_dir = project_dir / 'uploads'
    uploads_dir.mkdir(parents=True)
    source = uploads_dir / 'counts.csv'
    _write_expression_matrix(source)

    result = Convert10x(
        project_dir=str(project_dir),
        params={'source_path': str(source), 'input_format': 'expression_matrix'},
        progress_callback=lambda _percent, _message: None,
    ).run(str(source))

    output = result['output_adata']
    converted = anndata.read_h5ad(output)
    assert output.endswith('counts_imported.h5ad')
    assert result['summary']['input_format'] == 'expression_matrix'
    assert converted.shape == (2, 3)
    assert list(converted.obs['sample_id'].unique()) == ['counts']
    assert 'counts' in converted.layers
