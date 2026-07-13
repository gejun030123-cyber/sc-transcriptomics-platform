import io
import json
import zipfile


def _zip_bytes(names):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w') as archive:
        for name in names:
            archive.writestr(name, '')
    buffer.seek(0)
    return buffer


def test_batch_zip_endpoint_accepts_two_batches(test_project, monkeypatch):
    from app import create_app
    from config import Config
    from models import AnalysisTask
    from werkzeug.datastructures import MultiDict

    submitted = []
    import worker
    monkeypatch.setattr(worker, 'submit_task', lambda *args: submitted.append(args) or True)
    app = create_app()
    app.config['TESTING'] = True

    files = [
        ('batch_zip', (_zip_bytes(['filtered_feature_bc_matrix/matrix.mtx',
                                   'filtered_feature_bc_matrix/barcodes.tsv',
                                   'filtered_feature_bc_matrix/features.tsv']), 'control.zip')),
        ('batch_zip', (_zip_bytes(['filtered_feature_bc_matrix/matrix.mtx',
                                   'filtered_feature_bc_matrix/barcodes.tsv',
                                   'filtered_feature_bc_matrix/features.tsv']), 'treatment.zip')),
    ]
    with app.test_client() as client:
        response = client.post(
            f'/projects/{test_project}/upload/import-10x-batches',
            data=MultiDict(files + [
                ('batch_name', 'Control'),
                ('batch_name', 'Treatment'),
            ]),
            content_type='multipart/form-data',
        )

    assert response.status_code == 200
    payload = response.get_json()
    assert payload['batch_names'] == ['Control', 'Treatment']
    assert len(submitted) == 1
    task = AnalysisTask.get_by_id(payload['task_id'])
    params = json.loads(task.params_json)
    assert [item['batch_name'] for item in params['batch_sources']] == ['Control', 'Treatment']
    assert all(item['zip_path'].startswith(Config.uploads_dir(test_project)) for item in params['batch_sources'])


def test_batch_zip_extraction_rejects_path_traversal(tmp_path):
    malicious = _zip_bytes(['../escape/matrix.mtx'])
    from modules.convert_10x import _safe_extract_zip
    import pytest
    archive = tmp_path / 'bad.zip'
    archive.write_bytes(malicious.getvalue())

    with pytest.raises(ValueError):
        _safe_extract_zip(str(archive), str(tmp_path / 'extract'))
