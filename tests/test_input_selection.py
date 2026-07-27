from types import SimpleNamespace


def _task(module_name, path, task_id, finished_at):
    return SimpleNamespace(
        module_name=module_name,
        output_adata_path=str(path),
        id=task_id,
        finished_at=finished_at,
        started_at=finished_at,
        status='completed',
    )


def test_build_input_options_recommends_latest_valid_upstream(tmp_path):
    from routes.analysis import build_input_options

    qc_path = tmp_path / 'qc_output.h5ad'
    normalize_path = tmp_path / 'normalize_output.h5ad'
    qc_path.write_bytes(b'qc')
    normalize_path.write_bytes(b'normalize')
    tasks = [
        _task('qc', qc_path, 'task-qc', '2026-07-14T10:00:00'),
        _task('normalize', normalize_path, 'task-normalize', '2026-07-14T11:00:00'),
    ]

    options, recommended = build_input_options('hvg', tasks, [])

    assert recommended['module_name'] == 'normalize'
    assert recommended['path'] == str(normalize_path.resolve())
    assert options[0]['recommended'] is True
    assert '最新输出' in recommended['reason']


def test_build_input_options_uses_upload_for_first_step(tmp_path):
    from routes.analysis import build_input_options

    source = tmp_path / 'raw.h5ad'
    source.write_bytes(b'raw')
    options, recommended = build_input_options(
        'qc', [], [{'name': 'raw.h5ad', 'path': str(source)}]
    )

    assert recommended['source'] == 'upload'
    assert recommended['path'] == str(source.resolve())
    assert options[0]['recommended'] is True
