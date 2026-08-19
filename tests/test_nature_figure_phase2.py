"""Phase 2 regression tests for composer, new templates, export and gate."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _deg_fixture():
    return pd.DataFrame({
        'gene': [f'Gene{index:02d}' for index in range(30)],
        'mean_expression': np.geomspace(1, 10000, 30),
        'log2FC': np.r_[np.linspace(-3, -1.1, 8), np.zeros(14), np.linspace(1.1, 3, 8)],
        'padj': np.r_[np.geomspace(1e-5, 0.02, 8), np.full(14, 0.8),
                      np.geomspace(0.02, 1e-5, 8)],
    })


def test_director_renders_phase2_templates_with_fixed_encodings():
    import matplotlib.pyplot as plt
    from figure_engine import NatureFigureDirector

    director = NatureFigureDirector()
    figures = []
    ma = director.render(director.create_spec('ma', label_n=4), _deg_fixture())
    figures.append(ma)
    assert ma._nature_encodings['color'] == 'DEG direction'

    corr = np.array([[1.0, 0.91, 0.32], [0.91, 1.0, 0.28], [0.32, 0.28, 1.0]])
    correlation = director.render(
        director.create_spec('correlation', max_col_labels=3),
        {'correlation_matrix': corr, 'sample_labels': ['C1', 'C2', 'T1'],
         'groups': ['Control', 'Control', 'Treatment']},
    )
    figures.append(correlation)
    assert correlation._nature_encodings['color'] == 'correlation'

    scores = np.arange(24, dtype=float).reshape(4, 6)
    gsva = director.render(
        director.create_spec('gsva', max_row_labels=4, max_col_labels=6),
        {'scores': scores, 'pathways': ['P1', 'P2', 'P3', 'P4'],
         'samples': [f'S{index}' for index in range(6)]},
    )
    figures.append(gsva)
    assert 'pathway-wise z-score' in gsva._nature_encodings['color']

    upset = director.render(
        director.create_spec('upset', width='double', top_intersections=8),
        {'A': {'g1', 'g2', 'g3'}, 'B': {'g2', 'g3', 'g4'}, 'C': {'g3', 'g5'}},
    )
    figures.append(upset)
    assert upset._nature_encodings['matrix'] == 'set membership'

    wgcna = director.render(
        director.create_spec('wgcna', width='double'),
        {'correlation_matrix': [[0.8, -0.2], [-0.5, 0.62]],
         'pvalues': [[0.001, 0.4], [0.02, 0.008]],
         'module_labels': ['MEblue', 'MEbrown'], 'trait_labels': ['Age', 'Treatment']},
    )
    figures.append(wgcna)
    assert wgcna._nature_encodings['color'] == 'correlation'
    for figure in figures:
        plt.close(figure)


def test_composer_uses_subfigures_uniform_panel_labels_and_shared_legend():
    import matplotlib.pyplot as plt
    from figure_engine import FigurePanel, FigureSpec, NatureFigureComposer, NatureFigureDirector

    director = NatureFigureDirector()
    pca_spec = director.create_spec('pca', show_legend=True)
    panels = [
        FigurePanel(
            plot_type='pca', spec=pca_spec, row=0, column=0,
            data={'coordinates': [[-1, 0], [-0.8, 0.2], [1, 0], [0.9, -0.1]],
                  'groups': ['Control', 'Control', 'Treatment', 'Treatment'],
                  'explained_variance': [0.5, 0.2]},
        ),
        FigurePanel(plot_type='ma', data=_deg_fixture(), row=0, column=1,
                    spec=director.create_spec('ma', show_legend=True, label_n=0)),
    ]
    composite_spec = FigureSpec(
        plot_type='composite', width='double', height_mm=105,
        formats=('svg', 'pdf', 'png'),
    )
    figure = NatureFigureComposer(director).compose(
        panels, composite_spec, nrows=1, ncols=2, shared_legend=True)

    assert figure._nature_panel_labels == ['a', 'b']
    assert figure._nature_panel_grid == {'rows': 1, 'columns': 2}
    assert len(figure.subfigs) == 2
    assert len(figure.legends) == 1
    plt.close(figure)


def test_pdf_tiff_registration_and_strict_gate(tmp_path):
    import matplotlib.pyplot as plt
    from PIL import Image
    from figure_engine import (
        FigureReadinessError, FigureValidator, NatureFigureDirector, export_figure,
    )

    director = NatureFigureDirector()
    spec = director.create_spec(
        'ma', formats=('svg', 'pdf', 'png', 'tiff'), label_n=0,
    )
    figure = director.render(spec, _deg_fixture())
    paths, report = export_figure(figure, tmp_path / 'ma', spec)

    assert set(paths) == {'svg', 'pdf', 'png', 'tiff'}
    assert Path(paths['pdf']).read_bytes().startswith(b'%PDF-')
    with Image.open(paths['tiff']) as image:
        assert image.format == 'TIFF'
        assert image.info['dpi'][0] == pytest.approx(600, abs=1)
    assert report.ready

    figure.axes[0].text(0.5, 0.5, 'too small', fontsize=2)
    failed = FigureValidator().validate(figure, spec, export_paths=paths)
    assert not failed.ready
    with pytest.raises(FigureReadinessError):
        FigureValidator().require(failed)
    plt.close(figure)


def test_colour_vision_simulation_and_visual_signature(tmp_path):
    import matplotlib.pyplot as plt
    from PIL import Image
    from figure_engine import (
        NatureFigureDirector, audit_palette, compare_signatures,
        export_figure, figure_signature, simulate_image,
    )

    director = NatureFigureDirector()
    spec = director.create_spec('ma', formats=('png',), label_n=0)
    first = director.render(spec, _deg_fixture())
    second = director.render(spec, _deg_fixture())
    signature = figure_signature(first)
    comparison = compare_signatures(figure_signature(second), signature)
    assert comparison == {
        'pass': True, 'hash_distance': 0, 'maximum_distance': 8,
        'structural_match': True,
    }
    baseline = json.loads(
        (Path(__file__).parent / 'snapshots' / 'nature_phase2.json').read_text()
    )['ma_89mm']
    assert compare_signatures(signature, baseline, maximum_distance=8)['pass']
    paths, _ = export_figure(first, tmp_path / 'ma', spec)
    preview = simulate_image(paths['png'], tmp_path / 'ma_deuteranopia.png')
    with Image.open(preview) as image:
        assert image.size[0] > 0
    assert audit_palette(['#4C78A8', '#B65C5C'])['pass']
    plt.close(first)
    plt.close(second)


def test_visualization_modes_lock_nature_contract_and_accept_submission_formats():
    from modules.figure_style import normalize_visualization_params

    normalized = normalize_visualization_params({
        'figure_mode': 'nature_portfolio', 'theme': 'dark',
        'width_profile': 'double', 'nature_style': 'nature_communications',
        'static_formats': ['svg', 'pdf', 'png', 'tiff'],
    })
    assert normalized['theme'] == 'nature'
    assert normalized['bg_color'] == 'white'
    assert normalized['width_profile'] == 'double'
    assert normalized['static_formats'] == ('svg', 'pdf', 'png', 'tiff')
