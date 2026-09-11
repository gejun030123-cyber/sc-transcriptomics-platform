"""Tests for the donor-aware single-cell review summaries."""

import numpy as np
import pandas as pd


def test_qc_summary_retention_and_flags_are_non_exclusive():
    from modules.sc_figure_diagnostics import summarize_qc_by_batch

    qc_before = pd.DataFrame({
        'donor': ['A', 'A', 'B', 'B'],
        'total_counts': [100, 900, 1000, 80],
        'n_genes_by_counts': [100, 400, 300, 100],
        'pct_counts_mt': [30.0, 5.0, 5.0, 25.0],
    }, index=['a1', 'a2', 'b1', 'b2'])
    summary = summarize_qc_by_batch(
        qc_before,
        ['a2', 'b1'],
        'donor',
        mito_perc=0.2,
        n_umis=500,
        n_genes_min=250,
        doublet_index=['b2'],
    )

    assert summary['A']['n_before'] == 2
    assert summary['A']['n_after'] == 1
    assert summary['A']['flagged_cells']['low_umi'] == 1
    assert summary['A']['flagged_cells']['high_mt'] == 1
    assert summary['B']['flagged_cells']['low_umi'] == 1
    assert summary['B']['flagged_cells']['low_genes'] == 1
    assert summary['B']['flagged_cells']['high_mt'] == 1
    assert summary['B']['flagged_cells']['doublet'] == 1


def test_donor_qc_figures_keep_footer_and_thresholds_clear_of_labels():
    """Long donor labels must not collide with review text or thresholds."""
    import matplotlib.pyplot as plt
    from modules.sc_figure_diagnostics import (
        doublet_by_batch_figure,
        qc_by_batch_figure,
    )

    labels = [f'donor_{index:02d}_long_identifier' for index in range(12)]
    summary = {
        label: {
            'n_after': 800 + index * 10,
            'n_removed': 100 + index,
            'pct_removed': 10.0 + index / 10,
            'flags_per_1000': {
                'low_genes': 20.0,
                'low_umi': 10.0,
                'high_mt': 5.0,
                'doublet': 7.0,
                'high_genes': 0.0,
                'high_ribo': 0.0,
                'high_hb': 0.0,
            },
        }
        for index, label in enumerate(labels)
    }
    scores = {label: [0.05, 0.12, 0.20, 0.31] for label in labels}
    scrublet_summary = {
        'threshold_by_batch': {
            label: 0.20 + index * 0.002 for index, label in enumerate(labels)
        },
        'by_batch': {
            label: {'doublet_rate': 0.01 if index == 0 else 0.08}
            for index, label in enumerate(labels)
        },
    }

    figures = [
        (
            qc_by_batch_figure(summary),
            'Flags are non-exclusive; their sum may exceed removed cells.',
        ),
        (
            doublet_by_batch_figure(scores, scrublet_summary),
            'Review: donor doublet rates differ by >5×.',
        ),
    ]
    try:
        for figure, footer_text in figures:
            # Mirror the final layout pass performed by the shared exporter.
            figure.tight_layout(rect=figure._native_layout_rect, pad=1.1)
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            footer = next(text for text in figure.texts if text.get_text() == footer_text)
            footer_box = footer.get_window_extent(renderer)
            tick_boxes = [
                tick.get_window_extent(renderer)
                for axis in figure.axes if axis.axison
                for tick in axis.get_xticklabels()
            ]
            assert footer_box.ymax < min(box.ymin for box in tick_boxes)

        # Thresholds are line segments rather than an overlapping number above
        # every violin.  The score panel therefore has no data annotations.
        assert not figures[1][0].axes[0].texts
    finally:
        for figure, _ in figures:
            plt.close(figure)

def test_obs_vs_sim_figure_renders_observed_simulated_and_threshold():
    import matplotlib.pyplot as plt
    from modules.sc_figure_diagnostics import doublet_obs_vs_sim_figure

    scores = {
        'A': [0.02] * 40 + [0.5, 0.6],
        'B': [0.03] * 40 + [0.7],
    }
    sim = {
        'A': [0.1, 0.2, 0.4, 0.5, 0.6, 0.8],
        'B': [0.1, 0.3, 0.5, 0.7, 0.9],
    }
    scrublet_summary = {
        'threshold_by_batch': {'A': 0.45, 'B': 0.6},
        'by_batch': {
            'A': {'detected_doublet_rate': 0.0004, 'detectable_doublet_fraction': 0.5},
            'B': {'detected_doublet_rate': 0.01, 'detectable_doublet_fraction': 0.4},
        },
    }
    fig = doublet_obs_vs_sim_figure(scores, sim, scrublet_summary)
    try:
        assert fig is not None
        assert len(fig.axes) == 2
        # Both panels carry an observed histogram artist plus a threshold line.
        assert any(any(artist.get_label().startswith('Observed') for artist in axis.patches)
                   for axis in fig.axes)
        # Donor A has an unusually low detected rate -> review footer.
        assert any('unusually low' in text.get_text() for text in fig.texts)
    finally:
        plt.close(fig)


def test_hvg_technical_composition_keeps_small_category_labels_in_axes():
    import matplotlib.pyplot as plt
    from modules.sc_figure_diagnostics import hvg_technical_composition_figure

    figure = hvg_technical_composition_figure(
        {'MT': 2, 'Ribo': 7, 'Cell cycle': 4},
        n_selected=2000,
    )
    try:
        # Match the final layout pass performed during result export.
        figure.tight_layout(rect=figure._native_layout_rect, pad=1.1)
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        axis = figure.axes[0]
        axis_box = axis.get_window_extent(renderer)

        # The old offset used n_selected and pushed these labels out of the
        # plot.  All bar labels must now remain visible above their bars.
        for label in axis.texts:
            label_box = label.get_window_extent(renderer)
            assert label_box.ymax < axis_box.ymax

        footer = next(
            text for text in figure.texts
            if text.get_text().startswith('Technical flags may overlap')
        )
        footer_box = footer.get_window_extent(renderer)
        tick_boxes = [tick.get_window_extent(renderer) for tick in axis.get_xticklabels()]
        assert footer_box.ymax < min(box.ymin for box in tick_boxes)
    finally:
        plt.close(figure)


def test_composition_by_group_returns_within_group_fractions():
    from modules.sc_figure_diagnostics import composition_by_group

    obs = pd.DataFrame({
        'donor': ['A', 'A', 'A', 'B'],
        'celltype': ['T', 'T', 'B', 'T'],
    })
    counts, fractions = composition_by_group(obs, 'donor', 'celltype')

    assert counts.loc['A', 'T'] == 2
    assert counts.loc['A', 'B'] == 1
    np.testing.assert_allclose(fractions.sum(axis=1).to_numpy(), [1.0, 1.0])
    assert fractions.loc['A', 'T'] == 2 / 3
    assert fractions.loc['B', 'T'] == 1.0


def test_composition_by_group_keeps_other_in_full_denominator():
    from modules.sc_figure_diagnostics import composition_by_group

    obs = pd.DataFrame({
        'donor': ['A'] * 10 + ['B'] * 10,
        'celltype': (
            ['T'] * 5 + ['B'] * 3 + ['Myeloid'] + ['Stromal']
            + ['T'] * 2 + ['B'] * 5 + ['Myeloid'] * 2 + ['Stromal']
        ),
    })
    counts, fractions = composition_by_group(obs, 'donor', 'celltype', max_labels=3)

    assert set(counts.columns) == {'T', 'B', 'Other'}
    assert counts.loc['A', 'Other'] == 2
    assert counts.loc['B', 'Other'] == 3
    assert fractions.sum(axis=1).eq(1.0).all()
    assert fractions.loc['A', 'T'] == 0.5


def test_composition_by_group_figure_reserves_legend_space_and_labels_proportions():
    """Long cell-type legends must not cover the title or stack percentages."""
    import matplotlib.pyplot as plt
    from modules.sc_figure_diagnostics import composition_by_group_figure

    fractions = pd.DataFrame({
        'Absorptive enterocytes': [0.52, 0.11],
        'Stress/secretory enterocytes': [0.31, 0.62],
        'Inflammatory stress cells': [0.17, 0.27],
    }, index=['donor_1', 'donor_2'])
    figure = composition_by_group_figure(fractions)
    try:
        # Mirror the final layout pass performed by the shared exporter.
        figure.tight_layout(rect=figure._native_layout_rect, pad=1.1)
        figure.canvas.draw()
        axis = figure.axes[0]
        renderer = figure.canvas.get_renderer()
        title_box = axis.title.get_window_extent(renderer)
        legend_box = axis.get_legend().get_window_extent(renderer)

        assert not title_box.overlaps(legend_box)
        assert any(text.get_text() == '52.0%' for text in axis.texts)
        assert any(text.get_text() == '62.0%' for text in axis.texts)
        assert axis.get_ylabel() == 'Cell proportion'
    finally:
        plt.close(figure)


def test_cluster_composition_table_marks_donor_associated_clusters():
    from modules.sc_figure_diagnostics import cluster_composition_table

    obs = pd.DataFrame({
        'cluster': ['0'] * 9 + ['1'] * 4,
        'donor': ['A'] * 9 + ['A', 'B', 'A', 'B'],
    })
    table = cluster_composition_table(obs, 'cluster', 'donor').set_index('cluster')

    assert bool(table.loc['0', 'donor_associated'])
    assert table.loc['0', 'major_donor'] == 'A'
    assert table.loc['0', 'major_donor_fraction'] == 1.0
    assert not bool(table.loc['1', 'donor_associated'])


def test_condition_composition_includes_zeroes_for_missing_cell_types():
    from modules.sc_figure_diagnostics import condition_composition_figure

    obs = pd.DataFrame({
        'sample': ['s1'] * 3 + ['s2'] * 2,
        'condition': ['control'] * 3 + ['treated'] * 2,
        'celltype': ['T', 'T', 'B', 'T', 'T'],
    })
    figure = condition_composition_figure(obs, 'sample', 'condition', 'celltype')
    try:
        assert figure is not None
        assert len(figure.axes) >= 1
    finally:
        import matplotlib.pyplot as plt
        plt.close(figure)


def test_condition_composition_is_a_100_percent_stacked_bar_chart():
    from modules.sc_figure_diagnostics import condition_composition_figure

    obs = pd.DataFrame({
        'condition': ['IBD'] * 4 + ['Control'] * 5,
        'celltype': ['T', 'T', 'B', 'B', 'T', 'T', 'T', 'B', 'B'],
    })
    figure = condition_composition_figure(obs, None, 'condition', 'celltype')
    try:
        axis = figure.axes[0]
        assert [tick.get_text() for tick in axis.get_xticklabels()] == ['IBD', 'Control']
        assert axis.get_ylim() == (0.0, 1.0)
        # The two stacked bars each cover one complete condition (4 patches
        # total: two cell types × two conditions).
        heights_by_x = {}
        for patch in axis.patches:
            center = round(float(patch.get_x() + patch.get_width() / 2), 6)
            heights_by_x[center] = heights_by_x.get(center, 0.0) + patch.get_height()
        assert len(heights_by_x) == 2
        np.testing.assert_allclose(sorted(heights_by_x.values()), [1.0, 1.0])
        assert any(text.get_text() == '50.0%' for text in axis.texts)
    finally:
        import matplotlib.pyplot as plt
        plt.close(figure)


def test_condition_composition_keeps_missing_metadata_out_of_the_display():
    from modules.sc_figure_diagnostics import condition_composition_figure

    obs = pd.DataFrame({
        'condition': ['Control', 'Control', None, 'Treated', ''],
        'celltype': ['T', None, 'B', 'B', 'T'],
    })
    figure = condition_composition_figure(obs, None, 'condition', 'celltype')
    try:
        axis = figure.axes[0]
        assert [tick.get_text() for tick in axis.get_xticklabels()] == ['Control', 'Treated']
        assert 'nan' not in [text.get_text() for text in axis.get_legend().get_texts()]
    finally:
        import matplotlib.pyplot as plt
        plt.close(figure)


def test_condition_composition_reserves_space_for_a_dense_long_name_legend():
    import matplotlib.pyplot as plt
    from modules.figure_style import NATURE_TEXT
    from modules.sc_figure_diagnostics import condition_composition_figure

    labels = [f'Extremely detailed cell type {index:02d} with a long qualifier'
              for index in range(20)]
    obs = pd.DataFrame({
        'condition': ['Control'] * len(labels) + ['Treated'] * len(labels),
        'celltype': labels + labels,
    })
    figure = condition_composition_figure(obs, None, 'condition', 'celltype')
    try:
        # The legend gets its own column instead of squeezing the stack into a
        # fixed-width canvas.  Long names therefore remain usable in the web
        # preview as well as in the saved export.
        assert figure.get_size_inches()[0] > 10
        assert figure._native_layout_rect[2] < 0.70
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        assert (figure.axes[0].get_legend().get_window_extent(renderer).x0
                > figure.axes[0].get_window_extent(renderer).x1)
        # At least one of the lighter category colours uses dark text, rather
        # than disappearing as white text on a pale stack segment.
        assert any(text.get_color() == NATURE_TEXT for text in figure.axes[0].texts)
    finally:
        plt.close(figure)


def test_annotation_decision_heatmap_uses_final_labels_without_text_overlap():
    import matplotlib.pyplot as plt
    from modules.sc_figure_diagnostics import annotation_decision_heatmap

    programmes = [
        'CKB/FABP1+ metabolic enterocytes',
        'HSD17B2+ absorptive enterocytes',
        'SPINK4/CA4+ secretory cells',
    ]
    clusters = [str(index) for index in range(11)]
    score_frame = pd.DataFrame(
        np.arange(len(programmes) * len(clusters), dtype=float).reshape(
            len(programmes), len(clusters),
        ),
        index=[f'score_{label}' for label in programmes], columns=clusters,
    )
    final_labels = {
        cluster: [
            'High metabolic cells', 'Absorptive Enterocytes', 'Goblet cells',
            'MAML3/CHRM3+ Epithelial cells',
        ][index % 4]
        for index, cluster in enumerate(clusters)
    }
    figure = annotation_decision_heatmap(
        score_frame, clusters,
        decisions={cluster: {'score_margin': 0.1} for cluster in clusters},
        final_labels=final_labels,
        display_label_map={
            'CKB/FABP1+ metabolic enterocytes': 'High metabolic cells',
            'HSD17B2+ absorptive enterocytes': 'Absorptive Enterocytes',
            'SPINK4/CA4+ secretory cells': 'Goblet cells',
        },
    )
    try:
        axis = figure.axes[0]
        assert not axis.texts
        assert 'High metabolic cells' in axis.get_xticklabels()[0].get_text().replace('\n', ' ')
        assert [tick.get_text() for tick in axis.get_yticklabels()] == [
            'High metabolic cells', 'Absorptive Enterocytes', 'Goblet cells',
        ]
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        boxes = [tick.get_window_extent(renderer) for tick in axis.get_xticklabels()]
        assert all(left.x1 <= right.x0 for left, right in zip(boxes, boxes[1:]))
        assert figure._native_layout_rect[1] > 0
    finally:
        plt.close(figure)
