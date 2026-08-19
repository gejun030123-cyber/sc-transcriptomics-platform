"""GSVA/ssGSEA pathway-score heatmaps with fixed publication semantics."""

from figure_engine.spec import FigureSpec
from .heatmap import NatureHeatmap
from .common import wrap_term


class _NatureGeneSetHeatmap:
    plot_type = 'gsva'
    default_title = 'Gene-set variation analysis'

    def render(self, data, spec: FigureSpec, container=None):
        payload = dict(data or {})
        matrix = payload.get('matrix', payload.get('scores'))
        if matrix is None:
            raise ValueError(f'{self.__class__.__name__} 需要 pathway x sample score matrix')
        payload['matrix'] = matrix
        pathway_labels = payload.get(
            'pathway_labels', payload.get('gene_labels', payload.get('pathways')))
        if pathway_labels is not None:
            pathway_labels = [wrap_term(value, width=22, max_chars=54)
                              for value in pathway_labels]
        payload['gene_labels'] = pathway_labels
        payload['sample_labels'] = payload.get('sample_labels', payload.get('samples'))
        payload['annotations'] = payload.get('annotations', {})
        payload['colorbar_label'] = (
            'Pathway-wise z-score' if spec.score_scale == 'row' else 'Enrichment score'
        )
        heatmap_spec = spec.with_updates(
            plot_type='heatmap',
            title=spec.title or self.default_title,
            zscore='row' if spec.score_scale == 'row' else 'precomputed',
            max_row_labels=min(spec.max_row_labels, 30),
        )
        figure = NatureHeatmap().render(payload, heatmap_spec, container=container)
        if container is None:
            figure._nature_spec = spec.to_dict()
            figure._nature_spec_object = spec
            figure._nature_encodings = {
                'color': payload['colorbar_label'].lower(),
                'annotation': ', '.join(payload['annotations']),
            }
        return figure


class NatureGSVA(_NatureGeneSetHeatmap):
    plot_type = 'gsva'
    default_title = 'Gene-set variation analysis (GSVA)'


class NatureSSGSEA(_NatureGeneSetHeatmap):
    plot_type = 'ssgsea'
    default_title = 'Single-sample GSEA'
