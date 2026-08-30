"""Deterministic director for Nature figure specifications and renderers.

The director translates user-controlled scientific parameters into a
``FigureSpec``.  Visual constants remain owned by the style profiles and plot
templates, so callers cannot accidentally produce a publication figure with
arbitrary fonts, colours, or line widths.
"""

from __future__ import annotations

import math
from dataclasses import fields
from typing import Any, Mapping

from .spec import FigureSpec
from .templates import (
    NatureEnrichmentDotplot,
    NatureEnrichmentOverview,
    NatureEnrichmentBarplot,
    NatureEnrichmentChord,
    NatureEnrichmentCnetplot,
    NatureEnrichmentEmapplot,
    NatureCorrelationHeatmap,
    NatureDiagnostic,
    NatureEmbedding,
    NatureGSEA,
    NatureGSEARunning,
    NatureGSVA,
    NatureHeatmap,
    NatureMA,
    NaturePCA,
    NatureSSGSEA,
    NatureUpSet,
    NatureVolcano,
    NatureWGCNA,
)


_PLOT_ALIASES = {
    "pca": "pca",
    "volcano": "volcano",
    "heatmap": "heatmap",
    "gsea": "gsea",
    "ora": "enrichment_dotplot",
    "go": "enrichment_dotplot",
    "kegg": "enrichment_dotplot",
    "enrichment": "enrichment_dotplot",
    "enrichment_dotplot": "enrichment_dotplot",
    "enrichment_overview": "enrichment_overview",
    "multi_database_overview": "enrichment_overview",
    "pathway_overview": "enrichment_overview",
    "barplot": "enrichment_barplot",
    "enrichment_barplot": "enrichment_barplot",
    "chord": "enrichment_chord",
    "gochord": "enrichment_chord",
    "enrichment_chord": "enrichment_chord",
    "cnetplot": "enrichment_cnetplot",
    "enrichment_cnetplot": "enrichment_cnetplot",
    "emapplot": "enrichment_emapplot",
    "enrichment_emapplot": "enrichment_emapplot",
    "gsea_running": "gsea_running",
    "gsea_running_score": "gsea_running",
    "ma": "ma",
    "ma_plot": "ma",
    "correlation": "correlation_heatmap",
    "correlation_heatmap": "correlation_heatmap",
    "gsva": "gsva",
    "ssgsea": "ssgsea",
    "ssgsea_heatmap": "ssgsea",
    "upset": "upset",
    "wgcna": "wgcna",
    "module_trait": "wgcna",
    "embedding": "embedding",
    "umap": "embedding",
    "tsne": "embedding",
    # Bulk and single-cell diagnostics use one deterministic renderer. The
    # aliases remain explicit so callers can keep scientific plot names while
    # visual constants stay in the engine.
    "diagnostic": "diagnostic",
    "normalization": "diagnostic",
    "normalization_library": "diagnostic",
    "normalization_boxplot": "diagnostic",
    "normalization_pca": "diagnostic",
    "qc_overview": "diagnostic",
    "qc_pairs": "diagnostic",
    "qc_violin": "diagnostic",
    "variance": "diagnostic",
    "qq": "diagnostic",
    "trajectory": "diagnostic",
    "timecourse": "diagnostic",
    "boxplot": "diagnostic",
}


class NatureFigureDirector:
    """Choose a fixed template and construct a validated ``FigureSpec``."""

    def __init__(self) -> None:
        self._renderers = {
            "pca": NaturePCA(),
            "volcano": NatureVolcano(),
            "heatmap": NatureHeatmap(),
            "gsea": NatureGSEA(),
            "enrichment_dotplot": NatureEnrichmentDotplot(),
            "enrichment_overview": NatureEnrichmentOverview(),
            "enrichment_barplot": NatureEnrichmentBarplot(),
            "enrichment_chord": NatureEnrichmentChord(),
            "enrichment_cnetplot": NatureEnrichmentCnetplot(),
            "enrichment_emapplot": NatureEnrichmentEmapplot(),
            "gsea_running": NatureGSEARunning(),
            "ma": NatureMA(),
            "correlation_heatmap": NatureCorrelationHeatmap(),
            "gsva": NatureGSVA(),
            "ssgsea": NatureSSGSEA(),
            "upset": NatureUpSet(),
            "wgcna": NatureWGCNA(),
            "embedding": NatureEmbedding(),
            "diagnostic": NatureDiagnostic(),
        }

    @staticmethod
    def normalize_plot_type(plot_type: str) -> str:
        key = str(plot_type).strip().lower()
        try:
            return _PLOT_ALIASES[key]
        except KeyError as exc:
            supported = ", ".join(sorted(_PLOT_ALIASES))
            raise ValueError(
                f"Unsupported plot_type {plot_type!r}; expected one of: {supported}"
            ) from exc

    def create_spec(self, plot_type: str, **options: Any) -> FigureSpec:
        """Create a spec while rejecting unknown, visually arbitrary options."""

        normalized = self.normalize_plot_type(plot_type)
        allowed = {field.name for field in fields(FigureSpec)} - {"plot_type"}
        unknown = sorted(set(options) - allowed)
        if unknown:
            raise TypeError(
                "Unknown FigureSpec option(s): " + ", ".join(unknown)
            )
        return FigureSpec(plot_type=normalized, **options)

    def spec_from_params(
        self,
        plot_type: str,
        params: Mapping[str, Any] | None,
        **overrides: Any,
    ) -> FigureSpec:
        """Map platform parameters to the small scientific-control surface.

        ``params`` may contain the legacy ``_visualization`` dictionary.  Only
        semantically meaningful controls are accepted here; visual constants
        stay fixed in the selected style profile.
        """

        source = dict(params or {})
        visualization = dict(source.get("_visualization") or {})
        normalized = self.normalize_plot_type(plot_type)

        mode = (
            source.get("_figure_mode")
            or visualization.get("figure_mode")
            or (
                "nature_portfolio"
                if visualization.get("theme") == "nature"
                else "publication"
            )
        )
        width = (
            source.get("_figure_width")
            or visualization.get("width_profile")
            or ("double" if normalized in {"heatmap", "enrichment_overview"} else "single")
        )
        style = source.get("_nature_style") or visualization.get("nature_style") or "nature"

        options: dict[str, Any] = {
            "mode": mode,
            "width": width,
            "style": style,
        }
        requested_formats = visualization.get(
            "static_formats", visualization.get("export_formats"))
        if requested_formats:
            if isinstance(requested_formats, str):
                requested_formats = requested_formats.replace(',', ' ').split()
            options["formats"] = tuple(requested_formats)

        if normalized == "pca":
            options.update(
                confidence_ellipse=_as_bool(source.get("confidence_ellipse", False)),
                show_sample_labels=_as_bool(source.get("show_sample_labels", False)),
                outlier_labels=_as_tuple(source.get("outlier_labels", ())),
            )
        elif normalized == "volcano":
            if "log2fc_threshold" in source:
                fc_threshold = float(source["log2fc_threshold"])
            else:
                legacy_fc = max(float(source.get("fc_threshold", 2.0)), 1e-12)
                fc_threshold = math.log2(legacy_fc)
            options.update(
                fc_threshold=fc_threshold,
                fdr_threshold=float(source.get(
                    "fdr_threshold", source.get("pval_threshold", 0.05)
                )),
                label_n=int(source.get("label_n", source.get("top_n", 8))),
                label_strategy=str(source.get("label_strategy", "top_significant")),
                label_genes=_as_tuple(source.get("label_genes", ())),
            )
        elif normalized == "heatmap":
            options.update(
                zscore=str(source.get("zscore", "row")),
                row_cluster=_as_bool(source.get("row_cluster", True)),
                col_cluster=_as_bool(source.get("col_cluster", True)),
                cluster_method=str(source.get("cluster_method", "average")),
                distance_metric=str(
                    source.get("distance_metric", source.get("cluster_metric", "correlation"))
                ),
                max_row_labels=int(source.get("max_row_labels", 35)),
                max_col_labels=int(source.get("max_col_labels", 24)),
            )
        elif normalized in {
            "gsea", "enrichment_dotplot", "enrichment_overview", "enrichment_barplot",
            "enrichment_chord", "enrichment_cnetplot",
            "enrichment_emapplot", "gsea_running",
        }:
            options.update(
                top_n=int(source.get("pathway_number", source.get("top_n", 12))),
                term_wrap=int(source.get("term_wrap", source.get("term_width", 34))),
                term_max_chars=int(source.get("term_max_chars", 82)),
                max_genes=int(source.get("enrichment_max_genes", source.get("max_genes", 30))),
                similarity_threshold=float(source.get("similarity_threshold", 0.15)),
                redundancy_threshold=float(source.get("redundancy_threshold", 0.85)),
                running_term_n=int(source.get("running_term_n", 2)),
                pathway_selection=str(source.get("pathway_selection", "top")),
                pathway_terms=_as_tuple(source.get("target_pathways", source.get("pathway_terms", ()))),
                database_scope=_as_tuple(source.get("database_scope", source.get("databases", ()))),
                gene_label_strategy=str(source.get("gene_label_strategy", "all")),
                gene_labels=_as_tuple(source.get("target_genes", source.get("gene_labels", ()))),
                max_gene_labels=int(source.get("max_gene_labels", 30)),
            )
        elif normalized == "ma":
            options.update(
                fc_threshold=float(source.get("log2fc_threshold", 1.0)),
                fdr_threshold=float(source.get("fdr_threshold", source.get("pval_threshold", 0.05))),
                label_n=int(source.get("label_n", 6)),
                label_strategy=str(source.get("label_strategy", "top_significant")),
                label_genes=_as_tuple(source.get("label_genes", ())),
            )
        elif normalized == "correlation_heatmap":
            options.update(
                correlation_method=str(source.get("corr_method", source.get("correlation_method", "pearson"))),
                col_cluster=_as_bool(source.get("cluster", True)),
                mask_diagonal=_as_bool(source.get("mask_diagonal", True)),
                annotate_cells=_as_bool(source.get("annotate_cells", False)),
                max_col_labels=int(source.get("max_col_labels", 18)),
            )
        elif normalized in {"gsva", "ssgsea"}:
            options.update(
                score_scale=str(source.get("score_scale", "row")),
                row_cluster=_as_bool(source.get("row_cluster", True)),
                col_cluster=_as_bool(source.get("col_cluster", True)),
                max_row_labels=int(source.get("max_row_labels", 30)),
                max_col_labels=int(source.get("max_col_labels", 18)),
            )
        elif normalized == "upset":
            options.update(
                top_intersections=int(source.get("upset_top_n", source.get("top_intersections", 20))),
            )
        elif normalized == "wgcna":
            options.update(
                annotate_cells=_as_bool(source.get("annotate_cells", True)),
                max_row_labels=int(source.get("max_row_labels", 30)),
                max_col_labels=int(source.get("max_col_labels", 16)),
            )

        options.update(overrides)
        return self.create_spec(normalized, **options)

    def render(self, spec: FigureSpec, data: Any):
        """Render data with the single template registered for ``spec``."""

        normalized = self.normalize_plot_type(spec.plot_type)
        if normalized != spec.plot_type:
            spec = spec.with_updates(plot_type=normalized)
        renderer = self._renderers[normalized]
        return renderer.render(data, spec)

    def render_into(self, spec: FigureSpec, data: Any, container):
        """Render an editable panel into a Matplotlib SubFigure."""

        normalized = self.normalize_plot_type(spec.plot_type)
        if normalized != spec.plot_type:
            spec = spec.with_updates(plot_type=normalized)
        return self._renderers[normalized].render(data, spec, container=container)

    def create_and_render(
        self,
        plot_type: str,
        data: Any,
        params: Mapping[str, Any] | None = None,
        **overrides: Any,
    ):
        spec = self.spec_from_params(plot_type, params, **overrides)
        return spec, self.render(spec, data)


def _as_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _as_tuple(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return tuple(item.strip() for item in value.replace("\n", ",").replace(";", ",").split(",") if item.strip())
    return tuple(str(item).strip() for item in value if str(item).strip())
