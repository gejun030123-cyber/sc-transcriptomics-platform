"""Auditable, theme-prioritised GO BP/CC/MF enrichment views.

These templates intentionally render the three GO ontologies in a *single*
axis split into semantic sections.  They are used after ``sc_cell_go`` has
selected a fixed number of already FDR-significant rows per ontology.  The
renderer does not re-rank, re-test, or apply a theme filter: selection reason
and rank remain source-table fields produced by the analysis module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from figure_engine.spec import FigureSpec
from figure_engine.style import get_style
from .common import attach_contract, strip_term_id, wrap_term


GO_SECTIONS = (
    ("GO_BP", "Biological Process", "#168CA0"),
    ("GO_CC", "Cellular Component", "#E99816"),
    ("GO_MF", "Molecular Function", "#8C3196"),
)

# Three independently tested Human pathway libraries can use the same
# single-axis, divided-layout contract as GO.  The internal keys are
# canonicalised by ``_database_key`` so historical spelling/case variants in
# result tables do not make a section silently disappear.
PATHWAY_SECTIONS = (
    ("KEGG", "KEGG Human", "#168CA0"),
    ("REACTOME", "Reactome", "#E99816"),
    ("WIKIPATHWAYS", "WikiPathways Human", "#8C3196"),
)


def _first_column(frame, candidates):
    lookup = {str(column).strip().lower(): column for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def _database_key(value):
    text = str(value or "").upper().replace("-", "_")
    aliases = {
        "GO_BP": "GO_BP", "GO_BIOLOGICAL_PROCESS": "GO_BP", "BP": "GO_BP",
        "GO_CC": "GO_CC", "GO_CELLULAR_COMPONENT": "GO_CC", "CC": "GO_CC",
        "GO_MF": "GO_MF", "GO_MOLECULAR_FUNCTION": "GO_MF", "MF": "GO_MF",
    }
    return aliases.get(text, text)


def _ratio_value(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    if isinstance(value, str) and "/" in value:
        left, right = value.split("/", 1)
        try:
            denominator = float(right)
            return float(left) / denominator if denominator else np.nan
        except ValueError:
            return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _count_value(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    if isinstance(value, str) and "/" in value:
        value = value.split("/", 1)[0]
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _bubble_sizes(values, low=26.0, high=124.0):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return np.asarray([], dtype=float)
    finite = values[np.isfinite(values)]
    if not len(finite) or np.isclose(finite.min(), finite.max()):
        return np.full(len(values), (low + high) / 2.0)
    root = np.sqrt(np.clip(values, 0, None))
    scaled = (root - np.nanmin(root)) / (np.nanmax(root) - np.nanmin(root))
    return low + scaled * (high - low)


def _format_qvalue(value):
    value = float(value)
    return f"{value:.1e}" if value < .01 else f"{value:.3f}"


def _priority_frame(data, *, sections, view_label):
    """Normalise the small, explicit triptych source-table contract.

    A fully empty source remains valid so that the plot can truthfully retain
    all three ontology sections.  Non-empty rows require both term and FDR.
    """
    frame = pd.DataFrame(data).copy()
    if frame.empty and not len(frame.columns):
        frame = pd.DataFrame(columns=["Database", "Term", "Adjusted P-value", "Overlap"])
    database_column = _first_column(frame, ("Database", "database", "gene_set", "Gene_set"))
    term_column = _first_column(frame, ("Term", "Description", "pathway", "term"))
    fdr_column = _first_column(frame, (
        "Adjusted P-value", "Adjusted p-value", "Enrichment FDR", "FDR", "fdr", "padj", "qvalue",
    ))
    ratio_column = _first_column(frame, ("GeneRatio", "gene_ratio", "fraction", "Overlap"))
    count_column = _first_column(frame, ("Count", "count", "Gene Count", "gene_count", "num", "Overlap"))
    if database_column is None:
        raise ValueError(f"{view_label} 三分区图需要 Database/gene_set 字段")
    if term_column is None:
        raise ValueError(f"{view_label} 三分区图需要 Term/Description 字段")
    if fdr_column is None:
        raise ValueError(f"{view_label} 三分区图需要 FDR/Adjusted P-value 字段")
    if ratio_column is None:
        raise ValueError(f"{view_label} 气泡图需要 GeneRatio 或 Overlap 字段")
    if count_column is None:
        raise ValueError(f"{view_label} 三分区图需要 Count 或 Overlap 字段")

    frame["_database"] = frame[database_column].map(_database_key)
    frame["_term"] = frame[term_column].map(str)
    frame["_display_term"] = frame["_term"].map(strip_term_id)
    frame["_fdr"] = pd.to_numeric(frame[fdr_column], errors="coerce")
    frame["_ratio"] = frame[ratio_column].map(_ratio_value)
    frame["_count"] = frame[count_column].map(_count_value)
    rank_column = _first_column(frame, ("selection_rank", "Selection rank", "display_rank"))
    frame["_selection_rank"] = (
        pd.to_numeric(frame[rank_column], errors="coerce") if rank_column else np.nan
    )
    frame = frame.loc[frame["_database"].isin({item[0] for item in sections})].copy()
    frame = frame.replace([np.inf, -np.inf], np.nan)
    # The source selection has already enforced FDR significance.  This only
    # removes malformed rows rather than silently relaxing the result rule.
    return frame.dropna(subset=["_fdr", "_ratio", "_count"]).copy()


def _go_priority_frame(data):
    """Compatibility wrapper for the established GO triptych contract."""
    return _priority_frame(data, sections=GO_SECTIONS, view_label="GO")


def _pathway_triptych_frame(data):
    """Normalise KEGG/Reactome/WikiPathways Human triptych source rows."""
    return _priority_frame(data, sections=PATHWAY_SECTIONS, view_label="Human pathway")


def _section_layout(frame, *, sections, top_n):
    """Create deterministic y positions, section centres, and separators."""
    rows = []
    layout_sections = []
    boundaries = []
    cursor = 0.0
    for index, (database, label, color) in enumerate(sections):
        subset = frame.loc[frame["_database"].eq(database)].copy()
        subset = subset.sort_values(
            ["_selection_rank", "_fdr", "_ratio", "_count", "_display_term"],
            ascending=[True, True, False, False, True], na_position="last",
            kind="mergesort",
        ).head(top_n).reset_index(drop=True)
        if subset.empty:
            centre = cursor
            layout_sections.append({
                "database": database, "label": label, "color": color,
                "centre": centre, "empty": True,
            })
            cursor += 1.0
        else:
            # ``np.arange(start, start + n)`` can include an extra endpoint
            # for decimal starts such as 3.8.  Build the integer length first
            # so every selected term receives exactly one y coordinate.
            positions = cursor + np.arange(len(subset), dtype=float)
            subset["_y"] = positions
            rows.append(subset)
            layout_sections.append({
                "database": database, "label": label, "color": color,
                "centre": float(np.mean(positions)), "empty": False,
            })
            cursor += float(len(subset))
        if index < len(sections) - 1:
            boundaries.append(cursor - .5)
            cursor += .8
    display = pd.concat(rows, ignore_index=True) if rows else frame.iloc[0:0].copy()
    return display, layout_sections, boundaries, max(cursor - .3, 1.0)


def _legend_axis(fig, *, qvalues, counts, cmap, norm, style, profile):
    """Draw the reference-style combined q-value and Count legend."""
    ax = fig.add_axes([.79, .58, .18, .30])
    ax._nature_auxiliary = True
    ax.set_axis_off()
    ax.text(.5, .98, "Legend", ha="center", va="top", fontsize=profile.legend_font_pt,
            color=style.text)
    ax.text(.50, .84, "qvalue", ha="center", va="center", fontsize=profile.legend_font_pt,
            color=style.text)
    qvalues = np.asarray(qvalues, dtype=float)
    for ypos, value in zip(np.linspace(.70, .46, len(qvalues)), qvalues):
        ax.scatter(.25, ypos, s=24, color=cmap(norm(value)), edgecolors="none", clip_on=False)
        ax.text(.43, ypos, _format_qvalue(value), ha="left", va="center",
                fontsize=max(5.2, profile.legend_font_pt - .3), color=style.text)
    ax.text(.50, .33, "Count", ha="center", va="center", fontsize=profile.legend_font_pt,
            color=style.text)
    counts = np.asarray(counts, dtype=float)
    for ypos, value, size in zip(np.linspace(.22, .03, len(counts)), counts, _bubble_sizes(counts)):
        ax.scatter(.25, ypos, s=size, color="#333333", edgecolors="none", clip_on=False)
        ax.text(.43, ypos, str(int(round(value))), ha="left", va="center",
                fontsize=max(5.2, profile.legend_font_pt - .3), color=style.text)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return ax


class NatureGOPriorityDotplot:
    """Reference-style GO dot plot: one x axis and three labelled sections."""

    plot_type = "go_priority_dotplot"
    sections = GO_SECTIONS
    view_label = "GO"
    y_label = "GO Description"
    default_title = "GO Enrichment Plot"
    semantic_warning = (
        "每个 GO 本体按固定名额依次展示炎症、脂代谢、其他确认主题与 FDR Top 通路；"
        "统计和 FDR 均来自完整 GO 本体。",
    )

    def render(self, data, spec: FigureSpec, container=None):
        if container is not None:
            raise ValueError(f"{self.view_label} 三分区气泡图需要独立 Figure，不支持嵌入单一 panel")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap, Normalize

        style = get_style(spec.style)
        frame = _priority_frame(data, sections=self.sections, view_label=self.view_label)
        display, sections, boundaries, y_max = _section_layout(
            frame, sections=self.sections, top_n=spec.top_n,
        )
        with style.context(spec):
            profile = style.profile(spec.with_updates(plot_type=self.plot_type))
            fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
            cmap = LinearSegmentedColormap.from_list(
                "go_qvalue_priority", ("#FF120C", "#CE0052", "#7600A8", "#1F00EE"), N=256,
            )
            if display.empty:
                norm = Normalize(vmin=.001, vmax=.05, clip=True)
                qvalues = np.asarray([.001, .01, .025, .05])
                counts = np.asarray([1, 2, 3])
            else:
                q_min, q_max = float(display["_fdr"].min()), float(display["_fdr"].max())
                if np.isclose(q_min, q_max):
                    q_min = max(np.finfo(float).tiny, q_min * .5)
                    q_max = min(1.0, q_max * 1.5)
                norm = Normalize(vmin=q_min, vmax=q_max, clip=True)
                qvalues = np.linspace(q_min, q_max, 4)
                count_min, count_max = float(display["_count"].min()), float(display["_count"].max())
                counts = np.unique(np.round(np.linspace(count_min, count_max, 3))).astype(float)
                if len(counts) == 1:
                    counts = np.asarray([counts[0], counts[0] + 1, counts[0] + 2])

            if not display.empty:
                ax.scatter(
                    display["_ratio"], display["_y"], s=_bubble_sizes(display["_count"]),
                    c=display["_fdr"], cmap=cmap, norm=norm, edgecolors="white", linewidths=.45,
                    alpha=.98, zorder=3,
                )
                labels = [wrap_term(value, width=40, max_chars=90) for value in display["_display_term"]]
                ax.set_yticks(display["_y"], labels)
                ratio_max = max(float(display["_ratio"].max()), .01)
                ax.set_xlim(0, ratio_max * 1.16)
            else:
                ax.set_xlim(0, 1)
                ax.set_yticks([])

            for boundary in boundaries:
                ax.axhline(boundary, color=style.text, linewidth=.75, zorder=1)
            for section in sections:
                if section["empty"]:
                    ax.text(.50, section["centre"], "No FDR-significant pathways",
                            transform=ax.get_yaxis_transform(), ha="center", va="center",
                            color=style.muted_text, fontsize=profile.axis_font_pt)
                ax.text(1.08, section["centre"], section["label"], rotation=90,
                        transform=ax.get_yaxis_transform(), ha="center", va="center",
                        color=style.text, fontsize=profile.axis_font_pt, clip_on=False)
            ax.set_ylim(y_max, -.7)
            ax.set_xlabel("GeneRatio")
            ax.set_ylabel(self.y_label)
            ax.set_title(spec.title or self.default_title, pad=7)
            ax.tick_params(axis="y", length=0, pad=3, labelsize=max(5.4, profile.tick_font_pt - .1))
            style.apply_axis(ax, profile)
            ax.spines["top"].set_visible(True)
            ax.spines["right"].set_visible(True)
            ax.margins(x=0)
            fig.subplots_adjust(left=.43, right=.76, top=.91, bottom=.10)
            _legend_axis(fig, qvalues=qvalues, counts=counts, cmap=cmap, norm=norm,
                         style=style, profile=profile)

        fig._nature_panel_grid = {"nrows": 1, "ncols": 1, "panels": 3}
        return attach_contract(
            fig, spec.with_updates(plot_type=self.plot_type), style,
            semantic_warnings=self.semantic_warning,
            encodings={"x": "GeneRatio", "color": "FDR q value", "size": "overlap gene Count"},
        )


class NatureGOPriorityBarplot:
    """Reference-style GO bar plot with category-coloured terms and counts."""

    plot_type = "go_priority_barplot"
    sections = GO_SECTIONS
    view_label = "GO"
    y_label = "GO Description"
    default_title = "GO enrichment"
    semantic_warning = (
        "每个 GO 本体按固定名额依次展示炎症、脂代谢、其他确认主题与 FDR Top 通路；"
        "柱长和柱尾数字为重叠基因 Count，统计和 FDR 均来自完整 GO 本体。",
    )

    def render(self, data, spec: FigureSpec, container=None):
        if container is not None:
            raise ValueError(f"{self.view_label} 三分区柱状图需要独立 Figure，不支持嵌入单一 panel")
        import matplotlib.pyplot as plt

        style = get_style(spec.style)
        frame = _priority_frame(data, sections=self.sections, view_label=self.view_label)
        display, sections, boundaries, y_max = _section_layout(
            frame, sections=self.sections, top_n=spec.top_n,
        )
        color_by_database = {database: color for database, _label, color in self.sections}
        with style.context(spec):
            profile = style.profile(spec.with_updates(plot_type=self.plot_type))
            fig, ax = plt.subplots(figsize=profile.figsize, dpi=profile.dpi)
            if not display.empty:
                colors = [color_by_database.get(value, style.neutral_dark) for value in display["_database"]]
                ax.barh(display["_y"], display["_count"], color=colors, edgecolor="none",
                        height=.72, alpha=.98, zorder=2)
                labels = [wrap_term(value, width=44, max_chars=92) for value in display["_display_term"]]
                ax.set_yticks(display["_y"], labels)
                for tick, database in zip(ax.get_yticklabels(), display["_database"]):
                    tick.set_color(color_by_database.get(database, style.text))
                count_max = max(float(display["_count"].max()), 1.0)
                ax.set_xlim(0, count_max * 1.16)
                for y, count in zip(display["_y"], display["_count"]):
                    ax.text(count + count_max * .018, y, str(int(round(count))), ha="left", va="center",
                            color=style.text, fontsize=max(5.3, profile.legend_font_pt - .2))
            else:
                ax.set_xlim(0, 1)
                ax.set_yticks([])

            for boundary in boundaries:
                ax.axhline(boundary, color=style.text, linewidth=.75, zorder=1)
            for section in sections:
                if section["empty"]:
                    ax.text(.50, section["centre"], "No FDR-significant pathways",
                            transform=ax.get_yaxis_transform(), ha="center", va="center",
                            color=style.muted_text, fontsize=profile.axis_font_pt)
                ax.text(1.08, section["centre"], section["label"], rotation=90,
                        transform=ax.get_yaxis_transform(), ha="center", va="center",
                        color=section["color"], fontsize=profile.axis_font_pt, clip_on=False)
            ax.set_ylim(y_max, -.7)
            ax.set_xlabel("Gene Count")
            ax.set_ylabel(self.y_label)
            ax.set_title(spec.title or self.default_title, pad=7)
            ax.tick_params(axis="y", length=0, pad=3, labelsize=max(5.2, profile.tick_font_pt - .2))
            style.apply_axis(ax, profile)
            ax.spines["top"].set_visible(True)
            ax.spines["right"].set_visible(True)
            ax.margins(x=0)
            fig.subplots_adjust(left=.43, right=.76, top=.91, bottom=.10)

        fig._nature_panel_grid = {"nrows": 1, "ncols": 1, "panels": 3}
        return attach_contract(
            fig, spec.with_updates(plot_type=self.plot_type), style,
            semantic_warnings=self.semantic_warning,
            encodings={"x": "overlap gene Count", "color": f"{self.view_label} section", "label": "overlap gene Count"},
        )


class NaturePathwayTriptychDotplot(NatureGOPriorityDotplot):
    """KEGG/Reactome/WikiPathways Human dot plot on one divided axis."""

    plot_type = "pathway_triptych_dotplot"
    sections = PATHWAY_SECTIONS
    view_label = "Human pathway"
    y_label = "Pathway"
    default_title = "Human Pathway Enrichment Plot"
    semantic_warning = (
        "KEGG、Reactome 与 WikiPathways Human 各自保留完整库内独立 FDR；"
        "本图仅将各库已显著的 FDR Top 通路并列展示，不合并或重新校正 FDR。",
    )


class NaturePathwayTriptychBarplot(NatureGOPriorityBarplot):
    """KEGG/Reactome/WikiPathways Human count bar plot on one divided axis."""

    plot_type = "pathway_triptych_barplot"
    sections = PATHWAY_SECTIONS
    view_label = "Human pathway"
    y_label = "Pathway"
    default_title = "Human pathway enrichment"
    semantic_warning = (
        "KEGG、Reactome 与 WikiPathways Human 各自保留完整库内独立 FDR；"
        "柱长和柱尾数字为重叠基因 Count，本图不合并或重新校正 FDR。",
    )
