"""
Locus Viewer - Gene visualization component for RADS dashboard.
Creates interactive gene arrow diagrams similar to geneviewer.
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import polars as pl
from typing import Optional


# Teal/slate color palette inspired by defense system biology figures
# Primary palette: dark teal to light sage
TEAL_PALETTE = {
    "darkest": "#2d4a4a",      # Very dark teal
    "dark": "#3d5a5a",         # Dark teal
    "medium_dark": "#4a6670",  # Slate blue-teal
    "medium": "#5a7a7a",       # Medium teal
    "medium_light": "#6b9090", # Medium-light teal
    "light": "#7a9e9e",        # Light teal
    "lighter": "#8fb3b3",      # Lighter teal/sage
    "lightest": "#a8c4c4",     # Very light sage
    "pale": "#b5cece",         # Pale sage
    "accent": "#4a6670",       # Accent slate
}

# Color palettes for different annotation sources
DOMAIN_COLORS = {
    # InterProScan analysis types - teal/slate gradient
    "Pfam": "#2d4a4a",         # Darkest teal
    "TIGRFAM": "#3d5a5a",      # Dark teal
    "SMART": "#4a6670",        # Slate blue-teal
    "CDD": "#5a7a7a",          # Medium teal
    "Gene3D": "#6b9090",       # Medium-light teal
    "SUPERFAMILY": "#7a9e9e",  # Light teal
    "ProSiteProfiles": "#8fb3b3", # Lighter teal
    "PANTHER": "#a8c4c4",      # Light sage
    "Hamap": "#b5cece",        # Pale sage
    "default_interproscan": "#8a9a9a",
}

DEFENSE_SYSTEM_COLORS = {
    # DefenseFinder system types - distinctive colors for visibility
    # Using more saturated/contrasting colors to distinguish from other genes
    "RM": "#c0392b",           # Restriction-Modification (red)
    "CRISPR": "#2980b9",       # CRISPR-Cas (blue)
    "Abi": "#8e44ad",          # Abortive infection (purple)
    "TA": "#d35400",           # Toxin-Antitoxin (orange)
    "BREX": "#27ae60",         # BREX (green)
    "DISARM": "#16a085",       # DISARM (teal)
    "Gabija": "#e74c3c",       # Gabija (light red)
    "Hachiman": "#3498db",     # Hachiman (light blue)
    "Lamassu": "#9b59b6",      # Lamassu (light purple)
    "Lamassu-Fam": "#9b59b6",  # Lamassu family (light purple)
    "Thoeris": "#f39c12",      # Thoeris (yellow-orange)
    "Zorya": "#1abc9c",        # Zorya (turquoise)
    "Druantia": "#e67e22",     # Druantia (carrot orange)
    "Kiwa": "#2ecc71",         # Kiwa (emerald)
    "Wadjet": "#34495e",       # Wadjet (wet asphalt)
    "Septu": "#95a5a6",        # Septu (concrete)
    "RosmerTA": "#c0392b",     # RosmerTA (red - TA system)
    "MazEF": "#d35400",        # MazEF (orange - TA system)
    "PD-Lambda-1": "#8e44ad",  # PD-Lambda (purple)
    "Dodola": "#27ae60",       # Dodola (green)
    "AbiH": "#2980b9",         # AbiH (blue)
    "AbiC": "#3498db",         # AbiC (light blue)
    "AbiJ": "#1abc9c",         # AbiJ (turquoise)
    "AbiE": "#16a085",         # AbiE (teal)
    "PrrC": "#e74c3c",         # PrrC (light red)
    "RloC": "#f39c12",         # RloC (yellow-orange)
    "default_defense": "#e74c3c",  # Default: light red for visibility
}

# Special colors - teal theme
QUERY_HIT_COLOR = "#2d4a4a"      # Dark teal for query hits (recombinases)
DOWNSTREAM_COLOR = "#6b9090"     # Medium teal for co-transcribed downstream genes
NO_ANNOTATION_COLOR = "#d5e0e0"  # Light gray-teal for unannotated genes


def create_gene_arrow(
    x_start: float,
    x_end: float,
    y_center: float,
    strand: int,
    color: str,
    label: str = "",
    hover_text: str = "",
    arrow_height: float = 0.4,
    arrow_head_width: float = 0.15,
) -> dict:
    """
    Create a gene arrow shape for Plotly.

    Args:
        x_start: Start position of gene
        x_end: End position of gene
        y_center: Y coordinate for center of arrow
        strand: 1 for forward, -1 for reverse
        color: Fill color for the arrow
        label: Gene label
        hover_text: Text to display on hover
        arrow_height: Height of arrow body
        arrow_head_width: Relative width of arrow head

    Returns:
        Dictionary with shape and annotation data
    """
    gene_length = abs(x_end - x_start)
    head_length = min(gene_length * arrow_head_width, gene_length * 0.3)

    half_height = arrow_height / 2

    if strand >= 0:  # Forward strand (left to right arrow)
        # Arrow pointing right
        path = f"M {x_start},{y_center - half_height} " \
               f"L {x_end - head_length},{y_center - half_height} " \
               f"L {x_end - head_length},{y_center - half_height - 0.1} " \
               f"L {x_end},{y_center} " \
               f"L {x_end - head_length},{y_center + half_height + 0.1} " \
               f"L {x_end - head_length},{y_center + half_height} " \
               f"L {x_start},{y_center + half_height} Z"
    else:  # Reverse strand (right to left arrow)
        # Arrow pointing left
        path = f"M {x_end},{y_center - half_height} " \
               f"L {x_start + head_length},{y_center - half_height} " \
               f"L {x_start + head_length},{y_center - half_height - 0.1} " \
               f"L {x_start},{y_center} " \
               f"L {x_start + head_length},{y_center + half_height + 0.1} " \
               f"L {x_start + head_length},{y_center + half_height} " \
               f"L {x_end},{y_center + half_height} Z"

    return {
        "path": path,
        "color": color,
        "label": label,
        "hover_text": hover_text,
        "x_center": (x_start + x_end) / 2,
        "y_center": y_center,
    }


def create_locus_figure(
    orfs: pl.DataFrame,
    contig_id: str,
    hit_to_contig_mapping: Optional[pl.DataFrame] = None,
    interproscan: Optional[pl.DataFrame] = None,
    defensefinder_genes: Optional[pl.DataFrame] = None,
    downstream_orfs: Optional[list] = None,
    title: str = None,
    height: int = 300,
) -> go.Figure:
    """
    Create an interactive locus visualization figure.

    Args:
        orfs: DataFrame with orf_id, contig, start, end, strand columns
        contig_id: ID of contig to visualize
        hit_to_contig_mapping: Mapping from BLAST hits to contig ORF IDs
        interproscan: InterProScan results for domain coloring
        defensefinder_genes: DefenseFinder gene annotations
        downstream_orfs: List of downstream ORF IDs (co-transcribed)
        title: Plot title
        height: Figure height in pixels

    Returns:
        Plotly Figure object
    """
    # Filter ORFs for this contig
    contig_orfs = orfs.filter(pl.col("contig") == contig_id).sort("start")

    if len(contig_orfs) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No ORFs found for this contig",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    # Prepare data structures for coloring
    # Get query hit ORFs from the hit-to-contig mapping
    query_hit_orfs = set()
    if hit_to_contig_mapping is not None and len(hit_to_contig_mapping) > 0:
        # The mapping contains contig_orf_id which matches our orf_id format
        query_hit_orfs = set(hit_to_contig_mapping["contig_orf_id"].to_list())

    downstream_orf_set = set(downstream_orfs) if downstream_orfs else set()

    # Build InterProScan lookup
    ips_lookup = {}
    if interproscan is not None and len(interproscan) > 0:
        for row in interproscan.iter_rows(named=True):
            prot_id = row.get("protein_accession", "")
            if prot_id:
                if prot_id not in ips_lookup:
                    ips_lookup[prot_id] = []
                ips_lookup[prot_id].append({
                    "analysis": row.get("analysis", ""),
                    "signature_desc": row.get("signature_desc", ""),
                    "start": row.get("start", 0),
                    "stop": row.get("stop", 0),
                })

    # Build DefenseFinder lookup
    defense_lookup = {}
    if defensefinder_genes is not None and len(defensefinder_genes) > 0:
        for row in defensefinder_genes.iter_rows(named=True):
            hit_id = row.get("hit_id", "")
            if hit_id:
                defense_lookup[hit_id] = {
                    "gene_name": row.get("gene_name", ""),
                    "type": row.get("type", ""),
                    "subtype": row.get("subtype", ""),
                }

    # Create figure
    fig = go.Figure()

    # Calculate x-axis range
    min_pos = contig_orfs["start"].min()
    max_pos = contig_orfs["end"].max()
    x_range = [min_pos - (max_pos - min_pos) * 0.05, max_pos + (max_pos - min_pos) * 0.05]

    shapes = []
    annotations = []
    hover_traces = []

    y_center = 0.5

    for row in contig_orfs.iter_rows(named=True):
        orf_id = row["orf_id"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        # Determine color based on annotations (priority order)
        color = NO_ANNOTATION_COLOR
        annotation_text = "No annotation"

        # Check if it's a query hit (highest priority - red)
        if orf_id in query_hit_orfs:
            color = QUERY_HIT_COLOR
            annotation_text = "Query Hit (Recombinase)"

        # Check if it's a downstream co-transcribed gene
        elif orf_id in downstream_orf_set:
            color = DOWNSTREAM_COLOR
            annotation_text = "Co-transcribed downstream"

        # Check DefenseFinder annotations
        elif orf_id in defense_lookup:
            defense_info = defense_lookup[orf_id]
            defense_type = defense_info.get("type", "")
            color = DEFENSE_SYSTEM_COLORS.get(defense_type, DEFENSE_SYSTEM_COLORS["default_defense"])
            annotation_text = f"Defense: {defense_info.get('gene_name', '')} ({defense_type})"

        # Check InterProScan annotations
        elif orf_id in ips_lookup:
            domains = ips_lookup[orf_id]
            if domains:
                # Use the first domain's analysis type for coloring
                analysis = domains[0].get("analysis", "")
                color = DOMAIN_COLORS.get(analysis, DOMAIN_COLORS["default_interproscan"])
                desc = domains[0].get("signature_desc", "Unknown")
                annotation_text = f"{analysis}: {desc}"

        # Create hover text
        hover_text = f"<b>{orf_id}</b><br>" \
                    f"Position: {start:,} - {end:,}<br>" \
                    f"Strand: {'+' if strand > 0 else '-'}<br>" \
                    f"Length: {abs(end - start):,} bp<br>" \
                    f"Annotation: {annotation_text}"

        # Create arrow shape
        arrow = create_gene_arrow(
            x_start=start,
            x_end=end,
            y_center=y_center,
            strand=strand,
            color=color,
            label=orf_id.split("_")[-1] if "_" in orf_id else orf_id,
            hover_text=hover_text,
        )

        # Add shape
        shapes.append(
            dict(
                type="path",
                path=arrow["path"],
                fillcolor=color,
                line=dict(color="black", width=1),
                layer="above",
            )
        )

        # Add invisible scatter point for hover
        hover_traces.append(
            go.Scatter(
                x=[arrow["x_center"]],
                y=[y_center],
                mode="markers",
                marker=dict(size=20, opacity=0),
                hoverinfo="text",
                hovertext=hover_text,
                showlegend=False,
            )
        )

    # Add all hover traces
    for trace in hover_traces:
        fig.add_trace(trace)

    # Add shapes
    fig.update_layout(shapes=shapes)

    # Add scale bar
    scale_length = (max_pos - min_pos) / 5
    scale_text = f"{scale_length/1000:.1f} kb" if scale_length >= 1000 else f"{scale_length:.0f} bp"

    fig.add_shape(
        type="line",
        x0=min_pos, x1=min_pos + scale_length,
        y0=0.1, y1=0.1,
        line=dict(color="black", width=2),
    )
    fig.add_annotation(
        x=min_pos + scale_length/2,
        y=0.05,
        text=scale_text,
        showarrow=False,
        font=dict(size=10),
    )

    # Update layout
    fig.update_layout(
        title=title or f"Locus: {contig_id}",
        xaxis=dict(
            title="Position (bp)",
            range=x_range,
            showgrid=True,
            gridcolor="lightgray",
        ),
        yaxis=dict(
            range=[0, 1],
            showticklabels=False,
            showgrid=False,
        ),
        height=height,
        hovermode="closest",
        plot_bgcolor="white",
        showlegend=False,
    )

    return fig


def create_defense_locus_figure(
    orfs: pl.DataFrame,
    contig_id: str,
    defensefinder_genes: pl.DataFrame,
    defensefinder_systems: Optional[pl.DataFrame] = None,
    height: int = 350,
) -> go.Figure:
    """
    Create a specialized locus view for defense systems.
    Shows genes colored by defense system type with system boundaries.

    Args:
        orfs: DataFrame with ORF information
        contig_id: Contig to visualize
        defensefinder_genes: DefenseFinder gene annotations
        defensefinder_systems: DefenseFinder system information (for boundaries)
        height: Figure height

    Returns:
        Plotly Figure
    """
    # Filter ORFs for this contig
    contig_orfs = orfs.filter(pl.col("contig") == contig_id).sort("start")

    if len(contig_orfs) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No ORFs found for this contig",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    # Build defense gene lookup
    defense_lookup = {}
    if defensefinder_genes is not None and len(defensefinder_genes) > 0:
        for row in defensefinder_genes.iter_rows(named=True):
            hit_id = row.get("hit_id", "")
            if hit_id:
                defense_lookup[hit_id] = {
                    "gene_name": row.get("gene_name", ""),
                    "type": row.get("type", ""),
                    "subtype": row.get("subtype", ""),
                }

    fig = go.Figure()

    # Calculate x-axis range
    min_pos = contig_orfs["start"].min()
    max_pos = contig_orfs["end"].max()
    x_range = [min_pos - (max_pos - min_pos) * 0.05, max_pos + (max_pos - min_pos) * 0.05]

    shapes = []
    hover_traces = []
    y_center = 0.5

    # Track defense systems found
    defense_types_found = set()

    for row in contig_orfs.iter_rows(named=True):
        orf_id = row["orf_id"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        # Check DefenseFinder annotations
        if orf_id in defense_lookup:
            defense_info = defense_lookup[orf_id]
            defense_type = defense_info.get("type", "Unknown")
            defense_types_found.add(defense_type)
            color = DEFENSE_SYSTEM_COLORS.get(defense_type, DEFENSE_SYSTEM_COLORS["default_defense"])
            gene_name = defense_info.get("gene_name", "")
            subtype = defense_info.get("subtype", "")
            annotation_text = f"{gene_name} ({defense_type})"
        else:
            color = NO_ANNOTATION_COLOR
            annotation_text = "Not part of defense system"

        hover_text = f"<b>{orf_id}</b><br>" \
                    f"Position: {start:,} - {end:,}<br>" \
                    f"Strand: {'+' if strand > 0 else '-'}<br>" \
                    f"Length: {abs(end - start):,} bp<br>" \
                    f"{annotation_text}"

        # Create arrow shape
        arrow = create_gene_arrow(
            x_start=start,
            x_end=end,
            y_center=y_center,
            strand=strand,
            color=color,
        )

        shapes.append(
            dict(
                type="path",
                path=arrow["path"],
                fillcolor=color,
                line=dict(color="black", width=1),
                layer="above",
            )
        )

        hover_traces.append(
            go.Scatter(
                x=[arrow["x_center"]],
                y=[y_center],
                mode="markers",
                marker=dict(size=20, opacity=0),
                hoverinfo="text",
                hovertext=hover_text,
                showlegend=False,
            )
        )

    # Add hover traces
    for trace in hover_traces:
        fig.add_trace(trace)

    # Add legend traces for defense types found
    for dtype in sorted(defense_types_found):
        color = DEFENSE_SYSTEM_COLORS.get(dtype, DEFENSE_SYSTEM_COLORS["default_defense"])
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=15, color=color, symbol="square"),
                name=dtype,
                showlegend=True,
            )
        )

    # Add shapes
    fig.update_layout(shapes=shapes)

    # Update layout
    fig.update_layout(
        title=f"Defense Systems: {contig_id}",
        xaxis=dict(
            title="Position (bp)",
            range=x_range,
            showgrid=True,
            gridcolor="lightgray",
        ),
        yaxis=dict(
            range=[0, 1],
            showticklabels=False,
            showgrid=False,
        ),
        height=height,
        hovermode="closest",
        plot_bgcolor="white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
        ),
    )

    return fig


def get_contig_choices(orfs: pl.DataFrame) -> list[str]:
    """Get list of unique contig IDs for dropdown selection."""
    if orfs is None or len(orfs) == 0:
        return []
    return sorted(orfs["contig"].unique().to_list())


def get_contigs_with_hits(
    orfs: pl.DataFrame,
    blast_results: pl.DataFrame,
    cotranscription_mapping: Optional[pl.DataFrame] = None,
) -> list[str]:
    """
    Get list of contigs that contain query hits.

    Since contigs are extracted around BLAST hits, every contig by definition
    contains a query hit. This function returns all unique contigs, or uses
    the cotranscription mapping if provided to identify the specific hit ORFs.
    """
    if orfs is None or len(orfs) == 0:
        return []

    # All contigs have hits by definition (they were extracted because of BLAST hits)
    # Return all unique contigs
    return sorted(orfs["contig"].unique().to_list())


def get_contigs_with_defense(
    orfs: pl.DataFrame,
    defensefinder_genes: pl.DataFrame,
) -> list[str]:
    """Get list of contigs that contain defense system genes."""
    if orfs is None or defensefinder_genes is None or len(defensefinder_genes) == 0:
        return []

    defense_orfs = set(defensefinder_genes["hit_id"].to_list())

    contigs_with_defense = set()
    for row in orfs.iter_rows(named=True):
        if row["orf_id"] in defense_orfs:
            contigs_with_defense.add(row["contig"])

    return sorted(list(contigs_with_defense))


def create_color_legend_html() -> str:
    """Create HTML legend for gene colors."""
    legend_items = [
        (QUERY_HIT_COLOR, "Query Hit (Recombinase)"),
        (DOWNSTREAM_COLOR, "Co-transcribed Downstream"),
        ("#e74c3c", "Defense System Gene"),  # Representative defense color
        (NO_ANNOTATION_COLOR, "No Annotation"),
    ]

    html = '<div style="display: flex; flex-wrap: wrap; gap: 15px; margin: 10px 0;">'
    for color, label in legend_items:
        html += f'''
        <div style="display: flex; align-items: center; gap: 5px;">
            <div style="width: 20px; height: 12px; background-color: {color};
                        border: 1px solid black;"></div>
            <span style="font-size: 0.85rem;">{label}</span>
        </div>
        '''
    html += '</div>'

    return html
