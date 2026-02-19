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
    # DefenseFinder system types - navy blue color scheme
    # Matching the antiphage/defense system biology figure aesthetic
    "RM": "#1a3a4a",           # Restriction-Modification (darkest navy)
    "CRISPR": "#234b5e",       # CRISPR-Cas (dark navy)
    "Abi": "#2c5d73",          # Abortive infection (navy)
    "TA": "#1a3a4a",           # Toxin-Antitoxin (darkest navy)
    "BREX": "#356f88",         # BREX (medium navy)
    "DISARM": "#3e819d",       # DISARM (navy-teal)
    "Gabija": "#2c5d73",       # Gabija (navy)
    "Hachiman": "#234b5e",     # Hachiman (dark navy)
    "Lamassu": "#356f88",      # Lamassu (medium navy)
    "Lamassu-Fam": "#356f88",  # Lamassu family (medium navy)
    "Thoeris": "#3e819d",      # Thoeris (navy-teal)
    "Zorya": "#2c5d73",        # Zorya (navy)
    "Druantia": "#234b5e",     # Druantia (dark navy)
    "Kiwa": "#356f88",         # Kiwa (medium navy)
    "Wadjet": "#1a3a4a",       # Wadjet (darkest navy)
    "Septu": "#3e819d",        # Septu (navy-teal)
    "RosmerTA": "#1a3a4a",     # RosmerTA (darkest navy - TA system)
    "MazEF": "#234b5e",        # MazEF (dark navy - TA system)
    "PD-Lambda-1": "#2c5d73",  # PD-Lambda (navy)
    "Dodola": "#356f88",       # Dodola (medium navy)
    "AbiH": "#234b5e",         # AbiH (dark navy)
    "AbiC": "#2c5d73",         # AbiC (navy)
    "AbiJ": "#356f88",         # AbiJ (medium navy)
    "AbiE": "#3e819d",         # AbiE (navy-teal)
    "PrrC": "#1a3a4a",         # PrrC (darkest navy)
    "RloC": "#234b5e",         # RloC (dark navy)
    "default_defense": "#2c5d73",  # Default: navy blue
}

# Special colors for locus viewer gene categories
# Palette inspired by defense system biology figures (dark teal / slate / cool gray)
QUERY_HIT_COLOR = "#2d3d3d"      # Darkest teal-black for query hits (recombinases)
DOWNSTREAM_COLOR = "#5a8a7a"     # Medium green-teal for co-transcribed downstream genes
DEFENSE_COLOR = "#6a9eae"        # Light blue-teal for defense system genes
MGE_COLOR = "#4a5a6a"            # Dark slate-blue for mobile genetic elements
INTERPROSCAN_OTHER_COLOR = "#b0b0b0"  # Neutral gray for other InterProScan domains
NO_ANNOTATION_COLOR = "#ffffff"  # White for unannotated genes

# MGE-associated Pfam accessions for identifying mobile genetic elements
MGE_PFAM_ACCESSIONS = {
    "PF00665",   # rve - Integrase core domain
    "PF01609",   # DDE_Tnp_1
    "PF01527",   # HTH_Tnp_1
    "PF01526",   # DDE_Tnp_IS3 (IS3/IS911)
    "PF01610",   # DDE_Tnp_4 (IS4)
    "PF02371",   # Transposase_mut (MuDR)
    "PF02914",   # DDE_Tnp_IS1
    "PF03400",   # DDE_Tnp_IS21
    "PF05598",   # DUF772 (IS200/IS605)
    "PF00589",   # Phage_integrase
    "PF13009",   # DDE transposase variant
    "PF13276",   # DDE transposase variant
    "PF13358",   # DDE transposase variant
    "PF03389",   # MobA/MobL mobilization
    "PF01076",   # Relaxase/TraI conjugative
}

# MGE keyword patterns for description matching
MGE_KEYWORDS = [
    "transposase", "integrase", "recombinase", "conjugal transfer",
    "mobilization", "insertion sequence", "is element", "conjugative",
]


def is_mge_domain(signature_accession: str, signature_desc: str) -> bool:
    """Check if a domain annotation indicates a mobile genetic element."""
    if signature_accession and signature_accession in MGE_PFAM_ACCESSIONS:
        return True
    if signature_desc:
        desc_lower = signature_desc.lower()
        for keyword in MGE_KEYWORDS:
            if keyword in desc_lower:
                return True
    return False


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
    Create a boxarrow gene shape for Plotly (geneviewer style).

    Boxarrow style: a pentagon shape with a rectangular body and the end
    converging to a point. The arrow head stays within the same height as
    the body (no extension beyond). This matches the geneviewer R package
    boxarrow marker style.

    Shape for forward strand (+):
        _______________
       |               \\
       |                >
       |_______________/

    Shape for reverse strand (-):
         _______________
        /               |
       <                |
        \\_______________|

    Args:
        x_start: Start position of gene
        x_end: End position of gene
        y_center: Y coordinate for center of arrow
        strand: 1 for forward, -1 for reverse
        color: Fill color for the arrow
        label: Gene label
        hover_text: Text to display on hover
        arrow_height: Height of arrow body
        arrow_head_width: Relative width of arrow head (as fraction of gene length)

    Returns:
        Dictionary with shape and annotation data
    """
    gene_length = abs(x_end - x_start)
    head_length = min(gene_length * arrow_head_width, gene_length * 0.3)

    half_height = arrow_height / 2

    if strand >= 0:  # Forward strand (left to right arrow)
        # Pentagon boxarrow pointing right
        # 5 points: top-left -> top-right (before head) -> point -> bottom-right (before head) -> bottom-left
        path = f"M {x_start},{y_center - half_height} " \
               f"L {x_end - head_length},{y_center - half_height} " \
               f"L {x_end},{y_center} " \
               f"L {x_end - head_length},{y_center + half_height} " \
               f"L {x_start},{y_center + half_height} Z"
    else:  # Reverse strand (right to left arrow)
        # Pentagon boxarrow pointing left
        # 5 points: top-right -> top-left (after head) -> point -> bottom-left (after head) -> bottom-right
        path = f"M {x_end},{y_center - half_height} " \
               f"L {x_start + head_length},{y_center - half_height} " \
               f"L {x_start},{y_center} " \
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
                    "signature_accession": row.get("signature_accession", ""),
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
            if orf_id in ips_lookup:
                domains = ips_lookup[orf_id]
                seen = set()
                unique_descs = []
                for d in domains:
                    key = (d.get("analysis", ""), d.get("signature_desc", "Unknown"))
                    if key not in seen:
                        seen.add(key)
                        unique_descs.append(f"{key[0]}: {key[1]}")
                annotation_text = "<br>".join(unique_descs) if unique_descs else "Co-transcribed downstream"
            else:
                annotation_text = "Co-transcribed downstream"

        # Check DefenseFinder annotations
        elif orf_id in defense_lookup:
            defense_info = defense_lookup[orf_id]
            defense_type = defense_info.get("type", "")
            color = DEFENSE_COLOR
            annotation_text = f"Defense: {defense_info.get('gene_name', '')} ({defense_type})"

        # Check for MGE-associated domains (before general InterProScan)
        elif orf_id in ips_lookup and any(
            is_mge_domain(d.get("signature_accession", ""), d.get("signature_desc", ""))
            for d in ips_lookup[orf_id]
        ):
            # Find the matching MGE domain for annotation text
            mge_domain = next(
                d for d in ips_lookup[orf_id]
                if is_mge_domain(d.get("signature_accession", ""), d.get("signature_desc", ""))
            )
            color = MGE_COLOR
            desc = mge_domain.get("signature_desc", "Unknown")
            annotation_text = f"MGE: {desc}"

        # Check InterProScan annotations (other domains - gray)
        elif orf_id in ips_lookup:
            domains = ips_lookup[orf_id]
            if domains:
                color = INTERPROSCAN_OTHER_COLOR
                analysis = domains[0].get("analysis", "")
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
    hit_to_contig_mapping: Optional[pl.DataFrame] = None,
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
        hit_to_contig_mapping: Mapping of query hits to contig ORFs
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

    # Build query hit set
    query_hit_orfs = set()
    if hit_to_contig_mapping is not None and len(hit_to_contig_mapping) > 0:
        query_hit_orfs = set(hit_to_contig_mapping["contig_orf_id"].to_list())

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

        # Check query hits first (highest priority)
        if orf_id in query_hit_orfs:
            color = QUERY_HIT_COLOR
            annotation_text = "Query Hit (Recombinase)"
        # Check DefenseFinder annotations
        elif orf_id in defense_lookup:
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

    # Add query hit legend entry if any query hits are present
    if query_hit_orfs:
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=15, color=QUERY_HIT_COLOR, symbol="square"),
                name="Query Hit (Recombinase)",
                showlegend=True,
            )
        )

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
        (DEFENSE_COLOR, "Defense System Gene"),
        (MGE_COLOR, "Mobile Genetic Element"),
        (INTERPROSCAN_OTHER_COLOR, "Other Domain"),
        (NO_ANNOTATION_COLOR, "Unannotated"),
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


def get_contigs_with_mge(
    orfs: pl.DataFrame,
    interproscan: pl.DataFrame,
) -> list[str]:
    """Get contigs containing non-query MGE-associated genes."""
    if orfs is None or interproscan is None or len(interproscan) == 0:
        return []

    # Find ORF IDs with MGE-associated InterProScan hits
    mge_orfs = set()
    for row in interproscan.iter_rows(named=True):
        prot_id = row.get("protein_accession", "")
        sig_acc = row.get("signature_accession", "")
        sig_desc = row.get("signature_desc", "")
        if prot_id and is_mge_domain(sig_acc, sig_desc):
            mge_orfs.add(prot_id)

    # Map MGE ORFs to their contigs
    contigs_with_mge = set()
    for row in orfs.iter_rows(named=True):
        if row["orf_id"] in mge_orfs:
            contigs_with_mge.add(row["contig"])

    return sorted(list(contigs_with_mge))
