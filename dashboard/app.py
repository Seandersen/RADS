"""
RADS Results Explorer - PyShiny Dashboard
Interactive visualization of RADS pipeline results.
"""

from pathlib import Path
from shiny import App, reactive, render, ui
import plotly.express as px
import plotly.graph_objects as go
import polars as pl

# Project root is one level above the dashboard directory
PROJECT_ROOT = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"

# Import data loader
import sys
sys.path.insert(0, str(Path(__file__).parent))
from utils.data_loader import (
    find_results_dir,
    load_blast_results,
    load_genome_manifest,
    load_cotranscription_results,
    load_contig_orfs,
    load_interproscan_results,
    get_summary_stats,
    load_defensefinder_systems,
    load_defensefinder_genes,
    load_pipeline_metrics,
    load_hit_to_contig_mapping,
    load_defense_scores,
    load_binomial_results,
)
from utils.locus_viewer import (
    create_locus_figure,
    create_defense_locus_figure,
    get_contig_choices,
    get_contigs_with_hits,
    get_contigs_with_defense,
    get_contigs_with_mge,
    create_color_legend_html,
    is_mge_domain,
)


# Teal/slate color palette for charts (matching defense system biology theme)
TEAL_PALETTE = [
    "#2d4a4a",  # Darkest teal
    "#3d5a5a",  # Dark teal
    "#4a6670",  # Slate blue-teal
    "#5a7a7a",  # Medium teal
    "#6b9090",  # Medium-light teal
    "#7a9e9e",  # Light teal
    "#8fb3b3",  # Lighter teal/sage
    "#a8c4c4",  # Light sage
    "#b5cece",  # Pale sage
    "#c5dada",  # Very pale sage
]

# Set plotly default color sequence
px.defaults.color_discrete_sequence = TEAL_PALETTE

# Chart styling constants
CHART_COLORS = {
    "primary": "#2d4a4a",
    "secondary": "#4a6670",
    "tertiary": "#6b9090",
    "light": "#8fb3b3",
    "background": "#f5f8f8",
    "grid": "#e0e8e8",
}

# Custom CSS for teal/slate theme
CUSTOM_CSS = """
<style>
/* Sidebar styling */
.bslib-sidebar-layout > .sidebar {
    background-color: #2d4a4a !important;
    color: #ffffff !important;
}

.bslib-sidebar-layout > .sidebar .form-label,
.bslib-sidebar-layout > .sidebar label,
.bslib-sidebar-layout > .sidebar h4,
.bslib-sidebar-layout > .sidebar h5 {
    color: #ffffff !important;
}

.bslib-sidebar-layout > .sidebar hr {
    border-color: #4a6670 !important;
}

/* Value boxes */
.bslib-value-box {
    background-color: #3d5a5a !important;
    border: none !important;
}

.bslib-value-box .value-box-title {
    color: #b5cece !important;
}

.bslib-value-box .value-box-value {
    color: #ffffff !important;
}

.bslib-value-box .value-box-showcase {
    color: #8fb3b3 !important;
}

/* Nav tabs */
.nav-link.active {
    background-color: #4a6670 !important;
    color: #ffffff !important;
    border-color: #4a6670 !important;
}

.nav-link {
    color: #2d4a4a !important;
}

.nav-link:hover {
    color: #3d5a5a !important;
    border-color: #6b9090 !important;
}

/* Cards */
.card {
    border-color: #b5cece !important;
}

.card-header {
    background-color: #f5f8f8 !important;
    border-color: #b5cece !important;
    color: #2d4a4a !important;
}

/* Buttons */
.btn-primary, .btn-default {
    background-color: #4a6670 !important;
    border-color: #3d5a5a !important;
    color: #ffffff !important;
}

.btn-primary:hover, .btn-default:hover {
    background-color: #3d5a5a !important;
    border-color: #2d4a4a !important;
}

/* Sliders */
.irs--shiny .irs-bar {
    background: #4a6670 !important;
    border-color: #4a6670 !important;
}

.irs--shiny .irs-single {
    background-color: #3d5a5a !important;
}

.irs--shiny .irs-handle {
    border-color: #2d4a4a !important;
}

/* Select inputs */
.selectize-input {
    border-color: #b5cece !important;
}

.selectize-input.focus {
    border-color: #6b9090 !important;
    box-shadow: 0 0 0 0.2rem rgba(74, 102, 112, 0.25) !important;
}

/* Radio buttons */
.shiny-input-container .radio label input:checked + span::before {
    background-color: #4a6670 !important;
    border-color: #4a6670 !important;
}

/* Page background */
body {
    background-color: #f5f8f8 !important;
}

/* Main content area */
.bslib-sidebar-layout > .main {
    background-color: #ffffff !important;
}
</style>
"""


# UI Definition
app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h4("RADS Explorer"),
        ui.hr(),
        ui.input_select(
            "sample_select",
            "Select Sample",
            choices=[],
        ),
        ui.output_ui("sample_name_display"),
        ui.hr(),
        ui.h5("Filters"),
        ui.input_slider(
            "identity_filter",
            "Min % Identity",
            min=0,
            max=100,
            value=30,
        ),
        ui.input_slider(
            "length_filter",
            "Min Alignment Length",
            min=0,
            max=500,
            value=0,
        ),
        ui.hr(),
        ui.download_button("download_blast", "Download BLAST Results"),
        width=280,
    ),
    ui.navset_card_tab(
        ui.nav_panel(
            "Summary",
            ui.card(
                ui.card_header("Pipeline Results Breakdown"),
                ui.output_ui("pipeline_breakdown_plot"),
                full_screen=True,
            ),
            ui.layout_columns(
                ui.card(
                    ui.card_header("Pipeline Overview"),
                    ui.output_ui("pipeline_status"),
                ),
                ui.layout_columns(
                    ui.div(
                        ui.value_box(
                            "Total Genomes",
                            ui.output_text("stat_genomes"),
                            showcase=ui.HTML('<i class="fa-solid fa-dna" style="font-size: 2rem;"></i>'),
                            theme="primary",
                        ),
                        title="Count of genomes in the genome manifest",
                    ),
                    ui.div(
                        ui.value_box(
                            "BLAST Hits",
                            ui.output_text("stat_blast_hits"),
                            showcase=ui.HTML('<i class="fa-solid fa-bullseye" style="font-size: 2rem;"></i>'),
                            theme="success",
                        ),
                        title="Total rows in master BLAST results",
                    ),
                    ui.div(
                        ui.value_box(
                            "Genomes with Hits",
                            ui.output_text("stat_genomes_hits"),
                            showcase=ui.HTML('<i class="fa-solid fa-check" style="font-size: 2rem;"></i>'),
                            theme="info",
                        ),
                        title="Number of unique genomes with at least one BLAST hit",
                    ),
                    ui.div(
                        ui.value_box(
                            "Hits per Mb",
                            ui.output_text("stat_hits_per_mb"),
                            showcase=ui.HTML('<i class="fa-solid fa-chart-line" style="font-size: 2rem;"></i>'),
                            theme="secondary",
                        ),
                        title="BLAST hits divided by total input megabases",
                    ),
                    ui.div(
                        ui.value_box(
                            "Defense Systems",
                            ui.output_text("stat_defense_systems"),
                            showcase=ui.HTML('<i class="fa-solid fa-shield" style="font-size: 2rem;"></i>'),
                            theme="danger",
                        ),
                        title="Count of unique defense system IDs from DefenseFinder",
                    ),
                    ui.div(
                        ui.value_box(
                            "Co-transcribed Pairs",
                            ui.output_text("stat_cotx"),
                            showcase=ui.HTML('<i class="fa-solid fa-link" style="font-size: 2rem;"></i>'),
                            theme="warning",
                        ),
                        title="Number of co-transcribed gene pairs identified downstream of query hits",
                    ),
                    ui.div(
                        ui.value_box(
                            "Known Defense Systems per Contig",
                            ui.output_text("stat_discovery_per_contig"),
                            showcase=ui.HTML('<i class="fa-solid fa-magnifying-glass" style="font-size: 2rem;"></i>'),
                            theme="light",
                        ),
                        title="Ratio of unique defense systems to total contigs",
                    ),
                    ui.div(
                        ui.value_box(
                            "Genomes with Known Defense (%)",
                            ui.output_text("stat_discovery_per_genome"),
                            showcase=ui.HTML('<i class="fa-solid fa-bacteria" style="font-size: 2rem;"></i>'),
                            theme="dark",
                        ),
                        title="Proportion of genomes with at least one defense system",
                    ),
                    col_widths=[3, 3, 3, 3],
                ),
                col_widths=[4, 8],
            ),
        ),
        ui.nav_panel(
            "BLAST Results",
            ui.layout_columns(
                ui.card(
                    ui.card_header("Identity vs Length"),
                    ui.output_ui("blast_scatter_plot"),
                    full_screen=True,
                ),
                ui.card(
                    ui.card_header("Identity Distribution"),
                    ui.output_ui("identity_histogram"),
                    full_screen=True,
                ),
                col_widths=[6, 6],
            ),
            ui.card(
                ui.card_header("BLAST Hits Table"),
                ui.output_data_frame("blast_table"),
                full_screen=True,
            ),
        ),
        ui.nav_panel(
            "Contig Analysis",
            ui.layout_columns(
                ui.card(
                    ui.card_header("ORF Length Distribution"),
                    ui.output_ui("orf_length_plot"),
                    full_screen=True,
                ),
                ui.card(
                    ui.card_header("ORFs per Contig"),
                    ui.output_ui("orfs_per_contig_plot"),
                    full_screen=True,
                ),
                col_widths=[6, 6],
            ),
            ui.card(
                ui.card_header("ORF Details"),
                ui.output_data_frame("orf_table"),
                full_screen=True,
            ),
        ),
        ui.nav_panel(
            "Co-transcription",
            ui.card(
                ui.card_header("Co-transcribed Gene Pairs"),
                ui.output_data_frame("cotx_table"),
                full_screen=True,
            ),
            ui.card(
                ui.card_header("Distance Distribution"),
                ui.output_ui("cotx_distance_plot"),
                full_screen=True,
            ),
            ui.card(
                ui.card_header("Domain Annotations for Co-transcribed Genes"),
                ui.layout_columns(
                    ui.output_ui("cotx_domain_bar_chart"),
                    ui.output_data_frame("cotx_domain_table"),
                    col_widths=[6, 6],
                ),
                full_screen=True,
            ),
            ui.card(
                ui.card_header("DefenseFinder Hits in Co-transcribed Genes"),
                ui.output_ui("cotx_defense_summary_text"),
                ui.output_data_frame("cotx_defense_table"),
                full_screen=True,
            ),
            ui.card(
                ui.card_header("Defense Score Distribution"),
                ui.output_ui("cotx_defense_score_plot"),
                full_screen=True,
            ),
            ui.card(
                ui.card_header("Binomial Domain Enrichment"),
                ui.layout_columns(
                    ui.output_ui("binomial_bar_chart"),
                    ui.output_data_frame("binomial_table"),
                    col_widths=[6, 6],
                ),
                full_screen=True,
            ),
        ),
        ui.nav_panel(
            "Locus Viewer",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.h5("Contig Selection"),
                    ui.input_radio_buttons(
                        "locus_filter_type",
                        "Filter contigs by:",
                        choices={
                            "all": "All contigs",
                            "hits": "Contigs with query hits",
                        },
                        selected="hits",
                    ),
                    ui.hr(),
                    ui.layout_columns(
                        ui.input_action_button(
                            "locus_select_all",
                            "Select All Filtered",
                            class_="btn-sm",
                        ),
                        ui.input_action_button(
                            "locus_clear_all",
                            "Clear All",
                            class_="btn-sm",
                        ),
                        col_widths=[6, 6],
                    ),
                    ui.output_ui("locus_selected_count"),
                    ui.input_selectize(
                        "locus_contig_select",
                        "Select Contig(s):",
                        choices=[],
                        multiple=True,
                    ),
                    ui.hr(),
                    ui.h5("Locus Filters"),
                    ui.panel_conditional(
                        "input.locus_has_cotranscribed",
                        ui.input_slider(
                            "locus_defense_score_range",
                            "Defense Score Range",
                            min=0.0,
                            max=1.0,
                            value=[0.0, 1.0],
                            step=0.05,
                        ),
                        ui.input_slider(
                            "locus_binomial_pvalue",
                            "Binomial p-value threshold",
                            min=0.001,
                            max=0.1,
                            value=0.05,
                            step=0.001,
                        ),
                    ),
                    ui.input_checkbox(
                        "locus_has_cotranscribed",
                        "Only contigs with co-transcribed genes",
                        value=False,
                    ),
                    ui.input_checkbox(
                        "locus_has_defense",
                        "Only contigs with defense systems",
                        value=False,
                    ),
                    ui.input_checkbox(
                        "locus_has_mge",
                        "Only contigs with other MGEs",
                        value=False,
                    ),
                    ui.input_selectize(
                        "locus_defense_type_filter",
                        "Defense System Type:",
                        choices=[],
                        multiple=True,
                    ),
                    ui.input_action_button(
                        "locus_reset_filters",
                        "Reset Filters",
                    ),
                    ui.hr(),
                    ui.h6("Color Legend"),
                    ui.output_ui("locus_color_legend"),
                    width=280,
                ),
                ui.card(
                    ui.card_header(ui.output_text("locus_title")),
                    ui.output_ui("locus_viewer_plot"),
                    full_screen=True,
                ),
            ),
            ui.card(
                ui.card_header("ORF Details for Selected Contig"),
                ui.output_data_frame("locus_orf_table"),
                full_screen=True,
            ),
        ),
        ui.nav_panel(
            "Domain Annotations",
            ui.layout_columns(
                ui.card(
                    ui.card_header("Top Domains"),
                    ui.output_ui("domain_frequency_plot"),
                    full_screen=True,
                ),
                ui.card(
                    ui.card_header("Analysis Types"),
                    ui.output_ui("analysis_type_plot"),
                    full_screen=True,
                ),
                col_widths=[6, 6],
            ),
            ui.card(
                ui.card_header("Binomial Domain Enrichment (All Domains)"),
                ui.output_ui("domain_binomial_plot"),
                full_screen=True,
            ),
            ui.card(
                ui.card_header("InterProScan Results"),
                ui.output_data_frame("interproscan_table"),
                full_screen=True,
            ),
        ),
        ui.nav_panel(
            "DefenseFinder",
            # Row 1 — category overview
            ui.layout_columns(
                ui.card(
                    ui.card_header("Defense Systems by Antiphage Category"),
                    ui.output_ui("defense_category_plot"),
                    full_screen=True,
                ),
                ui.card(
                    ui.card_header("Category × System Type Breakdown"),
                    ui.output_ui("defense_category_sunburst"),
                    full_screen=True,
                ),
                col_widths=[6, 6],
            ),
            # Row 2 — per-category drilldown
            ui.card(
                ui.card_header("Individual Systems by Category"),
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_select(
                            "defense_category_select",
                            "Select Category:",
                            choices=[
                                "All",
                                "Abortive infection",
                                "CBASS",
                                "Nucleic acid restriction",
                                "CRISPR-Cas",
                                "Toxin-antitoxin",
                                "Retrons",
                                "tRNA degradation",
                                "Unknown mechanism",
                            ],
                            selected="All",
                        ),
                        width=240,
                    ),
                    ui.output_ui("defense_category_drilldown"),
                ),
                full_screen=True,
            ),
            # Row 3 — original type / subtype charts
            ui.layout_columns(
                ui.card(
                    ui.card_header("Defense System Types"),
                    ui.output_ui("defense_system_type_plot"),
                    full_screen=True,
                ),
                ui.card(
                    ui.card_header("Defense Systems by Subtype (Top 15)"),
                    ui.output_ui("defense_subtype_plot"),
                    full_screen=True,
                ),
                col_widths=[6, 6],
            ),
            ui.card(
                ui.card_header("Defense System Locus Viewer"),
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.input_selectize(
                            "defense_type_filter",
                            "Filter by Defense Type:",
                            choices=[],
                            multiple=True,
                        ),
                        ui.input_select(
                            "defense_contig_select",
                            "Select contig with defense system:",
                            choices=[],
                        ),
                        ui.output_ui("defense_system_info"),
                        width=280,
                    ),
                    ui.output_ui("defense_locus_plot"),
                ),
                full_screen=True,
            ),
            ui.card(
                ui.card_header("Defense Genes Table"),
                ui.output_data_frame("defense_genes_table"),
                full_screen=True,
            ),
        ),
    ),
    title="RADS Results Explorer",
    fillable=True,
)

# Wrap with head content for CSS
app_ui = ui.page_fluid(
    ui.tags.head(
        ui.HTML(CUSTOM_CSS),
        ui.HTML('<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'),
    ),
    app_ui,
)


def server(input, output, session):
    # Initialize sample selector
    @reactive.effect
    def _():
        samples = find_results_dir(str(RESULTS_DIR))
        if samples:
            ui.update_select("sample_select", choices=samples, selected=samples[0])

    # Reactive data loaders
    @reactive.calc
    def results_dir():
        sample = input.sample_select()
        if sample:
            return str(RESULTS_DIR / sample)
        return None

    @reactive.calc
    def blast_data():
        rdir = results_dir()
        if rdir:
            return load_blast_results(rdir)
        return None

    @reactive.calc
    def filtered_blast_data():
        df = blast_data()
        if df is None or len(df) == 0:
            return None

        min_identity = input.identity_filter()
        min_length = input.length_filter()

        filtered = df.filter(
            (df["pident"] >= min_identity) &
            (df["length"] >= min_length)
        )
        return filtered

    @reactive.calc
    def cotx_data():
        rdir = results_dir()
        if rdir:
            return load_cotranscription_results(rdir)
        return None

    @reactive.calc
    def orf_data():
        rdir = results_dir()
        if rdir:
            return load_contig_orfs(rdir)
        return None

    @reactive.calc
    def summary_stats():
        rdir = results_dir()
        if rdir:
            return get_summary_stats(rdir)
        return {}

    @reactive.calc
    def interproscan_data():
        rdir = results_dir()
        if rdir:
            return load_interproscan_results(rdir)
        return None

    @reactive.calc
    def defensefinder_systems():
        rdir = results_dir()
        if rdir:
            return load_defensefinder_systems(rdir)
        return None

    @reactive.calc
    def defensefinder_genes():
        rdir = results_dir()
        if rdir:
            return load_defensefinder_genes(rdir)
        return None

    @reactive.calc
    def defense_scores_data():
        rdir = results_dir()
        if rdir:
            return load_defense_scores(rdir)
        return None

    @reactive.calc
    def hit_to_contig_mapping():
        rdir = results_dir()
        if rdir:
            return load_hit_to_contig_mapping(rdir)
        return None

    @reactive.calc
    def pipeline_metrics():
        rdir = results_dir()
        if rdir:
            return load_pipeline_metrics(rdir)
        return None

    # Sample name display in sidebar
    @render.ui
    def sample_name_display():
        metrics = pipeline_metrics()
        sample = input.sample_select()
        if metrics and metrics.get("sample_name"):
            return ui.tags.div(
                ui.tags.strong("Sample: "),
                ui.tags.span(metrics["sample_name"], style="color: #0d6efd;"),
                style="margin-top: 5px; font-size: 0.9rem;"
            )
        elif sample:
            return ui.tags.div(
                ui.tags.strong("Sample: "),
                ui.tags.span(sample, style="color: #0d6efd;"),
                style="margin-top: 5px; font-size: 0.9rem;"
            )
        return None

    # Summary tab outputs
    @render.text
    def stat_genomes():
        stats = summary_stats()
        return str(stats.get("total_genomes", 0))

    @render.text
    def stat_blast_hits():
        stats = summary_stats()
        return str(stats.get("total_blast_hits", 0))

    @render.text
    def stat_genomes_hits():
        stats = summary_stats()
        return str(stats.get("genomes_with_hits", 0))

    @render.text
    def stat_cotx():
        stats = summary_stats()
        return str(stats.get("cotranscribed_pairs", 0))

    @render.text
    def stat_defense_systems():
        # First try to count from actual DefenseFinder data
        df = defensefinder_systems()
        if df is not None and len(df) > 0:
            # Count unique system IDs
            if "sys_id" in df.columns:
                return str(df["sys_id"].n_unique())
            return str(len(df))
        # Fallback to metrics file
        metrics = pipeline_metrics()
        if metrics and metrics.get("defense_systems", 0) > 0:
            return str(metrics.get("defense_systems", 0))
        return "0"

    @render.text
    def stat_hits_per_mb():
        metrics = pipeline_metrics()
        if metrics and metrics.get('hits_per_mb', 0) > 0:
            return f"{metrics.get('hits_per_mb', 0):.3f}"
        # Calculate from data if not in metrics
        stats = summary_stats()
        if stats.get("total_genomes", 0) > 0:
            # Rough estimate based on average genome size
            blast_hits = stats.get("total_blast_hits", 0)
            # Assume ~5 Mb average genome size
            total_mb = stats.get("total_genomes", 0) * 5
            if total_mb > 0:
                return f"{blast_hits / total_mb:.3f}"
        return "N/A"

    @render.text
    def stat_discovery_per_contig():
        # Calculate from actual data: defense systems per contig analyzed
        df = defensefinder_systems()
        stats = summary_stats()
        contigs = stats.get("total_contigs", 0)

        if df is not None and len(df) > 0 and contigs > 0:
            # Count unique systems
            if "sys_id" in df.columns:
                unique_systems = df["sys_id"].n_unique()
            else:
                unique_systems = len(df)
            rate = unique_systems / contigs
            return f"{rate:.3f}"

        # Fallback to metrics file
        metrics = pipeline_metrics()
        if metrics and metrics.get('discovery_rate_per_contig', 0) > 0:
            return f"{metrics.get('discovery_rate_per_contig', 0):.3f}"
        return "0.000"

    @render.text
    def stat_discovery_per_genome():
        # DefenseFinder runs with --db-type unordered on all_contigs.faa, so the
        # 'replicon' column is the input filename for every row — not genome IDs.
        # Use the metrics JSON which has the correct count from the pipeline.
        metrics = pipeline_metrics()
        if metrics and metrics.get('discovery_rate_per_genome', 0) > 0:
            return f"{metrics.get('discovery_rate_per_genome', 0) * 100:.1f}%"

        # Fallback: estimate from number of unique defense systems vs total genomes
        df = defensefinder_systems()
        stats = summary_stats()
        total_genomes = stats.get("total_genomes", 0)
        if df is not None and len(df) > 0 and total_genomes > 0:
            if "sys_id" in df.columns:
                unique_systems = df["sys_id"].n_unique()
                rate = min(unique_systems / total_genomes, 1.0)
                return f"~{rate * 100:.1f}%"

        return "0.0%"

    @render.ui
    def pipeline_status():
        stats = summary_stats()
        metrics = pipeline_metrics()

        items = [
            f"Genomes processed: {stats.get('total_genomes', 0)}",
            f"Total ORFs found: {stats.get('total_orfs', 0)}",
            f"Extracted contigs: {stats.get('total_contigs', 0)}",
            f"Domain annotations: {stats.get('domain_annotations', 0)}",
        ]

        # Add query info from metrics
        if metrics:
            query_file = metrics.get("query_file", "N/A")
            query_name = metrics.get("query_name", "N/A")
            items.append(f"Query file: {query_file}")
            items.append(f"Query sequence: {query_name}")
            items.append(f"Total input size: {metrics.get('total_input_mb', 0):.2f} Mb")

        return ui.tags.ul([ui.tags.li(item) for item in items])

    @render.ui
    def pipeline_breakdown_plot():
        stats = summary_stats()
        if not stats:
            return ui.p("No pipeline results available.")

        # Get defense system count same way the tile does
        df_def = defensefinder_systems()
        n_defense = 0
        if df_def is not None and len(df_def) > 0 and "sys_id" in df_def.columns:
            n_defense = df_def["sys_id"].n_unique()

        n_genomes = stats.get("total_genomes", 0)
        n_hits = stats.get("total_blast_hits", 0)
        n_genomes_hits = stats.get("genomes_with_hits", 0)
        n_genomes_no_hits = n_genomes - n_genomes_hits
        n_contigs = stats.get("total_contigs", 0)
        n_orfs = stats.get("total_orfs", 0)
        n_cotx = stats.get("cotranscribed_pairs", 0)
        n_domains = stats.get("domain_annotations", 0)

        # Nodes: Genomes -> Genomes with Hits / No Hits -> Contigs -> ORFs -> downstream outputs
        #   0: Genomes
        #   1: Genomes with Hits
        #   2: Genomes without Hits
        #   3: Extracted Contigs
        #   4: ORFs
        #   5: BLAST Hits
        #   6: Co-transcribed Pairs
        #   7: Defense Systems
        #   8: Domain Annotations
        node_labels = [
            f"Genomes ({n_genomes})",
            f"Genomes with Hits ({n_genomes_hits})",
            f"No Hits ({n_genomes_no_hits})",
            f"Extracted Contigs ({n_contigs})",
            f"ORFs ({n_orfs})",
            f"BLAST Hits ({n_hits})",
            f"Co-transcribed Pairs ({n_cotx})",
            f"Defense Systems ({n_defense})",
            f"Domain Annotations ({n_domains})",
        ]
        node_colors = [
            TEAL_PALETTE[0], TEAL_PALETTE[1], TEAL_PALETTE[7],
            TEAL_PALETTE[2], TEAL_PALETTE[3],
            TEAL_PALETTE[4], TEAL_PALETTE[5], TEAL_PALETTE[6], TEAL_PALETTE[8],
        ]

        # Links: source -> target with value
        sources = [0, 0, 1, 3, 4, 4, 4]
        targets = [1, 2, 3, 4, 5, 6, 7]
        values =  [n_genomes_hits, max(n_genomes_no_hits, 1), n_contigs or n_genomes_hits,
                   n_orfs or n_contigs, n_hits, n_cotx, n_defense]
        # Add domain annotations link from ORFs if available
        if n_domains > 0:
            sources.append(4)
            targets.append(8)
            values.append(n_domains)

        def hex_to_rgba(hex_color, alpha=0.4):
            h = hex_color.lstrip("#")
            r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
            return f"rgba({r},{g},{b},{alpha})"

        link_colors = [hex_to_rgba(node_colors[t]) for t in targets]

        import plotly.graph_objects as go
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=20,
                thickness=25,
                label=node_labels,
                color=node_colors,
            ),
            link=dict(
                source=sources,
                target=targets,
                value=values,
                color=link_colors,
            ),
        )])
        fig.update_layout(
            height=350,
            margin=dict(t=20, b=20, l=20, r=20),
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    # BLAST tab outputs
    @render.ui
    def blast_scatter_plot():
        df = filtered_blast_data()
        if df is None or len(df) == 0:
            return ui.p("No data available")

        pdf = df.to_pandas()
        fig = px.scatter(
            pdf,
            x="length",
            y="pident",
            color="evalue",
            color_continuous_scale="teal",
            hover_data=["query_id", "subject_id", "genome", "evalue"],
            labels={"length": "Alignment Length", "pident": "% Identity", "evalue": "E-value"},
        )
        fig.update_layout(
            height=400,
            plot_bgcolor=CHART_COLORS["background"],
            showlegend=False,
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.ui
    def identity_histogram():
        df = filtered_blast_data()
        if df is None or len(df) == 0:
            return ui.p("No data available")

        fig = px.histogram(
            df.to_pandas(),
            x="pident",
            nbins=20,
            labels={"pident": "% Identity"},
            color_discrete_sequence=[CHART_COLORS["primary"]],
        )
        fig.update_layout(height=400, plot_bgcolor=CHART_COLORS["background"])
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.data_frame
    def blast_table():
        df = filtered_blast_data()
        if df is None:
            return None
        return render.DataTable(df.to_pandas(), filters=True, height="400px")

    # Contig analysis outputs
    @render.ui
    def orf_length_plot():
        df = orf_data()
        if df is None or len(df) == 0:
            return ui.p("No data available")

        fig = px.histogram(
            df.to_pandas(),
            x="length",
            nbins=30,
            labels={"length": "ORF Length (bp)"},
            color_discrete_sequence=[CHART_COLORS["secondary"]],
        )
        fig.update_layout(height=350, plot_bgcolor=CHART_COLORS["background"])
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.ui
    def orfs_per_contig_plot():
        df = orf_data()
        if df is None or len(df) == 0:
            return ui.p("No data available")

        counts = df.group_by("contig").len()
        fig = px.bar(
            counts.to_pandas(),
            x="contig",
            y="len",
            labels={"contig": "Contig", "len": "ORF Count"},
            color_discrete_sequence=[CHART_COLORS["tertiary"]],
        )
        median_val = counts["len"].median()
        fig.add_hline(y=median_val, line_dash="dash", line_color="#3d5a5a",
                      annotation_text=f"Median: {median_val:.0f}",
                      annotation_position="top right")
        fig.update_layout(
            xaxis_tickangle=-45,
            margin=dict(b=100),
            height=350,
            plot_bgcolor=CHART_COLORS["background"],
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.data_frame
    def orf_table():
        df = orf_data()
        if df is None:
            return None
        return render.DataTable(df.to_pandas(), filters=True, height="400px")

    # Co-transcription outputs
    @render.data_frame
    def cotx_table():
        df = cotx_data()
        if df is None:
            return None
        # Join with defense scores to add defense_score column
        scores = defense_scores_data()
        if scores is not None and len(scores) > 0 and "downstream_orf" in scores.columns and "downstream_orf" in df.columns:
            score_cols = scores.select(["downstream_orf", "defense_score"])
            df = df.join(score_cols, on="downstream_orf", how="left")
        return render.DataTable(df.to_pandas(), filters=True, height="400px")

    @render.ui
    def cotx_distance_plot():
        df = cotx_data()
        if df is None or len(df) == 0:
            return ui.p("No co-transcribed pairs found")

        fig = px.histogram(
            df.to_pandas(),
            x="distance",
            nbins=20,
            labels={"distance": "Distance (bp)"},
            title="Distance between query and downstream ORF",
            color_discrete_sequence=[CHART_COLORS["tertiary"]],
        )
        fig.update_layout(height=300, plot_bgcolor=CHART_COLORS["background"])
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    # Co-transcription: Domain annotations for co-transcribed genes
    @render.ui
    def cotx_domain_bar_chart():
        cotx = cotx_data()
        ips = interproscan_data()
        if cotx is None or ips is None or len(cotx) == 0 or len(ips) == 0:
            return ui.p("No domain data for co-transcribed genes.")
        if "downstream_orf" not in cotx.columns or "protein_accession" not in ips.columns:
            return ui.p("Missing required columns.")
        downstream_ids = set(cotx["downstream_orf"].to_list())
        cotx_domains = ips.filter(pl.col("protein_accession").is_in(list(downstream_ids)))
        if len(cotx_domains) == 0:
            return ui.p("No domain annotations found for co-transcribed genes.")
        if "signature_desc" not in cotx_domains.columns:
            return ui.p("No signature descriptions available.")
        counts = (
            cotx_domains.filter(pl.col("signature_desc") != "-")
            .group_by("signature_desc").len()
            .sort("len", descending=True).head(15)
        )
        if len(counts) == 0:
            return ui.p("No domains found.")
        fig = px.bar(
            counts.to_pandas(), x="len", y="signature_desc", orientation="h",
            labels={"len": "Count", "signature_desc": "Domain"},
            color_discrete_sequence=[CHART_COLORS["secondary"]],
        )
        fig.update_layout(height=400, yaxis={"categoryorder": "total ascending"},
                          margin=dict(l=200), plot_bgcolor=CHART_COLORS["background"])
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.data_frame
    def cotx_domain_table():
        cotx = cotx_data()
        ips = interproscan_data()
        if cotx is None or ips is None or len(cotx) == 0 or len(ips) == 0:
            return None
        if "downstream_orf" not in cotx.columns or "protein_accession" not in ips.columns:
            return None
        downstream_ids = set(cotx["downstream_orf"].to_list())
        cotx_domains = ips.filter(pl.col("protein_accession").is_in(list(downstream_ids)))
        if len(cotx_domains) == 0:
            return None
        display_cols = [c for c in ["protein_accession", "analysis", "signature_accession",
                                     "signature_desc", "interpro_accession", "interpro_desc"]
                        if c in cotx_domains.columns]
        return render.DataTable(cotx_domains.select(display_cols).to_pandas(), filters=True, height="400px")

    # Co-transcription: DefenseFinder hits
    @render.ui
    def cotx_defense_summary_text():
        cotx = cotx_data()
        genes = defensefinder_genes()
        if cotx is None or genes is None or len(cotx) == 0 or len(genes) == 0:
            return ui.p("No DefenseFinder data for co-transcribed genes.")
        if "downstream_orf" not in cotx.columns or "hit_id" not in genes.columns:
            return ui.p("Missing required columns.")
        downstream_ids = set(cotx["downstream_orf"].to_list())
        gene_ids = set(genes["hit_id"].to_list())
        matches = downstream_ids & gene_ids
        return ui.p(f"{len(matches)} of {len(downstream_ids)} co-transcribed genes have DefenseFinder hits")

    @render.data_frame
    def cotx_defense_table():
        cotx = cotx_data()
        genes = defensefinder_genes()
        if cotx is None or genes is None or len(cotx) == 0 or len(genes) == 0:
            return None
        if "downstream_orf" not in cotx.columns or "hit_id" not in genes.columns:
            return None
        downstream_ids = list(set(cotx["downstream_orf"].to_list()))
        matching = genes.filter(pl.col("hit_id").is_in(downstream_ids))
        if len(matching) == 0:
            return None
        return render.DataTable(matching.to_pandas(), filters=True, height="300px")

    # Co-transcription: Defense score distribution
    @render.ui
    def cotx_defense_score_plot():
        cotx = cotx_data()
        scores = defense_scores_data()
        if cotx is None or scores is None or len(cotx) == 0 or len(scores) == 0:
            return ui.p("No defense scores for co-transcribed genes.")
        if "downstream_orf" not in cotx.columns or "downstream_orf" not in scores.columns:
            return ui.p("Missing required columns.")
        downstream_ids = list(set(cotx["downstream_orf"].to_list()))
        cotx_scores = scores.filter(pl.col("downstream_orf").is_in(downstream_ids))
        scored = cotx_scores.filter(pl.col("defense_score").is_not_null())
        if len(scored) == 0:
            return ui.p("No scored co-transcribed genes.")
        fig = px.histogram(
            scored.to_pandas(), x="defense_score", nbins=20,
            labels={"defense_score": "Defense Score"},
            color_discrete_sequence=[CHART_COLORS["primary"]],
        )
        fig.update_layout(height=350, plot_bgcolor=CHART_COLORS["background"],
                          xaxis=dict(range=[0, 1]))
        fig.add_vline(x=0.1, line_dash="dash", line_color="#6b9090",
                      annotation_text="Low threshold", annotation_position="top right")
        fig.add_vline(x=0.5, line_dash="dash", line_color="#3d5a5a",
                      annotation_text="High threshold", annotation_position="top right")
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    # Co-transcription: Binomial domain enrichment
    @reactive.calc
    def binomial_data():
        rdir = results_dir()
        if rdir:
            return load_binomial_results(rdir)
        return None

    @reactive.calc
    def cotx_binomial_data():
        """Filter binomial results to only domains found in co-transcribed genes."""
        df = binomial_data()
        if df is None or len(df) == 0:
            return None

        cotx = cotx_data()
        ips = interproscan_data()
        if cotx is None or ips is None or len(cotx) == 0 or len(ips) == 0:
            return None

        if "downstream_orf" not in cotx.columns or "protein_accession" not in ips.columns:
            return None

        # Get InterPro accessions from co-transcribed genes
        downstream_ids = set(cotx["downstream_orf"].to_list())
        cotx_domains = ips.filter(pl.col("protein_accession").is_in(list(downstream_ids)))

        if len(cotx_domains) == 0:
            return None

        # Collect InterPro accessions from co-transcribed gene annotations
        cotx_interpro_ids = set()
        if "interpro_accession" in cotx_domains.columns:
            cotx_interpro_ids = set(
                cotx_domains.filter(pl.col("interpro_accession") != "-")
                ["interpro_accession"].unique().to_list()
            )

        if not cotx_interpro_ids:
            return None

        # Filter binomial results to only co-transcribed gene domains
        id_col = "V13" if "V13" in df.columns else None
        if id_col is None:
            return df

        return df.filter(pl.col(id_col).is_in(list(cotx_interpro_ids)))

    @render.ui
    def binomial_bar_chart():
        df = cotx_binomial_data()
        if df is None or len(df) == 0:
            return ui.p("No binomial enrichment results for co-transcribed gene domains.")
        # Filter to significant domains
        if "p_adju" in df.columns:
            sig = df.filter(pl.col("p_adju") < 0.05)
        else:
            sig = df
        if len(sig) == 0:
            return ui.p("No significantly enriched domains found in co-transcribed genes.")
        # Sort by p_scaled descending
        if "p_scaled" in sig.columns:
            sig = sig.sort("p_scaled", descending=True).head(20)
        label_col = "V14" if "V14" in sig.columns else "V13" if "V13" in sig.columns else sig.columns[0]
        score_col = "p_scaled" if "p_scaled" in sig.columns else "p_adju"
        fig = px.bar(
            sig.to_pandas(), x=score_col, y=label_col, orientation="h",
            labels={score_col: "-10*log10(adj. p-value)", label_col: "Domain"},
            color_discrete_sequence=[CHART_COLORS["primary"]],
        )
        fig.update_layout(height=400, yaxis={"categoryorder": "total ascending"},
                          margin=dict(l=200), plot_bgcolor=CHART_COLORS["background"])
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.data_frame
    def binomial_table():
        df = cotx_binomial_data()
        if df is None or len(df) == 0:
            return None
        if "p_adju" in df.columns:
            df = df.filter(pl.col("p_adju") < 0.05)
        if "p_scaled" in df.columns:
            df = df.sort("p_scaled", descending=True)
        return render.DataTable(df.to_pandas(), filters=True, height="400px")

    # Domain Annotations tab outputs
    @render.ui
    def domain_frequency_plot():
        df = interproscan_data()
        if df is None or len(df) == 0:
            return ui.p("No InterProScan results available. Run pipeline with InterProScan enabled.")

        # Count domain frequencies by signature description
        if "signature_desc" not in df.columns:
            return ui.p("No domain annotations found")

        counts = (
            df.filter(df["signature_desc"] != "-")
            .group_by("signature_desc")
            .len()
            .sort("len", descending=True)
            .head(15)
        )

        if len(counts) == 0:
            return ui.p("No domain annotations found")

        fig = px.bar(
            counts.to_pandas(),
            x="len",
            y="signature_desc",
            orientation="h",
            labels={"len": "Count", "signature_desc": "Domain"},
            color_discrete_sequence=[CHART_COLORS["secondary"]],
        )
        fig.update_layout(
            height=400,
            yaxis={"categoryorder": "total ascending"},
            margin=dict(l=200),
            plot_bgcolor=CHART_COLORS["background"],
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.ui
    def analysis_type_plot():
        df = interproscan_data()
        if df is None or len(df) == 0:
            return ui.p("No InterProScan results available")

        if "analysis" not in df.columns:
            return ui.p("No analysis data found")

        counts = df.group_by("analysis").len().sort("len", descending=True)

        fig = px.pie(
            counts.to_pandas(),
            values="len",
            names="analysis",
            title="Annotations by Database",
            color_discrete_sequence=TEAL_PALETTE,
        )
        fig.update_layout(height=400)
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.ui
    def domain_binomial_plot():
        df = binomial_data()
        if df is None or len(df) == 0:
            return ui.p("No binomial enrichment results available. Run the binomial analysis pipeline first.")

        # Filter to significant domains
        if "p_adju" in df.columns:
            sig = df.filter(pl.col("p_adju") < 0.05)
        else:
            sig = df

        if len(sig) == 0:
            return ui.p("No significantly enriched domains found (adjusted p < 0.05).")

        # Sort by p_scaled descending and take top 20
        if "p_scaled" in sig.columns:
            sig = sig.sort("p_scaled", descending=True).head(20)

        label_col = "V14" if "V14" in sig.columns else "V13" if "V13" in sig.columns else sig.columns[0]
        score_col = "p_scaled" if "p_scaled" in sig.columns else "p_adju"

        fig = px.bar(
            sig.to_pandas(),
            x=score_col,
            y=label_col,
            orientation="h",
            labels={score_col: "-10*log10(adj. p-value)", label_col: "Domain"},
            color_discrete_sequence=[CHART_COLORS["secondary"]],
        )
        fig.update_layout(
            height=500,
            yaxis={"categoryorder": "total ascending"},
            margin=dict(l=250),
            plot_bgcolor=CHART_COLORS["background"],
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.data_frame
    def interproscan_table():
        df = interproscan_data()
        if df is None or len(df) == 0:
            return None

        # Select key columns for display
        display_cols = [
            col for col in ["protein_accession", "analysis", "signature_accession",
                           "signature_desc", "start", "stop", "interpro_accession",
                           "interpro_desc", "go_annotations"]
            if col in df.columns
        ]

        if not display_cols:
            return None

        result = df.select(display_cols)

        # Join with binomial data to add enrichment p-value
        binom = binomial_data()
        if (binom is not None and len(binom) > 0
                and "interpro_accession" in result.columns):
            id_col = "V13" if "V13" in binom.columns else None
            p_col = "p_adju" if "p_adju" in binom.columns else None
            if id_col and p_col:
                binom_lookup = binom.select([
                    pl.col(id_col).alias("interpro_accession"),
                    pl.col(p_col).alias("binomial_p_adjusted"),
                ]).unique(subset=["interpro_accession"])
                result = result.join(binom_lookup, on="interpro_accession", how="left")

        return render.DataTable(result.to_pandas(), filters=True, height="400px")

    # DefenseFinder tab outputs — category charts

    # Colour palette for the 8 fixed antiphage categories
    CATEGORY_COLORS = {
        "Abortive infection":       "#2d4a4a",
        "CBASS":                    "#3d6b6b",
        "Nucleic acid restriction": "#4a7d74",
        "CRISPR-Cas":               "#5a9090",
        "Toxin-antitoxin":          "#6b9e96",
        "Retrons":                  "#7aaeae",
        "tRNA degradation":         "#8fb3b3",
        "Unknown mechanism":        "#a8c4c4",
    }
    CATEGORY_ORDER = list(CATEGORY_COLORS.keys())

    @render.ui
    def defense_category_plot():
        df = defensefinder_systems()
        if df is None or len(df) == 0:
            return ui.p("No defense systems found. Run pipeline with DefenseFinder enabled.")
        if "antiphage_category" not in df.columns:
            return ui.p("Category data unavailable.")

        # Deduplicate on sys_id so each system counts once
        if "sys_id" in df.columns:
            deduped = df.unique(subset=["sys_id"])
        else:
            deduped = df

        counts = (
            deduped.group_by("antiphage_category")
            .len()
            .sort("len", descending=True)
        )

        pdf = counts.to_pandas()
        pdf["color"] = pdf["antiphage_category"].map(
            lambda c: CATEGORY_COLORS.get(c, "#a8c4c4")
        )

        fig = go.Figure(go.Bar(
            x=pdf["antiphage_category"],
            y=pdf["len"],
            marker_color=pdf["color"],
            hovertemplate="<b>%{x}</b><br>Systems: %{y}<extra></extra>",
        ))
        fig.update_layout(
            xaxis_title="Antiphage Category",
            yaxis_title="Number of Defense Systems",
            xaxis_tickangle=-30,
            height=420,
            margin=dict(b=120),
            plot_bgcolor=CHART_COLORS["background"],
            showlegend=False,
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.ui
    def defense_category_sunburst():
        df = defensefinder_systems()
        if df is None or len(df) == 0:
            return ui.p("No defense systems found.")
        if "antiphage_category" not in df.columns or "type" not in df.columns:
            return ui.p("Category data unavailable.")

        # Deduplicate on sys_id
        if "sys_id" in df.columns:
            deduped = df.unique(subset=["sys_id"])
        else:
            deduped = df

        counts = (
            deduped.group_by(["antiphage_category", "type"])
            .len()
            .sort("len", descending=True)
        )

        pdf = counts.to_pandas()
        fig = px.sunburst(
            pdf,
            path=["antiphage_category", "type"],
            values="len",
            color="antiphage_category",
            color_discrete_map=CATEGORY_COLORS,
            hover_data={"len": True},
        )
        fig.update_traces(
            hovertemplate="<b>%{label}</b><br>Count: %{value}<extra></extra>",
            insidetextorientation="radial",
        )
        fig.update_layout(
            height=420,
            margin=dict(t=10, b=10, l=10, r=10),
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.ui
    def defense_category_drilldown():
        df = defensefinder_systems()
        if df is None or len(df) == 0:
            return ui.p("No defense systems found.")
        if "antiphage_category" not in df.columns:
            return ui.p("Category data unavailable.")

        # Deduplicate on sys_id
        if "sys_id" in df.columns:
            deduped = df.unique(subset=["sys_id"])
        else:
            deduped = df

        selected = input.defense_category_select()
        if selected != "All":
            deduped = deduped.filter(pl.col("antiphage_category") == selected)

        if len(deduped) == 0:
            return ui.p(f"No systems found for category: {selected}")

        # Group by type (system family)
        counts = (
            deduped.group_by("type")
            .len()
            .sort("len", descending=True)
        )

        pdf = counts.to_pandas()

        # Pick a single colour for the selected category, cycle through palette otherwise
        if selected != "All":
            bar_color = CATEGORY_COLORS.get(selected, CHART_COLORS["secondary"])
            colors = [bar_color] * len(pdf)
        else:
            colors = [
                CATEGORY_COLORS.get(
                    deduped.filter(pl.col("type") == t)["antiphage_category"][0]
                    if len(deduped.filter(pl.col("type") == t)) > 0 else "",
                    CHART_COLORS["secondary"]
                )
                for t in pdf["type"]
            ]

        fig = go.Figure(go.Bar(
            y=pdf["type"],
            x=pdf["len"],
            orientation="h",
            marker_color=colors,
            hovertemplate="<b>%{y}</b><br>Systems: %{x}<extra></extra>",
        ))
        fig.update_layout(
            xaxis_title="Number of Defense Systems",
            yaxis_title="System Type",
            yaxis={"categoryorder": "total ascending"},
            height=max(350, len(pdf) * 24),
            margin=dict(l=160, r=20, t=20, b=40),
            plot_bgcolor=CHART_COLORS["background"],
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.ui
    def defense_system_type_plot():
        df = defensefinder_systems()
        if df is None or len(df) == 0:
            return ui.p("No defense systems found. Run pipeline with DefenseFinder enabled.")

        if "type" not in df.columns:
            return ui.p("No defense system type data found")

        counts = df.group_by("type").len().sort("len", descending=True)

        if len(counts) == 0:
            return ui.p("No defense systems found")

        fig = px.bar(
            counts.to_pandas(),
            x="type",
            y="len",
            labels={"type": "Defense System Type", "len": "Count"},
            title="Defense Systems by Type",
            color_discrete_sequence=[CHART_COLORS["primary"]],
        )
        fig.update_layout(
            xaxis_tickangle=-45,
            height=400,
            margin=dict(b=100),
            plot_bgcolor=CHART_COLORS["background"],
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.ui
    def defense_subtype_plot():
        df = defensefinder_systems()
        if df is None or len(df) == 0:
            return ui.p("No defense systems found")

        if "subtype" not in df.columns:
            return ui.p("No subtype data found")

        counts = df.group_by("subtype").len().sort("len", descending=True).head(15)

        if len(counts) == 0:
            return ui.p("No subtypes found")

        fig = px.bar(
            counts.to_pandas(),
            x="len",
            y="subtype",
            orientation="h",
            labels={"len": "Count", "subtype": "Subtype"},
            title="Top Defense System Subtypes",
            color_discrete_sequence=[CHART_COLORS["secondary"]],
        )
        fig.update_layout(
            height=400,
            yaxis={"categoryorder": "total ascending"},
            margin=dict(l=150),
            plot_bgcolor=CHART_COLORS["background"],
        )
        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))

    @render.data_frame
    def defense_genes_table():
        df = defensefinder_genes()
        if df is None or len(df) == 0:
            return None
        return render.DataTable(df.to_pandas(), filters=True, height="400px")

    # Download handler
    @render.download(filename="blast_results_filtered.csv")
    def download_blast():
        df = filtered_blast_data()
        if df is not None:
            return df.write_csv()
        return ""

    # ==================== Locus Viewer Tab ====================

    # Populate defense type selectize choices
    @reactive.effect
    def _update_defense_type_choices():
        df = defensefinder_genes()
        if df is not None and len(df) > 0 and "type" in df.columns:
            types = sorted(df["type"].unique().to_list())
            ui.update_selectize("locus_defense_type_filter", choices=types)

    # Reset filters button handler
    @reactive.effect
    @reactive.event(input.locus_reset_filters)
    def _reset_locus_filters():
        ui.update_slider("locus_defense_score_range", value=[0.0, 1.0])
        ui.update_slider("locus_binomial_pvalue", value=0.05)
        ui.update_checkbox("locus_has_cotranscribed", value=False)
        ui.update_checkbox("locus_has_defense", value=False)
        ui.update_checkbox("locus_has_mge", value=False)
        ui.update_selectize("locus_defense_type_filter", selected=[])

    # Select All Filtered button handler
    @reactive.effect
    @reactive.event(input.locus_select_all)
    def _select_all_contigs():
        contigs = filtered_contigs()
        if contigs:
            ui.update_selectize("locus_contig_select", selected=contigs)

    # Clear All button handler
    @reactive.effect
    @reactive.event(input.locus_clear_all)
    def _clear_all_contigs():
        ui.update_selectize("locus_contig_select", selected=[])

    # Selected count display
    @render.ui
    def locus_selected_count():
        selected = input.locus_contig_select()
        contigs = filtered_contigs()
        selected_list = list(selected) if isinstance(selected, (tuple, list)) else [selected] if selected else []
        selected_list = [c for c in selected_list if c and not c.startswith("No ")]
        total = len(contigs) if contigs else 0
        return ui.tags.small(
            f"{len(selected_list)} of {total} contigs selected",
            style="color: #8fb3b3; display: block; margin-bottom: 5px;"
        )

    @reactive.calc
    def filtered_contigs():
        """Apply all locus filters to produce a filtered contig list."""
        filter_type = input.locus_filter_type()
        orfs = orf_data()
        blast = blast_data()
        defense_genes_df = defensefinder_genes()
        interproscan = interproscan_data()
        scores = defense_scores_data()

        if orfs is None or len(orfs) == 0:
            return []

        # Step 1: Base contig list from radio button
        if filter_type == "hits":
            contigs = get_contigs_with_hits(orfs, blast)
        else:
            contigs = get_contig_choices(orfs)

        if not contigs:
            return []

        contig_set = set(contigs)

        # Step 2: Defense score range filter
        score_range = input.locus_defense_score_range()
        if scores is not None and len(scores) > 0 and (score_range[0] > 0.0 or score_range[1] < 1.0):
            scored = scores.filter(pl.col("defense_score").is_not_null())
            if len(scored) > 0 and "contig" in scored.columns:
                matching_contigs = set(
                    scored.filter(
                        (pl.col("defense_score") >= score_range[0]) &
                        (pl.col("defense_score") <= score_range[1])
                    )["contig"].unique().to_list()
                )
                contig_set &= matching_contigs

        # Step 2b: Binomial p-value filter (only when co-transcribed filter is active)
        if input.locus_has_cotranscribed():
            binom_threshold = input.locus_binomial_pvalue()
            binom = binomial_data()
            if binom is not None and len(binom) > 0 and "p_adju" in binom.columns:
                id_col = "V13" if "V13" in binom.columns else None
                if id_col is not None:
                    # Get significant domain IDs at the chosen threshold
                    sig_domains = set(
                        binom.filter(pl.col("p_adju") <= binom_threshold)
                        [id_col].to_list()
                    )
                    if sig_domains:
                        ips = interproscan_data()
                        if ips is not None and len(ips) > 0 and "interpro_accession" in ips.columns:
                            # Find ORFs that have significant domains
                            sig_orfs = set(
                                ips.filter(pl.col("interpro_accession").is_in(list(sig_domains)))
                                ["protein_accession"].unique().to_list()
                            )
                            # Map to contigs
                            if orfs is not None and sig_orfs:
                                sig_contigs = set(
                                    orfs.filter(pl.col("orf_id").is_in(list(sig_orfs)))
                                    ["contig"].unique().to_list()
                                )
                                contig_set &= sig_contigs

        # Step 3: Has co-transcribed genes filter
        if input.locus_has_cotranscribed():
            cotx = cotx_data()
            if cotx is not None and len(cotx) > 0 and "downstream_orf" in cotx.columns:
                # Get contigs that have downstream co-transcribed ORFs
                downstream_ids = set(cotx["downstream_orf"].to_list())
                cotx_contigs = set()
                for row in orfs.iter_rows(named=True):
                    if row["orf_id"] in downstream_ids:
                        cotx_contigs.add(row["contig"])
                contig_set &= cotx_contigs
            else:
                contig_set = set()

        # Step 4: Has defense systems filter
        if input.locus_has_defense():
            defense_contigs = set(get_contigs_with_defense(orfs, defense_genes_df))
            contig_set &= defense_contigs

        # Step 4: Has other MGEs filter
        if input.locus_has_mge():
            mge_contigs = set(get_contigs_with_mge(orfs, interproscan))
            contig_set &= mge_contigs

        # Step 5: Defense type filter
        selected_types = input.locus_defense_type_filter()
        if selected_types and len(selected_types) > 0 and defense_genes_df is not None and len(defense_genes_df) > 0:
            type_set = set(selected_types)
            # Find ORFs matching selected defense types
            matching_orf_ids = set()
            for row in defense_genes_df.iter_rows(named=True):
                if row.get("type", "") in type_set:
                    matching_orf_ids.add(row.get("hit_id", ""))
            # Map to contigs
            type_contigs = set()
            for row in orfs.iter_rows(named=True):
                if row["orf_id"] in matching_orf_ids:
                    type_contigs.add(row["contig"])
            contig_set &= type_contigs

        return sorted(list(contig_set))

    # Update contig selector based on filtered contigs
    @reactive.effect
    def _update_locus_contigs():
        contigs = filtered_contigs()

        if not contigs:
            filter_type = input.locus_filter_type()
            if filter_type == "hits":
                ui.update_selectize("locus_contig_select", choices=["No contigs with hits"])
            else:
                ui.update_selectize("locus_contig_select", choices=["No matching contigs"])
            return

        ui.update_selectize("locus_contig_select", choices=contigs, selected=[contigs[0]])

    @reactive.calc
    def selected_contig_orfs():
        """Get ORFs for the selected contig(s)."""
        orfs = orf_data()
        selected = input.locus_contig_select()

        if orfs is None or not selected:
            return None

        # Handle tuple/list from selectize multiple
        selected_list = list(selected) if isinstance(selected, (tuple, list)) else [selected]
        # Filter out placeholder values
        selected_list = [c for c in selected_list if c and not c.startswith("No ")]
        if not selected_list:
            return None

        return orfs.filter(pl.col("contig").is_in(selected_list))

    @render.text
    def locus_title():
        selected = input.locus_contig_select()
        if not selected:
            return "Locus Viewer"
        selected_list = list(selected) if isinstance(selected, (tuple, list)) else [selected]
        selected_list = [c for c in selected_list if c and not c.startswith("No ")]
        if not selected_list:
            return "Locus Viewer"
        if len(selected_list) == 1:
            return f"Locus View: {selected_list[0]}"
        return f"Locus View: {len(selected_list)} contigs selected"

    @render.ui
    def locus_color_legend():
        return ui.HTML(create_color_legend_html())

    @render.ui
    def locus_viewer_plot():
        orfs = orf_data()
        selected = input.locus_contig_select()

        if orfs is None or not selected:
            return ui.p("Select a contig to view its locus structure.")

        selected_list = list(selected) if isinstance(selected, (tuple, list)) else [selected]
        selected_list = [c for c in selected_list if c and not c.startswith("No ")]
        if not selected_list:
            return ui.p("Select a contig to view its locus structure.")

        # Get annotation data
        mapping = hit_to_contig_mapping()
        interproscan = interproscan_data()
        defense_genes = defensefinder_genes()
        cotx = cotx_data()

        # Get downstream ORF IDs
        downstream_orfs = None
        if cotx is not None and len(cotx) > 0:
            downstream_orfs = cotx["downstream_orf"].to_list()

        # Generate one plot per contig, stack vertically
        html_parts = []
        for i, contig_id in enumerate(selected_list):
            fig = create_locus_figure(
                orfs=orfs,
                contig_id=contig_id,
                hit_to_contig_mapping=mapping,
                interproscan=interproscan,
                defensefinder_genes=defense_genes,
                downstream_orfs=downstream_orfs,
                height=350,
            )
            html_parts.append(f'<div style="margin-bottom: 10px;">')
            html_parts.append(fig.to_html(include_plotlyjs=False, full_html=False))
            html_parts.append('</div>')

        return ui.HTML(
            f'<div style="max-height: 80vh; overflow-y: auto;">{"".join(html_parts)}</div>'
        )

    @render.data_frame
    def locus_orf_table():
        contig_orfs = selected_contig_orfs()
        if contig_orfs is None or len(contig_orfs) == 0:
            return None

        return render.DataTable(
            contig_orfs.to_pandas(),
            filters=True,
            height="300px"
        )

    # ==================== Defense Locus Viewer ====================

    # Populate defense type filter choices in DefenseFinder tab
    @reactive.effect
    def _update_defense_type_filter():
        df = defensefinder_genes()
        if df is not None and len(df) > 0 and "type" in df.columns:
            types = sorted(df["type"].unique().to_list())
            ui.update_selectize("defense_type_filter", choices=types)

    @reactive.effect
    def _update_defense_contigs():
        orfs = orf_data()
        defense_genes_df = defensefinder_genes()

        contigs = get_contigs_with_defense(orfs, defense_genes_df)

        # Apply defense type filter if selected
        selected_types = input.defense_type_filter()
        if selected_types and len(selected_types) > 0 and defense_genes_df is not None and len(defense_genes_df) > 0:
            type_set = set(selected_types)
            # Find ORFs matching selected defense types
            matching_orf_ids = set()
            for row in defense_genes_df.iter_rows(named=True):
                if row.get("type", "") in type_set:
                    matching_orf_ids.add(row.get("hit_id", ""))
            # Map to contigs
            if orfs is not None:
                type_contigs = set()
                for row in orfs.iter_rows(named=True):
                    if row["orf_id"] in matching_orf_ids:
                        type_contigs.add(row["contig"])
                contigs = sorted(list(set(contigs) & type_contigs))

        if contigs:
            ui.update_select("defense_contig_select", choices=contigs, selected=contigs[0])
        else:
            ui.update_select("defense_contig_select", choices=["No contigs with defense systems"])

    @render.ui
    def defense_system_info():
        contig_id = input.defense_contig_select()
        defense_genes = defensefinder_genes()

        if not contig_id or contig_id.startswith("No ") or defense_genes is None:
            return ui.p("No defense systems found.")

        # Count defense genes in this contig
        orfs = orf_data()
        if orfs is None:
            return ui.p("No data available.")

        defense_gene_ids = set(defense_genes["hit_id"].to_list())
        contig_orfs = orfs.filter(pl.col("contig") == contig_id)

        defense_count = sum(1 for row in contig_orfs.iter_rows(named=True)
                          if row["orf_id"] in defense_gene_ids)

        return ui.tags.div(
            ui.tags.p(f"Defense genes in contig: {defense_count}"),
            ui.tags.p(f"Total ORFs: {len(contig_orfs)}"),
        )

    @render.ui
    def defense_locus_plot():
        orfs = orf_data()
        contig_id = input.defense_contig_select()
        defense_genes = defensefinder_genes()
        defense_systems = defensefinder_systems()
        mapping = hit_to_contig_mapping()

        if orfs is None or not contig_id or contig_id.startswith("No "):
            return ui.p("Select a contig to view defense system locus.")

        if defense_genes is None or len(defense_genes) == 0:
            return ui.p("No defense genes found. Run DefenseFinder to annotate defense systems.")

        fig = create_defense_locus_figure(
            orfs=orfs,
            contig_id=contig_id,
            defensefinder_genes=defense_genes,
            defensefinder_systems=defense_systems,
            hit_to_contig_mapping=mapping,
            height=350,
        )

        return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))


app = App(app_ui, server)
