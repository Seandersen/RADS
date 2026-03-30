#!/usr/bin/env python3
"""
RADS Pipeline — Standalone HTML Report Generator
=================================================
Generates a self-contained interactive HTML report from pipeline results.
Requires an internet connection to load Plotly.js from CDN when viewing.

Usage:
    python workflow/scripts/generate_report.py \\
        --results results/SAMPLE_NAME \\
        --output  results/SAMPLE_NAME/report.html

Snakemake usage:
    See workflow/rules/report.smk
"""

import argparse
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import plotly.graph_objects as go
import plotly.express as px

# ---------------------------------------------------------------------------
# Color scheme (matches RADS PyShiny dashboard)
# ---------------------------------------------------------------------------
BG = "#f5f8f8"
TEAL = [
    "#2d4a4a", "#3d5a5a", "#4a6670", "#5a7a7a", "#6b9090",
    "#7a9e9e", "#8fb3b3", "#a8c4c4", "#b5cece", "#c5dada",
]
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
_BASE_LAYOUT = dict(
    plot_bgcolor=BG,
    paper_bgcolor="white",
    font=dict(family="system-ui, -apple-system, sans-serif", color="#2d4a4a", size=12),
    hoverlabel=dict(bgcolor="white", font_size=12),
)


def _layout(**overrides):
    """Merge _BASE_LAYOUT with per-chart overrides (including margin)."""
    return {**_BASE_LAYOUT, "margin": dict(t=50, b=60, l=70, r=30), **overrides}

# ---------------------------------------------------------------------------
# Antiphage category mapping  (mirrors dashboard/utils/data_loader.py)
# ---------------------------------------------------------------------------
_CAT_EXACT = {
    "Cas": "CRISPR-Cas", "CRISPR-Cas": "CRISPR-Cas",
    "CBASS": "CBASS", "Thoeris": "CBASS", "Pycsar": "CBASS",
    "RM": "Nucleic acid restriction", "BREX": "Nucleic acid restriction",
    "DISARM": "Nucleic acid restriction", "Dnd": "Nucleic acid restriction",
    "Dpd": "Nucleic acid restriction", "Wadjet": "Nucleic acid restriction",
    "Zorya": "Nucleic acid restriction", "Shedu": "Nucleic acid restriction",
    "NixI": "Nucleic acid restriction", "Nhi": "Nucleic acid restriction",
    "pAgo": "Nucleic acid restriction", "RADAR": "Nucleic acid restriction",
    "SspBCDE": "Nucleic acid restriction",
    "Retron": "Retrons",
    "PrrC": "tRNA degradation", "RloC": "tRNA degradation", "CapRel": "tRNA degradation",
    "DRT": "Toxin-antitoxin", "MazEF": "Toxin-antitoxin", "RexAB": "Toxin-antitoxin",
    "RnlAB": "Toxin-antitoxin", "SoFIC": "Toxin-antitoxin",
    "RosmerTA": "Toxin-antitoxin", "ShosTA": "Toxin-antitoxin",
    "Gabija": "Abortive infection", "Druantia": "Abortive infection",
    "Hachiman": "Abortive infection", "Lamassu-Fam": "Abortive infection",
    "Lamassu": "Abortive infection", "PARIS": "Abortive infection",
    "Paris": "Abortive infection", "Avs": "Abortive infection",
    "BstA": "Abortive infection", "Kiwa": "Abortive infection",
    "Lit": "Abortive infection", "Shango": "Abortive infection",
    "JukAB": "Abortive infection", "Septu": "Abortive infection",
    "SEFIR": "Abortive infection", "GasderMIN": "Abortive infection",
    "SpbK": "Abortive infection", "Stk2": "Abortive infection",
    "Pif": "Abortive infection", "DdmDE": "Abortive infection",
    "Dsr": "Abortive infection", "Viperin": "Abortive infection",
    "DarTG": "Abortive infection", "Detocs": "Abortive infection",
    "Borvo": "Abortive infection", "Abi": "Abortive infection",
    "Menshen": "Unknown mechanism", "Mokosh": "Unknown mechanism",
    "Aditi": "Unknown mechanism", "Dazbog": "Unknown mechanism",
    "Tiamat": "Unknown mechanism", "Dodola": "Unknown mechanism",
    "Eleos": "Unknown mechanism", "NLR": "Unknown mechanism",
    "Azaca": "Unknown mechanism", "Bunzi": "Unknown mechanism",
    "Uzume": "Unknown mechanism", "ISG15-like": "Unknown mechanism",
    "MADS": "Unknown mechanism",
}
_CAT_PREFIX = [
    ("Cas_Type_", "CRISPR-Cas"), ("CBASS_Type_", "CBASS"),
    ("RM_Type_", "Nucleic acid restriction"), ("Retron_", "Retrons"),
    ("Abi", "Abortive infection"), ("Gao_", "Unknown mechanism"),
    ("PD-Lambda-", "Unknown mechanism"), ("PD-T7-", "Unknown mechanism"),
    ("PD-T4-", "Abortive infection"), ("FS_", "Nucleic acid restriction"),
]


def _category(system_type: str, subtype: str = "") -> str:
    for field in (system_type, subtype):
        if not field or field == "Unknown":
            continue
        if field in _CAT_EXACT:
            return _CAT_EXACT[field]
        for prefix, cat in _CAT_PREFIX:
            if field.startswith(prefix):
                return cat
    return "Unknown mechanism"


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_metrics(results_dir: Path) -> dict:
    f = results_dir / "metrics" / "pipeline_metrics.json"
    return json.loads(f.read_text()) if f.exists() else {}


def load_blast(results_dir: Path) -> list[dict]:
    f = results_dir / "blast_results" / "master_blast.txt"
    if not f.exists():
        return []
    rows = []
    with open(f) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) == len(header):
                rows.append(dict(zip(header, parts)))
    return rows


def load_defensefinder(results_dir: Path) -> list[dict]:
    f = results_dir / "defensefinder" / "defense_finder_systems.tsv"
    if not f.exists():
        return []
    rows = []
    seen_sys = set()
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("replicon\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            model_fqn = parts[4] if len(parts) > 4 else ""
            fqn_parts = model_fqn.split("/")
            subtype = fqn_parts[-1] if fqn_parts else "Unknown"
            system_type = fqn_parts[-2] if len(fqn_parts) >= 2 else subtype
            sys_id = parts[5] if len(parts) > 5 else ""
            if sys_id in seen_sys:
                continue
            seen_sys.add(sys_id)
            rows.append({
                "sys_id": sys_id,
                "type": system_type,
                "subtype": subtype,
                "category": _category(system_type, subtype),
                "replicon": parts[0],
            })
    return rows


def load_cotranscription(results_dir: Path) -> list[dict]:
    f = results_dir / "cotranscription" / "cotranscribed_details.txt"
    if not f.exists():
        return []
    rows = []
    with open(f) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) == len(header):
                rows.append(dict(zip(header, parts)))
    return rows


def load_interproscan(results_dir: Path) -> list[dict]:
    f = results_dir / "interproscan_results.tsv"
    if not f.exists():
        return []
    rows = []
    with open(f) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 13:
                rows.append({
                    "protein": parts[0],
                    "analysis": parts[3],
                    "sig_acc": parts[4],
                    "sig_desc": parts[5],
                    "ipr_acc": parts[11],
                    "ipr_desc": parts[12],
                })
    return rows


def load_organism_names(results_dir: Path) -> dict[str, str]:
    """Parse genus from first FASTA header of each staged genome."""
    genomes_dir = results_dir / "genomes"
    organism_map = {}
    if not genomes_dir.exists():
        return organism_map
    for fna in genomes_dir.glob("*.fna"):
        genome_id = fna.stem
        try:
            with open(fna) as fh:
                first = fh.readline()
            if first.startswith(">"):
                parts = first[1:].split(None, 1)
                if len(parts) >= 2:
                    words = parts[1].strip().split()
                    genus = words[0] if words else "Unknown"
                    # Strip non-alpha characters (plasmid prefixes, etc.)
                    genus = re.sub(r"[^A-Za-z]", "", genus) or "Unknown"
                    organism_map[genome_id] = genus
                else:
                    organism_map[genome_id] = "Unknown"
        except Exception:
            organism_map[genome_id] = "Unknown"
    return organism_map


# ---------------------------------------------------------------------------
# Chart builders
# ---------------------------------------------------------------------------

def _fig_html(fig, first: bool = False) -> str:
    return fig.to_html(
        include_plotlyjs="cdn" if first else False,
        full_html=False,
        config={"displayModeBar": True, "displaylogo": False,
                "modeBarButtonsToRemove": ["select2d", "lasso2d"]},
    )


def chart_organism_breakdown(organism_map: dict) -> tuple[str, str]:
    if not organism_map:
        return "Organism Breakdown", "<p class='no-data'>No genome data available.</p>"
    counts = Counter(organism_map.values())
    top = counts.most_common(20)
    labels, vals = zip(*top)
    colors = (TEAL * 5)[:len(labels)]
    fig = go.Figure(go.Bar(
        y=list(labels), x=list(vals), orientation="h",
        marker_color=colors,
        hovertemplate="<b>%{y}</b><br>Genomes: %{x}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=40, l=160, r=30)),
        title="Input Genomes by Genus (top 20)",
        xaxis_title="Number of Genomes",
        yaxis=dict(categoryorder="total ascending"),
        height=max(350, len(labels) * 22),
    )
    note = (f"<p class='chart-note'>Genus parsed from NCBI FASTA headers. "
            f"{len(organism_map):,} total genomes, {len(counts):,} unique genera shown.</p>")
    return "Input Genomes by Genus", _fig_html(fig, first=True) + note


def chart_blast_identity(blast: list[dict]) -> tuple[str, str]:
    if not blast:
        return "BLAST Identity Distribution", "<p class='no-data'>No BLAST results available.</p>"
    pidents = [float(r["pident"]) for r in blast if "pident" in r]
    fig = go.Figure(go.Histogram(
        x=pidents, nbinsx=30,
        marker_color=TEAL[2], marker_line_color=TEAL[0], marker_line_width=0.5,
        hovertemplate="Identity: %{x:.1f}%<br>Count: %{y}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(),
        title="BLAST Hit Identity Distribution",
        xaxis_title="% Identity",
        yaxis_title="Number of Hits",
        height=340,
    )
    return "BLAST Identity Distribution", _fig_html(fig)


def chart_blast_hits_per_genome(blast: list[dict]) -> tuple[str, str]:
    if not blast:
        return "Hits per Genome", "<p class='no-data'>No BLAST results available.</p>"
    counts = Counter(r["genome"] for r in blast if "genome" in r)
    sorted_counts = sorted(counts.values(), reverse=True)
    fig = go.Figure(go.Histogram(
        x=sorted_counts, nbinsx=20,
        marker_color=TEAL[3], marker_line_color=TEAL[0], marker_line_width=0.5,
        hovertemplate="Hits per genome: %{x}<br>Genomes: %{y}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(),
        title="Distribution of BLAST Hits per Genome",
        xaxis_title="Number of Hits",
        yaxis_title="Number of Genomes",
        height=340,
    )
    note = (f"<p class='chart-note'>{len(counts):,} genomes with at least one hit. "
            f"Max {max(sorted_counts)} hits in a single genome.</p>")
    return "Hits per Genome", _fig_html(fig) + note


def chart_defense_categories(defense: list[dict]) -> tuple[str, str]:
    if not defense:
        return "Defense System Categories", "<p class='no-data'>No DefenseFinder results available.</p>"
    counts = Counter(r["category"] for r in defense)
    ordered = sorted(counts.items(), key=lambda x: -x[1])
    labels, vals = zip(*ordered)
    colors = [CATEGORY_COLORS.get(l, TEAL[4]) for l in labels]
    fig = go.Figure(go.Bar(
        x=list(labels), y=list(vals),
        marker_color=colors,
        hovertemplate="<b>%{x}</b><br>Systems: %{y}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=110, l=70, r=30)),
        title="Defense Systems by Antiphage Category",
        xaxis_title="Category",
        yaxis_title="Number of Systems",
        xaxis_tickangle=-25,
        height=380,
    )
    return "Defense System Categories", _fig_html(fig)


def chart_defense_sunburst(defense: list[dict]) -> tuple[str, str]:
    if not defense:
        return "Category × System Type", "<p class='no-data'>No DefenseFinder results available.</p>"
    # Build (category, type) counts
    ct_counts: dict[tuple, int] = Counter(
        (r["category"], r["type"]) for r in defense
    )
    parents, labels, values, colors = [], [], [], []
    cat_totals: dict[str, int] = Counter(r["category"] for r in defense)

    # Add category nodes
    for cat, total in cat_totals.items():
        parents.append("")
        labels.append(cat)
        values.append(total)
        colors.append(CATEGORY_COLORS.get(cat, TEAL[4]))

    # Add type nodes
    for (cat, typ), cnt in ct_counts.items():
        parents.append(cat)
        labels.append(typ)
        values.append(cnt)
        colors.append(CATEGORY_COLORS.get(cat, TEAL[4]))

    fig = go.Figure(go.Sunburst(
        parents=parents, labels=labels, values=values,
        marker=dict(colors=colors),
        hovertemplate="<b>%{label}</b><br>Systems: %{value}<extra></extra>",
        insidetextorientation="radial",
        branchvalues="total",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=10, l=10, r=10)),
        title="Category × System Type (click to expand)",
        height=440,
    )
    return "Category × System Type", _fig_html(fig)


def chart_defense_top_types(defense: list[dict]) -> tuple[str, str]:
    if not defense:
        return "Top System Types", "<p class='no-data'>No DefenseFinder results available.</p>"
    counts = Counter(r["type"] for r in defense)
    top = counts.most_common(15)
    labels, vals = zip(*top)
    colors = [
        CATEGORY_COLORS.get(
            next((r["category"] for r in defense if r["type"] == lbl), "Unknown mechanism"),
            TEAL[4]
        )
        for lbl in labels
    ]
    fig = go.Figure(go.Bar(
        y=list(labels), x=list(vals), orientation="h",
        marker_color=colors,
        hovertemplate="<b>%{y}</b><br>Systems: %{x}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=40, l=180, r=30)),
        title="Top 15 Defense System Types",
        xaxis_title="Number of Systems",
        yaxis=dict(categoryorder="total ascending"),
        height=420,
    )
    return "Top System Types", _fig_html(fig)


def chart_cotx_distance(cotx: list[dict]) -> tuple[str, str]:
    if not cotx:
        return "Co-transcription Gap Distance", "<p class='no-data'>No co-transcription results available.</p>"
    try:
        gaps = [int(r["gap_bp"]) for r in cotx if "gap_bp" in r]
    except Exception:
        gaps = []
    if not gaps:
        return "Co-transcription Gap Distance", "<p class='no-data'>Gap data unavailable.</p>"
    fig = go.Figure(go.Histogram(
        x=gaps, nbinsx=25,
        marker_color=TEAL[1], marker_line_color=TEAL[0], marker_line_width=0.5,
        hovertemplate="Gap: %{x} bp<br>Pairs: %{y}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(),
        title="Gap Distance Between Query Hit and Co-transcribed ORF",
        xaxis_title="Gap (bp)",
        yaxis_title="Number of Co-transcribed Pairs",
        height=320,
    )
    note = f"<p class='chart-note'>{len(gaps):,} co-transcribed pairs identified.</p>"
    return "Co-transcription Gap Distance", _fig_html(fig) + note


def chart_top_domains(ips: list[dict]) -> tuple[str, str]:
    if not ips:
        return "Top InterProScan Domains", "<p class='no-data'>No InterProScan results available.</p>"
    descs = [r["ipr_desc"] for r in ips if r.get("ipr_desc") and r["ipr_desc"] != "-"]
    if not descs:
        descs = [r["sig_desc"] for r in ips if r.get("sig_desc") and r["sig_desc"] != "-"]
    counts = Counter(descs)
    top = counts.most_common(20)
    if not top:
        return "Top InterProScan Domains", "<p class='no-data'>No domain descriptions found.</p>"
    labels, vals = zip(*top)
    colors = (TEAL * 4)[:len(labels)]
    fig = go.Figure(go.Bar(
        y=list(labels), x=list(vals), orientation="h",
        marker_color=colors,
        hovertemplate="<b>%{y}</b><br>Hits: %{x}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=40, l=330, r=30)),
        title="Top 20 InterPro Domains in Co-transcribed ORFs",
        xaxis_title="Number of Hits",
        yaxis=dict(categoryorder="total ascending"),
        height=520,
    )
    return "Top 20 InterPro Domains", _fig_html(fig)


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

def metric_card(label: str, value: str, subtitle: str = "") -> str:
    sub = f"<div class='card-sub'>{subtitle}</div>" if subtitle else ""
    return f"""
    <div class='metric-card'>
        <div class='card-label'>{label}</div>
        <div class='card-value'>{value}</div>
        {sub}
    </div>"""


def section(title: str, section_id: str, content: str) -> str:
    return f"""
    <section id='{section_id}'>
        <h2 class='section-title'>{title}</h2>
        <div class='section-body'>{content}</div>
    </section>"""


def chart_card(title: str, html: str) -> str:
    return f"""
    <div class='chart-card'>
        <div class='chart-inner'>{html}</div>
    </div>"""


def two_col(left: str, right: str) -> str:
    return f"<div class='two-col'>{left}{right}</div>"


CSS = """
:root {
    --teal-dark:  #2d4a4a;
    --teal-mid:   #4a6670;
    --teal-light: #8fb3b3;
    --teal-pale:  #c5dada;
    --bg:         #f5f8f8;
    --white:      #ffffff;
    --text:       #2d3a3a;
    --border:     #d0e0e0;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: system-ui, -apple-system, 'Segoe UI', sans-serif;
    background: var(--bg);
    color: var(--text);
    font-size: 14px;
    line-height: 1.5;
}
/* ---------- Header ---------- */
header {
    background: var(--teal-dark);
    color: white;
    padding: 0;
    position: sticky;
    top: 0;
    z-index: 100;
    box-shadow: 0 2px 8px rgba(0,0,0,0.2);
}
.header-inner {
    max-width: 1400px;
    margin: 0 auto;
    padding: 14px 24px;
    display: flex;
    align-items: center;
    justify-content: space-between;
}
.header-title { font-size: 1.35rem; font-weight: 700; letter-spacing: 0.5px; }
.header-meta  { font-size: 0.8rem; color: var(--teal-pale); text-align: right; }
/* ---------- Nav ---------- */
nav {
    background: var(--teal-mid);
    overflow-x: auto;
    white-space: nowrap;
}
.nav-inner {
    max-width: 1400px;
    margin: 0 auto;
    padding: 0 24px;
}
nav a {
    display: inline-block;
    color: #d0e8e8;
    text-decoration: none;
    padding: 10px 16px;
    font-size: 0.82rem;
    font-weight: 500;
    transition: background 0.15s, color 0.15s;
}
nav a:hover { background: rgba(255,255,255,0.15); color: white; }
/* ---------- Main layout ---------- */
main {
    max-width: 1400px;
    margin: 0 auto;
    padding: 32px 24px 64px;
}
/* ---------- Metric cards ---------- */
.metrics-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(170px, 1fr));
    gap: 14px;
    margin-bottom: 36px;
}
.metric-card {
    background: var(--teal-dark);
    color: white;
    border-radius: 10px;
    padding: 18px 16px;
    text-align: center;
}
.card-label { font-size: 0.72rem; color: var(--teal-pale); text-transform: uppercase; letter-spacing: 0.6px; margin-bottom: 8px; }
.card-value { font-size: 1.8rem; font-weight: 700; line-height: 1; }
.card-sub   { font-size: 0.72rem; color: var(--teal-light); margin-top: 4px; }
/* ---------- Sections ---------- */
section { margin-bottom: 48px; }
.section-title {
    font-size: 1.1rem;
    font-weight: 700;
    color: var(--teal-dark);
    border-left: 4px solid var(--teal-light);
    padding-left: 12px;
    margin-bottom: 18px;
}
.section-body { }
/* ---------- Chart cards ---------- */
.chart-card {
    background: var(--white);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 16px;
    margin-bottom: 18px;
    box-shadow: 0 1px 4px rgba(45,74,74,0.07);
}
.chart-inner { width: 100%; overflow-x: auto; }
.two-col {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 18px;
}
@media (max-width: 900px) { .two-col { grid-template-columns: 1fr; } }
/* ---------- Misc ---------- */
.chart-note {
    font-size: 0.78rem;
    color: #6b8888;
    margin-top: 6px;
    padding-left: 4px;
}
.no-data {
    color: #7a9090;
    font-style: italic;
    padding: 24px;
    text-align: center;
}
footer {
    text-align: center;
    padding: 24px;
    color: #7a9090;
    font-size: 0.78rem;
    border-top: 1px solid var(--border);
}
"""


def build_html(
    sample_name: str,
    generated_at: str,
    metrics: dict,
    organism_map: dict,
    blast: list[dict],
    defense: list[dict],
    cotx: list[dict],
    ips: list[dict],
) -> str:

    # ---- metric cards ----
    def fmt(val, decimals=0):
        if val is None:
            return "—"
        if isinstance(val, float):
            return f"{val:,.{decimals}f}"
        return f"{val:,}" if isinstance(val, int) else str(val)

    cards_html = "".join([
        metric_card("Total Genomes",       fmt(metrics.get("total_genomes")),               "input"),
        metric_card("Total Input",         fmt(metrics.get("total_input_mb"), 1) + " Mb",   "sequenced"),
        metric_card("BLAST Hits",          fmt(metrics.get("blast_hits")),                  "query matches"),
        metric_card("Hits per Mb",         fmt(metrics.get("hits_per_mb"), 3),               ""),
        metric_card("Defense Systems",     fmt(metrics.get("defense_systems")),              "DefenseFinder"),
        metric_card("Contigs Analyzed",    fmt(metrics.get("contigs_analyzed")),             "flanking regions"),
        metric_card("Co-transcribed Pairs",fmt(len(cotx)),                                   ""),
        metric_card("IPS Annotations",     fmt(len(ips)),                                    "domain hits"),
    ])

    # ---- charts (first=True marks where Plotly CDN is injected) ----
    first_chart = True

    def get_chart(fn, *args):
        nonlocal first_chart
        title, html = fn(*args)
        if "<p class='no-data'>" not in html and first_chart:
            first_chart = False
        return title, html

    org_title,    org_html    = get_chart(chart_organism_breakdown, organism_map)
    bi_title,     bi_html     = get_chart(chart_blast_identity,     blast)
    bhg_title,    bhg_html    = get_chart(chart_blast_hits_per_genome, blast)
    dc_title,     dc_html     = get_chart(chart_defense_categories, defense)
    ds_title,     ds_html     = get_chart(chart_defense_sunburst,   defense)
    dt_title,     dt_html     = get_chart(chart_defense_top_types,  defense)
    cx_title,     cx_html     = get_chart(chart_cotx_distance,      cotx)
    ips_title,    ips_html    = get_chart(chart_top_domains,         ips)

    # Note: Plotly CDN injected by organism chart (first real chart)
    # All subsequent charts have include_plotlyjs=False

    query_name = metrics.get("query_name", "Unknown")

    content = f"""
    <div class='metrics-grid'>{cards_html}</div>

    {section("Input Genomes", "genomes",
        chart_card(org_title, org_html)
    )}

    {section("BLAST Results", "blast",
        two_col(
            chart_card(bi_title,  bi_html),
            chart_card(bhg_title, bhg_html),
        )
    )}

    {section("Defense Systems (DefenseFinder)", "defense",
        chart_card(dc_title, dc_html) +
        two_col(
            chart_card(ds_title, ds_html),
            chart_card(dt_title, dt_html),
        )
    )}

    {section("Co-transcription Analysis", "cotranscription",
        chart_card(cx_title, cx_html)
    )}

    {section("InterProScan Domain Annotation", "interproscan",
        chart_card(ips_title, ips_html)
    )}
    """

    nav_links = "".join(
        f"<a href='#{s}'>{label}</a>"
        for s, label in [
            ("genomes",        "Input Genomes"),
            ("blast",          "BLAST Results"),
            ("defense",        "Defense Systems"),
            ("cotranscription","Co-transcription"),
            ("interproscan",   "InterProScan"),
        ]
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>RADS Report — {sample_name}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>{CSS}</style>
</head>
<body>
<header>
  <div class='header-inner'>
    <div>
      <div class='header-title'>RADS — Recombinase Associated Defense Search</div>
      <div style='font-size:0.8rem;color:var(--teal-pale);margin-top:2px;'>
        Sample: <b>{sample_name}</b> &nbsp;|&nbsp; Query: <b>{query_name}</b>
      </div>
    </div>
    <div class='header-meta'>
      Generated {generated_at}<br>
      defense-finder-models v2.0+
    </div>
  </div>
</header>
<nav><div class='nav-inner'>{nav_links}</div></nav>
<main>{content}</main>
<footer>
  RADS Pipeline Report &nbsp;·&nbsp; Generated {generated_at} &nbsp;·&nbsp;
  <a href='https://github.com/Seandersen/RADS' style='color:#8fb3b3;'>github.com/Seandersen/RADS</a>
</footer>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate RADS HTML report")
    parser.add_argument("--results", required=True, help="Path to sample results directory")
    parser.add_argument("--output",  required=True, help="Output HTML file path")
    args = parser.parse_args()

    results_dir = Path(args.results)
    out_path    = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading results from: {results_dir}")
    metrics      = load_metrics(results_dir)
    blast        = load_blast(results_dir)
    defense      = load_defensefinder(results_dir)
    cotx         = load_cotranscription(results_dir)
    ips          = load_interproscan(results_dir)
    organism_map = load_organism_names(results_dir)

    sample_name  = metrics.get("sample_name", results_dir.name)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"  Genomes:        {len(organism_map):,}")
    print(f"  BLAST hits:     {len(blast):,}")
    print(f"  Defense systems:{len(defense):,}")
    print(f"  Co-tx pairs:    {len(cotx):,}")
    print(f"  IPS annotations:{len(ips):,}")
    print("Building HTML...")

    html = build_html(
        sample_name=sample_name,
        generated_at=generated_at,
        metrics=metrics,
        organism_map=organism_map,
        blast=blast,
        defense=defense,
        cotx=cotx,
        ips=ips,
    )

    out_path.write_text(html, encoding="utf-8")
    print(f"Report written → {out_path}")


if __name__ == "__main__":
    main()
