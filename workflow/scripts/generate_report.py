#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RADS Pipeline — Standalone HTML Report Generator
=================================================
Generates a self-contained interactive HTML report from pipeline results.
Requires an internet connection to load Plotly.js from CDN when viewing.

Usage:
    python workflow/scripts/generate_report.py \\
        --results results/SAMPLE_NAME \\
        --output  results/SAMPLE_NAME/report.html \\
        [--include-locus-viewer]

Snakemake usage:
    See workflow/rules/report.smk
"""

import argparse
import csv
import json
import math
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
    "Anti-defense":             "#9bbfbf",
    "Unknown mechanism":        "#a8c4c4",
}

# Locus viewer colors (ported from dashboard/utils/locus_viewer.py)
QUERY_HIT_COLOR        = "#2d3d3d"
DOWNSTREAM_COLOR       = "#5a8a7a"
DEFENSE_COLOR          = "#6a9eae"
MGE_COLOR              = "#4a5a6a"
IPS_OTHER_COLOR        = "#b0b0b0"
NO_ANNOTATION_COLOR    = "#e8e8e8"

MGE_PFAM_ACCESSIONS = {
    "PF00665", "PF01609", "PF01527", "PF01526", "PF01610",
    "PF02371", "PF02914", "PF03400", "PF05598", "PF00589",
    "PF13009", "PF13276", "PF13358", "PF03389", "PF01076",
}
MGE_KEYWORDS = [
    "transposase", "integrase", "recombinase", "conjugal transfer",
    "mobilization", "insertion sequence", "is element", "conjugative",
]

_BASE_LAYOUT = dict(
    plot_bgcolor=BG,
    paper_bgcolor="white",
    font=dict(family="system-ui, -apple-system, sans-serif", color="#2d4a4a", size=12),
    hoverlabel=dict(bgcolor="white", font_size=12),
)


def _layout(**overrides):
    return {**_BASE_LAYOUT, "margin": dict(t=50, b=60, l=70, r=30), **overrides}


# ---------------------------------------------------------------------------
# Antiphage category mapping
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
    "RosmerTA": "Toxin-antitoxin", "ShosTA": "Toxin-antitoxin", "SanaTA": "Toxin-antitoxin",
    "Anti_RM": "Anti-defense",
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


def load_blast(results_dir: Path) -> list:
    # Columns: query_id, subject_id, length, nident, pident, evalue, genome
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


def load_defensefinder(results_dir: Path) -> list:
    """Load defense genes from defense_finder_genes.tsv, deduplicated by hit_id.

    Matches the counting method used by the figures script (defense_plots.py):
    one entry per unique gene hit, using the genes TSV (not systems TSV).
    """
    f = results_dir / "defensefinder" / "defense_finder_genes.tsv"
    if not f.exists():
        return []
    rows = []
    seen_hits = set()
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("replicon\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            hit_id = parts[1] if len(parts) > 1 else ""
            if hit_id == "hit_id" or hit_id in seen_hits:
                continue
            seen_hits.add(hit_id)
            model_fqn   = parts[4] if len(parts) > 4 else ""
            sys_id      = parts[5] if len(parts) > 5 else ""
            fqn_parts   = model_fqn.split("/")
            subtype     = fqn_parts[-1] if fqn_parts else "Unknown"
            system_type = fqn_parts[-2] if len(fqn_parts) >= 2 else subtype
            rows.append({
                "sys_id":   sys_id,
                "type":     system_type,
                "subtype":  subtype,
                "category": _category(system_type, subtype),
                "replicon": hit_id,
            })
    return rows


def load_defensefinder_genes(results_dir: Path) -> list:
    """Load defense_finder_genes.tsv (gene-level, not deduplicated)."""
    f = results_dir / "defensefinder" / "defense_finder_genes.tsv"
    if not f.exists():
        return []
    rows = []
    header = None
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            if len(parts) < len(header):
                continue
            row = dict(zip(header, parts))
            rows.append(row)
    return rows


def load_cotranscription(results_dir: Path) -> list:
    # Columns: blast_hit_id, hit_contig_orf, downstream_orf, strand, gap_bp
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


def load_downstream_orfs(results_dir: Path) -> set:
    f = results_dir / "cotranscription" / "downstream_orf_ids.txt"
    if not f.exists():
        return set()
    with open(f) as fh:
        return {line.strip() for line in fh if line.strip()}


def load_hit_to_contig_mapping(results_dir: Path) -> list:
    # Columns: blast_hit_id, contig_orf_id, genome, strand, genome_start, genome_stop
    f = results_dir / "cotranscription" / "hit_to_contig_mapping.tsv"
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


def load_interproscan(results_dir: Path) -> list:
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


def load_binomial(results_dir: Path) -> list:
    # Columns: V13 (ipr_acc), V14 (ipr_desc), p.value, p_adju, p_scaled, probability, n
    f = results_dir / "BinomialAnalysis.csv"
    if not f.exists():
        return []
    rows = []
    with open(f, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


def load_defense_scores(results_dir: Path) -> list:
    f = results_dir / "defense_scores.tsv"
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


def load_contig_orfs(results_dir: Path) -> list:
    """Parse contig ORF headers from contig_orfs/all_contigs.faa.
    Returns list of {orf_id, contig, start, end, strand}.
    """
    f = results_dir / "contig_orfs" / "all_contigs.faa"
    if not f.exists():
        return []
    orfs = []
    with open(f) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].split(" # ")
            if len(parts) < 4:
                continue
            full_id = parts[0].strip()
            try:
                rel_start = int(parts[1])
                rel_stop = int(parts[2])
                strand = int(parts[3])
            except ValueError:
                continue
            # ORF ID format: genome_start-end:._orfnum
            m = re.match(r"(.+?)_(\d+)-(\d+):\._(\d+)$", full_id)
            if m:
                contig_id = f"{m.group(1)}_{m.group(2)}-{m.group(3)}"
                orfs.append({
                    "orf_id": full_id,
                    "contig": contig_id,
                    "start": rel_start,
                    "end": rel_stop,
                    "strand": strand,
                })
    return orfs


def load_organism_names(results_dir: Path) -> dict:
    """Return genome_id → genus mapping.

    Priority:
    1. genus_manifest.tsv saved alongside results (small, always transferable)
    2. Scan genomes/ directory FASTA headers (only available if genomes were staged locally)
    """
    organism_map = {}

    # 1. Try compact manifest written by stage_genomes rule
    manifest = results_dir / "genus_manifest.tsv"
    if manifest.exists():
        with open(manifest) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("genome_id"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    organism_map[parts[0]] = parts[1]
        return organism_map

    # 2. Fall back to scanning genome FASTA files
    genomes_dir = results_dir / "genomes"
    if not genomes_dir.exists():
        # Try alternate layout where pipeline ran from a nested subdirectory
        sample_name = results_dir.name
        root = results_dir.parent.parent
        try:
            for sub in root.iterdir():
                if not sub.is_dir():
                    continue
                alt = sub / "results" / sample_name / "genomes"
                if alt.exists():
                    genomes_dir = alt
                    break
        except PermissionError:
            pass

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
                    genus = re.sub(r"[^A-Za-z]", "", genus) or "Unknown"
                    organism_map[genome_id] = genus
                else:
                    organism_map[genome_id] = "Unknown"
        except Exception:
            organism_map[genome_id] = "Unknown"
    return organism_map


# ---------------------------------------------------------------------------
# Chart helpers
# ---------------------------------------------------------------------------

def _fig_html(fig) -> str:
    """Serialize a Plotly figure to HTML fragment (no CDN — loaded in <head>)."""
    return fig.to_html(
        include_plotlyjs=False,
        full_html=False,
        config={"displayModeBar": True, "displaylogo": False,
                "modeBarButtonsToRemove": ["select2d", "lasso2d"]},
    )


def _no_data(msg: str = "No data available.") -> str:
    return f"<p class='no-data'>{msg}</p>"


# ---------------------------------------------------------------------------
# Chart builders — Summary
# ---------------------------------------------------------------------------

def chart_pipeline_funnel(metrics: dict, blast: list, defense: list) -> str:
    genomes = metrics.get("total_genomes", 0)
    blast_genomes = len({r["genome"] for r in blast if "genome" in r})
    contigs = metrics.get("contigs_analyzed", 0)
    def_systems = len(defense)

    stages = ["Input Genomes", "Genomes with BLAST Hits",
              "Contigs Analyzed", "Defense Systems"]
    values = [genomes, blast_genomes, contigs, def_systems]

    colors = [TEAL[0], TEAL[1], TEAL[2], TEAL[3]]
    fig = go.Figure(go.Funnel(
        y=stages, x=values,
        marker=dict(color=colors),
        textinfo="value+percent initial",
        hovertemplate="<b>%{y}</b><br>Count: %{x:,}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=30, l=200, r=80)),
        title="Pipeline Analysis Funnel",
        height=300,
    )
    return _fig_html(fig)


def _genus_bar(counts: Counter, title: str) -> str:
    top = counts.most_common(20)
    if not top:
        return _no_data("No genus data available.")
    labels, vals = zip(*top)
    colors = (TEAL * 5)[:len(labels)]
    fig = go.Figure(go.Bar(
        y=list(labels), x=list(vals), orientation="h",
        marker_color=colors,
        hovertemplate="<b>%{y}</b><br>Genomes: %{x}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=40, l=160, r=30)),
        title=title,
        xaxis_title="Number of Genomes",
        yaxis=dict(categoryorder="total ascending"),
        height=max(350, len(labels) * 22),
    )
    return _fig_html(fig)


def chart_organism_breakdown(organism_map: dict) -> str:
    if not organism_map:
        return _no_data("No genome data available.")
    counts = Counter(organism_map.values())
    return _genus_bar(
        counts,
        f"Input Genomes by Genus (top 20 of {len(counts):,} genera)",
    )


def chart_output_organism_breakdown(organism_map: dict, blast: list) -> str:
    """Genus breakdown limited to genomes that had at least one BLAST hit."""
    if not organism_map or not blast:
        return _no_data("No BLAST hit genomes available.")
    hit_genomes = {r["genome"] for r in blast if "genome" in r}
    hit_genera = [organism_map[g] for g in hit_genomes if g in organism_map]
    if not hit_genera:
        return _no_data("Genome IDs from BLAST results not found in genome directory.")
    counts = Counter(hit_genera)
    return _genus_bar(
        counts,
        f"Genomes with BLAST Hits by Genus (top 20 of {len(counts):,} genera)",
    )


# ---------------------------------------------------------------------------
# Chart builders — BLAST
# ---------------------------------------------------------------------------

def chart_blast_identity(blast: list) -> str:
    if not blast:
        return _no_data("No BLAST results available.")
    pidents = []
    for r in blast:
        try:
            pidents.append(float(r["pident"]))
        except (KeyError, ValueError):
            pass
    if not pidents:
        return _no_data("No identity data available.")
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
        height=320,
    )
    return _fig_html(fig)


def chart_blast_scatter(blast: list) -> str:
    if not blast:
        return _no_data("No BLAST results available.")
    pidents, lengths, evalues = [], [], []
    for r in blast:
        try:
            pidents.append(float(r["pident"]))
            lengths.append(int(r["length"]))
            evalues.append(float(r["evalue"]))
        except (KeyError, ValueError):
            pass
    if not pidents:
        return _no_data()
    # log10 evalue for color
    log_eval = [-math.log10(max(e, 1e-300)) for e in evalues]
    fig = go.Figure(go.Scatter(
        x=pidents, y=lengths,
        mode="markers",
        marker=dict(
            color=log_eval,
            colorscale=[[0, TEAL[8]], [1, TEAL[0]]],
            size=5, opacity=0.6,
            colorbar=dict(title="-log₁₀(E-value)", thickness=14),
            showscale=True,
        ),
        hovertemplate=(
            "Identity: %{x:.1f}%<br>"
            "Alignment length: %{y}<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=60, l=70, r=80)),
        title="BLAST Hits: Identity vs. Alignment Length",
        xaxis_title="% Identity",
        yaxis_title="Alignment Length (aa)",
        height=340,
    )
    return _fig_html(fig)


def chart_blast_hits_per_genome(blast: list) -> str:
    if not blast:
        return _no_data("No BLAST results available.")
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
        height=320,
    )
    return _fig_html(fig)


def chart_blast_table(blast: list, max_rows: int = 100) -> str:
    if not blast:
        return _no_data("No BLAST results available.")
    def sort_key(r):
        try:
            return float(r.get("pident", 0))
        except ValueError:
            return 0.0
    top = sorted(blast, key=sort_key, reverse=True)[:max_rows]
    cols = ["query_id", "subject_id", "pident", "length", "evalue", "genome"]
    return html_table(
        headers=["Query", "Subject", "% Identity", "Align. Len.", "E-value", "Genome"],
        rows=[[_esc(r.get(c, "")) for c in cols] for r in top],
        table_id="blast-table",
        caption=f"Top {len(top)} BLAST hits (sorted by % identity, highest first)",
    )


# ---------------------------------------------------------------------------
# Chart builders — Co-transcription
# ---------------------------------------------------------------------------

def chart_cotx_distance(cotx: list) -> str:
    if not cotx:
        return _no_data("No co-transcription results available.")
    try:
        gaps = [int(r["gap_bp"]) for r in cotx if "gap_bp" in r]
    except Exception:
        gaps = []
    if not gaps:
        return _no_data("Gap data unavailable.")
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
    return _fig_html(fig) + note


def chart_cotx_domains(ips: list, cotx: list) -> str:
    """Top InterPro domains found specifically in co-transcribed ORFs."""
    if not ips or not cotx:
        return _no_data("No co-transcription or InterProScan data available.")
    downstream_ids = {r["downstream_orf"] for r in cotx if "downstream_orf" in r}
    if not downstream_ids:
        return _no_data("No downstream ORF IDs found in co-transcription data.")
    cotx_ips = [r for r in ips if r.get("protein") in downstream_ids]
    if not cotx_ips:
        return _no_data("No InterProScan annotations for co-transcribed ORFs.")
    descs = [r["ipr_desc"] for r in cotx_ips
             if r.get("ipr_desc") and r["ipr_desc"] not in ("-", "")]
    if not descs:
        descs = [r["sig_desc"] for r in cotx_ips
                 if r.get("sig_desc") and r["sig_desc"] not in ("-", "")]
    if not descs:
        return _no_data("No domain descriptions in co-transcribed ORF annotations.")
    counts = Counter(descs)
    top = counts.most_common(20)
    labels, vals = zip(*top)
    colors = (TEAL * 4)[:len(labels)]
    fig = go.Figure(go.Bar(
        y=list(labels), x=list(vals), orientation="h",
        marker_color=colors,
        hovertemplate="<b>%{y}</b><br>Hits: %{x}<extra></extra>",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=40, l=330, r=30)),
        title="Top InterPro Domains in Co-transcribed ORFs (top 20)",
        xaxis_title="Number of Hits",
        yaxis=dict(categoryorder="total ascending"),
        height=520,
    )
    note = (f"<p class='chart-note'>{len(downstream_ids):,} downstream ORFs; "
            f"{len(cotx_ips):,} IPS annotations shown.</p>")
    return _fig_html(fig) + note


def chart_defense_scores(scores: list) -> str:
    if not scores:
        return _no_data("No defense scores available.")
    cols = list(scores[0].keys())
    score_col = next((c for c in cols if "score" in c.lower()), None)

    def parse_score(r):
        try:
            return float(r.get(score_col, 0)) if score_col else 0.0
        except ValueError:
            return 0.0

    if score_col is None:
        # Unknown format — just show all columns as-is
        return html_table(
            headers=cols,
            rows=[[_esc(r.get(c, "")) for c in cols] for r in scores[:100]],
            caption="Defense scores (first 100 rows)",
        )

    top50    = sorted(scores, key=parse_score, reverse=True)[:50]
    bottom50 = sorted(scores, key=parse_score)[:50]

    display_cols = cols  # show all columns from the TSV
    headers = [_esc(c) for c in display_cols]

    top_table = html_table(
        headers=headers,
        rows=[[_esc(r.get(c, "")) for c in display_cols] for r in top50],
        table_id="scores-top-table",
        caption=f"Top 50 by {score_col} (highest = most defense-associated)",
    )
    bottom_table = html_table(
        headers=headers,
        rows=[[_esc(r.get(c, "")) for c in display_cols] for r in bottom50],
        table_id="scores-bottom-table",
        caption=f"Bottom 50 by {score_col} (lowest scores)",
    )
    return (f"<h4 style='margin:12px 0 6px;color:var(--teal-dark)'>Top 50</h4>"
            + top_table
            + f"<h4 style='margin:18px 0 6px;color:var(--teal-dark)'>Bottom 50</h4>"
            + bottom_table)


# ---------------------------------------------------------------------------
# Chart builders — DefenseFinder
# ---------------------------------------------------------------------------

def chart_defense_categories(defense: list) -> str:
    if not defense:
        return _no_data("No DefenseFinder results available.")
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
    return _fig_html(fig)


def chart_defense_sunburst(defense: list) -> str:
    if not defense:
        return _no_data("No DefenseFinder results available.")
    ct_counts = Counter((r["category"], r["subtype"]) for r in defense)
    cat_totals = Counter(r["category"] for r in defense)
    ids, parents, labels, values, colors = [], [], [], [], []
    for cat, total in cat_totals.items():
        ids.append(f"cat:{cat}")
        parents.append("")
        labels.append(cat)
        values.append(total)
        colors.append(CATEGORY_COLORS.get(cat, TEAL[4]))
    for (cat, subtype), cnt in ct_counts.items():
        ids.append(f"subtype:{cat}/{subtype}")
        parents.append(f"cat:{cat}")
        labels.append(subtype)
        values.append(cnt)
        colors.append(CATEGORY_COLORS.get(cat, TEAL[4]))
    fig = go.Figure(go.Sunburst(
        ids=ids, parents=parents, labels=labels, values=values,
        marker=dict(colors=colors),
        hovertemplate="<b>%{label}</b><br>Systems: %{value}<extra></extra>",
        insidetextorientation="radial",
        branchvalues="total",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=10, l=10, r=10)),
        title="Category × DefenseFinder Type (click to expand)",
        height=440,
    )
    return _fig_html(fig)


def chart_defense_top_types(defense: list) -> str:
    if not defense:
        return _no_data("No DefenseFinder results available.")
    counts = Counter(r["type"] for r in defense)
    top = counts.most_common(15)
    labels, vals = zip(*top)
    colors = [
        CATEGORY_COLORS.get(
            next((r["category"] for r in defense if r["type"] == lbl), "Unknown mechanism"),
            TEAL[4],
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
    return _fig_html(fig)


def chart_defense_table(defense: list) -> str:
    if not defense:
        return _no_data("No DefenseFinder results available.")
    counts = Counter((r["type"], r["category"]) for r in defense)
    rows = sorted(counts.items(), key=lambda x: -x[1])
    return html_table(
        headers=["Type", "Category", "Count"],
        rows=[[_esc(typ), _esc(cat), str(cnt)] for (typ, cat), cnt in rows],
        table_id="defense-table",
        caption=f"{len(defense)} total systems across {len(rows)} unique types",
    )


# ---------------------------------------------------------------------------
# Chart builders — Domain Annotations
# ---------------------------------------------------------------------------

def chart_top_domains(ips: list) -> str:
    """Top 20 InterPro domains across all co-transcribed ORFs (contig IPS)."""
    if not ips:
        return _no_data("No InterProScan results available.")
    descs = [r["ipr_desc"] for r in ips if r.get("ipr_desc") and r["ipr_desc"] != "-"]
    if not descs:
        descs = [r["sig_desc"] for r in ips if r.get("sig_desc") and r["sig_desc"] != "-"]
    counts = Counter(descs)
    top = counts.most_common(20)
    if not top:
        return _no_data("No domain descriptions found.")
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
    return _fig_html(fig)


def chart_ips_analysis_pie(ips: list) -> str:
    """Breakdown of InterProScan analysis types."""
    if not ips:
        return _no_data("No InterProScan results available.")
    counts = Counter(r["analysis"] for r in ips if r.get("analysis"))
    if not counts:
        return _no_data()
    labels = list(counts.keys())
    vals = list(counts.values())
    colors = (TEAL * 4)[:len(labels)]
    fig = go.Figure(go.Pie(
        labels=labels, values=vals,
        marker=dict(colors=colors),
        hovertemplate="<b>%{label}</b><br>Hits: %{value:,}<br>%{percent}<extra></extra>",
        textinfo="label+percent",
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=30, l=30, r=30)),
        title="InterProScan Analysis Type Distribution",
        height=360,
        showlegend=False,
    )
    return _fig_html(fig)


def chart_binomial_bar(binomial: list) -> str:
    """Bar chart of top 20 significantly enriched domains from binomial analysis."""
    if not binomial:
        return _no_data("No binomial analysis results available.")
    # Filter significant (p_adju < 0.05) and sort
    sig = []
    for r in binomial:
        try:
            p_adj = float(r.get("p_adju", 1.0))
        except (ValueError, TypeError):
            p_adj = 1.0
        if p_adj < 0.05:
            sig.append(r)
    if not sig:
        return _no_data("No significantly enriched domains (p_adju < 0.05) found.")
    sig.sort(key=lambda r: float(r.get("p_adju", 1.0)))
    top = sig[:20]
    labels = [r.get("V14", r.get("interpro_desc", r.get("V13", "?"))) for r in top]
    try:
        p_vals = [-math.log10(max(float(r.get("p_adju", 1.0)), 1e-300)) for r in top]
    except (ValueError, TypeError):
        p_vals = [1.0] * len(top)
    colors = (TEAL * 4)[:len(top)]
    fig = go.Figure(go.Bar(
        y=labels, x=p_vals, orientation="h",
        marker_color=colors,
        customdata=[[r.get("V13", ""), r.get("p_adju", "")] for r in top],
        hovertemplate=(
            "<b>%{y}</b><br>"
            "Accession: %{customdata[0]}<br>"
            "Adj. p-value: %{customdata[1]}<br>"
            "-log₁₀(p_adj): %{x:.2f}<extra></extra>"
        ),
    ))
    fig.update_layout(
        **_layout(margin=dict(t=50, b=40, l=340, r=30)),
        title=f"Top Enriched Domains — Binomial Analysis ({len(sig)} significant)",
        xaxis_title="-log₁₀(adjusted p-value)",
        yaxis=dict(categoryorder="total ascending"),
        height=max(400, len(top) * 26),
    )
    note = (f"<p class='chart-note'>"
            f"{len(sig):,} / {len(binomial):,} domains significant at p_adj &lt; 0.05</p>")
    return _fig_html(fig) + note


def chart_binomial_table(binomial: list, max_rows: int = 200) -> str:
    if not binomial:
        return _no_data("No binomial analysis results available.")
    sig = []
    for r in binomial:
        try:
            p_adj = float(r.get("p_adju", 1.0))
        except (ValueError, TypeError):
            p_adj = 1.0
        if p_adj < 0.05:
            sig.append((p_adj, r))
    sig.sort(key=lambda x: x[0])
    rows_to_show = [r for _, r in sig[:max_rows]]
    cols = ["V13", "V14", "p.value", "p_adju", "n", "probability"]
    avail = [c for c in cols if any(c in r for r in rows_to_show)]
    headers_map = {
        "V13": "IPR Accession", "V14": "Description",
        "p.value": "p-value", "p_adju": "Adj. p-value",
        "n": "Count", "probability": "Probability",
    }
    return html_table(
        headers=[headers_map.get(c, c) for c in avail],
        rows=[[_esc(r.get(c, "")) for c in avail] for r in rows_to_show],
        table_id="binomial-table",
        caption=f"Significantly enriched domains (p_adj < 0.05) — top {len(rows_to_show)}",
    )


# ---------------------------------------------------------------------------
# Locus viewer (optional)
# ---------------------------------------------------------------------------

def _is_mge_domain(sig_acc: str, sig_desc: str) -> bool:
    if sig_acc and sig_acc in MGE_PFAM_ACCESSIONS:
        return True
    if sig_desc:
        dl = sig_desc.lower()
        return any(kw in dl for kw in MGE_KEYWORDS)
    return False


def _gene_arrow_path(x_start, x_end, y_center, strand,
                     arrow_height=0.4, arrow_head_frac=0.15) -> str:
    """Return SVG path string for a boxarrow gene shape."""
    gene_len = abs(x_end - x_start)
    head_len = min(gene_len * arrow_head_frac, gene_len * 0.3)
    hh = arrow_height / 2
    if strand >= 0:
        return (f"M {x_start},{y_center - hh} "
                f"L {x_end - head_len},{y_center - hh} "
                f"L {x_end},{y_center} "
                f"L {x_end - head_len},{y_center + hh} "
                f"L {x_start},{y_center + hh} Z")
    else:
        return (f"M {x_end},{y_center - hh} "
                f"L {x_start + head_len},{y_center - hh} "
                f"L {x_start},{y_center} "
                f"L {x_start + head_len},{y_center + hh} "
                f"L {x_end},{y_center + hh} Z")


def build_locus_figure(
    contig_id: str,
    orfs: list,          # filtered to this contig, sorted by start
    query_hit_orf_ids: set,
    downstream_orf_ids: set,
    ips_by_protein: dict,   # protein_id → list of {sig_acc, sig_desc, ipr_desc}
    defense_gene_ids: set,
) -> go.Figure:
    """Build a Plotly gene-arrow figure for one contig."""
    if not orfs:
        fig = go.Figure()
        fig.add_annotation(text="No ORFs found", xref="paper", yref="paper",
                           x=0.5, y=0.5, showarrow=False)
        return fig

    min_pos = min(o["start"] for o in orfs)
    max_pos = max(o["end"] for o in orfs)
    span = max(max_pos - min_pos, 1)
    x_range = [min_pos - span * 0.03, max_pos + span * 0.03]

    shapes = []
    hover_x, hover_y, hover_text = [], [], []
    y_center = 0.5

    for orf in orfs:
        orf_id = orf["orf_id"]
        start, end, strand = orf["start"], orf["end"], orf["strand"]

        if orf_id in query_hit_orf_ids:
            color = QUERY_HIT_COLOR
            ann = "Query Hit (Recombinase)"
        elif orf_id in downstream_orf_ids:
            color = DOWNSTREAM_COLOR
            domains = ips_by_protein.get(orf_id, [])
            descs = list({d.get("ipr_desc") or d.get("sig_desc", "?")
                          for d in domains if (d.get("ipr_desc") or d.get("sig_desc")) not in ("", "-", None)})
            ann = "<br>".join(descs[:3]) if descs else "Co-transcribed downstream"
        elif orf_id in defense_gene_ids:
            color = DEFENSE_COLOR
            ann = "Defense system gene"
        elif any(_is_mge_domain(d.get("sig_acc", ""), d.get("sig_desc", ""))
                 for d in ips_by_protein.get(orf_id, [])):
            color = MGE_COLOR
            mge = next(d for d in ips_by_protein[orf_id]
                       if _is_mge_domain(d.get("sig_acc", ""), d.get("sig_desc", "")))
            ann = f"MGE: {mge.get('sig_desc', 'mobile element')}"
        elif orf_id in ips_by_protein:
            color = IPS_OTHER_COLOR
            d0 = ips_by_protein[orf_id][0]
            ann = d0.get("ipr_desc") or d0.get("sig_desc") or "Annotated domain"
        else:
            color = NO_ANNOTATION_COLOR
            ann = "No annotation"

        path = _gene_arrow_path(start, end, y_center, strand)
        shapes.append(dict(type="path", path=path, fillcolor=color,
                           line=dict(color="#555", width=0.8), layer="above"))

        cx = (start + end) / 2
        hover_x.append(cx)
        hover_y.append(y_center)
        hover_text.append(
            f"<b>{orf_id}</b><br>"
            f"Position: {start:,}–{end:,}<br>"
            f"Strand: {'+' if strand >= 0 else '−'}<br>"
            f"Length: {abs(end - start):,} bp<br>"
            f"{ann}"
        )

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=hover_x, y=hover_y, mode="markers",
        marker=dict(size=18, opacity=0),
        hoverinfo="text", hovertext=hover_text, showlegend=False,
    ))

    # Scale bar
    scale_len = span / 5
    scale_txt = f"{scale_len/1000:.1f} kb" if scale_len >= 1000 else f"{scale_len:.0f} bp"
    shapes.append(dict(type="line",
                       x0=min_pos, x1=min_pos + scale_len,
                       y0=0.12, y1=0.12,
                       line=dict(color="#333", width=2)))
    fig.add_annotation(x=min_pos + scale_len / 2, y=0.06,
                       text=scale_txt, showarrow=False, font=dict(size=10))

    fig.update_layout(
        shapes=shapes,
        title=f"Locus: {contig_id}",
        xaxis=dict(title="Position (bp)", range=x_range,
                   showgrid=True, gridcolor="#e0e8e8"),
        yaxis=dict(range=[0, 1], showticklabels=False, showgrid=False),
        height=260,
        hovermode="closest",
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=False,
        margin=dict(t=40, b=50, l=50, r=20),
        font=dict(family="system-ui, sans-serif", size=11, color="#2d4a4a"),
    )
    return fig


def build_locus_data(
    contig_orfs: list,
    hit_mapping: list,
    downstream_orf_ids: set,
    ips: list,
    defense_genes: list,
    max_contigs: int = 250,
) -> dict:
    """
    Pre-render Plotly figures for each contig and return as {contig_id: json_str}.
    Limited to max_contigs to keep HTML size manageable.
    """
    # Index ORFs by contig
    orfs_by_contig = defaultdict(list)
    for o in contig_orfs:
        orfs_by_contig[o["contig"]].append(o)
    for lst in orfs_by_contig.values():
        lst.sort(key=lambda o: o["start"])

    # Set of query-hit ORF IDs
    query_hit_orf_ids = {r["contig_orf_id"] for r in hit_mapping if "contig_orf_id" in r}

    # IPS lookup: protein_id → list of annotation dicts
    ips_by_protein = defaultdict(list)
    for r in ips:
        pid = r.get("protein", "")
        if pid:
            ips_by_protein[pid].append({
                "sig_acc": r.get("sig_acc", ""),
                "sig_desc": r.get("sig_desc", ""),
                "ipr_desc": r.get("ipr_desc", ""),
            })

    # Defense gene IDs (hit_id column in genes TSV matches contig ORF IDs)
    defense_gene_id_set = set()
    for g in defense_genes:
        hid = g.get("hit_id", "")
        if hid:
            defense_gene_id_set.add(hid)

    # Prioritize contigs with query hits, then sort alphabetically
    contig_ids = sorted(orfs_by_contig.keys())
    with_hits = [c for c in contig_ids if any(
        o["orf_id"] in query_hit_orf_ids for o in orfs_by_contig[c]
    )]
    without_hits = [c for c in contig_ids if c not in set(with_hits)]
    ordered = with_hits + without_hits
    ordered = ordered[:max_contigs]

    locus_data = {}
    for contig_id in ordered:
        fig = build_locus_figure(
            contig_id=contig_id,
            orfs=orfs_by_contig[contig_id],
            query_hit_orf_ids=query_hit_orf_ids,
            downstream_orf_ids=downstream_orf_ids,
            ips_by_protein=ips_by_protein,
            defense_gene_ids=defense_gene_id_set,
        )
        locus_data[contig_id] = fig.to_json()
    return locus_data


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

def _esc(s: str) -> str:
    """HTML-escape a string."""
    return (str(s)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;"))


def metric_card(label: str, value: str, subtitle: str = "") -> str:
    sub = f"<div class='card-sub'>{subtitle}</div>" if subtitle else ""
    return f"<div class='metric-card'><div class='card-label'>{label}</div><div class='card-value'>{value}</div>{sub}</div>"


def chart_card(title: str, html: str) -> str:
    return f"<div class='chart-card'><h3 class='chart-title'>{title}</h3><div class='chart-inner'>{html}</div></div>"


def two_col(left: str, right: str) -> str:
    return f"<div class='two-col'>{left}{right}</div>"


def three_col(a: str, b: str, c: str) -> str:
    return f"<div class='three-col'>{a}{b}{c}</div>"


def html_table(headers: list, rows: list, table_id: str = "",
               caption: str = "", max_rows: int = 500) -> str:
    if not rows:
        return _no_data("No data to display.")
    id_attr = f" id='{table_id}'" if table_id else ""
    cap = f"<caption>{_esc(caption)}</caption>" if caption else ""
    th = "".join(f"<th>{h}</th>" for h in headers)
    displayed = rows[:max_rows]
    trs = "".join(
        "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"
        for row in displayed
    )
    note = ""
    if len(rows) > max_rows:
        note = f"<p class='chart-note'>Showing first {max_rows:,} of {len(rows):,} rows.</p>"
    return (f"<div class='table-wrap'><table{id_attr} class='data-table'>"
            f"{cap}<thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table></div>{note}")


# ---------------------------------------------------------------------------
# CSS + JS
# ---------------------------------------------------------------------------

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
header {
    background: var(--teal-dark);
    color: white;
    position: sticky;
    top: 0;
    z-index: 200;
    box-shadow: 0 2px 8px rgba(0,0,0,0.2);
}
.header-inner {
    max-width: 1400px; margin: 0 auto; padding: 14px 24px;
    display: flex; align-items: center; justify-content: space-between;
}
.header-title { font-size: 1.35rem; font-weight: 700; letter-spacing: 0.5px; }
.header-meta  { font-size: 0.8rem; color: var(--teal-pale); text-align: right; }

/* ---------- Tab nav ---------- */
.tab-bar {
    background: var(--teal-mid);
    overflow-x: auto; white-space: nowrap;
    position: sticky; top: 53px; z-index: 100;
}
.tab-bar-inner { max-width: 1400px; margin: 0 auto; padding: 0 24px; }
.tab-btn {
    display: inline-block; background: none; border: none; cursor: pointer;
    color: #d0e8e8; padding: 10px 16px;
    font-size: 0.82rem; font-weight: 500;
    font-family: inherit;
    transition: background 0.15s, color 0.15s;
    border-bottom: 3px solid transparent;
}
.tab-btn:hover { background: rgba(255,255,255,0.12); color: white; }
.tab-btn.active { color: white; border-bottom-color: var(--teal-pale); }

/* ---------- Tab panels ---------- */
.tab-panel { display: none; }
.tab-panel.active { display: block; }

/* ---------- Main layout ---------- */
main { max-width: 1400px; margin: 0 auto; padding: 28px 24px 64px; }

/* ---------- Metric cards ---------- */
.metrics-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(165px, 1fr));
    gap: 14px; margin-bottom: 32px;
}
.metric-card {
    background: var(--teal-dark); color: white;
    border-radius: 10px; padding: 16px 14px; text-align: center;
}
.card-label { font-size: 0.7rem; color: var(--teal-pale); text-transform: uppercase; letter-spacing: 0.6px; margin-bottom: 6px; }
.card-value { font-size: 1.75rem; font-weight: 700; line-height: 1; }
.card-sub   { font-size: 0.7rem; color: var(--teal-light); margin-top: 4px; }

/* ---------- Section headings ---------- */
.section-title {
    font-size: 1rem; font-weight: 700; color: var(--teal-dark);
    border-left: 4px solid var(--teal-light);
    padding-left: 10px; margin: 24px 0 14px;
}

/* ---------- Chart cards ---------- */
.chart-card {
    background: var(--white); border: 1px solid var(--border);
    border-radius: 10px; padding: 14px; margin-bottom: 16px;
    box-shadow: 0 1px 4px rgba(45,74,74,0.06);
}
.chart-title {
    font-size: 0.82rem; font-weight: 600; color: var(--teal-mid);
    text-transform: uppercase; letter-spacing: 0.4px; margin-bottom: 10px;
}
.chart-inner { width: 100%; overflow-x: auto; }
.two-col {
    display: grid; grid-template-columns: 1fr 1fr; gap: 16px;
}
.three-col {
    display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px;
}
@media (max-width: 1000px) { .three-col { grid-template-columns: 1fr 1fr; } }
@media (max-width: 700px)  { .two-col, .three-col { grid-template-columns: 1fr; } }

/* ---------- Tables ---------- */
.table-wrap { overflow-x: auto; margin-bottom: 6px; }
.data-table {
    width: 100%; border-collapse: collapse; font-size: 0.8rem;
}
.data-table caption {
    text-align: left; font-size: 0.75rem; color: #6b8888;
    padding: 0 0 6px; font-style: italic;
}
.data-table th {
    background: var(--teal-dark); color: white;
    padding: 7px 10px; text-align: left; font-size: 0.75rem;
    white-space: nowrap; position: sticky; top: 0;
}
.data-table td {
    padding: 5px 10px; border-bottom: 1px solid var(--border);
    max-width: 320px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
}
.data-table tr:hover td { background: #eef5f5; }

/* ---------- Locus viewer ---------- */
.locus-controls {
    display: flex; align-items: center; gap: 14px;
    margin-bottom: 14px; flex-wrap: wrap;
}
.locus-controls label { font-size: 0.82rem; font-weight: 600; color: var(--teal-dark); }
#contig-select {
    padding: 6px 10px; border: 1px solid var(--border);
    border-radius: 6px; font-size: 0.82rem; min-width: 260px;
    background: white;
}
.locus-nav { display: flex; gap: 8px; }
.locus-nav button {
    padding: 5px 14px; background: var(--teal-mid); color: white;
    border: none; border-radius: 6px; cursor: pointer; font-size: 0.8rem;
}
.locus-nav button:hover { background: var(--teal-dark); }
.locus-legend {
    display: flex; flex-wrap: wrap; gap: 12px;
    padding: 10px 14px; background: var(--white);
    border: 1px solid var(--border); border-radius: 8px; margin-bottom: 14px;
    font-size: 0.78rem;
}
.legend-item { display: flex; align-items: center; gap: 5px; }
.legend-swatch {
    width: 20px; height: 12px; border-radius: 2px;
    border: 1px solid rgba(0,0,0,0.15); flex-shrink: 0;
}
#locus-info { font-size: 0.8rem; color: #6b8888; margin-bottom: 8px; min-height: 18px; }

/* ---------- Misc ---------- */
.chart-note { font-size: 0.75rem; color: #6b8888; margin-top: 5px; padding-left: 4px; }
.no-data { color: #7a9090; font-style: italic; padding: 24px; text-align: center; }
footer {
    text-align: center; padding: 24px; color: #7a9090;
    font-size: 0.78rem; border-top: 1px solid var(--border);
}
"""

TAB_JS = """
(function() {
    function showTab(id) {
        document.querySelectorAll('.tab-panel').forEach(function(p) {
            p.classList.remove('active');
        });
        document.querySelectorAll('.tab-btn').forEach(function(b) {
            b.classList.remove('active');
        });
        var panel = document.getElementById('tab-' + id);
        var btn   = document.getElementById('btn-' + id);
        if (panel) panel.classList.add('active');
        if (btn)   btn.classList.add('active');
        // Trigger Plotly relayout so hidden charts resize properly
        if (typeof Plotly !== 'undefined') {
            var plots = panel ? panel.querySelectorAll('.plotly-graph-div') : [];
            plots.forEach(function(el) { Plotly.relayout(el, {}); });
        }
    }
    window.showTab = showTab;
    document.addEventListener('DOMContentLoaded', function() {
        var first = document.querySelector('.tab-btn');
        if (first) showTab(first.getAttribute('data-tab'));
    });
})();
"""

LOCUS_JS_TEMPLATE = """
(function() {
    var LOCUS_DATA = LOCUS_DATA_PLACEHOLDER;
    var contigIds  = CONTIG_IDS_PLACEHOLDER;
    var currentIdx = 0;

    function renderLocus(contigId) {
        var info = document.getElementById('locus-info');
        if (info) info.textContent = 'Contig: ' + contigId;
        var container = document.getElementById('locus-plot');
        if (!container) return;
        var figJson = LOCUS_DATA[contigId];
        if (!figJson) {
            Plotly.purge(container);
            container.innerHTML = '<p style="padding:20px;color:#7a9090">No data for this contig.</p>';
            return;
        }
        var fig = JSON.parse(figJson);
        Plotly.react(container, fig.data, fig.layout, {displayModeBar: true, displaylogo: false});
    }

    document.addEventListener('DOMContentLoaded', function() {
        var sel = document.getElementById('contig-select');
        if (sel) {
            sel.addEventListener('change', function() {
                currentIdx = contigIds.indexOf(this.value);
                renderLocus(this.value);
            });
            if (contigIds.length > 0) renderLocus(contigIds[0]);
        }
        var prevBtn = document.getElementById('locus-prev');
        var nextBtn = document.getElementById('locus-next');
        if (prevBtn) prevBtn.addEventListener('click', function() {
            if (currentIdx > 0) {
                currentIdx--;
                sel.value = contigIds[currentIdx];
                renderLocus(contigIds[currentIdx]);
            }
        });
        if (nextBtn) nextBtn.addEventListener('click', function() {
            if (currentIdx < contigIds.length - 1) {
                currentIdx++;
                sel.value = contigIds[currentIdx];
                renderLocus(contigIds[currentIdx]);
            }
        });
    });
})();
"""


# ---------------------------------------------------------------------------
# Full HTML assembly
# ---------------------------------------------------------------------------

def build_html(
    sample_name: str,
    generated_at: str,
    metrics: dict,
    organism_map: dict,
    blast: list,
    defense: list,
    cotx: list,
    ips: list,
    binomial: list,
    scores: list,
    locus_data=None,   # dict or None — None = locus viewer omitted
) -> str:

    def fmt(val, decimals=0):
        if val is None:
            return "—"
        if isinstance(val, float):
            return f"{val:,.{decimals}f}"
        return f"{val:,}" if isinstance(val, int) else str(val)

    # ---- Metric cards ----
    cards = "".join([
        metric_card("Total Genomes",     fmt(metrics.get("total_genomes")),             "input"),
        metric_card("Total Input",       fmt(metrics.get("total_input_mb"), 1) + " Mb", "sequenced"),
        metric_card("BLAST Hits",        fmt(metrics.get("blast_hits")),                "query matches"),
        metric_card("Hits per Mb",       fmt(metrics.get("hits_per_mb"), 3),            ""),
        metric_card("Defense Genes",      fmt(len(defense)),                             "DefenseFinder"),
        metric_card("Contigs Analyzed",  fmt(metrics.get("contigs_analyzed")),          "flanking regions"),
        metric_card("Co-tx Pairs",       fmt(len(cotx)),                                ""),
        metric_card("IPS Annotations",   fmt(len(ips)),                                 "domain hits"),
    ])

    query_name = metrics.get("query_name", "Unknown")

    # ============================================================
    # Tab 1 — Summary
    # ============================================================
    t_summary = f"""
    <div class='metrics-grid'>{cards}</div>
    {chart_card("Pipeline Analysis Funnel",
        chart_pipeline_funnel(metrics, blast, defense))}
    <div class='two-col'>
        {chart_card("Input Genomes by Genus", chart_organism_breakdown(organism_map))}
        {chart_card("Output Genomes by Genus", chart_output_organism_breakdown(organism_map, blast))}
    </div>
    {chart_card("Defense System Categories", chart_defense_categories(defense))}
    """

    # ============================================================
    # Tab 2 — BLAST Results
    # ============================================================
    t_blast = f"""
    <div class='two-col'>
        {chart_card("Identity Distribution", chart_blast_identity(blast))}
        {chart_card("Hits per Genome", chart_blast_hits_per_genome(blast))}
    </div>
    {chart_card("Identity vs. Alignment Length", chart_blast_scatter(blast))}
    <h3 class='section-title'>Top BLAST Hits</h3>
    <div class='chart-card'>
        <div class='chart-inner'>{chart_blast_table(blast)}</div>
    </div>
    """

    # ============================================================
    # Tab 3 — Co-transcription
    # ============================================================
    t_cotx = f"""
    <div class='two-col'>
        {chart_card("Gap Distance Distribution", chart_cotx_distance(cotx))}
        {chart_card("Domains in Co-transcribed ORFs", chart_cotx_domains(ips, cotx))}
    </div>
    <h3 class='section-title'>Defense Scores</h3>
    <div class='chart-card'>
        <div class='chart-inner'>{chart_defense_scores(scores)}</div>
    </div>
    """

    # ============================================================
    # Tab 4 — DefenseFinder
    # ============================================================
    t_defense = f"""
    {chart_card("Defense Categories", chart_defense_categories(defense))}
    <div class='two-col'>
        {chart_card("Category × System Type (Sunburst)", chart_defense_sunburst(defense))}
        {chart_card("Top System Types", chart_defense_top_types(defense))}
    </div>
    <h3 class='section-title'>All Defense Systems</h3>
    <div class='chart-card'>
        <div class='chart-inner'>{chart_defense_table(defense)}</div>
    </div>
    """

    # ============================================================
    # Tab 5 — Domain Annotations
    # ============================================================
    t_domains = f"""
    <div class='two-col'>
        {chart_card("Top 20 InterPro Domains (Co-transcribed ORFs)", chart_top_domains(ips))}
        {chart_card("InterProScan Analysis Types", chart_ips_analysis_pie(ips))}
    </div>
    <h3 class='section-title'>Binomial Domain Enrichment Analysis</h3>
    {chart_card("Enriched Domains (p<sub>adj</sub> &lt; 0.05)", chart_binomial_bar(binomial))}
    <div class='chart-card'>
        <div class='chart-inner'>{chart_binomial_table(binomial)}</div>
    </div>
    """

    # ============================================================
    # Tab 6 — Locus Viewer (optional)
    # ============================================================
    if locus_data is not None and locus_data:
        contig_ids = sorted(locus_data.keys())
        options = "".join(
            f"<option value='{_esc(c)}'>{_esc(c)}</option>" for c in contig_ids
        )
        locus_legend = f"""
        <div class='locus-legend'>
            <div class='legend-item'><div class='legend-swatch' style='background:{QUERY_HIT_COLOR}'></div>Query Hit</div>
            <div class='legend-item'><div class='legend-swatch' style='background:{DOWNSTREAM_COLOR}'></div>Co-transcribed</div>
            <div class='legend-item'><div class='legend-swatch' style='background:{DEFENSE_COLOR}'></div>Defense Gene</div>
            <div class='legend-item'><div class='legend-swatch' style='background:{MGE_COLOR}'></div>Mobile Genetic Element</div>
            <div class='legend-item'><div class='legend-swatch' style='background:{IPS_OTHER_COLOR}'></div>Other Domain</div>
            <div class='legend-item'><div class='legend-swatch' style='background:{NO_ANNOTATION_COLOR};border:1px solid #ccc'></div>No Annotation</div>
        </div>
        """
        t_locus = f"""
        <div class='locus-controls'>
            <label for='contig-select'>Select Contig:</label>
            <select id='contig-select'>{options}</select>
            <div class='locus-nav'>
                <button id='locus-prev'>&larr; Prev</button>
                <button id='locus-next'>Next &rarr;</button>
            </div>
            <span style='font-size:0.78rem;color:#6b8888'>{len(contig_ids):,} contigs</span>
        </div>
        {locus_legend}
        <div id='locus-info'></div>
        <div class='chart-card' style='padding:8px'>
            <div id='locus-plot'></div>
        </div>
        <p class='chart-note'>Gene arrows are colored by annotation priority:
            query hit &gt; co-transcribed &gt; defense gene &gt; MGE &gt; other domain &gt; unannotated.</p>
        """
        locus_js = (LOCUS_JS_TEMPLATE
                    .replace("LOCUS_DATA_PLACEHOLDER", json.dumps(locus_data))
                    .replace("CONTIG_IDS_PLACEHOLDER", json.dumps(contig_ids)))
    else:
        t_locus = None
        locus_js = ""

    # ============================================================
    # Build tab bar and panels
    # ============================================================
    tabs = [
        ("summary",    "Summary"),
        ("blast",      "BLAST Results"),
        ("cotx",       "Co-transcription"),
        ("defense",    "DefenseFinder"),
        ("domains",    "Domain Annotations"),
    ]
    if t_locus is not None:
        tabs.append(("locus", "Locus Viewer"))

    tab_contents = {
        "summary": t_summary,
        "blast":   t_blast,
        "cotx":    t_cotx,
        "defense": t_defense,
        "domains": t_domains,
    }
    if t_locus is not None:
        tab_contents["locus"] = t_locus

    tab_btns = "".join(
        f"<button class='tab-btn' id='btn-{tid}' data-tab='{tid}' onclick='showTab(\"{tid}\")'>{label}</button>"
        for tid, label in tabs
    )
    tab_panels = "".join(
        f"<div class='tab-panel' id='tab-{tid}'>{tab_contents[tid]}</div>"
        for tid, _ in tabs
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>RADS Report — {_esc(sample_name)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>{CSS}</style>
</head>
<body>
<header>
  <div class='header-inner'>
    <div>
      <div class='header-title'>RADS — Recombinase Associated Defense Search</div>
      <div style='font-size:0.8rem;color:var(--teal-pale);margin-top:2px;'>
        Sample: <b>{_esc(sample_name)}</b> &nbsp;|&nbsp; Query: <b>{_esc(query_name)}</b>
      </div>
    </div>
    <div class='header-meta'>
      Generated {_esc(generated_at)}<br>
      defense-finder-models v2.0+
    </div>
  </div>
</header>
<div class='tab-bar'><div class='tab-bar-inner'>{tab_btns}</div></div>
<main>
{tab_panels}
</main>
<footer>
  RADS Pipeline Report &nbsp;·&nbsp; Generated {_esc(generated_at)} &nbsp;·&nbsp;
  <a href='https://github.com/Seandersen/RADS' style='color:#8fb3b3;'>github.com/Seandersen/RADS</a>
</footer>
<script>{TAB_JS}</script>
{f'<script>{locus_js}</script>' if locus_js else ''}
</body>
</html>"""


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate RADS HTML report")
    parser.add_argument("--results", required=True, help="Path to sample results directory")
    parser.add_argument("--output",  required=True, help="Output HTML file path")
    parser.add_argument("--include-locus-viewer", action="store_true",
                        help="Include gene-arrow locus viewer (larger HTML, needs more memory)")
    parser.add_argument("--locus-max-contigs", type=int, default=250,
                        help="Max contigs to pre-render in locus viewer (default: 250)")
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
    binomial     = load_binomial(results_dir)
    scores       = load_defense_scores(results_dir)

    sample_name  = metrics.get("sample_name", results_dir.name)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"  Genomes:         {len(organism_map):,}")
    print(f"  BLAST hits:      {len(blast):,}")
    print(f"  Defense systems: {len(defense):,}")
    print(f"  Co-tx pairs:     {len(cotx):,}")
    print(f"  IPS annotations: {len(ips):,}")
    print(f"  Binomial rows:   {len(binomial):,}")
    print(f"  Defense scores:  {len(scores):,}")

    locus_data = None
    if args.include_locus_viewer:
        print("Building locus viewer data...")
        contig_orfs  = load_contig_orfs(results_dir)
        hit_mapping  = load_hit_to_contig_mapping(results_dir)
        downstream   = load_downstream_orfs(results_dir)
        def_genes    = load_defensefinder_genes(results_dir)
        print(f"  Contig ORFs:     {len(contig_orfs):,}")
        locus_data = build_locus_data(
            contig_orfs=contig_orfs,
            hit_mapping=hit_mapping,
            downstream_orf_ids=downstream,
            ips=ips,
            defense_genes=def_genes,
            max_contigs=args.locus_max_contigs,
        )
        print(f"  Locus contigs:   {len(locus_data):,}")

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
        binomial=binomial,
        scores=scores,
        locus_data=locus_data,
    )

    out_path.write_text(html, encoding="utf-8")
    size_kb = out_path.stat().st_size / 1024
    print(f"Report written → {out_path}  ({size_kb:,.0f} KB)")


if __name__ == "__main__":
    main()
