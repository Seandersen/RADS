"""
Generate a performance estimate figure for RADS pipeline.
Estimates are rough and based on typical HPC workloads; actual times
depend heavily on genome size, BLAST hit rate, and node speed.

Usage:
    python docs/generate_performance_estimate.py
Output:
    docs/performance_estimate.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Timing model (all times in seconds unless noted)
# ---------------------------------------------------------------------------
# Per-genome parallelisable work  (prodigal + diamond makedb + blast + extract steps)
T_PER_GENOME = 120          # seconds per genome

# Global sequential steps (scale weakly with genomes)
T_PRODIGAL_CONTIGS_PER_G = 0.05   # prodigal on extracted contigs
T_DEFENSEFINDER_BASE = 60
T_DEFENSEFINDER_PER_G = 0.3

# InterProScan on contig ORFs only (~100 proteins per genome on average)
# IPS scales roughly linearly with proteins; benefit from threads is limited
# Pfam only (faster):   ~1 min / 1000 proteins at 20 cores
# 4 apps (Pfam+SMART+CDD+SUPERFAMILY): ~5 min / 1000 proteins at 20 cores
CONTIGS_PROTEINS_PER_G = 100       # rough average ORFs in extracted flanking regions
IPS_PFAM_RATE   = 60  / 1000      # sec per protein, 20-core equivalent
IPS_4APPS_RATE  = 300 / 1000      # sec per protein, 20-core equivalent

# Binomial whole-genome InterProScan
# Assumes ~30% of genomes have BLAST hits; ~4000 proteins per genome
HIT_RATE          = 0.30
PROTEINS_PER_GENOME = 4000
BINOMIAL_PFAM_RATE  = IPS_PFAM_RATE
BINOMIAL_4APPS_RATE = IPS_4APPS_RATE

# IPS thread-scaling: diminishing returns above ~8 cores
# We model it as sqrt-ish: effective_rate = rate / min(cores/8, 2.5)
def ips_scale(cores):
    return min(cores / 8, 2.5)

def estimate_hours(n_genomes, cores, ips_4apps=False, binomial=False):
    """Return total estimated wall-clock hours for the pipeline."""
    # --- Parallelisable per-genome phase ---
    t_parallel = (T_PER_GENOME * n_genomes) / cores

    # --- Sequential global phases ---
    t_prodigal_contigs = T_PRODIGAL_CONTIGS_PER_G * n_genomes
    t_df = T_DEFENSEFINDER_BASE + T_DEFENSEFINDER_PER_G * n_genomes

    # InterProScan on contigs
    n_contig_proteins = CONTIGS_PROTEINS_PER_G * n_genomes
    ips_rate = IPS_4APPS_RATE if ips_4apps else IPS_PFAM_RATE
    t_ips = (n_contig_proteins * ips_rate) / ips_scale(cores)

    # Binomial whole-genome InterProScan
    t_binomial = 0
    if binomial:
        n_wg_proteins = HIT_RATE * n_genomes * PROTEINS_PER_GENOME
        t_binomial = (n_wg_proteins * ips_rate) / ips_scale(cores)

    total_seconds = t_parallel + t_prodigal_contigs + t_df + t_ips + t_binomial
    return total_seconds / 3600


def phase_breakdown(n_genomes, cores, ips_4apps=False, binomial=False):
    """Return dict of phase → hours for a stacked bar."""
    t_parallel = (T_PER_GENOME * n_genomes) / cores / 3600

    t_df = (T_DEFENSEFINDER_BASE + T_DEFENSEFINDER_PER_G * n_genomes) / 3600

    n_cp = CONTIGS_PROTEINS_PER_G * n_genomes
    ips_rate = IPS_4APPS_RATE if ips_4apps else IPS_PFAM_RATE
    t_ips = (n_cp * ips_rate) / ips_scale(cores) / 3600

    t_binomial = 0
    if binomial:
        n_wg = HIT_RATE * n_genomes * PROTEINS_PER_GENOME
        t_binomial = (n_wg * ips_rate) / ips_scale(cores) / 3600

    return {
        "Per-genome\n(translate, BLAST, extract)": t_parallel,
        "DefenseFinder": t_df,
        "InterProScan\n(contigs only)": t_ips,
        "Binomial\n(whole-genome IPS)": t_binomial,
    }


# ---------------------------------------------------------------------------
# Plot setup
# ---------------------------------------------------------------------------
genome_range = np.array([50, 100, 250, 500, 1000, 2500, 5000, 10000])
core_counts  = [4, 8, 16, 32, 64]
core_colors  = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4", "#9467bd"]

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.patch.set_facecolor("#f8f9fa")
for ax in axes:
    ax.set_facecolor("#f8f9fa")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# ---------------------------------------------------------------------------
# Panel 1 — Scaling with core count (Standard: Pfam only, no binomial)
# ---------------------------------------------------------------------------
ax1 = axes[0]
ax1.set_title("Core Scaling\n(Pfam-only IPS, no binomial)", fontsize=12, fontweight="bold")

for cores, color in zip(core_counts, core_colors):
    times = [estimate_hours(g, cores, ips_4apps=False, binomial=False)
             for g in genome_range]
    ax1.plot(genome_range, times, "o-", color=color, linewidth=2,
             markersize=5, label=f"{cores} cores")

# Shade the 20-core reference used in the paper
times_20 = [estimate_hours(g, 20, ips_4apps=False, binomial=False)
            for g in genome_range]
ax1.plot(genome_range, times_20, "s--", color="#555555", linewidth=1.5,
         markersize=4, label="20 cores (config default)", zorder=5)

ax1.set_xscale("log")
ax1.set_xlabel("Number of Genomes", fontsize=11)
ax1.set_ylabel("Estimated Wall-Clock Time (hours)", fontsize=11)
ax1.legend(fontsize=9, loc="upper left")
ax1.grid(True, alpha=0.3, linestyle="--")
ax1.set_xticks(genome_range)
ax1.set_xticklabels([str(g) for g in genome_range], rotation=45, ha="right", fontsize=8)

# ---------------------------------------------------------------------------
# Panel 2 — Configuration comparison (fixed 20 cores)
# ---------------------------------------------------------------------------
ax2 = axes[1]
ax2.set_title("Configuration Comparison\n(20 cores)", fontsize=12, fontweight="bold")

configs = [
    ("Minimal\n(no IPS, no binomial)",     False, False, "#aec7e8", "-"),
    ("Standard\n(Pfam only)",              False, False, "#1f77b4", "-"),
    ("Standard + Binomial\n(Pfam only)",   False, True,  "#1f77b4", "--"),
    ("Full IPS\n(4 apps, no binomial)",    True,  False, "#d62728", "-"),
    ("Full IPS + Binomial\n(4 apps)",      True,  True,  "#d62728", "--"),
]

for label, ips4, binom, color, ls in configs:
    if label.startswith("Minimal"):
        # Minimal: no IPS at all — only per-genome + DefenseFinder
        times = [(T_PER_GENOME * g / 20 +
                  T_DEFENSEFINDER_BASE + T_DEFENSEFINDER_PER_G * g) / 3600
                 for g in genome_range]
    else:
        times = [estimate_hours(g, 20, ips_4apps=ips4, binomial=binom)
                 for g in genome_range]
    ax2.plot(genome_range, times, linestyle=ls, color=color,
             linewidth=2, marker="o", markersize=4, label=label)

ax2.set_xscale("log")
ax2.set_xlabel("Number of Genomes", fontsize=11)
ax2.set_ylabel("Estimated Wall-Clock Time (hours)", fontsize=11)
ax2.legend(fontsize=8, loc="upper left")
ax2.grid(True, alpha=0.3, linestyle="--")
ax2.set_xticks(genome_range)
ax2.set_xticklabels([str(g) for g in genome_range], rotation=45, ha="right", fontsize=8)

# ---------------------------------------------------------------------------
# Panel 3 — Phase breakdown stacked bar at 20 cores, 5000 genomes
# ---------------------------------------------------------------------------
ax3 = axes[2]
ax3.set_title("Phase Breakdown\n(5,000 genomes, 20 cores)", fontsize=12, fontweight="bold")

bar_configs = [
    ("Minimal",                False, False),
    ("Standard\n(Pfam)",       False, False),
    ("Std +\nBinomial",        False, True),
    ("Full IPS\n(4 apps)",     True,  False),
    ("Full IPS\n+ Binomial",   True,  True),
]

phase_colors = {
    "Per-genome\n(translate, BLAST, extract)": "#4e79a7",
    "DefenseFinder":                            "#76b7b2",
    "InterProScan\n(contigs only)":             "#f28e2b",
    "Binomial\n(whole-genome IPS)":             "#e15759",
}

n_bars = len(bar_configs)
x = np.arange(n_bars)
bar_width = 0.6

bottoms = np.zeros(n_bars)
phase_names = list(phase_colors.keys())

for phase in phase_names:
    values = []
    for label, ips4, binom in bar_configs:
        if label == "Minimal":
            bd = phase_breakdown(5000, 20, ips_4apps=False, binomial=False)
            # In minimal, IPS and binomial are zero
            bd["InterProScan\n(contigs only)"] = 0
            bd["Binomial\n(whole-genome IPS)"] = 0
        else:
            bd = phase_breakdown(5000, 20, ips_4apps=ips4, binomial=binom)
        values.append(bd[phase])
    ax3.bar(x, values, bar_width, bottom=bottoms,
            color=phase_colors[phase], label=phase)
    bottoms += np.array(values)

ax3.set_xticks(x)
ax3.set_xticklabels([c[0] for c in bar_configs], fontsize=9)
ax3.set_ylabel("Estimated Wall-Clock Time (hours)", fontsize=11)
ax3.legend(fontsize=8, loc="upper left", bbox_to_anchor=(0, 1))
ax3.grid(True, alpha=0.3, linestyle="--", axis="y")

# ---------------------------------------------------------------------------
# Footer note
# ---------------------------------------------------------------------------
fig.text(
    0.5, -0.04,
    "Estimates assume ~120 s/genome for per-genome steps, ~100 contig ORFs/genome, "
    "30% BLAST hit rate, 4,000 proteins/genome for binomial IPS.\n"
    "InterProScan thread scaling is sub-linear; actual times vary with cluster speed, "
    "genome size, and hit rate. Use as order-of-magnitude guide only.",
    ha="center", fontsize=8, color="#555555", style="italic"
)

fig.suptitle("RADS Pipeline — Estimated Runtime by Configuration",
             fontsize=14, fontweight="bold", y=1.02)

plt.tight_layout()
out_path = "docs/performance_estimate.png"
plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
print(f"Saved → {out_path}")
