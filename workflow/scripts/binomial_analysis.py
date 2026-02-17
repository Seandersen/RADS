#!/usr/bin/env python3
"""
Binomial domain enrichment analysis for RADS pipeline.

Compares Pfam domain frequencies in extracted contigs vs. whole genomes
to identify domains that are significantly enriched near query sequences.

Ported from RADS.Rmd R binomial analysis logic.
"""

import argparse
import sys
import numpy as np
import polars as pl
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


def load_interproscan_tsv(path: str) -> pl.DataFrame:
    """Load InterProScan TSV (no header, standard 13-16 column format)."""
    base_columns = [
        "protein_accession", "md5", "length", "analysis",
        "signature_accession", "signature_desc", "start", "stop",
        "score", "status", "date", "interpro_accession", "interpro_desc",
    ]

    with open(path) as f:
        first_line = f.readline()
        if not first_line.strip():
            return pl.DataFrame()
        num_cols = len(first_line.strip().split("\t"))

    if num_cols == 13:
        columns = base_columns
    elif num_cols == 14:
        columns = base_columns + ["go_annotations"]
    elif num_cols >= 15:
        columns = base_columns + ["go_annotations", "pathways"]
    else:
        columns = [f"col_{i}" for i in range(num_cols)]

    return pl.read_csv(path, separator="\t", has_header=False, new_columns=columns[:num_cols])


def run_binomial_analysis(contigs_ips_path: str, genomes_ips_path: str, output_path: str):
    """Run binomial enrichment analysis comparing contig domains to genome domains."""
    contigs_ips = load_interproscan_tsv(contigs_ips_path)
    genomes_ips = load_interproscan_tsv(genomes_ips_path)

    if len(contigs_ips) == 0 or len(genomes_ips) == 0:
        print("Warning: Empty InterProScan input. Writing empty results.")
        pl.DataFrame({
            "V13": [], "V14": [], "p.value": [], "p_adju": [],
            "p_scaled": [], "probability": [], "n": [],
        }).write_csv(output_path)
        return

    # Filter to Pfam only
    if "analysis" in contigs_ips.columns:
        pfam_contigs = contigs_ips.filter(pl.col("analysis") == "Pfam")
    else:
        pfam_contigs = contigs_ips

    if "analysis" in genomes_ips.columns:
        pfam_genomes = genomes_ips.filter(pl.col("analysis") == "Pfam")
    else:
        pfam_genomes = genomes_ips

    if len(pfam_contigs) == 0 or len(pfam_genomes) == 0:
        print("Warning: No Pfam domains found. Writing empty results.")
        pl.DataFrame({
            "V13": [], "V14": [], "p.value": [], "p_adju": [],
            "p_scaled": [], "probability": [], "n": [],
        }).write_csv(output_path)
        return

    # Use interpro_accession (V13) and interpro_desc (V14)
    ipro_col = "interpro_accession" if "interpro_accession" in pfam_contigs.columns else pfam_contigs.columns[min(11, len(pfam_contigs.columns) - 1)]
    desc_col = "interpro_desc" if "interpro_desc" in pfam_contigs.columns else pfam_contigs.columns[min(12, len(pfam_contigs.columns) - 1)]

    # Get unique InterPro ID + description pairs from contigs
    unique_ids = (
        pfam_contigs
        .select([ipro_col, desc_col])
        .unique()
        .filter(pl.col(ipro_col) != "-")
        .filter(pl.col(ipro_col).is_not_null())
    )

    if len(unique_ids) == 0:
        print("Warning: No InterPro accessions in contigs. Writing empty results.")
        pl.DataFrame({
            "V13": [], "V14": [], "p.value": [], "p_adju": [],
            "p_scaled": [], "probability": [], "n": [],
        }).write_csv(output_path)
        return

    # Total Pfam domains in genomes
    domains_total = len(pfam_genomes)

    # Count genome domains by interpro_accession
    genome_ipro_col = "interpro_accession" if "interpro_accession" in pfam_genomes.columns else pfam_genomes.columns[min(11, len(pfam_genomes.columns) - 1)]
    genome_counts = pfam_genomes.group_by(genome_ipro_col).len().rename({"len": "genome_count"})

    # Count contig domains by interpro_accession
    contig_counts = pfam_contigs.group_by(ipro_col).len().rename({"len": "contig_count"})

    results = []
    for row in unique_ids.iter_rows(named=True):
        ipro_id = row[ipro_col]
        ipro_desc = row[desc_col]

        # Successes = count in contigs
        contig_match = contig_counts.filter(pl.col(ipro_col) == ipro_id)
        successes = contig_match["contig_count"][0] if len(contig_match) > 0 else 0

        # Trials = count of that domain in whole genomes
        genome_match = genome_counts.filter(pl.col(genome_ipro_col) == ipro_id)
        trials = genome_match["genome_count"][0] if len(genome_match) > 0 else 0

        if trials == 0 or domains_total == 0:
            continue

        # Probability = trials / total pfam domains in genomes
        probability = trials / domains_total

        # Binomial test
        try:
            result = binomtest(successes, trials, probability)
            p_value = result.pvalue
        except Exception:
            p_value = np.nan

        results.append({
            "V13": ipro_id,
            "V14": ipro_desc,
            "p.value": p_value,
            "probability": probability,
            "n": successes,
        })

    if not results:
        print("Warning: No valid binomial tests. Writing empty results.")
        pl.DataFrame({
            "V13": [], "V14": [], "p.value": [], "p_adju": [],
            "p_scaled": [], "probability": [], "n": [],
        }).write_csv(output_path)
        return

    results_df = pl.DataFrame(results)

    # BH correction
    p_values = results_df["p.value"].to_numpy()
    valid_mask = ~np.isnan(p_values)
    p_adju = np.full_like(p_values, np.nan)
    if valid_mask.sum() > 0:
        _, p_adju[valid_mask], _, _ = multipletests(p_values[valid_mask], method="fdr_bh")

    # Scaled p-value: -10 * log10(p_adju), cap zeros at smallest representable
    with np.errstate(divide="ignore", invalid="ignore"):
        # Replace exact 0 with smallest float to avoid inf
        p_adju_safe = np.where(p_adju == 0, np.finfo(float).tiny, p_adju)
        p_scaled = -10 * np.log10(p_adju_safe)
    p_scaled = np.where(np.isinf(p_scaled) | np.isnan(p_adju), np.nan, p_scaled)

    results_df = results_df.with_columns([
        pl.Series("p_adju", p_adju),
        pl.Series("p_scaled", p_scaled),
    ])

    # Sort by p_scaled descending
    results_df = results_df.sort("p_scaled", descending=True, nulls_last=True)

    results_df.write_csv(output_path)
    print(f"Binomial analysis complete: {len(results_df)} domains tested, "
          f"{results_df.filter(pl.col('p_adju') < 0.05).height} significant (p_adju < 0.05)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RADS Binomial Domain Enrichment Analysis")
    parser.add_argument("--contigs-interproscan", required=True,
                        help="InterProScan results for extracted contigs")
    parser.add_argument("--genomes-interproscan", required=True,
                        help="InterProScan results for whole genomes")
    parser.add_argument("--output", required=True,
                        help="Output CSV path for binomial results")
    args = parser.parse_args()

    run_binomial_analysis(args.contigs_interproscan, args.genomes_interproscan, args.output)
