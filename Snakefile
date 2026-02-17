# RADS - Recombinase Associated Defense Search
# Snakemake Pipeline
# ================================================

import os
from pathlib import Path

configfile: "config/config.yaml"

# Configuration
SAMPLE = config["sample_name"]
QUERY = config["query_file"]
UPSTREAM = config["upstream_nt"]
DOWNSTREAM = config["downstream_nt"]

# Determine genome source path
if config["download"]["enabled"]:
    GENOMES_PATH = f"results/{SAMPLE}/downloaded_genomes/ncbi_dataset/data"
else:
    GENOMES_PATH = config["genomes_path"]


# =============================================================================
# IMPORTANT: rule all MUST be the first rule to be the default target
# =============================================================================
def get_all_targets(wildcards=None):
    """Collect all target files, including optional binomial analysis."""
    targets = [
        # Phase 1 outputs (steps 1-4)
        f"results/{SAMPLE}/blast_results/master_blast.txt",
        # Phase 2 outputs (steps 5-7)
        f"results/{SAMPLE}/interproscan_results.tsv",
        f"results/{SAMPLE}/cotranscription/cotranscribed_sequences.faa",
        # Phase 3 outputs (DefenseFinder and metrics)
        f"results/{SAMPLE}/defensefinder/defense_finder_systems.tsv",
        f"results/{SAMPLE}/metrics/pipeline_metrics.json",
    ]
    # Defense scores (proximity/density to known defense systems)
    targets.append(f"results/{SAMPLE}/defense_scores.tsv")
    # Optional: Binomial domain enrichment analysis
    if config.get("binomial", {}).get("enabled", False):
        targets.append(f"results/{SAMPLE}/BinomialAnalysis.csv")
    return targets

rule all:
    """Final target rule - requests outputs from all phases."""
    input:
        get_all_targets,


# =============================================================================
# Checkpoint and helper functions for dynamic genome discovery
# =============================================================================
checkpoint discover_genomes:
    """Discover all genome directories after staging."""
    input:
        genomes_dir = f"results/{SAMPLE}/genomes"
    output:
        manifest = f"results/{SAMPLE}/genome_manifest.txt"
    run:
        genome_dir = Path(input.genomes_dir)
        genomes = [f.stem for f in genome_dir.glob("*.fna")]
        with open(output.manifest, "w") as f:
            for g in genomes:
                f.write(f"{g}\n")


def get_all_genomes(wildcards):
    """Get list of all genome IDs from checkpoint."""
    checkpoint_output = checkpoints.discover_genomes.get(**wildcards).output.manifest
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    return genomes


def get_all_blast_results(wildcards):
    """Aggregate function for all BLAST results."""
    genomes = get_all_genomes(wildcards)
    return expand(
        f"results/{SAMPLE}/blast_results/{{genome}}_blast.txt",
        genome=genomes
    )


def get_all_contigs(wildcards):
    """Aggregate function for all extracted contigs."""
    genomes = get_all_genomes(wildcards)
    return expand(
        f"results/{SAMPLE}/contigs/{{genome}}_contigs.fna",
        genome=genomes
    )


# =============================================================================
# Include rule modules
# =============================================================================
# Phase 1: Core pipeline
include: "workflow/rules/download_genomes.smk"
include: "workflow/rules/stage_genomes.smk"
include: "workflow/rules/translate.smk"
include: "workflow/rules/build_databases.smk"
include: "workflow/rules/blast_search.smk"
# Phase 2: Analysis
include: "workflow/rules/extract_contigs.smk"
include: "workflow/rules/process_orfs.smk"
include: "workflow/rules/cotranscription.smk"
# Phase 3: DefenseFinder and metrics
include: "workflow/rules/defensefinder.smk"
include: "workflow/rules/metrics.smk"
# Phase 4: Optional analyses
include: "workflow/rules/binomial_analysis.smk"
include: "workflow/rules/defense_score.smk"
