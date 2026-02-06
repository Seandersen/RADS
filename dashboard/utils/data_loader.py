"""
Data loader utilities for RADS dashboard.
Loads and parses pipeline output files.
"""

import json
import polars as pl
from pathlib import Path
from typing import Optional


def find_results_dir(base_path: str = "results") -> list[str]:
    """Find all sample result directories."""
    base = Path(base_path)
    if not base.exists():
        return []
    return [d.name for d in base.iterdir() if d.is_dir()]


def load_blast_results(results_dir: str) -> Optional[pl.DataFrame]:
    """
    Load master BLAST results file.

    Returns DataFrame with columns:
    - query_id, subject_id, length, nident, pident, evalue, genome
    """
    blast_file = Path(results_dir) / "blast_results" / "master_blast.txt"
    if not blast_file.exists():
        return None

    try:
        df = pl.read_csv(
            blast_file,
            separator="\t",
            has_header=True,
        )
        return df
    except Exception as e:
        print(f"Error loading BLAST results: {e}")
        return None


def load_genome_manifest(results_dir: str) -> list[str]:
    """Load list of processed genomes."""
    manifest_file = Path(results_dir) / "genome_manifest.txt"
    if not manifest_file.exists():
        return []

    with open(manifest_file) as f:
        return [line.strip() for line in f if line.strip()]


def load_cotranscription_results(results_dir: str) -> Optional[pl.DataFrame]:
    """
    Load co-transcription analysis results.

    Returns DataFrame with columns:
    - blast_hit_id, hit_contig_orf, downstream_orf, strand, distance
    """
    cotx_file = Path(results_dir) / "cotranscription" / "cotranscribed_details.txt"
    if not cotx_file.exists():
        return None

    try:
        df = pl.read_csv(
            cotx_file,
            separator="\t",
            has_header=True,
        )
        # Rename gap_bp to distance for consistency with dashboard code
        if "gap_bp" in df.columns:
            df = df.rename({"gap_bp": "distance"})
        return df
    except Exception as e:
        print(f"Error loading cotranscription results: {e}")
        return None


def load_interproscan_results(results_dir: str) -> Optional[pl.DataFrame]:
    """
    Load InterProScan domain annotation results.

    Standard InterProScan TSV format with columns for protein accession,
    sequence MD5, length, analysis, signature accession, signature description,
    start, stop, score, status, date, interpro accession, interpro description,
    and optionally GO annotations and pathways.
    """
    ips_file = Path(results_dir) / "interproscan_results.tsv"
    if not ips_file.exists():
        return None

    try:
        # Check if file is empty
        if ips_file.stat().st_size == 0:
            return pl.DataFrame({
                "protein_accession": [],
                "analysis": [],
                "signature_desc": [],
                "start": [],
                "stop": [],
            })

        # Read first line to detect number of columns
        with open(ips_file) as f:
            first_line = f.readline()
            num_cols = len(first_line.strip().split('\t'))

        # Define column names based on number of columns
        base_columns = [
            "protein_accession", "md5", "length", "analysis",
            "signature_accession", "signature_desc", "start", "stop",
            "score", "status", "date", "interpro_accession",
            "interpro_desc"
        ]

        if num_cols == 13:
            columns = base_columns
        elif num_cols == 14:
            columns = base_columns + ["go_annotations"]
        elif num_cols >= 15:
            columns = base_columns + ["go_annotations", "pathways"]
        else:
            columns = [f"col_{i}" for i in range(num_cols)]

        df = pl.read_csv(
            ips_file,
            separator="\t",
            has_header=False,
            new_columns=columns,
        )
        return df
    except Exception as e:
        print(f"Error loading InterProScan results: {e}")
        return None


def load_contig_orfs(results_dir: str) -> Optional[pl.DataFrame]:
    """
    Parse ORF information from contig protein FASTA headers.

    Prodigal header format: >contig_orf # start # end # strand # info
    """
    faa_file = Path(results_dir) / "contig_orfs" / "all_contigs.faa"
    if not faa_file.exists():
        return None

    try:
        orfs = []
        with open(faa_file) as f:
            for line in f:
                if line.startswith(">"):
                    # Parse prodigal header
                    parts = line[1:].split("#")
                    if len(parts) >= 4:
                        orf_id = parts[0].strip()
                        start = int(parts[1].strip())
                        end = int(parts[2].strip())
                        strand = int(parts[3].strip())

                        # Extract contig name from ORF ID
                        contig = orf_id.rsplit("_", 1)[0] if "_" in orf_id else orf_id

                        orfs.append({
                            "orf_id": orf_id,
                            "contig": contig,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "length": abs(end - start),
                        })

        return pl.DataFrame(orfs)
    except Exception as e:
        print(f"Error loading contig ORFs: {e}")
        return None


def get_summary_stats(results_dir: str) -> dict:
    """
    Calculate summary statistics from pipeline results.
    """
    stats = {
        "total_genomes": 0,
        "total_blast_hits": 0,
        "genomes_with_hits": 0,
        "total_contigs": 0,
        "total_orfs": 0,
        "cotranscribed_pairs": 0,
        "domain_annotations": 0,
    }

    # Genome count
    genomes = load_genome_manifest(results_dir)
    stats["total_genomes"] = len(genomes)

    # BLAST stats
    blast_df = load_blast_results(results_dir)
    if blast_df is not None and len(blast_df) > 0:
        stats["total_blast_hits"] = len(blast_df)
        stats["genomes_with_hits"] = blast_df["genome"].n_unique()

    # ORF stats
    orfs_df = load_contig_orfs(results_dir)
    if orfs_df is not None:
        stats["total_orfs"] = len(orfs_df)
        stats["total_contigs"] = orfs_df["contig"].n_unique()

    # Cotranscription stats
    cotx_df = load_cotranscription_results(results_dir)
    if cotx_df is not None:
        stats["cotranscribed_pairs"] = len(cotx_df)

    # InterProScan stats
    ips_df = load_interproscan_results(results_dir)
    if ips_df is not None:
        stats["domain_annotations"] = len(ips_df)

    return stats


def load_defensefinder_systems(results_dir: str) -> Optional[pl.DataFrame]:
    """
    Load DefenseFinder defense systems results.

    Handles the macsyfinder output format which includes:
    - Comment lines starting with #
    - Multiple header rows (repeated for each system block)
    - Columns: replicon, hit_id, gene_name, hit_pos, model_fqn, sys_id, etc.

    Returns DataFrame with columns for defense system information including
    type and subtype extracted from model_fqn.
    """
    systems_file = Path(results_dir) / "defensefinder" / "defense_finder_systems.tsv"
    if not systems_file.exists():
        return None

    try:
        # Check if file has content
        if systems_file.stat().st_size == 0:
            return pl.DataFrame({
                "sys_id": [],
                "type": [],
                "subtype": [],
                "hit_id": [],
                "gene_name": [],
            })

        # Read file manually to handle comments and multiple headers
        rows = []
        with open(systems_file) as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue
                # Skip header lines (they start with "replicon")
                if line.startswith("replicon\t"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 6:
                    # Extract type and subtype from model_fqn
                    # Format: defense-finder-models/DefenseFinder/Lamassu-Fam/Lamassu-Protease
                    # or: defense-finder-models/RM/RM/RM_Type_I
                    model_fqn = parts[4] if len(parts) > 4 else ""
                    fqn_parts = model_fqn.split("/")

                    # Extract type (second-to-last part) and subtype (last part)
                    if len(fqn_parts) >= 2:
                        subtype = fqn_parts[-1] if fqn_parts[-1] else "Unknown"
                        # Type is the category - simplify from the path
                        if len(fqn_parts) >= 3:
                            system_type = fqn_parts[-2]  # e.g., "Lamassu-Fam", "RM"
                        else:
                            system_type = subtype
                    else:
                        system_type = "Unknown"
                        subtype = "Unknown"

                    rows.append({
                        "replicon": parts[0],
                        "hit_id": parts[1],
                        "gene_name": parts[2],
                        "hit_pos": parts[3] if len(parts) > 3 else "",
                        "model_fqn": model_fqn,
                        "sys_id": parts[5] if len(parts) > 5 else "",
                        "type": system_type,
                        "subtype": subtype,
                        "hit_status": parts[8] if len(parts) > 8 else "",
                        "hit_score": parts[11] if len(parts) > 11 else "",
                    })

        if not rows:
            return pl.DataFrame({
                "sys_id": [],
                "type": [],
                "subtype": [],
                "hit_id": [],
                "gene_name": [],
            })

        return pl.DataFrame(rows)
    except Exception as e:
        print(f"Error loading DefenseFinder systems: {e}")
        return None


def load_defensefinder_genes(results_dir: str) -> Optional[pl.DataFrame]:
    """
    Load DefenseFinder defense genes results.

    Handles the macsyfinder output format (same as systems file).
    Returns DataFrame with columns for defense gene information.
    """
    genes_file = Path(results_dir) / "defensefinder" / "defense_finder_genes.tsv"
    if not genes_file.exists():
        return None

    try:
        if genes_file.stat().st_size == 0:
            return pl.DataFrame({
                "hit_id": [],
                "gene_name": [],
                "type": [],
                "subtype": [],
            })

        # Read file manually to handle comments and multiple headers
        rows = []
        with open(genes_file) as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue
                # Skip header lines
                if line.startswith("replicon\t") or line.startswith("hit_id\t"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 3:
                    # Extract type and subtype from model_fqn if present
                    model_fqn = parts[4] if len(parts) > 4 else ""
                    fqn_parts = model_fqn.split("/")

                    if len(fqn_parts) >= 2:
                        subtype = fqn_parts[-1] if fqn_parts[-1] else "Unknown"
                        if len(fqn_parts) >= 3:
                            system_type = fqn_parts[-2]
                        else:
                            system_type = subtype
                    else:
                        system_type = "Unknown"
                        subtype = "Unknown"

                    rows.append({
                        "hit_id": parts[1] if len(parts) > 1 else parts[0],
                        "gene_name": parts[2] if len(parts) > 2 else "",
                        "type": system_type,
                        "subtype": subtype,
                        "sys_id": parts[5] if len(parts) > 5 else "",
                    })

        if not rows:
            return pl.DataFrame({
                "hit_id": [],
                "gene_name": [],
                "type": [],
                "subtype": [],
            })

        return pl.DataFrame(rows)
    except Exception as e:
        print(f"Error loading DefenseFinder genes: {e}")
        return None


def load_hit_to_contig_mapping(results_dir: str) -> Optional[pl.DataFrame]:
    """
    Load the hit-to-contig mapping from cotranscription analysis.

    This maps BLAST hit IDs to their corresponding ORFs in extracted contigs.
    """
    mapping_file = Path(results_dir) / "cotranscription" / "hit_to_contig_mapping.tsv"
    if not mapping_file.exists():
        return None

    try:
        df = pl.read_csv(
            mapping_file,
            separator="\t",
            has_header=True,
        )
        return df
    except Exception as e:
        print(f"Error loading hit-to-contig mapping: {e}")
        return None


def load_pipeline_metrics(results_dir: str) -> Optional[dict]:
    """
    Load pipeline metrics JSON file.

    Returns dictionary with metrics including:
    - total_input_bases, total_input_mb
    - blast_hits, hits_per_mb
    - defense_systems, contigs_analyzed, discovery_rate_per_contig
    - total_genomes, discovery_rate_per_genome
    - query_file, query_name, sample_name
    """
    metrics_file = Path(results_dir) / "metrics" / "pipeline_metrics.json"
    if not metrics_file.exists():
        return None

    try:
        with open(metrics_file) as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading pipeline metrics: {e}")
        return None
