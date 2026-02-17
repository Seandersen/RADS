#!/usr/bin/env python3
"""
defense_score.py — Score co-transcribed genes for defense island proximity.

Reads existing RADS pipeline outputs and scores each co-transcribed gene
based on its proximity to DefenseFinder-identified defense genes and the
local density of defense genes. A low score means the gene is isolated
from known defense systems and would likely be missed by traditional
defense island searching.
"""

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict


def parse_prodigal_headers(faa_path):
    """Parse Prodigal FASTA headers to extract ORF coordinates.

    Returns dict: orf_id -> {contig, start, stop, strand}
    """
    orfs = {}
    with open(faa_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            parts = line[1:].split(" # ")
            if len(parts) < 4:
                continue
            orf_id = parts[0].strip()
            start = int(parts[1].strip())
            stop = int(parts[2].strip())
            strand = int(parts[3].strip())

            # Derive contig region from ORF ID
            # Format: ContigAcc_start-stop:._N
            contig = extract_contig(orf_id)

            orfs[orf_id] = {
                "contig": contig,
                "start": start,
                "stop": stop,
                "strand": strand,
            }
    return orfs


def extract_contig(orf_id):
    """Extract contig region identifier from an ORF ID.

    ORF ID format: ContigAcc_start-stop:._N
    Contig region: ContigAcc_start-stop
    """
    match = re.match(r"(.+?_\d+-\d+):\._\d+$", orf_id)
    if match:
        return match.group(1)
    # Fallback: everything before the last _N
    match = re.match(r"(.+)_\d+$", orf_id)
    if match:
        return match.group(1)
    return orf_id


def parse_cotranscribed(path):
    """Parse cotranscribed_details.txt.

    Returns list of dicts with keys: blast_hit_id, hit_contig_orf,
    downstream_orf, strand, gap_bp
    """
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def parse_defense_genes(path):
    """Parse defense_finder_genes.tsv.

    Handles both the simple TSV format (from empty-file fallback) and the
    real DefenseFinder/macsyfinder output which may have leading blank lines
    and repeated header rows. Extracts type/subtype from model_fqn if
    type/subtype columns are not present.

    Returns list of dicts with at least keys: hit_id, type, subtype.
    Returns None if the file is missing or empty (header-only),
    indicating DefenseFinder didn't produce results.
    """
    if not os.path.exists(path):
        return None

    rows = []
    header = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            # Skip repeated header lines
            if parts == header:
                continue
            row = dict(zip(header, parts))

            # Extract type/subtype from model_fqn if not already present
            if "type" not in row or "subtype" not in row:
                fqn = row.get("model_fqn", "")
                fqn_parts = fqn.split("/")
                if len(fqn_parts) >= 3:
                    row["type"] = fqn_parts[-2]
                    row["subtype"] = fqn_parts[-1]
                elif len(fqn_parts) >= 2:
                    row["type"] = fqn_parts[-1]
                    row["subtype"] = fqn_parts[-1]
                else:
                    row["type"] = ""
                    row["subtype"] = ""

            rows.append(row)

    if len(rows) == 0:
        return None

    return rows


def parse_interproscan(path):
    """Parse InterProScan results TSV (no header).

    Returns dict: orf_id -> set of InterPro accessions
    """
    domains = defaultdict(set)
    if not path or not os.path.exists(path):
        return domains

    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 12:
                continue
            orf_id = row[0]
            interpro_acc = row[11] if len(row) > 11 else ""
            if interpro_acc and interpro_acc != "-":
                domains[orf_id].add(interpro_acc)
    return domains



def midpoint(orf_info):
    """Calculate the midpoint of an ORF."""
    return (orf_info["start"] + orf_info["stop"]) / 2.0


def distance_between(orf_a, orf_b):
    """Calculate distance in bp between two ORFs (edge-to-edge, 0 if overlapping)."""
    if orf_a["stop"] < orf_b["start"]:
        return orf_b["start"] - orf_a["stop"]
    elif orf_b["stop"] < orf_a["start"]:
        return orf_a["start"] - orf_b["stop"]
    else:
        return 0  # overlapping


def compute_scores(cotranscribed, orf_coords, defense_genes, interpro_domains,
                   scale, window, w_prox, w_dens):
    """Compute defense scores for all co-transcribed genes.

    Returns list of output row dicts.
    """
    # Determine if DefenseFinder produced results
    df_missing = defense_genes is None

    # Index defense genes by contig
    defense_by_contig = defaultdict(list)
    if not df_missing:
        for dg in defense_genes:
            hit_id = dg["hit_id"]
            if hit_id in orf_coords:
                contig = orf_coords[hit_id]["contig"]
                defense_by_contig[contig].append({
                    "hit_id": hit_id,
                    "type": dg.get("type", ""),
                    "subtype": dg.get("subtype", ""),
                    "coords": orf_coords[hit_id],
                })

    # Set of defense gene IDs for quick lookup
    defense_ids = set()
    if not df_missing:
        for dg in defense_genes:
            defense_ids.add(dg["hit_id"])

    results = []
    for row in cotranscribed:
        downstream_orf = row["downstream_orf"]
        blast_hit_id = row["blast_hit_id"]

        # Get coordinates for the downstream orf
        if downstream_orf not in orf_coords:
            continue

        orf_info = orf_coords[downstream_orf]
        contig = orf_info["contig"]

        out = {
            "downstream_orf": downstream_orf,
            "blast_hit_id": blast_hit_id,
            "contig": contig,
            "orf_start": orf_info["start"],
            "orf_stop": orf_info["stop"],
            "strand": orf_info["strand"],
        }

        if df_missing:
            # DefenseFinder didn't run — mark all as NA
            out["nearest_defense_gene"] = "NA"
            out["nearest_defense_type"] = "NA"
            out["nearest_distance_bp"] = "NA"
            out["defense_genes_in_window"] = "NA"
            out["proximity_score"] = "NA"
            out["density_score"] = "NA"
            out["defense_score"] = "NA"
        else:
            # Check if this gene IS a defense gene
            is_defense = downstream_orf in defense_ids

            contig_defense = defense_by_contig.get(contig, [])

            if is_defense:
                nearest_dist = 0
                nearest_id = downstream_orf
                # Find type/subtype for this gene
                nearest_type = ""
                for dg in contig_defense:
                    if dg["hit_id"] == downstream_orf:
                        nearest_type = dg.get("type", "")
                        break
            elif contig_defense:
                nearest_dist = float("inf")
                nearest_id = "none"
                nearest_type = ""
                for dg in contig_defense:
                    d = distance_between(orf_info, dg["coords"])
                    if d < nearest_dist:
                        nearest_dist = d
                        nearest_id = dg["hit_id"]
                        nearest_type = dg.get("type", "")
            else:
                nearest_dist = float("inf")
                nearest_id = "none"
                nearest_type = ""

            # Count defense genes within window
            count_in_window = 0
            orf_mid = midpoint(orf_info)
            for dg in contig_defense:
                dg_mid = midpoint(dg["coords"])
                if abs(orf_mid - dg_mid) <= window / 2.0:
                    count_in_window += 1

            # Proximity score
            if nearest_dist == 0:
                proximity = 1.0
            elif nearest_dist == float("inf"):
                proximity = 0.0
            else:
                proximity = math.exp(-nearest_dist / scale)

            # Density score (cap at 10)
            max_genes = 10
            density = min(count_in_window / max_genes, 1.0)

            # Composite
            defense_score = (w_prox * proximity) + (w_dens * density)

            out["nearest_defense_gene"] = nearest_id
            out["nearest_defense_type"] = nearest_type if nearest_type else "none"
            out["nearest_distance_bp"] = int(nearest_dist) if nearest_dist != float("inf") else 0
            out["defense_genes_in_window"] = count_in_window
            out["proximity_score"] = round(proximity, 4)
            out["density_score"] = round(density, 4)
            out["defense_score"] = round(defense_score, 4)

        # Domain annotations
        orf_domains = interpro_domains.get(downstream_orf, set())
        out["interpro_domains"] = ";".join(sorted(orf_domains)) if orf_domains else ""

        results.append(out)

    return results


OUTPUT_COLUMNS = [
    "downstream_orf",
    "blast_hit_id",
    "contig",
    "orf_start",
    "orf_stop",
    "strand",
    "nearest_defense_gene",
    "nearest_defense_type",
    "nearest_distance_bp",
    "defense_genes_in_window",
    "proximity_score",
    "density_score",
    "defense_score",
    "interpro_domains",
]


def main():
    parser = argparse.ArgumentParser(
        description="Score co-transcribed genes for defense island proximity."
    )
    parser.add_argument(
        "--cotranscribed", required=True,
        help="Path to cotranscribed_details.txt"
    )
    parser.add_argument(
        "--defense-genes", required=True,
        help="Path to defense_finder_genes.tsv"
    )
    parser.add_argument(
        "--contigs-faa", required=True,
        help="Path to all_contigs.faa (Prodigal output)"
    )
    parser.add_argument(
        "--interproscan", default=None,
        help="Path to interproscan_results.tsv (optional)"
    )
    parser.add_argument(
        "--output", default="defense_scores.tsv",
        help="Output TSV path (default: defense_scores.tsv)"
    )
    parser.add_argument(
        "--window", type=int, default=10000,
        help="Window size in bp for density calculation (default: 10000)"
    )
    parser.add_argument(
        "--scale", type=int, default=2000,
        help="Exponential decay scale in bp for proximity (default: 2000)"
    )
    parser.add_argument(
        "--w-proximity", type=float, default=0.6,
        help="Weight for proximity score (default: 0.6)"
    )
    parser.add_argument(
        "--w-density", type=float, default=0.4,
        help="Weight for density score (default: 0.4)"
    )

    args = parser.parse_args()

    # Parse inputs
    cotranscribed = parse_cotranscribed(args.cotranscribed)
    if not cotranscribed:
        # Empty input — write header-only output
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
            writer.writeheader()
        print(f"No co-transcribed genes found. Wrote header-only {args.output}")
        sys.exit(0)

    orf_coords = parse_prodigal_headers(args.contigs_faa)
    defense_genes = parse_defense_genes(args.defense_genes)
    interpro_domains = parse_interproscan(args.interproscan)

    if defense_genes is None:
        print("WARNING: DefenseFinder results missing or empty — scoring columns set to NA",
              file=sys.stderr)

    # Compute scores
    results = compute_scores(
        cotranscribed, orf_coords, defense_genes, interpro_domains,
        args.scale, args.window,
        args.w_proximity, args.w_density,
    )

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"Wrote {len(results)} scored genes to {args.output}")

    # Summary statistics
    if defense_genes is not None and results:
        scores = [r["defense_score"] for r in results if r["defense_score"] != "NA"]
        if scores:
            print(f"  Score range: {min(scores):.4f} - {max(scores):.4f}")
            print(f"  Mean score:  {sum(scores)/len(scores):.4f}")
            low = sum(1 for s in scores if s < 0.1)
            high = sum(1 for s in scores if s >= 0.5)
            print(f"  Low (<0.1):  {low}  |  High (>=0.5): {high}")



if __name__ == "__main__":
    main()
