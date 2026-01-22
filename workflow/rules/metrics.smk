# Metrics Calculation
# ===================
# Calculate pipeline metrics including hits per megabase and discovery rates

import json


rule calculate_genome_sizes:
    """
    Calculate total genome sizes from staged genome files.
    Sums total bases from all input genomes.
    """
    input:
        genomes_dir = f"results/{SAMPLE}/genomes",
        manifest = f"results/{SAMPLE}/genome_manifest.txt"
    output:
        size_file = f"results/{SAMPLE}/metrics/total_genome_size.txt"
    log:
        f"logs/{SAMPLE}/calculate_genome_sizes.log"
    run:
        from pathlib import Path

        genomes_path = Path(input.genomes_dir)
        total_bases = 0
        genome_count = 0

        with open(log[0], "w") as logfile:
            for fna_file in genomes_path.glob("*.fna"):
                file_bases = 0
                with open(fna_file) as f:
                    for line in f:
                        if not line.startswith(">"):
                            file_bases += len(line.strip())
                total_bases += file_bases
                genome_count += 1
                logfile.write(f"{fna_file.name}: {file_bases} bases\n")

            logfile.write(f"\nTotal genomes: {genome_count}\n")
            logfile.write(f"Total bases: {total_bases}\n")

        # Write output
        Path(output.size_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output.size_file, "w") as f:
            f.write(f"{total_bases}\n")


rule calculate_metrics:
    """
    Calculate all pipeline metrics and write to JSON.

    Metrics include:
    - EF_B0058 hits per megabase
    - Discovery rate per contig
    - Discovery rate per genome
    - Query information
    """
    input:
        genome_size = f"results/{SAMPLE}/metrics/total_genome_size.txt",
        blast_results = f"results/{SAMPLE}/blast_results/master_blast.txt",
        manifest = f"results/{SAMPLE}/genome_manifest.txt",
        defense_systems = f"results/{SAMPLE}/defensefinder/defense_finder_systems.tsv",
        contigs = f"results/{SAMPLE}/all_contigs_filtered.fna"
    output:
        metrics = f"results/{SAMPLE}/metrics/pipeline_metrics.json"
    params:
        query_file = config["query_file"],
        sample_name = config["sample_name"]
    log:
        f"logs/{SAMPLE}/calculate_metrics.log"
    run:
        import json
        from pathlib import Path

        with open(log[0], "w") as logfile:
            # Read total genome size
            with open(input.genome_size) as f:
                total_bases = int(f.read().strip())
            total_mb = total_bases / 1_000_000

            # Count BLAST hits
            blast_hits = 0
            blast_file = Path(input.blast_results)
            if blast_file.exists() and blast_file.stat().st_size > 0:
                with open(blast_file) as f:
                    for i, line in enumerate(f):
                        if i > 0:  # Skip header
                            blast_hits += 1

            # Count genomes
            with open(input.manifest) as f:
                total_genomes = len([line for line in f if line.strip()])

            # Count defense systems (skip header)
            defense_systems = 0
            defense_file = Path(input.defense_systems)
            if defense_file.exists() and defense_file.stat().st_size > 0:
                with open(defense_file) as f:
                    for i, line in enumerate(f):
                        if i > 0 and line.strip():  # Skip header
                            defense_systems += 1

            # Count contigs
            contigs_analyzed = 0
            contigs_file = Path(input.contigs)
            if contigs_file.exists():
                with open(contigs_file) as f:
                    for line in f:
                        if line.startswith(">"):
                            contigs_analyzed += 1

            # Extract query name from FASTA header
            query_name = "Unknown"
            query_file = Path(params.query_file)
            if query_file.exists():
                with open(query_file) as f:
                    first_line = f.readline()
                    if first_line.startswith(">"):
                        # Extract ID from header (first word after >)
                        query_name = first_line[1:].split()[0].strip()

            # Calculate metrics
            hits_per_mb = blast_hits / total_mb if total_mb > 0 else 0
            discovery_rate_per_contig = defense_systems / contigs_analyzed if contigs_analyzed > 0 else 0
            discovery_rate_per_genome = defense_systems / total_genomes if total_genomes > 0 else 0

            # Build metrics dictionary
            metrics = {
                "total_input_bases": total_bases,
                "total_input_mb": round(total_mb, 2),
                "blast_hits": blast_hits,
                "hits_per_mb": round(hits_per_mb, 3),
                "defense_systems": defense_systems,
                "contigs_analyzed": contigs_analyzed,
                "discovery_rate_per_contig": round(discovery_rate_per_contig, 3),
                "total_genomes": total_genomes,
                "discovery_rate_per_genome": round(discovery_rate_per_genome, 3),
                "query_file": str(params.query_file),
                "query_name": query_name,
                "sample_name": str(params.sample_name)
            }

            # Log metrics
            logfile.write("Pipeline Metrics Summary:\n")
            logfile.write("=" * 40 + "\n")
            for key, value in metrics.items():
                logfile.write(f"{key}: {value}\n")

            # Write JSON output
            with open(output.metrics, "w") as f:
                json.dump(metrics, f, indent=2)
