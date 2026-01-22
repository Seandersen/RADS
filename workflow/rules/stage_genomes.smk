# Stage genomes - collect and organize input genomes
# ==================================================

import os
from pathlib import Path


def get_genome_source(wildcards):
    """Determine input based on whether download is enabled."""
    if config["download"]["enabled"]:
        return f"results/{SAMPLE}/downloaded_genomes/ncbi_dataset/data"
    else:
        return config["genomes_path"]


rule stage_genomes:
    """
    Stage genomes by copying/linking them to the working directory.
    Handles both downloaded and local genome sources.
    Corresponds to RADS.sh mvgenomes() function.
    """
    input:
        source_dir = get_genome_source
    output:
        genomes_dir = directory(f"results/{SAMPLE}/genomes")
    log:
        f"logs/{SAMPLE}/stage_genomes.log"
    benchmark:
        f"logs/{SAMPLE}/benchmarks/stage_genomes.txt"
    run:
        import shutil
        from pathlib import Path

        source = Path(input.source_dir)
        dest = Path(output.genomes_dir)
        dest.mkdir(parents=True, exist_ok=True)

        with open(log[0], "w") as logfile:
            genome_count = 0

            # Handle NCBI dataset structure (directories with .fna files)
            for genome_dir in source.iterdir():
                if genome_dir.is_dir():
                    for fna_file in genome_dir.glob("*.fna"):
                        # Use directory name as genome ID
                        genome_id = genome_dir.name
                        dest_file = dest / f"{genome_id}.fna"
                        shutil.copy2(fna_file, dest_file)
                        logfile.write(f"Staged: {fna_file} -> {dest_file}\n")
                        genome_count += 1
                elif genome_dir.suffix == ".fna":
                    # Handle flat directory structure
                    dest_file = dest / genome_dir.name
                    shutil.copy2(genome_dir, dest_file)
                    logfile.write(f"Staged: {genome_dir} -> {dest_file}\n")
                    genome_count += 1

            logfile.write(f"\nTotal genomes staged: {genome_count}\n")
