# Binomial Domain Enrichment Analysis
# Compares Pfam domain frequencies in RADS contigs vs. whole genomes
#
# When binomial.whole_genome_interproscan is not pre-supplied, the whole-genome
# InterProScan step is parallelized:
#   1. split_proteins_for_binomial  (checkpoint) — seqkit splits cleaned FAA into chunks
#   2. run_interproscan_chunk       (scatter)     — one SLURM job per chunk
#   3. aggregate_interproscan_chunks (gather)     — cat all chunk TSVs into final TSV


import glob as _glob
from pathlib import Path


def get_whole_genome_interproscan(wildcards):
    """Determine source of whole-genome InterProScan results."""
    user_path = config.get("binomial", {}).get("whole_genome_interproscan", "")
    if user_path:
        return user_path
    return f"results/{SAMPLE}/binomial/interproscan_wholegenomes.tsv"


# =============================================================================
# Step 1: Extract and translate genomes that have BLAST hits
# =============================================================================

rule extract_genomes_with_hits:
    """Combine genome FASTA files for genomes that have BLAST hits."""
    input:
        blast=f"results/{SAMPLE}/blast_results/master_blast.txt",
        genomes_dir=f"results/{SAMPLE}/genomes",
    output:
        combined=f"results/{SAMPLE}/binomial/genomes_with_hits.fna",
    run:
        import polars as pl

        blast = pl.read_csv(input.blast, separator="\t", has_header=True)
        genome_ids = blast["genome"].unique().to_list()

        genomes_dir = Path(input.genomes_dir)
        with open(output.combined, "w") as out:
            for gid in genome_ids:
                fna = genomes_dir / f"{gid}.fna"
                if fna.exists():
                    with open(fna) as f:
                        out.write(f.read())


rule translate_genomes_for_binomial:
    """Translate combined genomes with Prodigal to get protein sequences."""
    input:
        fna=f"results/{SAMPLE}/binomial/genomes_with_hits.fna",
    output:
        faa=f"results/{SAMPLE}/binomial/genomes_with_hits.faa",
        genes=f"results/{SAMPLE}/binomial/genomes_with_hits_prodigal.txt",
    params:
        mode=config["prodigal"]["mode"],
    log:
        f"logs/{SAMPLE}/binomial/translate_genomes.log",
    shell:
        """
        prodigal \
            -p {params.mode} \
            -i {input.fna} \
            -o {output.genes} \
            -a {output.faa} \
            2>&1 | tee {log}
        """


rule clean_proteins_for_binomial:
    """Strip stop codon asterisks from Prodigal output for InterProScan compatibility."""
    input:
        faa=f"results/{SAMPLE}/binomial/genomes_with_hits.faa",
    output:
        cleaned=f"results/{SAMPLE}/binomial/genomes_with_hits_cleaned.faa",
    log:
        f"logs/{SAMPLE}/binomial/clean_proteins.log",
    shell:
        """
        sed 's/\\*//g' {input.faa} > {output.cleaned} 2>> {log}
        echo "Prepared $(grep -c '>' {output.cleaned} || echo 0) sequences for InterProScan" >> {log}
        """


# =============================================================================
# Step 2: Split cleaned proteins into chunks for parallel InterProScan
# =============================================================================

checkpoint split_proteins_for_binomial:
    """
    Split the cleaned whole-genome protein file into equal-sized chunks so
    each chunk can be submitted as its own SLURM job.

    seqkit split2 names chunks: {stem}.part_001.faa, {stem}.part_002.faa, …
    If the input is empty a single empty placeholder chunk is created so the
    downstream gather rule always has something to aggregate.
    """
    input:
        faa=f"results/{SAMPLE}/binomial/genomes_with_hits_cleaned.faa",
    output:
        chunk_dir=directory(f"results/{SAMPLE}/binomial/chunks"),
    params:
        chunk_size=config.get("binomial", {}).get("interproscan_chunk_size", 50000),
    log:
        f"logs/{SAMPLE}/binomial/split_proteins.log",
    shell:
        """
        mkdir -p {output.chunk_dir}

        if [ -s {input.faa} ]; then
            seqkit split2 \
                --by-size {params.chunk_size} \
                --out-dir {output.chunk_dir} \
                {input.faa} \
                2>&1 | tee {log}
            echo "Split into $(ls {output.chunk_dir}/*.faa | wc -l) chunks of up to {params.chunk_size} sequences" >> {log}
        else
            # No proteins — create an empty placeholder so the DAG can proceed
            touch {output.chunk_dir}/empty.faa
            echo "Input was empty; created placeholder chunk" | tee {log}
        fi
        """


# =============================================================================
# Helpers: discover chunks after the checkpoint resolves
# =============================================================================

def _get_chunk_stems(wildcards):
    """Return the file stems of all .faa chunks produced by the checkpoint."""
    chunk_dir = checkpoints.split_proteins_for_binomial.get(**wildcards).output.chunk_dir
    faa_files = sorted(_glob.glob(f"{chunk_dir}/*.faa"))
    return [Path(f).stem for f in faa_files]


def get_interproscan_chunk_inputs(wildcards):
    """Return paths to all chunk FAA files (used by the scatter rule)."""
    stems = _get_chunk_stems(wildcards)
    return expand(
        f"results/{SAMPLE}/binomial/chunks/{{chunk}}.faa",
        chunk=stems,
    )


def get_interproscan_chunk_results(wildcards):
    """Return paths to all chunk TSV outputs (used by the gather rule)."""
    stems = _get_chunk_stems(wildcards)
    return expand(
        f"results/{SAMPLE}/binomial/chunks/{{chunk}}.tsv",
        chunk=stems,
    )


# =============================================================================
# Step 3: Run InterProScan on each chunk (scatter — one SLURM job per chunk)
# =============================================================================

rule run_interproscan_chunk:
    """
    Run InterProScan on a single protein chunk.
    One job is submitted per chunk, so all chunks run in parallel on SLURM.
    """
    input:
        faa=f"results/{SAMPLE}/binomial/chunks/{{chunk}}.faa",
    output:
        tsv=f"results/{SAMPLE}/binomial/chunks/{{chunk}}.tsv",
    params:
        interproscan_path=config.get("interproscan", {}).get("path", "interproscan.sh"),
        applications=config.get("interproscan", {}).get("applications", "Pfam"),
    threads:
        config.get("interproscan", {}).get("threads", 8)
    log:
        f"logs/{SAMPLE}/binomial/interproscan_chunk_{{chunk}}.log",
    shell:
        """
        INTERPROSCAN_PATH=$(eval echo {params.interproscan_path})

        if [ -s {input.faa} ] && [ -x "$INTERPROSCAN_PATH" ]; then
            bash "$INTERPROSCAN_PATH" \
                -i {input.faa} \
                -f tsv \
                -o {output.tsv} \
                -appl {params.applications} \
                -cpu {threads} \
                --disable-precalc \
                2>&1 | tee {log}
        else
            touch {output.tsv}
            echo "Skipped: empty chunk or InterProScan not found at $INTERPROSCAN_PATH" | tee {log}
        fi
        """


# =============================================================================
# Step 4: Gather — concatenate all chunk TSVs into the final whole-genome TSV
# =============================================================================

rule aggregate_interproscan_chunks:
    """
    Concatenate all per-chunk InterProScan TSV files into a single whole-genome
    InterProScan results file for use in the binomial analysis.
    """
    input:
        chunks=get_interproscan_chunk_results,
    output:
        tsv=f"results/{SAMPLE}/binomial/interproscan_wholegenomes.tsv",
    log:
        f"logs/{SAMPLE}/binomial/aggregate_interproscan.log",
    run:
        with open(output.tsv, "w") as out_fh, open(log[0], "w") as logf:
            total_lines = 0
            for tsv_file in sorted(input.chunks):
                with open(tsv_file) as fh:
                    for line in fh:
                        out_fh.write(line)
                        total_lines += 1
            logf.write(
                f"Aggregated {total_lines} lines from {len(input.chunks)} chunks\n"
            )


# =============================================================================
# Step 5: Binomial enrichment analysis (uses whole-genome TSV from above)
# =============================================================================

rule run_binomial_analysis:
    """Run binomial domain enrichment analysis."""
    input:
        contigs_ips=f"results/{SAMPLE}/interproscan_results.tsv",
        genomes_ips=get_whole_genome_interproscan,
    output:
        csv=f"results/{SAMPLE}/BinomialAnalysis.csv",
    log:
        f"logs/{SAMPLE}/binomial/binomial_analysis.log",
    shell:
        """
        python workflow/scripts/binomial_analysis.py \
            --contigs-interproscan {input.contigs_ips} \
            --genomes-interproscan {input.genomes_ips} \
            --output {output.csv} \
            2>&1 | tee {log}
        """
