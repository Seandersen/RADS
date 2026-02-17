# Binomial Domain Enrichment Analysis
# Compares Pfam domain frequencies in RADS contigs vs. whole genomes


def get_whole_genome_interproscan(wildcards):
    """Determine source of whole-genome InterProScan results."""
    user_path = config.get("binomial", {}).get("whole_genome_interproscan", "")
    if user_path:
        return user_path
    return f"results/{SAMPLE}/binomial/interproscan_wholegenomes.tsv"


rule extract_genomes_with_hits:
    """Combine genome FASTA files for genomes that have BLAST hits."""
    input:
        blast=f"results/{SAMPLE}/blast_results/master_blast.txt",
        genomes_dir=f"results/{SAMPLE}/genomes",
    output:
        combined=f"results/{SAMPLE}/binomial/genomes_with_hits.fna",
    run:
        import polars as pl
        from pathlib import Path

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


rule run_whole_genome_interproscan:
    """Run InterProScan on translated proteins from genomes with BLAST hits."""
    input:
        faa=f"results/{SAMPLE}/binomial/genomes_with_hits_cleaned.faa",
    output:
        tsv=f"results/{SAMPLE}/binomial/interproscan_wholegenomes.tsv",
    params:
        interproscan_path=config.get("interproscan", {}).get("path", "interproscan.sh"),
        applications=config.get("interproscan", {}).get("applications", "Pfam"),
    threads: config.get("diamond", {}).get("threads", 8)
    log:
        f"logs/{SAMPLE}/binomial/interproscan_wholegenomes.log",
    shell:
        """
        {params.interproscan_path} \
            -i {input.faa} \
            -o {output.tsv} \
            -f TSV \
            -appl {params.applications} \
            --cpu {threads} \
            --disable-precalc \
            2>&1 | tee {log}
        """


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
