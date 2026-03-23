# Download genomes from NCBI using datasets CLI or efetch
# =======================================================

rule download_genomes:
    """
    Download genomes from NCBI.
    Supports two modes:
    1. Taxon-based download using datasets CLI (default)
    2. Accession file download using efetch (when accession_file is set)
    """
    output:
        directory(f"results/{SAMPLE}/downloaded_genomes/ncbi_dataset/data")
    params:
        taxon = config["download"].get("taxon", ""),
        source = config["download"].get("source", "refseq"),
        assembly_level = config["download"].get("assembly_level", "complete"),
        max_genomes = config["download"].get("max_genomes", 0),
        random_sample = config["download"].get("random_sample", False),
        accession_file = config["download"].get("accession_file", ""),
        outdir = f"results/{SAMPLE}/downloaded_genomes"
    log:
        f"logs/{SAMPLE}/download_genomes.log"
    benchmark:
        f"logs/{SAMPLE}/benchmarks/download_genomes.txt"
    shell:
        """
        mkdir -p {params.outdir}/ncbi_dataset/data

        if [ -n "{params.accession_file}" ] && [ -f "{params.accession_file}" ]; then
            # Mode 2: Download by accession file using efetch
            echo "Downloading from accession file: {params.accession_file}" > {log}

            total=$(wc -l < "{params.accession_file}")
            count=0

            while IFS= read -r acc || [ -n "$acc" ]; do
                if [ -z "$acc" ]; then continue; fi
                count=$((count + 1))
                echo "[$count/$total] Downloading $acc..." >> {log}

                # Create directory for this accession
                acc_dir="{params.outdir}/ncbi_dataset/data/$acc"
                mkdir -p "$acc_dir"

                # Download using efetch
                efetch -db nuccore -id "$acc" -format fasta > "$acc_dir/$acc.fna" 2>> {log}

                # Check if download succeeded
                if [ ! -s "$acc_dir/$acc.fna" ]; then
                    echo "WARNING: Failed to download $acc" >> {log}
                    rm -rf "$acc_dir"
                fi

                # Rate limiting to avoid NCBI throttling
                sleep 0.4
            done < "{params.accession_file}"

            echo "Download complete. Successfully downloaded $count sequences." >> {log}

        else
            # Mode 1: Taxon-based download using datasets CLI
            if [ {params.max_genomes} -gt 0 ]; then
                # Get list of accessions (avoid SIGPIPE by saving to temp file first)
                datasets summary genome taxon "{params.taxon}" \
                    --assembly-source {params.source} \
                    --assembly-level {params.assembly_level} \
                    --as-json-lines > {params.outdir}/all_records.jsonl 2>> {log} || true

                # Extract N accessions (random or first N)
                if [ "{params.random_sample}" = "True" ] || [ "{params.random_sample}" = "true" ]; then
                    shuf -n {params.max_genomes} {params.outdir}/all_records.jsonl | \
                        jq -r '.accession' > {params.outdir}/accessions.txt
                    echo "Random sample of {params.max_genomes} genomes selected" >> {log}
                else
                    head -n {params.max_genomes} {params.outdir}/all_records.jsonl | \
                        jq -r '.accession' > {params.outdir}/accessions.txt
                    echo "First {params.max_genomes} genomes selected" >> {log}
                fi

                rm -f {params.outdir}/all_records.jsonl

                echo "Downloading $(wc -l < {params.outdir}/accessions.txt) genomes..." >> {log}

                # Download by accession list
                datasets download genome accession \
                    --inputfile {params.outdir}/accessions.txt \
                    --include genome \
                    --filename {params.outdir}/genomes.zip \
                    2>&1 | tee -a {log}

                rm -f {params.outdir}/accessions.txt
            else
                # No limit - download all matching
                datasets download genome taxon "{params.taxon}" \
                    --assembly-source {params.source} \
                    --assembly-level {params.assembly_level} \
                    --include genome \
                    --filename {params.outdir}/genomes.zip \
                    2>&1 | tee {log}
            fi

            unzip -o {params.outdir}/genomes.zip -d {params.outdir} >> {log} 2>&1
            rm -f {params.outdir}/genomes.zip
        fi
        """


# Helper rule to skip download when using local genomes
rule skip_download:
    """
    Placeholder rule when download is disabled.
    Creates a marker file indicating local genomes should be used.
    """
    output:
        touch(f"results/{SAMPLE}/.using_local_genomes")
