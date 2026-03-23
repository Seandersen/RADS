# BLAST search using Diamond
# ==========================

rule blast_search:
    """
    Run Diamond BLASTP against a genome database.
    Corresponds to RADS.sh blast() function.
    Handles empty database files gracefully.
    """
    input:
        db = f"results/{SAMPLE}/diamond_dbs/{{genome}}.dmnd",
        query = QUERY
    output:
        hits = f"results/{SAMPLE}/blast_results/{{genome}}_blast.txt"
    params:
        identity = config["diamond"]["identity"],
        max_target_seqs = config["diamond"]["max_target_seqs"],
        block_size = config["diamond"].get("block_size", 0)
    threads:
        config["diamond"]["threads"]
    log:
        f"logs/{SAMPLE}/blast/{{genome}}.log"
    benchmark:
        f"logs/{SAMPLE}/benchmarks/blast_{{genome}}.txt"
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        # Check if database is valid (not empty marker file)
        if [ -s {input.db} ] && file {input.db} | grep -q "data"; then
            # Build max_target_seqs option only if > 0
            MAX_TARGETS_OPT=""
            if [ {params.max_target_seqs} -gt 0 ]; then
                MAX_TARGETS_OPT="--max-target-seqs {params.max_target_seqs}"
            fi

            # Build block_size option only if > 0
            BLOCK_SIZE_OPT=""
            if [ {params.block_size} -gt 0 ]; then
                BLOCK_SIZE_OPT="--block-size {params.block_size}"
            fi

            diamond blastp \
                -d {input.db} \
                --query {input.query} \
                --threads {threads} \
                --out {output.hits} \
                --outfmt 6 qseqid sseqid length nident pident evalue \
                $MAX_TARGETS_OPT \
                $BLOCK_SIZE_OPT \
                --id {params.identity} \
                2>&1 | tee {log}
        else
            echo "Empty or invalid database, creating empty output" | tee {log}
            touch {output.hits}
        fi
        """


rule aggregate_blast_results:
    """
    Combine all BLAST results into a master file.
    Adds genome ID as a column for downstream analysis.
    """
    input:
        blast_files = get_all_blast_results
    output:
        master = f"results/{SAMPLE}/blast_results/master_blast.txt"
    log:
        f"logs/{SAMPLE}/aggregate_blast.log"
    run:
        with open(output.master, "w") as out, open(log[0], "w") as logfile:
            # Write header
            out.write("query_id\tsubject_id\tlength\tnident\tpident\tevalue\tgenome\n")

            total_hits = 0
            for blast_file in input.blast_files:
                # Extract genome ID from filename
                genome_id = Path(blast_file).stem.replace("_blast", "")

                with open(blast_file) as f:
                    for line in f:
                        if line.strip():
                            out.write(f"{line.strip()}\t{genome_id}\n")
                            total_hits += 1

                logfile.write(f"Processed: {blast_file}\n")

            logfile.write(f"\nTotal hits aggregated: {total_hits}\n")
