# Process ORFs in extracted contigs
# ==================================

rule translate_contigs:
    """
    Run Prodigal on concatenated contigs to find ORFs.
    """
    input:
        contigs = f"results/{SAMPLE}/all_contigs_filtered.fna"
    output:
        proteins = f"results/{SAMPLE}/contig_orfs/all_contigs.faa",
        genes = f"results/{SAMPLE}/contig_orfs/all_contigs_prodigal.txt"
    params:
        mode = config["prodigal"]["mode"]
    log:
        f"logs/{SAMPLE}/translate_contigs.log"
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            prodigal \
                -p {params.mode} \
                -i {input.contigs} \
                -o {output.genes} \
                -a {output.proteins} \
                2>&1 | tee {log}
        else
            touch {output.proteins} {output.genes}
            echo "Empty input contigs, creating empty outputs" >> {log}
        fi
        """


rule prepare_interproscan_input:
    """
    Prepare protein file for InterProScan by removing stop codon asterisks.
    """
    input:
        proteins = f"results/{SAMPLE}/contig_orfs/all_contigs.faa"
    output:
        cleaned = f"results/{SAMPLE}/contig_orfs/interproscan_input.faa"
    log:
        f"logs/{SAMPLE}/prepare_interproscan.log"
    shell:
        """
        sed 's/\\*//g' {input.proteins} > {output.cleaned} 2>> {log}
        echo "Prepared $(grep -c '>' {output.cleaned} || echo 0) sequences for InterProScan" >> {log}
        """


rule run_interproscan:
    """
    Run InterProScan for domain annotation.
    Requires InterProScan to be installed manually.
    """
    input:
        proteins = f"results/{SAMPLE}/contig_orfs/interproscan_input.faa"
    output:
        results = f"results/{SAMPLE}/interproscan_results.tsv"
    params:
        interproscan_path = config["interproscan"]["path"],
        enabled = config["interproscan"]["enabled"]
    log:
        f"logs/{SAMPLE}/interproscan.log"
    conda:
        "../envs/interproscan.yaml"
    shell:
        """
        # Expand tilde in path
        INTERPROSCAN_PATH=$(eval echo {params.interproscan_path})

        if [ "{params.enabled}" = "True" ] || [ "{params.enabled}" = "true" ]; then
            if [ -s {input.proteins} ] && [ -x "$INTERPROSCAN_PATH" ]; then
                bash "$INTERPROSCAN_PATH" \
                    -i {input.proteins} \
                    -f tsv \
                    -o {output.results} \
                    -exclappl PRINTS \
                    2>&1 | tee {log}
            else
                touch {output.results}
                echo "InterProScan not found at $INTERPROSCAN_PATH or empty input, creating empty output" >> {log}
            fi
        else
            touch {output.results}
            echo "InterProScan disabled in config" >> {log}
        fi
        """


# Alternative rule when InterProScan is disabled
rule skip_interproscan:
    """
    Create placeholder output when InterProScan is disabled.
    """
    input:
        proteins = f"results/{SAMPLE}/contig_orfs/interproscan_input.faa"
    output:
        placeholder = f"results/{SAMPLE}/.interproscan_skipped"
    shell:
        """
        touch {output.placeholder}
        """
