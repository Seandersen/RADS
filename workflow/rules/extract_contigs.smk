# Extract contigs around BLAST hits
# ==================================

rule extract_orf_ids:
    """
    Extract ORF IDs from BLAST results for a genome.
    """
    input:
        blast = f"results/{SAMPLE}/blast_results/{{genome}}_blast.txt"
    output:
        orf_ids = f"results/{SAMPLE}/orf_extraction/{{genome}}_orf_ids.txt"
    log:
        f"logs/{SAMPLE}/extract_orf/{{genome}}.log"
    shell:
        """
        if [ -s {input.blast} ]; then
            cut -f2 {input.blast} > {output.orf_ids} 2>> {log}
        else
            touch {output.orf_ids}
            echo "Empty BLAST file, creating empty output" >> {log}
        fi
        """


rule extract_coordinates:
    """
    Extract coordinates from ORF headers in translated proteins.
    """
    input:
        orf_ids = f"results/{SAMPLE}/orf_extraction/{{genome}}_orf_ids.txt",
        proteins = f"results/{SAMPLE}/translated/{{genome}}.faa"
    output:
        coordinates = f"results/{SAMPLE}/orf_extraction/{{genome}}_coordinates.txt",
        formatted = f"results/{SAMPLE}/orf_extraction/{{genome}}_formatted_coords.txt",
        final = f"results/{SAMPLE}/orf_extraction/{{genome}}_final_coords.txt"
    params:
        upstream = UPSTREAM,
        downstream = DOWNSTREAM
    log:
        f"logs/{SAMPLE}/extract_coords/{{genome}}.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        if [ -s {input.orf_ids} ]; then
            # Extract full headers for matching ORFs
            seqkit grep -f {input.orf_ids} {input.proteins} | seqkit seq -n > {output.coordinates} 2>> {log}

            # Parse coordinates from prodigal headers (format: >contig_1 # start # end # strand # info)
            awk 'BEGIN{{OFS="\\t"; FS="#"}} {{
                up=$2; down=$3;
                gsub(/^[ \\t]+|[ \\t]+$/, "", up);
                gsub(/^[ \\t]+|[ \\t]+$/, "", down);
                if (up>down) print down, up; else print up, down
            }}' {output.coordinates} > {output.formatted} 2>> {log}

            # Add flanking regions
            awk -v upint={params.upstream} -v downint={params.downstream} \
                'BEGIN{{OFS="\\t"}} {{
                    up=$1-upint; down=$2+downint;
                    if (up<0) up=0;
                    print up, down
                }}' {output.formatted} > {output.final} 2>> {log}
        else
            touch {output.coordinates} {output.formatted} {output.final}
            echo "No ORF IDs found, creating empty outputs" >> {log}
        fi
        """


rule create_bed_file:
    """
    Create BED file for contig extraction.
    """
    input:
        blast = f"results/{SAMPLE}/blast_results/{{genome}}_blast.txt",
        final_coords = f"results/{SAMPLE}/orf_extraction/{{genome}}_final_coords.txt"
    output:
        contigs_list = f"results/{SAMPLE}/orf_extraction/{{genome}}_contigs_list.txt",
        bed = f"results/{SAMPLE}/bed_files/{{genome}}.bed"
    log:
        f"logs/{SAMPLE}/bed_file/{{genome}}.log"
    shell:
        """
        if [ -s {input.blast} ]; then
            # Extract contig names from subject IDs (format: contig_orf -> contig)
            cut -f2 {input.blast} | sed 's/_[0-9]*$//' > {output.contigs_list} 2>> {log}

            # Create BED file: contig, start, end
            paste {output.contigs_list} {input.final_coords} > {output.bed} 2>> {log}
        else
            touch {output.contigs_list} {output.bed}
            echo "No BLAST hits, creating empty outputs" >> {log}
        fi
        """


rule extract_contig_sequences:
    """
    Extract contig sequences using seqkit subseq.
    """
    input:
        bed = f"results/{SAMPLE}/bed_files/{{genome}}.bed",
        genome = f"results/{SAMPLE}/genomes/{{genome}}.fna"
    output:
        contigs = f"results/{SAMPLE}/contigs/{{genome}}_contigs.fna"
    log:
        f"logs/{SAMPLE}/extract_contigs/{{genome}}.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        if [ -s {input.bed} ]; then
            seqkit subseq --bed {input.bed} {input.genome} > {output.contigs} 2>> {log}
        else
            touch {output.contigs}
            echo "Empty BED file, creating empty output" >> {log}
        fi
        """


rule concatenate_contigs:
    """
    Concatenate all extracted contigs into a single file.
    """
    input:
        contigs = get_all_contigs
    output:
        concatenated = f"results/{SAMPLE}/all_contigs.fna",
        filtered = f"results/{SAMPLE}/all_contigs_filtered.fna"
    log:
        f"logs/{SAMPLE}/concatenate_contigs.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        # Concatenate all contig files
        cat {input.contigs} > {output.concatenated} 2>> {log}

        # Filter to keep only sequences >= 50bp
        seqkit seq -m 50 {output.concatenated} > {output.filtered} 2>> {log}

        echo "Concatenated $(grep -c '>' {output.concatenated} || echo 0) contigs" >> {log}
        echo "Filtered to $(grep -c '>' {output.filtered} || echo 0) contigs (>=50bp)" >> {log}
        """
