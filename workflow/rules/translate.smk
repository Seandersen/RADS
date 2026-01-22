# Translate genomes using Prodigal
# =================================

rule translate_genome:
    """
    Translate a single genome using Prodigal.
    Corresponds to RADS.sh translate() function.
    """
    input:
        genome = f"results/{SAMPLE}/genomes/{{genome}}.fna",
        manifest = f"results/{SAMPLE}/genome_manifest.txt"
    output:
        proteins = f"results/{SAMPLE}/translated/{{genome}}.faa",
        genes = f"results/{SAMPLE}/translated/{{genome}}_prodigal.txt"
    params:
        mode = config["prodigal"]["mode"]
    log:
        f"logs/{SAMPLE}/translate/{{genome}}.log"
    benchmark:
        f"logs/{SAMPLE}/benchmarks/translate_{{genome}}.txt"
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal \
            -p {params.mode} \
            -i {input.genome} \
            -o {output.genes} \
            -a {output.proteins} \
            2>&1 | tee {log}
        """


def get_all_translated(wildcards):
    """Get all translated protein files."""
    genomes = get_all_genomes(wildcards)
    return expand(
        f"results/{SAMPLE}/translated/{{genome}}.faa",
        genome=genomes
    )


rule translate_all:
    """Aggregate rule to ensure all genomes are translated."""
    input:
        get_all_translated
    output:
        touch(f"results/{SAMPLE}/.translation_complete")
