# Build Diamond databases from translated proteins
# =================================================

rule build_database:
    """
    Build a Diamond database from translated proteins.
    Corresponds to RADS.sh makedbs() function.
    Handles empty input files gracefully.
    """
    input:
        proteins = f"results/{SAMPLE}/translated/{{genome}}.faa"
    output:
        db = f"results/{SAMPLE}/diamond_dbs/{{genome}}.dmnd"
    params:
        threads = config["diamond"]["threads"]
    log:
        f"logs/{SAMPLE}/diamond_db/{{genome}}.log"
    benchmark:
        f"logs/{SAMPLE}/benchmarks/diamond_db_{{genome}}.txt"
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        if [ -s {input.proteins} ]; then
            diamond makedb \
                --in {input.proteins} \
                --db {output.db} \
                --threads {params.threads} \
                2>&1 | tee {log}
        else
            echo "Empty input file, creating empty database marker" | tee {log}
            touch {output.db}
        fi
        """


def get_all_databases(wildcards):
    """Get all Diamond database files."""
    genomes = get_all_genomes(wildcards)
    return expand(
        f"results/{SAMPLE}/diamond_dbs/{{genome}}.dmnd",
        genome=genomes
    )


rule build_all_databases:
    """Aggregate rule to ensure all databases are built."""
    input:
        get_all_databases
    output:
        touch(f"results/{SAMPLE}/.databases_complete")
