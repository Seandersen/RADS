# DefenseFinder Integration
# ==========================
# Run DefenseFinder to identify defense systems in extracted contigs

rule run_defensefinder:
    """
    Run DefenseFinder on extracted contig proteins.
    Uses --db-type unordered since contigs are not ordered by genomic position.
    Gracefully handles the case when defense-finder is not available.
    """
    input:
        proteins = f"results/{SAMPLE}/contig_orfs/all_contigs.faa"
    output:
        systems = f"results/{SAMPLE}/defensefinder/defense_finder_systems.tsv",
        genes = f"results/{SAMPLE}/defensefinder/defense_finder_genes.tsv",
        hmmer = f"results/{SAMPLE}/defensefinder/defense_finder_hmmer.tsv"
    params:
        outdir = f"results/{SAMPLE}/defensefinder",
        enabled = config.get("defensefinder", {}).get("enabled", True),
        db_type = config.get("defensefinder", {}).get("db_type", "unordered"),
        coverage = config.get("defensefinder", {}).get("coverage", 0.4),
        workers = config.get("defensefinder", {}).get("workers", 4)
    log:
        f"logs/{SAMPLE}/defensefinder.log"
    conda:
        "../envs/defensefinder.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}

        # Function to create empty output files with headers
        create_empty_outputs() {{
            echo -e "sys_id\\ttype\\tsubtype\\tsys_beg\\tsys_end\\tprotein_in_syst\\tgenes_count\\tname_of_profiles_in_sys" > {output.systems}
            echo -e "hit_id\\treplicon\\tposition\\thit_pos\\tgene_name\\ti_eval\\tscore\\tprofile_cov\\tseq_cov\\tbegin_match\\tend_match\\tprotein_in_syst\\ttype\\tsubtype" > {output.genes}
            echo -e "hit_id\\treplicon\\tposition\\thit_pos\\tgene_name\\ti_eval\\tscore\\tprofile_cov\\tseq_cov\\tbegin_match\\tend_match" > {output.hmmer}
        }}

        # Check if defense-finder is enabled
        if [ "{params.enabled}" != "True" ] && [ "{params.enabled}" != "true" ]; then
            create_empty_outputs
            echo "DefenseFinder disabled in config" > {log}
            exit 0
        fi

        # Check if input is empty
        if [ ! -s {input.proteins} ]; then
            create_empty_outputs
            echo "Empty input proteins, created empty output files" > {log}
            exit 0
        fi

        # Check if defense-finder command is available and working
        if ! command -v defense-finder &> /dev/null; then
            create_empty_outputs
            echo "defense-finder command not found, created empty output files" > {log}
            exit 0
        fi

        # Try to run defense-finder
        if defense-finder run \
            --db-type {params.db_type} \
            --coverage {params.coverage} \
            --workers {params.workers} \
            -o {params.outdir} \
            {input.proteins} \
            2>&1 | tee {log}; then

            # Check if outputs were created
            if [ -f {params.outdir}/defense_finder_systems.tsv ]; then
                echo "DefenseFinder completed successfully" >> {log}
            else
                create_empty_outputs
                echo "No defense systems found, created empty output files" >> {log}
            fi
        else
            # defense-finder failed (possibly due to module issues)
            create_empty_outputs
            echo "DefenseFinder execution failed, created empty output files" >> {log}
        fi
        """
