# DefenseFinder Integration
# ==========================
# Run DefenseFinder to identify defense systems in extracted contigs

rule run_defensefinder:
    """
    Run DefenseFinder on extracted contig proteins.
    Uses --db-type unordered since contigs are not ordered by genomic position.
    Uses pixi defensefinder environment to avoid model version conflicts.
    Handles DefenseFinder post-processing bug by copying results from /tmp/defense-finder.
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
        coverage = config.get("defensefinder", {}).get("coverage", 0.4)
    threads:
        config.get("defensefinder", {}).get("workers", 4)
    log:
        f"logs/{SAMPLE}/defensefinder.log"
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

        # Function to copy results from /tmp/defense-finder (handles post-processing bug)
        copy_from_tmp() {{
            if [ -f /tmp/defense-finder/all_systems.tsv ]; then
                echo "Copying results from /tmp/defense-finder..." >> {log}
                cp /tmp/defense-finder/all_systems.tsv {output.systems}
                # Create genes file from all_systems.tsv
                grep -v "^#" /tmp/defense-finder/all_systems.tsv > {output.genes} 2>/dev/null || touch {output.genes}
                # Create hmmer file with reformatted columns
                echo -e "hit_id\\treplicon\\tgene_name\\ti_eval\\tscore\\tprofile_cov\\tseq_cov\\tbegin_match\\tend_match" > {output.hmmer}
                grep -v "^#" /tmp/defense-finder/all_systems.tsv | awk -F'\\t' 'NR>1 {{print $2"\\t"$1"\\t"$3"\\t"$11"\\t"$12"\\t"$13"\\t"$14"\\t"$15"\\t"$16}}' >> {output.hmmer} 2>/dev/null || true
                return 0
            fi
            return 1
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

        # Clean up any previous /tmp/defense-finder results
        rm -rf /tmp/defense-finder 2>/dev/null || true

        # Run defense-finder using pixi environment to avoid model version conflicts
        # The pixi defensefinder environment has compatible versions of defense-finder and models
        if pixi run -e defensefinder defense-finder run \
            --db-type {params.db_type} \
            --workers {threads} \
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
            # defense-finder command failed - check if results are in /tmp/defense-finder
            # This handles a bug where post-processing fails but results exist
            echo "DefenseFinder command returned error, checking for results in /tmp/defense-finder..." >> {log}
            if copy_from_tmp; then
                echo "Successfully recovered results from /tmp/defense-finder" >> {log}
            else
                create_empty_outputs
                echo "DefenseFinder execution failed and no results found, created empty output files" >> {log}
            fi
        fi
        """
