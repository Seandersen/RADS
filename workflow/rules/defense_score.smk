# Defense Score Calculation
# Scores co-transcribed genes for defense island proximity


rule calculate_defense_scores:
    """Score co-transcribed genes based on proximity to defense systems."""
    input:
        cotranscribed=f"results/{SAMPLE}/cotranscription/cotranscribed_details.txt",
        defense_genes=f"results/{SAMPLE}/defensefinder/defense_finder_genes.tsv",
        contigs_faa=f"results/{SAMPLE}/contig_orfs/all_contigs.faa",
        interproscan=f"results/{SAMPLE}/interproscan_results.tsv",
    output:
        scores=f"results/{SAMPLE}/defense_scores.tsv",
    log:
        f"logs/{SAMPLE}/defense_scores.log",
    shell:
        """
        python defense_score.py \
            --cotranscribed {input.cotranscribed} \
            --defense-genes {input.defense_genes} \
            --contigs-faa {input.contigs_faa} \
            --interproscan {input.interproscan} \
            --output {output.scores} \
            2>&1 | tee {log}
        """
