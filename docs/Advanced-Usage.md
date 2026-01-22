# Advanced Usage

This guide covers advanced features and customization options for RADS.

## Running Specific Pipeline Steps

### Target Specific Outputs

```bash
# Run only up to BLAST
pixi run snakemake results/{sample}/blast_results/master_blast.txt --cores 8

# Run only DefenseFinder
pixi run snakemake results/{sample}/defensefinder/defense_finder_systems.tsv --cores 4

# Run multiple specific targets
pixi run snakemake \
    results/{sample}/blast_results/master_blast.txt \
    results/{sample}/cotranscription/cotranscribed_sequences.faa \
    --cores 8
```

### Skip Specific Steps

```yaml
# Disable optional components in config/config.yaml
interproscan:
  enabled: false

defensefinder:
  enabled: false
```

### Force Re-run Specific Rules

```bash
# Re-run a specific rule and all downstream
pixi run snakemake --forcerun blast_search --cores 8

# Re-run everything
pixi run snakemake --forceall --cores 8
```

## Working with Multiple Samples

### Batch Processing

Create a script to run multiple samples:

```bash
#!/bin/bash
# run_all_samples.sh

samples=("sample1" "sample2" "sample3")

for sample in "${samples[@]}"; do
    echo "Processing $sample..."
    pixi run snakemake \
        --configfile config/${sample}_config.yaml \
        --cores 8 \
        || echo "Failed: $sample"
done
```

### Parallel Sample Processing

```bash
# Using GNU parallel
parallel -j 2 'pixi run snakemake \
    --configfile config/{}_config.yaml --cores 4' ::: sample1 sample2 sample3
```

### Comparing Results Across Samples

```python
# compare_samples.py
import polars as pl
from pathlib import Path

samples = ["sample1", "sample2", "sample3"]
all_results = []

for sample in samples:
    blast_file = Path(f"results/{sample}/blast_results/master_blast.txt")
    if blast_file.exists():
        df = pl.read_csv(blast_file, separator="\t")
        df = df.with_columns(pl.lit(sample).alias("sample"))
        all_results.append(df)

combined = pl.concat(all_results)
print(combined.group_by("sample").len())
```

## Custom Query Sequences

### Multiple Query Proteins

The pipeline supports multi-FASTA queries:

```fasta
>protein1
MKTQPIKVN...
>protein2
MLIVGQRPK...
>protein3
MSDFTYHNK...
```

### Query from GenBank

```bash
# Download protein sequence from NCBI
pixi run efetch -db protein -id "WP_123456789.1" -format fasta > resources/my_query.fa
```

### Query Preprocessing

```bash
# Remove problematic characters
pixi run seqkit seq -w 0 resources/raw_query.fa > resources/clean_query.fa
```

## Customizing BLAST Parameters

### E-value Threshold

Edit `workflow/rules/blast_search.smk`:

```python
rule blast_search:
    shell:
        """
        diamond blastp \
            -d {input.db} \
            -q {params.query} \
            --evalue 1e-10 \  # Add e-value threshold
            ...
        """
```

### Additional Output Columns

```python
# Add more BLAST output columns
-outfmt "6 qseqid sseqid length nident pident evalue bitscore qcovhsp"
```

### Use BLASTp Instead of Diamond

For more sensitive searches:

```python
rule blast_search:
    conda:
        "../envs/blast.yaml"
    shell:
        """
        blastp \
            -db {input.db} \
            -query {params.query} \
            -outfmt "6 qseqid sseqid length nident pident evalue" \
            -num_threads {threads} \
            -out {output}
        """
```

## Customizing Flanking Region Extraction

### Asymmetric Flanking Regions

```yaml
# config/config.yaml
upstream_nt: 10000    # More upstream
downstream_nt: 3000   # Less downstream
```

### Variable Flanking by Hit

Modify `workflow/rules/extract_contigs.smk` to use different flanking sizes based on hit properties.

### Include Hit Sequence in Output

By default, flanking regions include the hit. To exclude:

```python
# In extract_contigs.smk, modify coordinate calculation:
start = max(0, hit_start - upstream_nt)  # Start before hit
end = min(genome_len, hit_end + downstream_nt)  # End after hit
```

## Custom Annotations

### Adding New Annotation Tools

Create a new rule file `workflow/rules/my_annotation.smk`:

```python
rule run_my_tool:
    input:
        proteins = f"results/{SAMPLE}/contig_orfs/all_contigs.faa"
    output:
        results = f"results/{SAMPLE}/my_annotation/results.tsv"
    log:
        f"logs/{SAMPLE}/my_annotation.log"
    shell:
        """
        my_tool -i {input.proteins} -o {output.results} 2>&1 | tee {log}
        """
```

Add to Snakefile:
```python
include: "workflow/rules/my_annotation.smk"

rule all:
    input:
        ...
        f"results/{SAMPLE}/my_annotation/results.tsv",
```

### Custom Domain Database

```python
rule custom_hmmscan:
    input:
        proteins = f"results/{SAMPLE}/contig_orfs/all_contigs.faa",
        db = "resources/my_custom.hmm"
    output:
        f"results/{SAMPLE}/custom_hmmscan/results.txt"
    shell:
        """
        hmmscan --tblout {output} {input.db} {input.proteins}
        """
```

## Cluster/HPC Execution

### SLURM Configuration

Create `profiles/slurm/config.yaml`:

```yaml
executor: slurm
jobs: 100
default-resources:
  slurm_partition: compute
  mem_mb: 4000
  runtime: 60
  cpus_per_task: 4
```

Run:
```bash
pixi run snakemake --profile profiles/slurm
```

### PBS/Torque

```yaml
executor: cluster-generic
cluster-generic-submit-cmd: qsub
default-resources:
  nodes: 1
  ppn: 4
  walltime: "01:00:00"
```

### Cloud Execution (AWS Batch)

```yaml
executor: aws-batch
aws-batch-queue: my-queue
aws-batch-job-role: arn:aws:iam::123456789:role/batch-role
```

## Workflow Visualization

### Generate DAG

```bash
# Full DAG
pixi run snakemake --dag | dot -Tsvg > dag.svg

# Simplified rulegraph
pixi run snakemake --rulegraph | dot -Tpng > rulegraph.png

# File graph
pixi run snakemake --filegraph | dot -Tpdf > filegraph.pdf
```

### Generate Report

```bash
pixi run snakemake --report report.html
```

## Performance Tuning

### Resource Allocation

Add resources to rules:

```python
rule translate_genome:
    resources:
        mem_mb=4000,
        runtime=30
    threads: 2
```

### Caching

Enable between-workflow caching:

```bash
pixi run snakemake --cache --cores 8
```

### Benchmarking

Add benchmarks to track performance:

```python
rule blast_search:
    benchmark:
        "benchmarks/{sample}/blast_search.txt"
```

## Integration with Other Tools

### Export to Galaxy

```bash
# Create Galaxy-compatible output structure
mkdir -p galaxy_export
cp results/{sample}/blast_results/master_blast.txt galaxy_export/
cp results/{sample}/contig_orfs/all_contigs.faa galaxy_export/
```

### Export to Anvi'o

```bash
# Convert to Anvi'o-compatible format
pixi run python scripts/export_anvio.py results/{sample}
```

### Integration with Nextflow

For Nextflow compatibility, outputs follow standard naming conventions that can be imported into Nextflow workflows.

## Extending the Dashboard

### Add Custom Tab

In `dashboard/app.py`:

```python
# Add to app_ui navset_card_tab:
ui.nav_panel(
    "My Analysis",
    ui.card(
        ui.card_header("Custom Visualization"),
        ui.output_ui("my_custom_viz"),
    ),
),

# Add to server function:
@render.ui
def my_custom_viz():
    # Load your data
    df = load_my_data(results_dir())
    if df is None:
        return ui.p("No data")

    # Create visualization
    fig = px.bar(df.to_pandas(), x="category", y="count")
    return ui.HTML(fig.to_html(include_plotlyjs="cdn", full_html=False))
```

### Add Custom Data Loader

In `dashboard/utils/data_loader.py`:

```python
def load_my_data(results_dir: str) -> Optional[pl.DataFrame]:
    """Load custom analysis results."""
    file_path = Path(results_dir) / "my_analysis" / "results.tsv"
    if not file_path.exists():
        return None

    return pl.read_csv(file_path, separator="\t")
```

## Development and Testing

### Run Tests

```bash
# Use test configuration
pixi run snakemake --configfile config/test.yaml --cores 4

# Dry run to check syntax
pixi run snakemake -n
```

### Lint Snakefiles

```bash
pixi run snakefmt workflow/rules/
pixi run snakefmt Snakefile
```

### Debug Mode

```bash
# Maximum verbosity
pixi run snakemake --cores 4 -p --verbose --debug-dag

# Print shell commands
pixi run snakemake --cores 4 --printshellcmds
```

## Next Steps

- [[Pipeline-Overview]] - Understanding the pipeline
- [[Troubleshooting]] - Solving common issues
- [Snakemake Documentation](https://snakemake.readthedocs.io/) - Official Snakemake docs
