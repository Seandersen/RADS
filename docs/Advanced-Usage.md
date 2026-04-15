# Advanced Usage

## Running Specific Pipeline Steps

### Target Specific Outputs

```bash
# Run only up to BLAST
pixi run snakemake results/{sample}/blast_results/master_blast.txt --cores 8

# Run only DefenseFinder
pixi run snakemake results/{sample}/defensefinder/defense_finder_systems.tsv --cores 4

# Run binomial analysis (requires interproscan and co-transcription to be complete)
pixi run snakemake results/{sample}/BinomialAnalysis.csv --cores 8
```

### Force Re-run Specific Rules

```bash
# Re-run a specific rule and all downstream steps
pixi run snakemake --forcerun blast_search --cores 8

# Re-run everything
pixi run snakemake --forceall --cores 8
```

---

## Generating the HTML Report

After the pipeline completes, generate a self-contained shareable report:

```bash
pixi run python workflow/scripts/generate_report.py \
    --results results/my_analysis \
    --output  results/my_analysis/report.html
```

Add `--include-locus-viewer` to embed interactive gene-arrow diagrams (larger file):

```bash
pixi run python workflow/scripts/generate_report.py \
    --results results/my_analysis \
    --output  results/my_analysis/report.html \
    --include-locus-viewer \
    --locus-max-contigs 200
```

The report opens in any browser with no server required. See [Dashboard Guide](Dashboard-Guide.md#static-html-report) for a full description of the report tabs.

---

## HPC / SLURM Execution

> **Experimental.** See [SLURM Usage](SLURM-Usage.md) for setup, configuration, and known limitations.

RADS includes a SLURM profile at `profiles/slurm/` that submits each rule as an independent batch job. Quick start:

```bash
pixi install        # one-time, on login node
pixi run slurm      # submit pipeline
```

Full documentation including cluster-specific settings, parallelization details, and binomial chunk configuration is in [SLURM-Usage.md](SLURM-Usage.md).

---

## Working with Multiple Samples

### Batch Processing

```bash
#!/bin/bash
for sample in sample1 sample2 sample3; do
    echo "Processing $sample..."
    pixi run snakemake \
        --configfile config/${sample}_config.yaml \
        --cores 8 \
        || echo "Failed: $sample"
done
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

---

## Custom Query Sequences

### Multiple Query Proteins

The pipeline supports multi-FASTA queries:

```fasta
>protein1
MKTQPIKVN...
>protein2
MLIVGQRPK...
```

### Download Query from NCBI

```bash
pixi run efetch -db protein -id "WP_123456789.1" -format fasta > resources/my_query.fa
```

---

## Customizing BLAST Parameters

### E-value Threshold

Edit `workflow/rules/blast_search.smk` to add an e-value cutoff:

```python
diamond blastp \
    -d {input.db} \
    --query {input.query} \
    --evalue 1e-10 \
    ...
```

### Additional Output Columns

```
--outfmt 6 qseqid sseqid length nident pident evalue bitscore qcovhsp
```

---

## Customizing Flanking Region Extraction

```yaml
# config/config.yaml
upstream_nt: 10000    # more upstream
downstream_nt: 3000   # less downstream
```

---

## Adding Custom Annotation Rules

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

Add to `Snakefile`:
```python
include: "workflow/rules/my_annotation.smk"

rule all:
    input:
        ...
        f"results/{SAMPLE}/my_annotation/results.tsv",
```

---

## Performance Tuning

### Runtime Estimates

Key factors affecting runtime:

- **InterProScan dominates at scale.** Using `applications: "Pfam"` instead of all databases reduces IPS time by 5–10×.
- **Binomial analysis adds cost.** It runs InterProScan on whole-genome proteomes of all genomes with BLAST hits. Providing a precomputed `whole_genome_interproscan` TSV skips this step entirely.
- **Core scaling has diminishing returns** above ~16–32 cores because InterProScan and DefenseFinder are largely sequential.
- **`max_genomes`** is the fastest way to reduce runtime for exploratory runs.

### Precomputed Whole-Genome InterProScan

If you have InterProScan results for your genomes from a previous run:

```yaml
binomial:
  whole_genome_interproscan: "/path/to/precomputed_ips.tsv"
```

---

## Workflow Visualization

```bash
# Full DAG
pixi run snakemake --dag | dot -Tsvg > dag.svg

# Simplified rule graph
pixi run snakemake --rulegraph | dot -Tpng > rulegraph.png
```

---

## Debug Mode

```bash
# Maximum verbosity
pixi run snakemake --cores 4 -p --verbose --debug-dag

# Print shell commands as they run
pixi run snakemake --cores 4 --printshellcmds
```

---

## Next Steps

- [Pipeline Overview](Pipeline-Overview.md) — understanding the analysis steps
- [SLURM Usage](SLURM-Usage.md) — HPC cluster submission (experimental)
- [Troubleshooting](Troubleshooting.md) — solving common issues
- [Snakemake Documentation](https://snakemake.readthedocs.io/) — official Snakemake docs
