# RADS Wiki - Recombinase Associated Defense Search

Welcome to the RADS documentation! This wiki provides comprehensive documentation for the Snakemake implementation of RADS.

## What is RADS?

RADS (Recombinase Associated Defense Search) is a bioinformatics pipeline for extracting and analyzing genomic loci surrounding a query protein sequence. Originally developed for identifying defense systems near the EF_B0058 serine recombinase, RADS can be used to study any genomic neighborhood around proteins of interest.

## Key Features

- **Automated Genome Processing** - Download genomes from NCBI or use local files
- **Parallel Analysis** - Process hundreds of genomes simultaneously
- **Defense System Detection** - Optional DefenseFinder integration
- **Domain Annotation** - Optional InterProScan integration
- **Interactive Dashboard** - Visualize and explore results
- **Co-transcription Analysis** - Identify nearby co-transcribed genes
- **Reproducible Environments** - Managed with Pixi/Conda

## Quick Links

| Page | Description |
|------|-------------|
| [[Installation]] | Setup instructions for all platforms |
| [[Configuration]] | Configuration file reference |
| [[Pipeline-Overview]] | Detailed pipeline steps |
| [[Dashboard-Guide]] | Using the results dashboard |
| [[Troubleshooting]] | Common issues and solutions |
| [[Advanced-Usage]] | Advanced features and customization |

## Quick Start

```bash
# 1. Install Pixi
curl -fsSL https://pixi.sh/install.sh | bash

# 2. Clone repository
git clone https://github.com/Seandersen/RADS.git
cd RADS
git checkout snakemake-pipeline

# 3. Install dependencies
pixi install

# 4. Configure (edit config/config.yaml)
# 5. Run pipeline
pixi run snakemake --cores 8

# 6. View results
pixi run dashboard
```

## Pipeline Output

The pipeline generates:

| Output | Description |
|--------|-------------|
| `master_blast.txt` | Combined BLAST results from all genomes |
| `all_contigs.fna` | Extracted flanking sequences |
| `all_contigs.faa` | Predicted ORFs from flanking regions |
| `cotranscribed_sequences.faa` | Downstream co-transcribed proteins |
| `interproscan_results.tsv` | Domain annotations (optional) |
| `defense_finder_systems.tsv` | Defense systems (optional) |
| `pipeline_metrics.json` | Summary statistics |

## Support

- **Issues**: [GitHub Issues](https://github.com/Seandersen/RADS/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Seandersen/RADS/discussions)

## Citation

If you use RADS in your research, please cite:

> Andersen SE, Kirsch JM, Hesselberth JR, Duerkop BA. RADS: Recombinase Associated Defense Search. [Publication details pending]
