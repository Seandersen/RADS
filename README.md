# RADS - Recombinase Associated Defense Search

A bioinformatics pipeline for extracting and analyzing genomic loci surrounding a query protein sequence.

RADS was developed by Shelby E Andersen in collaboration with Joshua M Kirsch, Jay R Hesselberth, and Breck A Duerkop as a tool for identifying defense systems near recombinase genes. RADS can be used to extract any sized genomic neighborhood around proteins of interest.

## Features

- **Automated genome download** from NCBI using accession lists
- **Parallel processing** of hundreds of genomes
- **Defense system detection** via DefenseFinder integration
- **Domain annotation** via InterProScan
- **Co-transcription analysis** to identify nearby genes
- **Interactive dashboard** for exploring results
- **Shareable HTML report** for distributing results without a server
- **Reproducible environments** via Pixi or Conda

## Documentation

Full documentation is available in the [docs/](docs/) folder:

| Guide | Description |
|-------|-------------|
| [Home](docs/Home.md) | Wiki overview and quick links |
| [Installation](docs/Installation.md) | Setup instructions for all platforms |
| [Configuration](docs/Configuration.md) | Configuration file reference |
| [Pipeline Overview](docs/Pipeline-Overview.md) | Detailed pipeline steps |
| [Dashboard Guide](docs/Dashboard-Guide.md) | Using the results explorer |
| [Troubleshooting](docs/Troubleshooting.md) | Common issues and solutions |
| [Advanced Usage](docs/Advanced-Usage.md) | Customization and advanced features |
| [SLURM Usage](docs/SLURM-Usage.md) | HPC cluster execution (experimental) |

## Quick Start

### 1. Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

> **Conda/Mamba users:** If Pixi is not available on your system, see [Installation](docs/Installation.md) for the Conda-based setup.

### 2. Clone and Setup

```bash
git clone https://github.com/Seandersen/RADS.git
cd RADS
git checkout snakemake-pipeline
pixi install
```

### 3. Configure

Edit `config/config.yaml`:

```yaml
sample_name: "my_analysis"

download:
  enabled: true
  accession_file: "my_accessions.txt"  # One accession per line

query_file: "resources/EFB0058.fa"

upstream_nt: 5000
downstream_nt: 5000

diamond:
  identity: 30
  threads: 8
```

### 4. Create Accession File

One NCBI nucleotide accession per line:

```
NC_000913.3
NZ_CP088158.1
CP016471.1
```

### 5. Run Pipeline

```bash
# Dry run (preview steps)
pixi run dry-run

# Full pipeline
pixi run snakemake --cores 8
```

### 6. Generate Shareable Report

```bash
pixi run python workflow/scripts/generate_report.py \
    --results results/my_analysis \
    --output  results/my_analysis/report.html
```

Open `report.html` in any browser — no server required. Add `--include-locus-viewer` to embed gene-arrow diagrams (larger file).

### 7. Explore Results Interactively

```bash
# Local machine
pixi run dashboard

# HPC (then access via browser — see Dashboard Guide)
pixi run dashboard-hpc
```

Open `http://localhost:8000` in your browser.

## Pipeline Overview

```
Genomes → Translate → BLAST → Extract Contigs → Annotate
                                      ↓
                              Co-transcription
                              DefenseFinder
                              InterProScan
```

| Step | Description |
|------|-------------|
| Download | Fetch genomes from NCBI |
| Translate | Predict ORFs with Prodigal |
| BLAST | Search query vs all genomes (Diamond) |
| Extract | Get flanking regions around hits |
| Annotate | Domain and defense system analysis |

## Output Structure

```
results/{sample}/
├── blast_results/master_blast.txt     # Combined BLAST results
├── all_contigs_filtered.fna           # Extracted flanking regions
├── contig_orfs/all_contigs.faa        # Predicted ORFs
├── cotranscription/                   # Co-transcribed genes
├── defensefinder/                     # Defense systems
├── interproscan_results.tsv           # Domain annotations
├── BinomialAnalysis.csv               # Enriched domains
├── defense_scores.tsv                 # Defense association scores
├── metrics/pipeline_metrics.json      # Summary statistics
└── report.html                        # Shareable HTML report
```

## Useful Commands

```bash
# Preview workflow (no execution)
pixi run dry-run

# Resume after failure
pixi run snakemake --cores 8 --rerun-incomplete

# Run specific step
pixi run snakemake results/{sample}/blast_results/master_blast.txt --cores 4

# Clean results
pixi run clean
```

## Requirements

- **OS**: macOS (Intel/Apple Silicon), Linux
- **Disk**: ~10 GB for dependencies + genome data
- **Memory**: 8 GB+ recommended

### Optional (separate installation)

- **InterProScan**: For domain annotation (~15 GB) — see [Installation](docs/Installation.md)

## Citation

If you use RADS in your research, please cite:

> Andersen SE, Kirsch JM, Hesselberth JR, Duerkop BA. RADS: Recombinase Associated Defense Search. [Publication details pending]

## Support

- **Documentation**: [docs/](docs/)
- **Issues**: [GitHub Issues](https://github.com/Seandersen/RADS/issues)

## License

MIT License - See [LICENSE](LICENSE) for details.

---

## Legacy Bash Pipeline

The original bash script (`RADS.sh`) is available on the `main` branch. See the main branch README for legacy documentation.
