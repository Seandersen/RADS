# RADS - Recombinase Associated Defense Search

A bioinformatics pipeline for extracting and analyzing genomic loci surrounding a query protein sequence.

RADS was developed by Shelby E Andersen in collaboration with Joshua M Kirsch, Jay R Hesselberth, and Breck A Duerkop as a tool for identifying defense systems near recombinase genes. RADS can be used to extract any sized genomic neighborhood around proteins of interest.

## Features

- **Automated genome download** from NCBI using accession lists
- **Parallel processing** of hundreds of genomes
- **Defense system detection** via DefenseFinder integration
- **Domain annotation** via InterProScan (optional)
- **Co-transcription analysis** to identify nearby genes
- **Interactive dashboard** for exploring results
- **Reproducible environments** via Pixi/Conda

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
| [Advanced Usage](docs/Advanced-Usage.md) | Customization and HPC |

## Quick Start

### 1. Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 2. Clone and Setup

```bash
git clone https://github.com/Seandersen/RADS.git
cd RADS
git checkout snakemake-pipeline
pixi install
```

### Alternative: Using Mamba/Conda

If you prefer mamba (recommended over conda for faster dependency resolution):

```bash
# Install mamba if not already installed
conda install -n base -c conda-forge mamba

# Create environment using mamba
mamba env create -f environment.yaml
conda activate rads

# Run pipeline with conda environments
snakemake --cores 8 --use-conda
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

interproscan:
  enabled: false  # Requires separate installation

defensefinder:
  enabled: true
```

### 4. Create Accession File

Create a file with NCBI nucleotide accessions (one per line):

```
NC_000913.3
NZ_CP088158.1
CP016471.1
```

### 5. Run Pipeline

```bash
# Full pipeline
pixi run snakemake --cores 8

# Dry run (preview)
pixi run snakemake -n
```

### 6. Explore Results

```bash
pixi run dashboard
# Open http://localhost:8080
```

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
| BLAST | Search query vs all genomes |
| Extract | Get flanking regions around hits |
| Annotate | Domain and defense system analysis |

## Output Structure

```
results/{sample}/
├── blast_results/master_blast.txt     # Combined BLAST results
├── all_contigs.fna                    # Extracted flanking regions
├── contig_orfs/all_contigs.faa        # Predicted ORFs
├── cotranscription/                   # Co-transcribed genes
├── defensefinder/                     # Defense systems
├── interproscan_results.tsv           # Domain annotations
└── metrics/pipeline_metrics.json      # Summary statistics
```

## Dashboard

The interactive dashboard provides:

- Summary statistics and metrics
- BLAST result visualization (scatter plots, histograms)
- Contig and ORF analysis
- Co-transcription pair tables
- Domain annotation charts
- Defense system breakdown

![Dashboard Screenshot](docs/images/dashboard-screenshot.png)

## Useful Commands

```bash
# Visualize workflow
pixi run snakemake --dag | dot -Tsvg > dag.svg

# Run specific step
pixi run snakemake results/{sample}/blast_results/master_blast.txt --cores 4

# Resume after failure
pixi run snakemake --cores 8 --rerun-incomplete

# Clean results
pixi run clean
```

## Requirements

- **OS**: macOS (Intel/Apple Silicon), Linux
- **Disk**: ~10GB for dependencies + genome data
- **Memory**: 8GB+ recommended

### Optional

- **InterProScan**: For domain annotation (~15GB)
- **DefenseFinder database**: For defense system detection

## Citation

If you use RADS in your research, please cite:

> Andersen SE, Kirsch JM, Hesselberth JR, Duerkop BA. RADS: Recombinase Associated Defense Search. [Publication details pending]

## Support

- **Documentation**: [docs/](docs/)
- **Issues**: [GitHub Issues](https://github.com/Seandersen/RADS/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Seandersen/RADS/discussions)

## License

MIT License - See [LICENSE](LICENSE) for details.

---

## Legacy Bash Pipeline

The original bash script (`RADS.sh`) is available on the `main` branch. See the main branch README for legacy documentation.
