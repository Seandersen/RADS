# Installation Guide

This guide covers installing RADS on local machines and HPC clusters.

## Prerequisites

- **OS**: macOS (Intel/Apple Silicon), Linux (x86_64)
- **Internet connection**: Required for downloading genomes from NCBI
- **Disk space**: ~10 GB for pipeline dependencies, plus space for genomes

---

## Method 1: Pixi (Recommended)

[Pixi](https://pixi.sh) manages all dependencies in a single step and works on most systems, including many HPC clusters.

### Step 1: Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

After installation, restart your terminal or run:
```bash
source ~/.bashrc   # or ~/.zshrc on macOS
```

### Step 2: Clone the Repository

```bash
git clone https://github.com/Seandersen/RADS.git
cd RADS
git checkout snakemake-pipeline
```

### Step 3: Install Dependencies

```bash
pixi install
```

This installs all required packages including Snakemake, Diamond, Prodigal, SeqKit, NCBI datasets CLI, and the dashboard dependencies (Shiny, Plotly, Polars).

### Step 4: Verify

```bash
pixi run snakemake --version
pixi run dry-run
```

---

## Method 2: Conda / Mamba

Use this method if Pixi is not available on your system (common on some HPC clusters).
Mamba is strongly recommended over Conda for faster dependency resolution.

### Step 1: Install Mamba

**Option A — Miniforge (recommended, includes mamba):**

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

**Option B — Add mamba to existing Conda:**

```bash
conda install -n base -c conda-forge mamba
```

### Step 2: Clone and Create Environment

```bash
git clone https://github.com/Seandersen/RADS.git
cd RADS
git checkout snakemake-pipeline

mamba env create -f environment.yaml
conda activate rads
```

### Step 3: Install DefenseFinder

DefenseFinder is not in the Conda environment file and must be installed separately:

```bash
conda activate rads
pip install mdmparis-defense-finder
defense-finder update
```

### Step 4: Run the Pipeline

Run **without** `--use-conda` so Snakemake uses your active `rads` environment directly:

```bash
conda activate rads
snakemake --cores 8
```

> **Why not `--use-conda`?** Snakemake's per-rule conda environments may not have access to the DefenseFinder models installed above. Running against the `rads` environment directly is simpler and more reliable.

### HPC note

On shared HPC clusters, configure Conda to copy files rather than symlink (required on network filesystems like NFS or Lustre):

```bash
conda config --set always_copy true
conda config --set channel_priority strict
```

---

## Optional: InterProScan

InterProScan provides detailed domain annotations but requires a separate manual installation (~15 GB).

### Download and Setup

```bash
mkdir my_interproscan && cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz
tar -pxvzf interproscan-5.69-101.0-64-bit.tar.gz
cd interproscan-5.69-101.0
python3 setup.py -f interproscan.properties
```

### Configure RADS

```yaml
# config/config.yaml
interproscan:
  enabled: true
  path: "/absolute/path/to/interproscan-5.69-101.0/interproscan.sh"
  applications: "Pfam"
```

**Recommended: Pfam only.** Some InterProScan analyses (MobiDBLite, Panther, ProSiteProfiles) have binary dependencies that frequently fail on HPC systems.

Safe analyses to add: `Pfam`, `CDD`, `TIGRFAM`

**Platform notes:**
- Linux x86_64: Full InterProScan support
- macOS Intel: Works but slower
- macOS Apple Silicon (M1/M2/M3): Requires Rosetta 2 or Docker
- HPC: Use Pfam-only to avoid binary compatibility issues

---

## Verifying Your Installation

```bash
# Dry run — checks config and previews steps without executing
pixi run dry-run
# or (conda):
snakemake -n
```

---

## Troubleshooting Installation

### Conda solver too slow

```bash
# Install libmamba solver
conda install -n base conda-libmamba-solver
conda config --set solver libmamba

# Or use mamba directly
mamba env create -f environment.yaml
```

### DefenseFinder module error

If you see `ModuleNotFoundError: No module named 'macsypy'`:

```bash
pip install macsyfinder mdmparis-defense-finder --force-reinstall
```

### Pixi not available on compute nodes

If your HPC does not export the login-node `PATH` to batch jobs, add Pixi explicitly:

```bash
export PATH="$HOME/.pixi/bin:$PATH"
```

Or add this to your job script header.

---

## Next Steps

- [Configuration](Configuration.md) — set up your analysis
- [Pipeline Overview](Pipeline-Overview.md) — understand the steps
- [Advanced Usage](Advanced-Usage.md) — HPC and customization
