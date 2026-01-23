# Installation Guide

This guide covers installing RADS on different platforms.

## Prerequisites

- **Operating System**: macOS (Intel/Apple Silicon), Linux (x86_64)
- **Internet connection**: Required for downloading genomes from NCBI
- **Disk space**: ~10GB for pipeline dependencies, plus space for genomes

## Method 1: Pixi (Recommended)

[Pixi](https://pixi.sh) is a fast, cross-platform package manager that handles all dependencies automatically.

### Step 1: Install Pixi

**macOS / Linux:**
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

**Windows (PowerShell):**
```powershell
iwr -useb https://pixi.sh/install.ps1 | iex
```

After installation, restart your terminal or run:
```bash
source ~/.bashrc  # or ~/.zshrc on macOS
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

This installs all required packages:
- Snakemake (workflow management)
- Diamond (BLAST searches)
- Prodigal (ORF prediction)
- SeqKit (sequence manipulation)
- NCBI datasets CLI (genome downloads)
- DefenseFinder (defense system detection)
- Dashboard dependencies (Shiny, Plotly, Polars)

### Step 4: Verify Installation

```bash
# Check Snakemake
pixi run snakemake --version

# Check Diamond
pixi run diamond --version

# Run dry-run to verify workflow
pixi run dry-run
```

## Method 2: Mamba (Recommended) or Conda

Mamba is strongly recommended over conda for faster dependency resolution. Conda's default solver can be extremely slow with complex bioinformatics environments.

### Step 1: Install Mamba

**Option A: Install Miniforge (includes mamba)**

Download [Miniforge](https://github.com/conda-forge/miniforge) - this is the recommended approach:
```bash
# macOS/Linux
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

**Option B: Add mamba to existing conda**
```bash
conda install -n base -c conda-forge mamba
```

### Step 2: Create Environment

```bash
git clone https://github.com/Seandersen/RADS.git
cd RADS
git checkout snakemake-pipeline

# Using mamba (recommended - much faster)
mamba env create -f environment.yaml
conda activate rads

# Or using conda (slower)
conda env create -f environment.yaml
conda activate rads
```

### Step 3: Run Pipeline

```bash
snakemake --cores 8 --use-conda
```

**Note**: When using `--use-conda`, Snakemake creates isolated environments for each rule. Mamba will be used automatically if available, significantly speeding up environment creation.

## Method 3: HPC / Supercomputer Installation

When running on shared computing clusters or supercomputers, you may encounter issues with network filesystems. Follow these additional steps:

### Configure Conda for Network Filesystems

Network filesystems (NFS, CIFS, Lustre) often don't support symbolic links. Configure conda to copy files instead:

```bash
conda config --set always_copy true
conda config --set channel_priority strict
```

### Use Local Storage for Conda Environments

Snakemake creates conda environments for each rule. Store these on local disk (not network storage) to avoid symlink errors:

```bash
# Create a local directory for conda environments
mkdir -p /tmp/$USER/snakemake_conda

# Run with --conda-prefix pointing to local storage
snakemake --cores 20 --use-conda --conda-prefix /tmp/$USER/snakemake_conda
```

Common local storage paths on HPC systems:
- `/tmp/$USER/`
- `/scratch/$USER/`
- `/local/$USER/`
- `$HOME/.snakemake/conda` (if home is on local disk)

### Example HPC Run Command

```bash
snakemake --cores 20 --use-conda --conda-prefix /tmp/$USER/snakemake_conda
```

### Troubleshooting HPC Issues

**Error: `[Errno 95] Operation not supported: 'cacert.pem'`**

This indicates symlink issues on network filesystem. Solution:
```bash
conda config --set always_copy true
rm -rf .snakemake/conda/*  # Clear failed environments
snakemake --cores 20 --use-conda --conda-prefix /tmp/$USER/snakemake_conda
```

## Optional: InterProScan Setup

InterProScan provides domain annotations but requires separate installation due to its size (~15GB).

### Download InterProScan

```bash
mkdir my_interproscan && cd my_interproscan

# Download (choose appropriate version)
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz

# Extract
tar -pxvzf interproscan-5.69-101.0-64-bit.tar.gz
cd interproscan-5.69-101.0

# Setup
python3 setup.py -f interproscan.properties
```

### Configure RADS

Update `config/config.yaml`:

```yaml
interproscan:
  enabled: true
  path: "/absolute/path/to/interproscan-5.69-101.0/interproscan.sh"
  applications: "Pfam"  # Recommended: Pfam only for reliability
```

### Recommended: Use Pfam Only

Some InterProScan analyses (MobiDBLite, Panther, ProSiteProfiles) have binary dependencies that frequently fail on HPC systems. We recommend running only Pfam analysis:

```yaml
interproscan:
  enabled: true
  path: "/path/to/interproscan.sh"
  applications: "Pfam"  # Most reliable
```

If you need more analyses, these are generally safe to add:
```yaml
  applications: "Pfam,CDD,TIGRFAM"
```

**Avoid these analyses** (often have binary/dependency issues):
- MobiDBLite
- Panther
- ProSiteProfiles

### Platform Notes

- **Linux x86_64**: Full InterProScan support
- **macOS Intel**: Works but slower
- **macOS Apple Silicon (M1/M2/M3)**: InterProScan requires Rosetta 2 or Docker
- **HPC/Clusters**: Use Pfam-only mode to avoid binary compatibility issues
- **Alternative**: Disable InterProScan and annotate sequences manually

## Optional: DefenseFinder Database

DefenseFinder requires its database to be downloaded on first run:

```bash
pixi run defense-finder update
```

This downloads the MacSyFinder models for defense system detection.

## Verifying Your Installation

Run a test analysis with the included test data:

```bash
# Create a small test accession file
echo "NC_000913.3" > test_accessions.txt

# Update config for test
# Edit config/config.yaml:
#   sample_name: "installation_test"
#   download:
#     enabled: true
#     accession_file: "test_accessions.txt"

# Run dry-run
pixi run snakemake -n

# If dry-run succeeds, run full test
pixi run snakemake --cores 4
```

## Troubleshooting Installation

### Pixi Installation Fails

```bash
# Clear Pixi cache
rm -rf ~/.pixi

# Reinstall
curl -fsSL https://pixi.sh/install.sh | bash
```

### Conda Solver Too Slow

Use Mamba instead (strongly recommended):
```bash
# Option 1: Install mamba directly
conda install -n base -c conda-forge mamba

# Then use mamba instead of conda
mamba env create -f environment.yaml

# Option 2: Use libmamba solver with conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

### DefenseFinder Module Error

If you see "ModuleNotFoundError: No module named 'macsypy'":
```bash
# Reinstall DefenseFinder dependencies
pixi run pip install macsyfinder defense-finder --force-reinstall
```

### Permission Denied Errors

```bash
# Fix script permissions
chmod +x RADS.sh
chmod +x workflow/scripts/*.py
```

## Running the Dashboard Without Pixi

If pixi is not available (e.g., on HPC systems), you can run the dashboard using conda/pip:

### Install Dashboard Dependencies

```bash
# Option 1: Create a new conda environment
conda create -n rads-dashboard python=3.10 -c conda-forge -y
conda activate rads-dashboard
pip install shiny polars plotly pandas pyarrow

# Option 2: Install in existing environment
pip install shiny polars plotly pandas pyarrow
```

### Run the Dashboard

```bash
cd /path/to/RADS
shiny run dashboard/app.py --port 8000
```

### Remote Access (SSH Tunnel)

When running on a remote server/HPC, create an SSH tunnel to access the dashboard:

**On your local machine**, open a new terminal:
```bash
ssh -L 8000:localhost:8000 username@remote-server
```

Then open `http://localhost:8000` in your local browser.

### Keep Dashboard Running (screen/tmux)

To keep the dashboard running after disconnecting:

```bash
# Start a screen session
screen -S dashboard
shiny run dashboard/app.py --port 8000

# Detach: press Ctrl+A, then D
# Reattach later: screen -r dashboard
```

See [[Dashboard-Guide]] for detailed dashboard usage.

## Next Steps

- [[Configuration]] - Learn how to configure the pipeline
- [[Pipeline-Overview]] - Understand the analysis steps
