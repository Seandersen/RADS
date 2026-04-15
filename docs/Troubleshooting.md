# Troubleshooting Guide

---

## Installation Issues

### Pixi Won't Install

```bash
# Alternative installation method
wget -qO- https://pixi.sh/install.sh | bash
```

### Conda Solver Too Slow

Use Mamba instead:
```bash
conda install -n base -c conda-forge mamba
mamba env create -f environment.yaml

# Or enable the libmamba solver for conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

### DefenseFinder Module Error

**Symptom:** `ModuleNotFoundError: No module named 'macsypy'`

```bash
# Reinstall
pip install macsyfinder mdmparis-defense-finder --force-reinstall

# Or disable if not needed
# config/config.yaml:
defensefinder:
  enabled: false
```

---

## Pipeline Execution Issues

### SyntaxError When Running Scripts

**Symptom:** `SyntaxError: Non-ASCII character` or `SyntaxError: invalid syntax` on f-strings

**Cause:** The system `python` command points to Python 2.

**Solution:** Always use `pixi run python` or `python3`:

```bash
# Wrong:
python workflow/scripts/generate_report.py ...

# Correct:
pixi run python workflow/scripts/generate_report.py ...
# or:
python3 workflow/scripts/generate_report.py ...
```

### IncompleteFilesException

**Symptom:**
```
IncompleteFilesException: The files below seem to be incomplete...
```

**Cause:** A previous run was interrupted before a rule finished writing its output. Snakemake marks the output as incomplete.

**Solution:**
```bash
# Clear the incomplete marker
rm -rf .snakemake/incomplete/

# Then re-run normally
pixi run snakemake --cores 8
# or add --rerun-incomplete to be explicit:
pixi run snakemake --cores 8 --rerun-incomplete
```

### "Nothing to be done"

**Symptom:** Snakemake says nothing to do but expected outputs are missing.

1. Check if outputs exist:
   ```bash
   ls -la results/{sample}/blast_results/master_blast.txt
   ```
2. Clear Snakemake metadata and retry:
   ```bash
   rm -rf .snakemake/metadata
   pixi run snakemake --cores 8
   ```
3. Force regeneration:
   ```bash
   pixi run snakemake --forceall --cores 8
   ```

### Genome Download Fails

**Valid accession formats:**
```
NC_000913.3
NZ_CP088158.1
CP016471.1
GCF_000001405.40
```

**Test a single accession:**
```bash
pixi run datasets summary genome accession NC_000913.3
pixi run datasets download genome accession NC_000913.3 --filename test.zip
```

**Use local genomes instead:**
```yaml
# config/config.yaml
download:
  enabled: false
genomes_path: "path/to/local/genomes"
```

### BLAST Returns No Hits

1. **Check query format** — must be amino acid FASTA:
   ```bash
   head resources/query.fa
   # Expected: >protein_name followed by amino acid sequence
   ```
2. **Lower identity threshold:**
   ```yaml
   diamond:
     identity: 20   # default is 30
   ```
3. **Verify databases were built:**
   ```bash
   ls -lh results/{sample}/diamond_dbs/*.dmnd
   ```

### InterProScan Fails

1. **Verify installation:**
   ```bash
   /path/to/interproscan.sh --version
   java -version   # needs Java 11+
   ```
2. **Use absolute path in config** (not relative)
3. **Switch to Pfam-only** if binary errors occur (MobiDB, Panther, ProSiteProfiles):
   ```yaml
   interproscan:
     applications: "Pfam"
   ```
4. **Disable if not needed:**
   ```yaml
   interproscan:
     enabled: false
   ```

**Rerun after fixing:**
```bash
rm -f results/{sample}/interproscan_results.tsv
pixi run snakemake results/{sample}/interproscan_results.tsv --cores 8
```

### Out of Memory

```bash
# Reduce parallelism
pixi run snakemake --cores 2

# Limit genome count in config
download:
  max_genomes: 50
```

### Disk Space Issues

```bash
# Remove downloaded genomes after staging (they're copied to results/{sample}/genomes/)
rm -rf results/{sample}/downloaded_genomes

# Clean results and logs entirely
pixi run clean

# Clean Snakemake logs
rm -rf .snakemake/log/*
```

---

## Dashboard Issues

### Dashboard Won't Start

```bash
# Check dependencies
pixi run python -c "import shiny; import plotly; import polars"

# Check for port conflicts
lsof -i :8000

# Kill existing process and retry
pkill -f "shiny run"
pixi run dashboard
```

### No Samples in Dropdown

```bash
ls results/
# Should show completed sample directories
```

### Plots Not Rendering

- Verify files exist: `ls results/{sample}/blast_results/master_blast.txt`
- Check for content: `wc -l results/{sample}/blast_results/master_blast.txt`
- Try refreshing the browser and checking the browser console for errors

### Missing Defense Score or Binomial Plots

- Defense scores require both DefenseFinder **and** co-transcription to complete first
- Binomial enrichment requires `binomial.enabled: true` in config and InterProScan to be enabled
- Check that `results/{sample}/defense_scores.tsv` and `results/{sample}/BinomialAnalysis.csv` exist

### Dashboard Access from HPC

See [Dashboard Guide](Dashboard-Guide.md) for:
- Open OnDemand node proxy setup
- SSH tunnel setup
- `pixi run dashboard-hpc` (use instead of `pixi run dashboard` on remote servers)

---

## Configuration Issues

### Config File Not Found

```bash
pwd   # must be in the RADS directory
ls config/config.yaml

# Or specify explicitly
pixi run snakemake --configfile /absolute/path/to/config.yaml --cores 8
```

### YAML Syntax Errors

**Common causes:** tabs instead of spaces, missing colons, unquoted special characters.

```bash
pixi run python -c "import yaml; yaml.safe_load(open('config/config.yaml'))"
```

---

## HPC / Network Filesystem Issues

### Symlink Errors

**Symptom:** `[Errno 95] Operation not supported: 'cacert.pem'`

**Cause:** Network filesystems (NFS, Lustre) don't support symbolic links used by Conda.

```bash
conda config --set always_copy true
rm -rf .snakemake/conda/*
```

### Pixi Not Found on Compute Nodes

```bash
export PATH="$HOME/.pixi/bin:$PATH"
```

Add this to your job script header, or see [SLURM Usage](SLURM-Usage.md) for profile-level solutions.

---

## Common Error Reference

| Error | Cause | Solution |
|-------|-------|----------|
| `MissingOutputException` | Rule didn't create expected output | Check rule log, verify input exists |
| `IncompleteFilesException` | Previous run interrupted | `rm -rf .snakemake/incomplete/` then rerun |
| `MissingInputException` | Required input file missing | Run prerequisite rules first |
| `command not found: diamond` | Tool not in PATH | Use `pixi run snakemake` |
| `SyntaxError: Non-ASCII character` | Running Python 2 | Use `pixi run python` or `python3` |
| `[Errno 95] Operation not supported` | Network filesystem symlink issue | `conda config --set always_copy true` |
| `FileNotFoundError: bin/panther/epa-ng` | InterProScan binary missing | Use `applications: "Pfam"` in config |
| `ModuleNotFoundError: No module named 'statsmodels'` | Missing dependency | Use `pixi run python`; statsmodels is in the pixi environment |

---

## Getting Help

Run with verbose output to capture more detail:
```bash
pixi run snakemake --cores 4 --verbose --printshellcmds
```

Check logs:
```bash
ls .snakemake/log/
cat logs/{sample}/blast_search.log
```

When reporting issues, include: the command run, full error traceback, config file, and relevant log files.

File issues at: https://github.com/Seandersen/RADS/issues
