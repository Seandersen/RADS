# Troubleshooting Guide

This page covers common issues and their solutions.

## Installation Issues

### Pixi Won't Install

**Symptom:** `curl: command not found` or download fails

**Solution:**
```bash
# On macOS, install curl via Homebrew
brew install curl

# Alternative installation method
wget -qO- https://pixi.sh/install.sh | bash
```

### Dependencies Fail to Install

**Symptom:** `pixi install` shows solver errors

**Solution:**
```bash
# Clear cache and retry
rm -rf ~/.pixi/cache
pixi install

# If specific package fails, try pinning version
# Edit pixi.toml and specify exact version
```

### DefenseFinder Module Error

**Symptom:** `ModuleNotFoundError: No module named 'macsypy'`

**Cause:** Python version incompatibility (Python 3.13)

**Solution:**
```bash
# Option 1: Disable DefenseFinder in config
# Edit config/config.yaml:
defensefinder:
  enabled: false

# Option 2: The pipeline handles this gracefully
# Empty output files will be created automatically
```

## Pipeline Execution Issues

### "Nothing to be done"

**Symptom:** Snakemake says nothing to do, but outputs are missing

**Solutions:**

1. **Check if outputs exist:**
   ```bash
   ls -la results/{sample}/blast_results/master_blast.txt
   ```

2. **Force regeneration:**
   ```bash
   pixi run snakemake --forceall --cores 8
   ```

3. **Clear Snakemake metadata:**
   ```bash
   rm -rf .snakemake/metadata
   pixi run snakemake --cores 8
   ```

4. **Explicitly target outputs:**
   ```bash
   pixi run snakemake results/{sample}/metrics/pipeline_metrics.json --cores 8
   ```

### Genome Download Fails

**Symptom:** `download_genomes` rule fails

**Solutions:**

1. **Check accession format:**
   ```bash
   # Valid formats:
   NC_000913.3
   NZ_CP088158.1
   CP016471.1
   GCF_000001405.40
   ```

2. **Check internet connection:**
   ```bash
   pixi run datasets summary genome accession NC_000913.3
   ```

3. **Download manually:**
   ```bash
   pixi run datasets download genome accession NC_000913.3 --filename test.zip
   ```

4. **Use local genomes instead:**
   ```yaml
   # config/config.yaml
   download:
     enabled: false
   genomes_path: "path/to/local/genomes"
   ```

### BLAST Search Returns No Hits

**Symptom:** `master_blast.txt` is empty or has only header

**Solutions:**

1. **Check query format:**
   ```bash
   head resources/query.fa
   # Should be amino acid FASTA:
   # >protein_name
   # MKTQPIKVN...
   ```

2. **Lower identity threshold:**
   ```yaml
   diamond:
     identity: 20  # Lower from default 30
   ```

3. **Verify databases built correctly:**
   ```bash
   ls -la results/{sample}/diamond_dbs/*.dmnd
   ```

4. **Test BLAST manually:**
   ```bash
   pixi run diamond blastp -d results/{sample}/diamond_dbs/genome1.dmnd \
       -q resources/query.fa -o test_blast.txt
   ```

### InterProScan Fails

**Symptom:** Rule fails with InterProScan errors

**Solutions:**

1. **Verify InterProScan installation:**
   ```bash
   /path/to/interproscan.sh --version
   ```

2. **Check Java version:**
   ```bash
   java -version  # Needs Java 11+
   ```

3. **Disable if not needed:**
   ```yaml
   interproscan:
     enabled: false
   ```

4. **Check path is absolute:**
   ```yaml
   interproscan:
     path: "/absolute/path/to/interproscan.sh"  # Not relative
   ```

### Out of Memory

**Symptom:** Jobs killed with memory errors

**Solutions:**

1. **Reduce parallelism:**
   ```bash
   pixi run snakemake --cores 2  # Fewer parallel jobs
   ```

2. **Process fewer genomes:**
   ```yaml
   download:
     max_genomes: 50  # Limit initial run
   ```

3. **Increase system swap:**
   ```bash
   # Linux
   sudo fallocate -l 8G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile
   ```

### Disk Space Issues

**Symptom:** "No space left on device"

**Solutions:**

1. **Clean intermediate files:**
   ```bash
   pixi run clean  # Removes results and logs
   ```

2. **Remove downloaded genomes after staging:**
   ```bash
   rm -rf results/{sample}/downloaded_genomes
   ```

3. **Clean Snakemake logs:**
   ```bash
   rm -rf .snakemake/log/*
   ```

## Dashboard Issues

### Dashboard Won't Start

**Symptom:** Error when running `pixi run dashboard`

**Solutions:**

1. **Check dependencies:**
   ```bash
   pixi run python -c "import shiny; import plotly; import polars"
   ```

2. **Check port availability:**
   ```bash
   # Kill existing process
   pkill -f "shiny run"

   # Use different port
   pixi run shiny run dashboard/app.py --port 9000
   ```

3. **Check for Python errors:**
   ```bash
   pixi run python dashboard/app.py
   ```

### No Samples in Dropdown

**Symptom:** Sample selector is empty

**Solution:**
```bash
# Ensure results directory exists with completed analyses
ls results/
# Should show sample directories like: efb0058_hits/
```

### Plots Not Rendering

**Symptom:** Blank charts or "No data available"

**Solutions:**

1. **Check data files exist:**
   ```bash
   ls results/{sample}/blast_results/master_blast.txt
   ```

2. **Verify data format:**
   ```bash
   head results/{sample}/blast_results/master_blast.txt
   ```

3. **Test data loading:**
   ```bash
   pixi run python -c "
   import sys
   sys.path.insert(0, 'dashboard')
   from utils.data_loader import load_blast_results
   df = load_blast_results('results/{sample}')
   print(len(df) if df is not None else 'No data')
   "
   ```

## Configuration Issues

### Config File Not Found

**Symptom:** `FileNotFoundError: config/config.yaml`

**Solution:**
```bash
# Ensure you're in the RADS directory
pwd  # Should be /path/to/RADS

# Check config exists
ls config/config.yaml

# Or specify config explicitly
pixi run snakemake --configfile /absolute/path/to/config.yaml --cores 8
```

### YAML Syntax Errors

**Symptom:** `yaml.scanner.ScannerError`

**Common causes:**
- Tabs instead of spaces (use 2 spaces)
- Missing colons
- Unquoted special characters

**Solution:** Validate YAML:
```bash
pixi run python -c "import yaml; yaml.safe_load(open('config/config.yaml'))"
```

## Getting Help

### Debug Mode

Run with verbose output:
```bash
pixi run snakemake --cores 4 --verbose --printshellcmds
```

### Check Logs

```bash
# Snakemake logs
ls .snakemake/log/

# Rule-specific logs
cat logs/{sample}/blast_search.log
cat logs/{sample}/defensefinder.log
```

### Report Issues

When reporting issues, include:

1. **Command run:**
   ```bash
   pixi run snakemake --cores 8
   ```

2. **Error message:** (full traceback)

3. **Config file:** (sanitized)

4. **System info:**
   ```bash
   uname -a
   pixi run snakemake --version
   pixi run python --version
   ```

5. **Log files:** Attach relevant logs from `logs/` directory

File issues at: https://github.com/Seandersen/RADS/issues

## Common Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| `MissingOutputException` | Rule didn't create expected output | Check rule log, verify input exists |
| `ProtectedOutputException` | Trying to overwrite protected file | Use `--forceall` or delete output |
| `MissingInputException` | Required input file missing | Run prerequisite rules first |
| `WorkflowError: Conda environment file cannot be found` | Missing env file | Check workflow/envs/ exists |
| `command not found: diamond` | Tool not in PATH | Run via `pixi run snakemake` |
