# SLURM / HPC Cluster Usage

> **Status: Experimental.** The SLURM profile has been used successfully on RMACC Alpine but has not been validated on other clusters. If you run into issues, please open a [GitHub Issue](https://github.com/Seandersen/RADS/issues) with your cluster details.

RADS ships with a SLURM profile at `profiles/slurm/` that uses the Snakemake 8 native SLURM executor (`snakemake-executor-plugin-slurm`) to submit each rule as an independent batch job. This allows per-genome steps to run in parallel across hundreds of cluster nodes simultaneously.

---

## Quick Start

```bash
# 1. Install the pixi environment on the login node (one time)
pixi install

# 2. Submit the pipeline
pixi run slurm
```

`pixi run slurm` expands to:
```bash
mkdir -p logs/slurm && snakemake --profile profiles/slurm
```

Snakemake stays running on the login node as the orchestrator; all compute happens on cluster nodes via `sbatch`.

---

## Installing Pixi on HPC

Most HPC systems do not ship pixi. Install it to your home directory:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
# Follow the prompt to add ~/.pixi/bin to your PATH in ~/.bashrc
source ~/.bashrc
```

Verify pixi is available on a compute node:
```bash
srun --pty bash -c "which pixi"
```

If your cluster does not source `~/.bashrc` for batch jobs, you have two options:

**Option A — `--export=ALL` (default for the RADS profile)**

The profile sets `slurm_extra="--export=ALL"`, which copies the full login-node environment (including `PATH`) into every job. This works on most clusters.

**Option B — explicit PATH in the profile**

If `--export=ALL` is blocked by your cluster policy, open `profiles/slurm/config.yaml` and replace the `slurm_extra` line:

```yaml
default-resources:
  - slurm_extra="--export=NONE --export=PATH=/home/YOUR_USER/.pixi/bin:$PATH"
```

> The `run_defensefinder` rule calls `pixi run -e defensefinder ...` directly. If pixi is not on PATH on compute nodes, that rule will fail. Confirm with `srun` (above) before a full run.

---

## Cluster-Specific Settings

Open `profiles/slurm/config.yaml` and adjust the `default-resources` block for your cluster:

```yaml
default-resources:
  - mem_mb=4000
  - runtime=60
  - slurm_extra="--export=ALL"
  # Uncomment and set these if your cluster requires them:
  # - slurm_partition=compute
  # - slurm_account=your_account
```

To target a specific partition for all jobs:
```yaml
  - slurm_partition=high_mem
```

To override resources for a single rule (e.g., give InterProScan more time):
```yaml
set-resources:
  run_interproscan:
    mem_mb: 64000
    runtime: 960
```

---

## What Gets Parallelized

| Rule | Parallelism |
|---|---|
| `translate_genome` | One job per genome |
| `build_database` | One job per genome |
| `blast_search` | One job per genome — 8 CPUs, 16 GB each |
| `extract_orf_ids`, `extract_coordinates`, `create_bed_file`, `extract_contig_sequences` | One job per genome |
| `run_interproscan` | Single job, 8 CPUs, 32 GB |
| `run_defensefinder` | Single job, 8 CPUs, 16 GB |
| `run_interproscan_chunk` *(binomial)* | One job per ~50k-protein chunk |

Per-genome rules scale linearly with genome count and are all submitted simultaneously (up to the `jobs: 100` cap in the profile).

---

## Binomial Analysis: Parallel InterProScan

When binomial analysis is enabled, RADS splits the whole-genome protein file into chunks and runs each chunk as its own SLURM job:

```
genomes_with_hits_cleaned.faa
        │
        ▼  split_proteins_for_binomial (checkpoint)
        │  seqkit split2 --by-size 50000
        │
   ┌────┴──────────────────────────┐
   │  chunk_001.faa  chunk_002.faa │  ← one SLURM job each
   └────┬──────────────────────────┘
        │  run_interproscan_chunk (scatter)
        ▼
   aggregate_interproscan_chunks (gather)
        ▼
   run_binomial_analysis
```

Configure chunk size in `config/config.yaml`:

```yaml
binomial:
  interproscan_chunk_size: 50000   # proteins per job
```

To skip whole-genome InterProScan entirely by supplying a precomputed file:

```yaml
binomial:
  whole_genome_interproscan: "/path/to/precomputed_ips.tsv"
```

---

## InterProScan on HPC

InterProScan is not managed by pixi — it requires a separate manual installation (see [Installation](Installation.md)). The path in `config/config.yaml` must be accessible from compute nodes:

```yaml
interproscan:
  path: "/scratch/shared/interproscan-5.76-107.0/interproscan.sh"
```

Use a shared filesystem path (not local scratch) so all nodes can reach it.

---

## Monitoring Jobs

```bash
# Watch job queue
watch squeue -u $USER

# Per-job stdout/stderr (named by SLURM job ID)
ls logs/slurm/

# Per-rule logs
ls logs/{sample_name}/
```

---

## Resuming After Failure

```bash
pixi run slurm
```

Snakemake skips completed outputs and resubmits only what is missing or incomplete. If partial outputs exist:

```bash
mkdir -p logs/slurm && snakemake --profile profiles/slurm --rerun-incomplete
```

---

## Dry Run

Preview which jobs will be submitted without actually submitting:

```bash
pixi run dry-run
# or, to see SLURM resource annotations:
snakemake --profile profiles/slurm -n
```

---

## Adjusting Concurrency

The default cap is 100 simultaneous jobs. Change `jobs:` in `profiles/slurm/config.yaml`:

```yaml
jobs: 50    # conservative
jobs: 200   # if your allocation allows
```
