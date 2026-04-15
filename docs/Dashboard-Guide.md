# Dashboard Guide

The RADS Results Explorer is an interactive web dashboard for visualizing and exploring pipeline results. For sharing results without running a server, see [Static HTML Report](#static-html-report) below.

---

## Starting the Dashboard

### Local Machine

```bash
pixi run dashboard
```

Then open `http://localhost:8000` in your browser.

**Custom port:**
```bash
pixi run shiny run dashboard/app.py --port 9000
```

### HPC / Remote Server

Use `dashboard-hpc` (disables auto-reload, which is not needed on a server):

```bash
pixi run dashboard-hpc
```

Then access the dashboard using one of the methods below depending on your cluster setup.

#### Option A — Open OnDemand (OOD) Node Proxy

If your cluster provides Open OnDemand (e.g., RMACC Alpine, OSC, TACC), you can access the dashboard directly through the browser without an SSH tunnel.

1. Open an OOD terminal (Interactive Apps → Clusters → Shell)
2. Request a compute node and start the dashboard:
   ```bash
   salloc --nodes=1 --ntasks=1 --time=2:00:00
   pixi run dashboard-hpc
   ```
3. Note the hostname of your compute node:
   ```bash
   hostname   # e.g., c3cpu-a7-u26-3
   ```
4. In your browser, navigate to:
   ```
   https://<your-ood-url>/rnode/<hostname>/8000/
   ```
   For example, on RMACC Alpine:
   ```
   https://ondemand-rmacc.rc.colorado.edu/rnode/c3cpu-a7-u26-3/8000/
   ```

> **Note:** The dashboard must be running on the compute node whose hostname appears in the URL. Login nodes are typically firewalled and will not work.

#### Option B — SSH Tunnel

When Open OnDemand is not available, create an SSH tunnel from your local machine:

**Step 1:** Start the dashboard on the remote server:
```bash
pixi run dashboard-hpc
# or (conda):
shiny run dashboard/app.py --host 0.0.0.0 --port 8000
```

**Step 2:** On your **local machine**, open a new terminal and create the tunnel:
```bash
ssh -L 8000:localhost:8000 username@remote-server
```

**Step 3:** Open `http://localhost:8000` in your local browser.

#### Keeping the Dashboard Running (screen/tmux)

```bash
screen -S dashboard
pixi run dashboard-hpc
# Detach: Ctrl+A then D
# Reattach: screen -r dashboard
```

---

## Static HTML Report

For sharing results without running a server, generate a self-contained HTML file:

```bash
pixi run python workflow/scripts/generate_report.py \
    --results results/my_analysis \
    --output  results/my_analysis/report.html
```

The report opens in any browser with no server or internet connection required (charts load from CDN on first open).

### Optional: Locus Viewer

Add `--include-locus-viewer` to embed interactive gene-arrow diagrams for all contigs. This increases file size substantially (tens of MB for large runs) but requires no additional tools to view:

```bash
pixi run python workflow/scripts/generate_report.py \
    --results results/my_analysis \
    --output  results/my_analysis/report.html \
    --include-locus-viewer \
    --locus-max-contigs 200   # default 250; reduce if file is too large
```

### Report Contents

The HTML report includes five tabs:

| Tab | Contents |
|-----|----------|
| **Summary** | Metric cards, pipeline funnel, input/output genomes by genus, defense category bar |
| **BLAST Results** | Identity histogram, hits-per-genome histogram, scatter plot, top hits table |
| **Co-transcription** | Gap distance histogram, top domains in co-transcribed ORFs, defense scores (top 50 / bottom 50) |
| **DefenseFinder** | Category bar, category × type sunburst, top system types, full systems table |
| **Domain Annotations** | Top IPS domains, analysis-type pie, binomial enrichment bar and table |
| **Locus Viewer** *(optional)* | Gene-arrow diagrams per contig, navigable by dropdown or Prev/Next |

---

## Dashboard Overview

The live dashboard has a sidebar for global controls and seven tabbed panels.

### Sidebar Controls

| Control | Description |
|---------|-------------|
| **Select Sample** | Choose which analysis to view |
| **Min % Identity** | Filter BLAST hits by minimum percent identity |
| **Min Alignment Length** | Filter BLAST hits by minimum alignment length |
| **Download BLAST Results** | Export filtered results as CSV |

All visualizations update in real-time as filters are adjusted.

---

## Dashboard Tabs

### 1. Summary

High-level overview: metric cards, Sankey diagram of data flow, key statistics.

### 2. BLAST Results

- Identity vs. length scatter plot
- Identity distribution histogram
- Sortable, filterable hits table

**Columns:** query_id, subject_id, length, nident, pident, evalue, genome

### 3. Contig Analysis

- ORF length distribution histogram
- ORFs per contig bar chart
- ORF details table

### 4. Co-transcription

Analysis of genes immediately downstream of query hits on the same strand.

- Co-transcribed gene pairs table
- Distance distribution histogram
- Domain annotations for co-transcribed genes
- DefenseFinder hits in co-transcribed genes
- Defense score distribution
- Binomial domain enrichment (requires binomial analysis to be enabled)

### 5. Locus Viewer

Interactive gene-arrow diagrams for genomic loci around query hits.

**Sidebar filters:** contig scope, defense score range, binomial p-value threshold, MGE filter, defense type

**Color legend:**

| Color | Description |
|-------|-------------|
| Dark teal-black | Query hits (recombinases) |
| Medium green-teal | Co-transcribed downstream genes |
| Light blue-teal | DefenseFinder defense genes |
| Dark slate-blue | Mobile genetic elements |
| Gray | Other InterProScan domains |
| White | Unannotated genes |

Hover any gene arrow to see ORF ID, coordinates, strand, length, and annotation details.

### 6. Domain Annotations

- Top domains across all co-transcribed ORFs
- Analysis-type breakdown pie chart
- Full binomial enrichment results (sorted by adjusted p-value)
- Complete InterProScan results table

### 7. DefenseFinder

- Defense system categories bar chart
- Subtype breakdown
- Systems and genes tables
- Defense locus viewer (contigs filtered to those with defense annotations)

---

## Sample Selection

The **Select Sample** dropdown shows all completed analyses in the `results/` directory. To add a new sample, run the pipeline with a new `sample_name`, then refresh the dashboard and select it.

---

## Troubleshooting

**Dashboard won't start:**
```bash
pixi run python -c "import shiny; import plotly; import polars"
lsof -i :8000   # check for port conflicts
```

**No data displayed:**
- Verify pipeline completed: `ls results/{sample}/blast_results/master_blast.txt`
- Check the file has content: `wc -l results/{sample}/blast_results/master_blast.txt`

**Missing defense score or binomial plots:**
- Defense scores require both DefenseFinder and co-transcription to complete
- Binomial enrichment requires `BinomialAnalysis.csv` — check that `binomial.enabled: true` in config

**Slow performance:** Use identity/length filters in the sidebar to reduce the displayed dataset.

---

## Next Steps

- [Advanced Usage](Advanced-Usage.md) — customizing the pipeline
- [Troubleshooting](Troubleshooting.md) — common issues
