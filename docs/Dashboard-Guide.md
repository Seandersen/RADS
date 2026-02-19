# Dashboard Guide

The RADS Results Explorer is an interactive web dashboard for visualizing and exploring pipeline results. The dashboard features a teal/slate color theme inspired by defense system biology figures.

## Starting the Dashboard

### With Pixi (Local Machine)

```bash
pixi run dashboard
```

Then open http://localhost:8000 in your web browser.

**Custom port:**
```bash
pixi run shiny run dashboard/app.py --port 9000
```

### Without Pixi (HPC / Remote Servers)

If pixi isn't available, install dependencies with conda/pip:

```bash
# Create environment (one-time setup)
conda create -n rads-dashboard python=3.10 -c conda-forge -y
conda activate rads-dashboard
pip install shiny polars plotly pandas pyarrow

# Run dashboard
shiny run dashboard/app.py --port 8000
```

### Remote Access via SSH Tunnel

When running on a remote server (HPC, cloud, etc.), you need an SSH tunnel to access the dashboard in your local browser.

**Step 1:** Start the dashboard on the remote server:
```bash
shiny run dashboard/app.py --port 8000
```

**Step 2:** On your **local machine**, create an SSH tunnel (new terminal):
```bash
ssh -L 8000:localhost:8000 username@remote-server
```

**Step 3:** Open http://localhost:8000 in your local browser.

### Keep Dashboard Running (screen/tmux)

To keep the dashboard running after disconnecting from SSH:

```bash
# Start a screen session
screen -S dashboard

# Run the dashboard
shiny run dashboard/app.py --port 8000

# Detach from screen: press Ctrl+A, then D
# The dashboard keeps running in the background

# Reattach later:
screen -r dashboard
```

## Dashboard Overview

The dashboard has a sidebar for global controls and seven tabbed panels for different views. The interface uses a coordinated teal/slate color palette throughout.

### Sidebar Controls

| Control | Description |
|---------|-------------|
| **Select Sample** | Choose which analysis to view |
| **Min % Identity** | Filter BLAST hits by minimum percent identity |
| **Min Alignment Length** | Filter BLAST hits by minimum alignment length |
| **Download BLAST Results** | Export filtered results as CSV |

All visualizations update in real-time as filters are adjusted.

## Dashboard Tabs

### 1. Summary Tab

Provides a high-level overview of the pipeline run.

#### Value Boxes

| Metric | Description |
|--------|-------------|
| Total Genomes | Number of genomes processed |
| BLAST Hits | Total query matches found (after filtering) |
| Genomes with Hits | Genomes containing at least one hit |
| Hits per Mb | BLAST hit density (hits per megabase of input sequence) |
| Defense Systems | Defense systems detected by DefenseFinder |
| Co-transcribed Pairs | Downstream gene pairs identified |
| Discovery Rate (per contig) | Known defense systems per contig with a hit |
| Discovery Rate (per genome) | Fraction of hit-containing genomes with a defense system |

#### Pipeline Overview

Text summary of key pipeline statistics, including total genomes, ORFs, contigs, annotations, and query information.

#### Pipeline Results Breakdown

A Sankey diagram showing how data flows through the pipeline:
- Genomes → Genomes with Hits / No Hits
- Hits → Contigs → ORFs
- ORFs → BLAST Hits / Co-transcribed Pairs / Defense Systems / Domain Annotations

### 2. BLAST Results Tab

Detailed exploration of Diamond BLAST search results.

#### Identity vs Length Scatter Plot

- **X-axis**: Alignment length
- **Y-axis**: Percent identity
- **Color**: Source genome
- **Hover**: Query ID, subject ID, E-value

Use this to identify:
- High-confidence hits (high identity, long alignment)
- Partial matches (high identity, short alignment)
- Distant homologs (low identity, long alignment)

#### Identity Distribution Histogram

Shows the distribution of percent identity values across all hits.

#### BLAST Hits Table

Interactive table with all BLAST results:
- Sortable columns
- Filterable by any column
- Full-screen mode available

Columns:
| Column | Description |
|--------|-------------|
| query_id | Query sequence name |
| subject_id | Hit protein ID |
| length | Alignment length |
| nident | Identical positions |
| pident | Percent identity |
| evalue | E-value |
| genome | Source genome |

### 3. Contig Analysis Tab

Analysis of extracted flanking regions.

#### ORF Length Distribution

Histogram showing the size distribution of predicted ORFs in extracted contigs.

#### ORFs per Contig

Bar chart showing how many ORFs were predicted in each contig, with a median line.

#### ORF Details Table

Interactive table with ORF information:
- ORF ID
- Parent contig
- Start/end positions
- Strand
- Length

### 4. Co-transcription Tab

Analysis of genes immediately downstream of query hits that are likely co-transcribed, along with statistical enrichment and defense association scoring.

#### Co-transcribed Gene Pairs Table

Shows genes immediately downstream of query hits on the same strand and within the distance threshold:

| Column | Description |
|--------|-------------|
| blast_hit_id | Original BLAST hit identifier |
| hit_contig_orf | Query hit ORF in extracted contig |
| downstream_orf | Co-transcribed downstream gene |
| strand | Strand orientation (+1 / -1) |
| distance | Intergenic gap between genes (bp) |

#### Distance Distribution

Histogram of intergenic distances for all co-transcribed pairs.

#### Domain Annotations for Co-transcribed Genes

Horizontal bar chart and table showing which InterProScan domains are found in co-transcribed genes, along with annotation source (Pfam, TIGRFAM, CDD, etc.).

#### DefenseFinder Hits in Co-transcribed Genes

Summary and table of co-transcribed genes that are also annotated as defense genes by DefenseFinder. These represent cases where the query hit is directly adjacent to a known defense system.

#### Defense Score Distribution

Plot of defense scores for all co-transcribed genes. The defense score quantifies how associated a co-transcribed gene is with nearby defense systems based on spatial proximity and local defense gene density. Low scores indicate genes that are isolated from known defense islands and would be missed by traditional proximity-based detection.

See [Pipeline Overview](Pipeline-Overview.md#step-12-calculate-defense-scores) for details on how the score is calculated.

#### Binomial Domain Enrichment

Bar chart and table of Pfam domains that are statistically enriched in the extracted contigs compared to whole-genome background frequencies. Domains are ranked by adjusted p-value (Benjamini-Hochberg correction). Only available when the binomial analysis step is enabled in the pipeline.

### 5. Locus Viewer Tab

Interactive gene arrow diagrams for visualizing genomic loci around query hits.

#### Sidebar Controls

| Control | Description |
|---------|-------------|
| **Filter contigs by** | Choose viewing scope |
| - All contigs | Show all extracted contigs |
| - Contigs with query hits | Show contigs containing recombinase hits |
| - Contigs with defense systems | Show contigs with DefenseFinder annotations |
| **Select / Deselect All** | Quickly select or clear the contig list |
| **Select contig(s)** | Choose one or more specific contigs to visualize |
| **Defense Score Range** | Show only contigs whose co-transcribed genes fall within the selected score range (0–1) |
| **Binomial p-value threshold** | Filter contigs to those with at least one co-transcribed gene whose top domain meets the enrichment threshold |
| **Only co-transcribed genes** | Highlight contigs that have at least one co-transcribed downstream gene |
| **Only defense systems** | Highlight contigs that have at least one DefenseFinder defense gene |
| **Only MGEs** | Highlight contigs with mobile genetic element domain annotations |
| **Defense type filter** | Filter to contigs containing specific defense system types |
| **Reset Filters** | Clear all locus filters back to defaults |

#### Locus Visualization

Genes are drawn as directional arrows:
- **Position**: Horizontal placement reflects genomic coordinates
- **Direction**: Arrow orientation shows strand (→ forward, ← reverse)
- **Color**: Indicates gene annotation type (see color legend below)

Multiple contigs can be selected and are stacked vertically, each in its own panel.

#### Color Legend

| Color | Description |
|-------|-------------|
| **Dark teal-black** | Query hits (recombinases) — highest priority |
| **Medium green-teal** | Co-transcribed downstream genes |
| **Light blue-teal** | DefenseFinder defense system genes |
| **Dark slate-blue** | Mobile genetic elements (transposases, integrases, etc.) |
| **Gray** | Other InterProScan domain annotations |
| **White** | Unannotated genes |

#### Hover Annotations

Hovering over any gene arrow shows:
- ORF ID
- Genomic coordinates (start–end)
- Strand
- Gene length
- **Annotation**: For co-transcribed genes, lists all InterProScan domain hits found for that gene (database and description, one per line). Falls back to "Co-transcribed downstream" if no domain annotations are available. For other gene types, shows the relevant annotation label.

#### ORF Details Table

Below the locus diagram, a table lists all ORFs in the selected contig with their coordinates, strand, and length.

### 6. Domain Annotations Tab

InterProScan results for all ORFs in extracted contigs (requires InterProScan to be enabled).

#### Top Domains Chart

Horizontal bar chart of the most frequently detected protein domains across all contigs.

#### Analysis Types Pie Chart

Breakdown of annotations by database source (Pfam, CDD, TIGRFAM, SMART, Gene3D, SUPERFAMILY, PANTHER, ProSiteProfiles, Hamap, etc.).

#### Binomial Domain Enrichment (All Domains)

Bar chart of all statistically enriched Pfam domains from the binomial analysis, sorted by adjusted p-value. Complements the Co-transcription tab view, which focuses on domains found specifically in co-transcribed genes.

#### InterProScan Results Table

Full annotation table with:
- Protein accession
- Analysis database
- Signature accession and description
- Domain coordinates
- InterPro accession
- GO annotations

### 7. DefenseFinder Tab

Defense system detection results (requires DefenseFinder to be enabled).

#### Defense System Types

Bar chart of detected defense system categories (RM, Abi, CRISPR-Cas, TA, BREX, DISARM, etc.).

#### Defense Systems by Subtype

Detailed breakdown by system subtype.

#### Defense Systems Table

Full details of detected systems:
- System ID
- Type and subtype
- Position range on contig
- Gene count

#### Defense Genes Table

Individual genes contributing to each defense system, with gene name, system type and subtype, and contig location.

#### Defense System Locus Viewer

A specialized locus viewer focused on contigs containing defense systems:

1. **Filter by defense type**: Optionally restrict the contig list to a specific defense system type
2. **Select a contig**: Choose from contigs that have at least one DefenseFinder defense system
3. **View system info**: See the defense system type and subtype for the selected contig
4. **Visualize the locus**: Gene arrow diagram with genes color-coded by annotation

**Color coding in the defense locus viewer:**

| Color | Description |
|-------|-------------|
| **Dark teal-black** | Query hits (recombinases) — same color as in the main locus viewer |
| **Navy blue gradient** | DefenseFinder defense genes, shaded by system type (see table below) |
| **White** | Genes not part of a defense system |

Defense system gene colors:
| System Type | Color |
|-------------|-------|
| RM (Restriction-Modification) | Darkest navy |
| CRISPR-Cas | Dark navy |
| Abi (Abortive infection) | Navy |
| TA (Toxin-Antitoxin) | Darkest navy |
| BREX | Medium navy |
| DISARM | Navy-teal |
| Other systems | Navy gradient |

The query hit coloring in the defense locus viewer uses the same `QUERY_HIT_COLOR` as the main locus viewer, so it is easy to identify the recombinase in the context of each defense system.

## Interactive Features

### Filtering

Use sidebar sliders to filter BLAST results:
- Drag **Min % Identity** to show only high-confidence hits
- Adjust **Min Alignment Length** to filter short alignments

All visualizations update in real-time.

### Full-Screen Mode

Click the expand icon on any card to view charts in full-screen mode.

### Data Export

Click **Download BLAST Results** to export filtered BLAST hits as a CSV file.

### Table Interactions

- Click column headers to sort
- Use filter boxes to search within columns
- Scroll horizontally for wide tables

## Sample Selection

The **Select Sample** dropdown shows all available results in the `results/` directory.

To analyze a new sample:
1. Run the pipeline with a new `sample_name`
2. Refresh the dashboard
3. Select the new sample from the dropdown

## Troubleshooting

### Dashboard Won't Start

```bash
# Check Shiny is installed
pixi run python -c "import shiny; print(shiny.__version__)"

# Check for port conflicts
lsof -i :8000
```

### No Data Displayed

- Verify pipeline completed successfully
- Check that `results/{sample}/` contains expected files
- Ensure `master_blast.txt` exists and has data

### Plots Not Rendering

- Check Plotly is installed: `pixi run python -c "import plotly"`
- Try refreshing the browser
- Check browser console for JavaScript errors

### Missing Defense Score or Binomial Plots

- Defense scores require DefenseFinder and co-transcription analysis to complete first
- Binomial enrichment requires `binomial.enabled: true` in `config/config.yaml` and InterProScan to be enabled
- Check that `defense_scores.tsv` and `BinomialAnalysis.csv` exist in the results directory

### Slow Performance

For large datasets:
- Use identity/length filters to reduce displayed data
- Consider running on a subset of genomes first

## Customizing the Dashboard

The dashboard code is organized in `dashboard/`:

```
dashboard/
├── app.py                 # Main application and UI
└── utils/
    ├── data_loader.py     # Data loading functions
    └── locus_viewer.py    # Locus visualization components
```

### Color Theme

The dashboard uses a teal/slate color palette defined in `app.py`. The locus viewer color constants are in `locus_viewer.py`:

```python
QUERY_HIT_COLOR = "#2d3d3d"          # Dark teal-black for recombinases
DOWNSTREAM_COLOR = "#5a8a7a"         # Green-teal for co-transcribed genes
DEFENSE_COLOR = "#6a9eae"            # Light blue-teal for defense genes
MGE_COLOR = "#4a5a6a"                # Dark slate-blue for MGEs
INTERPROSCAN_OTHER_COLOR = "#b0b0b0" # Gray for other domain annotations
NO_ANNOTATION_COLOR = "#ffffff"      # White for unannotated genes
```

To modify colors, edit these constants in `locus_viewer.py` and the `TEAL_PALETTE` / `CUSTOM_CSS` variables in `app.py`.

### Add New Visualizations

Edit the UI in `app_ui` and add corresponding server functions:

```python
# In app_ui, add a new card:
ui.card(
    ui.card_header("My Custom Plot"),
    ui.output_ui("my_custom_plot"),
)

# In server, add the render function:
@render.ui
def my_custom_plot():
    df = filtered_blast_data()
    if df is None:
        return ui.p("No data")
    fig = px.scatter(df.to_pandas(), x="pident", y="evalue",
                     color_discrete_sequence=TEAL_PALETTE)
    return ui.HTML(fig.to_html(include_plotlyjs=False, full_html=False))
```

### Customize Locus Viewer

Edit `dashboard/utils/locus_viewer.py` to modify:
- `DOMAIN_COLORS` - Colors for InterProScan domain annotations
- `DEFENSE_SYSTEM_COLORS` - Colors for DefenseFinder system types
- `QUERY_HIT_COLOR` - Color for recombinase query hits (shared by both locus viewers)
- `DOWNSTREAM_COLOR` - Color for co-transcribed genes
- `MGE_COLOR` - Color for mobile genetic element genes
- `create_gene_arrow()` - Gene arrow shape and styling
- `create_locus_figure()` - Main locus diagram
- `create_defense_locus_figure()` - Defense system locus diagram

### Modify Data Loading

Edit `dashboard/utils/data_loader.py` to add new data sources or modify parsing:
- `load_blast_results()` - BLAST output parsing
- `load_contig_orfs()` - ORF information from Prodigal
- `load_hit_to_contig_mapping()` - Maps BLAST hits to contig ORFs
- `load_defensefinder_systems()` - Defense system results
- `load_interproscan_results()` - Domain annotations
- `load_defense_scores()` - Defense scores for co-transcribed genes
- `load_binomial_results()` - Binomial domain enrichment results

## Next Steps

- [[Advanced-Usage]] - Customizing the pipeline
- [[Troubleshooting]] - Common issues
