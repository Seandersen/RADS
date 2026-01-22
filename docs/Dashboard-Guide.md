# Dashboard Guide

The RADS Results Explorer is an interactive web dashboard for visualizing and exploring pipeline results. The dashboard features a teal/slate color theme inspired by defense system biology figures.

## Starting the Dashboard

```bash
pixi run dashboard
```

Then open http://localhost:8000 in your web browser.

**Custom port:**
```bash
pixi run shiny run dashboard/app.py --port 9000
```

## Dashboard Overview

The dashboard has a sidebar for controls and tabbed panels for different views. The interface uses a coordinated teal/slate color palette throughout.

### Sidebar Controls

| Control | Description |
|---------|-------------|
| **Select Sample** | Choose which analysis to view |
| **Min % Identity** | Filter BLAST hits by identity |
| **Min Alignment Length** | Filter by alignment length |
| **Download BLAST Results** | Export filtered results as CSV |

## Dashboard Tabs

### 1. Summary Tab

Provides an overview of the pipeline results.

#### Value Boxes

| Metric | Description |
|--------|-------------|
| Total Genomes | Number of genomes processed |
| BLAST Hits | Total query matches found |
| Genomes with Hits | Genomes containing at least one hit |
| Co-transcribed Pairs | Downstream gene pairs identified |
| Defense Systems | Defense systems detected (if DefenseFinder enabled) |
| Hits per Mb | BLAST hit density (hits per megabase) |
| Discovery Rate (per contig) | Defense systems per contig |
| Discovery Rate (per genome) | Defense systems per genome |

#### Pipeline Overview

Shows:
- Genomes processed
- Total ORFs found
- Extracted contigs
- Domain annotations
- Query information
- Total input size

#### Hits per Genome Chart

Bar chart showing the distribution of BLAST hits across genomes.

### 2. BLAST Results Tab

Detailed exploration of BLAST search results.

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

Bar chart showing how many ORFs were predicted in each contig.

#### ORF Details Table

Interactive table with ORF information:
- ORF ID
- Parent contig
- Start/end positions
- Strand
- Length

### 4. Co-transcription Tab

Analysis of potentially co-transcribed gene pairs.

#### Co-transcribed Gene Pairs Table

Shows genes immediately downstream of query hits:

| Column | Description |
|--------|-------------|
| blast_hit_id | Original BLAST hit |
| hit_contig_orf | Hit ORF in extracted contig |
| downstream_orf | Co-transcribed gene |
| strand | Strand orientation |
| distance | Gap between genes (bp) |

#### Distance Distribution

Histogram of intergenic distances for co-transcribed pairs.

### 5. Locus Viewer Tab

Interactive gene arrow diagrams for visualizing genomic loci around query hits.

#### Contig Selection

Use the sidebar controls to filter and select contigs:

| Control | Description |
|---------|-------------|
| **Filter contigs by** | Choose viewing mode |
| - All contigs | Show all extracted contigs |
| - Contigs with query hits | Show contigs containing recombinase hits |
| - Contigs with defense systems | Show contigs with DefenseFinder annotations |
| **Select contig** | Choose specific contig to visualize |

#### Locus Visualization

The locus viewer displays genes as directional arrows showing:
- **Gene position**: Horizontal placement indicates genomic coordinates
- **Gene direction**: Arrow direction shows strand orientation (→ forward, ← reverse)
- **Gene function**: Color coding indicates annotation type

#### Color Legend

| Color | Description |
|-------|-------------|
| **Dark teal** | Query hits (recombinases) |
| **Medium teal** | Co-transcribed downstream genes |
| **Light gray-teal** | Unannotated genes |
| **Teal gradient** | Domain annotations (varies by database) |

#### Gene Information

Hover over any gene arrow to see detailed information:
- ORF ID
- Start/end coordinates
- Strand orientation
- Domain annotations (if available)
- Defense system association (if applicable)

#### ORF Details Table

Below the locus diagram, a table displays all ORFs in the selected contig with:
- ORF ID
- Contig name
- Start/end positions
- Strand
- Length

### 6. Domain Annotations Tab

InterProScan results (if enabled).

#### Top Domains Chart

Horizontal bar chart of the most frequently detected protein domains.

#### Analysis Types Pie Chart

Breakdown of annotations by database source (Pfam, CDD, SMART, etc.).

#### InterProScan Results Table

Full annotation table with:
- Protein accession
- Analysis database
- Signature accession/description
- Domain coordinates
- InterPro accession
- GO annotations

### 7. DefenseFinder Tab

Defense system detection results (if enabled).

#### Defense System Types

Bar chart of detected defense system categories (e.g., RM, Abi, CRISPR-Cas).

#### Defense Systems by Subtype

Detailed breakdown by system subtype.

#### Defense Systems Table

Full details of detected systems:
- System ID
- Type and subtype
- Position range
- Proteins involved
- Gene count

#### Defense Genes Table

Individual genes contributing to defense systems.

#### Defense System Locus Viewer

A specialized locus viewer for contigs containing defense systems:

1. **Select contig**: Choose from contigs with detected defense systems
2. **View defense system info**: See system type and subtype details
3. **Visualize locus**: Gene arrow diagram with defense genes highlighted

Defense system genes are color-coded by system type:
| System Type | Color |
|-------------|-------|
| RM (Restriction-Modification) | Darkest teal |
| CRISPR-Cas | Dark teal |
| Abi (Abortive infection) | Slate teal |
| TA (Toxin-Antitoxin) | Medium teal |
| BREX | Medium-light teal |
| DISARM | Light teal |
| Other systems | Teal gradient |

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
lsof -i :8080
```

### No Data Displayed

- Verify pipeline completed successfully
- Check that `results/{sample}/` contains expected files
- Ensure `master_blast.txt` exists and has data

### Plots Not Rendering

- Check Plotly is installed: `pixi run python -c "import plotly"`
- Try refreshing the browser
- Check browser console for JavaScript errors

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

The dashboard uses a teal/slate color palette defined in `app.py`:

```python
TEAL_PALETTE = [
    "#2d4a4a",  # Darkest teal
    "#3d5a5a",  # Dark teal
    "#4a6670",  # Slate blue-teal
    "#5a7a7a",  # Medium teal
    "#6b9090",  # Medium-light teal
    "#7a9e9e",  # Light teal
    "#8fb3b3",  # Lighter teal/sage
    "#a8c4c4",  # Light sage
    "#b5cece",  # Pale sage
]
```

To modify colors, edit the `TEAL_PALETTE`, `CHART_COLORS`, and `CUSTOM_CSS` variables in `app.py`, and the color constants in `locus_viewer.py`.

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
    # Create your plot using the teal palette
    fig = px.scatter(df.to_pandas(), x="pident", y="evalue",
                     color_discrete_sequence=TEAL_PALETTE)
    return ui.HTML(fig.to_html(include_plotlyjs="cdn", full_html=False))
```

### Customize Locus Viewer

Edit `dashboard/utils/locus_viewer.py` to modify:
- `DOMAIN_COLORS` - Colors for InterProScan domain annotations
- `DEFENSE_SYSTEM_COLORS` - Colors for DefenseFinder system types
- `QUERY_HIT_COLOR` - Color for recombinase query hits
- `DOWNSTREAM_COLOR` - Color for co-transcribed genes
- `create_gene_arrow()` - Gene arrow shape and styling
- `create_locus_figure()` - Overall locus diagram layout

### Modify Data Loading

Edit `dashboard/utils/data_loader.py` to add new data sources or modify parsing:
- `load_blast_results()` - BLAST output parsing
- `load_contig_orfs()` - ORF information from Prodigal
- `load_hit_to_contig_mapping()` - Maps BLAST hits to contig ORFs
- `load_defensefinder_systems()` - Defense system results
- `load_interproscan_results()` - Domain annotations

## Next Steps

- [[Advanced-Usage]] - Customizing the pipeline
- [[Troubleshooting]] - Common issues
