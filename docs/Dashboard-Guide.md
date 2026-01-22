# Dashboard Guide

The RADS Results Explorer is an interactive web dashboard for visualizing and exploring pipeline results.

## Starting the Dashboard

```bash
pixi run dashboard
```

Then open http://localhost:8080 in your web browser.

**Custom port:**
```bash
pixi run shiny run dashboard/app.py --port 9000
```

## Dashboard Overview

The dashboard has a sidebar for controls and tabbed panels for different views.

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

### 5. Domain Annotations Tab

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

### 6. DefenseFinder Tab

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

The dashboard code is in `dashboard/app.py`. You can modify:

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
    # Create your plot
    fig = px.scatter(df.to_pandas(), x="pident", y="evalue")
    return ui.HTML(fig.to_html(include_plotlyjs="cdn", full_html=False))
```

### Modify Data Loading

Edit `dashboard/utils/data_loader.py` to add new data sources or modify parsing.

## Next Steps

- [[Advanced-Usage]] - Customizing the pipeline
- [[Troubleshooting]] - Common issues
