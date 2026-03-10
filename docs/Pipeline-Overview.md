# Pipeline Overview

This page explains each step of the RADS Snakemake pipeline in detail.

## Pipeline Architecture

```
┌─────────────────┐     ┌─────────────────┐
│  NCBI Download  │ OR  │  Local Genomes  │
└────────┬────────┘     └────────┬────────┘
         │                       │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │   Stage Genomes       │
         │   (organize .fna)     │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │  Translate Genomes    │
         │    (Prodigal)         │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │  Build BLAST DBs      │
         │    (Diamond)          │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │   BLAST Search        │
         │  Query vs Genomes     │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │  Extract Contigs      │
         │ (flanking regions)    │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │  Translate Contigs    │
         │   (ORF prediction)    │
         └───────────┬───────────┘
                     │
         ┌───────────┼───────────┐
         ▼           ▼           ▼
┌─────────────┐ ┌─────────────┐ ┌─────────────┐
│InterProScan │ │Co-transcr.  │ │DefenseFinder│
│  (domains)  │ │  Analysis   │ │  (defense)  │
└──────┬──────┘ └──────┬──────┘ └──────┬──────┘
       │               │               │
       │          ┌────┴────┐          │
       │          ▼         ▼          │
       │   ┌──────────┐ ┌──────────┐  │
       │   │ Defense  │ │ Binomial │  │
       │   │  Score   │ │ Analysis │  │
       │   └──────────┘ └──────────┘  │
       │               │               │
       └───────────────┼───────────────┘
                       ▼
              ┌─────────────────┐
              │ Pipeline Metrics│
              └─────────────────┘
```

## Phase 1: Genome Preparation

### Step 1: Download Genomes

**Rule:** `download_genomes`
**Location:** `workflow/rules/download_genomes.smk`

Downloads genome sequences from NCBI using the accession list.

**Input:**
- Accession file (one accession per line)

**Output:**
- `results/{sample}/downloaded_genomes/ncbi_dataset/data/`

**Tools:** NCBI datasets CLI, efetch

**Notes:**
- Supports both RefSeq (NZ_*, NC_*) and GenBank (CP*, GCA_*) accessions
- Handles failed downloads gracefully
- Can be skipped if using local genomes

### Step 2: Stage Genomes

**Rule:** `stage_genomes`
**Location:** `workflow/rules/stage_genomes.smk`

Organizes genome files into a standardized structure and creates the genome manifest.

**Input:**
- Downloaded or local genome directories

**Output:**
- `results/{sample}/genomes/*.fna` - Organized genome files
- `results/{sample}/genome_manifest.txt` - List of all genome IDs

### Step 3: Translate Genomes

**Rule:** `translate_genome`
**Location:** `workflow/rules/translate.smk`

Predicts ORFs in each genome using Prodigal.

**Input:**
- `results/{sample}/genomes/{genome}.fna`

**Output:**
- `results/{sample}/translated/{genome}.faa` - Protein sequences

**Tools:** Prodigal

**Parameters:**
- `mode`: "meta" (default) or "single"
- Meta mode is recommended for diverse genome sets

### Step 4: Build BLAST Databases

**Rule:** `build_database`
**Location:** `workflow/rules/build_databases.smk`

Creates Diamond BLAST databases for each genome.

**Input:**
- `results/{sample}/translated/{genome}.faa`

**Output:**
- `results/{sample}/diamond_dbs/{genome}.dmnd`

**Tools:** Diamond makedb

## Phase 2: BLAST and Extraction

### Step 5: BLAST Search

**Rule:** `blast_search`
**Location:** `workflow/rules/blast_search.smk`

Searches query sequence against all genome databases.

**Input:**
- Query FASTA file
- Diamond databases

**Output:**
- `results/{sample}/blast_results/{genome}_blast.txt` - Per-genome results
- `results/{sample}/blast_results/master_blast.txt` - Combined results

**Tools:** Diamond blastp

**Output Columns:**
| Column | Description |
|--------|-------------|
| query_id | Query sequence identifier |
| subject_id | Hit protein identifier |
| length | Alignment length |
| nident | Number of identical matches |
| pident | Percent identity |
| evalue | E-value |
| genome | Source genome |

### Step 6: Extract Contigs

**Rule:** `extract_contigs`
**Location:** `workflow/rules/extract_contigs.smk`

Extracts flanking genomic regions around each BLAST hit.

**Input:**
- BLAST results
- Original genome sequences

**Output:**
- `results/{sample}/bed_files/{genome}.bed` - Extraction coordinates
- `results/{sample}/contigs/{genome}_contigs.fna` - Extracted sequences
- `results/{sample}/all_contigs.fna` - Combined contigs
- `results/{sample}/all_contigs_filtered.fna` - Filtered (non-empty) contigs

**Tools:** SeqKit subseq

**Process:**
1. Parse BLAST hits to get hit coordinates
2. Calculate upstream/downstream boundaries using `upstream_nt` and `downstream_nt` config values
3. Handle chromosome boundaries (no wraparound)
4. Extract sequences using BED coordinates

### Step 7: Translate Contigs

**Rule:** `translate_contigs`
**Location:** `workflow/rules/process_orfs.smk`

Predicts ORFs in extracted contigs using Prodigal.

**Input:**
- `results/{sample}/all_contigs_filtered.fna`

**Output:**
- `results/{sample}/contig_orfs/all_contigs.faa` - Predicted proteins
- `results/{sample}/contig_orfs/all_contigs_prodigal.txt` - Prodigal output

**Tools:** Prodigal (meta mode)

## Phase 3: Annotation and Analysis

### Step 8: InterProScan (Optional)

**Rule:** `run_interproscan`
**Location:** `workflow/rules/process_orfs.smk`

Annotates protein domains in contig ORFs using InterProScan.

**Input:**
- `results/{sample}/contig_orfs/interproscan_input.faa`

**Output:**
- `results/{sample}/interproscan_results.tsv`

**Tools:** InterProScan

**Notes:**
- Requires separate InterProScan installation; set path in `config.yaml`
- Creates an empty file if disabled or unavailable
- Excludes the PRINTS database for speed
- Results are used by the locus viewer, domain annotations tab, co-transcription tab, and binomial analysis

### Step 9: Co-transcription Analysis

**Rules:** `map_blast_hits_to_contig_orfs`, `identify_downstream_orfs`, `extract_cotranscribed_sequences`
**Location:** `workflow/rules/cotranscription.smk`

Identifies genes immediately downstream of query hits on the same strand that are likely co-transcribed.

**Process:**
1. Map original BLAST hit coordinates to contig ORFs (coordinate matching, ±3 bp tolerance)
2. For each mapped hit, find the immediately adjacent ORF on the same strand
3. Retain pairs where the intergenic gap is within `distance_threshold`
4. Extract protein sequences of co-transcribed ORFs

**Output:**
- `results/{sample}/cotranscription/hit_to_contig_mapping.tsv` - BLAST hit → contig ORF mapping
- `results/{sample}/cotranscription/downstream_orf_ids.txt` - IDs of co-transcribed ORFs
- `results/{sample}/cotranscription/cotranscribed_details.txt` - Full pair details
- `results/{sample}/cotranscription/cotranscribed_sequences.faa` - Co-transcribed protein sequences

**Parameters:**
- `distance_threshold`: Maximum intergenic gap in bp (default: 100)

### Step 10: DefenseFinder (Optional)

**Rule:** `run_defensefinder`
**Location:** `workflow/rules/defensefinder.smk`

Detects bacterial defense systems in contig ORFs using DefenseFinder (MacSyFinder).

**Input:**
- `results/{sample}/contig_orfs/all_contigs.faa`

**Output:**
- `results/{sample}/defensefinder/defense_finder_systems.tsv`
- `results/{sample}/defensefinder/defense_finder_genes.tsv`
- `results/{sample}/defensefinder/defense_finder_hmmer.tsv`

**Tools:** DefenseFinder, MacSyFinder

**Notes:**
- Uses "unordered" mode since contigs are extracted fragments, not full chromosomes
- Gracefully handles missing dependencies

## Phase 4: Scoring and Enrichment

### Step 11: Calculate Defense Scores

**Rule:** `calculate_defense_scores`
**Location:** `workflow/rules/defense_score.smk`
**Script:** `defense_score.py`

Scores each co-transcribed downstream gene based on its spatial association with known defense systems. This is run independently of the binomial analysis and does not require it.

**Inputs:**
- `results/{sample}/cotranscription/cotranscribed_details.txt`
- `results/{sample}/defensefinder/defense_finder_genes.tsv`
- `results/{sample}/contig_orfs/all_contigs.faa`
- `results/{sample}/interproscan_results.tsv` (if available)

**Output:**
- `results/{sample}/defense_scores.tsv`

**Scoring logic:**

The defense score quantifies how closely associated a co-transcribed gene is with known defense systems on the same contig. It combines two components:

1. **Proximity**: Distance (in genes) from the co-transcribed ORF to the nearest DefenseFinder-annotated gene on the same contig
2. **Local density**: Fraction of genes within a sliding window that are defense genes

A **low score** indicates that the co-transcribed gene is spatially isolated from known defense systems — these are candidates that would be missed by traditional defense island detection methods. A **high score** indicates the gene sits within or immediately adjacent to a known defense cluster.

Scores are reported in the range 0–1 and visualized as a distribution in the Co-transcription tab of the dashboard.

### Step 12: Binomial Domain Enrichment Analysis (Optional)

**Rules:** `extract_genomes_with_hits`, `translate_genomes_for_binomial`, `clean_proteins_for_binomial`, `run_whole_genome_interproscan`, `run_binomial_analysis`
**Location:** `workflow/rules/binomial_analysis.smk`
**Script:** `workflow/scripts/binomial_analysis.py`

Identifies Pfam domains that are statistically enriched in the extracted contigs (flanking regions around query hits) compared to the background frequency in the full genomes that contained hits.

**Enable in config:**
```yaml
binomial:
  enabled: true
```

**Process:**
1. Collect all genomes that produced at least one BLAST hit
2. Combine and translate all proteins from those genomes (Prodigal)
3. Run InterProScan on the whole-genome protein set (or provide a precomputed TSV via `whole_genome_interproscan` in config)
4. Filter both contig and whole-genome annotations to Pfam domains only
5. For each Pfam domain observed in the contigs:
   - Count occurrences in contigs (*k*)
   - Count total occurrences in whole genomes (*n*)
   - Estimate background rate (*p* = whole-genome Pfam frequency)
   - Compute one-sided binomial test: P(X ≥ k | n, p)
6. Apply Benjamini-Hochberg multiple testing correction

**Output:**
- `results/{sample}/binomial/genomes_with_hits.fna` - Combined genome sequences
- `results/{sample}/binomial/genomes_with_hits.faa` - Translated proteins
- `results/{sample}/binomial/interproscan_wholegenomes.tsv` - Whole-genome domain annotations
- `results/{sample}/BinomialAnalysis.csv` - Enrichment results

**BinomialAnalysis.csv columns:**
| Column | Description |
|--------|-------------|
| interpro_id | Pfam domain accession |
| description | Domain description |
| contig_count | Occurrences in extracted contigs |
| genome_count | Occurrences in whole genomes |
| p_value | Binomial test p-value |
| adjusted_p_value | Benjamini-Hochberg corrected p-value |

Results are displayed in both the Co-transcription tab (focused on co-transcribed genes) and the Domain Annotations tab (all contig ORFs) of the dashboard.

### Step 13: Calculate Metrics

**Rule:** `calculate_metrics`
**Location:** `workflow/rules/metrics.smk`

Calculates summary statistics for the pipeline run.

**Output:**
- `results/{sample}/metrics/pipeline_metrics.json`

**Metrics:**
```json
{
  "total_input_bases": 222146374,
  "total_input_mb": 222.15,
  "blast_hits": 370,
  "hits_per_mb": 1.666,
  "defense_systems": 0,
  "contigs_analyzed": 370,
  "discovery_rate_per_contig": 0.0,
  "total_genomes": 222,
  "discovery_rate_per_genome": 0.0,
  "query_file": "resources/EFB0058.fa",
  "query_name": "EFB_0058",
  "sample_name": "efb0058_hits_test"
}
```

## Output Directory Structure

```
results/{sample_name}/
├── downloaded_genomes/          # Raw NCBI downloads
│   └── ncbi_dataset/data/
├── genomes/                     # Organized .fna files
│   ├── genome1.fna
│   └── genome2.fna
├── translated/                  # Protein sequences
│   ├── genome1.faa
│   └── genome2.faa
├── diamond_dbs/                 # BLAST databases
│   ├── genome1.dmnd
│   └── genome2.dmnd
├── blast_results/
│   ├── genome1_blast.txt
│   ├── genome2_blast.txt
│   └── master_blast.txt         # Combined results
├── bed_files/                   # Extraction coordinates
├── contigs/                     # Extracted flanking regions
├── all_contigs.fna              # Combined contigs
├── all_contigs_filtered.fna     # Non-empty contigs
├── contig_orfs/
│   ├── all_contigs.faa          # ORFs from contigs
│   └── interproscan_input.faa   # Cleaned for InterProScan
├── interproscan_results.tsv     # Domain annotations (contig ORFs)
├── cotranscription/
│   ├── hit_to_contig_mapping.tsv
│   ├── downstream_orf_ids.txt
│   ├── cotranscribed_details.txt
│   └── cotranscribed_sequences.faa
├── defensefinder/
│   ├── defense_finder_systems.tsv
│   ├── defense_finder_genes.tsv
│   └── defense_finder_hmmer.tsv
├── defense_scores.tsv           # Defense association scores for co-transcribed genes
├── binomial/                    # (if binomial.enabled: true)
│   ├── genomes_with_hits.fna
│   ├── genomes_with_hits.faa
│   └── interproscan_wholegenomes.tsv
├── BinomialAnalysis.csv         # Pfam enrichment results (if binomial enabled)
├── metrics/
│   ├── total_genome_size.txt
│   └── pipeline_metrics.json
└── genome_manifest.txt          # List of processed genomes
```

## Running Individual Steps

You can run specific pipeline steps by naming their output files:

```bash
# Just download genomes
pixi run snakemake results/{sample}/downloaded_genomes/ncbi_dataset/data --cores 4

# Run up to BLAST
pixi run snakemake results/{sample}/blast_results/master_blast.txt --cores 8

# Run DefenseFinder specifically
pixi run snakemake results/{sample}/defensefinder/defense_finder_systems.tsv --cores 4

# Run defense scoring (requires co-transcription and DefenseFinder to be complete)
pixi run snakemake results/{sample}/defense_scores.tsv --cores 4

# Run binomial analysis (requires interproscan and co-transcription; binomial must be enabled)
pixi run snakemake results/{sample}/BinomialAnalysis.csv --cores 8
```

## Parallelization

The pipeline automatically parallelizes:

- Genome translation (per genome)
- Database building (per genome)
- BLAST searches (per genome)
- Contig extraction (per genome)

Use `--cores N` to control parallelism:

```bash
pixi run snakemake --cores 16  # Use 16 cores
```

For estimated runtimes by genome count, core count, and configuration, see the [[Advanced-Usage#runtime-estimates|Performance Tuning]] section of Advanced Usage.

## Next Steps

- [[Dashboard-Guide]] - Explore results interactively
- [[Troubleshooting]] - Common pipeline issues
