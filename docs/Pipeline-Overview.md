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
2. Calculate upstream/downstream boundaries
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

Annotates protein domains using InterProScan.

**Input:**
- `results/{sample}/contig_orfs/interproscan_input.faa`

**Output:**
- `results/{sample}/interproscan_results.tsv`

**Tools:** InterProScan

**Notes:**
- Requires separate InterProScan installation
- Creates empty file if disabled or unavailable
- Excludes PRINTS database (faster)

### Step 9: Co-transcription Analysis

**Rules:** `map_blast_hits_to_contig_orfs`, `identify_downstream_orfs`, `extract_cotranscribed_sequences`
**Location:** `workflow/rules/cotranscription.smk`

Identifies genes immediately downstream of query hits that may be co-transcribed.

**Process:**
1. Map original BLAST hits to contig ORFs by matching coordinates
2. For each hit, find the next ORF on the same strand
3. Check if within distance threshold
4. Extract sequences of co-transcribed ORFs

**Output:**
- `results/{sample}/cotranscription/hit_to_contig_mapping.tsv`
- `results/{sample}/cotranscription/downstream_orf_ids.txt`
- `results/{sample}/cotranscription/cotranscribed_details.txt`
- `results/{sample}/cotranscription/cotranscribed_sequences.faa`

**Parameters:**
- `distance_threshold`: Maximum bp between genes (default: 100)

### Step 10: DefenseFinder (Optional)

**Rule:** `run_defensefinder`
**Location:** `workflow/rules/defensefinder.smk`

Detects bacterial defense systems using DefenseFinder.

**Input:**
- `results/{sample}/contig_orfs/all_contigs.faa`

**Output:**
- `results/{sample}/defensefinder/defense_finder_systems.tsv`
- `results/{sample}/defensefinder/defense_finder_genes.tsv`
- `results/{sample}/defensefinder/defense_finder_hmmer.tsv`

**Tools:** DefenseFinder, MacSyFinder

**Notes:**
- Uses "unordered" mode since contigs aren't in genomic order
- Gracefully handles missing dependencies

### Step 11: Calculate Metrics

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
├── interproscan_results.tsv     # Domain annotations
├── cotranscription/
│   ├── hit_to_contig_mapping.tsv
│   ├── downstream_orf_ids.txt
│   ├── cotranscribed_details.txt
│   └── cotranscribed_sequences.faa
├── defensefinder/
│   ├── defense_finder_systems.tsv
│   ├── defense_finder_genes.tsv
│   └── defense_finder_hmmer.tsv
├── metrics/
│   ├── total_genome_size.txt
│   └── pipeline_metrics.json
└── genome_manifest.txt          # List of processed genomes
```

## Running Individual Steps

You can run specific pipeline steps:

```bash
# Just download genomes
pixi run snakemake results/{sample}/downloaded_genomes/ncbi_dataset/data --cores 4

# Run up to BLAST
pixi run snakemake results/{sample}/blast_results/master_blast.txt --cores 8

# Run DefenseFinder specifically
pixi run snakemake results/{sample}/defensefinder/defense_finder_systems.tsv --cores 4
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

## Next Steps

- [[Dashboard-Guide]] - Explore results interactively
- [[Troubleshooting]] - Common pipeline issues
