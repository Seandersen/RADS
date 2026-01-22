# Configuration Reference

RADS is configured through the `config/config.yaml` file. This page documents all available options.

## Configuration File Location

```
RADS/
├── config/
│   ├── config.yaml     # Main configuration (edit this)
│   └── test.yaml       # Test configuration example
```

## Complete Configuration Reference

```yaml
# ===========================
# RADS Pipeline Configuration
# ===========================

# Sample name (used for output directories)
# Results will be saved to results/{sample_name}/
sample_name: "my_analysis"

# ===========================
# Input Sources
# ===========================

# Option 1: Download genomes from NCBI
download:
  enabled: true                      # Set true to download from NCBI
  accession_file: "accessions.txt"   # File with nucleotide accessions (one per line)

  # Alternative: Taxon-based download (used if accession_file is not set)
  taxon: "Bacillota"                 # NCBI taxon query
  source: "refseq"                   # refseq or genbank
  assembly_level: "complete"         # complete, chromosome, scaffold, contig
  max_genomes: 0                     # Limit number of genomes (0 = no limit)

# Option 2: Use local genomes (when download.enabled = false)
# Should point to directory containing genome directories with .fna files
genomes_path: "path/to/genomes/ncbi_dataset/data"

# ===========================
# Query Settings
# ===========================

# Query sequence file (amino acid FASTA format)
query_file: "resources/EFB0058.fa"

# Flanking region to extract around each hit
upstream_nt: 5000      # Nucleotides upstream of hit
downstream_nt: 5000    # Nucleotides downstream of hit

# ===========================
# Analysis Settings
# ===========================

# Prodigal ORF prediction
prodigal:
  mode: "meta"    # "meta" for metagenomes/mixed samples
                  # "single" for single genome (requires training)

# Diamond BLAST settings
diamond:
  identity: 30           # Minimum percent identity threshold
  threads: 8             # Number of parallel threads
  max_target_seqs: 0     # Max hits per query (0 = no limit)

# Co-transcription analysis
cotranscription:
  identity_threshold: 95    # % identity for full-length query matching
  distance_threshold: 100   # Max bp distance for co-transcription

# ===========================
# Optional Features
# ===========================

# InterProScan domain annotation
# Requires separate InterProScan installation
interproscan:
  enabled: true
  path: "/path/to/interproscan-5.69-101.0/interproscan.sh"

# DefenseFinder defense system detection
defensefinder:
  enabled: true
  db_type: "unordered"    # Use "unordered" for contigs
  coverage: 0.4           # Minimum coverage threshold
  workers: 4              # Number of parallel workers
```

## Configuration Options Explained

### Sample Name

```yaml
sample_name: "my_analysis"
```

- Used to name output directories
- Results saved to `results/{sample_name}/`
- Use descriptive names without spaces

### Genome Sources

#### Option A: Download from NCBI (Recommended)

```yaml
download:
  enabled: true
  accession_file: "my_accessions.txt"
```

The accession file should contain one NCBI nucleotide accession per line:
```
NC_000913.3
NZ_CP088158.1
CP016471.1
```

#### Option B: Taxon-Based Download

```yaml
download:
  enabled: true
  taxon: "Bacillota"
  source: "refseq"
  assembly_level: "complete"
  max_genomes: 100
```

#### Option C: Local Genomes

```yaml
download:
  enabled: false
genomes_path: "/path/to/my_genomes/ncbi_dataset/data"
```

Expected structure:
```
genomes_path/
├── GCF_000001/
│   └── GCF_000001.fna
├── GCF_000002/
│   └── GCF_000002.fna
```

### Query Sequence

```yaml
query_file: "resources/my_query.fa"
```

- Must be amino acid FASTA format
- Can contain single or multiple sequences
- Example format:
```fasta
>EF_B0058
MKTQPIKVN...
```

### Flanking Region Size

```yaml
upstream_nt: 5000
downstream_nt: 5000
```

- Controls how much sequence to extract around hits
- Larger values capture more context but increase file sizes
- Recommended: 5000-10000 bp for defense system analysis

### BLAST Settings

```yaml
diamond:
  identity: 30
  threads: 8
  max_target_seqs: 0
```

| Setting | Description | Recommended |
|---------|-------------|-------------|
| `identity` | Minimum % identity | 30-40 for distant homologs |
| `threads` | Parallel threads | Match CPU cores |
| `max_target_seqs` | Hits per query | 0 (unlimited) |

### Co-transcription Settings

```yaml
cotranscription:
  identity_threshold: 95
  distance_threshold: 100
```

- `distance_threshold`: Maximum intergenic distance (bp) to consider genes co-transcribed
- Lower values are more stringent

### InterProScan

```yaml
interproscan:
  enabled: true
  path: "/absolute/path/to/interproscan.sh"
```

- Set `enabled: false` if InterProScan is not installed
- Path must be absolute, not relative

### DefenseFinder

```yaml
defensefinder:
  enabled: true
  db_type: "unordered"
  coverage: 0.4
  workers: 4
```

- `db_type`: Use "unordered" since extracted contigs aren't in genomic order
- `coverage`: Minimum model coverage (0.4 = 40%)
- Set `enabled: false` if DefenseFinder has issues

## Example Configurations

### Minimal Configuration

```yaml
sample_name: "quick_test"
download:
  enabled: true
  accession_file: "test_accessions.txt"
query_file: "resources/EFB0058.fa"
upstream_nt: 5000
downstream_nt: 5000
diamond:
  identity: 30
  threads: 4
interproscan:
  enabled: false
defensefinder:
  enabled: false
```

### Full Analysis

```yaml
sample_name: "full_analysis"
download:
  enabled: true
  accession_file: "all_accessions.txt"
query_file: "resources/EFB0058.fa"
upstream_nt: 10000
downstream_nt: 10000
prodigal:
  mode: "meta"
diamond:
  identity: 25
  threads: 16
cotranscription:
  identity_threshold: 95
  distance_threshold: 150
interproscan:
  enabled: true
  path: "/opt/interproscan/interproscan.sh"
defensefinder:
  enabled: true
  db_type: "unordered"
  coverage: 0.4
  workers: 8
```

## Using Multiple Configurations

You can create different config files for different analyses:

```bash
# Run with specific config
pixi run snakemake --configfile config/my_special_config.yaml --cores 8

# Override specific values
pixi run snakemake --config sample_name="test123" --cores 8
```

## Next Steps

- [[Pipeline-Overview]] - Learn what each step does
- [[Troubleshooting]] - Common configuration issues
