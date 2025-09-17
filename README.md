# RADS
Recombinase Associated Defense Search: a package for extracting loci surrounding a query


RADS was developed by Shelby E Andersen in collaboration with Joshua M Kirsch, Jay R Hesselberth, and Breck A Duerkop as a command line tool for extractng loci surrounding the EF_B0058 serine recombinase for the purpose of identifying surrounding ORFS to find antiphage defense systems. RADS may now be used to extract any sized locus surrounding an ORF of interest.

# Installation

```{bash}
##Create and activate a conda environment for RADS
conda create -n RADS
conda activate RADS

##Install diamond blast into the RADS conda environment
##If using a non-silicon chip Mac or another linux distribution
conda install -c bioconda -c conda-forge diamond
##If using an apple silicon chip Mac
brew install diamond

##Install seqkit into the RADS conda environment
conda install seqkit

##Install prodigal into the RADS conda environment
##If using a Mac
brew install prodigal
##If using another linux distribution
conda install prodigal

##Pull RADS from GitHub
git clone https://github.com/Seandersen/RADS.git

##Navigate into the RADS directory
cd RADS/

##Install interproscan into the RADS directory
mkdir my_interproscan
cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz.md5
tar -pxvzf interproscan-5.69-101.0-*-bit.tar.gz
##If errors with "Cannot utime: operation not permitted", run:
#sudo tar -pxvzf interproscan-5.69-101.0-*-bit.tar.gz

python3 setup.py -f interproscan.properties
```

If you wish to test RADS, test genomic information has been included in the download.
```{bash}
## conda activate RADS
./RADS.sh -g testgenomes/ncbi_dataset/data
```

# Quick Start
Due to its original intended purpose, RADS will default with 5000 nucleotides up- and down-stream of the queried ORF. EFB0058.fa, included in the download, is the default query. 

```{bash}
##conda activate RADS
./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. (Required. Often ncbi_dataset/data/)] -n [sample name to append to all file and directory names (Optional)]
```

For help in the command line...
```{bash}
##conda activate RADS
./RADS.sh -h
```

To change the defaults, simply supply the corresponding options.
```{bash}
-u	integer for number of nucleotides upstream of ORF to extract. Default=5000 if not specified
-d	integer for number of nucleotides downstream of ORF to extract. Default=5000 if not specified
-q	amino acid fasta (.faa/.fa) file for ORF(s) to use as query. Default=EFB0058.fa if not specified
-n  sample name to be added to all output files and directories
```

If you want to run the categorization performed in [pub link], utilize the RMarkdown file RADSCat.rmd

# Output File Information
RADS will generate eight output directories with various data, and one output text file in the original directory. Information on those directories and files can be found here.

| File | Data Contained |
| --- | --- |
| genomes/ | genomic .fna files parsed from the input genomes directory supplied by -g |
| genomestranslated/ | genomes translated into amino acid fasta (.faa) files by prodigal |
| diamonddbs/ | Diamond blast databases generated from amino acid fasta files to be used to blast against |
| blast_results_30 | 30% amino acid identity hits of query supplied by -q or EF_B0058 (default) |
| master.txt | list of all blast results concatenated from all genomes |
| EF_B0058_ORFs/ | list of ORF IDs and coordinates containing query (EFB0058 default) |
| EF_B0058_flanks/ | files of flanks coordinates. Size defaults to 5000nt up- and down-stream. Can be changed by providing integers with options -u and -d |
| bedfiles/ | bed files used by seqkit for contig extraction |
| allcontigsconcatenated_.fna | all RADS output contigs as nucleic acids in a single .fna file |
| allcontigsconcatenated_.txt | all RADS output contigs' ORF coordinates |
| allcontigsconcatenated_.faa | all RADS output contigs translated to protein sequences |
|interproscaninput_.faa | allcontigsconcatenated_.faa without * for use by interproscan |
| interproscaninut_.txt | output of interproscan - all available domain data for proteins in RADS contigs |

# Common Issues and How to Fix Them
| Issue | Solution |
| --- | --- |
| Running frozen during BLASTing | Cancel your run with ctrl + c, comment out chunks before the BLAST process (lines X-XX) (to avoid the laborious re-translating genomes step), re-start your run with the same sample name/query/etc. RADS was built to avoid re-BLASTing files that have already been completed |

