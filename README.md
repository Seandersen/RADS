# RADS
Recombinase Associated Defense Search: a package for extracting loci surrounding a query


RADS was developed by Shelby E Andersen in collaboration with Joshua M Kirsch, Jay R Hesselberth, and Breck A Duerkop as a command line tool for extractng loci surrounding the EF_B0058 serine recombinase for the purpose of identifying surrounding ORFS to find antiphage defense systems. RADS may now be used to extract any sized locus surrounding an ORF of interest.

# Installation

```{bash}
conda create -n RADS
conda install diamond
git clone [link]
```

# Quick Start
Due to its original intended purpose, RADS will default with 5000 nucleotides up- and down-stream of the queried ORF. EFB0058.fa, included in the download, is the default query. 

```{bash}
./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. (Required. Often ncbi_dataset/data/)]
```

For help in the command line...
```{bash}
./RADS.sh -h
```

To change the defaults, simply supply the corresponding options.
```{bash}
-u	integer for number of nucleotides upstream of ORF to extract. Default=5000 if not specified
-d	integer for number of nucleotides downstream of ORF to extract. Default=5000 if not specified
-q	amino acid fasta (.faa/.fa) file for ORF(s) to use as query. Default=EFB0058.fa if not specified
```
