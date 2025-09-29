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
./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. (Required. Often ncbi_dataset/data/)] -n [sample name to appended to all file and directory names. (Optional)] -s [step to run individually. (Optional)]
```

For help in the command line...
```{bash}
##conda activate RADS
./RADS.sh -h
```

To change the defaults, simply supply the corresponding options.
```{bash}
-h|help		  print help message and exit
-u|upint	  integer for number of nucleotides upstream of ORF to extract. Default=5000 if not specified
-d|downint	integer for number of nucleotides downstream of ORF to extract. Default=5000 if not specified
-q|query	  amino acid fasta (.faa/.fa) file for ORF(s) to use as query. Default=EFB0058.fa if not specified
-n|name		  sample name to be added to all output files and directories
-s|step		  run an individual step of the pipeline. Options:
         			  1|mvgenomes --> move genomes into genomes_
         			  2|translate --> translate genomic sequences into amino acid sequences
				        3|makedbs --> make diamond databases from amino acid sequences
				        4|blast --> blast for query in diamond databases
				        5|extractcontigs --> extract contigs of specified size (default 10kb) around query
				        6|contigorfprocessing --> search for ORFs and run InterProScan on contigs around query
				        7|cotranscription --> search for putatively cotranscribed genes around query
```

If you want to run the binomial analysis performed in [pub link], utilize the RMarkdown file RADS.rmd

# Output File Information
RADS will generate eight output directories with various data, and one output text file in the original directory. Information on those directories and files can be found here.

|Step |File | Data Contained |
| --- | --- | ---|
| 1: mvgenomes | genomes_${samplename}/ | genomic .fna files parsed from the input genomes directory supplied by -g |
| 2:translate | genomes_translated_${samplename}/ | genomes translated into amino acid fasta (.faa) files by prodigal |
| 3: makedbs | diamonddbs_${samplename}/ | Diamond blast databases generated from amino acid fasta files to be used to blast against |
| 4:blast | blast_results_30_${samplename} | 30% amino acid identity hits of query supplied by -q or EF_B0058 (default) |
| 4:blast | master${samplename}.txt | list of all blast results concatenated from all genomes |
| 5:extractcontigs | EFB0058_ORFs_${samplename}/ | list of ORF IDs and coordinates containing query (EFB0058 default) |
| 5:extractcontigs | EFB0058_flanks_${samplename}/ | files of flanks coordinates. Size defaults to 5000nt up- and down-stream. Can be changed by providing integers with options -u and -d |
| 5:extractcontigs | bedfiles_${samplename}/ | bed files used by seqkit for contig extraction |
| 5:extractcontigs | allcontigsconcatenated_${samplename}.fna | all RADS output contigs as nucleic acids in a single .fna file |
| 6:contigorfprocessing | allcontigsconcatenated_${samplename}.txt | all RADS output contigs' ORF coordinates |
| 6:contigorfprocessing | allcontigsconcatenated_${samplename}.faa | all RADS output contigs translated to protein sequences |
| 6:contigorfprocessing |interproscaninput_${samplename}.faa | allcontigsconcatenated_.faa without * for use by interproscan |
| 6:contigorfprocessing | interproscaninut_${samplename}.txt | output of interproscan - all available domain data for proteins in RADS contigs |
| 7:cotranscription | fulllengthqueryORFS.txt | query ORFs with the flank coordinates in the IDs |
| 7:cotranscription | EFB0058_cotxORFS_${samplename}/ | contains file manipulations of obtaining ORFs downstream of query |
| 7:cotranscription | EFB0058_cotxORFS_$samplename/finallists/ | contains list of likely cotranscribed ORFs for each input genome |

# Common Issues and How to Fix Them
| Issue | Solution |
| --- | --- |
| Running frozen during BLASTing | This is typically caused by running out of RAM due to a large number of files being searched. Cancel your run with ctrl + c, Run again isolating steps of the pipeline using -s starting with [-s blast] or [-s 4]. RADS is incorporated with an argument that ensures completed files will be ignored.|
| Errors running individual steps | Make sure you have provided the properly named files for all parts of the pipeline upstream. Check source code or File Info (above) to ensure you have properly provided and named all necessary files |

