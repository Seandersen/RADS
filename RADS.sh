upint=5000
downint=5000
query=EFB0058.fa
genomespath=""
samplename=""

while getopts ":u:d:q:g:n:h" opt; do
	case $opt in
		u)
			upint="$OPTARG"
			;;
		d)
			downint="$OPTARG"
			;;
		q)
			query="$OPTARG"
			;;
		g)
			genomespath="$OPTARG"
			;;
		h)	echo "usage:  ./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. (Required. Often ncbi_dataset/data/)] -n [sample name to appended to all file and directory names. (Optional)]"
			echo ""
			echo "[REQUIRED OPTIONS]"
			echo "-------------------------------------------------------------"
			echo "-g	path to directory containing directories with genomes information. Each genome should have its own directory containing .fna files"
			echo ""
			echo "[OPTIONAL OPTIONS]"
			echo "------------------------------------------------------------"
			echo "-u	integer for number of nucleotides upstream of ORF to extract. Default=5000 if not specified"
			echo "-d	integer for number of nucleotides downstream of ORF to extract. Default=5000 if not specified"
			echo "-q	amino acid fasta (.faa/.fa) file for ORF(s) to use as query. Default=EFB0058.fa if not specified"
			echo "-h	print help message and exit"
			echo "-n	sample name to be added to all output files and directories"
			exit 1
			;;
		n)
			if [ -n "$OPTARG" ]; then  # Check if an argument was provided to -n
                                samplename="$OPTARG"
                        fi
                        ;;
		\?)
			echo "option -$OPTARG requires an argument. Usage: bash ./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. Often ncbi_dataset/data/)" >&2
			exit 1
			;;
		:)
			echo "option -$OPTARG requires an argument. Usage: bash ./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. Often ncbi_dataset/data/)" >&2
			exit 1
			;;
	esac
done

if [ -z "$genomespath" ]; then
	echo "Error: Option -g (genomes path) is required." >&2
	exit 1
fi

if [ -z "$upint" ]; then
	echo "Error: Option -u (nucleotides upstream) is required as an integer." >&2
	exit 1
fi

if [ -z "$downint" ]; then
	echo "Error: option -d (nucleotides downstream) is required as an integer." >&2
	exit 1
fi
echo "Thank you for using RADS! RADS will now run with the following parameters:"
echo "upstream set to: $upint"
echo "downstream set to: $downint"
echo "query set to: $query"
echo "genomes path set to: $genomespath"
echo "sample name set to: $samplename"

mkdir genomes_${samplename}/
for i in $(ls $genomespath)
do
cat $genomespath/${i} > genomes_${samplename}/${i}.fna
done

mkdir genomes_translated_${samplename}/
cd genomes_${samplename}/
files=( ./*.fna )
for i in "${files[@]}"
do
prodigal -p meta  -i ${i} -o ../genomes_translated_${samplename}/${i}prodigal.txt -a ../genomes_translated_${samplename}/${i}translated.faa
done
cd ../

mkdir diamonddbs_${samplename}/

cd genomes_translated_${samplename}/
files=( ./*.faa )
for i in "${files[@]}"
diamond makedb --in ../00_QuerySequences/firmicutes_genomes_orfs_translations.faa --db diamonddbs_${samplename}/firmicutes_genomes_translated.db --threads 30


mkdir blast_results_30_${samplename}/
cd diamonddbs_${samplename}/
files=( *.dmnd )
for i in "${files[@]}"
do
if [ ! -f ../blast_results_30_$samplename/EFB0058_blast_${i}.txt ]; then
	diamond blastp -d ${i} --query ../$query --threads 40 --out ../blast_results_30_$samplename/EFB0058_blast_${i}.txt --outfmt 6 qseqid sseqid length nident --max-target-seqs 0 --id 30
fi
done
cd ../

for i in $(ls blast_results_30_${samplename}/)
do
cat blast_results_30_${samplename}/${i} >> master$samplename.txt
done

mkdir EFB0058_ORFs_$samplename/
mkdir EFB0058_flanks_$samplename/
mkdir bedfiles_$samplename/
mkdir contigs_$samplename/
mkdir EFB0058_formatted_coordinates_$samplename/
mkdir EFB0058_final_coordinates_$samplename/


for i in $(ls genomes_$samplename/)
do
  if [ -s blast_results_30_$samplename/EFB0058_blast_${i}translated.faa.db.dmnd.txt ]; then
    echo ${i}
    cut -f2 blast_results_30_$samplename/EFB0058_blast_${i}translated.faa.db.dmnd.txt > EFB0058_ORFS_$samplename/${i}_EFB0058_ORFs.txt
    echo "ORF IDs extracted"
    seqkit grep -f EFB0058_ORFS_$samplename/${i}_EFB0058_ORFs.txt genomes_translated_$samplename/${i}translated.faa | seqkit seq -n > EFB0058_ORFS_$samplename/${i}_EFB0058_hits_coordinates.txt
    echo "Coordinates extracted"
    awk 'BEGIN{OFS="\t"} {F="#"} {up=$3; down=$5; if (up>down) print down, up; else print up, down}' EFB0058_ORFS_$samplename/${i}_EFB0058_hits_coordinates.txt > EFB0058_formatted_coordinates_$samplename/${i}_EFB0058_formatted_coordinates.txt
    awk -v upint=$upint -v downint=$downint 'BEGIN{OFS="\t"} {up=$1-upint; down=$2+downint; if (up>0) print up, down; else print 0, down}' EFB0058_formatted_coordinates_$samplename/${i}_EFB0058_formatted_coordinates.txt > EFB0058_final_coordinates_$samplename/${i}_EFB0058_final_coordinates.txt
    echo "Coordinates formatted"
    cut -f 2 blast_results_30_$samplename/EFB0058_blast_${i}translated.faa.db.dmnd.txt | cut -f 1,2 -d "_" > EFB0058_flanks_$samplename/${i}_EFB0058_contigs.txt
    echo "Contig names extracted"
    paste EFB0058_flanks_$samplename/${i}_EFB0058_contigs.txt EFB0058_final_coordinates_$samplename/${i}_EFB0058_final_coordinates.txt > bedfiles_$samplename/${i}.bed
    echo "Bed file created"
    seqkit subseq --bed bedfiles_$samplename/${i}.bed genomes_$samplename/${i} > contigs_$samplename/${i}_EFB0058_contigs.fna
    echo "Contigs created for ${i} with flanks $upint up and $downint down"
  else
    echo "File ${i} is empty, skipping..."
  fi
done

##Translate all RADS contigs and peform domain prediction with interproscan (interproscan must be installed into the RADS directory)
files=(contigs_$samplename/*.fna)
for i in "${files[@]}"
do
cat ${i} >> allcontigsconcatenated_$samplename.fna
done

seqkit seq -m 50 allcontigsconcatenated_$samplename.fna > allcontigsconcatenated_filtered_$samplename.fna
prodigal -p meta -i allcontigsconcatenated_filtered_$samplename.fna -o allcontigsconcatenated_$samplename.txt -a allcontigsconcatenated_$samplename.faa
sed 's/*//g' allcontigsconcatenated_$samplename.faa > interposcaninput_$samplename.faa
bash my_interproscan/interproscan-*/interproscan.sh -i interposcaninput_$samplename.faa


##Generate a list of ORFs that are likely co-transcribed with query ORF
#Pull ORF ID, start, stop, and strand for 0058 homologues
for i in $(ls genomes_$samplename/)
do
awk 'BEGIN{OFS="\t"} {F="#"} {ORF=$1; start=$3; stop=$5; strand=$7; print ORF, start, stop, strand}' EFB0058_ORFS_$samplename/${i}_EFB0058_hits_coordinates.txt > EF0058_ORFS_$samplename/${i}_EFB0058_hits_coordinatesandstrand.txt
done

#Re-blast to get correct full query ID
diamond makedb --in allcontigsconcatenated_$samplename.faa --db diamonddbs_${samplename}/allcontigsconcatenated_$samplename.db --threads 30
diamond blastp -d diamonddbs_${samplename}/allcontigsconcatenated_$samplename.db --query $query --threads 40 --out fulllengthqueryORFS.txt --outfmt 6 qseqid sseqid length nident --max-target-seq 0 --id 30

mkdir EFB0058_cotxORFS_$samplename/

cut -f 2 fulllengthqueryORFS.txt > query_ORFS_long.txt
seqkit grep -f query_ORFS_long.txt allcontigsconcatenated_$samplename.faa | seqkit seq -n > EFB0058_ORFS_$samplename/EFB0058_hits_coordinatesfull.txt
awk 'BEGIN{OFS="\t"} {F="#"} {ORF=$1; start=$3; stop=$5; strand=$7; print ORF, start, stop, strand}' EFB0058_ORFS_$samplename/EFB0058_hits_coordinatesfull.txt > EFB0058_hits_coordinatesandstrand.txt

##generate up/downstream orf number determined by strand for 0058 homologue
while IFS= read -r line; do

	ORF=$(echo $line | awk 'BEGIN{OFS="\t"} {F=" "} {print $1}')
        strand=$(echo $line | awk 'BEGIN{OFS="\t"} {F=" "} {print $4}')

        prefix=$(echo $ORF | sed 's/:\..*//')
        number=$(echo $ORF | sed 's/*:\.._//')

	if ((strand>0)); then
               	new_number=$((number + 1))
       	elif ((strand<0)); then
               	new_number=$((number -1))
       	fi

	downstreamORF=${prefix}:.__${new_number}
	print downstreamORF >> EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs.txt
done < fulllengthqueryORFS.txt

##check whether the downstream ORF is on the same strand and starts within 100bp of query
mkdir EFB0058_cotxORFS_$samplename/finallists/
for i in $(ls genomes_$samplename)
do
seqkit grep -f EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs.txt genomes_translated_$samplename/${i}translated.faa | seqkit seq -n > EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs_fullheader.txt
awk'BEGIN{OFS="\t"} {F="#"} {ORF=$1; start=$3; stop=$5; strand=$7; print ORF, start, stop, strand}' EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs_fullheader.txt > EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs_coordinates.txt
paste EFB0058_ORFS_$samplename/${i}_EFB0058_hits_coordinatesandstrand.txt EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs_coordinates.txt > EFB0058_cotxORFS_$samplename/${i}_EFB0058anddownstream.txt
awk'BEGIN{OFS="\t"} {EFB0058ORF=$1; EFB0058stop=$3; EFB0058strand=$4; downstreamORF=$5; downstreamstart=$6; downstreamstrand=$8; if ((EFB0058strand==downstreamstrand) && (($3-$6) <=100) && (($3-$6)>=-100)) print downstreamORF}' EFB0058_cotxORFS_$samplename/${i}_EFB0058anddownstream.txt > EFB0058_cotxORFS/finallists/${i}_EFB0058_cotxORFs.txt
done
