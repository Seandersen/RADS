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
cat $genomespath/${i}/*.fna > genomes_${samplename}/${i}.fna
done

mkdir genomes_translated_${samplename}/

for i in $(ls genomes_$samplename/)
do
prodigal -i genomes_$samplename/${i} -o genomes_translated_${samplename}/${i}prodigal.txt -a genomes_translated_${samplename}/${i}translated.faa
done

mkdir diamonddbs_${samplename}/

cd genomes_translated_${samplename}/
for i in $(ls ./*.faa)
do
diamond makedb --in ${i} --db ../diamonddbs_${samplename}/${i}.db --threads 30
done
cd ../

mkdir blast_results_30_${samplename}/
cd diamonddbs_${samplename}/
for i in $(ls ./)
do
diamond blastp -d ${i} --query ../$query --threads 40 --out ../blast_results_30_$samplename/EFB0058_blast_${i}.txt --outfmt 6 qseqid sseqid length nident --max-target-seqs 0 --id 30
done
cd ../

cat blast_results_30_${samplename}/*.txt > master$samplename.txt

while IFS=$'\t' read -r pattern filename; do
	echo "Pattern:$pattern"
	echo "Filename:$filename"
	filepath="genomes/$filename"
	seqkit grep -p $pattern $filepath >> EFB0058_contigs/parsed_$filename
done < master$samplename.txt

mkdir EFB0058_ORFs_$samplename/
mkdir EFB0058_flanks_$samplename/
mkdir bedfiles_$samplename/
mkdir contigs_$samplename/

for i in $(ls genomes_$samplename/)
do
echo ${i}
cut -f2 blast_results_30_$samplename/EFB0058_blast_${i}translated.faa.db.dmnd.txt > EFB0058_ORFS_$samplename/${i}_EFB0058_ORFs.txt
echo “ORF IDs extracted”
seqkit grep -f EFB0058_ORFS_$samplename/${i}_EFB0058_ORFs.txt genomes_translated_$samplename/${i}translated.faa | seqkit seq -n > EFB0058_ORFS_$samplename/${i}_EFB0058_hits_coordinates.txt
echo “Coordinates extracted”
awk 'BEGIN{OFS="\t"} {F="#"} {up=$3; down=$5; if (up>down) print down, up; else print up, down}' EFB0058_ORFS_$samplename/${i}_EFB0058_hits_coordinates.txt > EFB0058_flanks_$samplename/${i}_EFB0058hits_flanks.txt
awk 'BEGIN{OFS="\t"} {up=$1-$upint; down=$2+$downint; if (up>0) print up, down; else print 0, down}' EFB0058_flanks_$samplename/${i}_EFB0058hits_flanks.txt > EFB0058_flanks_$samplename/${i}_EFB0058hits_flanks_forbed.txt
echo “coordinates formatted”
cut -f 2 blast_results_30_$samplename/EFB0058_blast_${i}translated.faa.db.dmnd.txt | cut -f 1 -d "_" > EFB0058_flanks_$samplename/${i}_EFB0058_contigs.txt
echo “Contig names extracted”
paste EFB0058_flanks_$samplename/${i}_EFB0058_contigs.txt EFB0058_flanks_$samplename/${i}_EFB0058hits_flanks_forbed.txt > bedfiles_$samplename/${i}.bed
echo “bed file created”
seqkit subseq --bed bedfiles_$samplename/${i}.bed genomes_$samplename/${i} > contigs_$samplename/${i}_EFB0058_contigs.fna
echo “Contigs created for ${i} with flanks $upint up and $downint down”
done

##Translate all RADS contigs and peform domain prediction with interproscan (interproscan must be installed into the RADS directory)
cat contigs_$samplename/*.fna > allcontigsconcatenated_$samplename.fna
prodigal -i allcontigsconcatenated_$samplename.fna -o allcontigsconcatenated_$samplename.txt -a allcontigsconcatenated_$samplename.faa
sed 's/*//g' allcontigsconcatenated_$samplename.faa > interposcaninput_$samplename.faa
./my_interproscan/interproscan-5.69-101.0/interproscan.sh -i interposcaninput_$samplename.faa
