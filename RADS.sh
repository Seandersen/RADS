#!/bin/sh

upint=5000
downint=5000
query=EFB0058.fa
genomespath=""
samplename=""
run_all=true
step=""

while getopts ":u:d:q:g:n:s:h" opt; do
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
		h)
			echo "usage:  ./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. (Required. Often ncbi_dataset/data/)] -n [sample name to appended to all file and directory names. (Optional)] -s [step to run individually. (Optional)]"
			echo ""
			echo "[REQUIRED OPTIONS]"
			echo "-------------------------------------------------------------"
			echo "-g|--genomespath	path to directory containing directories with genomes information. Each genome should have its own directory containing .fna files"
			echo ""
			echo "[OPTIONAL ARGUMENTS]"
			echo "------------------------------------------------------------"
			echo "-h|--help		print help message and exit"
            		echo "-u|--upint	integer for number of nucleotides upstream of ORF to extract. Default=5000 if not specified"
			echo "-d|--downint	integer for number of nucleotides downstream of ORF to extract. Default=5000 if not specified"
			echo "-q|--query	amino acid fasta (.faa/.fa) file for ORF(s) to use as query. Default=EFB0058.fa if not specified"
			echo "-n|--name		sample name to be added to all output files and directories"
			echo "-s|--step		run an individual step of the pipeline. Options:"
            echo "         		1|mvgenomes --> move genomes into genomes_${samplename}"
            echo "         		2|translate --> translate genomic sequences into amino acid sequences"
			echo "				3|makedbs --> make diamond databases from amino acid sequences"
			echo "				4|blast --> blast for query in diamond databases"
			echo "				5|extractcontigs --> extract contigs of specified size (default 10kb) around query"
			echo "				6|contigorfprocessing --> search for ORFs and run InterProScan on contigs around query"
			echo "				7|cotranscription --> search for putatively cotranscribed genes around query"
            		exit 1
			;;
		n)
			if [ -n "$OPTARG" ]; then  # Check if an argument was provided to -n
                                samplename="$OPTARG"
                        fi
                        ;;
        	s)
            		step="$OPTARG"
            		run_all=false
            		shift 2
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

mvgenomes( ) {
    mkdir genomes_${samplename}/
    for i in $(ls $genomespath)
    do
        cat $genomespath/${i} > genomes_${samplename}/${i}.fna
    done
    echo "genomes located"
}

translate( ){
    mkdir genomes_translated_${samplename}/
    cd genomes_${samplename}/
    files=( ./*.fna )
    for i in "${files[@]}"
    do
        prodigal \
        -p meta \
        -i ${i} \
        -o ../genomes_translated_${samplename}/${i}prodigal.txt \
        -a ../genomes_translated_${samplename}/${i}translated.faa
    done
    cd ../
    echo "genomes translated"
}

makedbs( ){
    mkdir diamonddbs_${samplename}/
    cd genomes_translated_${samplename}/
    files=( ./*.faa )
    for i in "${files[@]}"
    do
        diamond makedb \
        --in ${i} \
        --db ../diamonddbs_${samplename}/${i}.db \
        --threads 30
    done
    cd ../
}

blast( ){
    mkdir blast_results_30_${samplename}/
    cd diamonddbs_${samplename}
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
}

extractcontigs ( ){
    mkdir EFB0058_ORFS_$samplename/
    mkdir EFB0058_flanks_$samplename/
    mkdir bedfiles_$samplename/
    mkdir contigs_$samplename/
    mkdir EFB0058_formatted_coordinates_$samplename/
    mkdir EFB0058_final_coordinates_$samplename/

	max_jobs=20
	count=0
    for i in $(ls genomes_$samplename/)
    do
	{
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
		} &

		((count++))
		 if ((count % max_jobs == 0)); then
		 	wait
		fi
    done
	wait
}

contigorfprocessing( ){
    files=(contigs_$samplename/*.fna)
    for i in "${files[@]}"
    do
        cat ${i} >> allcontigsconcatenated_$samplename.fna
    done

    seqkit seq -m 50 allcontigsconcatenated_$samplename.fna > allcontigsconcatenated_filtered_$samplename.fna
    prodigal \
    -p meta \
    -i allcontigsconcatenated_filtered_$samplename.fna \
    -o allcontigsconcatenated_$samplename.txt \
    -a allcontigsconcatenated_$samplename.faa
    sed 's/*//g' allcontigsconcatenated_$samplename.faa > interposcaninput_$samplename.faa
    bash my_interproscan/interproscan-*/interproscan.sh -i interposcaninput_$samplename.faa
}

cotranscription ( ){
    ## Generate a list of ORFs that are likely co-transcribed with query ORF
    # Pull ORF ID, start, stop, and strand for 0058 homologues
    for i in $(ls genomes_$samplename/); do
        awk 'BEGIN{OFS="\t"; FS="#"} {ORF=$1; start=$3; stop=$5; strand=$7; print ORF, start, stop, strand}' \
        EFB0058_ORFS_$samplename/${i}_EFB0058_hits_coordinates.txt \
        > EFB0058_ORFS_$samplename/${i}_EFB0058_hits_coordinatesandstrand.txt
    done
    echo "Query Homolog ORFs, coordinates, and strand extracted"

    # Re-blast to get correct full query ID
    diamond makedb --in allcontigsconcatenated_$samplename.faa \
                --db diamonddbs_${samplename}/allcontigsconcatenated_$samplename.db \
                --threads 30
    diamond blastp -d diamonddbs_${samplename}/allcontigsconcatenated_$samplename.db \
                --query $query --threads 40 \
                --out fulllengthqueryORFS.txt \
                --outfmt 6 qseqid sseqid length nident --id 95

    mkdir -p EFB0058_cotxORFS_$samplename/

    cut -f 2 fulllengthqueryORFS.txt > query_ORFS_long.txt
    seqkit grep -f query_ORFS_long.txt allcontigsconcatenated_$samplename.faa \
        | seqkit seq -n \
        > EFB0058_ORFS_$samplename/EFB0058_hits_coordinatesfull.txt

    awk 'BEGIN{OFS="\t"; FS="#"} {ORF=$1; start=$2; stop=$3; strand=$4; print ORF, start, stop, strand}' \
        EFB0058_ORFS_$samplename/EFB0058_hits_coordinatesfull.txt \
        > EFB0058_hits_coordinatesandstrand.txt

    ## Generate up/downstream ORF number determined by strand
    mkdir EFB0058_cotxORFS_${samplename}

    while IFS= read -r line; do
        ORF=$(echo $line | awk 'BEGIN{OFS="\t"} {print $1}')
        strand=$(echo $line | awk 'BEGIN{OFS="\t"} {F=" "} {print $4}')

        prefix=$(echo $ORF | sed 's/:\..*//')                                       # part before first :
        number=$(echo $ORF | sed -E 's/.*_([0-9]+)$/\1/')     # extract trailing number

        if (( $strand > 0 )); then
            new_number=$((number + 1))
        elif (( $strand < 0 )); then
            new_number=$((number - 1))
        else
            continue
        fi

        downstreamORF="${prefix}:._${new_number}"
        echo "$downstreamORF" >> EFB0058_cotxORFS_${samplename}/all_downstreamORFs.txt
    done < EFB0058_hits_coordinatesandstrand.txt

    ## Check whether the downstream ORF is on the same strand and starts within 100bp of query
    mkdir -p EFB0058_cotxORFS_$samplename/finallists/

    for i in $(ls genomes_$samplename); do
        seqkit grep -f EFB0058_cotxORFS_${samplename}/all_downstreamORFs.txt allcontigsconcatenated_$samplename.faa \
            | seqkit seq -n \
            > EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs_fullheader.txt

        awk 'BEGIN{OFS="\t"; FS="#"} {ORF=$1; start=$2; stop=$3; strand=$4; print ORF, start, stop, strand}' \
            EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs_fullheader.txt \
            > EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs_coordinates.txt

        paste EFB0058_hits_coordinatesandstrand.txt \
            EFB0058_cotxORFS_$samplename/${i}_EFB0058_downstreamORFs_coordinates.txt \
            > EFB0058_cotxORFS_$samplename/${i}_EFB0058anddownstream.txt

        awk 'BEGIN{OFS="\t"} {
                 EFB0058ORF=$1; EFB0058stop=$3; EFB0058strand=$4;
                 downstreamORF=$5; downstreamstart=$6; downstreamstrand=$8;
                if ((EFB0058strand==downstreamstrand) && (($3-$6) <=100) && (($3-$6)>=-100))
                    print downstreamORF
            }' EFB0058_cotxORFS_$samplename/${i}_EFB0058anddownstream.txt \
            > EFB0058_cotxORFS_$samplename/finallists/${i}_EFB0058_cotxORFs.txt
    done
}

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

if $run_all; then
    echo "Thank you for using RADS! RADS will now run with the following parameters:"
    echo "upstream set to: $upint"
    echo "downstream set to: $downint"
    echo "query set to: $query"
    echo "genomes path set to: $genomespath"
    echo "sample name set to: $samplename"
    echo "running entire pipeline..."
    mvgenomes
    translate
    makedbs
    extractcontigs
    contigorfprocessing
    cotranscription
else
    case "$step" in
        1|mvgenomes)
		    echo "Thank you for using RADS! RADS will now move input genomes to genomes_${samplename}/"
		    mvgenomes ;;
        2|translate)
		    echo "Thank you for using RADS! RADS will now translate input genomes..."
        	translate ;;
	    3|makedbs)
		    echo "Thank you for using RADS! RADS will now make Diamond databases for input genomes..."
		    makedbs ;;
        4|blast)
		    echo "Thank you for using RADS! RADS will now BLAST for your query sequence..."
		    echo "query set to: $query"
		    blast;;
        5|extractcontigs)
		    echo "Thank you for using RADS! RADS will now extract contigs surrounding your query sequence..."
		    echo "upstream set to: $upint"
		    echo "downstream set to: $downint"
		    extractcontigs;;
        6|contigorfprocessing)
		    echo "Thank you for using RADS! RADS will now process your contigs..."
		    contigorfprocessing;;
        7|cotranscription)
		    echo "Thank you for using RADS! RADS will now search for ORFs cotranscribed with your query"
		    cotranscription;;
        *) echo "Error: invalid step argument";;
    esac
fi
