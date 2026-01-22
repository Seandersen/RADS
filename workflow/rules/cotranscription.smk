# Co-transcription analysis
# =========================
# Uses coordinate-based mapping to find ORFs downstream of BLAST hits

import re
from pathlib import Path


rule map_blast_hits_to_contig_orfs:
    """
    Map original BLAST hits to their corresponding ORFs in the extracted contigs
    by matching genomic coordinates.
    """
    input:
        master_blast = f"results/{SAMPLE}/blast_results/master_blast.txt",
        contig_orfs = f"results/{SAMPLE}/contig_orfs/all_contigs.faa",
        translated_dir = f"results/{SAMPLE}/translated"
    output:
        mapping = f"results/{SAMPLE}/cotranscription/hit_to_contig_mapping.tsv"
    log:
        f"logs/{SAMPLE}/map_hits_to_contigs.log"
    run:
        import os

        with open(log[0], 'w') as logfile:
            logfile.write("Starting coordinate-based mapping...\n")

            # Parse master BLAST results
            blast_hits = {}  # genome -> list of (subject_id, orf_num)
            with open(input.master_blast) as f:
                next(f)  # skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 7:
                        subject_id = parts[1]  # e.g., NZ_CP088158.1_2443
                        genome = parts[6]
                        # Extract ORF number from subject_id
                        match = re.match(r'(.+)_(\d+)$', subject_id)
                        if match:
                            blast_hits.setdefault(genome, []).append({
                                'subject_id': subject_id,
                                'orf_num': int(match.group(2))
                            })

            logfile.write(f"Found BLAST hits in {len(blast_hits)} genomes\n")

            # Get coordinates for each BLAST hit from translated genome files
            hit_coords = {}  # subject_id -> (start, stop, strand)
            for genome, hits in blast_hits.items():
                translated_file = os.path.join(input.translated_dir, f"{genome}.faa")
                if os.path.exists(translated_file):
                    orf_nums_needed = {h['orf_num'] for h in hits}
                    with open(translated_file) as f:
                        for line in f:
                            if line.startswith('>'):
                                # Parse header: >NZ_CP088158.1_2443 # 2649704 # 2650306 # -1 # ...
                                parts = line[1:].split(' # ')
                                if len(parts) >= 4:
                                    orf_id = parts[0].strip()
                                    match = re.match(r'(.+)_(\d+)$', orf_id)
                                    if match and int(match.group(2)) in orf_nums_needed:
                                        start = int(parts[1])
                                        stop = int(parts[2])
                                        strand = int(parts[3])
                                        hit_coords[orf_id] = (start, stop, strand)

            logfile.write(f"Retrieved coordinates for {len(hit_coords)} BLAST hits\n")

            # Parse contig ORF headers to get their coordinates
            # Format: >NZ_CP088158.1_2644705-2655306:._9 # 5000 # 5602 # -1 # ...
            contig_orfs_info = {}  # contig_id -> list of (orf_id, rel_start, rel_stop, strand)
            contig_regions = {}  # contig_id -> (genome_start, genome_end)

            with open(input.contig_orfs) as f:
                for line in f:
                    if line.startswith('>'):
                        parts = line[1:].split(' # ')
                        if len(parts) >= 4:
                            full_orf_id = parts[0].strip()
                            rel_start = int(parts[1])
                            rel_stop = int(parts[2])
                            strand = int(parts[3])

                            # Parse contig region from ORF ID
                            # Format: NZ_CP088158.1_2644705-2655306:._9
                            match = re.match(r'(.+?)_(\d+)-(\d+):\._(\d+)$', full_orf_id)
                            if match:
                                genome_id = match.group(1)
                                region_start = int(match.group(2))
                                region_end = int(match.group(3))
                                orf_num = int(match.group(4))

                                contig_id = f"{genome_id}_{region_start}-{region_end}"
                                contig_regions[contig_id] = (region_start, region_end, genome_id)
                                contig_orfs_info.setdefault(contig_id, []).append({
                                    'full_id': full_orf_id,
                                    'orf_num': orf_num,
                                    'rel_start': rel_start,
                                    'rel_stop': rel_stop,
                                    'strand': strand,
                                    'genome_start': region_start + rel_start,
                                    'genome_stop': region_start + rel_stop
                                })

            logfile.write(f"Parsed {len(contig_orfs_info)} contigs with ORF information\n")

            # Map BLAST hits to contig ORFs by matching coordinates
            with open(output.mapping, 'w') as out:
                out.write("blast_hit_id\tcontig_orf_id\tgenome\tstrand\tgenome_start\tgenome_stop\n")

                for blast_hit_id, (hit_start, hit_stop, hit_strand) in hit_coords.items():
                    # Find genome from blast_hit_id
                    match = re.match(r'(.+)_\d+$', blast_hit_id)
                    if not match:
                        continue
                    genome_id = match.group(1)

                    # Find matching contig ORF
                    found = False
                    for contig_id, orfs in contig_orfs_info.items():
                        if contig_id.startswith(genome_id):
                            region_start = contig_regions[contig_id][0]
                            for orf in orfs:
                                # Check if coordinates match (with small tolerance)
                                if abs(orf['genome_start'] - hit_start) <= 3 and \
                                   abs(orf['genome_stop'] - hit_stop) <= 3 and \
                                   orf['strand'] == hit_strand:
                                    out.write(f"{blast_hit_id}\t{orf['full_id']}\t{genome_id}\t{hit_strand}\t{hit_start}\t{hit_stop}\n")
                                    logfile.write(f"Mapped {blast_hit_id} -> {orf['full_id']}\n")
                                    found = True
                                    break
                        if found:
                            break

                    if not found:
                        logfile.write(f"No contig match for {blast_hit_id} at {hit_start}-{hit_stop}\n")


rule identify_downstream_orfs:
    """
    Find ORFs immediately downstream of mapped BLAST hits based on strand.
    For plus strand: downstream = next higher ORF number
    For minus strand: downstream = next lower ORF number (towards 3' end)
    """
    input:
        mapping = f"results/{SAMPLE}/cotranscription/hit_to_contig_mapping.tsv",
        contig_orfs = f"results/{SAMPLE}/contig_orfs/all_contigs.faa"
    output:
        downstream = f"results/{SAMPLE}/cotranscription/downstream_orf_ids.txt",
        details = f"results/{SAMPLE}/cotranscription/cotranscribed_details.txt"
    params:
        distance = config["cotranscription"]["distance_threshold"]
    log:
        f"logs/{SAMPLE}/identify_downstream.log"
    run:
        with open(log[0], 'w') as logfile:
            # Parse all contig ORF info
            contig_orfs_by_region = {}  # contig_region -> sorted list of orf info

            with open(input.contig_orfs) as f:
                for line in f:
                    if line.startswith('>'):
                        parts = line[1:].split(' # ')
                        if len(parts) >= 4:
                            full_orf_id = parts[0].strip()
                            rel_start = int(parts[1])
                            rel_stop = int(parts[2])
                            strand = int(parts[3])

                            match = re.match(r'(.+?)_(\d+)-(\d+):\._(\d+)$', full_orf_id)
                            if match:
                                contig_region = f"{match.group(1)}_{match.group(2)}-{match.group(3)}"
                                orf_num = int(match.group(4))

                                contig_orfs_by_region.setdefault(contig_region, []).append({
                                    'full_id': full_orf_id,
                                    'orf_num': orf_num,
                                    'rel_start': rel_start,
                                    'rel_stop': rel_stop,
                                    'strand': strand
                                })

            # Sort ORFs by position within each contig
            for region in contig_orfs_by_region:
                contig_orfs_by_region[region].sort(key=lambda x: x['rel_start'])

            # Read mapping and find downstream ORFs
            downstream_orfs = []

            with open(input.mapping) as f:
                next(f)  # skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        blast_hit_id = parts[0]
                        contig_orf_id = parts[1]
                        hit_strand = int(parts[3])

                        # Parse contig region from ORF ID
                        match = re.match(r'(.+?)_(\d+)-(\d+):\._(\d+)$', contig_orf_id)
                        if not match:
                            continue

                        contig_region = f"{match.group(1)}_{match.group(2)}-{match.group(3)}"
                        hit_orf_num = int(match.group(4))

                        if contig_region not in contig_orfs_by_region:
                            continue

                        orfs = contig_orfs_by_region[contig_region]

                        # Find the hit ORF and determine downstream
                        hit_orf = None
                        hit_idx = None
                        for i, orf in enumerate(orfs):
                            if orf['orf_num'] == hit_orf_num:
                                hit_orf = orf
                                hit_idx = i
                                break

                        if hit_orf is None:
                            logfile.write(f"Could not find hit ORF {hit_orf_num} in {contig_region}\n")
                            continue

                        # Downstream depends on strand
                        # Plus strand (+1): downstream is next ORF (higher position)
                        # Minus strand (-1): downstream is previous ORF (lower position)
                        if hit_strand > 0:
                            # Plus strand - downstream is next in list
                            if hit_idx + 1 < len(orfs):
                                downstream_orf = orfs[hit_idx + 1]
                                gap = downstream_orf['rel_start'] - hit_orf['rel_stop']
                        else:
                            # Minus strand - downstream is previous in list
                            if hit_idx - 1 >= 0:
                                downstream_orf = orfs[hit_idx - 1]
                                gap = hit_orf['rel_start'] - downstream_orf['rel_stop']

                        # Check distance threshold and same strand
                        if 'downstream_orf' in dir() and downstream_orf is not None:
                            if downstream_orf['strand'] == hit_strand and abs(gap) <= params.distance:
                                downstream_orfs.append({
                                    'hit_id': blast_hit_id,
                                    'hit_contig_orf': contig_orf_id,
                                    'downstream_orf': downstream_orf['full_id'],
                                    'strand': hit_strand,
                                    'gap': gap
                                })
                                logfile.write(f"{contig_orf_id} (strand {hit_strand}) -> {downstream_orf['full_id']} (gap: {gap}bp)\n")

                        # Reset for next iteration
                        downstream_orf = None

            # Write outputs
            with open(output.downstream, 'w') as out:
                for item in downstream_orfs:
                    out.write(f"{item['downstream_orf']}\n")

            with open(output.details, 'w') as out:
                out.write("blast_hit_id\thit_contig_orf\tdownstream_orf\tstrand\tgap_bp\n")
                for item in downstream_orfs:
                    out.write(f"{item['hit_id']}\t{item['hit_contig_orf']}\t{item['downstream_orf']}\t{item['strand']}\t{item['gap']}\n")

            logfile.write(f"\nFound {len(downstream_orfs)} co-transcribed downstream ORFs\n")


rule extract_cotranscribed_sequences:
    """
    Extract protein sequences for co-transcribed ORFs.
    """
    input:
        downstream = f"results/{SAMPLE}/cotranscription/downstream_orf_ids.txt",
        proteins = f"results/{SAMPLE}/contig_orfs/all_contigs.faa"
    output:
        sequences = f"results/{SAMPLE}/cotranscription/cotranscribed_sequences.faa"
    log:
        f"logs/{SAMPLE}/extract_cotx_seqs.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        if [ -s {input.downstream} ]; then
            seqkit grep -f {input.downstream} {input.proteins} > {output.sequences} 2>> {log}
            echo "Extracted $(grep -c '>' {output.sequences} || echo 0) sequences" >> {log}
        else
            touch {output.sequences}
            echo "No co-transcribed ORFs found" >> {log}
        fi
        """
