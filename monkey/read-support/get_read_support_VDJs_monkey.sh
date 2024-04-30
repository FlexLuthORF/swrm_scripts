#!/bin/bash
set -e -x

scratch=$PWD

function run_make_ref_masked {
    monkey_mask_ref=/home/zmvanw01/ref/monkey/Child_17thApril2024_bothhap_IG_masked_unique.fasta
    samtools faidx ${monkey_mask_ref}
}

function run_map_ccs_to_pers {
   
    outd=${scratch}/read_support/${sample}
    mkdir -p ${outd}/ccs_to_pers

    # Convert PacBio HiFi reads to FASTA
    samtools view ${ccs_bam} | awk '{ print ">"$1"\n"$10 }' > ${outd}/ccs_to_pers/reads.fasta
    samtools faidx ${outd}/ccs_to_pers/reads.fasta

    # Create personalized reference
    samtools view -F 0x100 -F 0x800 "${assemblies_fasta}" | awk '{print ">"$1"\n"$10}' > "${outd}/ccs_to_pers/contigs.fasta"

    samtools faidx ${outd}/ccs_to_pers/contigs.fasta

    pers_contigs=${outd}/ccs_to_pers/contigs.fasta

    cat ${monkey_mask_ref} ${pers_contigs} > ${outd}/ccs_to_pers/pers_ref.fasta

    # Uniqify the headers in pers_ref.fasta using awk
    #awk '/^>/ {++counter; print ">seq_" counter; next} {print}' ${outd}/ccs_to_pers/pers_ref.fasta > ${outd}/ccs_to_pers/pers_ref_uniqified.fasta
    #mv ${outd}/ccs_to_pers/pers_ref_uniqified.fasta ${outd}/ccs_to_pers/pers_ref.fasta

    #echo "Uniqified FASTA file: ${outd}/ccs_to_pers/pers_ref.fasta"

    samtools faidx ${outd}/ccs_to_pers/pers_ref.fasta

    minimap2 -ax map-hifi --secondary=no -t 10 -L ${outd}/ccs_to_pers/pers_ref.fasta ${outd}/ccs_to_pers/reads.fasta > ${outd}/ccs_to_pers/output.sam
    samtools view -Sbh ${outd}/ccs_to_pers/output.sam > ${outd}/ccs_to_pers/output.bam
    samtools sort -@ 10 ${outd}/ccs_to_pers/output.bam -o ${outd}/ccs_to_pers/output.sorted.bam
    samtools index ${outd}/ccs_to_pers/output.sorted.bam
    rm -f ${outd}/ccs_to_pers/output.sam
}

function get_read_support_vdj3 {
    base_outd="${scratch}/read_support/${sample}/imported_genes"
    bam_file="${scratch}/read_support/${sample}/ccs_to_pers/output.sorted.bam"
    ref="${scratch}/read_support/${sample}/ccs_to_pers/pers_ref.fasta"

    # Create the base output directory if it doesn't exist
    mkdir -p "$base_outd"

    if [ ! -f "${bam_file}.bai" ]; then
        samtools index "$bam_file"
    fi

    for gene_type in "igh" "igk" "igl"
    do
        case "$gene_type" in
            "igh")
                import_out="$igh_digger"
                ;;
            "igk")
                import_out="$igk_digger"
                ;;
            "igl")
                import_out="$igl_digger"
                ;;
        esac
        mkdir -p "${base_outd}/${gene_type}"  # Create the directory for import_out

        if [[ -f "$import_out" ]]; then
            tmp_file="${scratch}/read_support/${sample}/tmp/$(basename "$import_out")_read_support.tmp"
            mkdir -p "$(dirname "$tmp_file")"
            echo "Total_Positions,Average_Coverage,Mismatched_Positions,Matched_Positions,Position_Mismatches,Position_Matches,Percent_Accuracy,Positions_With_At_Least_10x_Coverage,Fully_Spanning_Reads,Fully_Spanning_Reads_100%_Match" > "$tmp_file"
            
            header=$(head -n 1 "$import_out")
            IFS=',' read -ra header_cols <<< "$header"
            for i in "${!header_cols[@]}"; do
                case "${header_cols[$i]}" in
                    "contig") contig_col=$i ;;
                    "start") start_col=$i ;;
                    "end") end_col=$i ;;
                    "ASC_match") gene_col=$i ;;
                esac
            done

            tmp_counts="${tmp_file}_counts"
            > "$tmp_counts"

            tail -n +2 "$import_out" | while IFS=, read -ra line
            do
                contig="${line[$contig_col]}"
                start=$(echo "${line[$start_col]}" | awk '{printf "%.0f", $1}')
                end=$(echo "${line[$end_col]}" | awk '{printf "%.0f", $1}')
                
                gene="${line[$gene_col]}"
                region="${contig}:${start}-${end}"

                contig_filename=$(echo "$contig" | tr '/' '_')
                tmp_bam="${base_outd}/${gene_type}/${contig_filename}_${start}_${end}.bam"
                mkdir -p "$(dirname "$tmp_bam")"
                samtools view -F 0x100 -F 0x800 -b "$bam_file" -o "$tmp_bam" -U "/dev/null" "${contig}:${start}-${end}"
                samtools index "$tmp_bam"

                samtools mpileup -f "$ref" -r "$region" "$tmp_bam" | \
                awk -v total_positions="$((end - start + 1))" -v sample="$sample" \
                'BEGIN {
                    total_reads=0; mismatched_positions=0; matched_positions=0; positions_with_10x=0;
                    mismatch_list=""; match_list="";
                }
                {
                    total_reads += length($5);
                    mismatches = length(gensub(/[.,]/, "", "g", $5));
                    matches = length(gensub(/[^.,]/, "", "g", $5));
                    mismatch_list = (mismatch_list == "" ? mismatches : mismatch_list ":" mismatches);
                    match_list = (match_list == "" ? matches : match_list ":" matches);

                    coverage = length($5);
                    if (coverage >= 10) {
                        positions_with_10x++;
                    }

                    mismatch_rate = mismatches / coverage;
                    match_rate = matches / coverage;

                    if (mismatch_rate > 0.2) {
                        mismatched_positions++;
                    }

                    if (match_rate > 0.8) {
                        matched_positions++;
                    }
                }
                END {
                    avg_reads_per_position = (total_positions > 0) ? total_reads / total_positions : 0;
                    percent_accuracy = (matched_positions / total_positions) * 100;
                    print total_positions, avg_reads_per_position, mismatched_positions, matched_positions, mismatch_list, match_list, percent_accuracy, positions_with_10x;
                }' OFS=',' >> "${tmp_file}_awk_out"

                python /home/zmvanw01/git_repos/swrm_scripts/monkey/read-support/match_subsequences.py "$tmp_bam" "$contig" "$start" "$end" "$gene" "$import_out" > "${tmp_file}_py_out"
                echo "gene is $gene. awk_out is $(cat "${tmp_file}_awk_out"). py_out is $(cat "${tmp_file}_py_out")"
                paste -d ',' "${tmp_file}_awk_out" "${tmp_file}_py_out" >> "$tmp_file"
                rm "${tmp_file}_awk_out" "${tmp_file}_py_out"
                rm "$tmp_bam" "${tmp_bam}.bai"
            done

            if [[ -f "$import_out" && -f "$tmp_file" ]]; then
                combined_file="${scratch}/read_support/${sample}/output/${gene_type}/$(basename "$import_out" .csv)_combined.csv"
                final_output="${scratch}/read_support/${sample}/output/${gene_type}/$(basename "$import_out" .csv)_with_read_support.csv"
                mkdir -p "$(dirname "$combined_file")"
                python -c "import pandas as pd; df1 = pd.read_csv('$import_out'); df2 = pd.read_csv('$tmp_file'); combined = pd.concat([df1, df2], axis=1); combined.to_csv('$final_output', index=False)"
                #rm -f "$tmp_file"
                

               # awk -v sample="$sample" 'BEGIN{FS=OFS=","} {
                #    if (NR == 1) {
                 #       print;
                  #  } else {
                   #     $2 = sample;
                    #    $3 = sample;
                     #   print;
                  #  }
               # }' "$combined_file" > "$final_output"

               # rm -f "$combined_file"
            fi
        fi
    done
}

sample=$1
assemblies_fasta=$2
igh_digger=$3
igk_digger=$4
igl_digger=$5
ccs_bam=$6
monkey_mask_ref=/home/zmvanw01/ref/monkey/Child_17thApril2024_bothhap_IG_masked_unique.fasta
mkdir -p ${scratch}/read_support/${sample}

# Process digger option files using Python
for digger in "$igh_digger" "$igk_digger" "$igl_digger"
do
    new_file_path="${scratch}/read_support/${sample}/$(basename $digger)_modified.csv"
    python -c "import pandas as pd; df = pd.read_csv('$digger'); df['notes'] = df['notes'].str.replace(',', ';', regex=False) if 'notes' in df.columns else df; df.to_csv('$new_file_path', index=False)"
done

# Update the digger variables to point to the new files
igh_digger="${scratch}/read_support/${sample}/$(basename $igh_digger)_modified.csv"
igk_digger="${scratch}/read_support/${sample}/$(basename $igk_digger)_modified.csv"
igl_digger="${scratch}/read_support/${sample}/$(basename $igl_digger)_modified.csv"

#run_make_ref_masked
run_map_ccs_to_pers
get_read_support_vdj3
