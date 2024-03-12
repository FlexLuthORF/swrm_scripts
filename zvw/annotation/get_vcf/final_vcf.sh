#!/bin/bash
set -e -x

# Function to genotype SV regions
function genotype_SV_regions {
    local bam_path="$1"
    local SV_regions_1bp="$2"
    local outd="$3"

    while read -r sv_region; do
        sv_name=$(echo "$sv_region" | cut -f4)
        grep -w "$sv_name" "$SV_regions_1bp" > "${sv_name}.bed"

        # Check if an entry for the sample and sv_name already exists
        if grep -q -P "^$sample\t$sv_name\t" "${outd}/SV_genotype_results.txt"; then
            echo "Entry for $sample and $sv_name already exists. Skipping appending step."
            continue
        fi

        mpileup_output=$(samtools mpileup -l "${sv_name}.bed" -f "${reffn}" "$bam_path" | head -1)
        genotype=$(echo "$mpileup_output" | awk '{print $5}')
        asterisk_count=$(echo "$genotype" | tr -cd '*' | wc -c)
        total_count=$(echo "$genotype" | wc -c)

        if [ "$asterisk_count" -eq 0 ]; then
            genotype_label="0/0"
        elif [ "$asterisk_count" -eq "$((total_count-1))" ]; then
            genotype_label="1/1"
        else
            genotype_label="0/1"
        fi

        echo -e "$sample\t$sv_name\t$genotype_label" >> "${outd}/SV_genotype_results.txt"
    done < "$SV_regions_1bp"
}


# Function to process VCF files
function process_vcf {
    local bam_file="$1"
    local sample="$2"
    local outd="$3"
    local reffn="$4"
    local num_threads="$5"
    local SV_regions_entire="$6"
    local changeg="$7"
    local anno_config_file="$8"
    local vcfanno="$9"

    sample_outd="${outd}/${sample}"
    mkdir -p "${sample_outd}"

    of="${sample_outd}/${sample}"
    bcftools mpileup -B -a QS -Ou -f "${reffn}" \
        --threads "${num_threads}" "$bam_file" | \
        bcftools call -m -Oz -o "${of}.vcf.gz"
    bcftools index "${of}.vcf.gz"

    bcftools mpileup -B -a QS -Ou -f "${reffn}" \
        --threads "${num_threads}" "$bam_file" | \
        bcftools call -m -Ov -o "${of}.vcf"

    mkdir -p "${sample_outd}/change_to_hemi"
    mkdir -p "${outd}/annotated_vcfs/${sample}"

    sv_regions_input="${SV_regions_entire}"
    while read -r sv_region; do
        sv_name=$(echo "$sv_region" | cut -f4)
        sv_genotype=$(grep -P "^$sample\t$sv_name\t" "${outd}/SV_genotype_results.txt" | cut -f3)

        if [ "$sv_genotype" == "0/1" ] || [ "$sv_genotype" == "1/0" ]; then
            grep -w "$sv_name" "${SV_regions_entire}" > "${sample_outd}/${sv_name}.bed"
            sv_regions_input="${sample_outd}/${sv_name}.bed"
            break
        fi
    done < "$SV_regions_entire"

    output_vcf="${sample_outd}/change_to_hemi/${sample}.vcf"

    if [ "$sv_regions_input" != "${SV_regions_entire}" ]; then
        python "${changeg}" "${of}.vcf" "${outd}/SV_genotype_results.txt" \
            "${sv_regions_input}" \
            "${output_vcf}"
        bgzip -c "${output_vcf}" > "${output_vcf}.gz"
        bcftools index "${output_vcf}.gz"
        "${vcfanno}" "${anno_config_file}" "${output_vcf}.gz" \
            > "${outd}/annotated_vcfs/${sample}/${sample}_annotated.vcf"
    else
        "${vcfanno}" "${anno_config_file}" "${of}.vcf.gz" > "${outd}/annotated_vcfs/${sample}/${sample}_annotated.vcf"
    fi
}

# Main script
scratch="$PWD"
outd="${scratch}/geno_analysis/per_samp"
mkdir -p "${outd}"

bam_file="$1"
sample=$(basename "$bam_file" .bam | cut -d '_' -f 1)
reffn="$2"
num_threads="$3"
SV_regions_entire="$4"
SV_regions_1bp="$5"
changeg="/home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/get_vcf/vcf_processing.py"
anno_config_file="$6"
vcfanno="$7"

samtools addreplacerg -r ID:"${sample}" -r SM:"${sample}" \
    -o "${scratch}/$sample/${sample}.editRG.bam" "${bam_file}"
samtools index "${scratch}/$sample/${sample}.editRG.bam"

bam_path="${scratch}/$sample/${sample}.editRG.bam"

genotype_SV_regions "$bam_path" "$SV_regions_1bp" "$outd"
process_vcf "$bam_path" "$sample" "$outd" "$reffn" "$num_threads" "$SV_regions_entire" "$changeg" "$anno_config_file" "$vcfanno"