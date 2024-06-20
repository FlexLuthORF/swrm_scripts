#!/bin/bash
set -e -x

function run_vcf2mt {
    scratch=$PWD
    vcf_dir=/home/egenge01/projects/CW50/multisamp_april2024
    python $PWD/vcf_to_matrix2_mod.py \
        ${vcf_dir}/multisample_final_snps_common_filt.vcf \
        $PWD/177_donors_names_only.txt \
        ${scratch}/multisamp_april2024/snps_matrix_biallelic.txt \
        biallelic
}

run_vcf2mt