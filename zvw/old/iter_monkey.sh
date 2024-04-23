#!/bin/bash

set -e -x

user=$(whoami)

input_file=$1

monkey_mask_ref=/home/zmvanw01/ref/monkey/Child_17thApril2024_bothhap_IG_masked_unique.fasta

while read -r sample assemblies_fasta igh_digger igk_digger igl_digger ccs_bam; do
    outdir=$PWD/monkey_processing/${sample}
    mkdir -p $outdir
    mkdir -p $PWD/logs

    sbatch --time=88:00:00 -p compute -o $PWD/logs/${sample}_job.txt --wrap="python /home/zmvanw01/git_repos/swrm_scripts/monkey/read-support/read-support_monkey.py ${sample} ${assemblies_fasta} ${igh_digger} ${igk_digger} ${igl_digger} ${ccs_bam}"

    count=$(squeue | grep $user | wc -l)

    while [ ${count} -gt 18 ]; do
        sleep 1s
        count=$(squeue | grep $user | wc -l)
    done
done < "$input_file"