#!/bin/bash
set -e -x

user=$(whoami)

input_file=$1
monkey_mask_ref=/home/zmvanw01/ref/monkey/Child_17thApril2024_bothhap_IG_masked_unique.fasta
cat $input_file | while read sample assemblies_fasta igh_digger igk_digger igl_digger ccs_bam monkey_mask_ref
do
    outdir=$PWD/monkey_processing/${sample}
    mkdir -p $outdir
    mkdir -p $PWD/logs
    sbatch --time=88:00:00 -p compute -o $PWD/logs/${sample}_job.txt --wrap="bash /home/zmvanw01/git_repos/swrm_scripts/monkey/read-support/get_read_support_VDJs_monkey.sh ${sample} ${assemblies_fasta} ${igh_digger} ${igk_digger} ${igl_digger} ${ccs_bam}"

    count=`squeue | grep $user | wc -l`
    while [ ${count} -gt 18 ]
    do
        sleep 1s
        count=`squeue | grep $user | wc -l`
    done
done