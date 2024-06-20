#!/bin/bash
set -e -x

user=$(whoami)
input_file=$1
cat $input_file | while read sample ccs
do
    outdir=$PWD/ccs_cov/${sample}
    mkdir -p $PWD/ccs_cov/${sample}
    #REMOVE PATH
    sbatch --time=88:00:00 -p compute -o ${outdir}/job.txt --wrap="sh /home/zmvanw01/git_repos/swrm_scripts/zvw/cov.sh ${sample} ${ccs}"
    count=`squeue | grep $user | wc -l`
    while [ ${count} -gt 15 ]
    do
	sleep 1s
	count=`squeue | grep $user | wc -l`
    done
done
