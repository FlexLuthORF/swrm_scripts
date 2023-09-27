#!/bin/bash
set -e -x

cat all_IGK_cov_filtered_samples.txt | grep -v relationship | cut -f1,9 | while read sample ccs
do
    outdir=$PWD/run_hifiasm/${sample}
    if [ ! -d "${outdir}" ]; then
        mkdir -p $PWD/run_hifiasm/${sample}
	mkdir -p $PWD/run_hifiasm/${sample}/job
	sbatch --time=88:00:00 -c 1 -p compute -o ${outdir}/job/job.txt --wrap="sh pipeline.sh ${outdir} ${ccs} 10"
    else
        echo "Directory ${outdir} already exists. Skipping..."
    fi
    count=`squeue | grep egenge01 | wc -l`
    while [ ${count} -gt 50 ]
    do
        sleep 1s
        count=`squeue | grep egenge01 | wc -l`
    done
done
