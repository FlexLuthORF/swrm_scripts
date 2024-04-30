#!/bin/bash
set -e -x

scratch=$PWD

reffn=/immune_receptor_genomics/current/reference.fasta
refbed=/home/zmvanw01/test-beds/sorted_region.bed

function run_get_ccs_cov {
    mkdir -p ${scratch}/cov_analysis
    outd=${scratch}/cov_analysis
    tail -n +2 swapped.csv | cut -d',' -f1 | while read sample
    do
        bam_path=${scratch}/${sample}/output.sorted.bam
        bamtocov --regions ${refbed} --report ${outd}/${sample}_stats.tsv ${bam_path} > ${outd}/${sample}_cov.bed
    done
}

run_get_ccs_cov
