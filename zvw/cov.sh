#!/bin/bash
set -e -x

scratch=$PWD

reffn=/home/egenge01/projects/IGL_ref_mod/reference_ready/modified_reference_renamed.fasta
refbed=/home/egenge01/projects/CW61/regions_for_cov_windows.bed

function run_get_ccs_cov {
    mkdir -p ${scratch}/cov_analysis/ccs_stats
    outd=${scratch}/cov_analysis/ccs_stats
    dd=/home/watsonlab/project/CW65/Seq_CW65_2-3-12-plex/minimap
    tail -n +2 swapped.csv | cut -d',' -f1 | while read sample
    do
        bam_path=$dd/${sample}/output.sorted.bam
        bamtocov --regions ${refbed} --report ${outd}/${sample}_stats.tsv ${bam_path} > ${outd}/${sample}_cov.bed
    done
}

run_get_ccs_cov
