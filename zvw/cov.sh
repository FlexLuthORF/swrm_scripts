#!/bin/bash
set -e -x

scratch=$PWD
sample=$1
ccs=$2
reffn=/home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta
#refbed=/home/zmvanw01/test-beds/sorted_region.bed
#refbed=/home/zmvanw01/projects/t_Parks/parks.bed
#refbed=/home/zmvanw01/240520-coverage.bed
refbed=/home/zmvanw01/IG_loci.bed

function run_get_ccs_cov {
    mkdir -p ${scratch}/ccs_cov
    outd=${scratch}/ccs_cov
    bam_path=${scratch}/ccs_cov/${sample}/ccs_to_ref.sorted.bam
    #bamtocounts --coords ${refbed} ${bam_path} > ${outd}/${sample}/${sample}_ccs-cov_counts.bed
    bamtocov --regions ${refbed} --report ${outd}/$sample/${sample}_stats.tsv ${bam_path} > ${outd}/$sample/${sample}_cov-cov.bed
}
function map_ccs_to_ref {
    dir=$scratch/ccs_cov
    outdir=${scratch}/ccs_cov/${sample}
    mkdir -p $outdir
    samtools view ${ccs} | awk '{ print ">"$1"\n"$10 }' > ${outdir}/reads.fasta
    samtools faidx ${outdir}/reads.fasta
    minimap2 -x map-hifi --secondary=no -t 11 -L -a ${reffn} ${outdir}/reads.fasta > ${outdir}/${sample}.sam
    samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam
    samtools sort -@ 10 ${outdir}/${sample}.bam -o ${outdir}/ccs_to_ref.sorted.bam
    samtools index ${outdir}/ccs_to_ref.sorted.bam
}
map_ccs_to_ref
run_get_ccs_cov
