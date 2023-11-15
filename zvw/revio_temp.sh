#!/bin/bash

data_path=$1
scratch=$2

/home/zmvanw01/anaconda3/envs/IGv2/bin/samtools view ${data_path} | awk '{ print ">"$1"\n"$10 }' > ${scratch}/ccs_reads.fasta
fasta=${scratch}/ccs_reads.fasta
/home/egenge01/minimap2/minimap2 -ax map-hifi -t 12 -L /home/egenge01/projects/IGK/data/make_franken/franken_ref/reference.fasta ${scratch}/ccs_reads.fasta > ${scratch}/output.sam
/home/zmvanw01/anaconda3/envs/IGv2/bin/samtools view -Sbh ${scratch}/output.sam > ${scratch}/output.bam
/home/zmvanw01/anaconda3/envs/IGv2/bin/samtools sort -@ 12 ${scratch}/output.bam -o ${scratch}/output.sorted.bam
/home/zmvanw01/anaconda3/envs/IGv2/bin/samtools index ${scratch}/output.sorted.bam
rm -f ${scratch}/output.sam
rm -f ${scratch}/output.bam
