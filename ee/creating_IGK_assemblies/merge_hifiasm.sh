#!/bin/bash
set -e -x

reffn=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/reference.fasta

dir=/home/egenge01/projects/IGK/run_all_samples/run_hifiasm
cat 40_samples.txt | grep HG00136 | while read sample 
do
    if [ ! -s ${dir}/${sample}/merged_bam/alg_asm20_to_ref/${sample}.sorted.bam ]
    then
	mkdir -p ${dir}/${sample}/merged_bam
	mkdir -p ${dir}/${sample}/merged_bam/alg_asm20_to_ref
	outdir=${dir}/${sample}/merged_bam/alg_asm20_to_ref
	mkdir -p ${outdir}/jobs
	# 	sbatch --time=72:00:00 -c 1 -p compute -o ${outdir}/jobs/align_job.txt --wrap="minimap2 -x asm20 -t 10 -L -a ${reffn} \
	# 	${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta > ${outdir}/${sample}.sam
	# 	 wait
	# 	 samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam 
	# 	 samtools sort -@ 10 ${outdir}/${sample}.bam -o ${outdir}/${sample}.sorted.bam
	# 	 samtools index ${outdir}/${sample}.sorted.bam
	# 	 rm -f ${outdir}/${sample}.sam
	# 	 rm -f ${outdir}/${sample}.bam"
	# fi
    	sbatch --time=72:00:00 -c 1 -p compute -o ${outdir}/jobs/align_job.txt --wrap="
         source activate IGv2
          samtools merge -f ${dir}/${sample}/merged_bam/merged.bam ${dir}/${sample}/break_at_soft_clip/2/1_asm20_hifi_asm_to_ref.sorted.bam \
    	  ${dir}/${sample}/break_at_soft_clip/2/2_asm20_hifi_asm_to_ref.sorted.bam
           wait
           samtools sort ${dir}/${sample}/merged_bam/merged.bam -o ${dir}/${sample}/merged_bam/merged.sorted.bam
            wait
            samtools index ${dir}/${sample}/merged_bam/merged.sorted.bam
            wait
              samtools fasta --reference ${reffn} ${dir}/${sample}/merged_bam/merged.bam > ${dir}/${sample}/merged_bam/merged_all_reads.fasta
             wait
             conda deactivate
              source activate seqkit-env
             seqkit rmdup --by-seq ${dir}/${sample}/merged_bam/merged_all_reads.fasta \
     	       -o ${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta"
    fi
done
