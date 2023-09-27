#!/bin/bash
set +e -x

cur_dir=/home/egenge01/projects/IGK/run_all_samples/curated_contigs
reffn=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/reference.fasta

function run_stitch {
    cat 40_samples.txt | grep HG03018 | while read sample
    do
	#    if [ ! -s ${dir}/${sample}/merged_bam/alg_asm20_to_ref/${sample}.sorted.bam ]
	#   then
	#	if [ ! -s $PWD/curated_24_samples/merged_alignments/${sample}/${sample}_merged.sorted.bam ]
	#	then
	for k in prox_contigs #dist_contigs
	do
	    for i in 2 #1 2
	    do
		for min_overlap in 300 100 # 200 400 500
		do
		    for max_err in 25 10 #75 100 1
		    do
			fasta=${cur_dir}/${k}/${sample}/hap${i}.fasta
			mkdir -p ${cur_dir}/${k}/${sample}/blastn_res/hap${i}/min_ovlp_${min_overlap}/max_err_${max_err}
			outd=${cur_dir}/${k}/${sample}/blastn_res/hap${i}/min_ovlp_${min_overlap}/max_err_${max_err}
			blastn -query ${fasta} \
			    -subject ${fasta} \
			    -outfmt "6 length pident nident mismatch gapopen gaps qseqid qstart qend qlen sseqid sstart send slen sstrand" \
			    > ${outd}/blast.txt
			wait
			touch ${outd}/fosmids_to_ignore.txt
			python /home/egenge01/bioinformatics/Fosmids/python/merge_fosmids.py \
			    ${outd}/blast.txt \
			    ${outd}/blast_edited.txt \
			    ${outd}/fosmid_groups.txt \
			    ${outd}/fosmids_to_ignore.txt \
			    ${outd}/fosmids_to_merge.txt \
			    ${fasta} \
			    ${outd}/merge_fasta.fasta \
			    ${min_overlap} \
			    ${max_err}
			wait
			mkdir -p ${outd}/alignment
	                sbatch --time=72:00:00 -c 1 -p gpu -o ${outd}/alignment/align_job.txt --wrap="minimap2 -x asm20 \
    -t 10 -L -a ${reffn} \
    ${outd}/merge_fasta.fasta > ${outd}/alignment/${sample}_${min_overlap}_${max_err}.sam
     wait
     samtools view -Sbh ${outd}/alignment/${sample}_${min_overlap}_${max_err}.sam > ${outd}/alignment/${sample}_${min_overlap}_${max_err}.bam 
     samtools sort -@ 10 ${outd}/alignment/${sample}_${min_overlap}_${max_err}.bam -o ${outd}/alignment/${sample}_${min_overlap}_${max_err}.sorted.bam
     samtools index ${outd}/alignment/${sample}_${min_overlap}_${max_err}.sorted.bam
     rm -f ${outd}/alignment/${sample}_${min_overlap}_${max_err}.sam
     rm -f ${outd}/alignment/${sample}_${min_overlap}_${max_err}.bam"
		    done
		done
	    done
	done
    done
}

function run_cp_checked_contigs {
    contig_dir=/home/egenge01/projects/IGK/run_all_samples/curated_contigs
#    cat 40_samples.txt | grep HG03018 | while read sample
    cat 40_samples.txt | grep HG03018 | while read sample
    do
#	if [ ! -s $PWD/curated_24_samples/merged_alignments/${sample}/${sample}_merged.sorted.bam ]
#	then 
	mkdir -p $PWD/curated_24_samples/proximal/${sample}
	mkdir -p $PWD/curated_24_samples/distal/${sample}
	for k in prox_contigs dist_contigs
	do
	    for i in 1 2
	    do
		if [[ "$sample" == "NA20351" && "$k" == "dist_contigs" && "$i" == 1 ]]; then
		    min_overlap=200
		    max_err=25
		elif [[ "$sample" == "HG02700" && "$k" == "dist_contigs" && "$i" == 2 ]]; then
		    min_overlap=100
		    max_err=25
		elif [[ "$sample" == "HG02513" && "$k" == "dist_contigs" && "$i" == 2 ]]; then
		    min_overlap=100
		    max_err=50
		elif [[ "$sample" == "HG02433" && "$k" == "prox_contigs" && "$i" == 2 ]]; then
		    min_overlap=100
		    max_err=25
		elif [[ "$sample" == "NA19735" && "$k" == "prox_contigs" && "$i" == 2 ]]; then
		    min_overlap=100
		    max_err=25
		elif [[ "$sample" == "HG03786" && "$k" == "prox_contigs" && "$i" == 2 ]]; then
		    min_overlap=100
		    max_err=25
		elif [[ "$sample" == "HG03196" ]]; then
		    min_overlap=100
		    max_err=25
		elif [[ "$sample" == "HG02391" ]]; then
		    min_overlap=200
		    max_err=25
		else
		    min_overlap=300
		    max_err=25
		fi
		if [[ "$k" == "prox_contigs" ]]; then
		    jj=proximal
		elif [[ "$k" == "dist_contigs" ]]; then
		    jj=distal
		fi
		cp ${contig_dir}/${k}/${sample}/blastn_res/hap${i}/min_ovlp_${min_overlap}/max_err_${max_err}/alignment/${sample}_${min_overlap}_${max_err}.sorted.bam $PWD/curated_24_samples/${jj}/${sample}/${sample}_hap${i}.bam
		cp ${contig_dir}/${k}/${sample}/blastn_res/hap${i}/min_ovlp_${min_overlap}/max_err_${max_err}/alignment/${sample}_${min_overlap}_${max_err}.sorted.bam.bai $PWD/curated_24_samples/${jj}/${sample}/${sample}_hap${i}.bam.bai
		#		cp ${contig_dir}/dist_contigs/${sample}/blastn_res/hap${i}/min_ovlp_${min_overlap}/max_err_${max_err}/alignment/${sample}_${min_overlap}_${max_err}.sorted.bam $PWD/curated_24_samples/distal/${sample}/${sample}_hap${i}.bam
		#		cp ${contig_dir}/dist_contigs/${sample}/blastn_res/hap${i}/min_ovlp_${min_overlap}/max_err_${max_err}/alignment/${sample}_${min_overlap}_${max_err}.sorted.bam.bai $PWD/curated_24_samples/distal/${sample}/${sample}_hap${i}.bam.bai
		mkdir -p $PWD/curated_24_samples/merged_alignments/${sample}
		#cp $PWD/curated_24_samples/proximal/${sample}/${sample}_hap1.bam $PWD/curated_24_samples/proximal/${sample}/${sample}_prox_hap1.bam
		#cp $PWD/curated_24_samples/proximal/${sample}/${sample}_hap2.bam $PWD/curated_24_samples/proximal/${sample}/${sample}_prox_hap2.bam
		#cp $PWD/curated_24_samples/distal/${sample}/${sample}_hap1.bam $PWD/curated_24_samples/distal/${sample}/${sample}_dist_hap1.bam
		#cp $PWD/curated_24_samples/distal/${sample}/${sample}_hap2.bam $PWD/curated_24_samples/distal/${sample}/${sample}_dist_hap2.bam
		cp $PWD/curated_24_samples/proximal/${sample}/${sample}_hap1.bam $PWD/curated_24_samples/proximal/${sample}/hap1.bam
		cp $PWD/curated_24_samples/proximal/${sample}/${sample}_hap2.bam $PWD/curated_24_samples/proximal/${sample}/hap2.bam
		cp $PWD/curated_24_samples/distal/${sample}/${sample}_hap1.bam $PWD/curated_24_samples/distal/${sample}/hap1.bam
		cp $PWD/curated_24_samples/distal/${sample}/${sample}_hap2.bam $PWD/curated_24_samples/distal/${sample}/hap2.bam
		wait
		samtools merge -f -r $PWD/curated_24_samples/merged_alignments/${sample}/${sample}_merged.bam \
		    $PWD/curated_24_samples/proximal/${sample}/hap1.bam $PWD/curated_24_samples/distal/${sample}/hap1.bam \
		    $PWD/curated_24_samples/proximal/${sample}/hap2.bam $PWD/curated_24_samples/distal/${sample}/hap2.bam
		wait
		samtools sort $PWD/curated_24_samples/merged_alignments/${sample}/${sample}_merged.bam \
		    -o $PWD/curated_24_samples/merged_alignments/${sample}/${sample}_merged.sorted.bam
		samtools index $PWD/curated_24_samples/merged_alignments/${sample}/${sample}_merged.sorted.bam
	    done	
	done	
    done
}

#run_stitch
run_cp_checked_contigs
