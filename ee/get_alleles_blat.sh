#!/bin/bash
set -e -x

scratch=$PWD

IMGT_alleles=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/alleles.fasta
gene_coords=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/gene_coords.bed
#gene_coords=/home/egenge01/projects/IGK/data/make_franken/make_annotations/IGK_merged_genes_l2_test_WL_script.bed
reffn=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/reference.fasta
#extract_seq=/home/egenge01/bioinformatics/python/extract_sequence_from_bam_EEmod2.py
extract_seq=/home/egenge01/bioinformatics/python/extract_sequence_from_bam_EEmod3.py

function run_get_alleles {
    # conda activate vcflib_env
    cat 36_samp_groups.txt | grep -v NA19467 | grep -v NA18522 | grep -v HG01925 | grep -v HG02018 | while read sample pop
    do
	#	mkdir -p $PWD/alleles_from_contigs/${sample}
	#	mkdir -p $PWD/curated_contig_fasta/${sample}
	#bam_path=$PWD/curated_24_samples/merged_alignments/${sample}/${sample}_merged.sorted.bam
	#samtools fasta -t ${bam_path} | awk '/^>/{sub(">", "> "++i " ")}1' > $PWD/curated_contig_fasta/${sample}/${sample}_h1h2_contigs.fasta
	#	python ${extract_seq} ${gene_coords} ${bam_path} > $PWD/alleles_from_contigs/${sample}/${sample}_allele_seqs.fasta
	#	samtools faidx $PWD/alleles_from_contigs/${sample}/${sample}_allele_seqs.fasta
	mkdir -p $PWD/alleles_from_contigs_blat/${sample}
	blat_out=$PWD/alleles_from_contigs_blat/${sample}/${sample}_IMGT_blat_res.txt
	/home/egenge01/blat/blat ${IMGT_alleles} $PWD/alleles_from_contigs/${sample}/${sample}_allele_seqs.fasta \
	    -out=blast8 ${blat_out}
	wait
#	Create a new directory for Rscript output
	mkdir -p $PWD/alleles_from_contigs/parse_blat/${sample}
        
        # Run your R script for this sample with input and output file names
	Rscript parse_blast_output10.R ${blat_out} $sample $PWD/alleles_from_contigs/parse_blat/${sample}
    done
}

function run_concatenate {
#    if [ -f "$PWD/alleles_from_contigs/parse_blat/concatenate/all_samples.txt" ]; then
#	rm "$PWD/alleles_from_contigs/parse_blat/concatenate/all_samples.txt"
#    fi
    cat 36_samp_groups.txt | grep -v NA19467 | grep -v NA18522 | grep -v HG01925 | grep -v HG02018 | while read sample pop
    do
	mkdir -p $PWD/alleles_from_contigs/parse_blat/concatenate
	cat $PWD/alleles_from_contigs/parse_blat/${sample}/${sample}_IGK.vs.IMGT_consolidated.txt \
	    | grep -v "sample_haplotype" >> $PWD/alleles_from_contigs/parse_blat/concatenate/all_samples_updated_orig.txt
    done
}


function run_concatenate_previous {
    if [ -f "$PWD/alleles_from_contigs/parse_blat/concatenate/all_samples.txt" ]; then
	rm "$PWD/alleles_from_contigs/parse_blat/concatenate/all_samples.txt"
    fi
    cat 36_samp_groups.txt | grep -v NA19467 | grep -v NA18522 | grep -v HG01925 | grep -v HG02018 | grep | HG02433 | while read sample pop
    do
	mkdir -p $PWD/alleles_from_contigs/parse_blat/concatenate
	cat $PWD/alleles_from_contigs/parse_blat/${sample}/${sample}_IGK.vs.IMGT_consolidated.txt | awk -F'\t' 'NR>1 {print $7, $8, $9, $6, $4, $5, $3, $10, $11}' \
	    | grep -v "sample_haplotype" >> $PWD/alleles_from_contigs/parse_blat/concatenate/all_samples.txt
    done
}



#run_get_alleles
run_concatenate
#run_concatenate_previous
