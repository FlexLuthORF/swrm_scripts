#!/bin/bash


sample_id=$1
hap=$2
contig=$3
gene=$4
allele=$5
allele_id=$6
pop=$7

allele_grep=$(echo "$allele" | cut -f2 -d'*')

read_support_span=$(sh personal_contig_based_hifi_support.sh ${sample_id} ${contig} ${gene} ${allele})
hifi_support=$(cat <path_to_vdjbase>/VDJbase-Genomics/results/IGHC/gene_lists/ccs-bam_imported_genes.csv | tr ',' '\t' | grep ${sample_id} | awk -v gene_name="$gene" '{if ($6==gene_name) print $0}' | cut -f8 | awk -v allele_hit="$allele" '{if ($1==allele_hit) print $0}' | wc -l )


cat <path_to_vdjbase>/VDJbase-Genomics/results/IGHC/gene_lists/ccs-bam_imported_genes.csv | tr ',' '\t' | grep ${sample_id} | awk -v gene_name="$gene" '{if ($6==gene_name) print $0}' | awk -v allele_hit="$allele" '{if ($5!=allele_hit) print $0}' | awk -F'\t' '{print ">"$5"\n"$59}' > ${sample_id}_non_supporting_alleles.fasta

alternate_allele_counts=$(python ~/allele_accuracy_check/test_by_clustering/fasta_unique_sequence_count.py ${sample_id}_non_supporting_alleles.fasta | sed ':a;N;$!ba;s/\n/,/g')

rm ${sample_id}_non_supporting_alleles.fasta

cat ~/projects/ighc/gene_coords/IGHC_T2T_gene_coords.bed | awk -v gene_name="$gene" '{if ($4==gene_name) print $1"\t"$2"\t"$3}' | while read chr start end
do
    total_full_cov_reads=$(python ~/projects/ighc/scripts/block_extract.py ~/projects/ighc/final_IGenotyper_data/${sample_id}/alignments/alignment_to_ighc/ccs_to_ref_phased.sorted.minimap2_realigned.grouped.bam ${chr} ${start} ${end} random | grep ">" | wc -l)


    no_support=$((total_full_cov_reads-hifi_support))

    echo -e "$sample_id\t$hap\t$contig\t$gene\t$allele\t$allele_id\t$pop\t$hifi_support\t$no_support\t$alternate_allele_counts\t$total_full_cov_reads\t$read_support_span"

done
