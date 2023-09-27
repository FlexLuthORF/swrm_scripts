#!/bin/bash
set -e -x
scratch=$PWD

#IMGT_alleles=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/alleles.fasta
gene_coords=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/gene_coords.bed
reffn=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/reference.fasta
extract_seq=/home/egenge01/bioinformatics/python/extract_sequence_from_bam_EEmod3.py
changeg=/home/egenge01/bioinformatics/python/change_genotypes.py
vcfanno=/home/egenge01/vcfanno-0.3.5/vcfanno
get_dens=$PWD/calculate_SNV_density.py


function run_fixed_persample {
    mkdir -p $PWD/geno_analysis/per_samp_updated
    outd=$PWD/geno_analysis/per_samp_updated
#    cp $PWD/geno_analysis/single_samp/SV_genotype_results.txt ${outd}
    # cat $PWD/36_samples_noAMR.txt | while read sample
    # do
    # 	mkdir -p ${outd}/${sample}
    # 	of=${outd}/${sample}/${sample}
    # 	bcftools mpileup -Oz -f ${reffn} \
    # 	    -B -a QS --threads 10 \
    # 	    -b $PWD/variant_analysis/bam_paths_noAMR.fofn.txt \
    # 	    -s ${sample} \
    # 	    --regions-file $PWD/regions/IGK_V_nodropout_updated.bed \
    # 	    -o ${of}.vcf.gz
    # 	bcftools index ${of}.vcf.gz
    # 	bcftools call -R $PWD/regions/IGK_V_nodropout_updated.bed \
    # 	    -m -Oz \
    # 	    -o ${of}_call.vcf.gz \
    # 	    ${of}.vcf.gz
    # 	bcftools index ${of}_call.vcf.gz
    # 	bcftools norm -a -m- -Ov -o ${of}_call_norm.vcf \
    # 	    ${of}_call.vcf.gz
    # 	bgzip -c ${of}_call_norm.vcf > ${of}_call_norm.vcf.gz
    # 	bcftools index ${of}_call_norm.vcf.gz
    # done
    cat ${outd}/SV_genotype_results.txt | grep 'IGKV1-NL1_HZ' |  while read sample pop SV_GT
    do
	of=${outd}/${sample}/${sample}
	mkdir -p ${outd}/${sample}/change_to_hemi
    	python ${changeg} ${of}_call_norm.vcf \
    	    $PWD/IGKV1-NL1_SV.bed \
    	    ${outd}/${sample}/change_to_hemi/${sample}.vcf
    done
    mkdir -p ${outd}/pre_merged_vcfs
    mkdir -p ${outd}/merged_vcfs
    cat ${outd}/SV_genotype_results.txt | grep 'IGKV1-NL1_HZ' |  while read sample pop SV_GT; do
    	cp ${outd}/${sample}/change_to_hemi/${sample}.vcf ${outd}/pre_merged_vcfs
    done    
    cat ${outd}/SV_genotype_results.txt | grep -v 'IGKV1-NL1_HZ' |  while read sample pop SV_GT; do
    	cp ${outd}/${sample}/${sample}_call_norm.vcf ${outd}/pre_merged_vcfs
    done
    find "${outd}/pre_merged_vcfs" -name "*.vcf" > "${outd}/pre_merged_vcfs/vcfs.fofn.txt"
	cat ${outd}/pre_merged_vcfs/vcfs.fofn.txt | while read file; do 
    	    bgzip -c ${file} > ${file}.gz
    	    bcftools index ${file}.gz
    	    echo ${file}.gz >> ${outd}/pre_merged_vcfs/zipped_vcfs.fofn.txt
	done
	bcftools merge -m both -l ${outd}/pre_merged_vcfs/zipped_vcfs.fofn.txt \
    	    -Ov -o ${outd}/merged_vcfs/merged.vcf
	wait
	bgzip -c ${outd}/merged_vcfs/merged.vcf > ${outd}/merged_vcfs/merged.vcf.gz
	bcftools index ${outd}/merged_vcfs/merged.vcf.gz
	
	${vcfanno} $PWD/variant_analysis/config.toml ${outd}/merged_vcfs/merged.vcf \
    	    > ${outd}/merged_vcfs/merged_anno.vcf
	bgzip -c ${outd}/merged_vcfs/merged_anno.vcf > ${outd}/merged_vcfs/merged_anno.vcf.gz
	bcftools index ${outd}/merged_vcfs/merged_anno.vcf.gz
	bcftools norm -a -m- -Ov -o ${outd}/merged_vcfs/merged_anno_norm.vcf \
	    ${outd}/merged_vcfs/merged_anno.vcf.gz
	bgzip -c ${outd}/merged_vcfs/merged_anno_norm.vcf > ${outd}/merged_vcfs/merged_anno_norm.vcf.gz
	bcftools index ${outd}/merged_vcfs/merged_anno_norm.vcf.gz
	bcftools +fill-tags ${outd}/merged_vcfs/merged_anno_norm.vcf.gz \
	    -Ov -o ${outd}/merged_vcfs/merged_anno_norm_filltags.vcf \
	    -- -S $PWD/36_noAMR_samp_groups.txt -t AF,AC,AC_Hemi,AC_Hom,AC_Het,MAF,NS
	bgzip -c ${outd}/merged_vcfs/merged_anno_norm_filltags.vcf > ${outd}/merged_vcfs/merged_anno_norm_filltags.vcf.gz
	bcftools index ${outd}/merged_vcfs/merged_anno_norm_filltags.vcf.gz
	bcftools view -i 'INFO/AC > 0.0' ${outd}/merged_vcfs/merged_anno_norm_filltags.vcf.gz \
	    -Oz -o ${outd}/merged_vcfs/AC_GT0.vcf.gz
	bcftools view -v snps ${outd}/merged_vcfs/AC_GT0.vcf.gz -Ov -o ${outd}/merged_vcfs/AC_GT0.snps.vcf
	bcftools view -v indels ${outd}/merged_vcfs/AC_GT0.vcf.gz -Ov -o ${outd}/merged_vcfs/AC_GT0.indels.vcf
	bgzip -c ${outd}/merged_vcfs/AC_GT0.snps.vcf > ${outd}/merged_vcfs/AC_GT0.snps.vcf.gz
	bcftools index ${outd}/merged_vcfs/AC_GT0.snps.vcf.gz
}

function run_do_filter_merged_final {
    outd=$PWD/geno_analysis/per_samp_updated
    bcftools filter -i'AC>1' -Ov -o ${outd}/merged_vcfs/AC_greaterthan1.snps.vcf ${outd}/merged_vcfs/AC_GT0.snps.vcf.gz
    bgzip -c ${outd}/merged_vcfs/AC_greaterthan1.snps.vcf > ${outd}/merged_vcfs/AC_greaterthan1.snps.vcf.gz
    bcftools index ${outd}/merged_vcfs/AC_greaterthan1.snps.vcf.gz
}

function run_get_per_sample_and_AF_files {
    outd=$PWD/geno_analysis/per_samp_updated
    #    inf=${outd}/merged_vcfs/AC_GT0.snps.vcf.gz
    inf=${outd}/merged_vcfs/AC_greaterthan1.snps.vcf
    # mkdir -p $PWD/geno_analysis/filtering/SNV_densities
    # outd=$PWD/geno_analysis/filtering/SNV_densities
    #    bedd=$PWD/regions
    #    for region in rss lpart1 introns gene_coords
    #    do
    #     	awk '{ total_length += $3 - $2 } END { print total_length }' ${bedd}/${region}_modified.bed \
    #     	    > ${outd}/${region}_length.txt
    #    done
    #    awk '{ total_length += $3 - $2 } END { print total_length }' ${bedd}/IGKV_nodropout_intergenicOnly.bed \
    #     	> ${outd}/intergenic_length.txt
    
    #    output_file="${outd}/merged_vcfs/per_sample.txt"
    output_file1="${outd}/merged_vcfs/AC1_filt_per_sample.txt"
    output_file2="${outd}/merged_vcfs/AC1_filt_AF_table.txt"
    if [ -e "$output_file1" ]; then
     	rm "$output_file1"
    fi
    if [ -e "$output_file2" ]; then
	rm "$output_file2"
    fi
    echo -e "Sample\tCHROM\tPOS\tGT\tGene\tIntrons\tLPart1\tRSS\tIntergenic" > "$output_file1"    
    bcftools query -l ${inf}.gz | while read S; do
    	bcftools view --samples $S -Ou ${inf}.gz |
	bcftools query -f "$S\t%CHROM\t%POS\t[%GT]\t%INFO/Gene\t%INFO/Introns\t%INFO/LPart1\t%INFO/RSS\t%INFO/Intergenic\n" &&
    	echo
    done >> "${output_file1}"
    
    bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/NS_AFR\t%INFO/NS_EAS\t%INFO/NS_EUR\t%INFO/NS_SAS\t%INFO/NS\t%INFO/MAF_AFR\t%INFO/MAF_EAS\t%INFO/MAF_EUR\t%INFO/MAF_SAS\t%INFO/MAF\t%INFO/AC_AFR\t%INFO/AC_EAS\t%INFO/AC_EUR\t%INFO/AC_SAS\t%INFO/AC_Hom_AFR\t%INFO/AC_Hom_EAS\t%INFO/AC_Hom_EUR\t%INFO/AC_Hom_SAS\t%INFO/AC_Het_AFR\t%INFO/AC_Het_EAS\t%INFO/AC_Het_EUR\t%INFO/AC_Het_SAS\t%INFO/AC_Hemi_AFR\t%INFO/AC_Hemi_EAS\t%INFO/AC_Hemi_EUR\t%INFO/AC_Hemi_SAS\t%INFO/AF_AFR\t%INFO/AF_EAS\t%INFO/AF_EUR\t%INFO/AF_SAS\t%INFO/AN\t%INFO/AC\t%INFO/Gene\t%INFO/Introns\t%INFO/LPart1\t%INFO/RSS\n' ${inf} \
	> "${output_file2}"
    


#    echo -e "Sample\tCHROM\tPOS\tGT\tAC\tGene\tIntrons\tLPart1\tRSS\tIntergenic\tAF_AFR\tAF_EAS\tAF_EUR\tAF_SAS\tAC_AFR\tAC_EAS\tAC_EUR\tAC_SAS" > "$output_file"
#	bcftools query -f "$S\t%CHROM\t%POS\t[%GT]\t%INFO/AC\t%INFO/Gene\t%INFO/Introns\t%INFO/LPart1\t%INFO/RSS\t%INFO/Intergenic\t%INFO/AF_AFR\t%INFO/AF_EAS\t%INFO/AF_EUR\t%INFO/AF_SAS\t%INFO/AC_AFR\t%INFO/AC_EAS\t%INFO/AC_EUR\t%INFO/AC_SAS\n" &&
}

function run_subset_vcf {
    mkdir -p $PWD/geno_analysis/per_samp_updated/subset_snps
    bedd=$PWD/regions
    ind=$PWD/geno_analysis/filtering
    
    bedtools subtract -a $PWD/regions/gene_coords_modified.bed \
	-b /home/egenge01/projects/IGK/data/make_franken/make_annotations/putative_pseudogenes.bed \
	> $PWD/regions/gene_coords_nopseudo_modified.bed
    bedtools intersect -a $PWD/regions/gene_coords_nopseudo_modified.bed \
	-b /home/egenge01/projects/IGK/data/make_franken/make_annotations/regions/franken_proximal_no_gap.bed \
	> $PWD/regions/gene_coords_nopseudo_prox_modified.bed
    bedtools intersect -a $PWD/regions/gene_coords_nopseudo_modified.bed \
	-b /home/egenge01/projects/IGK/data/make_franken/make_annotations/regions/franken_distal_no_gap.bed \
	> $PWD/regions/gene_coords_nopseudo_dist_modified.bed
    for region in gene_coords_nopseudo_prox gene_coords_nopseudo_dist #gene_coords_nopseudo rss lpart1 introns gene_coords
    do
	bedtools intersect -header -a ${ind}/AC_filter.snps.vcf -b ${bedd}/${region}_modified.bed \
	    > $PWD/geno_analysis/filtering/subset_snps/AC_filter.snps_${region}.vcf
    done
}

function run_prox_dist_density {
script=/home/egenge01/projects/IGK/run_all_samples/geno_analysis/filtering/calc_SNV_density_from_per_sampletxt.py
#sliding_script=/home/egenge01/projects/IGK/run_all_samples/sliding_window_SNVs.sh


}

function run_do_vcf_persamp {
#     mkdir -p $PWD/geno_analysis/single_samp
#     outd=$PWD/geno_analysis/single_samp
#     outf=${outd}/new_ms_vcf_update2.vcf
#     bcftools mpileup -Ov -f ${reffn} \
#     	-B -I --threads 10 \
# 	-a QS \
#     	-b $PWD/variant_analysis/bam_paths_noAMR.fofn.txt \
#     	--regions-file $PWD/regions/IGK_CJV_nodropout.bed \
#     	-o ${outf}
#     #	-s HG00136 \
#     bgzip -c ${outf} > ${outf}.gz
#     bcftools index ${outf}.gz
#     mkdir -p $PWD/geno_analysis/single_samp/norm_option
#     nod=$PWD/geno_analysis/single_samp/norm_option
#     bcftools call -R $PWD/regions/IGK_CJV_nodropout.bed \
#     	-mv -Oz \
#     	-o ${nod}/ms_called.vcf.gz \
#     	${outf}.gz
#     bcftools index ${nod}/ms_called.vcf.gz 
#     bcftools norm -a -m- -Ov -o ${nod}/ms_called_norm.vcf \
# 	${nod}/ms_called.vcf.gz 
#     bgzip -c ${nod}/ms_called_norm.vcf > ${nod}/ms_called_norm.vcf.gz
#     #	-S $PWD/36_samples_noAMR.txt \
    
#     	${outd}/multisample_multiallelic_all.vcf.gz
#     wait
#     bgzip -c ${outd}/multisample_consensus_variantOnly_called.vcf \
#      	> ${outd}/multisample_consensus_variantOnly_called.vcf.gz
#     tabix ${outd}/multisample_consensus_variantOnly_called.vcf.gz
#     bcftools index ${outd}/multisample_consensus_variantOnly_called.vcf.gz
#     bcftools +fill-tags ${outd}/multisample_consensus_variantOnly_called.vcf.gz \
#      	-Ov -o ${outd}/ms_consensus_AF_HWE.vcf -- -S $PWD/36_noAMR_samp_groups.txt -t AF,AC,AC_Hemi,AC_Hom,AC_Het,MAF,NS

    output_file="${outd}/SV_genotype_results.txt"
    # Check if the output file exists and delete it if it does
    if [ -e "$output_file" ]; then
        rm "$output_file"
    fi
    cat $PWD/36_noAMR_samp_groups.txt | while read sample pop
    do
    	bam_path=$PWD/variant_analysis/${sample}/${sample}_editRG2.bam
    	mpileup_output=$(samtools mpileup -l $PWD/IGKV1-NL1_region.bed "$bam_path" -t DP | head -1)
    	genotype=$(echo "$mpileup_output" | awk '{print $5}')
    	asterisk_count=$(echo "$genotype" | tr -cd '*' | wc -c)
    	# Decide genotype label based on asterisk count
    	if [ "$asterisk_count" -eq 0 ]; then
            genotype_label="IGKV1-NL1_INS.INS"
    	elif [ "$asterisk_count" -eq 1 ]; then
            genotype_label="IGKV1-NL1_HZ"
    	elif [ "$asterisk_count" -eq 2 ]; then
            genotype_label="IGKV1-NL1_DEL.DEL"
    	else
            genotype_label="unknown"
    	fi
    	echo -e "$sample\t$pop\t$genotype_label" >> "${outd}/SV_genotype_results.txt"
    done
    # mkdir -p ${outd}/split_new
    # outdir=${outd}/split_new
    # bgzip -c ${outd}/ms_consensus_AF_HWE.vcf > ${outd}/ms_consensus_AF_HWE.vcf.gz
    # bcftools index ${outd}/ms_consensus_AF_HWE.vcf.gz
    # bcftools +split -Ov -o ${outdir} \
    # 	${outd}/ms_consensus_AF_HWE.vcf.gz
    # mkdir -p ${outd}/change_to_hemi_new
    # cat ${outd}/SV_genotype_results.txt | grep 'IGKV1-NL1_HZ' |  while read sample pop SV_GT
    # do
    # 	python ${changeg} ${outd}/split/${sample}.vcf \
    # 	    $PWD/IGKV1-NL1_SV.bed \
    # 	    ${outd}/change_to_hemi/${sample}.vcf
    # done
    # mkdir -p ${outd}/pre_merged_vcfs
    # mkdir -p ${outd}/merged_vcfs
    # cat ${outd}/SV_genotype_results.txt | grep 'IGKV1-NL1_HZ' |  while read sample pop SV_GT; do
    # 	cp ${outd}/change_to_hemi/${sample}.vcf ${outd}/pre_merged_vcfs
    # done    
    # cat ${outd}/SV_genotype_results.txt | grep -v 'IGKV1-NL1_HZ' |  while read sample pop SV_GT; do
    # 	cp ${outd}/split/${sample}.vcf ${outd}/pre_merged_vcfs
    # done
    # find "${outd}/pre_merged_vcfs" -name "*.vcf" > "${outd}/pre_merged_vcfs/vcfs.fofn.txt"
    # cat ${outd}/pre_merged_vcfs/vcfs.fofn.txt | while read file; do 
    # 	bgzip -c ${file} > ${file}.gz
    # 	bcftools index ${file}.gz
    # 	echo ${file}.gz >> ${outd}/pre_merged_vcfs/zipped_vcfs.fofn.txt
    # done
    # bcftools merge -m both -l ${outd}/pre_merged_vcfs/zipped_vcfs.fofn.txt \
    # 	-Ov -o ${outd}/merged_vcfs/merged_change_hemi.vcf
    # wait
    # bgzip -c ${outd}/merged_vcfs/merged_change_hemi.vcf > ${outd}/merged_vcfs/merged_change_hemi.vcf.gz
    # bcftools index ${outd}/merged_vcfs/merged_change_hemi.vcf.gz

    # mkdir -p ${outd}/annotated_ms_vcf
    # ${vcfanno} $PWD/variant_analysis/config.toml ${outd}/merged_vcfs/merged_change_hemi.vcf \
    # 	> ${outd}/annotated_ms_vcf/merged_change_hemi_anno.vcf
    
    # mkdir -p ${outd}/filtering
    # ind2=${outd}/annotated_ms_vcf
    # outdir2=${outd}/filtering
    # #with -any I will split the multiallelic variants (SNPs+INDELs) into several records
    # bgzip ${outd}/annotated_ms_vcf/merged_change_hemi_anno.vcf -c \
    # 	> ${outd}/annotated_ms_vcf/merged_change_hemi_anno.vcf.gz
    # bcftools index ${outd}/annotated_ms_vcf/merged_change_hemi_anno.vcf.gz
    # bcftools norm -m -any ${ind2}/merged_change_hemi_anno.vcf.gz -o ${outd}/merged_change_hemi_anno.norm.vcf.gz -Oz
    # #Filter variants to include only those that occur in >1 sample (use AC > 1)
    # bcftools view -i 'INFO/AC > 1.0' ${outd}/merged_change_hemi_anno.norm.vcf.gz -Oz -o ${outdir2}/AC_filter.vcf.gz
    # #split into SNVs vs INDELs, see if you can include all homozygous REF for each
    # bcftools view -v snps ${outdir2}/AC_filter.vcf.gz -Ov -o ${outdir2}/AC_filter.snps.vcf
    # bcftools view -v indels ${outdir2}/AC_filter.vcf.gz -Ov -o ${outdir2}/AC_filter.indels.vcf
    # bgzip -c ${outdir2}/AC_filter.snps.vcf > ${outdir2}/AC_filter.snps.vcf.gz
    # bcftools index ${outdir2}/AC_filter.snps.vcf.gz
}

function run_add_AF_HWE {
#bcftools +fill-tags your_input.bcf -Ob -o output_filled.bcf -- -S sample-groups.txt -t AF,HWE
    bcftools +fill-tags $PWD/geno_analysis/multisample_multiallelic_all_called.vcf.gz \
	-Ov -o $PWD/geno_analysis/ms_ma_AF_HWE.vcf -- -S $PWD/36_noAMR_samp_groups.txt -t AF,AC,AC_Hemi,AC_Hom,AC_Het,MAF,NS
}

function run_gt_SV {
    output_file="$PWD/variant_analysis/SV_genotype_results.txt"
    # Check if the output file exists and delete it if it does
    if [ -e "$output_file" ]; then
        rm "$output_file"
    fi
    cat $PWD/36_noAMR_samp_groups.txt | while read sample pop
    do
	bam_path=$PWD/variant_analysis/${sample}/${sample}_editRG2.bam
	mpileup_output=$(samtools mpileup -l $PWD/IGKV1-NL1_region.bed "$bam_path" -t DP | head -1)
	genotype=$(echo "$mpileup_output" | awk '{print $5}')
	asterisk_count=$(echo "$genotype" | tr -cd '*' | wc -c)
	# Decide genotype label based on asterisk count
	if [ "$asterisk_count" -eq 0 ]; then
            genotype_label="IGKV1-NL1_INS.INS"
	elif [ "$asterisk_count" -eq 1 ]; then
            genotype_label="IGKV1-NL1_HZ"
	elif [ "$asterisk_count" -eq 2 ]; then
            genotype_label="IGKV1-NL1_DEL.DEL"
	else
            genotype_label="unknown"
	fi
	echo -e "$sample\t$pop\t$genotype_label" >> "$PWD/variant_analysis/SV_genotype_results.txt"
    done

}

function run_vcfanno {
    mkdir -p $PWD/geno_analysis/annotated_ms_vcf
    ${vcfanno} $PWD/variant_analysis/config.toml $PWD/geno_analysis/merged_vcfs/merged_change_hemi.vcf \
	> $PWD/geno_analysis/annotated_ms_vcf/merged_change_hemi_anno.vcf
}

function run_filter_getSNP {
    mkdir -p $PWD/geno_analysis/filtering
    ind=$PWD/geno_analysis/annotated_ms_vcf
    outd=$PWD/geno_analysis/filtering
    #with -any I will split the multiallelic variants (SNPs+INDELs) into several records
    bgzip $PWD/geno_analysis/annotated_ms_vcf/merged_change_hemi_anno.vcf -c \
	> $PWD/geno_analysis/annotated_ms_vcf/merged_change_hemi_anno.vcf.gz
    bcftools index $PWD/geno_analysis/annotated_ms_vcf/merged_change_hemi_anno.vcf.gz
    bcftools norm -m -any ${ind}/merged_change_hemi_anno.vcf.gz -o ${outd}/merged_change_hemi_anno.norm.vcf.gz -Oz
    #Filter variants to include only those that occur in >1 sample (use AC > 1)
    bcftools view -i 'INFO/AC > 1.0' ${outd}/merged_change_hemi_anno.norm.vcf.gz -Oz -o ${outd}/AC_filter.vcf.gz
    #split into SNVs vs INDELs, see if you can include all homozygous REF for each
    bcftools view -v snps ${outd}/AC_filter.vcf.gz -Ov -o ${outd}/AC_filter.snps.vcf
    bcftools view -v indels ${outd}/AC_filter.vcf.gz -Ov -o ${outd}/AC_filter.indels.vcf
    bgzip -c ${outd}/AC_filter.snps.vcf > ${outd}/AC_filter.snps.vcf.gz
    bcftools index ${outd}/AC_filter.snps.vcf.gz
}

function run_make_AF_table {
    mkdir -p $PWD/geno_analysis/filtering/AF_data
    ind=$PWD/geno_analysis/filtering
    outd=$PWD/geno_analysis/filtering/AF_data
    for ft in snps #indels
    do
	bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/NS_AFR\t%INFO/NS_EAS\t%INFO/NS_EUR\t%INFO/NS_SAS\t%INFO/NS\t%INFO/MAF_AFR\t%INFO/MAF_EAS\t%INFO/MAF_EUR\t%INFO/MAF_SAS\t%INFO/MAF\t%INFO/AC_AFR\t%INFO/AC_EAS\t%INFO/AC_EUR\t%INFO/AC_SAS\t%INFO/AC_Hom_AFR\t%INFO/AC_Hom_EAS\t%INFO/AC_Hom_EUR\t%INFO/AC_Hom_SAS\t%INFO/AC_Het_AFR\t%INFO/AC_Het_EAS\t%INFO/AC_Het_EUR\t%INFO/AC_Het_SAS\t%INFO/AC_Hemi_AFR\t%INFO/AC_Hemi_EAS\t%INFO/AC_Hemi_EUR\t%INFO/AC_Hemi_SAS\t%INFO/AF_AFR\t%INFO/AF_EAS\t%INFO/AF_EUR\t%INFO/AF_SAS\t%INFO/AN\t%INFO/AC\t%INFO/Gene\t%INFO/Introns\t%INFO/LPart1\t%INFO/RSS\n' ${ind}/AC_filter.${ft}.vcf \
	    > ${outd}/AC_filter.${ft}.table.txt
    done
}

function get_SNV_density {
    outd=$PWD/geno_analysis/filtering
    #inf=${ind}/AC_filter.snps.vcf
    # mkdir -p $PWD/geno_analysis/filtering/SNV_densities
    # outd=$PWD/geno_analysis/filtering/SNV_densities
    bedd=$PWD/regions
    for region in rss lpart1 introns gene_coords
    do
     	awk '{ total_length += $3 - $2 } END { print total_length }' ${bedd}/${region}_modified.bed \
     	    > ${outd}/${region}_length.txt
    done
    awk '{ total_length += $3 - $2 } END { print total_length }' ${bedd}/IGKV_nodropout_intergenicOnly.bed \
     	> ${outd}/intergenic_length.txt

    # output_file="$PWD/geno_analysis/filtering/per_sample.txt"

    # if [ -e "$output_file" ]; then
    # 	rm "$output_file"
    # fi
    
    # echo -e "Sample\tCHROM\tPOS\tGT\tGene\tIntrons\tLPart1\tRSS\tIntergenic" > "$output_file"
    
    # bcftools query -l ${inf}.gz | while read S; do
    # 	bcftools view --samples $S -Ou ${inf}.gz |
    # 	bcftools query -f "$S\t%CHROM\t%POS\t[%GT]\t%INFO/Gene\t%INFO/Introns\t%INFO/LPart1\t%INFO/RSS\t%INFO/Intergenic\n" &&
    # 	echo
    # done >> "$output_file"
}

function run_subset_vcf {
    mkdir -p $PWD/geno_analysis/filtering/subset_snps
    bedd=$PWD/regions
    ind=$PWD/geno_analysis/filtering
    
    bedtools subtract -a $PWD/regions/gene_coords_modified.bed \
	-b /home/egenge01/projects/IGK/data/make_franken/make_annotations/putative_pseudogenes.bed \
	> $PWD/regions/gene_coords_nopseudo_modified.bed
    bedtools intersect -a $PWD/regions/gene_coords_nopseudo_modified.bed \
	-b /home/egenge01/projects/IGK/data/make_franken/make_annotations/regions/franken_proximal_no_gap.bed \
	> $PWD/regions/gene_coords_nopseudo_prox_modified.bed
    bedtools intersect -a $PWD/regions/gene_coords_nopseudo_modified.bed \
	-b /home/egenge01/projects/IGK/data/make_franken/make_annotations/regions/franken_distal_no_gap.bed \
	> $PWD/regions/gene_coords_nopseudo_dist_modified.bed
    for region in gene_coords_nopseudo_prox gene_coords_nopseudo_dist #gene_coords_nopseudo rss lpart1 introns gene_coords
    do
	bedtools intersect -header -a ${ind}/AC_filter.snps.vcf -b ${bedd}/${region}_modified.bed \
	    > $PWD/geno_analysis/filtering/subset_snps/AC_filter.snps_${region}.vcf
    done
}

#run_gt_SV
#run_do_vcf
#run_add_AF_HWE
#run_splitvcf
#run_change_geno_hemiz
#run_merge_vcfs
#run_vcfanno
#run_filter_getSNP
#run_calc_af_population
#run_make_AF_table
#get_SNV_density
#run_do_vcf_persamp
#run_subset_vcf
#run_fixed_persample
#run_do_filter_merged_final
run_get_per_sample_and_AF_files
