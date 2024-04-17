# parsing VDJbase imported_genes.csv
sample=$1
contig_name=$2
allele_name=$4
gene=$3
scratch=/home/u0jana01/scratch_hifi_read_support_${sample}

mkdir -p ${scratch}


# loading alignment files
ccs_bam=~/projects/ighc/final_IGenotyper_data/${sample}/alignments/alignment_to_ighc/ccs_to_ref_phased.sorted.minimap2_realigned.grouped.bam

if [ "$gene" == "IGHG4" ] || [ "$gene" == "IGHG4D" ]; then
    contig_bam=~/VDJbase-Genomics/results/IGHC/hifiasm-assemblies/T2T/${sample}/${sample}/merged.bam

    cat ~/VDJbase-Genomics/results/IGHC/gene_lists/hifiasm_imported_genes.csv | awk -F',' -v sample=${sample} '$3 == sample' | awk -F',' -v gene_name=${gene} '{if ($6==gene_name) print $0}' | awk -F',' -v c_n=${contig_name} '{if ($5==c_n) print $0}'| tr ',' '\t' | cut -f8,17 | sort | uniq | awk '{print ">"$1"\n"$2}' > ${scratch}/allele_seq.fasta


    cat ~/VDJbase-Genomics_new/VDJbase-Genomics/results/IGHC/gene_lists/hifiasm_all_genes.csv | awk -F',' -v sample=${sample} '$3 == sample' | awk -F',' -v gene_name=${gene} '{if ($6==gene_name) print $0}' | awk -F',' -v c_n=${contig_name} '{if ($4==c_n) print $0}'| tr ',' '\t'| awk -F'\t' '{print $4"\t"$27"\t"$28"\n"$4"\t"$25"\t"$26"\n"$4"\t"$23"\t"$24"\n"$4"\t"$21"\t"$22"\n"$4"\t"$75"\t"$76"\n"$4"\t"$73"\t"$74"\n"$4"\t"$87"\t"$88"\n"$4"\t"$85"\t"$86"\n"$4"\t"$83"\t"$84}' | awk '{if ($3!=$2) print $0"\t"$3-$2+1}' > ${scratch}/personalized_aln_bound.bed 

    cat ~/VDJbase-Genomics_new/VDJbase-Genomics/results/IGHC/gene_lists/hifiasm_imported_genes.csv | grep ${sample} |awk -F',' '{print $5}'| sort | uniq > ${scratch}/contig_extract

else
    contig_bam=~/projects/ighc/final_IGenotyper_data/${sample}/alignments/alignment_to_ighc/contigs_to_ref_phased.sorted.minimap2_realigned.grouped_mapq_edited_sftcc.bam

    cat ~/VDJbase-Genomics/results/IGHC/gene_lists/IGenotyper_imported_genes.csv | awk -F',' -v sample=${sample} '{if ($3 ==  sample) print $0}' | awk -F',' -v gene=${gene} '{ if ($6 == gene) print $0}' | awk -F',' -v contig_name=${contig_name} '{ if ($5 == contig_name) print $0}' | tr ',' '\t' | cut -f8,17 | sort | uniq | awk '{ print ">"$1"\n"$2 }' > ${scratch}/allele_seq.fasta


    cat ~/VDJbase-Genomics_new/VDJbase-Genomics/results/IGHC/gene_lists/IGenotyper_all_genes.csv | awk -F',' -v sample=${sample} '{if ($3 ==  sample) print $0}' | awk -F',' -v gene=${gene} '{ if ($6 == gene) print $0}' | awk -F',' -v contig_name=${contig_name} '{ if ($4 == contig_name) print $0}' | tr ',' '\t' | awk -F'\t' '{print $4"\t"$27"\t"$28"\n"$4"\t"$25"\t"$26"\n"$4"\t"$23"\t"$24"\n"$4"\t"$21"\t"$22"\n"$4"\t"$75"\t"$76"\n"$4"\t"$73"\t"$74"\n"$4"\t"$87"\t"$88"\n"$4"\t"$85"\t"$86"\n"$4"\t"$83"\t"$84}' | awk '{if ($3!=$2) print $0"\t"$3-$2+1}' > ${scratch}/personalized_aln_bound.bed

    cat ~/VDJbase-Genomics_new/VDJbase-Genomics/results/IGHC/gene_lists/IGenotyper_imported_genes.csv | grep ${sample} | awk -F',' '{print $5}' | sort | uniq > ${scratch}/contig_extract

fi



# extracting all the hifi reads mapped to ighc
samtools view ${ccs_bam} ighc | awk '{print ">"$1"\n"$10}' > ${scratch}/ighc_hifi_reads.fasta

# extracting all the contigs representing the ighc locus:
samtools view -F 3884 ${contig_bam} ighc | awk '{print ">"$1"\n"$10}' > ${scratch}/ighc_contigs_temp.fasta

seqkit rmdup -n ${scratch}/ighc_contigs_temp.fasta -o ${scratch}/ighc_contigs_temp_1.fasta
cat ${scratch}/ighc_contigs_temp_1.fasta | awk '/^>/ {printf("%s%s\n",(NR==1)?"":"\n",$0);next} {printf("%s",$0)} END {printf("\n")}' > ${scratch}/ighc_contigs.fasta
samtools faidx ${scratch}/ighc_contigs.fasta 

rm ${scratch}/ighc_contigs_temp.fasta ${scratch}/ighc_contigs_temp_1.fasta

# alignment to personal ref
minimap2 -ax map-hifi ${scratch}/ighc_contigs.fasta ${scratch}/ighc_hifi_reads.fasta > ${scratch}/aln.sam
samtools view -b ${scratch}/aln.sam > ${scratch}/aln.bam
samtools sort ${scratch}/aln.bam > ${scratch}/aln_sorted.bam
samtools index ${scratch}/aln_sorted.bam

cat ${scratch}/ighc_contigs.fasta | awk -v header=${contig_name} '$1 == ">"header {print; getline; print}' > ${scratch}/contig_seq.fasta
samtools faidx ${scratch}/contig_seq.fasta

# whole gene alignment
#minimap2 -ax map-hifi ${scratch}/contig_seq.fasta ${scratch}/allele_seq_exons.fasta > ${scratch}/spliced_aln.sam
minimap2 -ax map-hifi ${scratch}/contig_seq.fasta ${scratch}/allele_seq.fasta --MD > ${scratch}/non-spliced_aln.sam
samtools view -b ${scratch}/non-spliced_aln.sam > ${scratch}/non-spliced_aln.bam
samtools sort ${scratch}/non-spliced_aln.bam > ${scratch}/non-spliced_aln_sorted.bam
samtools index ${scratch}/non-spliced_aln_sorted.bam


# creating tring_depth_ref_alt
echo -e "read_depth=\tmatch_depth=\talt_depth=" > ${scratch}/string_depth_match_alt


## tackling spliced alignment and exon specific position-based-average accuracy | no introns
cat ${scratch}/personalized_aln_bound.bed | while read chr start end len
do
    bam_file=${scratch}/aln_sorted.bam
    match_bases=0
    total_bases=$(($end-$start+1))
    matches=$(samtools mpileup -f ${scratch}/ighc_contigs.fasta -r ${chr}:${start}-${end} ${bam_file} | \
            awk -v min_depth=1\
            'BEGIN {mismatch=0} \
            {total=length($5); \
            if (total >= min_depth) { \
                gsub(/[ACTGNactgn*#]/, "X", $5); \
                gsub(/[+-][0-9]+[ACGTNacgtn]+/, "X", $5); \
                mismatch=gsub("X", "", $5); \
                if ((mismatch/total) <= 0.20) print $2}}' | awk -F'\t' -v start=${start} -v end=${end} '{if ($1 >= start && $1 <= end) print $0}'| wc -l)
        match_bases=$((match_bases + matches))

    [ -e "${scratch}/summary_over_exons" ] && echo -e "$match_bases\t$total_bases" >> ${scratch}/summary_over_exons || echo -e "$match_bases\t$total_bases" > ${scratch}/summary_over_exons
   
    
    [ -e "${scratch}/string_depth_match_alt" ] && samtools mpileup -f ${scratch}/ighc_contigs.fasta -r ${chr}:${start}-${end} ${bam_file} | \
            awk -v min_depth=0\
            'BEGIN {mismatch=0} \
            {total=length($5); \
            if (total >= min_depth) { \
                gsub(/[ACTGNactgn*#]/, "X", $5); \
                gsub(/[+-][0-9]+[ACGTNacgtn]+/, "X", $5); \
                mismatch=gsub("X", "", $5); \
                print total "\t" (total-mismatch) "\t" mismatch}}' >> ${scratch}/string_depth_match_alt || samtools mpileup -f ${scratch}/ighc_contigs.fasta -r ${chr}:${start}-${end} ${bam_file} | \
            awk -v min_depth=0\
            'BEGIN {mismatch=0} \
            {total=length($5); \
            if (total >= min_depth) { \
                gsub(/[ACTGNactgn*#]/, "X", $5); \
                gsub(/[+-][0-9]+[ACGTNacgtn]+/, "X", $5); \
                mismatch=gsub("X", "", $5); \
                print total "\t" (total-mismatch) "\t" mismatch}}' > ${scratch}/string_depth_match_alt

 
done

base_accuracy_over_exons=$(awk -F'\t' '{total_match+=$1; total_entry+=$2} END {printf "%.3f\n", ((total_match) * 100) / total_entry}' ${scratch}/summary_over_exons)

depth=$(cat ${scratch}/string_depth_match_alt | cut -f1 | tr '\n' ',')
match=$(cat ${scratch}/string_depth_match_alt | cut -f2 | tr '\n' ',')
alt=$(cat ${scratch}/string_depth_match_alt | cut -f3 | tr '\n' ',')

# samtools non-spliced non-spliced gene specific position-based-average accuracy | introns are there
samtools view -F 3884 ${scratch}/non-spliced_aln_sorted.bam | awk '{print $3"\t"$4"\t"$4 + length($10) - 1}' | while read chr start end
do
    bam_file=${scratch}/aln_sorted.bam
    mismatch_bases=0
    total_bases=$(($end-$start))

    mismatches=$(samtools mpileup -f ${scratch}/ighc_contigs.fasta -r ${chr}:${start}-${end} ${bam_file} | \
            awk \
            'BEGIN {mismatch=0} \
            {total=length($5); \
            gsub(/[ACTGNactgn*#]/, "X", $5); \
            gsub(/[+-][0-9]+[ACGTNacgtn]+/, "X", $5); \
            mismatch=gsub("X", "", $5); \
            if ((mismatch/total) > 0.20) print $2}' | wc -l)
        mismatch_bases=$((mismatch_bases + mismatches))

    base_accuracy=$(echo "scale=3; ($total_bases-$mismatch_bases)*100/$total_bases" | bc)


    count=$(samtools view -F 3884 ${scratch}/aln_sorted.bam ${chr}:${start}-${end} | cut -f1,6 | while read name cigar; do total_matches=$(echo ${cigar} | grep -oP '\d+(?=M)' | paste -sd+ - | bc); total_length=$(echo ${cigar} | grep -oP '\d+' | paste -sd+ - | bc); percentage=$(echo "scale=6; $total_matches / $total_length * 100" | bc); echo "$percentage"; done | awk '{if ($1>99.50) print $1}' | wc -l)
    hifi_support=${count}

    total_hifi=$(samtools view -F 3884 ${scratch}/aln_sorted.bam ${chr}:${start}-${end} | cut -f1,6 | wc -l)

    echo -e "$base_accuracy\t$base_accuracy_over_exons\t$depth\t$match\t$alt"
done

rm -r ${scratch}
