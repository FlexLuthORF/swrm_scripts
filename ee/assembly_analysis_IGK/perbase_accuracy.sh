$ cat perbase_accuracy.sh
#!/bin/bash
set -e -x

scratch=$PWD
samp_paths="/home/egenge01/projects/IGK/run_all_samples/36_samples.only.txt"

# Ensure ${samp_paths} is defined
if [ -z "${samp_paths}" ]; then
    echo "samp_paths is not set. Please provide a file with sample names."
    exit 1
fi

# Make output directory and define output file path
mkdir -p ${scratch}/map_to_pers_accuracy
outf=${scratch}/map_to_pers_accuracy/acc_output.txt


# Initialize the output file with headers
echo -e "Sample\tTotal_Bases\tMismatch_Bases" > $outf

# Iterate over samples from ${samp_paths}
cat ${samp_paths} | while read sample
do
    bam=${scratch}/map_to_pers_ref/${sample}/output.sorted.bam
    fasta=${scratch}/personal_refs/${sample}/ref.fasta
    
    # Validate if files exist
    if [[ ! -f "$bam" ]] || [[ ! -f "$fasta" ]]; then
        echo "BAM or FASTA file for sample ${sample} does not exist!"
        continue
    fi
    
    # Extract contigs not starting with "chr" or "igh"
    contigs=$(samtools view -H $bam | grep "@SQ" | awk '{if ($2 !~ /SN:chr/ && $2 !~ /SN:igh/) print $2}' | sed 's/SN://g')
    
    # Initialize count variables
    total_bases=0
    mismatch_bases=0
    
    # Iterate through valid contigs
    for contig in ${contigs}; do
	echo "${contig}"
        # Get total bases in contig
	contig_length=$(samtools faidx $fasta $contig | grep -v ">" | wc -c)
        total_bases=$((total_bases + contig_length))
        
        # Get bases with >25% mismatches, considering substitutions and indels
        samtools mpileup -f $fasta -r $contig $bam | head
	mismatches=$(samtools mpileup -f $fasta -r $contig $bam | \
            awk -v contig=$contig \
            'BEGIN {mismatch=0} \
            {total=length($5); \
            gsub(/[ACTGNactgn*#]/, "X", $5); \
            gsub(/[+-][0-9]+[ACGTNacgtn]+/, "X", $5); \
            mismatch=gsub("X", "", $5); \
            if ((mismatch/total) > 0.25) print contig, $2}' | wc -l)
        mismatch_bases=$((mismatch_bases + mismatches))
	echo "mismatch_bases is ${mismatch_bases}"
    done
    
    # Write results to output file
    echo -e "${sample}\t${total_bases}\t${mismatch_bases}" >> $outf
done
