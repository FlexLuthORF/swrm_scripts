function run_IG {
    scratch=$PWD
    export SJOB_DEFALLOC=NONE
    mkdir -p ${scratch}/run_igenotyper_bed-changes
    mkdir -p ${scratch}/IG_jobs_bed-changes
    
    cat /home/egenge01/projects/CW50/filtered_2023-12-2023_fofn.revioSamples.txt | while read sample data_path
    #do
     #   dir1="/home/egenge01/projects/CW50/run_igenotyper_bed-changes/${sample}"
     #   if [ -d "$dir1" ]; then
     #       echo "Directory $dir exists. Skipping ${sample}..."
     #       continue
     #   fi
     #   file1="/home/egenge01/projects/CW50/run_igenotyper2/$sample/alignments/contigs_to_ref_phased.sorted.bam"
     #   file2="/home/egenge01/projects/CW50/run_igenotyper2/$sample/alignments/ccs_to_ref_phased.sorted.bam"
     #   
     #   if [ -e "$file1" ] || [ -e "$file2" ]; then
     #       echo "Either $file1 or $file2 exists for sample $sample. Skipping..."
     #       continue
        fi
        samtools index ${data_path}
        wait
        sbatch --time=72:00:00 -p compute -o ${scratch}/IG_jobs3/${sample}_phase.txt \
            IG phase \
            --sample ${sample} \
            --threads 11 \
            ${data_path} \
            ${scratch}/run_igenotyper_bed-changes/${sample}
        
        sbatch --time=24:00:00 -p compute -o IG_jobs3/${sample}_assembly.txt \
            IG assembly \
            --threads 11 \
            ${scratch}/run_igenotyper_bed-changes/${sample}

        sbatch --time=8:00:00 -p compute -o IG_jobs3/${sample}_detect.txt \
            IG detect \
            ${scratch}/run_igenotyper_bed-changes/${sample}

    done
}