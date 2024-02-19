function run_IG {
    scratch=$PWD
    export SJOB_DEFALLOC=NONE

    # Check and create directories
    mkdir -p ${scratch}/run_igenotyper_bed-changes || { echo "Failed to create run_igenotyper_bed-changes directory"; exit 1; }
    mkdir -p ${scratch}/IG_jobs_bed-changes || { echo "Failed to create IG_jobs_bed-changes directory"; exit 1; }

    # Check if the input file exists and readable
    input_file="/home/zmvanw01/test-beds/newbams/newer_fofn.txt"
    if [[ ! -f "$input_file" ]] || [[ ! -r "$input_file" ]]; then
        echo "Input file does not exist or is not readable: $input_file"
        exit 1
    fi

    # Processing each sample
    cat $input_file | while read sample data_path
    do
        if [[ -z "$sample" ]] || [[ -z "$data_path" ]]; then
            echo "Empty sample or data path found. Skipping..."
            continue
        fi

        # Create a script for each sample
        script_path="${scratch}/IG_jobs_bed-changes/${sample}_job.sh"
        echo "#!/bin/bash" > $script_path || { echo "Failed to create script for sample: $sample"; continue; }

        # Adding commands to the script
        {
            #echo "samtools index ${data_path}"
            echo "IG phase --sample ${sample} --threads 11 ${data_path} ${scratch}/run_igenotyper_bed-changes/${sample}"
            echo "echo phase completed"
            echo "IG assembly --threads 11 ${scratch}/run_igenotyper_bed-changes/${sample}"
            echo "echo assembly completed"
            echo "IG detect ${scratch}/run_igenotyper_bed-changes/${sample}"
            echo "echo detect completed"
        } >> $script_path || { echo "Failed to write commands for sample: $sample"; continue; }

        # Submit the script as a job
        sbatch -p compute -o ${scratch}/IG_jobs_bed-changes/${sample}_job.txt $script_path  || echo "Failed to submit job for sample: $sample"
    done
}

# Execute the function
run_IG
