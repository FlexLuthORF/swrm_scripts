import os
import subprocess
import argparse
from multiprocessing import Pool

def process_sample(sample_data):
    subdir, output_base = sample_data
    for file in os.listdir(subdir):
        if file.endswith('.bam') and not file.endswith('.bam.pbi'):
            sample_id = os.path.basename(subdir)
            bam_path = os.path.join(subdir, file)
            output_dir = os.path.join(output_base, sample_id)

            # Create output directory if not exists
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # Run samtools index
            samtools_cmd = "samtools index {} -@ 10".format(bam_path)
            subprocess.call(samtools_cmd, shell=True)

            # Run IG phase
            output_path = os.path.join(output_dir, "{}.output".format(sample_id))
            ig_phase_cmd = "IG phase {} {} --threads 10".format(bam_path, output_path)
            subprocess.call(ig_phase_cmd, shell=True)

def main():
    parser = argparse.ArgumentParser(description='Process .bam files in parallel.')
    parser.add_argument('root_dir', help='Root directory containing subdirectories with .bam files')
    args = parser.parse_args()

    root_dir = args.root_dir
    print("root dir is "+root_dir)
    output_base = os.path.join(root_dir, 'run_igenotyper')
    print("output base is "+output_base)
    sample_dirs = [(os.path.join(root_dir, subdir), output_base) for subdir in os.listdir(root_dir)
               if os.path.isdir(os.path.join(root_dir, subdir)) and subdir != "run_igenotyper"
               and os.path.join(root_dir, subdir) != output_base]
    print(sample_dirs)
    pool = Pool()
    pool.map(process_sample, sample_dirs)
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()

