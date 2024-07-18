#!/bin/bash

# Directory where the output file will be saved
output_dir="."

# Output file name
output_file="${output_dir}/cov_merged.tsv"

# Remove the output file if it already exists
rm -f "${output_file}"

# Loop over each file matching the pattern
find . -type f -name "*stats*" | while read -r stats_file; do
  # Extract the sample ID from the file name (before the first underscore)
  file_name=$(basename "${stats_file}")
  sample_id=${file_name%%_stats*}

  # Append the sample ID as a new column (7th column) to each row of the file and write to the output file
  awk -v id="${sample_id}" '{print $0 "\t" id}' "${stats_file}" >> "${output_file}"
done
sed -i '/^interval/d' cov_merged.tsv
