#!/usr/bin/env bash

#SBATCH --job-name=format_hairpin_txt
#SBATCH --time=01:00:00
#SBATCH --mem=5gb
#SBATCH --output=/home/osdominguez/output/format_HAIRPIN_txt_%a_%A.out
#SBATCH --error=/home/osdominguez/output/format_HAIRPIN_txt_%a_%A.err

cd "/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files"

file1="assess_no_.txt"
file2="hairpin_phens.txt"
output_file="hairpin.txt"

boot_n=100

# Check if all input files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
    echo "Error: One or more input files do not exist."
    exit 1
fi

# Initialize the output file
> "$output_file"

# Read each line from file1
while IFS= read -r line1 || [ -n "$line1" ]; do
    # Read each line from file2
    while IFS= read -r line2 || [ -n "$line2" ]; do
            echo "${line1} ${line2} F 0" >> "$output_file"
            for (( i=1; i <= ${boot_n}; ++i ))
            do
                echo "${line1} ${line2} T ${i}" >> "$output_file"
            done
    done < "$file2"
done < "$file1"