#!/usr/bin/env bash

#SBATCH --job-name=format_merge_txt
#SBATCH --time=01:00:00
#SBATCH --mem=5gb
#SBATCH --output=/home/osdominguez/output/formatting/format_merge_txt_%a_%A.out
#SBATCH --error=/home/osdominguez/output/formatting/format_merge_txt_%a_%A.err

cd "/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files"

file1="phens.txt"
file2="assess_no_.txt"
file3="pops.txt"

output_file="merge.txt"

boot_n=100

# Check if all input files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ] || [ ! -f "$file3" ]; then
    echo "Error: One or more input files do not exist."
    exit 1
fi

# Initialize the output file
> "$output_file"

while IFS= read -r line1 || [ -n "$line1" ]; do
    while IFS= read -r line2 || [ -n "$line2" ]; do
        while IFS= read -r line3 || [ -n "$line3" ]; do
            echo "${line1} ${line2} ${line3} ${boot_n}" >> "$output_file"
        done < "$file3"
    done < "$file2"
done < "$file1"