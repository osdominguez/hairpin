#!/usr/bin/env bash

#SBATCH --job-name=format_tables_txt
#SBATCH --time=01:00:00
#SBATCH --mem=5gb
#SBATCH --output=/home/osdominguez/output/formatting/format_tables_txt_%a_%A.out
#SBATCH --error=/home/osdominguez/output/formatting/format_tables_txt_%a_%A.err

cd "/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files"

file1="phens.txt"
file2="pops.txt"
output_file="tables.txt"

boot_n=100

# Check if all input files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
    echo "Error: One or more input files do not exist."
    exit 1
fi

# Initialize the output file
> "$output_file"

while IFS= read -r line1 || [ -n "$line1" ]; do
    while IFS= read -r line2 || [ -n "$line2" ]; do
            echo "${line1} ${line2}" >> "$output_file"
    done < "$file2"
done < "$file1"