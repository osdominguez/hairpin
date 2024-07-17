#!/usr/bin/env bash

#SBATCH --job-name=format_txt
#SBATCH --time=01:00:00
#SBATCH --mem=5gb
#SBATCH --output=/home/osdominguez/output/formatting/format_HAIRPIN_txt_%a_%A.out
#SBATCH --error=/home/osdominguez/output/formatting/format_HAIRPIN_txt_%a_%A.err

cd "/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files"

file1="assess_no_.txt"
file2="pcs.txt"
file3="phens.txt"
output_file="linearity.txt"

# Check if all input files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ] || [ ! -f "$file3" ]; then
    echo "Error: One or more input files do not exist."
    exit 1
fi

# Initialize the output file
> "$output_file"

# Read each line from file1
while IFS= read -r line1 || [ -n "$line1" ]; do
    # Read each line from file2
    while IFS= read -r line2 || [ -n "$line2" ]; do
        # Read each line from file3
        while IFS= read -r line3 || [ -n "$line3" ]; do
            echo "$line1 $line2 $line3" >> "$output_file"
        done < "$file3"
    done < "$file2"
done < "$file1"

echo "Combinations generated and saved in ${output_file}."