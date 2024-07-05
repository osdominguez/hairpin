#!/usr/bin/env bash

#SBATCH --job-name=format_txt
#SBATCH --time=01:00:00
#SBATCH --mem=5gb
#SBATCH --output=/home/osdominguez/output/format_HAIRPIN_txt_%a_%A.out
#SBATCH --error=/home/osdominguez/output/format_HAIRPIN_txt_%a_%A.err

awk /gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/format_txt.awk > /gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/test.txt