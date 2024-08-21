#!/usr/bin/env bash

#SBATCH --job-name=hairpin_analysis
#SBATCH --time=01:00:00
#SBATCH --mem=75gb
#SBATCH --output=/home/osdominguez/output/hairpin_analysis/hairpin_merge_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_analysis/hairpin_merge_%a_%A.err

module load gcc/12.1.0
module load R/4.3.1

CODE="/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/merge_CI.R"

Rscript ${CODE} $1 $2 $3 $4 $5