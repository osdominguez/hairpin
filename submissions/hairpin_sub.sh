#!/usr/bin/env bash

#SBATCH --job-name=hairpin_analysis
#SBATCH --time=12:00:00
#SBATCH --mem=10gb
#SBATCH --output=/home/osdominguez/output/hairpin_analysis/hairpin_analysis_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_analysis/hairpin_analysis_%a_%A.err

module load gcc/12.1.0
module load R/4.3.1

CODE="/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/hairpin.R"

Rscript ${CODE} $1 $2 $3 $4 $5