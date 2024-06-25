#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_PGS_tables
#SBATCH --time=6:00:00
#SBATCH --mem=6gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/PGS_tables_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/PGS_tables_%a_%A.err

module load gcc/12.1.0
module load R/4.3.1

Rscript /gpfs/data/ukb-share/dahl/ophelia/hairpin/code/tables.R
