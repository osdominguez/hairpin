#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_sim
#SBATCH --time=20:00:00
#SBATCH --mem=60gb
#SBATCH --output=/home/osdominguez/output/hairpin_sim_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_sim_%a_%A.err

SCRIPT=/gpfs/data/ukb-share/dahl/ophelia/hairpin/sims/sim.R

module load gcc/12.1.0
module load R/4.3.1

Rscript ${SCRIPT}