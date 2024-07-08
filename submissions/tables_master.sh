#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_PGS_tables_master
#SBATCH --time=00:30:00
#SBATCH --mem=250mb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/PGS_tables_master_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/PGS_tables_master_%a_%A.err

SCRIPT_PATH='/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions/individual/tables_sub.sh'

sbatch ${SCRIPT_PATH} "height"
sbatch ${SCRIPT_PATH} "BMI"
sbatch ${SCRIPT_PATH} "LDL"
sbatch ${SCRIPT_PATH} "WHRadjBMI"
sbatch ${SCRIPT_PATH} "EA3"
sbatch ${SCRIPT_PATH} "EA4"