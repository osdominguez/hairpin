#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_PGS_tables_master
#SBATCH --time=04:00:00
#SBATCH --mem=20gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/PGS_tables_master_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/PGS_tables_master_%a_%A.err
#SBATCH --array=1-8%2

SCRIPT=/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/tables.R
TXT_PATH=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/tables.txt

readarray -t parms < <(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR == row {for(i=1; i<=NF; i++) print $i}' "${TXT_PATH}")

module load gcc/12.1.0
module load R/4.3.1

Rscript ${SCRIPT} ${parms[0]} ${parms[1]}
