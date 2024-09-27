#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_merge
#SBATCH --time=12:00:00
#SBATCH --mem=40gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/merge_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/merge_%a_%A.err
#SBATCH --array=1-60%2

SCRIPT="/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/merge.R"
TXT_PATH=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/merge.txt

[[ ( ${SLURM_ARRAY_TASK_ID} -gt $(awk 'END{print NR}' ${TXT_PATH}) ) ]] && { echo "slurm ID greater"; exit 0; }

module load gcc/12.1.0
module load R/4.3.1

readarray -t parms < <(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR == row {for(i=1; i<=NF; i++) print $i}' "${TXT_PATH}")

Rscript ${SCRIPT} ${parms[0]} ${parms[1]} ${parms[2]} ${parms[3]}