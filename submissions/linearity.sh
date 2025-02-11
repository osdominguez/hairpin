#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_linearity
#SBATCH --time=12:00:00
#SBATCH --mem=20gb
#SBATCH --output=/home/osdominguez/output/linearity/linearity_%A_%a.out
#SBATCH --error=/home/osdominguez/output/linearity/linearity_%A_%a.err
#SBATCH --array=1-60%5

TXT_PATH=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/linearity.txt

[[ ( ${SLURM_ARRAY_TASK_ID} -gt $(awk 'END{print NR}' ${TXT_PATH}) ) ]] && { echo "slurm ID greater"; exit 0; }

module load gcc/12.1.0
module load R/4.3.1

readarray -t parms < <(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR == row {for(i=1; i<=NF; i++) print $i}' "${TXT_PATH}")

Rscript /gpfs/data/ukb-share/dahl/ophelia/hairpin/code/linearity.R ${parms[0]} ${parms[1]} ${parms[2]} 