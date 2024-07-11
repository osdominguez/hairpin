#!/usr/bin/env bash

#SBATCH --job-name=hairpin_master
#SBATCH --time=12:00:00
#SBATCH --mem=20gb
#SBATCH --output=/home/osdominguez/output/hairpin_analysis/hairpin_sub_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_analysis/hairpin_sub_%a_%A.err
#SBATCH --array=1-1414%101

TXT_PATH=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/hairpin.txt
CODE="/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/hairpin.R"
OUT_DIR="/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting"
TMP_DIR="/scratch/osdominguez/temp_boot"

readarray -t parms < <(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR == row {for(i=1; i<=NF; i++) print $i}' "${TXT_PATH}")

module load gcc/12.1.0
module load R/4.3.1

[[ ! ${OUT_DIR}/${parms[2]}_${parms[0]}_base.table ]] && { exit 1; }
[[ ! ${TMP_DIR}/${parms[2]}_${parms[0]}_bootstrap_${parms[5]}.table ]] && { exit 1; }

Rscript ${CODE} ${parms[1]} ${parms[2]} ${parms[3]} ${parms[4]} ${parms[5]} ${parms[0]}