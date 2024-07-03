#!/usr/bin/env bash

#SBATCH --job-name=hairpin_master
#SBATCH --time=24:00:00
#SBATCH --mem=250mb
#SBATCH --output=/home/osdominguez/output/hairpin_analysis/hairpin_sub_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_analysis/hairpin_sub_%a_%A.err

HEIGHT_PATH="/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno"
BMI_PATH="/gpfs/data/ukb-share/extracted_phenotypes/BMI/BMI674178.pheno"
LDL_PATH="/gpfs/data/ukb-share/extracted_phenotypes/LDL/LDL674178.pheno"
WHRadjBMI_PATH="/gpfs/data/ukb-share/extracted_phenotypes/WHRadjBMI/WHRadjBMI.pheno"
EA3_PATH="/gpfs/data/ukb-share/extracted_phenotypes/EA3/EA3.pheno"
EA4_PATH="/gpfs/data/ukb-share/extracted_phenotypes/EA4/EA4.pheno"

HAIRPIN_PATH="/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions/individual/hairpin_sub_test.sh"
MERGE_PATH="/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions/individual/CIMerge_sub.sh"

TEMP_DIR="/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/temp"

# scripts   phen_path   Pheno name(for file)   ID as it will show in R   bootstrap?   n-boots

N_BOOTS=$1

hairpin_func() {
    sbatch ${HAIRPIN_PATH} $1 $2 $3 "F" ${N_boots}
    for ((i = 1 ; i < ${N_BOOTS} ; i++)); do
        sbatch ${HAIRPIN_PATH} $1 $2 $3 "T" $i
    done
    sbatch --wait ${HAIRPIN_PATH} $1 $2 $3 "T" $i
    sbatch ${MERGE_PATH} $1 $2 $3 $i
    rm ${TEMP_DIR}/${2}*
}

hairpin_func ${HEIGHT_PATH} "height" "X50"
sbatch ${CODE_PATH} ${HEIGHT_PATH} "height" "X50"
sbatch ${CODE_PATH} ${BMI_PATH} "BMI" "X21001"
sbatch ${CODE_PATH} ${LDL_PATH} "LDL" "X30780"
sbatch ${CODE_PATH} ${WHRadjBMI_PATH} "WHRadjBMI" "WHRadjBMI"
sbatch ${CODE_PATH} ${EA3_PATH} "EA3" "EA3"
sbatch ${CODE_PATH} ${EA4_PATH} "EA4" "EA4"
