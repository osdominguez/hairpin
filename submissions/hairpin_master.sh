#!/usr/bin/env bash

#SBATCH --job-name=hairpin_master
#SBATCH --time=00:10:00
#SBATCH --mem=500mb
#SBATCH --output=/home/osdominguez/output/hairpin_analysis/hairpin_sub_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_analysis/hairpin_sub_%a_%A.err

HEIGHT_PATH="/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno"
BMI_PATH="/gpfs/data/ukb-share/extracted_phenotypes/BMI/BMI674178.pheno"
LDL_PATH="/gpfs/data/ukb-share/extracted_phenotypes/LDL/LDL674178.pheno"
WHRadjBMI_PATH="/gpfs/data/ukb-share/extracted_phenotypes/WHRadjBMI/WHRadjBMI.pheno"
EA3_PATH="/gpfs/data/ukb-share/extracted_phenotypes/EA3/EA3.pheno"
EA4_PATH="/gpfs/data/ukb-share/extracted_phenotypes/EA4/EA4.pheno"

CODE_PATH="/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions/individual/hairpin_sub.sh"

# scripts   phen_path   Pheno name(for file)   ID as it will show in R   bootstrap?   n-boots

N_BOOTS=$1

hairpin_func() {
    sbatch ${CODE_PATH} $1 $2 $3 "F" ${N_boots}
    for ((i = 1 ; i < ${N_BOOTS} ; i++)); do
        sbatch ${CODE_PATH} $1 $2 $3 "T" $i
    done
    sbatch --wait ${CODE_PATH} $1 $2 $3 "T" $i
    #sbatch ${new script}  all the parameters + n_boots
}

sbatch ${CODE_PATH} ${HEIGHT_PATH} "height" "X50" "T" ${N_BOOTS}
sbatch ${CODE_PATH} ${BMI_PATH} "BMI" "X21001" "T" 100
sbatch ${CODE_PATH} ${LDL_PATH} "LDL" "X30780" "T" 100
sbatch ${CODE_PATH} ${WHRadjBMI_PATH} "WHRadjBMI" "WHRadjBMI" "T" 100
sbatch ${CODE_PATH} ${EA3_PATH} "EA3" "EA3" "T" 100
sbatch ${CODE_PATH} ${EA4_PATH} "EA4" "EA4" "T" 100
