#!/usr/bin/env bash

#SBATCH --job-name=hairpin_master
#SBATCH --time=00:10:00
#SBATCH --mem=500mb
#SBATCH --output=/home/osdominguez/output/hairpin_analysis/hairpin_sub_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_analysis/hairpin_sub_%a_%A.err

HEIGHT_PATH="/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno"
BMI_PATH="/gpfs/data/ukb-share/extracted_phenotypes/BMI/BMI674178.pheno"
LDL_PATH="/gpfs/data/ukb-share/extracted_phenotypes/LDL/LDL674178.pheno"
CODE_PATH="/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions/hairpin_sub.sh"

# scripts   phen_path   Pheno name(for file)   ID as it will show in R   bootstrap?   n-boots

sbatch ${CODE_PATH} ${HEIGHT_PATH} "height" "X50" "T" 100
#sbatch ${CODE_PATH} ${BMI_PATH} BMI T 100
#sbatch ${CODE_PATH} ${LDL_PATH} LDL T 100