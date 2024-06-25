#!/usr/bin/env bash

#SBATCH --job-name=hairpin_sub
#SBATCH --time=00:10:00
#SBATCH --mem=500mb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/hairpin_sub_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/hairpin_sub_%a_%A.err

HEIGHT_PATH="/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno"
BMI_PATH="/gpfs/data/ukb-share/extracted_phenotypes/BMI/BMI674178.pheno"
LDL_PATH="/gpfs/data/ukb-share/extracted_phenotypes/LDL/LDL674178.pheno"

sbatch hairpin_sub.sh ${HEIGHT_PATH} height T 100
sbatch hairpin_sub.sh ${BMI_PATH} BMI T 100
sbatch hairpin_sub.sh ${LDL_PATH} LDL T 100