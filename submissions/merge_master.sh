#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_merge
#SBATCH --time=12:00:00
#SBATCH --mem=75gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/linearity_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/linearity_%a_%A.err
#SBATCH --array=1-6%1

module load gcc/12.1.0
module load R/4.3.1

CODE="/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/merge_CI.R"

case ${SLURM_ARRAY_TASK_ID}} in

  1)
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno height X50 100
    ;;

  2)
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/BMI/BMI674178.pheno BMI X21001 100
    ;;

  3)
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/LDL/LDL674178.pheno LDL X30780 100
    ;;

  4)
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/WHRadjBMI/WHRadjBMI.pheno WHRadjBMI WHRadjBMI 100
    ;;
  5)
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/EA3/EA3.pheno EA3 EA3 100
    ;;

  6)
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/EA4/EA4.pheno EA4 EA4 100
    ;;
    
  *)
    echo -n "unknown"
    ;;