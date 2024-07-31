#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_merge
#SBATCH --time=12:00:00
#SBATCH --mem=75gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/merge_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/merge_%a_%A.err
#SBATCH --array=1-10%1

module load gcc/12.1.0
module load R/4.3.1

CODE="/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/merge_CI.R"

case ${SLURM_ARRAY_TASK_ID} in

  1)
    echo "height"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno height X50 100
    ;;

  2)
    echo "BMI"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/BMI/BMI674178.pheno BMI X21001 100
    ;;

  3)
    echo "LDL"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/LDL/LDL674178.pheno LDL X30780 100
    ;;

  4)
    echo "WHRadjBMI"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/WHRadjBMI/WHRadjBMI.pheno WHRadjBMI WHRadjBMI 100
    ;;
  5)
    echo "EA3"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/EA3/EA3.pheno EA3 EA3 100
    ;;

  6)
    echo "EA4"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/EA4/EA4.pheno EA4 EA4 100
    ;;

  7)
    echo "days_phys"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/days_phys/days_phys.pheno days_phys X884.0.0 100
    ;;

  8)
    echo "T2D"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/T2D/T2D.pheno T2D X2443.0.0 100
    ;;

  9)
    echo "height_on_EA4"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/EA4/EA4.pheno height_on_EA4 EA4 100
    ;;

  10)
    echo "height_on_BMI"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/BMI/BMI674178.pheno height_on_BMI X21001 100
    ;;

  *)
    echo -n "unknown"
    exit 0
    ;;
esac