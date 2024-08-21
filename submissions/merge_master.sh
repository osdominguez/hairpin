#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_merge
#SBATCH --time=12:00:00
#SBATCH --mem=40gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/merge_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/merge_%a_%A.err
#SBATCH --array=1-11%1

module load gcc/12.1.0
module load R/4.3.1

CODE="/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/merge_CI.R"

case ${SLURM_ARRAY_TASK_ID} in

  1)
    echo "EA4_NUKB_WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/EA4/EA4.pheno EA4_NUKB EA4 100 WBRT
    ;;

  2)
    echo "EA4_NUKB_WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/EA4/EA4.pheno EA4_NUKB EA4 100 WEUR
    ;;

  3)
    echo "EA3_NUKB_WBRT"
    Rscript ${CODE} //gpfs/data/ukb-share/extracted_phenotypes/EA3/EA3.pheno EA3_NUKB EA3 100 WBRT
    ;;

  4)
    echo "EA3_NUKB_WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/EA3/EA3.pheno EA3_NUKB EA3 100 WEUR
    ;;

  5)
    echo "GPpsy WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/GPpsy/GPpsy_observed.pheno GPpsy GPpsy 100 WBRT
    ;;

  6)
    echo "GPpsy WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/GPpsy/GPpsy_observed.pheno GPpsy GPpsy 100 WEUR
    ;;

  7)
    echo "DepAll"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/DepAll/DepAll_observed.pheno DepAll DepAll 100 WBRT
    ;;

  8)
    echo "ICD10Dep"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/ICD10Dep/ICD10Dep_observed.pheno ICD10Dep ICD10Dep 100 WBRT
    ;;

  9)
    echo "LifetimeMDD"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/LifetimeMDD/LifetimeMDD_observed.pheno LifetimeMDD LifetimeMDD 100 WBRT
    ;;

  10)
    echo "Psypsy"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/Psypsy/Psypsy_observed.pheno Psypsy Psypsy 100 WBRT
    ;;

  11)
    echo "DepAll"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/DepAll/DepAll_observed.pheno DepAll DepAll 100 WBRT
    ;;

  *)
    echo -n "unknown"
    exit 0
    ;;
esac