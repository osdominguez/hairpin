#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_merge
#SBATCH --time=12:00:00
#SBATCH --mem=40gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/merge_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/merge_%a_%A.err
#SBATCH --array=1-20%1

module load gcc/12.1.0
module load R/4.3.1

CODE="/gpfs/data/ukb-share/dahl/ophelia/hairpin/code/merge_CI.R"

case ${SLURM_ARRAY_TASK_ID} in

  1)
    echo "GPpsy_imputed WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/GPpsy/GPpsy_imputed.pheno GPpsy_imputed GPpsy 100 WBRT
    ;;

  2)
    echo "ICD10Dep_imputed WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/ICD10Dep/ICD10Dep_imputed.pheno ICD10Dep_imputed ICD10Dep 100 WBRT
    ;;

  3)
    echo "LifetimeMDD_imputed WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/LifetimeMDD/LifetimeMDD_imputed.pheno LifetimeMDD_imputed LifetimeMDD 100 WBRT
    ;;

  4)
    echo "Psypsy_imputed WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/Psypsy/Psypsy_imputed.pheno Psypsy_imputed Psypsy 100 WBRT
    ;;

  5)
    echo "DepAll_imputed WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/DepAll/DepAll_imputed.pheno DepAll_imputed DepAll 100 WBRT
    ;;

  6)
    echo "DepAll WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/DepAll/DepAll_observed.pheno DepAll DepAll 100 WBRT 
    ;;

  7)
    echo "Psypsy WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/Psypsy/Psypsy_observed.pheno Psypsy Psypsy 100 WBRT 
    ;;

  8)
    echo "LifetimeMDD WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/LifetimeMDD/LifetimeMDD_observed.pheno LifetimeMDD LifetimeMDD 100 WBRT
    ;;

  9)
    echo "ICD10Dep WBRT"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/ICD10Dep/ICD10Dep_observed.pheno ICD10Dep ICD10Dep 100 WBRT
    ;;

  10)
    echo "GPpsy WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/GPpsy/GPpsy_observed.pheno GPpsy GPpsy 100 WEUR 
    ;;
  
  11)
    echo "GPpsy_imputed WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/GPpsy/GPpsy_imputed.pheno GPpsy_imputed GPpsy 100 WEUR
    ;;

  12)
    echo "ICD10Dep_imputed WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/ICD10Dep/ICD10Dep_imputed.pheno ICD10Dep_imputed ICD10Dep 100 WEUR
    ;;

  13)
    echo "LifetimeMDD_imputed WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/LifetimeMDD/LifetimeMDD_imputed.pheno LifetimeMDD_imputed LifetimeMDD 100 WEUR
    ;;

  14)
    echo "Psypsy_imputed WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/Psypsy/Psypsy_imputed.pheno Psypsy_imputed Psypsy 100 WEUR
    ;;

  15)
    echo "DepAll_imputed WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/DepAll/DepAll_imputed.pheno DepAll_imputed DepAll 100 WEUR
    ;;

  16)
    echo "DepAll WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/DepAll/DepAll_observed.pheno DepAll DepAll 100 WEUR
    ;;

  17)
    echo "Psypsy WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/Psypsy/Psypsy_observed.pheno Psypsy Psypsy 100 WEUR
    ;;

  18)
    echo "LifetimeMDD WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/LifetimeMDD/LifetimeMDD_observed.pheno LifetimeMDD LifetimeMDD 100 WEUR
    ;;

  19)
    echo "ICD10Dep WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/ICD10Dep/ICD10Dep_observed.pheno ICD10Dep ICD10Dep 100 WEUR
    ;;

  20)
    echo "GPpsy WEUR"
    Rscript ${CODE} /gpfs/data/ukb-share/extracted_phenotypes/GPpsy/GPpsy_observed.pheno GPpsy GPpsy 100 WEUR
    ;;
  
  *)
    echo -n "unknown"
    exit 0
    ;;

esac