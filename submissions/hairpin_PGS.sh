#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_PGS
#SBATCH --time=12:00:00
#SBATCH --mem=10gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/hairpin_PGS_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/hairpin_PGS_%a_%A.err

module load gcc/12.1.0
module load R/4.3.1

COV_FILE=/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC_age.pheno
OUT_DIR=/scratch/osdominguez/prs_hairpin_outputs
PHEN_DIR=/gpfs/data/ukb-share/extracted_phenotypes

COVS=FID,IID,31-0.0,34-0.0

GEN_FILE=$1
SUM_FILE=$2 
IS_BIN=$3 
PHEN_FILE=$4
COL_STAT=$5
COL_SNPID=$6 
COL_CHR=$7
COL_BP=$8 
COL_A1=$9
COL_A2=${10} 
COL_P=${11} 
PGS_THRESH=${12} 
STAT_TYPE=${13}
PHEN_NAME=${14} 
PC_N=${15}
DIR=${16}
MAF=${17} 

i=1
while [ $i -le $PC_N ]
do
  COVS+=",22009-0.${i}"
  i=$(( $i + 1 ))
done

#really low p-values go fast so 12 hours. Be careful for some summary stats as they filter out some snps, we want ALL snps since that's what hairpin wants.
Rscript /gpfs/data/ukb-share/dahl/ophelia/PRSice.R \
    --prsice /gpfs/data/ukb-share/dahl/ophelia/PRSice_linux \
    --base ${SUM_FILE} \
    --target ${GEN_FILE} \
    --binary-target ${IS_BIN} \
    --pheno ${PHEN_DIR}/${PHEN_FILE}  \
    --cov ${COV_FILE} \
    --cov-col ${COVS} \
    --clump-p 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --stat ${COL_STAT} \
    --snp ${COL_SNPID} \
    --chr ${COL_CHR} \
    --bp ${COL_BP} \
    --A1 ${COL_A1} \
    --A2 ${COL_A2} \
    --pvalue ${COL_P} \
    --fastscore \
    --bar-levels ${PGS_THRESH} \
    --no-full \
    --base-info INFO:0.8 \
    --base-maf ${MAF} \
    --"${STAT_TYPE}" \
    --out ${OUT_DIR}/${DIR}/${PHEN_NAME}_prs_${PC_N}pc_${PGS_THRESH}pval

chgrp cri-ukb_share ${OUT_DIR}/${DIR}/*
chmod g+rx ${OUT_DIR}/${DIR}/*
