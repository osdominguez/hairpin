#!/usr/bin/env bash

#SBATCH --job-name=Hairpin_PGS
#SBATCH --time=12:00:00
#SBATCH --mem=10gb
#SBATCH --output=/home/osdominguez/output/hairpin_PGS/hairpin_PGS_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_PGS/hairpin_PGS_%a_%A.err
#SBATCH --array=1-12600%300

module load gcc/12.1.0
module load R/4.3.1

COV_FILE=/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno
OUT_DIR=/scratch/osdominguez/prs_hairpin_outputs

TXT_PATH=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/combinations.txt
PHENO_DIR="/gpfs/data/ukb-share/extracted_phenotypes"
BED_DIR="/scratch/osdominguez/bgen_files/white_brit_unrelated"
SUM_DIR="/gpfs/data/ukb-share/dahl/ophelia/hairpin/sum_stats/nodups_sumstats"

readarray -t parms < <(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR == row {for(i=1; i<=NF; i++) print $i}' "${TXT_PATH}")

as=${parms[0]}
pc_n=${parms[1]}
p_value=${parms[2]}
sum=${parms[3]}
binary=${parms[4]}
ppath=${parms[5]}
scol=${parms[6]}
id=${parms[7]}
chr=${parms[8]}
bp=${parms[9]}
a1=${parms[10]}
a2=${parms[11]}
pcol=${parms[12]}
stat=${parms[13]}
pname=${parms[14]}
maf=${parms[15]}
dir=${parms[16]}
bfile=${parms[17]}

[[ ! /scratch/osdominguez/prs_hairpin_outputs/${dir}/${pname}_prs_${as}${pc_n}pc_${p_value}pval.best ]] && { exit 1; }

if [[ ${as} == "as_" ]]; then
    COVS=FID,IID,X31.0.0,X34.0.0,X54.0.0,X22000-0.0,age2
else
    COVS=FID,IID,X31.0.0,X34.0.0,X54.0.0,age2
fi

i=1
while [ $i -le ${pc_n} ]
do
    COVS+=",X22009.0.${i}"
    i=$(( $i + 1 ))
done

echo ${#parms[@]}

#really low p-values go fast so 12 hours. Be careful for some summary stats as they filter out some snps, we want ALL snps since that's what hairpin wants.
Rscript /gpfs/data/ukb-share/dahl/ophelia/PRSice.R \
    --prsice /gpfs/data/ukb-share/dahl/ophelia/PRSice_linux \
    --base ${SUM_DIR}/${sum} \
    --target ${BED_DIR}/${bfile} \
    --binary-target ${binary} \
    --pheno ${PHENO_DIR}/${ppath}  \
    --cov ${COV_FILE} \
    --cov-col ${COVS} \
    --clump-p 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --stat ${scol} \
    --snp ${id} \
    --chr ${chr} \
    --bp ${bp} \
    --A1 ${a1} \
    --A2 ${a2} \
    --pvalue ${pcol} \
    --fastscore \
    --bar-levels ${p_value} \
    --no-full \
    --base-info INFO:0.8 \
    --base-maf ${maf} \
    --"${stat}" \
    --out ${OUT_DIR}/${dir}/${pname}_prs_${as}${pc_n}pc_${p_value}pval

chgrp cri-ukb_share ${OUT_DIR}/${dir}/*
chmod g+rx ${OUT_DIR}/${dir}/*