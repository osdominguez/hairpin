#!/usr/bin/env bash

#SBATCH --job-name=test_submit_pgs
#SBATCH --time=01:00:00
#SBATCH --mem=1gb
#SBATCH --output=/home/osdominguez/output/SUBMIT_HAIRPIN_PGS_%a_%A.out
#SBATCH --error=/home/osdominguez/output/SUBMIT_HAIRPIN_PGS_%a_%A.err

TXT_DIR=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files
PHENO_DIR=/gpfs/data/ukb-share/extracted_phenotypes
BED_DIR=/scratch/osdominguez/bgen_files/white_brit_unrelated
SUM_DIR=/gpfs/data/ukb-share/dahl/ophelia/hairpin/sum_stats/nodups_sumstats
SCRIPT_PATH="/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions/hairpin_PGS.sh"

mapfile -t pvals < ${TXT_DIR}/pvalues.txt

mapfile -t PC_number < ${TXT_DIR}/pcs.txt

sub_phenotype() {

n=0

for pc_n in "${PC_number[@]}"; do

    for p_value in "${pvals[@]}"; do 
    	
      #run odd chromosome PRS
       [ ! -f /scratch/osdominguez/prs_hairpin_outputs/odd/${12}_prs_${pc_n}pc_${p_value}pval.best ] && \
        sbatch ${SCRIPT_PATH} \
          ${BED_DIR}/ukb_imp_chr0O_v3_whitebrit_unrelated_QC \
          ${SUM_DIR}/$1 \
          $2 \
          ${PHEN_DIR}/$3 \
          $4 \
          $5 \
          $6 \
          $7 \
          $8 \
          $9 \
          ${10} \
          ${p_value} \
          ${11} \
          ${12} \
          ${pc_n} \
          "odd" \
          ${13}

      #run even chromosome PRS
        [ ! -f /scratch/osdominguez/prs_hairpin_outputs/even/${12}_prs_${pc_n}pc_${p_value}pval.best ] && \
        sbatch ${SCRIPT_PATH} \
          ${BED_DIR}/ukb_imp_chr0E_v3_whitebrit_unrelated_QC \
          ${SUM_DIR}/$1 \
          $2 \
          ${PHEN_DIR}/$3 \
          $4 \
          $5 \
          $6 \
          $7 \
          $8 \
          $9 \
          ${10} \
          ${p_value} \
          ${11} \
          ${12} \
          ${pc_n} \
          "even" \
          ${13}  
      
      #run all chromosome PRS
        [ ! -f /scratch/osdominguez/prs_hairpin_outputs/all/${12}_prs_${pc_n}pc_${p_value}pval.best ] && \
        sbatch ${SCRIPT_PATH} \
          ${BED_DIR}/ukb_imp_chr0E_v3_whitebrit_unrelated_QC \
          ${SUM_DIR}/$1 \
          $2 \
          ${PHEN_DIR}/$3 \
          $4 \
          $5 \
          $6 \
          $7 \
          $8 \
          $9 \
          ${10} \
          ${p_value} \
          ${11} \
          ${12} \
          ${pc_n} \
          "all" \
          ${13}  

    done

n=$((n + 1))

done

}

#function_call	sumstat_file_name is_binary pheno_file_path (relative to extracted phehnotypes) stat_column_name snpID_column_name   chr_col_name   bp_col_name   A1_col   A2_col   p_val_col   stat   pheno_name   mafname_mafscore 

sub_phenotype GIANT_HEIGHT_new_R.csv F Height/Height674178.pheno BETA RSID CHR POS EFFECT_ALLELE OTHER_ALLELLE P beta height EFFECT_ALLELE_FREQ:0.05

sub_phenotype SNP_gwas_mc_merge_nogc.tbl.uniq F BMI/BMI674178.pheno b SNP not_avail not_avail A1 A2 p beta BMI Freq1.Hapmap:0.05 

sub_phenotype jointGwasMc_LDL.txt F LDL/LDL674178.pheno beta rsid not_avail not_avail A1 A2 P-value beta LDL Freq.A1.1000G.EUR:0.05

sub_phenotype GIANT_2015_WHRadjBMI_COMBINED_EUR.txt F WHRadjBMI/WHRadjBMI.pheno b MarkerName not_vail not_vail Allele1 Allele2 p beta WHRadjBMI FreqAllele1HapMapCEU:0.05

sub_phenotype EA4_additive_excl_23andMe.txt F EA4/EA4.pheno Beta rsID Chr not_avail Effect_allele Other_allele P beta EA4 EAF_HRC:0.05

sub_phenotype EA4_additive_excl_23andMe.txt F EA3/EA3.pheno Beta rsID Chr not_avail Effect_allele Other_allele P beta EA3 EAF_HRC:0.05