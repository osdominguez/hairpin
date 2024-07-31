#!/usr/bin/env bash

#SBATCH --job-name=hairpin_reformat
#SBATCH --time=00:10:00
#SBATCH --mem=1gb
#SBATCH --output=/home/osdominguez/output/formatting/hairpin_reformat_%a_%A.out
#SBATCH --error=/home/osdominguez/output/formatting/hairpin_reformat_%a_%A.err

hairpin_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin

txt_dir=${hairpin_dir}/txt_files
sh_dir=${hairpin_dir}/submissions
sum_dir=${hairpin_dir}/sum_stats

# need some way of specifying bed files
# we ahve this with the one file having prefix, bedfile name, and the other stuff right?
# Also REALLY need to have checks for files existing in this

# put not_avail if the column is not available
phen_name=${1} # phen name
phen_id=${2} # phen ID (as it will show in .pheno file)
phen_path=${3} # Full path to pheno
sum_file=${4} # Name of the summary statistic file
binary=${5} # Is phenotype binary?
phen_rel_path=${6} # Relative path of the phenotype
stat_col=${7} # Name of the ranking stat column
ID_col=${8} # snpID column
chr_col=${9} # CHR column
bp_col=${10} # Base pair
a1=${11} # effect allele
a2=${12} # other allele
pcol=${13} # pvalue column
stat=${14} # stat type (should be beta)
maf=${15} # allele frequency (think this is effect allele) then put :0.05

echo "adding new phenotype to hairpin workflow: ${phen_name}..."

pinfo_1="${phen_path} ${phen_name} ${phen_id}"
pinfo_2="${sum_file} ${binary} ${phen_rel_path} ${stat_col} ${ID_col} ${chr_col} ${bp_col} ${a1} ${a2} ${pcol} ${stat} ${phen_name} ${maf}"

# Add new run specs
printf "\n${pinfo_1}" >> ${txt_dir}/hairpin_phens.txt
printf "\n${phen_name}" >> ${txt_dir}/phens.txt
printf "\n${pinfo_2}" >> ${txt_dir}/PGSphens.txt
printf '\nsbatch ${SCRIPT_PATH}'" \"${phen_name}\" " >> ${sh_dir}/tables_master.sh

echo "formatting txt files..."
#Format the txt files
rm ${txt_dir}/combinations.txt
sbatch --wait ${txt_dir}/format_PGS_txt.sh
rm ${txt_dir}/hairpin.txt
sbatch --wait ${txt_dir}/format_hairpin_txt.sh
rm ${txt_dir}/linearity.txt
sbatch --wait ${txt_dir}/format_linearity_txt.sh
echo "successfully formatted txt files"
echo "reformatting scripts..."
sbatch --wait ${sh_dir}/reformat_scripts.sh ${phen_name} ${phen_path} ${phen_id}
echo "successfully formatted scripts"
# For different bed files, ignore for now
#echo new_bed_info >> ${txt_dir}/sep.txt
echo "succesfully added phenotype"