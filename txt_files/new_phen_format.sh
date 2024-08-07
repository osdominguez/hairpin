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
new=${1}
phen_name=${2} # phen name
phen_id=${3} # phen ID (as it will show in .pheno file)
phen_path=${4} # Full path to pheno
sum_file=${5} # Name of the summary statistic file
binary=${6} # Is phenotype binary?
phen_rel_path=${7} # Relative path of the phenotype
stat_col=${8} # Name of the ranking stat column
ID_col=${9} # snpID column
chr_col=${10} # CHR column
bp_col=${11} # Base pair
a1=${12} # effect allele
a2=${13} # other allele
pcol=${14} # pvalue column
stat=${15} # stat type (should be beta)
maf=${16} # allele frequency (think this is effect allele) then put :0.05

if [[ ${new} == 1 ]]; then
    echo "adding new phenotype to hairpin workflow: ${phen_name}..."

    pinfo_1="${phen_path} ${phen_name} ${phen_id}"
    pinfo_2="${sum_file} ${binary} ${phen_rel_path} ${stat_col} ${ID_col} ${chr_col} ${bp_col} ${a1} ${a2} ${pcol} ${stat} ${phen_name} ${maf}"

    # Add new run specs
    printf "\n${pinfo_1}" >> ${txt_dir}/hairpin_phens.txt
    printf "\n${phen_name}" >> ${txt_dir}/phens.txt
    printf "\n${pinfo_2}" >> ${txt_dir}/PGSphens.txt
    printf '\nsbatch ${SCRIPT_PATH}'" \"${phen_name}\" " >> ${sh_dir}/tables_master.sh
fi

echo "formatting txt files..."
#Format the txt files
rm ${txt_dir}/combinations.txt
sbatch --wait ${txt_dir}/format_PGS_txt.sh
rm ${txt_dir}/hairpin.txt
sbatch --wait ${txt_dir}/format_hairpin_txt.sh
rm ${txt_dir}/linearity.txt
sbatch --wait ${txt_dir}/format_linearity_txt.sh
echo "successfully formatted txt files"

if [[ ${new} == 1 ]]; then
    echo "reformatting scripts..."
    sbatch --wait ${sh_dir}/reformat_scripts.sh ${phen_name} ${phen_path} ${phen_id}

    echo "successfully formatted scripts"
fi

# For different bed files, ignore for now
#echo new_bed_info >> ${txt_dir}/sep.txt
echo "succesfully added phenotype"