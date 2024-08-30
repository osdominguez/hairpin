#!/usr/bin/env bash

#SBATCH --job-name=hairpin_full
#SBATCH --time=24:00:00
#SBATCH --mem=10mb
#SBATCH --output=/home/osdominguez/output/hairpin_full_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_full_%a_%A.err

hairpin_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin

txt_dir=${hairpin_dir}/txt_files
sh_dir=${hairpin_dir}/submissions
sum_dir=${hairpin_dir}/sum_stats

# need some way of specifying bed files
    # we ahve this with the one file having prefix, bedfile name, and the other stuff right?
    # Also REALLY need to have checks for files existing in this

new=${1}

if [[ ${new} == 1 ]]; then
  
    # put not_avail if the column is not available
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

    sbatch --wait ${txt_dir}/new_phen_format.sh \
        ${new} \
        ${phen_name} \
        ${phen_id} \
        ${phen_path} \
        ${sum_file} \
        ${binary} \
        ${phen_rel_path} \
        ${stat_col} \
        ${ID_col} \
        ${chr_col} \
        ${bp_col} \
        ${a1} \
        ${a2} \
        ${pcol} \
        ${stat} \
        ${maf}

else
    sbatch --wait ${txt_dir}/new_phen_format.sh ${new}
fi

# Run hairpin workflow
echo "running hairpin PGS..."
sbatch --wait ${sh_dir}/hairpin_PGS.sh
echo "finished running hairpin PGS"

echo "running tabling for hairpin PGS..."
sbatch --wait ${sh_dir}/tables.sh
echo "finished tabling hairpin PGS"

echo "running hairpin..."
sbatch --wait ${sh_dir}/hairpin.sh
echo "finished running hairpin"

echo "combining hairpin bootstaps..."
sbatch --wait ${sh_dir}/merge.sh
echo "finished combining hairpin bootstraps"

echo "testing linearity..."
sbatch --wait ${sh_dir}/linearity.sh
echo "finished testing linearity"

echo "hairpin workflow completed"