#!/usr/bin/env bash

#SBATCH --job-name=hairpin_full
#SBATCH --time=24:00:00
#SBATCH --mem=1gb
#SBATCH --output=/home/osdominguez/output/hairpin_full_%a_%A.out
#SBATCH --error=/home/osdominguez/output/hairpin_full_%a_%A.err

txt_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files
sh_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions
sum_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/sum_stats

# need some way of specifying bed files
    # we ahve this with the one file having prefix, bedfile name, and the other stuff right?
    # Also REALLY need to have checks for files existing in this

new=${1}
if [[ ${new} ]]; then
    #Things we need with variable assignments:
    #new=0 # run with new phen/file?
    
    # example run for adding EA3:
    # sbatch workflow.sh 1 \
    #    EA3 \
    #    EA3 \
    #    /gpfs/data/ukb-share/extracted_phenotypes/EA3/EA3.pheno \
    #    EA4_additive_excl_23andMe.txt \
    #    F \
    #    EA3/EA3.pheno \
    #    Beta \
    #    rsID \
    #    Chr \
    #    not_avail \
    #    Effect_allele \
    #    Other_allele \
    #    P \
    #    beta \
    #    EAF_HRC:0.05
  

    phen_name=${2} # phen name
    phen_id=${3} # phen ID (as it will show in .pheno file)
    phen_path=${4}
    sum_file=${5}
    binary=${6}
    phen_rel_path=${7}
    stat_col=${8}
    ID_col=${9}
    chr_col=${10}
    bp_col=${11}
    a1=${12}
    a2=${13}
    pcol=${14}
    stat=${15}
    maf=${16}

    pinfo_1="${phen_path} ${phen_name} ${phen_id}"
    pinfo_2="${sum_file} ${binary} {phen_rel_path} ${stat_col} ${ID_col} ${chr_col} ${bp_col} ${a1} ${a2} ${pcol} ${stat} ${phen_name} ${maf}"

    # Add new run specs
    echo ${pinfo_1} >> ${txt_dir}/hairpin_phens.txt
    echo ${phen_name} >> ${txt_dir}/phens.txt
    echo ${pinfo_2} >> ${txt_dir}/PGS_phens.txt
    # For different bed files, ignore for now
    echo new_bed_info >> ${txt_dir}/sep.txt
fi

#Format the txt files
sbatch --wait ${txt_dir}/format_hairpin_txt.sh
sbatch --wait ${txt_dir}/format_lienarity_txt.sh
sbatch --wait ${txt_dir}/format_PGS_txt.sh

# Run hairpin workflow
sbatch --wait ${sh_dir}/hairpin_PGS_master.sh
sbatch --wait ${sh_dir}/tables_master.sh
sbatch --wait ${sh_dir}/hairpin_master.sh
sbatch --wait ${sh_dir}/merge_master.sh
sbatch --wait ${sh_dir}/linearity.sh
