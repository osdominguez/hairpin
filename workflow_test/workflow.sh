#!/usr/bin/env bash

txt_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files
sh_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions

# need some way of specifying bed files
    # we ahve this with the one file having prefix, bedfile name, and the other stuff right?
    # Also REALLY need to have checks for files existing in this

#Things we need with variable assignments:
new=0 # run with new phen/file?
phen_name=height # phen name
phen_id=X50 # phen ID (as it will show in .pheno file)
phen_path=/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno
#EA4_additive_excl_23andMe.txt F EA3/EA3.pheno Beta rsID Chr not_avail Effect_allele Other_allele P beta EA3 EAF_HRC:0.05

# Add new run specs
echo new_phen_info1 >> ${txt_dir}/hairpin_phens.txt
echo new_phen_info2 >> ${txt_dir}/phens.txt
echo new_phen_info3 >> ${txt_dir}/PGS_phens.txt
# For different bed files, ignore for now
echo new_bed_info >> ${txt_dir}/sep.txt

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
