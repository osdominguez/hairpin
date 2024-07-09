#!/usr/bin/env bash

txt_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files
sh_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions
sum_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/sum_stats

# need some way of specifying bed files
    # we ahve this with the one file having prefix, bedfile name, and the other stuff right?
    # Also REALLY need to have checks for files existing in this

#Things we need with variable assignments:
new=0 # run with new phen/file?
phen_name=height # phen name
phen_id=X50 # phen ID (as it will show in .pheno file)
phen_path=/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno
sum_file=EA4_additive_excl_23andMe.txt
binary=F
phen_rel_path=EA3/EA3.pheno
stat_col=Beta
ID_col=rsID
chr_col=Chr
bp_col=not_avail
a1=Effect_allele
a2=Other_allele
pcol=P
stat=beta
maf=EAF_HRC:0.05

pinfo_1="${phen_path} ${phen_name} ${phen_id}"
pinfo_2="${sum_file} ${binary} {phen_rel_path} ${stat_col} ${ID_col} ${chr_col} ${bp_col} ${a1} ${a2} ${pcol} ${stat} ${phen_name} ${maf}"

# Add new run specs
echo ${pinfo_1} >> ${txt_dir}/hairpin_phens.txt
echo ${phen_name} >> ${txt_dir}/phens.txt
echo ${pinfo_2} >> ${txt_dir}/PGS_phens.txt
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
