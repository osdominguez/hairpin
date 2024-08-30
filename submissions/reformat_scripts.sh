#!/usr/bin/env bash

#SBATCH --job-name=reformat_scripts
#SBATCH --time=00:10:00
#SBATCH --mem=250mb
#SBATCH --output=/home/osdominguez/output/formatting/reformat_%a_%A.out
#SBATCH --error=/home/osdominguez/output/formatting/reformat_%a_%A.err

txt_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files
sub_dir=/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions

phen_name=${1}
phen_path=${2}
phen_id=${3}

reform_array() {

    input_file=${1}
    output_file=${2}
    max_array=${3}

    if [ ! -f "$input_file" ]; then
        echo "Error: Input file '$input_file' not found."
        exit 1
    fi

    n_row=$(awk 'END { print NR }' ${input_file})

    echo ${n_row}

    new_line="#SBATCH --array=1-${n_row}%${max_array}"

    sed -i '8s/.*/'"$new_line"'/' "$output_file"

    echo "Line 8 replaced successfully with n_row = $n_row."
}

reform_array "${txt_dir}/combinations.txt" "${sub_dir}/hairpin_PGS.sh" 200
reform_array "${txt_dir}/hairpin.txt" "${sub_dir}/hairpin.sh" 101
reform_array "${txt_dir}/linearity.txt" "${sub_dir}/linearity.sh" 5
reform_array "${txt_dir}/tables.txt" "${sub_dir}/tables.sh" 2
reform_array "${txt_dir}/merge.txt" "${sub_dir}/merge.sh" 2