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

reform_merge () {
    
    input_file=/gpfs/data/ukb-share/dahl/ophelia/hairpin/submissions/merge_master.sh

    pattern='*)'  # Pattern to match '*)' with optional leading whitespace

    # Find line number of the pattern '*)'
    line_number=$(grep -n "$pattern" "$input_file" | cut -d':' -f1)

    if [[ -n "$line_number" ]]; then
        # Extract the number to update (assuming it's the line before '*)')
        prev_line_number=$((line_number - 5))
        number=$(sed -n "${prev_line_number}p" "$input_file" | grep -oP '[0-9]+')

        if [[ -n "$number" ]]; then
            # Increment the number
            new_number=$((number + 1))

            # Construct the new line to add
            new_line="\  ${new_number})"
            new_line+="\n    echo \"${phen_name}"\"
            new_line+="\n    Rscript \${CODE} ${phen_path} ${phen_name} ${phen_id} 100"
            new_line+="\n    ;;\n"

            # Insert the new line into the script file
            sed -i "${line_number}i ${new_line}" "$input_file"

            echo "New line added successfully."
        else
            echo "Error: Number extraction failed."
        fi
    else
        echo "Error: Pattern not found in the script file."
    fi

}


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

reform_merge
reform_array "${txt_dir}/combinations.txt" "${sub_dir}/hairpin_PGS_master.sh" 300
reform_array "${txt_dir}/hairpin.txt" "${sub_dir}/hairpin_master.sh" 101
reform_array "${txt_dir}/phens.txt" "${sub_dir}/merge_master.sh" 1
reform_array "${txt_dir}/linearity.txt" "${sub_dir}/linearity.sh" 40
reform_array "${txt_dir}/tables.txt" "${sub_dir}/tables.sh" 2