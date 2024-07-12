#!/usr/bin/env bash

#SBATCH --job-name=reformat_scripts
#SBATCH --time=00:10:00
#SBATCH --mem=250mb
#SBATCH --output=/home/osdominguez/output/formatting/reformat_%a_%A.out
#SBATCH --error=/home/osdominguez/output/format/reformat_%A.err

phen_name=${1}
phen_path=${2}
phen_id=${3}

reform_merge () {
    
    input_file=/gpfs/data/ukb-share/dahl/ophelia/hairpin/workflow_test/workflow.sh

    # Search for the line containing "#SBATCH --array=1-<number>%1"
    line=$(grep "#SBATCH --array=1-[0-9]%1" "$input_file")

    if [[ -n "$line" ]]; then
        # Extract the number after the dash
        number=$(echo "$line" | grep -oP '(?<=1-)[0-9]+')

        if [[ -n "$number" ]]; then
            # Increment the number
            new_number=$((number + 1))

            # Replace the old number with the new number in the line
            new_line=$(echo "$line" | sed "s/1-$number%1/1-$new_number%1/")

            # Replace the line in the file
            sed -i "s|$line|$new_line|g" "$input_file"

            echo "Number replaced successfully."
        else
            echo "Error: Number extraction failed."
        fi
    else
        echo "Error: Line not found in the input file."
    fi

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
            new_line+="\n  echo \"days_phys\""
            new_line+="\n    Rscript \${CODE} ${phen_path} ${phen_name} ${phen_id}"
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


