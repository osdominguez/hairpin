#!/usr/bin/awk -f

# Define variables to hold file names
BEGIN {
    pcs ="/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/pcs.txt"
    pvals ="/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/pvalues.txt"
    phens ="/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/PGSphens.txt"
    
    # Skip the first three arguments (which are file names)
    ARGC = 4
}

# Read pcs into an array
FNR == NR {
    pcs_lines[NR] = $0
    pcs_count = NR
    next
}

# Read pvals into an array
FNR > pcs_count && FNR <= pcs_count + pvals_count {
    pvals_lines[FNR - pcs_count] = $0
    pvals_count = FNR - pcs_count
    next
}

# Read phens into an array
FNR > pcs_count + pvals_count && FNR == NR {
    for (i = 1; i <= NF; i++) {
        phens_lines[FNR - pcs_count - pvals_count, i] = $i
    }
    phens_count = FNR - pcs_count - pvals_count
    next
}

# After reading all files, generate combinations
END {
    for (i = 1; i <= pcs_count; i++) {
        for (j = 1; j <= pvals_count; j++) {
            for (k = 1; k <= phens_count; k++) {
                printf "%s %s", pcs_lines[i], pvals_lines[j]
                for (l = 1; l <= NF; l++) {
                    printf " %s", phens_lines[k, l]
                }
                printf "\n"
            }
        }
    }
}