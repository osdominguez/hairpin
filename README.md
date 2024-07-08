# Most of the code here was made by Shevaugn originally or very heavily inspired by Renée

# Covars path:

/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno

# Covar Fields (for PRSice) :

## Without smoking:
FID IID X31.0.0 X34.0.0 X54.0.0 X22000.0.0 X22009.0.1 X22009.0.2 X22009.0.3 X22009.0.4 X22009.0.5 X22009.0.6 X22009.0.7 X22009.0.8 X22009.0.9 X22009.0.10 X22009.0.11 X22009.0.12 X22009.0.13 X22009.0.14 X22009.0.15 X22009.0.16 X22009.0.17 X22009.0.18 X22009.0.19 X22009.0.20 X22009.0.21 X22009.0.22 X22009.0.23 X22009.0.24 X22009.0.25 X22009.0.26 X22009.0.27 X22009.0.28 X22009.0.29 X22009.0.30 X22009.0.31 X22009.0.32 X22009.0.33 X22009.0.34 X22009.0.35 X22009.0.36 X22009.0.37 X22009.0.38 X22009.0.39 X22009.0.40 age2

They have X's inf front of the numerical ones due to the script I used being written in R

# Pheno paths:

/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno

# PGS Specs:

In hairpin/txt_files there are txt files for both the pcs and pvalues. These are new line delimited files that all the code reads in and is relative to. To add a new number of pcs or p-value threshold, just add a new line with the desired value. For formatting purposes please keep the numerically ascending order.
Then rerun the formatting scripts found in the txt_files directory