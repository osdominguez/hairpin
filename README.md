# Most of the code here was made by Shevaugn originally or very heavily inspired by Renée

# Covars path:

## Without smoking

/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC_age.pheno
/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC_sex.pheno

## With Smoking

/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC_smokingA.pheno

# Covar Fields (for PRSice) :

## Without smoking:
FID IID 31-0.0 34-0.0 22009-0.1 22009-0.2 22009-0.3 22009-0.4 22009-0.5 22009-0.6 22009-0.7 22009-0.8 22009-0.9 22009-0.10 22009-0.11 22009-0.12 22009-0.13 22009-0.14 22009-0.15 22009-0.16 22009-0.17 22009-0.18 22009-0.19 22009-0.20 22009-0.21 22009-0.22 22009-0.23 22009-0.24 22009-0.25 22009-0.26 22009-0.27 22009-0.28 22009-0.29 22009-0.30 22009-0.31 22009-0.32 22009-0.33 22009-0.34 22009-0.35 22009-0.36 22009-0.37 22009-0.38 22009-0.39 22009-0.40

## With smoking: 
FID IID 31-0.0 34-0.0 22009-0.1 22009-0.2 22009-0.3 22009-0.4 22009-0.5 22009-0.6 22009-0.7 22009-0.8 22009-0.9 22009-0.10 22009-0.11 22009-0.12 22009-0.13 22009-0.14 22009-0.15 22009-0.16 22009-0.17 22009-0.18 22009-0.19 22009-0.20 22009-0.21 22009-0.22 22009-0.23 22009-0.24 22009-0.25 22009-0.26 22009-0.27 22009-0.28 22009-0.29 22009-0.30 22009-0.31 22009-0.32 22009-0.33 22009-0.34 22009-0.35 22009-0.36 22009-0.37 22009-0.38 22009-0.39 22009-0.40 20116

# Pheno paths:

## Height:
/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno

# PGS Specs:

In hairpin/txt_files there are txt files for both the pcs and pvalues. These are new line delimited files that all the code reads in and is relative to. To add a new number of pcs or p-value threshold, just add a new line with the desired value. For formatting purposes please keep the numerically ascending order.