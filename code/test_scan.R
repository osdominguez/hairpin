library(dplyr) 
library(data.table)

txt_path <- file.path('/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files')

pcs <- scan(file.path(txt_path, '/pcs.txt'), what = integer())
pvals <- scan(file.path(txt_path, '/pvalues.txt'), what = character())



print(pcs)
print(pvals)