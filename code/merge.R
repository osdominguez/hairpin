
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=4) {
  stop("Four arguments must be supplied", call.=FALSE)
}

phen_name <- toString(args[1])
as <- toString(args[2])
pop <- toString(args[3])
boot_n <- as.numeric(args[4])

out_dir <- paste0("/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/", pop, "/")
tmp_dir <- paste0("/scratch/osdominguez/temp_boot/", pop, "/")
boot_df <- data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("phenotype", "threshold", "pc_num", "r2", "theta_eo", "replicate")))) 

for (i in 1:boot_n) {
  boot_file <- paste0(tmp_dir, phen_name, "_", as, "_bootstrap_", i, ".table")
  if (!file.exists(boot_file)) {
	  next
  }
	
  boot_i <- read.table(file = boot_file, header = TRUE)
  
  boot_df <- rbind(boot_df, boot_i)
}

write.table(boot_df,  file = paste0(out_dir, phen_name, "_", as, "_bootstrap.table"), row.names = F, quote = F) 
