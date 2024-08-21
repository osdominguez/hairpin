
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=5) {
  stop("Five arguments must be supplied", call.=FALSE)
}

phen_path <- toString(args[1])
phen_name <- toString(args[2])
phen_id <- toString(args[3])
boot_n <- as.numeric(args[4])
pop <- toString(args[5])

for (as in c("as", "noas")) {
  out_dir <- paste0("/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/", pop, "/")
  tmp_dir <- paste0("/scratch/osdominguez/temp_boot/", pop, "/")
  boot_df <- data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("phenotype", "threshold", "pc_num", "r2", "theta_eo", "replicate")))) 

  for (i in 1:boot_n) {
    boot_i <- read.table(file = paste0(tmp_dir, phen_name, "_", as, "_bootstrap_", i, ".table"), header = TRUE)
    
    boot_df <- rbind(boot_df, boot_i)
  }

  write.table(boot_df,  file = paste0(out_dir, phen_name, "_", as, "_bootstrap.table"), row.names = F, quote = F) 

  confint_table <- boot_df %>% 
    group_by(phenotype, pc_num, threshold) %>% 
    summarise(r2_max = mean(r2, na.rm = TRUE) + (1.96 * sd(r2, na.rm = TRUE)), 
              r2_min = mean(r2,  na.rm = TRUE) - (1.96 * sd(r2, na.rm = TRUE)),
              thetaeo_max = mean(theta_eo, na.rm = TRUE) + (1.96 * sd(theta_eo, na.rm = TRUE)),
              thetaeo_min = mean(theta_eo, na.rm = TRUE) - (1.96 * sd(theta_eo, na.rm = TRUE)),
              num_reps = n())

  write.table(confint_table, file = paste0(out_dir, phen_name, "_", as, "_confint.table"), row.names = F, quote = F)
}