rm(list=ls())

library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=1) {
    stop("Only one argument may be supplied", call.=FALSE)
} else {
    phen_name <- toString(args[1])
}

txt_path <- file.path('/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files')
table_dir <- "/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/"

hairpin_df <- read.table(paste0(table_dir, phen_name, "_base.table"), header = TRUE)
boot_df <- read.table(paste0(table_dir, phen_name, "_bootstrap.table"), header = TRUE)
cof_df <- read.table(paste0(table_dir, phen_name, "_cofint.table"), header = TRUE)


get_resid <- function(full_boot, boot_n, phen, pc){
  # make an empty dataframe that just contains all pvalue thresholds
  uniq_p <- unique(full_boot$threshold)
  eps_df <- as.data.frame(uniq_p)
  colnames(eps_df) <- "threshold"
  
  # for every bootstrap replicate, find the resiuduals for each pvalue threshold
  for (i in 1:boot_n) {
    
    # filter down to the specific hairpin bootstrap replicate
    boot_rep <- full_boot %>% filter(replicate == i & pc_num == pc & phenotype == phen) 
    best_p <- subset(boot_rep, threshold <= 1e-13 & threshold >= 1e-18)
    
    # run the gls linear model for that hairpin
    lin_mod <- lm(formula=theta_eo~0 + r2, data = best_p)
    
    # create empty dataframe to store results for that bootstrap replicate
    rep_resid <- data.frame(matrix(ncol=2, nrow=0, dimnames=list(NULL, c("threshold", paste0("residuals_",i) ))))
    
    # finding residuals for each threshold
    for (pval in uniq_p) {    
      
      # find the y (theta_eo) and the x (r2) 
      y <- boot_rep$theta_eo[boot_rep$threshold == pval]
      x <- boot_rep$r2[boot_rep$threshold == pval]
      
      # find the lines estimation of y 
      y_hat <- predict(lin_mod, data.frame(r2 = x))[[1]]
      
      # find the residual, y - y_hat
      residual = y - y_hat
      
      # create a new row where the column name is the replicate id and another column indicates the threshold
      column_name <- paste0("residual_", i)
      new_row <- data.frame(threshold = pval)
      new_row[, column_name] <- residual
      
      # add the new row into the dataframe 
      rep_resid <- rbind(rep_resid, new_row)
    }
    # merge that replicates total residuals and the epsilon_df
    eps_df <- merge(eps_df, rep_resid, by = "threshold")
  } 
  # make episilon go from least to greatest by threshold 
  eps_df <- eps_df[order(eps_df$threshold), ]
  
  # return the epsilon dataframe
  return(eps_df)
}

