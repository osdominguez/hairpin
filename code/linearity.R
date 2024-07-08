rm(list=ls())

library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=4) {
    stop("Only four arguments may be supplied", call.=FALSE)
} else {
    phen_name <- toString(args[1])
    pc_n <- as.numeric(args[2])
    as <- toString(args[3])
    boot_num <- as.numeric(args[4])
}

txt_path <- file.path('/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files')
table_dir <- "/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/"

hairpin_df <- read.table(paste0(table_dir, phen_name, "_", as, "_base.table"), header = TRUE)
boot_df <- read.table(paste0(table_dir, phen_name, "_", as, "_bootstrap.table"), header = TRUE)
cof_df <- read.table(paste0(table_dir, phen_name, "_", as, "_cofint.table"), header = TRUE)

hairpin_df <- hairpin_df %>% filter(pc_num == pc_n & phenotype == phen_name)
boot_df <- boot_df %>% filter(pc_num == pc_n & phenotype == phen_name)

hairpin_df <- hairpin_df[order(hairpin_df$threshold), ]

get_epsilon <- function(full_boot, boot_n, phen){
  # make an empty dataframe that just contains all pvalue thresholds
  uniq_p <- unique(full_boot$threshold)
  eps_df <- as.data.frame(uniq_p)
  colnames(eps_df) <- "threshold"
  
  # for every bootstrap replicate, find the resiuduals for each pvalue threshold
  for (i in 1:boot_n) {
    
    # filter down to the specific hairpin bootstrap replicate
    boot_rep <- full_boot %>% filter(replicate == i)
    # ask about why we're doing this here...
    # is it the GLS assumption? Or like we really don't want to use the small values because the whole purpose of GLS is to figure out which ones (of the small values) are bad? 
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

  # make column names the thresholds

  # return the epsilon dataframe
  return(eps_df)
}

get_beta <- function(eps, hairpin) {
  
  eps <- t(eps)
  colnames(eps) <- eps[1, ]
  eps <- eps[-1, ]

  
  omega <- cov(eps)

  # To make sure the matrix has nrow(omega) linearly independent bases
  omega <- omega + diag(1e-12, nrow(omega))

  thresholds_to_keep <- as.list(hairpin$threshold) 
  omega <- omega[rownames(omega) %in% thresholds_to_keep, colnames(omega) %in% thresholds_to_keep]

  X <- as.vector(hairpin$r2)
  Y <- as.vector(hairpin$theta_eo)

  beta <- solve(t(X)%*%solve(omega)%*%X)%*%t(X)%*%solve(omega)%*%Y

}

# function that finds z score and pvalue for each threshold in the dataframe
get_zvals <- function(full_bootstrap_ols,full_bootstrap_gls, pheno, epsilon_df){
  
  ##### Setting up the function ####
  
  # create empty datafraem containing distances to append on to
  distances_df <-  data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("threshold", "distance", "a_o", "pc_num", "replicate", "phenotype"))))
  
  # create list of pvalue thresholds present for current dataframe to iterate through 
  thresholds <- as.list(unique(full_bootstrap_ols$threshold))   
  
  ##### getting distances for every threshold for every replicate ####
  
  # for each hairpin (dictated by number of bootstraps present), find the distance for each threshold
  for (i in 1:boot_num) {
    
    # filter down to the specific hairpin bootstrap replicate
    boot_replicate_gls <- full_bootstrap_gls %>% filter(replicate == i) 
    
    beta <- get_beta(epsilon_df, boot_replicate_gls)
  
    # for each threshold find the distance 
    for (pval in thresholds) {
      
      # the distance between a point in a line is found using the formula | ax_0 + by_0 + c | / sqrt(a^2 + b^2)
      # where point P is P=(x_0, y_0) and line L is ax + by + c = 0. 
      
      # to find point P we need an x and y. The x is r2 and y is the theta_eo associated with the threshold. 
      # Found using code below. 
      x <- boot_replicate_gls$r2[boot_replicate_gls$threshold == pval]
      y <- boot_replicate_gls$theta_eo[boot_replicate_gls$threshold == pval]
      
      
      # linear model essentially gives y = mx + b where b is the intercept, m is beta, and y is prediction.
      # to get y = mx + b into standard form for distance formula just substract y from both sides. Thus, b 
      # is -1 (cause its multiplied by y), c is b (the intercept), and a is m (the beta). 
      a = beta
      b = -1
      c = 0 
      
      # plug it into distance formula
      distance = abs(a * x + b * y + c ) / sqrt(a^2 + b^2)
      
      # find the y-value for the line
      line_x <- boot_replicate_gls$r2[boot_replicate_gls$threshold == pval]
      line_y <- 0 + beta * line_x
      
      if ((!is.na(y) & !is.na(line_y)) & (y <= line_y)){
        distance <- -distance
      } 
      
      # create a new row to add the data into
      new_row <- data.frame(threshold = pval,
                            distance = distance, 
                            a_o = case_when(
                              !is.na(y) & !is.na(line_y) & (y <= line_y) ~ "below",
                              !is.na(y) & !is.na(line_y) & (y > line_y) ~ "above"
                            ),
                            pc_num = pc,
                            replicate = i, 
                            phenotype = pheno) 
      
      # add it into the dataframe 
      distances_df <- rbind(distances_df, new_row)
    }
  }
  
  ##### getting the z-score ######
  
  # group by threshold and calculate the z score
  z_table <- distances_df %>% 
    group_by(threshold) %>% 
    summarise(z_score_s = (mean(distance, na.rm = TRUE) - 0 ) / sd(distance, na.rm = TRUE) ,
              z_score_sn =  (mean(distance, na.rm = TRUE) - 0 ) / ( sd(distance, na.rm = TRUE) / sqrt(boot_num)))
  
  # use the z score to calculate the pvalue then make a categorical column
  z_table <- z_table %>% 
    mutate(p_value = pnorm(q=z_score_s,mean = 0, sd = 1, lower.tail=FALSE),
           on_line = case_when(
             p_value <= 0.05 ~ "off",
             p_value > 0.05 ~ "on"
           ), 
           pc_num = pc, 
           phenotype = pheno) 
  
  df_list <- list(distances_df, z_table)
  
  return(df_list)
  
}

epsilon <- get_residuals(full_boot, boot_num, pheno)

beta <- get_beta(epsilon, hairpin_df)

test <- get_zvals(boot_df, boot_df, pheno, pc)

ztable <- test[[2]] 

full_df <- merge(hairpin_df, ztable, by = c("threshold", "pc_num"))

write.table(ztable, file = paste0("/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/linearity/hairpin_", pheno, pc, "_", as, "_", "gls.png"), row.names = F, quote = F)

plot <- ggplot(full_df, aes(x=r2, y=theta_eo, colour = on_line, label = threshold)) + 
  geom_line(colour = "black") +
  geom_point() +
  geom_text(hjust = 1, vjust = 0) +
  theme_bw()  +
  labs(title = paste0("Hairpin ", pheno, " plot ", pc, " PC (gls)"), y = "theta even/odd", x = "R2")  +
  scale_color_manual(values = c("on" = "black", "off" = "red")) + 
  geom_abline(intercept = 0, slope = beta, color = "blue", linewidth = .5)


ggsave(filename = paste0("/gpfs/data/ukb-share/dahl/ophelia/hairpin/plots/hairpin_", pheno, pc, "_", as, "_", "gls.png"), plot = plot, width = 6, height = 4)