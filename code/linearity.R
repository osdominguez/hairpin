rm(list=ls())

library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Three arguments must be supplied", call.=FALSE)
} 

phen <- toString(args[1])
as <- toString(args[2])
pop <- toString(args[3])

table_dir <- "/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/"

blank_table <- data.frame(matrix(ncol=10, nrow=0, dimnames=list(NULL, c("phenotype", "pc_num", "threshold", "sd_dist", "se_dist", "mean_dist", "sd_teo", "z_jk", "pz", "on_line"))))

get_epsilon <- function(full_boot) {
  # Returns a dataframe with the residuals of each bootstrap run at each threshold for a phenotype and pc of interest.
  
  uniq_p <- unique(full_boot$threshold)
  eps_df <- as.data.frame(uniq_p)
  colnames(eps_df) <- "threshold"
  
  # for every bootstrap replicate, find the resiuduals for each pvalue threshold
  for (i in unique(full_boot$replicate)) {
    
      # filter down to the specific hairpin bootstrap replicate
      boot_rep <- full_boot %>% filter(replicate == i)
      
      # create empty dataframe to store results for that bootstrap replicate
      rep_resid <- data.frame(matrix(ncol=2, nrow=0, dimnames=list(NULL, c("threshold", paste0("residuals_",i) ))))
      
      # run the ols linear model for that hairpin
      lin_mod <- lm(formula=theta_eo~0 + r2, data = boot_rep)
      
      # finding residuals for each threshold
      for (pval in uniq_p) {    
        
        # find the y (theta_eo) and the x (r2) 
        y <- boot_rep$theta_eo[boot_rep$threshold == pval]
        x <- boot_rep$r2[boot_rep$threshold == pval]
        
        # find the line estimation of y 
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
  
  # order epsilon by threshold
  eps_df <- eps_df[order(eps_df$threshold), ]
  
  eps_df <- na.omit(eps_df)

  return(eps_df)
}

get_omega <- function(eps) {
  
  # Obtain omega from the given epsilons
  
  eps <- t(eps)
  colnames(eps) <- eps[1, ]
  eps <- eps[-1, ]
  
  eps <- as.matrix(eps)
  
  omega <- cov(eps)
  
  # To make sure the matrix has nrow(omega) linearly independent bases
  omega <- omega + diag(1e-12, nrow(omega))
  
  return(omega)
}

get_beta <- function(omega, hairpin) {
  
  X <- as.vector(hairpin$r2)
  Y <- as.vector(hairpin$theta_eo)
  
  beta <- solve(t(X)%*%solve(omega)%*%X)%*%t(X)%*%solve(omega)%*%Y
  
  return(as.numeric(beta))
  
}

get_dist <- function(hair, pmax, omega_j, rep = 0, tail) {
  
  # return the distances of the points on the gls line for a given hairpin
    
  distances_df <-  data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("threshold", "distance", "a_o", "replicate", "theta_eo"))))
  
  hair <- hair[order(hair$threshold), ]
  
  hair_j <- hair %>% filter(threshold <= pmax)
  
  beta_j <- get_beta(omega_j, hair_j)
  
  rm(hair_j)
  
  line_j <- beta_j * hair$r2 
  
  for (k in 1:length(hair$threshold)) {
    
    thresh <- hair$threshold[k]
    teo <- hair$theta_eo[k]
    
    if (tail == "two") {
      dist_jk <- line_j[k] - teo
    } else {
      dist_jk <- abs(line_j[k] - teo)
    }
    
    new_row <- data.frame(threshold = thresh,
                          dist = dist_jk, 
                          a_o = case_when(
                            (dist_jk < 0) ~ "below",
                            (dist_jk > 0) ~ "above",
                            .default = "on"
                          ),
                          replicate = rep,
                          theta_eo = teo) 
    
    distances_df <- rbind(distances_df, new_row)
    
  }
  
  return(distances_df)
}

get_zvals <- function(hair, boot, method, tail, pmax = NULL, prev = NULL) {
  
  # returns the ztable for the given hairpin, assuming the hair and boot dfs
  # have been filtered for the specific phenotype and pc of interest
  
  if (!(method %in% c("forward", "backward"))) {
    stop('ERROR: please enter a valid mathod')
  } 
  
  distances_df <-  data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("threshold", "distance", "a_o", "replicate", "theta_eo"))))
  
  reps <- unique(boot$replicate)
  n_boots <- length(reps)
  
  threshs <- unique(hair$threshold)
  threshs <- threshs[order(threshs)]
  
  if (is.null(pmax)) {
    pmax = case_when(
      method == "backward" ~ max(hair$threshold),
      method == "forward" ~ min(hair$threshold)
    )
  }
  
  boot_j <- boot %>% filter(threshold <= pmax)
  hair_j <- hair %>% filter(threshold <= pmax)
  
  eps_j <- get_epsilon(boot_j)
  
  omega_j <- get_omega(eps_j)
  
  hairbeta_j <- get_beta(omega_j, hair_j)
  
  for (i in reps) {
    
    boot_i <- boot %>% filter(replicate == i)

    dists_i <- get_dist(boot_i, pmax, omega_j, i, tail)
    
    distances_df <- rbind(distances_df, dists_i)
  }
  
  dists_hair <- get_dist(hair, pmax, omega_j, 0, tail) %>% select(-replicate)
  
  dists_hair$beta <- rep(hairbeta_j, nrow(dists_hair))
  
  z_table <- distances_df %>% 
    group_by(threshold) %>% 
    summarise(mean_dist = mean(dist, na.rm = TRUE),
              sd_dist = sqrt(sum((dist - mean_dist)^2)/(n_boots - 1)),
              se_dist = sqrt(sd_dist),
              sd_teo = sd(theta_eo, na.rm = TRUE))
  
  rm(distances_df)  
  
  z_table <- merge(z_table, dists_hair)
  
  rm(dists_hair)
  
  z_table <- z_table %>% mutate(z_jk = (mean_dist - 0)/sd_dist)
  
  if (tail == "two") {
    z_table <- z_table %>% mutate(pz = pnorm(z_jk, mean = 0, sd = 1, lower.tail = TRUE))
  } else {
    z_table <- z_table %>% mutate(pz = pnorm(z_jk, mean = 0, sd = 1, lower.tail = FALSE))
  }
  
  z_table <- z_table %>% mutate(on_line = case_when(pz < 0.05 ~ "off",
                                                    pz > 0.05 ~ "on"))
  
  offs <- (z_table %>% filter(on_line == "off"))$threshold
  
  if (method == "backward") {
    if (any(offs <= pmax) || hairbeta_j < 0) {
      if (pmax == min(threshs)) {
        return(z_table)
      } else {
        pmax_n <- threshs[which(threshs == pmax) - 1]
        return(get_zvals(hair, boot, method = method, tail = tail, pmax = pmax_n))
      } 
    } else {
      return(z_table)
    }
  } else {
    if (any(offs <= pmax)) {
      return(prev) 
    } else {  
      if (pmax == max(threshs)) {
        return(z_table)
      } else {
        pmax_n <- threshs[which(threshs == pmax) + 1]
        return(get_zvals(hair, boot, method = method, tail = tail, pmax = pmax_n, ztable))
      }
    }
  }
}

ztable <- function(phen_name, as, pop, method = "backward", tail = "two") {
  
  # returns the full dataframe of the given hairpin tested at gls
  
  hair_file <- paste0(table_dir, pop, "/", phen_name, "_", as, "_base.table")
  boot_file <- paste0(table_dir, pop, "/", phen_name, "_", as, "_bootstrap.table")

  if (!file.exists(hair_file) || !file.exists(boot_file)) {
    stop('ERROR: either the bootstrap or hairpin table does not exist')
  }

  hairpin_df <- read.table(hair_file, header = TRUE)
  boot_df <- read.table(boot_file, header = TRUE)

  hairpin_df <- na.omit(hairpin_df) %>% filter(phenotype == phen_name)
  boot_df <- na.omit(boot_df) %>% filter(phenotype == phen_name)
  
  f_df <- blank_table 
  
  print(paste0("Testing hairpin linearity for ", phen_name, " ", as, " in ", pop, "..."))
  
  for (pc in unique(unique(hairpin_df$pc_num))) {
    
    t_df <- get_zvals(hairpin_df %>% filter(pc_num == pc), 
                      boot_df %>% filter(pc_num == pc),
                      method,
                      tail)
    
    t_df$phenotype = rep(phen_name, nrow(t_df))
    t_df$pc_num = rep(pc, nrow(t_df))
    
    f_df <- rbind(f_df, t_df)
    
    print(paste0("finished testing pc ", pc))
  }
  
  print("Finished testing linearity!")
  
  f_df <- f_df[order(f_df$pc_num, f_df$threshold), ]
  f_df <- merge(f_df, hairpin_df)
  
  return(f_df)
  
}

ztab_phen <- ztable(phen, as, pop)

if (nrow(ztab_phen) == 0) {
	exit("ERROR: empty ztable")
} else {
	ztab_file <- paste0(table_dir, pop, "/", "ztab_", phen, "_", as, ".table")
	
	print("writing ztable...")

	write.table(ztab_phen, file = ztab_file, row.names = F, quote = F)
	
	print("successfully wrote ztable!")
}

