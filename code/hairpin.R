library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Six arguments must be supplied", call.=FALSE)
} 

phen_path <- toString(args[1])
phen_name <- toString(args[2])
phen_id <- toString(args[3])
boot <- as.logical(args[4])
run_n <- as.numeric(args[5])
as <- toString(args[6])

# this function takes in the covariate dataframe (df), PRS dataframe (df), pheno_name (string), and p value (string)
calc_r2 <- function(covarpc_df, prs_df, pheno_df, pheno_id, pval) {
  
  # use the FID & IID columns to merge the dataframes
  full_df <- merge(covarpc_df, prs_df, by = c("FID", "IID")) 
  full_df <- merge(full_df, pheno_df, by =  c("FID", "IID"))
  
  # get rid of the FID and IID columns. This allows the use of the period in the following lm function. 
  # also get rid of the odd and the even PRS scores, we only use the "all" score here
  full_df <- full_df %>% select(-c(FID, IID))
  full_df <- full_df[, !grepl("even|odd", colnames(full_df))] 
  
  # Regress based on the ID of the phenotype (as shown up in the .pheno file)
  full_mod <- lm(paste0(pheno_id," ~ ."), data = full_df, na.action = na.omit) 
  full_r2 <- summary(full_mod)$r.squared
  
  null_mod <- lm(paste0(pheno_id," ~ . - all_",pval), data = full_df, na.action = na.omit) 
  null_r2 <- summary(null_mod)$r.squared
  
  r2 <- full_r2 - null_r2
  
  return(r2)
}

calc_thetaeo <- function(covarpc_df, prs_df, pval) {
  
  # use the FID & IID columns to merge the dataframes
  full_df <- merge(covarpc_df, prs_df, by = c("FID", "IID"))
  full_df <- full_df %>%  select(-c("FID","IID")) 
  
  # make even and odd datafranes
  even_df <- full_df[, !grepl("all|odd", colnames(full_df))]
  odd_df <- full_df[, !grepl("all|even", colnames(full_df))]
  
  # run linear models and get residuals for each function
  even_lm <- lm(paste0("even_", pval," ~ . "), data = even_df, na.action = na.omit)
  even_resid <- resid(even_lm)
  
  odd_lm <- lm(paste0("odd_",pval," ~ . "), data = odd_df, na.action = na.omit)
  odd_resid <- resid(odd_lm)
  
  # calculate theta even-odd
  theta_eo_num <- cor(odd_resid, even_resid)
  
  return(theta_eo_num)
  
}

# this function bootstraps the hairpin for a specific phenotype
bootstrap_hairpin <- function(boot_num, pheno_filename, pheno_name, pheno_id, as) {
  
  hairpin_output <- make_hairpin(pheno_filename, pheno_name, pheno_id, as, boot_num) %>% mutate(replicate = boot_num)
  
  write.table(hairpin_output,  file = paste0(tmp_dir, pheno_name, "_", as, "_bootstrap_", boot_num, ".table"), row.names = F, quote = F)
}

# this function generates the necessary column names for a specific number of pcs
# R doesn't like the "-" in some of the headers so it changes them to "."
pc_header <- function(n_pcs) {
  cols <- c("FID", "IID", "X31.0.0", "X34.0.0", "X54.0.0", "age2")
  if(n_pcs <= 0) {
    return(cols)
  }
  for (i in 1:n_pcs) {
    cols <- append(cols, paste0("X22009.0.", i))
  }
  return(cols)
}

# bootstrap function
make_hairpin <- function(pheno_path, pheno_name, pheno_id, as, boot_seed = 0) {
  
  #### Setup for the Hairpin ####
  
  #initialize empty hairpin dataframe
  
  ## give bootsrap argument 
  
  hairpin_df <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("phenotype", "threshold", "pc_num", "r2", "theta_eo"))))  
  
  #read in phenotype file 
  pheno_df <- read.table(pheno_path, header = TRUE)
  
  #### Actually Creating the Hairpin ####
  
  # for each number of principal components, read in the PC file then work on each p-value threshold 
  for (pc in pcs) {
    
    #read in PC table and the PRS table that contains all the scores for that PC  
    pc_df <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno", header=TRUE, sep = " ") 
    
    pc_df <- pc_df %>% select(pc_header(pc))
    
    prs_df <- read.table(paste0("/scratch/osdominguez/tables/", pheno_name, pc, as, ".table"), header=TRUE, sep = " ")
    
    # if a bootstrap is being called for resample the dataframe with a set seed 
    if (boot_seed != 0) {
      # set the seed and resample the dataframe
      set.seed(boot_seed)
      prs_df <- prs_df[sample(nrow(prs_df), replace = TRUE), ]     
    }
    
    for (pval in pvals_list) {
      
      # select the FID and IID of participants. Also select the column names with the pvalue of interest in it
      prs_subset <- prs_df %>% select(c(FID, IID, contains(pval)))       
      
      # if we have the even, odd, and all columns go forward with getting r2 and thetaeo
      if (any(grepl("even", colnames(prs_subset))) && any(grepl("odd", colnames(prs_subset))) && any(grepl("all", colnames(prs_subset))))  {
        
        # call the hairpin functions
        thetaeo <- calc_thetaeo(pc_df,prs_subset,pval)
        r_2 <- calc_r2(pc_df, prs_subset, pheno_df, phen_id, pval)
        
        #create the new row to the hairpin dataframe
        new_row <- data.frame(phenotype = pheno_name,
                              threshold = pval,
                              pc_num = pc, 
                              r2 = r_2,
                              theta_eo = thetaeo)
      } else {
        # create the new row to the hairpin dataframe 
        new_row <- data.frame(phenotype = pheno_name,
                              threshold = pval,
                              pc_num = pc, 
                              r2 = NA,
                              theta_eo = NA)
        
      }
      
      # add the new row to the dataframe
      hairpin_df <- rbind(hairpin_df, new_row)
      
    }
  }
  
  # return the completed hairpin dataframe
  return(hairpin_df)
  
}

txt_path <- file.path('/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files')
out_dir <- "/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/"
tmp_dir <- "/scratch/osdominguez/temp_boot/"

pcs <- scan(file.path(txt_path, '/pcs.txt'), what = integer())
pvals_list <- scan(file.path(txt_path, '/pvalues.txt'), what = character())

if (boot) {
  bootstrap_hairpin(run_n, phen_path, phen_name, phen_id, as)
} else {
  base <- make_hairpin(phen_path, phen_name, phen_id, as)
  write.table(base, file = paste0(out_dir, phen_name, "_", as, "_base.table"), row.names = F, quote = F)
}
