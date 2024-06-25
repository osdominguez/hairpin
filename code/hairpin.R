library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)<2) {
    stop("At least two argument must be supplied", call.=FALSE)
} else if (length(args)==2) {
    boot <- FALSE
} else if (length(args) == 3) {
    boot <- as.logical(args[3])
    boot_n <- 100
} else {
    boot <- as.logical(args[3])
    boot_n <- as.numeric(args[4])
}

phen_path <- args[1]
phen_name <- args[2]

# this function takes in the covariate dataframe (df), PRS dataframe (df), pheno_name (string), and p value (string)
calc_r2 <- function(covarpc_df, prs_df, pheno_df, pheno_name, pval) {
  
    # use the FID & IID columns to merge the dataframes
    full_df <- merge(covarpc_df, prs_df, by = c("FID", "IID")) 
    full_df <- merge(full_df, pheno_df, by =  c("FID", "IID"))

    # get rid of the FID and IID columns. This allows the use of the period in the following lm function. 
    # also get rid of the odd and the even PRS scores, we only use the "all" score here
    full_df <- full_df %>% select(-c(FID, IID))
    full_df <- full_df[, !grepl("even|odd", colnames(full_df))] 

    # run linear model for the full model and get its R^2
    full_mod <- lm(paste0(pheno_name," ~ ."), data = full_df, na.action = na.omit) 
    full_r2 <- summary(full_mod)$r.squared

    # run linear model for the null model and get its R^2
    null_mod <- lm(paste0(pheno_name," ~ . - all_",pval), data = full_df, na.action = na.omit) 
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
bootstrap_hairpin <- function(boot_num, pheno_filename, pheno_name) {

                                #### Setup for the Bootstrap ####
    # initialize final bootstrap table 
    bootstrap_hairpin_df <- data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("phenotype", "threshold", "pc_num", "r2", "theta_eo", "replicate")))) 

                                    #### Actually doing the bootsrap ####



    # for each dataframe, make the hairpin and append it to the final dataframe
    for (i in 1:boot_num){
    
        # call the hairpin function
        hairpin_output <- make_hairpin(pheno_filename, pheno_name, i) %>% mutate(replicate = i)
        
        # append
        bootstrap_hairpin_df <- rbind(bootstrap_hairpin_df, hairpin_output)
        
    } 

  write.table(bootstrap_hairpin_df,  file = paste0("/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting",pheno_name,"_bootstrap.table"), row.names = F, quote = F)

  # finding the confidence intervals for each threshold 
  confint_table <- bootstrap_hairpin_df %>% 
    group_by(phenotype, pc_num, threshold) %>% 
    summarise(r2_max = mean(r2, na.rm = TRUE) + (1.96 * sd(r2, na.rm = TRUE)), 
              r2_min = mean(r2,  na.rm = TRUE) - (1.96 * sd(r2, na.rm = TRUE)),
              thetaeo_max = mean(theta_eo, na.rm = TRUE) + (1.96 * sd(theta_eo, na.rm = TRUE)),
              thetaeo_min = mean(theta_eo, na.rm = TRUE) - (1.96 * sd(theta_eo, na.rm = TRUE)),
              num_reps = n()) 
  
  # return the final dataframe
  return(confint_table)
}

# this function generates the necessary column names for a specific number of pcs
pc_header <- function(n_pcs) {
    cols <- c("FID", "IID", "31-0.0", "34-0.0")
    for (i in 1:n_pcs) {
        cols <- append(cols, paste0("22009-0.", i))
    }
    return(cols)
}

# bootstrap function
make_hairpin <- function(pheno_path, pheno_name, boot_seed = 0) {

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
        pc_df <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC_age.pheno", header=TRUE, sep = " ") 
        pc_df <- pc_df %>% select(pc_header(pc))

        prs_df <- read.table(paste0("/scratch/osdominguez/tables", pheno_name, pc, ".table"), header=TRUE, sep = " ")

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
            r_2 <- calc_r2(pc_df, prs_subset, pheno_df, pheno_name, pval)
            
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

pcs <- scan(file.path(txt_path, '/pcs.txt'), what = integer())
pvals_list <- scan(file.path(txt_path, '/pvalues.txt'), what = character())

base <- make_hairpin(phen_path, phen_name)
write.table(base, file = paste0(out_dir, phen_name, "_base.table"), row.names = F, quote = F) 

if (boot) {
    boot_hairpin <- bootstrap_hairpin(boot_n, phen_path, phen_name)
    write.table(boot, file = paste0(out_dir, phen_name, "_confint.table"), row.names = F, quote = F) 
}