library(dplyr)

list_vector <- function(l) {
	res <- c()
	for (i in 1:length(l)) {
		res <- c(res, l[[i]])
	}
	return(res)
}

# this function generates the necessary column names for a specific number of pcs
# R doesn't like the "-" in some of the headers so it changes them to "."
pc_header <- function(n_pcs) {
  cols <- c("FID", "IID", "X31.0.0", "X21003.0.0", "X22000.0.0", "age2")
  if (as == "as") {
    cols <- c(cols, "X54.0.0")
  }
  if(n_pcs > 0) {
    for (i in 1:n_pcs) {
      cols <- append(cols, paste0("X22009.0.", i))
    }
  }
  return(cols)
}

get_ctypes <- function(df, header) {
	cols <- colnames(df)
  	cols_i <- which(cols %in% header)
  	ctypes <- lapply(df[1, ], class)
  	ctypes <- list_vector(ctypes)
  	ctypes[setdiff(1:ncol(df), cols_i)] <- "NULL"
	return(ctypes)
}

# This file is designed to make the hairpin more memory efficient in an array context

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Seven arguments must be supplied", call.=FALSE)
} 

phen_path <- toString(args[1])
phen_name <- toString(args[2])
phen_id <- toString(args[3])
boot <- as.logical(args[4])
run_n <- as.numeric(args[5])
as <- toString(args[6])
pop <- toString(args[7])

txt_path <- file.path('/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files')
out_dir <- paste0("/gpfs/data/ukb-share/dahl/ophelia/hairpin/plotting/", pop,"/")
tmp_dir <- paste0("/scratch/osdominguez/temp_boot/", pop,"/")

dir.create(out_dir, showWarnings = FALSE)
dir.create(tmp_dir, showWarnings = FALSE)

print('scanning in hairpin run txt files...')

pcs <- scan(file.path(txt_path, '/pcs.txt'), what = integer())
pvals_list <- scan(file.path(txt_path, '/pvalues.txt'), what = character())

print('finished!')

if (boot) {
  hfile <- paste0(tmp_dir, phen_name, "_", as, "_bootstrap_", run_n, ".table")
} else {
  hfile = paste0(out_dir, phen_name, "_", as, "_base.table")
}

if (file.exists(hfile)) {
  hairpin_df <- read.table(file = hfile, header = TRUE)
  pcs <- setdiff(pcs, as.numeric(unique(hairpin_df$pc_num)))
  rm(hairpin_df)
} else {
  if (boot) {
    hairpin_df <- data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("phenotype", "threshold", "pc_num", "r2", "theta_eo", "replicate"))))
  } else {
    hairpin_df <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("phenotype", "threshold", "pc_num", "r2", "theta_eo"))))  
  }
    write.table(hairpin_df, hfile, row.names = F, quote = F)
}

#read in phenotype file 

print('reading in phenotype file...')

pheno_df <- read.table(phen_path, header = TRUE)

print('successfully read in phenotype!')
#### Actually Creating the Hairpin ####

# for each number of principal components, read in the PC file then work on each p-value threshold 
for (pc in pcs) {
  
  hairpin_df <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("phenotype", "threshold", "pc_num", "r2", "theta_eo")))) 

  #read in PC table and the PRS table that contains all the scores for that PC  
  pc_df <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno", header=TRUE, sep = " ", nrows = 1) 
  
  ctypes_pc <- get_ctypes(pc_df, pc_header(pc))

  rm(pc_df)

  pc_df <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno", header=TRUE, sep = " ", colClasses = ctypes_pc)

  rm(ctypes_pc)

  prs_table <- paste0("/scratch/osdominguez/tables/", pop, "/", phen_name, pc, as, ".table")

  if (!file.exists(prs_table)) {
    next
  }
  
  for (pval in pvals_list) {
    
    print(paste('running hairpin for', pval, 'pvalue and', pc, 'pcs...'))
	tryCatch( {
    prs_df <- read.table(prs_table, header=TRUE, sep = " ", nrows = 1)
    ctypes_prs <- get_ctypes(prs_df, c("FID", "IID", paste0(c("even", "odd", "all"), paste0("_", pval))))

    prs_subset <- read.table(prs_table, header=TRUE, sep = " ", colClasses = ctypes_prs)
    
    # if a bootstrap is being called for resample the dataframe with a set seed 
    if (run_n != 0) {
      # set the seed and resample the dataframe
      set.seed(run_n)
      prs_subset <- prs_subset[sample(nrow(prs_subset), replace = TRUE), ]     
    }
    
    rm(ctypes_prs)
    rm(prs_df)

    # if we have the even, odd, and all columns go forward with getting r2 and thetaeo
    if (any(grepl("even", colnames(prs_subset))) && any(grepl("odd", colnames(prs_subset))) && any(grepl("all", colnames(prs_subset))))  {
            
      # use the FID & IID columns to merge the dataframes
      full_df <- merge(pc_df, prs_subset)
      full_df <- merge(full_df, pheno_df)
      full_df <- full_df %>%  select(-c("FID","IID")) 
      
      full_df <- full_df[, !grepl("even|odd", colnames(full_df))] 

      # Regress based on the ID of the phenotype (as shown up in the .pheno file)
      full_mod <- lm(paste0(phen_id," ~ ."), data = full_df, na.action = na.omit) 
      full_r2 <- summary(full_mod)$r.squared
      rm(full_mod)

      null_mod <- lm(paste0(phen_id," ~ . - all_",pval), data = full_df, na.action = na.omit) 
      null_r2 <- summary(null_mod)$r.squared
      rm(null_mod)
      
      r_2 <- (full_r2 - null_r2)/(1 - null_r2)
      rm(full_r2)
      rm(null_r2)
      rm(full_df)

      # make even and odd datafranes
      full_df <- merge(pc_df, prs_subset)
      full_df <- full_df %>%  select(-c("FID","IID"))
      even_df <- full_df[, !grepl("all|odd", colnames(full_df))]
      odd_df <- full_df[, !grepl("all|even", colnames(full_df))]
      rm(full_df)
      rm(prs_subset)

      # run linear models and get residuals for each function
      even_lm <- lm(paste0("even_", pval," ~ . "), data = even_df, na.action = na.omit)
      even_resid <- resid(even_lm)
      rm(even_lm)

      odd_lm <- lm(paste0("odd_",pval," ~ . "), data = odd_df, na.action = na.omit)
      odd_resid <- resid(odd_lm)
      rm(odd_lm)
      
      # calculate theta even-odd
      thetaeo <- cor(odd_resid, even_resid)
      
      rm(even_resid)
      rm(odd_resid)

      #create the new row to the hairpin dataframe
      new_row <- data.frame(phenotype = phen_name,
                            threshold = pval,
                            pc_num = pc, 
                            r2 = r_2,
                            theta_eo = thetaeo)
      rm(thetaeo)
      rm(r_2)
     

    } else {
      rm(prs_subset)
      # create the new row to the hairpin dataframe 
      new_row <- data.frame(phenotype = phen_name,
                            threshold = pval,
                            pc_num = pc, 
                            r2 = NA,
                            theta_eo = NA)
      
    }

    # add the new row to the dataframe
    hairpin_df <- rbind(hairpin_df, new_row) 
    print('success!')
  	}, error = function(msg) {
		print(paste('hairpin failed with error', msg))
		print("continuing")
	})
  }

  if (boot) {
    hairpin_df <- hairpin_df %>% mutate(replicate = run_n)
    write.table(hairpin_df,  file = paste0(tmp_dir, phen_name, "_", as, "_bootstrap_", run_n, ".table"), row.names = F, col.names = F, quote = F, append=TRUE)
  } else {
    write.table(hairpin_df, file = paste0(out_dir, phen_name, "_", as, "_base.table"), row.names = F, col.names = F, quote = F, append=TRUE)
  }
}

if (boot) {
  write(paste0(pop,"/",phen_name,"_",as,"_bootstrap_",run_n,".table"),file="/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files/fin_hairpin.txt",append=TRUE)
}
