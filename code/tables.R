rm(list=ls())

library(dplyr) 
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=2) {
    stop("Only two arguments may be supplied", call.=FALSE)
} 

phen <- toString(args[1])
pop <- toString(args[2])

txt_path <- file.path('/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files')

pcs <- scan(file.path(txt_path, '/pcs.txt'), what = character())
pvals <- scan(file.path(txt_path, '/pvalues.txt'), what = character())
assess <- c("as", "noas")

# creates a function that has all the scores for every ID by PC number and stores the dataframes in a list
make_prs_df <- function(pheno) {

    for (as in assess) {
        # for every PC, add all the PRS score columns into one dataframe
        for (pc in pcs) {
            
            f <- paste0("/scratch/osdominguez/tables/", pop, "/", paste0(pheno,pc,as),".table")

            if (file.exists(f)) {
                # If table already exists read it in
                ids_df <- read.table(f, header=TRUE, sep =" ")
                alr_pvals <- colnames(ids_df)
            } else {
                # making new ID dataframe for each PC
                ids_df <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/white_british/whitebrit_unrelated.pheno", header=TRUE, sep = " ")
                alr_pvals <- c(0) 
                # for every pvalue add its PGS column  to the dataframe above
            }

            for (pval in pvals) {
                if (!(paste0("all_", pval) %in% alr_pvals)) {
                    odd_path <- paste0("/scratch/osdominguez/prs_hairpin_outputs/", pop, "/odd/", pheno,"_prs_",as,"_",pc,"pc_",pval,"pval.best")
                    even_path <- paste0("/scratch/osdominguez/prs_hairpin_outputs/", pop, "/even/", pheno,"_prs_",as,"_",pc,"pc_",pval,"pval.best")
                    all_path <- paste0("/scratch/osdominguez/prs_hairpin_outputs/", pop, "/all/", pheno,"_prs_",as,"_",pc,"pc_",pval,"pval.best")

                    o_exists <- try(read.table(odd_path, header=TRUE, sep = " ")) 
                    e_exists <- try(read.table(even_path, header=TRUE, sep= " ")) 
                    a_exists <- try(read.table(all_path, header=TRUE, sep = " ")) 

                    if (!is(o_exists, "try-error") & !is(e_exists, "try-error") & !is(a_exists, "try-error")) {

                        pval_odd <- read.table(odd_path, header=TRUE, sep = " ")
                        pval_even <- read.table(even_path, header=TRUE, sep = " ")
                        pval_all <- read.table(all_path, header = TRUE, sep = " ") 

                        names(pval_odd)[names(pval_odd) == 'PRS'] <- paste0("odd_",pval)
                        names(pval_even)[names(pval_even) == 'PRS'] <- paste0("even_",pval)     
                        names(pval_all)[names(pval_all) == 'PRS'] <- paste0("all_",pval)

                        ids_df <- merge(ids_df, pval_odd, by = c("FID","IID"), all = TRUE)
                        ids_df <- merge(ids_df, pval_even, by = c("FID","IID"), all = TRUE) 
                        ids_df <- merge(ids_df, pval_all, by = c("FID","IID"), all = TRUE) 

                        ids_df <- select(ids_df, -contains("In_Regression"))
                    }
                }    
            }

            write.table(ids_df, file = paste0("/scratch/osdominguez/tables/", pop, "/", paste0(pheno,pc,as),".table"), row.names = F, quote = F)
            
        }
    }
}


make_prs_df(phen)
