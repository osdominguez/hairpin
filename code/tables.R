library(dplyr) 
library(data.table)

txt_path <- file.path('/gpfs/data/ukb-share/dahl/ophelia/hairpin/txt_files')

pcs <- scan(file.path(txt_path, '/pcs.txt'), what = character())
pvals <- scan(file.path(txt_path, '/pvalues.txt'), what = character())

# creates a function that has all the scores for every ID by PC number and stores the dataframes in a list
make_prs_df <- function(pheno){

    # for every PC, add all the PRS score columns into one dataframe
    for (pc in pcs){

        # making new ID dataframe for each PC
        ids_df <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/white_british/whitebrit_unrelated.pheno", header=TRUE, sep = " ") 
        # for every pvalue add its PGS column  to the dataframe above
        
        for (pval in pvals) {

            odd_path <- paste0("/scratch/osdominguez/prs_hairpin_outputs/odd/", pheno,"_prs_",pc,"pc_",pval,"pval.best")
            even_path <- paste0("/scratch/osdominguez/prs_hairpin_outputs/even/", pheno,"_prs_",pc,"pc_",pval,"pval.best")
            all_path <- paste0("/scratch/osdominguez/prs_hairpin_outputs/all/", pheno,"_prs_",pc,"pc_",pval,"pval.best")

            o_exists <- try(read.table(odd_path, header=TRUE, sep = " ")) 
            e_exists <- try(read.table(even_path, header=TRUE, sep= " ")) 
            a_exists <- try(read.table(all_path, header=TRUE, sep = " ")) 

            if (!is(o_exists, "try-error") & !is(e_exists, "try-error") & !is(a_exists, "try-error")){

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

        write.table(ids_df, file = paste0("/scratch/osdominguez/tables/", paste0(pheno,pc),".table"), row.names = F, quote = F)
        
    }
}


make_prs_df("height")
make_prs_df("BMI")
make_prs_df("LDL") 
#make_prs_df("edu_years", prs_dfs)
