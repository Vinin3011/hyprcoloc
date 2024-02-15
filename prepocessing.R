library(hyprcoloc)
library(dplyr)
library(data.table)
library(purrr)

# create array of paths
path_array <- array(list())

ms_paths <- c("Harmonized GWAS results/single_GWAS/pmid21833088_MS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid24076602_MS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid27386562_MS_eur.tsv.gz")
path_array[[1]] <- ms_paths

als_paths <- c("Harmonized GWAS results/single_GWAS/pmid27455348_ALS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid28931804_ALS_eur-eas.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid29566793_ALS.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid29566793_ALS_eur.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid34873335_ALS_eur.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid34873335_ALS_eur-eas.tsv.gz")
path_array[[2]] <- als_paths

t2d_paths <-c("Harmonized GWAS results/single_GWAS/pmid30054458_T2D_eur-sas.tsv.gz")
path_array[[3]] <- t2d_paths

pd_paths <-c("Harmonized GWAS results/single_GWAS/pmid33111402_PD_eur.tsv.gz",
             "Harmonized GWAS results/single_GWAS/pmid33987465_PD_eur.tsv.gz",
             "Harmonized GWAS results/single_GWAS/pmid34594039_PD_meta-eas-eur.tsv.gz",
             "Harmonized GWAS results/single_GWAS/pmid34594039_PD_nonmeta-eas.tsv.gz")
path_array[[4]] <- pd_paths

ra_paths <-c("Harmonized GWAS results/single_GWAS/pmid33310728_RA_eas-eur.tsv.gz")
path_array[[5]] <- ra_paths

dlb_paths <-c("Harmonized GWAS results/pmid33589841_DLB_eur.tsv.gz")
path_array[[6]] <- dlb_paths

scz_paths <-c("Harmonized GWAS results/single_GWAS/pmid34594039_SCZ_meta-eas-eur.tsv.gz",
             "Harmonized GWAS results/single_GWAS/pmid34594039_SCZ_nonmeta-eas.tsv.gz")
path_array[[7]] <- scz_paths

c3_paths <-c("Harmonized GWAS results/single_GWAS/pmid34648354_C3_eur.tsv.gz")
path_array[[8]] <- c3_paths

crp_paths <-c("Harmonized GWAS results/single_GWAS/pmid34648354_CRP_eur.tsv.gz")
path_array[[9]] <- crp_paths

aph_paths <-c("Harmonized GWAS results/single_GWAS/pmid34737426_APH_eur.tsv.gz")
path_array[[10]] <- aph_paths

# print path array
path_array

# Create lists to fill with betas and ses df
betas_df_list <- list()
ses_df_list <- list()


# Fill list
for (i in seq_along(path_array)) {
  sublist <- path_array[[i]]  # Get the sublist
  
  # Get trait name
  trait <- extract_trait(sublist[1])
  print(paste("Creating betas and ses dataframe for: ", trait))
  
  # attempt to create betas df from sublist and continue if an error occurs
  result_betas_df <- tryCatch(create_betas_df(sublist), error = function(e){
    print(paste("create_betas_df produced the following error for trait ", trait, ":",e))
    print("continue with next trait...")
    return(NULL)  # Return NULL to assign an empty dataframe
  })
  
  # Skip the rest of the loop if an error occurred
  if(is.null(result_betas_df)){
    next  
  }
  
  # attempt to create betas df from sublist and continue if an error occurs
  result_ses_df <- tryCatch(create_ses_df(sublist), error = function(e){
    print(paste("create_betas_df produced the following error for trait ", trait, ":",e))
    print("continue with next trait...")
    return(NULL)  # Return NULL to assign an empty dataframe
  })
  
  # Skip the rest of the loop if an error occurred
  if(is.null(result_ses_df)){
    next  
  }
  
  # Store the result dataframe along with the trait name in the results list
  betas_df_list[[trait]] <- result_betas_df
  ses_df_list[[trait]] <- result_ses_df
}

# create subsets
betas_test_list <- betas_df_list[1:4]
ses_test_list <- ses_df_list[1:4]

# Merge the betas data frames based on the "SNP" column
merged_betas <- reduce(betas_test_list, inner_join, by = "SNP")

# Merge the ses data frames based on the "SNP" column
merged_ses <- reduce(ses_test_list, inner_join, by = "SNP")

# remove zero values to avoid error in hyprcoloc
merged_ses <- merged_ses[merged_ses$MS != 0, ]

# merged_ses_first1000 <- head(merged_ses_nonzero, 1000)

# Convert into matrices ------------------------------------

# betas 

# Extract row names from the first column
rownames <- merged_betas[, 1]
# Remove the first column before converting to matrix
merged_betas <- merged_betas[, -1]
# Convert dataframe to matrix and set row names
merged_betas_matrix <- as.matrix(merged_betas)
rownames(merged_betas_matrix) <- rownames


# ses 

# Extract row names from the first column
rownames <- merged_ses[, 1]
# Remove the first column before converting to matrix
merged_ses <- merged_ses[, -1]
# Convert dataframe to matrix and set row names
merged_ses_matrix <- as.matrix(merged_ses)
rownames(merged_ses_matrix) <- rownames

