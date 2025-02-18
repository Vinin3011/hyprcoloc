library(hyprcoloc)
library(dplyr)
library(data.table)
library(purrr)


traits_of_interest <- list("ALS_80610_29566793-GCST005647-EFO_0000253")

# cluster traits by region?
region_traits <- FALSE
complete_names <- TRUE

# paths for cytokines and single GWAS
cytokines_directory <- "Harmonized GWAS results/cytokines_new"
single_GWAS_directory <- "Harmonized GWAS results/full_unfiltered"
metabolite_directory <- "Harmonized GWAS results/metabolites"
cytokines_paths <- list.files(cytokines_directory, full.names = TRUE)
single_GWAS_paths <- list.files(single_GWAS_directory, full.names = TRUE)
metabolite_paths <- list.files(metabolite_directory, full.names = TRUE)

# group the paths by traits in lists
single_GWAS_traits_list <- group_paths_by_trait(single_GWAS_paths, complete_names = TRUE)
cytokine_traits_list <- group_paths_by_trait(cytokines_paths, fullnames = region_traits)
metabolite_traits_list <- group_paths_by_trait(metabolite_paths, fullnames = region_traits)

all_traits_list <- c(single_GWAS_traits_list, cytokine_traits_list, metabolite_traits_list)

# get paths of interest
paths_of_interest <- get_paths_of_interest(all_traits_list, traits_of_interest)

# Create lists to fill with betas and ses df
betas_df_list <- list()
ses_df_list <- list()


# Create ses and betas df for every trait of interest
for (i in seq_along(paths_of_interest)) {
  sublist <- paths_of_interest[[i]]  # Get the sublist
  
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

# Merge the betas data frames based on the "SNP" column
merged_betas <- reduce(betas_df_list, inner_join, by = "SNP")

# Merge the ses data frames based on the "SNP" column
merged_ses <- reduce(ses_df_list, inner_join, by = "SNP")

# Get the column names except the first one (assuming it's "SNP")
cols <- names(merged_ses)[-1]

# Remove rows with zero values in any column except the first one
merged_ses_nonzero <- merged_ses[rowSums(merged_ses[cols] == 0) == 0, ]
merged_betas_nonzero <- merged_betas[rowSums(merged_betas[cols] == 0) == 0, ]


# Convert into matrices ------------------------------------

# betas 

# Extract row names from the first column
rownames <- merged_betas_nonzero[, 1]
# Remove the first column before converting to matrix
merged_betas_nonzero <- merged_betas_nonzero[, -1]
# Convert dataframe to matrix and set row names
merged_betas_matrix <- as.matrix(merged_betas_nonzero)
rownames(merged_betas_matrix) <- rownames


# ses 

# Extract row names from the first column
rownames <- merged_ses_nonzero[, 1]
# Remove the first column before converting to matrix
merged_ses_nonzero <- merged_ses_nonzero[, -1]
# Convert dataframe to matrix and set row names
merged_ses_matrix <- as.matrix(merged_ses_nonzero)
rownames(merged_ses_matrix) <- rownames

