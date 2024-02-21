library(hyprcoloc)
library(dplyr)
library(data.table)
library(purrr)

# All trait groups of interest
trait_groups_of_interest <- list(
  list("MCP1", "RANTES", "eotaxin", "IL6"),
  list("IL10", "RA", "ALS"),
  list("CRP", "MCSF", "MS", "SCGFb"),
  list("CTACK", "IL18", "MIP1b", "IL8", "T2D")
  )

# cluster traits by region?
region_traits <- FALSE

# paths for cytokines and single GWAS
cytokines_directory <- "Harmonized GWAS results/cytokines"
single_GWAS_directory <- "Harmonized GWAS results/single_GWAS"
metabolite_directory <- "Harmonized GWAS results/metabolites"
cytokines_paths <- list.files(cytokines_directory, full.names = TRUE)
single_GWAS_paths <- list.files(single_GWAS_directory, full.names = TRUE)
metabolite_paths <- list.files(metabolite_directory, full.names = TRUE)

# group the paths by traits in lists
single_GWAS_traits_list <- group_paths_by_trait(single_GWAS_paths, fullnames = region_traits)
cytokine_traits_list <- group_paths_by_trait(cytokines_paths, fullnames = region_traits)
metabolite_traits_list <- group_paths_by_trait(metabolite_paths, fullnames = region_traits)

all_traits_list <- c(single_GWAS_traits_list, cytokine_traits_list, metabolite_traits_list)

# List of all betas and ses data frames for all groups of interest
merged_df_list <- list()

for (i in seq_along(trait_groups_of_interest)) {
  # get group of traits
  traits_of_interest <- trait_groups_of_interest[[i]]
  # String of trait names
  traits_group <- paste(traits_of_interest, collapse = ", ")
  print(paste("Creating merged betas and ses dataframes for: ", traits_group))
  
  # get paths of interest
  paths_of_interest <- get_paths_of_interest(all_traits_list, traits_of_interest)
  
  # create the merged data frames
  current_merged_dfs <- create_merged_betas_for_traits(paths_of_interest)
  
  # append to list of all groups
  merged_df_list[[traits_group]] <- current_merged_dfs
}

# List of all hyprcoloc results
res_list <- list()

for (i in seq_along(merged_df_list)) {
  # Get list of data frames
  current_dfs <- merged_df_list[[i]]
  
  # extract data frames
  merged_betas_matrix <- current_dfs[["betas"]]
  merged_ses_matrix <- current_dfs[["ses"]]
 
  # Arrange arguments for analysis function
  traits <- colnames(merged_betas_matrix)
  rsid <- rownames(merged_betas_matrix)
 
  # run hyprcoloc analysis with extracted data frames
  res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, snpscores = TRUE)

  sensitvity_plot <- sensitivity.plot(merged_betas_matrix, merged_ses_matrix, trait.names = traits, snp.id=rsid, 
                   reg.thresh = c(0.6,0.7,0.8,0.9), align.thresh = c(0.6,0.7,0.8,0.9), prior.c = c(0.02, 0.01, 0.005), equal.thresholds = FALSE, similarity.matrix = TRUE)
  
  output_list <- list(
    traits = traits,
    res = res,
    sens_plot = sensitvity_plot
  )

  # group
  traits_group <- paste(traits, collapse = ", ")
  
  # Add to result list 
  res_list[[traits_group]] <- output_list
}






