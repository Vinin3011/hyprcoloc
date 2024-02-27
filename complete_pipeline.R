library(hyprcoloc)
library(dplyr)
library(data.table)
library(purrr)

trait_groups_of_interest <- list(
  list("bNGF_27989323-GCST004421-EFO_0008035", "pmid33987465_PD_eur_full", "ALS_80610_29566793-GCST005647-EFO_0000253", "pmid21833088_MS_eur_full", "RA_311292_33310728-GCST90013534-EFO_0000685", "eotaxin_27989323-GCST004460-EFO_0008122", "TNFa_27989323-GCST004426-EFO_0004684", "T2D_659316_30054458-GCST006867-EFO_0001360",
       "pmid24816252_asparagine_eur_pval", "MCP1_27989323-GCST004438-EFO_0004749")
  # ,list("eotaxin_27989323-GCST004460-EFO_0008122", "MCP1_27989323-GCST004438-EFO_0004749", "IL6_27989323-GCST004446-EFO_0004810", "IFNg_27989323-GCST004456-EFO_0008165")
)
# cluster traits by region?
region_traits <- FALSE
complete_name <- TRUE
new_gwas <- TRUE

# paths for cytokines and single GWAS
cytokines_directory <- "Harmonized GWAS results/cytokines_new"
single_GWAS_directory <- "Harmonized GWAS results/full_unfiltered"
metabolite_directory <- "Harmonized GWAS results/metabolites"
cytokines_paths <- list.files(cytokines_directory, full.names = TRUE)
single_GWAS_paths <- list.files(single_GWAS_directory, full.names = TRUE)
metabolite_paths <- list.files(metabolite_directory, full.names = TRUE)

# group the paths by traits in lists
single_GWAS_traits_list <- group_paths_by_trait(single_GWAS_paths, fullnames = region_traits, complete_names = complete_name)
cytokine_traits_list <- group_paths_by_trait(cytokines_paths, fullnames = region_traits, complete_names = complete_name)
metabolite_traits_list <- group_paths_by_trait(metabolite_paths, fullnames = region_traits, complete_names = complete_name)

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
  current_merged_dfs <- create_merged_dfs_for_traits(paths_of_interest, complete_name = complete_name, new_gwas = new_gwas)
  
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

trait_dataframe <- data.frame(traits = traits)
rownames(trait_dataframe) <- traits

similarity_matrices <- list()
for (entry in res_list) {
  matrix <- entry$sens_plot[[2]]
  print(matrix)
  similarity_matrices <- c(similarity_matrices, list(matrix))
}

# Create an empty list to store data frames
similarity_data_frames <- list()

# Iterate over each matrix in the list
for (matrix in similarity_matrices) {
  # Convert the matrix into a data frame
  matrix_df <- as.data.frame(matrix)
  
  # Append the data frame to the list
  similarity_data_frames <- c(similarity_data_frames, list(matrix_df))
}

# Convert trait_dataframe to a matrix to preserve row names
trait_matrix <- as.matrix(trait_dataframe)

# Iterate over each similarity matrix in similarity_data_frames
for(i in seq_along(similarity_data_frames)) {
  # Subset the current similarity matrix to match the row names of trait_matrix
  subset_similarity <- similarity_data_frames[[i]][rownames(trait_matrix), ]
  
  # Combine trait_matrix and subset_similarity using cbind
  trait_matrix <- cbind(trait_matrix, subset_similarity)
  
}

merged_df <- as.data.frame(trait_matrix)
#merged_df$traits <- NULL 
merged_df <- merged_df[!rowSums(is.na(merged_df)) == ncol(merged_df), ]

# Remove non-numeric columns
merged_df_numeric <- merged_df[, sapply(merged_df, is.numeric)]

my_palette <- c("white","blue", "darkblue")

# Plot heatmap for the current merged dataframe
heatmap(as.matrix(merged_df_numeric), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette,
        main = paste("Heatmap for Similarity Matrix", i))

