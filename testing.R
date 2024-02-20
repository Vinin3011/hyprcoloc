
# test functions

ms_betas_unique2 <- create_betas_df(ms_paths) # betas_df function

# -----------------------------------------------
# Check if two dataframes are identical
identical <- identical(ms_betas_unique, ms_betas_unique2)

if (identical) {
  print("The data frames are identical.")
} else {
  print("The data frames are not identical.")
}

# -----------------------------------------------
# Check if values in SNP are unique and print the non-unique values if present

# Get counts of each SNP value ()
snp_counts <- table(dataframe_of_interest$SNP)

# Filter for non-unique values
non_unique_values <- names(snp_counts)[snp_counts > 1]

if (length(non_unique_values) == 0) {
  print("All values in the 'SNP' column are unique.")
} else {
  print("Non-unique values in the 'SNP' column along with their counts:")
  print(data.frame(SNP = non_unique_values, Count = snp_counts[non_unique_values]))
}

# ----------------------------------------------- 
cytokines_directory <- "Harmonized GWAS results/cytokines"

# Retrieve paths of all .tsv.gz files in the directory
cytokines_paths <- list.files(cytokines_directory, full.names = TRUE)

# Print the file paths
print(cytokines_paths)

test_cytokine <- read.table(gzfile(cytokines_paths[1]), header = TRUE, sep = "\t")


# -------------------------------------------------
# read metabolite tsv files
test_betas_list <- list()

for (i in seq_along(metabolite_traits_list)){
  sublist <- metabolite_traits_list[[i]]
  
  # Get trait name
  trait <- extract_trait(sublist[1])
  print(paste("Creating betas and ses dataframe for: ", trait))
  
  result_betas_df <- tryCatch(create_betas_df(sublist), error = function(e){
    print(paste("create_betas_df produced the following error for trait ", trait, ":",e))
    print("continue with next trait...")
    return(NULL)  # Return NULL to assign an empty dataframe
  })
  
  # Skip the rest of the loop if an error occurred
  if(is.null(result_betas_df)){
    next  
  }
  
  test_betas_list[[trait]] <- result_betas_df
}


# -------------------------------------------------

# cleanup
rm(identical, ms_betas_unique2)
rm(file_paths)
rm(test_betas_list)


