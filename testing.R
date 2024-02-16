
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
# cleanup
rm(identical, ms_betas_unique2)
rm(file_paths)


