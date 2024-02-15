library(hyprcoloc)
library(dplyr)
library(data.table)

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

path_array


# MS
# betas table creation

# Read each .gz file and bind the resulting data frames together
ms_betas <- lapply(ms_paths, function(ms_paths) {
  # Read the data and explicitly specify the column types to ensure consistency
  data <- read.table(gzfile(ms_paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(beta_effect = "numeric"))
  new_column_name <- as.character(data$trait[1])
  data %>%
    select(SNP, beta_effect) %>%
    rename(!!new_column_name := beta_effect)
}) %>% bind_rows() # combine into single data frame

# remove NA values
ms_betas <- ms_betas[ms_betas$SNP != "NA", ]

# remove duplicates
ms_betas_unique <- ms_betas[!duplicated(ms_betas$SNP), ]

# Check if values in SNP are unique and print the non-unique values if present

# Get counts of each SNP value
snp_counts <- table(ms_betas_unique$SNP)

# Filter for non-unique values
non_unique_values <- names(snp_counts)[snp_counts > 1]

if (length(non_unique_values) == 0) {
  print("All values in the 'SNP' column are unique.")
} else {
  print("Non-unique values in the 'SNP' column along with their counts:")
  print(data.frame(SNP = non_unique_values, Count = snp_counts[non_unique_values]))
}

# se table creation
# Read each .gz file and bind the resulting data frames together
ms_ses <- lapply(ms_paths, function(ms_paths) {
  # Read the data and explicitly specify the column types to ensure consistency
  data <- read.table(gzfile(ms_paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(SE = "numeric"))
  new_column_name <- as.character(data$trait[1])
  data %>%
    select(SNP, SE) %>%
    rename(!!new_column_name := SE)
}) %>% bind_rows() # combine into single data frame

# remove NA values
ms_ses <- ms_ses[ms_ses$SNP != "NA", ]

# remove duplicates
ms_ses_unique <- ms_ses[!duplicated(ms_ses$SNP), ]

# Check if values in SNP are unique and print the non-unique values if present

# Get counts of each SNP value
snp_counts <- table(ms_ses_unique$SNP)

# Filter for non-unique values
non_unique_values <- names(snp_counts)[snp_counts > 1]

if (length(non_unique_values) == 0) {
  print("All values in the 'SNP' column are unique.")
} else {
  print("Non-unique values in the 'SNP' column along with their counts:")
  print(data.frame(SNP = non_unique_values, Count = snp_counts[non_unique_values]))
}

# ALS
# betas table creation
# Read each .gz file and bind the resulting data frames together
als_betas <- lapply(als_paths, function(als_paths) {
  # Read the data and explicitly specify the column types to ensure consistency
  data <- read.table(gzfile(als_paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(beta_effect = "numeric"))
  new_column_name <- as.character(data$trait[1])
  data %>%
    select(SNP, beta_effect) %>%
    rename(!!new_column_name := beta_effect)
}) %>% bind_rows() # combine into single data frame

# remove NA values
als_betas <- als_betas[als_betas$SNP != "NA", ]

# remove duplicates
als_betas_unique <- als_betas[!duplicated(als_betas$SNP), ]

# Check if values in SNP are unique and print the non-unique values if present

# Get counts of each SNP value
snp_counts <- table(als_betas_unique$SNP)

# Filter for non-unique values
non_unique_values <- names(snp_counts)[snp_counts > 1]

if (length(non_unique_values) == 0) {
  print("All values in the 'SNP' column are unique.")
} else {
  print("Non-unique values in the 'SNP' column along with their counts:")
  print(data.frame(SNP = non_unique_values, Count = snp_counts[non_unique_values]))
}

# se table creation
# Read each .gz file and bind the resulting data frames together
als_ses <- lapply(als_paths, function(als_paths) {
  # Read the data and explicitly specify the column types to ensure consistency
  data <- read.table(gzfile(als_paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(SE = "numeric"))
  new_column_name <- as.character(data$trait[1])
  data %>%
    select(SNP, SE) %>%
    rename(!!new_column_name := SE)
}) %>% bind_rows() # combine into single data frame

# remove NA values
als_ses <- als_ses[als_ses$SNP != "NA", ]

# remove duplicates
als_ses_unique <- als_ses[!duplicated(als_ses$SNP), ]

# Check if values in SNP are unique and print the non-unique values if present

# Get counts of each SNP value
snp_counts <- table(als_ses_unique$SNP)

# Filter for non-unique values
non_unique_values <- names(snp_counts)[snp_counts > 1]

if (length(non_unique_values) == 0) {
  print("All values in the 'SNP' column are unique.")
} else {
  print("Non-unique values in the 'SNP' column along with their counts:")
  print(data.frame(SNP = non_unique_values, Count = snp_counts[non_unique_values]))
}


# Merge the betas data frames based on the "SNP" column
merged_betas <- merge(ms_betas_unique, als_betas_unique, by = "SNP")

merged_betas_first1000 <- head(merged_betas, 1000)

# Merge the betas data frames based on the "SNP" column
merged_ses <- merge(ms_ses_unique, als_ses_unique, by = "SNP")

merged_ses_nonzero <- merged_ses[merged_ses$MS != 0, ]
merged_ses_first1000 <- head(merged_ses_nonzero, 1000)

# Convert into matrix

# betas

# Extract row names from the first column
rownames <- merged_betas_first1000[, 1]

# Remove the first column before converting to matrix
merged_betas_first1000 <- merged_betas_first1000[, -1]

# Convert dataframe to matrix and set row names
merged_betas_first1000_matrix <- as.matrix(merged_betas_first1000)
rownames(merged_betas_first1000_matrix) <- rownames

# se

# Extract row names from the first column
rownames <- merged_ses_first1000[, 1]

# Remove the first column before converting to matrix
merged_ses_first1000 <- merged_ses_first1000[, -1]

# Convert dataframe to matrix and set row names
merged_ses_first1000_matrix <- as.matrix(merged_ses_first1000)
rownames(merged_ses_first1000_matrix) <- rownames