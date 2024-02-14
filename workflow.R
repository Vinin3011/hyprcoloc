library(hyprcoloc)
library(dplyr)
library(data.table)

#data1 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid21833088_MS_eur.tsv.gz"), header = TRUE, sep = "\t")
#data2 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid24076602_MS_eur.tsv.gz"), header = TRUE, sep = "\t")
#data3 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid27386562_MS_eur.tsv.gz"), header = TRUE, sep = "\t")

ms_paths <- c("Harmonized GWAS results/single_GWAS/pmid21833088_MS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid24076602_MS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid27386562_MS_eur.tsv.gz")

als_paths <- c("Harmonized GWAS results/single_GWAS/pmid27455348_ALS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid28931804_ALS_eur-eas.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid29566793_ALS.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid29566793_ALS_eur.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid34873335_ALS_eur.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid34873335_ALS_eur-eas.tsv.gz")

t2d_paths <-c()
pd_paths <-c()

#datax <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid21833088_MS_eur.tsv.gz"))


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

# Merge the betas data frames based on the "SNP" column
merged_ses <- merge(ms_ses_unique, als_ses_unique, by = "SNP")

# hyprcoloc workflow
traits <- paste0("T", 1:dim(merged_betas)[2])
rsid <- rownames(merged_betas)
res <- hyprcoloc(merged_betas, merged_ses, trait.names=traits, snp.id=rsid)
res
