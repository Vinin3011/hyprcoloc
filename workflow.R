library(hyprcoloc)
library(dplyr)
library(data.table)
hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)

#data1 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid21833088_MS_eur.tsv.gz"), header = TRUE, sep = "\t")
#data2 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid24076602_MS_eur.tsv.gz"), header = TRUE, sep = "\t")
#data3 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid27386562_MS_eur.tsv.gz"), header = TRUE, sep = "\t")

ms_paths <- c("Harmonized GWAS results/single_GWAS/pmid21833088_MS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid24076602_MS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid27386562_MS_eur.tsv.gz")

als_paths <-c("Harmonized GWAS results/single_GWAS/pmid27455348_ALS_eur.tsv.gz", 
              "Harmonized GWAS results/single_GWAS/pmid28931804_ALS_eur-eas.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid29566793_ALS.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid29566793_ALS_eur.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid34873335_ALS_eur.tsv.gz",
              "Harmonized GWAS results/single_GWAS/pmid34873335_ALS_eur-eas.tsv.gz")
t2d_paths <-c()
pd_paths <-c()



# Read each .gz file and bind the resulting data frames together
ms_betas <- lapply(ms_paths, function(ms_paths) {
  # Read the data and explicitly specify the column types to ensure consistency
  data <- read.table(gzfile(ms_paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(beta_effect = "numeric")) 
  new_column_name <- as.character(data$trait[1])
  data %>%
    select(SNP, beta_effect) %>%
    rename(!!new_column_name := beta_effect)
}) %>% bind_rows()


# Get counts of each SNP value
snp_counts <- table(ms_betas$SNP)

# Filter for non-unique values
non_unique_values <- names(snp_counts)[snp_counts > 1]

if (length(non_unique_values) == 0) {
  print("All values in the 'SNP' column are unique.")
} else {
  print("Non-unique values in the 'SNP' column along with their counts:")
  print(data.frame(SNP = non_unique_values, Count = snp_counts[non_unique_values]))
}
