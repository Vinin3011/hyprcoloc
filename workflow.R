library(hyprcoloc)
library(dplyr)
hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)

#data1 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid21833088_MS_eur.tsv.gz"), header = TRUE, sep = "\t")
#data2 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid24076602_MS_eur.tsv.gz"), header = TRUE, sep = "\t")
#data3 <- read.table(gzfile("Harmonized GWAS results/single_GWAS/pmid27386562_MS_eur.tsv.gz"), header = TRUE, sep = "\t")

ms_paths <- c("Harmonized GWAS results/single_GWAS/pmid21833088_MS_eur.tsv.gz", "Harmonized GWAS results/single_GWAS/pmid24076602_MS_eur.tsv.gz", "Harmonized GWAS results/single_GWAS/pmid27386562_MS_eur.tsv.gz")

# Read each .gz file and bind the resulting data frames together
ms_betas <- lapply(ms_paths, function(ms_paths) {
  # Read the data and explicitly specify the column types to ensure consistency
  data <- read.table(gzfile(ms_paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(beta_effect = "numeric")) 
  new_column_name <- as.character(data$trait[1])
  data %>%
    select(SNP, beta_effect) %>%
    rename(!!new_column_name := beta_effect)
}) %>% bind_rows()






