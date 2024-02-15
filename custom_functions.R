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