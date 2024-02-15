# Function to extract betas from GWAS tsv.gz data defined by their paths
create_betas_df <- function(paths){
  # Read each .gz file and bind the resulting data frames together
  output <- lapply(paths, function(paths) {
    # Read the data and explicitly specify the column types to ensure consistency
    data <- read.table(gzfile(paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(beta_effect = "numeric"))
    new_column_name <- as.character(data$trait[1])
    data %>%
      select(SNP, beta_effect) %>%
      rename(!!new_column_name := beta_effect)
  }) %>% bind_rows() # combine into single data frame
  
  # remove NA values
  output <- output[output$SNP != "NA", ]
  
  # remove duplicates
  output <- output[!duplicated(output$SNP), ]
  
  return(output)
}

# -------------------------------------------------------------------

# Function to extract ses from GWAS tsv.gz data defined by their paths
create_ses_df <- function(paths){
  # Read each .gz file and bind the resulting data frames together
  output <- lapply(paths, function(paths) {
    # Read the data and explicitly specify the column types to ensure consistency
    data <- read.table(gzfile(paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(SE = "numeric"))
    new_column_name <- as.character(data$trait[1])
    data %>%
      select(SNP, SE) %>%
      rename(!!new_column_name := SE)
  }) %>% bind_rows() # combine into single data frame
  
  # remove NA values
  output <- output[output$SNP != "NA", ]
  
  # remove duplicates
  output <- output[!duplicated(output$SNP), ]
  
  return(output)
}

# -------------------------------------------------------------------

# Function to extract the trait from a standard GWAS file path name
extract_trait <- function(file_path) {
  # Extract the trait abbreviation using regular expression
  trait <- gsub(".*pmid[0-9]+_([A-Za-z0-9]+)_.*", "\\1", file_path)
  return(trait)
}
