# Function to group paths together
group_paths_by_trait <- function(path_list){
  output_list <- list()
  
  for (i in seq_along(path_list)){
    trait_path <- path_list[i]
    # get trait name
    trait <- extract_trait(trait_path)
    # see if already present in list and create new entry or extend existing entry
    if(trait %in% names(output_list)){
      existing_traits <- output_list[[trait]]
      output_list[[trait]] <- c(existing_traits, trait_path)
    } else {
      output_list[[trait]] <- trait_path
    }
  }
  
  return(output_list)
}

# ----------------------------------------------------------------------

# Get subset paths list specified by traits from traits_to_paths_dictionary
get_paths_of_interest <- function(trait_path_dict, list_of_traits){
  output_list <- list()
  for (i in seq_along(list_of_traits)){
    trait <- list_of_traits[[i]] # get trait
    output_list[length(output_list)+1] <- trait_path_dict[trait] # Append path(s) specified by trait
  }
  
  return(output_list)
}

# ----------------------------------------------------------------------

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
  trait <- gsub(".*pmid[0-9]+_([A-Za-z0-9]+).*", "\\1", file_path)
  return(trait)
}
