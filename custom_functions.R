# Function to create from a list of paths corresponding to GWAS studies of traits of interest
# a merged betas and merged ses data frame

create_merged_dfs_for_traits <- function(paths_of_interest, fullname = FALSE, complete_name = FALSE, new_gwas = FALSE){
  # Create lists to fill with betas and ses df
  betas_df_list <- list()
  ses_df_list <- list()
  traitnames_list <- list()
  
  
  # Create ses and betas df for every trait of interest
  for (i in seq_along(paths_of_interest)) {
    sublist <- paths_of_interest[[i]]  # Get the sublist
    
    # Get trait name
    trait <- extract_trait(sublist[1], fullname, complete_name)
    print(paste("Creating betas and ses dataframe for: ", trait))
    # Store trait in list
    traitnames_list[length(traitnames_list)+1] <- trait
    
    # attempt to create betas df from sublist and continue if an error occurs
    result_betas_df <- tryCatch(create_betas_df(sublist, new_gwas = new_gwas), error = function(e){
      print(paste("create_betas_df produced the following error for trait ", trait, ":",e))
      print("continue with next trait...")
      return(NULL)  # Return NULL to assign an empty dataframe
    })
    
    # Skip the rest of the loop if an error occurred
    if(is.null(result_betas_df)){
      next  
    }
    
    # attempt to create betas df from sublist and continue if an error occurs
    result_ses_df <- tryCatch(create_ses_df(sublist, new_gwas = new_gwas), error = function(e){
      print(paste("create_betas_df produced the following error for trait ", trait, ":",e))
      print("continue with next trait...")
      return(NULL)  # Return NULL to assign an empty dataframe
    })
    
    # Skip the rest of the loop if an error occurred
    if(is.null(result_ses_df)){
      next  
    }
    
    # Store the result dataframe along with the trait name in the results list
    betas_df_list[[trait]] <- result_betas_df
    ses_df_list[[trait]] <- result_ses_df
  }
  
  # Merge the betas data frames based on the "SNP" column
  merged_betas <- reduce(betas_df_list, inner_join, by = "SNP")
  
  # Merge the ses data frames based on the "SNP" column
  merged_ses <- reduce(ses_df_list, inner_join, by = "SNP")
  
  # Get the column names except the first one (assuming it's "SNP")
  cols <- names(merged_ses)[-1]
  
  # Remove rows with zero values in any column except the first one
  merged_ses_nonzero <- merged_ses[rowSums(merged_ses[cols] == 0) == 0, ]
  merged_betas_nonzero <- merged_betas[rowSums(merged_betas[cols] == 0) == 0, ]
  
  
  # Convert into matrices ------------------------------------
  
  # betas 
  
  # Extract row names from the first column
  rownames <- merged_betas_nonzero[, 1]
  # Remove the first column before converting to matrix
  merged_betas_nonzero <- merged_betas_nonzero[, -1]
  # Convert dataframe to matrix and set row names
  merged_betas_matrix <- as.matrix(merged_betas_nonzero)
  rownames(merged_betas_matrix) <- rownames
  
  
  # ses 
  
  # Extract row names from the first column
  rownames <- merged_ses_nonzero[, 1]
  # Remove the first column before converting to matrix
  merged_ses_nonzero <- merged_ses_nonzero[, -1]
  # Convert dataframe to matrix and set row names
  merged_ses_matrix <- as.matrix(merged_ses_nonzero)
  rownames(merged_ses_matrix) <- rownames
  
  # output a list with traitnames and both merged matrices (betas and ses)
  output <- list(
    traits = traitnames_list,
    ses = merged_ses_matrix,
    betas = merged_betas_matrix
  )
  
  return(output)
}

# -------------------------------------------------------------------



# Function to group paths together
group_paths_by_trait <- function(path_list, fullnames = FALSE, complete_names = FALSE){
  output_list <- list()
  
  for (i in seq_along(path_list)){
    trait_path <- path_list[i]
    # get trait name
    trait <- extract_trait(trait_path, fullname = fullnames, complete_name = complete_names)
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
create_betas_df <- function(paths, new_gwas = FALSE){
  # Read each .gz file and bind the resulting data frames together
  output <- lapply(paths, function(paths) {
    # Read the data and explicitly specify the column types to ensure consistency
    data <- read.table(gzfile(paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(beta_effect = "numeric"))
    # If new GWAS data. Retrieve trait from path, since there is no trait column
    if(new_gwas){
      new_column_name <- extract_trait(paths, new_gwas = TRUE)
    }else{
      new_column_name <- as.character(data$trait[1])
    }
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
create_ses_df <- function(paths, new_gwas = FALSE){
  # Read each .gz file and bind the resulting data frames together
  output <- lapply(paths, function(paths) {
    # Read the data and explicitly specify the column types to ensure consistency
    data <- read.table(gzfile(paths), header = TRUE, sep = "\t",na.strings = "-NA", stringsAsFactors = FALSE, colClasses = c(SE = "numeric"))
    # If new GWAS data. Retrieve trait from path, since there is no trait column
    if(new_gwas){
      new_column_name <- extract_trait(paths, new_gwas = TRUE)
    }else{
      new_column_name <- as.character(data$trait[1])
    }
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
extract_trait <- function(file_path, fullname = FALSE, complete_name = FALSE, new_gwas = FALSE) {
  # Get the basename of the path
  base_name <- basename(as.character(file_path))
  
  # Specify if we want to be more specific for the trait name
  pat <- ".*pmid[0-9]+_([A-Za-z0-9]+).*"
  if(fullname){
    pat <- ".*pmid[0-9]+_([^\\.]+)\\..*" # capture everything after pmid to first dot
  }
   else if(complete_name){
    pat <- ".*/([^/.]+)\\..*" # capture everything until dot
  }
  else if(new_gwas){
    # If pathbase starts with pmid just keep specific pattern
    # Else just take the start of the string until firts underscore
    if(!startsWith(base_name,"pmid")){
      pat <- "^([^_]+)_.*"
    }
  }
  
  # Extract the trait abbreviation using regular expression
  trait <- gsub(pat, "\\1", base_name)
  return(trait)
}

