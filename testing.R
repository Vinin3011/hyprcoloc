ms_betas_unique2 <- create_betas_df(ms_paths)

identical <- identical(ms_betas_unique, ms_betas_unique2)

if (identical) {
  print("The data frames are identical.")
} else {
  print("The data frames are not identical.")
}

# cleanup
rm(identical, ms_betas_unique2)
