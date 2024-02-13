install.packages("devtools")
library(devtools)
install_github("cnfoley/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE)
library(hyprcoloc)
browseVignettes("hyprcoloc")

