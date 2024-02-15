# hyprcoloc workflow
traits <- colnames(merged_betas_matrix)
rsid <- rownames(merged_betas_matrix)
res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid)


# 2.1.1 Labelling a trait as either continuous or binary in analyses

binary.traits <- rep(1, ncol(merged_betas_matrix))
res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, binary.outcomes = binary.traits)

# 2.1.2 Choosing a subset of traits to analyze

res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, trait.subset = c("MS","ALS"))

# 2.2 Computing a credible set of snps for each cluster of colocalized traits

res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, snpscores = TRUE)

res[[2]][[1]]

#cred.sets(res, value = 0.95)

#> [[1]]
#> rs11591147 
#>          1 
#> 
#> [[2]]
#>  rs12117612   rs7532349  rs11206481  rs12145624  rs12724445  rs12126037 
#> 0.419677651 0.419677651 0.025324539 0.022969148 0.022969148 0.022332481 
#>  rs13375783   rs1544909 
#> 0.014513746 0.009991906 

# 3.1.0.1 Conditionally uniform configuration priors
#fehlender parameter reg.steps = 2

prior.options = c(1e-4, 1e-10, 1e-20, 1e-25, 1e-100);
for(i in prior.options){
  res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, 
                   uniform.priors = TRUE, prior.1 = i);
  print(paste0("prior.1 = ",i));
  print(res);
}

