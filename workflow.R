# hyprcoloc workflow
traits <- colnames(merged_betas_matrix)
rsid <- rownames(merged_betas_matrix)

res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid)


# 2.1.1 Labelling a trait as either continuous or binary in analyses

binary.traits <- rep(0, ncol(merged_betas_matrix))
res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, binary.outcomes = binary.traits)

# 2.1.2 Choosing a subset of traits to analyze

res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, trait.subset = c("MS","ALS"))

# 2.2 Computing a credible set of snps for each cluster of colocalized traits

traits <- colnames(merged_betas_matrix)
rsid <- rownames(merged_betas_matrix)
res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, snpscores = TRUE)
res[[1]]

cred.sets(res, value = 0.95);

# > cred.sets(res, value = 0.95);
# [[1]]
# rs9633340  rs11265206   rs4326616 rs149783572 rs183012206  rs75431404  rs10908724   rs4656826 
# 0.02998592  0.02852620  0.02788553  0.02640320  0.02583718  0.02474243  0.02429282  0.02407639 
# rs12410729  rs79508328  rs79433881   rs4656237   rs4399156  rs12142553  rs12118628  rs11265187 
# 0.02407639  0.02407639  0.02356134  0.02246796  0.02243733  0.02243733  0.02243733  0.02195802 
# rs12035635   rs4551589   rs4290055  rs74571338  rs74355389   rs4656230   rs4656232   rs4656820 
# 0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403 
# rs78879802  rs11265185  rs11265186  rs10908716  rs10908717   rs4656821  rs12131860  rs12144843 
# 0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403 
# rs12144848  rs12042314   rs4492615  rs11265193  rs75906889   rs4128725  rs12145616  rs12129890 
# 0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403  0.02183403 
# rs74520022  rs76546294 
# 0.02183403  0.01933677 
# 
# [[2]]
# rs72831623 rs13412535 rs68066031 
# 0.7570090  0.1352653  0.1026782 

# 3.1.0.1 Conditionally uniform configuration priors

prior.options = c(1e-4, 1e-10, 1e-20, 1e-25, 1e-100);
for(i in prior.options){
  res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, 
                   uniform.priors = TRUE, prior.1 = i, reg.steps = 2);
  print(paste0("prior.1 = ",i));
  print(res);
}

# [1] "prior.1 = 1e-04"
# $results
# iteration        traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp
# 1         1 MCP1, eotaxin         0.9970         1.000     rs9633340                      0.030
# 2         2     IL6, IFNg         0.9375         0.947    rs72831623                      0.757
# dropped_trait
# 1            NA
# 2            NA
# 
# attr(,"class")
# [1] "hyprcoloc"
# [1] "prior.1 = 1e-10"
# $results
# iteration        traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp
# 1         1 MCP1, eotaxin         0.9961        0.9991     rs9633340                       0.03
# 2         2          None             NA        0.0000          <NA>                         NA
# dropped_trait
# 1          <NA>
#   2           IL6
# 
# attr(,"class")
# [1] "hyprcoloc"
# [1] "prior.1 = 1e-20"
# $results
# iteration traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp dropped_trait
# 1         1   None             NA             0            NA                         NA          MCP1
# 2         2   None             NA             0            NA                         NA           IL6
# 3         3   None             NA             0            NA                         NA       eotaxin
# 
# attr(,"class")
# [1] "hyprcoloc"
# [1] "prior.1 = 1e-25"
# $results
# iteration traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp dropped_trait
# 1         1   None             NA             0            NA                         NA          MCP1
# 2         2   None             NA             0            NA                         NA           IL6
# 3         3   None             NA             0            NA                         NA       eotaxin
# 
# attr(,"class")
# [1] "hyprcoloc"
# [1] "prior.1 = 1e-100"
# $results
# iteration traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp dropped_trait
# 1         1   None             NA             0            NA                         NA          MCP1
# 2         2   None             NA             0            NA                         NA           IL6
# 3         3   None             NA             0            NA                         NA       eotaxin


# 3.1.0.2 Variant specific configuration priors 

prior.1 = 1e-4;
prior.c.options = c(0.05, 0.02, 0.01, 0.005);
for(i in prior.c.options){
  res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid,
                   uniform.priors = FALSE, prior.1 = prior.1, prior.c = i);
  print(c(paste0("prior.1 = ",prior.1), paste0("prior.c = ",i)));
  print(res);
}

# /*
# [1] "prior.1 = 1e-04" "prior.c = 0.05" 
# $results
# iteration        traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp
# 1         1 MCP1, eotaxin         0.9552         1.000     rs9633340                      0.030
# 2         2     IL6, IFNg         0.8538         0.986    rs72831623                      0.757
# dropped_trait
# 1            NA
# 2            NA
# 
# attr(,"class")
# [1] "hyprcoloc"
# [1] "prior.1 = 1e-04" "prior.c = 0.02" 
# $results
# iteration        traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp
# 1         1 MCP1, eotaxin         0.8950        0.9999     rs9633340                      0.030
# 2         2     IL6, IFNg         0.6962        0.9658    rs72831623                      0.757
# dropped_trait
# 1            NA
# 2            NA
# 
# attr(,"class")
# [1] "hyprcoloc"
# [1] "prior.1 = 1e-04" "prior.c = 0.01" 
# $results
# iteration        traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp
# 1         1 MCP1, eotaxin         0.8100        0.9998     rs9633340                      0.030
# 2         2     IL6, IFNg         0.5263        0.9338    rs72831623                      0.757
# dropped_trait
# 1            NA
# 2            NA
# 
# attr(,"class")
# [1] "hyprcoloc"
# [1] "prior.1 = 1e-04" "prior.c = 0.005"
# $results
# iteration        traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp
# 1         1 MCP1, eotaxin         0.6807        0.9997     rs9633340                       0.03
# 2         2          None             NA        0.8759          <NA>                         NA
# dropped_trait
# 1        <NA>
# 2         IL6

# 3.1.0.3 Evidential strength and the regional and alignment thresholds

res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, uniform.priors = FALSE,
                 reg.thresh = 0.95, align.thresh = 0.95)

# $results
# iteration traits posterior_prob regional_prob candidate_snp posterior_explained_by_snp dropped_trait
# 1         1   None             NA        0.9999            NA                         NA       eotaxin
# 2         2   None             NA        0.0286            NA                         NA          MCP1
# 3         3   None             NA        0.9658            NA                         NA           IL6

# 3.1.0.4 Analysis protocol: a one stop sensitivity assessment
sensitivity.plot(merged_betas_matrix, merged_ses_matrix, trait.names = traits, snp.id=rsid, 
                 reg.thresh = c(0.6,0.7,0.8,0.9), align.thresh = c(0.6,0.7,0.8,0.9), prior.c = c(0.02, 0.01, 0.005), equal.thresholds = FALSE)

sensitivity.plot(merged_betas_matrix, merged_ses_matrix, trait.names = traits, snp.id=rsid, 
                 reg.thresh = c(0.6,0.7,0.8,0.9), align.thresh = c(0.6,0.7,0.8,0.9), prior.c = c(0.02, 0.01, 0.005), equal.thresholds = FALSE, similarity.matrix = TRUE)
res

# 3.2.1 Assessing differences between the Bayesian divisive clustering criteria
res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, bb.selection = "regional")
print(res)

res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, bb.selection = "align")
print(res)

res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, bb.selection = FALSE)
print(res)


traits <- colnames(merged_betas_matrix)
rsid <- rownames(merged_betas_matrix)
res <- hyprcoloc(merged_betas_matrix, merged_ses_matrix, trait.names=traits, snp.id=rsid, snpscores = TRUE)
res[[1]]
sensitivity.plot(merged_betas_matrix, merged_ses_matrix, trait.names = traits, snp.id=rsid, 
                 reg.thresh = c(0.6,0.7,0.8,0.9), align.thresh = c(0.6,0.7,0.8,0.9), prior.c = c(0.02, 0.01, 0.005), equal.thresholds = FALSE)

