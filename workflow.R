# hyprcoloc workflow
traits <- paste0("T", 1:dim(merged_betas_first1000_matrix)[2])
rsid <- rownames(merged_betas_first1000_matrix)
res <- hyprcoloc(merged_betas_first1000_matrix, merged_ses_first1000_matrix, trait.names=traits, snp.id=rsid)
res


betas <- hyprcoloc::test.betas
ses <- hyprcoloc::test.ses

#hallo
