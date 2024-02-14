# hyprcoloc workflow
traits <- paste0("T", 1:dim(merged_betas)[2])
rsid <- rownames(merged_betas)
res <- hyprcoloc(merged_betas, merged_ses, trait.names=traits, snp.id=rsid)
res
