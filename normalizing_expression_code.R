library("edgeR")

#all samples should be apart of the same group
#make sure that your columns showing presence or absene of the gene is NOT included in the df that is put through this code!
num_samples <- length(colnames(exp_dat_all_reps))
group <- rep(1, num_samples)

#normalize all replicates and samples
exp_list <- DGEList(counts=exp_dat_all_reps, group=group)
exp_normalized <- calcNormFactors(exp_list)
normalized_exp_df <- cpm(exp_normalized, normalized.lib.sizes=FALSE)


