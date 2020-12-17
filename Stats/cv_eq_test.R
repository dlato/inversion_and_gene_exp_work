#test to see if coef of var are different between groups
#https://cran.r-project.org/web/packages/cvequality/vignettes/how_to_test_CVs.html

library(cvequality)
library(dplyr)

#read in dat
g_dat <- read.table("Stats/inversions_gene_exp_wtest_data.csv", sep = "\t")
summary(g_dat)

#CV test between inver and non-inver (all blocks)
# null = no difference in CV between groups
#Feltz and Miller's test (1996)
with(g_dat, 
     asymptotic_test(norm_exp,
                     inversion))
#Krishnamoorthy and Lee's modified signed-likelihood ratio test (2014)
with(g_dat, 
     mslr_test(nr = 1e4, 
               norm_exp,
               inversion))

#CV test between rev_comp and non-rev_comp (all blocks)
# null = no difference in CV between groups
#Feltz and Miller's test (1996)
with(g_dat, 
     asymptotic_test(norm_exp,
                     rev_comp))
#Krishnamoorthy and Lee's modified signed-likelihood ratio test (2014)
with(g_dat, 
     mslr_test(nr = 1e4, 
               norm_exp,
               rev_comp))

#CV test between SIG inver and non-SIG inver
sig_d <- g_dat %>% filter(!is.na(block_w_pvalue)) %>%
                  mutate(sig = ifelse(block_w_pvalue <= 0.05, "yes", "no"))
# null = no difference in CV between groups
#Feltz and Miller's test (1996)
with(sig_d, 
     asymptotic_test(norm_exp,
                     sig))
#Krishnamoorthy and Lee's modified signed-likelihood ratio test (2014)
with(sig_d, 
     mslr_test(nr = 1e4, 
               norm_exp,
               sig))

#CV test between rev_comp and non-rev_comp (ATCC blocks)
atcc_d <- g_dat %>% filter(strain == "ATCC")
# null = no difference in CV between groups
#Feltz and Miller's test (1996)
with(atcc_d, 
     asymptotic_test(norm_exp,
                     rev_comp))
#Krishnamoorthy and Lee's modified signed-likelihood ratio test (2014)
with(atcc_d, 
     mslr_test(nr = 1e4, 
               norm_exp,
               rev_comp))
