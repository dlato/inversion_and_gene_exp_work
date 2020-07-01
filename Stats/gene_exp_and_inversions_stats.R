#test script for inversions and gene expression analysis
library(tidyverse)
library(dplyr)
library(tidyr)
library(broom)
library(rstatix)
library(ggpubr)

#read in the data
gei_dat <- read.csv("Sample_final_df.csv", header = TRUE)

############
# ALL inversions
############
#Wilcox signed-rank test
# to see if there is a difference btwn inversions and non-inverted segments
#null hyp = there is no difference in gene expression btwn inverted and non-inverted segments
wilcox.test(gei_dat$norm_exp[gei_dat$inversion == 1], gei_dat$norm_exp[gei_dat$inversion == 0])
t.test(gei_dat$norm_exp[gei_dat$inversion == 1], gei_dat$norm_exp[gei_dat$inversion == 0])

############
# by block inversions
############
uniq_block_w_pvalue <- c()
uniq_block_t_pvalue <- c()
uniq_block_t_stat <- c()
uniq_block_t_confinf <- c()
for (i in unique(gei_dat$block)) {
  tmp_df <- gei_dat[which(gei_dat$block == i),]
  if (length(unique(tmp_df$rev_comp)) > 1) {
    #wilcox test
    wilcox_test <- wilcox.test(tmp_df$norm_exp[tmp_df$rev_comp == 1], tmp_df$norm_exp[tmp_df$rev_comp == 0])
    uniq_block_w_pvalue <- c(uniq_block_w_pvalue,wilcox_test$p.value)
    #t test
    if (length(which(tmp_df$rev_comp == 1)) > 1 & length(which(tmp_df$rev_comp == 0)) > 1){
      t_test <- t.test(tmp_df$norm_exp[tmp_df$rev_comp == 1], tmp_df$norm_exp[tmp_df$rev_comp == 0])
      uniq_block_t_pvalue <- c(uniq_block_t_pvalue,t_test$pvalue)
      uniq_block_t_stat <- c(uniq_block_t_stat,t_test$statistic)
      uniq_block_t_confinf <- c(uniq_block_t_confinf,t_test$conf.inf)
    } else {
      uniq_block_t_pvalue <- c(uniq_block_t_pvalue,"NA")
      uniq_block_t_stat <- c(uniq_block_t_stat,"NA")
      uniq_block_t_confinf <- c(uniq_block_t_confinf,"NA")
    }
  } else {
    #not able to do tests
    uniq_block_w_pvalue <- c(uniq_block_w_pvalue,"NA")
    uniq_block_t_pvalue <- c(uniq_block_t_pvalue,"NA")
    uniq_block_t_stat <- c(uniq_block_t_stat,"NA")
    uniq_block_t_confinf <- c(uniq_block_t_confinf,"NA")
  }
}

#add test results to df
count = 1
block_w_pvalue <- vector(mode='numeric', length = length(gei_dat$block))
block_t_pvalue <- vector(mode='numeric', length = length(gei_dat$block))
block_t_stat <- vector(mode='numeric', length = length(gei_dat$block))
block_t_coninf <- vector(mode='numeric', length = length(gei_dat$block))
for (b in unique(gei_dat$block)) {
  bloc_loc <- which(gei_dat$block == b)
  block_w_pvalue[bloc_loc] <- uniq_block_w_pvalue[count] 
  block_t_pvalue[bloc_loc] <- uniq_block_t_pvalue[count] 
  block_t_stat[bloc_loc] <- uniq_block_t_stat[count] 
  block_t_coninf[bloc_loc] <- uniq_block_t_confinf[count] 
  count <- count + 1
}
gei_dat['block_w_pvalue'] <- block_w_pvalue
gei_dat['block_t_pvalue'] <- block_t_pvalue
gei_dat['block_t_stat'] <- block_t_stat
gei_dat['block_t_confinf'] <- block_t_coninf





#t.test(gei_dat$norm_exp[gei_dat$block == "Block800",], gei_dat$norm_exp[gei_dat$block == "Block800",])
#
#test_inver <- gei_dat %>% filter(inversion == 1)
#test_revcomp <- gei_dat %>% filter(rev_comp == 1)
#
#gei_dat %>% filter(block == "Block800") %>%
#  filter(rev_comp == 1)
#gei_dat %>% filter(block == "Block800") %>%
#  filter(rev_comp == 0)
#
#gei_dat %>% select(block, inversion, rev_comp, norm_exp) %>%
#  group_by(block, rev_comp)
#
#gei_dat %>%
#  select(filter(rev_comp == 1),block, norm_exp, inversion)
#
#
#gei_dat %>%
#  filter(block == "Block800" && block == "Block799") %>%
#  group_by(block) %>%
#  do(tidy(t.test(norm_exp ~ rev_comp, data = .)))
#
##inverted rows we want
#gei_dat[gei_dat$block %in% gei_dat$block[gei_dat$rev_comp == 1],]
#gei_dat %>%
#  filter(rev_comp == 1)
#