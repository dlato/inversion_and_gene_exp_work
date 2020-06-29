#test script for inversions and gene expression analysis
library(tidyverse)
library(dplyr)
library(tidyr)
library(broom)

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
t.test(gei_dat$norm_exp[gei_dat$block == "Block800",], gei_dat$norm_exp[gei_dat$block == "Block800",])

test_inver <- gei_dat %>% filter(inversion == 1)
test_revcomp <- gei_dat %>% filter(rev_comp == 1)

gei_dat %>% filter(block == "Block800") %>%
  filter(rev_comp == 1)
gei_dat %>% filter(block == "Block800") %>%
  filter(rev_comp == 0)

gei_dat %>% select(block, inversion, rev_comp, norm_exp) %>%
  group_by(block, rev_comp)

gei_dat %>%
  select(filter(rev_comp == 1),block, norm_exp, inversion)


gei_dat %>%
  filter(block == "Block800" && block == "Block799") %>%
  group_by(block) %>%
  do(tidy(t.test(norm_exp ~ rev_comp, data = .)))

#inverted rows we want
gei_dat[gei_dat$block %in% gei_dat$block[gei_dat$rev_comp == 1],]
gei_dat %>%
  filter(rev_comp == 1)
