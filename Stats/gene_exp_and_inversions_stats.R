#test script for inversions and gene expression analysis
library(tidyverse)
library(dplyr)
library(tidyr)

#read in the data
gei_dat <- read.csv("Sample_final_df.csv", header = TRUE)

#Wilcox signed-rank test
# to see if there is a difference btwn inversions and non-inverted segments
gei_dat %>% select(norm_exp, inversion) %>% 
  group_by(inversion) %>%
  summarise(norm_exp = list(norm_exp)) %>% 
  spread(inversion, norm_exp) %>%
  mutate(p_value = t.test(unlist(0), unlist(1))$p.value,
         t_value = t.test(unlist(0), unlist(1))$statistic)
  
  
  wilcox.test(gei_dat$norm_exp, gei_dat$norm_exp, paired=TRUE)
