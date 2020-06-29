#test script for inversions and gene expression analysis
library(tidyverse)
library(dplyr)
library(tidyr)

#read in the data
gei_dat <- read.csv("Sample_final_df.csv", header = TRUE)

#Wilcox signed-rank test
# to see if there is a difference btwn inversions and non-inverted segments
#null hyp = there is no difference in gene expression btwn inverted and non-inverted segments
wilcox.test(gei_dat$norm_exp[gei_dat$inversion == 1], gei_dat$norm_exp[gei_dat$inversion == 0])
t.test(gei_dat$norm_exp[gei_dat$inversion == 1], gei_dat$norm_exp[gei_dat$inversion == 0])
