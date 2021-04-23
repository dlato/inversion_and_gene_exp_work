#HNS and Fis permutation test

library(regioneR)
library(dplyr)

#######
# HNS
#######

#read in data
hns_d <- read.table("./Stats/hns_binding_all_datasets.csv")
#make proper bed file format
chrom <- rep("chr1", length(nrow(hns_d)))
hns_d <- cbind(chrom,hns_d)
head(hns_d)
summary(hns_d)
levels(hns_d$data_set)
#split into multiple dataframes
hns_ALL <- hns_d[,-2]

hns_g <- hns_d %>%
         filter(data_set == "G") %>%
         select(c("chrom","start","end"))
hns_h1 <- hns_d %>%
         filter(data_set == "H1") %>%
         select(c("chrom","start","end"))
hns_h2 <- hns_d %>%
         filter(data_set == "H2") %>%
         select(c("chrom","start","end"))
hns_h3 <- hns_d %>%
         filter(data_set == "H3") %>%
         select(c("chrom","start","end"))
hns_h4 <- hns_d %>%
         filter(data_set == "H4") %>%
         select(c("chrom","start","end"))
hns_h5 <- hns_d %>%
         filter(data_set == "H5") %>%
         select(c("chrom","start","end"))
hns_h6 <- hns_d %>%
         filter(data_set == "H6") %>%
         select(c("chrom","start","end"))
hns_l <- hns_d %>%
         filter(data_set == "L") %>%
         select(c("chrom","start","end"))
hns_o <- hns_d %>%
         filter(data_set == "O") %>%
         select(c("chrom","start","end"))
hns_u <- hns_d %>%
         filter(data_set == "U") %>%
         select(c("chrom","start","end"))
summary(hns_g)

rownames(hns_d) <- hns_d[,1]
#convert to proper object
hns_reg <- toGRanges(hns_ALL)
