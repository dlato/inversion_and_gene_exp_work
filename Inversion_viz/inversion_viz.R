#script for creating inversion visualization
#library(tidyverse)
library(dplyr)
library(tidyr)
library("stringr")
library(scales) # needed for oob parameter
library(viridis)
#library(tximport)
library(readr)
library(reshape2)
library(data.table)
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("GenomicRanges")
library(ggforce)
#set theme
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            #change size of facet header text
            theme(strip.text = element_text(size =10.49)) +
            theme(plot.title = element_text(hjust = 0.5, size = 18),
                  panel.background = element_rect(fill = "white", colour = NA),
                  panel.grid.major = element_line(colour = "grey90", size = 0.2),
                  panel.grid.minor = element_line(colour = "grey98", size = 0.5),
                  panel.spacing = unit(0.25, "lines"),
                  axis.text=element_text(size=18),
                  axis.title = element_text(size = 18),
                  #plot margins
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                  #for second legend on y-axis
                  axis.text.y.right = element_text(size=18),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  #                  legend.title = element_blank(),
                  legend.text = element_text(size = 18),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facest
                  panel.spacing.x=unit(0, "lines"),
                  #                  legend.key = element_blank(),
                  #                  legend.background=element_blank(),
                  #                  legend.position="none")
                  legend.position="none")
)



print("###############################################################################")
print("INVERSION VISUALIZATION PARALLEL SETS")
print("###############################################################################")
print("read in block info file")
#file_name <- "inversion_block_info_all.txt"
file_name <- "Inversion_viz/inversion_block_info_all_ATCC_revcomp.txt"
block_inf <- read.table(file_name,sep = "\t", header = FALSE)
colnames(block_inf) <- c("block","strain","start","end","rev_comp","inversion")
print("make column of midpoint of each block")
summary(block_inf)
block_inf <- within(block_inf, midpoint <- (end + start) /2)
colnames(block_inf)[colnames(block_inf) == "midpoint"] <- "midpoint"
head(block_inf)
#get the genome pos into Mbp
block_inf$midpoint <- as.numeric(block_inf$midpoint)
block_inf$midpoint <- as.numeric(block_inf$midpoint / 1000000)
head(block_inf)
bi_dat <- block_inf %>% select(block,strain,midpoint)
bi_dat <-  spread(bi_dat, strain, midpoint)
head(bi_dat)
#bi_dat <- bi_dat[1:10,]
bi_dat <- bi_dat[order(bi_dat[,5],decreasing = FALSE),]
bi_dat
summary(bi_dat)

print("test parallel sets plot")
ps_dat <- bi_dat %>%
  gather_set_data(2:5)
head(ps_dat)
summary(ps_dat)

ps <- (ggplot(data =ps_dat, aes(x, id = id, split = y, value = 1))
  #+ geom_parallel_sets(aes(fill = U00096000 ))
  + geom_parallel_sets(aes(fill = U00096 ))
  + xlab("Strain") 
  + ylab("Genomic Position")
  + coord_flip()
  + scale_x_discrete(expand = c(0,0))
)
#pdf("Inversion_viz/parallel_sets_first_10_blocks.pdf", width = 15, height = 7)
#pdf("Inversion_viz/parallel_sets_all.pdf", width = 15, height = 7)
pdf("Inversion_viz/parallel_sets_all_ATCC_revcomp.pdf", width = 15, height = 7)
ps
dev.off()
