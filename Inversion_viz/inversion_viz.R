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
                  text = element_text(size = 24),
                  axis.text=element_text(),
                  axis.title = element_text(),
                  #plot margins
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                  #for second legend on y-axis
                  axis.text.y.right = element_text(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  #                  legend.title = element_blank(),
                  legend.text = element_text(),
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
file_name <- "inversion_block_info_all.txt"
#file_name <- "Inversion_viz/inversion_block_info_all_ATCC_revcomp.txt"
block_inf <- read.table(file_name,sep = "\t", header = FALSE)
colnames(block_inf) <- c("block","strain","start","end","rev_comp","inversion")
summary(block_inf)
#block_inf$strain <- factor(block_inf$strain, c("U00096", "NC_010473", "NZ_CP009273", "NZ_CP009072"))
block_inf$strain <- factor(block_inf$strain, c("U00096000", "NC_010473", "NZ_CP009273", "CP009072"))
print("make column of midpoint of each block")
summary(block_inf)
block_inf <- within(block_inf, midpoint <- (end + start) /2)
colnames(block_inf)[colnames(block_inf) == "midpoint"] <- "midpoint"
head(block_inf)
#get the genome pos into Mbp
block_inf$midpoint <- as.numeric(block_inf$midpoint)
#block_inf$midpoint <- as.numeric(block_inf$midpoint / 1000000)
head(block_inf)

tmpd <- block_inf %>% filter(strain == "NZ_CP009072")
sum(tmpd$rev_comp)
length(tmpd$rev_comp)


##test rev comp the non-revcomp data
#manual_rev <- block_inf %>% filter(strain == "CP009072") %>%
#              mutate(midpoint, as.numeric(abs(midpoint - 5130767)))
#colnames(manual_rev)[8] <- "rev_midpoint"
#manual_rev <- manual_rev %>% select(!midpoint)
#head(manual_rev)
#colnames(manual_rev)[7] <- "midpoint"
#sum(manual_rev$rev_comp)
#length(manual_rev$rev_comp)
#
#other_t <- block_inf %>% filter(strain != "CP009072")
#head(other_t)
#block_inf <- rbind(manual_rev,other_t)

#making new column with inversion info for each block
inver <- unique(block_inf %>% select(block, inversion))
inversion <- inver$inversion
head(inver)
bi_dat <- block_inf %>% select(block,strain,midpoint)
bi_dat <-  spread(bi_dat, strain, midpoint)
#adding inversion col to data for colouring (hopefully)
bi_dat <- cbind(bi_dat, inversion)
bi_dat$inversion <- as.factor(bi_dat$inversion)
head(bi_dat)
#bi_dat <- bi_dat[c(1:10,1140:1150),]
bi_dat <- bi_dat[order(bi_dat[,5],decreasing = FALSE),]
bi_dat
summary(bi_dat)


print("test parallel sets plot")
ps_dat <- bi_dat %>%
  gather_set_data(2:5)
#ps_dat$x <- factor(ps_dat$x, rev(c("U00096", "NC_010473", "NZ_CP009273", "NZ_CP009072")))
#levels(ps_dat$x) <- list("K12 MG1655"="U00096", "K12 DH10B"="NC_010473", "BW25113"="NZ_CP009273", "ATCC"="NZ_CP009072")
ps_dat$x <- factor(ps_dat$x, rev(c("U00096000", "NC_010473", "NZ_CP009273", "CP009072")))
levels(ps_dat$x) <- list("K12 MG1655"="U00096000", "K12 DH10B"="NC_010473", "BW25113"="NZ_CP009273", "ATCC"="CP009072")
ps_dat$x <- factor(ps_dat$x, rev(c("K12 MG1655", "K12 DH10B", "BW25113", "ATCC")))
head(ps_dat)
summary(ps_dat)

ecolT <- substitute(italic(ecoli)~strain, list(ecoli="E.coli",strain="Strain"))

ps <- (ggplot(data =ps_dat, aes(x, id = id, split = y, value = 1))
#  + geom_parallel_sets(aes(fill = U00096000))
  + geom_parallel_sets(aes(fill = inversion))
  + scale_fill_manual(values = c("#2E294E","#BEBEBE"))
#  + geom_parallel_sets(aes(fill = U00096 ))
  + xlab(ecolT) 
  + ylab("Genomic Position")
  + coord_flip()
  + scale_x_discrete(expand = c(0,0))
)
ps
#pdf("Inversion_viz/parallel_sets_maual_revcomp.pdf", width = 15, height = 7)
pdf("Inversion_viz/parallel_sets_all.pdf", width = 15, height = 7)
#pdf("Inversion_viz/parallel_sets_all_ATCC_revcomp.pdf", width = 15, height = 7)
ps
dev.off()


#
#
##checking to see what the difference is btwn rev comp and non rev comp blocks
#rev_ATCC1 <- block_inf
#rev_ATCC1 <- rev_ATCC1[order(rev_ATCC1[,1],decreasing = FALSE),]
#
#
#tmp_rev <- rev_ATCC1 %>% filter(strain == "NZ_CP009273" | strain == "NC_010473")
#levels(tmp_rev$block) <- sub("^Block", "", levels(tmp_rev$block))
#levels(tmp_rev$block) <- sub("\\.mafft", "", levels(tmp_rev$block))
#tmp_rev$block <- as.numeric(tmp_rev$block)
#head(tmp_rev)
#rev_ATCC <- rev_ATCC1 %>% filter(strain == "NZ_CP009273")
#rev_ATCC <- rev_ATCC[order(rev_ATCC[,3],decreasing = FALSE),]
#rev_ATCC
#
#print("read in block info file")
#file_name <- "inversion_block_info_all.txt"
#block_inf <- read.table(file_name,sep = "\t", header = FALSE)
#colnames(block_inf) <- c("block","strain","start","end","rev_comp","inversion")
#summary(block_inf)
#block_inf$strain <- factor(block_inf$strain, c("U00096000", "NC_010473", "NZ_CP009273", "CP009072"))
#print("make column of midpoint of each block")
#summary(block_inf)
#block_inf <- within(block_inf, midpoint <- (end + start) /2)
#colnames(block_inf)[colnames(block_inf) == "midpoint"] <- "midpoint"
#head(block_inf)
##get the genome pos into Mbp
#block_inf$midpoint <- as.numeric(block_inf$midpoint)
#block_inf$midpoint <- as.numeric(block_inf$midpoint / 1000000)
#head(block_inf)
#nonrev1 <- block_inf
#nonrev1 <- nonrev1[order(nonrev1[,1],decreasing = FALSE),]
#
#
#manual_rev <- nonrev1 %>% filter(strain == "CP009072") %>%
#              mutate(start, abs(start - 5130767)) %>%
#              mutate(end, abs(end - 5130767))
#colnames(manual_rev)[8] <- "rev_start"
#colnames(manual_rev)[9] <- "rev_end"
#head(manual_rev)
#
#nonrev <- nonrev1 %>% filter(strain == "NZ_CP009273")
#nonrev <- nonrev[order(nonrev[,3],decreasing = FALSE),]
#
#
#sum(rev_ATCC$rev_comp)
#sum(nonrev$rev_comp)
#rev_ATCC[1:15,]
#nonrev[1:15,]
#nonrev1 %>% filter(block == "Block7.mafft")
#nonrev1 %>% filter(block == "Block58.mafft")
#rev_ATCC1 %>% filter(block == "Block773.mafft")
#