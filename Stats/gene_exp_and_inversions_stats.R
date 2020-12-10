#test script for inversions and gene expression analysis
#library(tidyverse)
library(dplyr)
library(plyr)
library(tidyr)
library("stringr")
#library(broom)
#library(rstatix)
#library(ggpubr)
library(DESeq2)
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
                  #plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                  #for second legend on y-axis
                  axis.text.y.right = element_text(size=18),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 18),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facest
                  panel.spacing.x=unit(0, "lines"),
#                  legend.key = element_blank(),
                  legend.background=element_blank(),
                  #legend background
                  legend.key = element_rect(fill = NA),
#                  legend.position="none")
                  legend.position="top")
)


#read in the data
#gbk_start and gbk_end are the starts and ends for each gene in the block
# start and end are the starts and ends for the WHOLE block
#gei_dat <- read.csv("../Sample_final_df.csv", header = TRUE)
#getting command line arguments
args <- commandArgs(TRUE)
file_name <- as.character(args[1])
gei_dat <- read.csv(file_name, header = TRUE)

print("make column of length of each block")
gei_dat <- within(gei_dat, block_len <- end - start)
colnames(gei_dat)[colnames(gei_dat) == "block_len"] <- "block_len"
print("make column of midpoint of each block")
gei_dat <- within(gei_dat, midpoint <- (end + start) /2)
colnames(gei_dat)[colnames(gei_dat) == "midpoint"] <- "midpoint"
head(gei_dat)
bac_name <- as.character(args[5])
replicon <- as.character(args[6])

print("summary of block lengths:")
print("non-inversion")
summary(gei_dat$block_len[which(gei_dat$inversion == 0)])
print("inversion")
summary(gei_dat$block_len[which(gei_dat$inversion == 1)])

print("summary of gei_dat")
summary(gei_dat)
print("K12MG")
gei_dat %>% filter(strain == "K12MG") %>% filter(rev_comp == 1)
print("K12DH")
gei_dat %>% filter(strain == "K12DH") %>% filter(rev_comp == 1)
print("-------------------------------------------------------")
gei_dat %>% filter(block == "Block800")
print("-------------------------------------------------------")
gei_dat %>% filter(block == "Block799")
print("-------------------------------------------------------")
print("BW25113")
gei_dat %>% filter(strain == "BW25113") %>% filter(rev_comp == 1)
print("summary of gei_dat k12 only")
tmp_k12 <- gei_dat %>% filter(strain == "K12MG")
inver_tmp_k12 <- tmp_k12
summary(tmp_k12)


options("scipen"=100, "digits"=10)

print("################################################################################")
print("#ORIGIN SCALING AND BIDIRECTIONALITY MIDPOINT                                            ")
print("################################################################################")
#first scaling things to the origin (if necessary)
max_pos <- as.numeric(args[2])
print("max_pos")
max_pos
oriC_pos <- as.numeric(args[3])
print("oriC")
oriC_pos
terminus <- as.numeric(args[4])
print("ter")
terminus
new_pos <- gei_dat$midpoint
tmp_pos <- gei_dat$midpoint
print("MIN POS")
min(gei_dat$midpoint)

if (bac_name == "E.coli" | replicon == "pSymA") {
  to_shift_ter <- max_pos - oriC_pos
  shifted_ter <-terminus + to_shift_ter
  terminus <- shifted_ter
}
print("shifted ter")
terminus

if (replicon == "pSymB") {
  shifted_ter <- terminus - oriC_pos
  terminus <- shifted_ter
}

if (bac_name == "E.coli" | replicon == "pSymA" | replicon == "pSymB")
{
  for(i in 1:length(tmp_pos)) {
    if (tmp_pos[i] >= oriC_pos) {
      new_pos[i] <- tmp_pos[i] - oriC_pos
    } else {
      tmp_end <- max_pos - oriC_pos
      new_pos[i] <- tmp_pos[i] + tmp_end
    }
  }
  tmp_pos <- new_pos
}


#now accounting for the bidirectionality. if things are between the start pos and
#the terminus then they will stay as the same position. If not, then they will be
#changed to a new position starting at 1 and going to the terminus
new_pos2 <- tmp_pos
#also have to account for bidirectional replicaion in the strand, in
#the left replichore a complemented gene (1) is actually on the leading
#strand. so for the left replichore, all 0 -> 1, and 1 -> 0
new_strand <- gei_dat$gbk_strand
if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name ==
"S.meliloti") {
  for(i in 1:length(tmp_pos)) {
    #left replichore
    if (tmp_pos[i] > terminus) {
      new_pos2[i] <- max_pos - tmp_pos[i]
      # making sure the strand column accounts for bidirectional rep
      if (gei_dat$gbk_strand[i] == 0) {
        new_strand[i] <- 1
      } else {
        new_strand[i] <- 0
      }
    } else {
    }
  }
  tmp_pos <- new_pos2

  print("max tmp_pos")
  max(tmp_pos)
}


if (bac_name == "Streptomyces") {
  for(i in 1:length(tmp_pos)) {
    # right replichore
    if (tmp_pos[i] >= oriC_pos) {
      new_pos[i] <- tmp_pos[i] - oriC_pos
    }#if btwn origin and end of genome
    # left replichore
    if (tmp_pos[i] <= oriC_pos) {
      new_pos[i] <- -1 * (oriC_pos - tmp_pos[i])
      # making sure the strand column accounts for bidirectional rep
      if (gei_dat$gbk_strand[i] == 0) {
        new_strand[i] <- 1
      } else {
        new_strand[i] <- 0
      }
    }#if btwn origin and beginning of genome
    if (tmp_pos[i] == oriC_pos) {
      new_pos[i] <- 0
    }#if equal to origin
  }
  tmp_pos <- new_pos
}


gei_dat$gbk_strand <- new_strand

gei_dat$midpoint <- tmp_pos
#gei_dat <- as.data.frame(cbind(gei_dat$block, gei_dat$gene, gei_dat$sec, tmp_pos, gei_dat$dS, gei_dat$dN, gei_dat$omega, gei_dat$sec_len))
#colnames(gei_dat) <- c("block","end","gbk_end","gbk_gene_id","gbk_locus_tag","gbk_midpoint","gbk_old_locus_tag","gbk_start","gbk_strand","gene_id","inversion","locus_tag","norm_exp","rev_comp","start","strain","taxa","gene_name")
head(gei_dat)
max(gei_dat$midpoint)
min(gei_dat$midpoint)

print("################################################################################")
print("#ORIGIN SCALING AND BIDIRECTIONALITY gbk_midpoint                                            ")
print("################################################################################")
#first scaling things to the origin (if necessary)
max_pos <- as.numeric(args[2])
print("max_pos")
max_pos
oriC_pos <- as.numeric(args[3])
print("oriC")
oriC_pos
terminus <- as.numeric(args[4])
print("ter")
terminus
new_pos <- gei_dat$gbk_midpoint
tmp_pos <- gei_dat$gbk_midpoint
print("MIN POS")
min(gei_dat$gbk_midpoint)

if (bac_name == "E.coli" | replicon == "pSymA") {
  to_shift_ter <- max_pos - oriC_pos
  shifted_ter <-terminus + to_shift_ter
  terminus <- shifted_ter
}
print("shifted ter")
terminus

if (replicon == "pSymB") {
  shifted_ter <- terminus - oriC_pos
  terminus <- shifted_ter
}

if (bac_name == "E.coli" | replicon == "pSymA" | replicon == "pSymB")
{
  for(i in 1:length(tmp_pos)) {
    if (tmp_pos[i] >= oriC_pos) {
      new_pos[i] <- tmp_pos[i] - oriC_pos
    } else {
      tmp_end <- max_pos - oriC_pos
      new_pos[i] <- tmp_pos[i] + tmp_end
    }
  }
  tmp_pos <- new_pos
}


print("now for bidir")
#now accounting for the bidirectionality. if things are between the start pos and
#the terminus then they will stay as the same position. If not, then they will be
#changed to a new position starting at 1 and going to the terminus
new_pos2 <- tmp_pos
#also have to account for bidirectional replicaion in the strand, in
#the left replichore a complemented gene (1) is actually on the leading
#strand. so for the left replichore, all 0 -> 1, and 1 -> 0
new_strand <- gei_dat$gbk_strand
if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name ==
"S.meliloti") {
  for(i in 1:length(tmp_pos)) {
    #left replichore
    if (tmp_pos[i] > terminus) {
      new_pos2[i] <- max_pos - tmp_pos[i]
      # making sure the strand column accounts for bidirectional rep
      if (gei_dat$gbk_strand[i] == 0) {
        new_strand[i] <- 1
      } else {
        new_strand[i] <- 0
      }
    } else {
    }
  }
  tmp_pos <- new_pos2

  print("max tmp_pos")
  max(tmp_pos)
}


if (bac_name == "Streptomyces") {
  for(i in 1:length(tmp_pos)) {
    # right replichore
    if (tmp_pos[i] >= oriC_pos) {
      new_pos[i] <- tmp_pos[i] - oriC_pos
    }#if btwn origin and end of genome
    # left replichore
    if (tmp_pos[i] <= oriC_pos) {
      new_pos[i] <- -1 * (oriC_pos - tmp_pos[i])
      # making sure the strand column accounts for bidirectional rep
      if (gei_dat$gbk_strand[i] == 0) {
        new_strand[i] <- 1
      } else {
        new_strand[i] <- 0
      }
    }#if btwn origin and beginning of genome
    if (tmp_pos[i] == oriC_pos) {
      new_pos[i] <- 0
    }#if equal to origin
  }
  tmp_pos <- new_pos
}


gei_dat$gbk_strand <- new_strand

gei_dat$gbk_midpoint <- tmp_pos
#gei_dat <- as.data.frame(cbind(gei_dat$block, gei_dat$gene, gei_dat$sec, tmp_pos, gei_dat$dS, gei_dat$dN, gei_dat$omega, gei_dat$sec_len))
#colnames(gei_dat) <- c("block","end","gbk_end","gbk_gene_id","gbk_locus_tag","gbk_midpoint","gbk_old_locus_tag","gbk_start","gbk_strand","gene_id","inversion","locus_tag","norm_exp","rev_comp","start","strain","taxa","gene_name")
head(gei_dat)
summary(gei_dat)
max(gei_dat$gbk_midpoint)
min(gei_dat$gbk_midpoint)

print("############    ")
print("# Total number of inversions identified")
print("############    ")
atcc_inver <- gei_dat %>% filter(inversion == 1)
length(unique(atcc_inver$block))
atcc_inver <- gei_dat %>% filter(strain == "ATCC") %>% filter(inversion == 1)
atc_b <- unique(atcc_inver$block)
atcc_inver <- gei_dat %>% filter(strain == "K12DH") %>% filter(inversion == 1)
dh_b <- unique(atcc_inver$block)
print("############    ")
print("# Total number of inversions identified in ATCC")
print("############    ")
length(setdiff(atc_b,dh_b))
print("############    ")
print("# Total number of inversions identified in DH")
print("############    ")
length(setdiff(dh_b,atc_b))
print("intersect")
length(intersect(atc_b, dh_b))

print("############    ")
print("# ALL inversions=1")
print("############    ")
print("#Wilcox signed-rank test                                                                     ")
print("# to see if there is a difference btwn inversions and non-inverted segments                  ")
print("#null hyp = there is no difference in gene expression btwn inverted and non-inverted segments")
wilcox.test(gei_dat$norm_exp[gei_dat$inversion == 1], gei_dat$norm_exp[gei_dat$inversion == 0])
t.test(gei_dat$norm_exp[gei_dat$inversion == 1], gei_dat$norm_exp[gei_dat$inversion == 0])
warnings()
print("############    ")
print("# ALL rev_comp=1")
print("############    ")
wilcox.test(gei_dat$norm_exp[gei_dat$rev_comp == 1], gei_dat$norm_exp[gei_dat$rev_comp == 0])
t.test(gei_dat$norm_exp[gei_dat$rev_comp == 1], gei_dat$norm_exp[gei_dat$rev_comp == 0])
warnings()

print("############         ")
print("# by block inversions")
print("#Wilcox signed-rank test                                                                     ")
print("############         ")
print("blocks with at least one rev comp taxa")
tmp_df <- gei_dat[which(gei_dat$rev_comp == 1),]
print(unique(tmp_df$block))
tmp_df <- gei_dat[which(gei_dat$block == "Block800"),]
print(tmp_df)
print("################################################")
uniq_block_w_pvalue <- c()
uniq_block_t_pvalue <- c()
uniq_block_t_stat <- c()
uniq_block_t_confinf1 <- c()
uniq_block_t_confinf2 <- c()
uniq_avg_exp_invert <- c()
uniq_avg_exp_noninvert <- c()
uniq_avg_len_invert <- c()
uniq_avg_len_noninvert <- c()
uniq_block_tax_num <- c()
for (i in unique(gei_dat$block)) {
  tmp_df <- gei_dat[which(gei_dat$block == i),]
print("--------------------")
print(tmp_df)
  if (length(unique(tmp_df$rev_comp)) > 1) {
    #mean experssion for each group
    i_avg <- mean(tmp_df$norm_exp[tmp_df$rev_comp == 1])
    ni_avg <- mean(tmp_df$norm_exp[tmp_df$rev_comp == 0])
    #mean block length for each group
    i_len <- mean(tmp_df$block_len[tmp_df$rev_comp == 1])
    ni_len <- mean(tmp_df$block_len[tmp_df$rev_comp == 0])
    tnum <- length(unique(tmp_df$strain))
print("means")
print(i_len)
print(ni_len)
print("tnum")
print(tnum)
    uniq_avg_exp_invert <- c(uniq_avg_exp_invert,i_avg)
    uniq_avg_exp_noninvert <- c(uniq_avg_exp_noninvert,ni_avg)
    uniq_avg_len_invert <- c(uniq_avg_len_invert,i_len)
    uniq_avg_len_noninvert <- c(uniq_avg_len_noninvert,ni_len)
    uniq_block_tax_num <- c(uniq_block_tax_num,tnum)
    #wilcox test
    wilcox_test <- wilcox.test(tmp_df$norm_exp[tmp_df$rev_comp == 1], tmp_df$norm_exp[tmp_df$rev_comp == 0])
    uniq_block_w_pvalue <- c(uniq_block_w_pvalue,wilcox_test$p.value)
    #t test
    if (length(which(tmp_df$rev_comp == 1)) > 1 && length(which(tmp_df$rev_comp == 0)) > 1){
      t_test <- t.test(tmp_df$norm_exp[tmp_df$rev_comp == 1], tmp_df$norm_exp[tmp_df$rev_comp == 0])
      uniq_block_t_pvalue <- c(uniq_block_t_pvalue,t_test$p.value)
      uniq_block_t_stat <- c(uniq_block_t_stat,t_test$statistic)
      uniq_block_t_confinf1 <- c(uniq_block_t_confinf1,t_test$conf.int[1])
      uniq_block_t_confinf2 <- c(uniq_block_t_confinf2,t_test$conf.int[2])
    } else {
      uniq_block_t_pvalue <- c(uniq_block_t_pvalue,"NA")
      uniq_block_t_stat <- c(uniq_block_t_stat,"NA")
      uniq_block_t_confinf1 <- c(uniq_block_t_confinf1,"NA")
      uniq_block_t_confinf2 <- c(uniq_block_t_confinf2,"NA")
    }
  } else {
    tnum <- length(unique(tmp_df$strain))
    uniq_block_tax_num <- c(uniq_block_tax_num,tnum)
    #not able to do tests
    uniq_block_w_pvalue <- c(uniq_block_w_pvalue,"NA")
    uniq_block_t_pvalue <- c(uniq_block_t_pvalue,"NA")
    uniq_block_t_stat <- c(uniq_block_t_stat,"NA")
    uniq_block_t_confinf1 <- c(uniq_block_t_confinf1,"NA")
    uniq_block_t_confinf2 <- c(uniq_block_t_confinf2,"NA")
    uniq_avg_exp_invert <- c(uniq_avg_exp_invert,"NA")
    uniq_avg_exp_noninvert <- c(uniq_avg_exp_noninvert,"NA")
    uniq_avg_len_invert <- c(uniq_avg_len_invert,"NA")
    uniq_avg_len_noninvert <- c(uniq_avg_len_noninvert,"NA")
  }
}
warnings()

#add test results to df
count = 1
block_tax_num <- vector(mode='numeric', length = length(gei_dat$block))
block_w_pvalue <- vector(mode='numeric', length = length(gei_dat$block))
block_t_pvalue <- vector(mode='numeric', length = length(gei_dat$block))
block_t_stat <- vector(mode='numeric', length = length(gei_dat$block))
block_t_coninf1 <- vector(mode='numeric', length = length(gei_dat$block))
block_t_coninf2 <- vector(mode='numeric', length = length(gei_dat$block))
block_avg_exp_invert <- vector(mode='numeric', length = length(gei_dat$block))
block_avg_exp_noninvert <- vector(mode='numeric', length = length(gei_dat$block))
block_avg_len_invert <- vector(mode='numeric', length = length(gei_dat$block))
block_avg_len_noninvert <- vector(mode='numeric', length = length(gei_dat$block))
for (b in unique(gei_dat$block)) {
  bloc_loc <- which(gei_dat$block == b)
  block_tax_num[bloc_loc] <- uniq_block_tax_num[count] 
  block_w_pvalue[bloc_loc] <- uniq_block_w_pvalue[count] 
  block_t_pvalue[bloc_loc] <- uniq_block_t_pvalue[count] 
  block_t_stat[bloc_loc] <- uniq_block_t_stat[count] 
  block_t_coninf1[bloc_loc] <- uniq_block_t_confinf1[count] 
  block_t_coninf2[bloc_loc] <- uniq_block_t_confinf2[count] 
  block_avg_exp_invert[bloc_loc] <- uniq_avg_exp_invert[count] 
  block_avg_exp_noninvert[bloc_loc] <- uniq_avg_exp_noninvert[count] 
  block_avg_len_invert[bloc_loc] <- uniq_avg_len_invert[count] 
  block_avg_len_noninvert[bloc_loc] <- uniq_avg_len_noninvert[count] 
  count <- count + 1
}
gei_dat['block_tax_num'] <- block_tax_num
gei_dat['block_w_pvalue'] <- block_w_pvalue
gei_dat['block_t_pvalue'] <- block_t_pvalue
gei_dat['block_t_stat'] <- block_t_stat
gei_dat['block_t_confinf1'] <- block_t_coninf1
gei_dat['block_t_confinf2'] <- block_t_coninf2
gei_dat['block_t_confinf2'] <- block_t_coninf2
gei_dat['block_avg_exp_invert'] <- block_avg_exp_invert
gei_dat['block_avg_exp_noninvert'] <- block_avg_exp_noninvert
gei_dat['block_avg_len_invert'] <- block_avg_len_invert
gei_dat['block_avg_len_noninvert'] <- block_avg_len_noninvert

head(gei_dat)
print("#####################")
print("ONLY USING BLOCKS WITH ALL 4 TAXA")
print("#####################")
gei_dat <- gei_dat %>% filter(block_tax_num == 4)
print("SAVED DATA TO FILE")
write.table(gei_dat, 'inversions_gene_exp_wtest_data.csv', sep = "\t")
bp_dat <- gei_dat
inver_dat_bidir <- gei_dat

print("make df with just block info")
block_df <- subset(gei_dat,select = c("block","start","end","midpoint","gbk_start","gbk_end","gbk_midpoint","block_w_pvalue","block_t_pvalue","block_t_stat","block_t_confinf1","block_t_confinf2", "block_avg_exp_invert", "block_avg_exp_noninvert","block_avg_len_invert", "block_avg_len_noninvert","inversion","rev_comp","strain"))
#block_df_uniq <- unique(block_df)
block_df_uniq <- block_df

#flip so ggplot can use
block_df_w <- melt(block_df_uniq,
        # ID variables - all the variables to keep but not split apart
        # on
    id.vars=c("block", "start", "end","midpoint","gbk_start","gbk_end","gbk_midpoint","block_t_stat","block_t_confinf1","block_t_confinf2","block_avg_exp_invert", "block_avg_exp_noninvert","block_avg_len_invert", "block_avg_len_noninvert","inversion","rev_comp", "strain"),
        # The source columns
    measure.vars=c("block_w_pvalue", "block_t_pvalue"),
        # Name of the destination column that will identify the
        # original
        # column that the measurement came from
    variable.name="test",
    value.name="pvalue"
)
block_df_w$pvalue <- as.numeric(block_df_w$pvalue)
summary(block_df_w)
print("non zero pvals") 
complete_block_df <- block_df_w[which(block_df_w$pvalue != "NA"),]
summary(complete_block_df)
complete_block_df <- complete_block_df %>% filter(test == "block_w_pvalue")
head(complete_block_df)
print("SAVED complete DATA TO FILE")
write.table(complete_block_df, 'inversions_gene_exp_complete_data.csv', sep = "\t")
n_occur <- data.frame(table(complete_block_df$gbk_midpoint))
n_occur[n_occur$Freq > 2,]
complete_block_df[complete_block_df$gbk_midpoint %in% n_occur$Var1[n_occur$Freq > 2],]
#block_df_w$block_t_stat <- as.numeric(block_df_w$block_t_stat)
#block_df_w$block_t_confinf <- as.numeric(block_df_w$block_t_confinf)
complete_block_df$pvalue <- as.numeric(complete_block_df$pvalue)
#str(block_df_w)
#head(complete.cases(block_df_w$pvalue))
#complete_block_df <- block_df_w[complete.cases(block_df_w), ]
print("TEST mult sig")
all_sig <- c()
for (i in unique(complete_block_df$block)) {
  tmp_df <- complete_block_df[which(complete_block_df$block == i),]
  if (all(tmp_df$pvalue <= 0.05)) {
    all_sig <- c(all_sig,"yes")
  } else {
    all_sig <- c(all_sig, "no")
  }#else

}#for
#add test results to df
count = 1
allSig <- vector(mode='numeric', length = length(complete_block_df$block))
for (b in unique(complete_block_df$block)) {
  bloc_loc <- which(complete_block_df$block == b)
  allSig[bloc_loc] <- all_sig[count] 
  count <- count + 1
}#for
complete_block_df['allSig'] <- allSig

print("number of SIG blocks with all sig")
td <- complete_block_df[which(complete_block_df$allSig == "yes"),]
length(unique(td$block))


#split_block <- split(complete_block_df, complete_block_df$block)
#split_block
#split_df <- lapply(split_block, function(x) if (x$pvalue <= 0.05)
#{print("all_sig")} else {print("not all_sig")})
#split_df
#
#complete_tibble <- complete_block_df %>%
#group_by(block) %>%
#mutate(if (pvalue<=0.05){all_sig = "yes"}
#    else {all_sig = "no"})
#complete_tibble[1]
#warnings()


print("complete HEAD")
head(complete_block_df)
cor_dat <- complete_block_df
cor_dat <- within(cor_dat, non_bidir_midpoint <- (gbk_end + gbk_start) /2)
summary(cor_dat)
#complete_block_df[which(complete_block_df$pvalue <= 0.05),]
#plot pvalues
p <- ggplot(complete_block_df, aes(x=test, y=pvalue)) + 
  geom_boxplot()
pdf("pvalue_box_plots.pdf")
p
dev.off()

print("########################")
print("INFO ON BLOCKS: diff gene exp between inversions and non-inversions")
print("########################")
print("total number of blocks")
length(unique(gei_dat$block))
print("number of inverted blocks")
tmp <- gei_dat[which(gei_dat$rev_comp ==1),]
length(unique(tmp$block))
print("percent of inverted blocks")
(length(unique(tmp$block)) / length(unique(gei_dat$block))) * 100
print("number of blocks that were tested")
length(unique(complete_block_df$block))
print("percent of blocks that were tested")
(length(unique(complete_block_df$block)) / length(unique(gei_dat$block))) * 100
print("number of SIGNIFICANT blocks")
tmp <- complete_block_df[which(complete_block_df$pvalue <= 0.05),]
length(unique(tmp$block))
print("percent of SIG tested blocks")
(length(unique(tmp$block)) / length(unique(complete_block_df$block))) * 100

print("gene exp averages in sig blocks")
tmp2 <- subset(tmp,select = c("block","pvalue","block_avg_exp_invert", "block_avg_exp_noninvert","block_avg_len_invert", "block_avg_len_noninvert","rev_comp","inversion"))
df <- tmp2 %>%
  mutate(exp_reg = ifelse(block_avg_exp_invert > block_avg_exp_noninvert, "up", "down"))
df <- df %>%
  mutate(len_reg = ifelse(block_avg_len_invert > block_avg_len_noninvert, "long", "short"))
print("####################")
print("GENE EXP UP/DOWN INFO")
print("####################")
print("number of SIG blocks with inversion gene exp > noninversions gene exp")
up <- df[which(df$exp_reg == "up"),]
length(unique(up$block))
print("percent of SIG blocks with inversions exp > noninversion exp")
(length(unique(up$block)) / length(unique(df$block))) *100
print("number of SIG blocks with inversion gene exp < noninversions gene exp")
up <- df[which(df$exp_reg == "down"),]
length(unique(up$block))
print("percent of SIG blocks with inversions exp < noninversion exp")
(length(unique(up$block)) / length(unique(df$block))) *100
print("####################")
print("BLOCK LENGTH INFO")
print("####################")
print("number of SIG blocks with inversion gene len > noninversions gene len")
up <- df[which(df$len_reg == "long"),]
length(unique(up$block))
print("percent of SIG blocks with inversions len > noninversion len")
(length(unique(up$block)) / length(unique(df$block))) *100
print("number of SIG blocks with inversion gene len < noninversions gene len")
up <- df[which(df$len_reg == "short"),]
length(unique(up$block))
print("percent of SIG blocks with inversions len < noninversion len")
(length(unique(up$block)) / length(unique(df$block))) *100

#make df with sig diff in gene exp btwn inversions
blocks_new <- complete_block_df %>%
    mutate(sig = if_else(pvalue <= 0.05, 'yes', 'no'))
#complete_block_df[which(complete_block_df$pvalue <=0.05),]
#complete_block_df$sig <- which(complete_block_df$pvalue <=0.05)
blocks_new$class <- rep("BlockX",length(blocks_new$block))
sig_blocks <- blocks_new[which(blocks_new$sig == "yes"),]
print("head sig_blocks")
head(sig_blocks)
summary(sig_blocks)
reg_df <- sig_blocks %>%
  mutate(exp_reg = ifelse(block_avg_exp_invert > block_avg_exp_noninvert, "up", "down"))
head(reg_df)
pdf("scatter_plot_test.pdf")
ggplot(reg_df, aes(x=gbk_midpoint, y=strain, color=exp_reg))+
#  geom_jitter(aes(tt, val), data = df, colour = I("red"), 
#               position = position_jitter(width = 0.05)) +
#  geom_point(size = 3) +
  geom_point(size = 2, alpha=0.4)
#  geom_errorbar(aes(ymin=val-sd, ymax=val+sd), width = 0.01, size = 1)
dev.off()

print("###################")
print("df with avg inversion info, only K12")
print("###################")
full_blocks_new_df <- blocks_new
#k12_df <- reg_df %>%
k12_df <- blocks_new %>%
          filter(strain == "K12MG") %>%
        gather(class, avg_exp, block_avg_exp_invert:block_avg_exp_noninvert, factor_key = TRUE)
k12_df$avg_exp <- as.numeric(k12_df$avg_exp)
head(k12_df)
summary(k12_df) 
print("###################")
print("2 data pts per block")
print("###################")
#I think midpoint is for the block and gbk_midpoint is for each gene
k12_df <- k12_df %>% select(block, midpoint,sig,class,avg_exp)
k12_df <- unique(k12_df)
head(k12_df)
k12_df[which(k12_df$block == "Block62"),]
k12_df_non_sig <- k12_df %>% filter(sig == "no")
head(k12_df_non_sig)
k12_df_sig <- k12_df %>% filter(sig == "yes")
head(k12_df_sig)
k12MG_df <- k12_df
#adjust genomic position to be in Mbp
k12_df_sig$midpoint = k12_df_sig$midpoint / 1000000
k12_df_non_sig$midpoint = k12_df_non_sig$midpoint / 1000000
print("TEST HEAD")
head(k12_df_sig)
k12_df_sig %>% filter(block == "Block164")
override.shape <- c(16, 16, 17)
override.linetype <- c(0,1, 3)
print("GRAPH IS ONLY WITH INVERTED BLOCKS THAT WERE ABLE TO BE TESTED")
#cols <- c("block_avg_exp_noninvert" = "#3d405b","block_avg_exp_invert" = "#e07a5f")
p <- (ggplot(k12_df_sig, aes(x=midpoint, y=avg_exp, color=class,shape= class))
#  geom_jitter(aes(tt, val), data = df, colour = I("red"), 
   #non-sig pts in light grey, un-filled
   + geom_point(data = k12_df_non_sig,aes(x=midpoint,y =avg_exp,colour = "#BEBEBE",shape=class), alpha=0.5 ,show.legend = FALSE, size=2.5 )
   #sig pts
   + geom_point(data = k12_df_sig,aes(x=midpoint,y = avg_exp, color = class), size=2.5)
   + geom_smooth(data=k12_df_sig,aes(x=midpoint,y = avg_exp,color=class,linetype=class),span=0.5, method = "loess")
   + scale_linetype_manual(values=c("solid", "dotted"))
#   + scale_color_manual(values = c("#BEBEBE","#e07a5f","#3d405b"),limits = "Inverted Sequences", "Non-inverted Sequences")
#   + scale_color_manual(values = c("#BEBEBE","#e07a5f" = "block_avg_exp_invert" ,"#3d405b"= "block_avg_exp_noninvert" ),limits = "block_avg_exp_invert", "block_avg_exp_noninvert",limits =c("Inversion","No Inversion"))
   + scale_color_manual(values = c("#BEBEBE","#e07a5f","#3d405b"),labels =c("Non-significant","Inverted","Non-inverted"))
##               position = position_jitter(width = 0.05)) +
##  geom_point(size = 3) +
#   + geom_point(size = 2, alpha=0.4)
##  geom_errorbar(aes(ymin=val-sd, ymax=val+sd), width = 0.01, size = 1)
   + scale_y_continuous(trans='log10')
   + scale_x_continuous(limits = c(0, 3))
   + labs(title= "",x = "Distance from the Origin of Replication (Mbp)", y = "Average Gene Expression (CPM)") 
   #colour for legend
   + guides(colour = guide_legend(override.aes = list(shape = override.shape, linetype = override.linetype, fill=NA,unit(3,"line"))))
   + scale_shape(guide = FALSE)
   + scale_linetype(guide = FALSE)
#   + guides(color=guide_legend(override.aes=list(fill=NA)))
#   + guides(fill = FALSE)
#   + scale_colour_manual(values = cols, limits = "Inverted Sequences", "Non-inverted Sequences")
#   + scale_y_continuous(trans='log10',labels = function(x) ifelse(x ==0, "0", x),breaks=c(0.0001,0.001,0.01,0.1, 1, 10,100))
)
pdf("genome_pos_inversions_k12.pdf")
p
dev.off()

#print("##############")
#print("boxplot test data")
#print("##############")
#head(bp_dat)
#bp_dat$rev_comp <- as.factor(bp_dat$rev_comp)
#sig_blocks <- k12_df_sig$block
#bp_sig <- bp_dat %>% filter(block == sig_blocks) 
#head(bp_sig)
#tail(bp_sig)
#
#p <- (ggplot(bp_sig, aes(x=block, y=norm_exp, color=rev_comp))
#   + geom_boxplot() 
###  geom_jitter(aes(tt, val), data = df, colour = I("red"), 
##   #non-sig pts in light grey, un-filled
##   + geom_point(data = k12_df_non_sig,aes(x=midpoint,y = avg_exp, color = "#BEBEBE"), alpha=0.7, shape=1)
##   #sig pts
##   + geom_point(data = k12_df_sig,aes(x=midpoint,y = avg_exp, color = class))
###   + geom_smooth(data=k12_df_sig,aes(x=midpoint,y = avg_exp,color =class),span=0.5, method = "loess")
##   + scale_color_manual(values = c("#BEBEBE","#e07a5f","#3d405b"))
####               position = position_jitter(width = 0.05)) +
####  geom_point(size = 3) +
###   + geom_point(size = 2, alpha=0.4)
####  geom_errorbar(aes(ymin=val-sd, ymax=val+sd), width = 0.01, size = 1)
#   + scale_y_continuous(trans='log10')
##   + scale_y_continuous(trans='log10',labels = function(x) ifelse(x ==0, "0", x),breaks=c(0.0001,0.001,0.01,0.1, 1, 10,100))
#)
#pdf("genome_pos_inversions_k12_boxplot.pdf")
#p
#dev.off()

print("###################")
print("# info about sig inversion blocks")
print("###################")
head(blocks_new)
blocks_new$length <- blocks_new$end - blocks_new$start
head(blocks_new)
print("diff btwn BLOCK LENGTH of inversions that have sig gene exp or not?")
wilcox.test(blocks_new$length[blocks_new$sig == "yes"], blocks_new$length[blocks_new$sig == "no"])
print("mean block length: sig blocks")
mean(blocks_new$length[blocks_new$sig == "yes"])
print("mean block length: non-sig blocks")
mean(blocks_new$length[blocks_new$sig == "no"])
print("diff btwn POSITION of inversions that have sig gene exp or not?")
wilcox.test(blocks_new$midpoint[blocks_new$sig == "yes"], blocks_new$midpoint[blocks_new$sig == "no"])
sub_blocks_new <- subset(blocks_new, select = c("sig","length","block","midpoint"))
sub_blocks_new <- unique(sub_blocks_new)
print("unique df")
head(sub_blocks_new)
tail(sub_blocks_new)
wilcox.test(sub_blocks_new$length[sub_blocks_new$sig == "yes"], sub_blocks_new$length[sub_blocks_new$sig == "no"])


print("###################")
print("# Distance from the origin info:")
print("###################")
print("# LOGISTIC REGRESSION on rev_comp and dist (midpoint)")
print("###################")
#glm with logit link
log_reg_rev_mid = glm(rev_comp ~ midpoint, family=binomial, gei_dat, control = list(maxit=1000))
print(summary(log_reg_rev_mid))
print("###################")
print("# LOGISTIC REGRESSION on rev_comp and dist (gbk_midpoint)")
print("###################")
#glm with logit link
log_reg_rev_gbk= glm(rev_comp ~ gbk_midpoint, family=binomial, gei_dat, control = list(maxit=1000))
print(summary(log_reg_rev_gbk))
print("###################")
print("# LOGISTIC REGRESSION on inversions and dist (midpoint)")
print("###################")
#glm with logit link
log_reg_inv_mid = glm(inversion ~ midpoint, family=binomial, gei_dat, control = list(maxit=1000))
print(summary(log_reg_inv_mid))
print("###################")
print("# LOGISTIC REGRESSION on inversions and dist (gbk_midpoint)")
print("###################")
#glm with logit link
log_reg_inv_gbk= glm(inversion ~ gbk_midpoint, family=binomial, gei_dat, control = list(maxit=1000))
print(summary(log_reg_inv_gbk))
print("#################################################################")
print("SIG diff in gene exp blocks")
print("###################")
print("# LOGISTIC REGRESSION on SIG rev_comp and dist (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_rev_mid = glm(rev_comp ~ midpoint, family=binomial, sig_blocks, control = list(maxit=1000))
print(summary(sig_log_reg_rev_mid))
print("###################")
print("# LOGISTIC REGRESSION on SIG rev_comp and dist (gbk_midpoint)")
print("###################")
#glm with logit link
sig_log_reg_rev_gbk= glm(rev_comp ~ gbk_midpoint, family=binomial, sig_blocks, control = list(maxit=1000))
print(summary(sig_log_reg_rev_gbk))
print("###################")
print("# LOGISTIC REGRESSION on all blocks sig and dist (gbk_midpoint)")
print("###################")
sig_df <- blocks_new %>% 
    mutate(sig = recode(sig, 
                      "yes" = "1", 
                      "no" = "0"))
sig_df$sig <- as.integer(as.character(sig_df$sig))
summary(sig_df)
#glm with logit link
sig_log_reg_rev_gbk= glm(sig ~ gbk_midpoint, family=binomial, sig_df, control = list(maxit=1000))
print(summary(sig_log_reg_rev_gbk))

print("###################")
print("# LOGISTIC REGRESSION on all blocks sig and dist (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_rev_gbk= glm(sig ~ midpoint, family=binomial, sig_df, control = list(maxit=1000))
print(summary(sig_log_reg_rev_gbk))

print("###################")
print("# LOGISTIC REGRESSION on inverted blocks sig and dist (gbk_midpoint)")
print("###################")
sig_df <- sig_df %>% filter(inversion == 1)
#glm with logit link
sig_log_reg_rev_gbk= glm(sig ~ gbk_midpoint, family=binomial, sig_df, control = list(maxit=1000))
print(summary(sig_log_reg_rev_gbk))

print("###################")
print("# LOGISTIC REGRESSION on inverted blocks sig and dist (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_rev_gbk= glm(sig ~ midpoint, family=binomial, sig_df, control = list(maxit=1000))
print(summary(sig_log_reg_rev_gbk))




print("###################")
print("# LOGISTIC REGRESSION on STRAINS")
print("# ALL BLOCKS")
print("###################")
print("K12MG")
tmp_bd <- gei_dat %>% filter(strain == "K12MG")
head(tmp_bd)
print("###################")
print("# LOGISTIC REGRESSION on K12MG rev_comp (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_inv_gbk= glm(rev_comp ~ midpoint, family=binomial, tmp_bd, control = list(maxit=1000))
print(summary(sig_log_reg_inv_gbk))
print("###################")
print("# LOGISTIC REGRESSION on K12MG inversion (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_inv_gbk= glm(inversion ~ midpoint, family=binomial, tmp_bd, control = list(maxit=1000))
print(summary(sig_log_reg_inv_gbk))

print("K12DH")
tmp_bd <- gei_dat %>% filter(strain == "K12DH")
head(tmp_bd)
print("###################")
print("# LOGISTIC REGRESSION on K12DH rev_comp (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_inv_gbk= glm(rev_comp ~ midpoint, family=binomial, tmp_bd, control = list(maxit=1000))
print(summary(sig_log_reg_inv_gbk))
print("###################")
print("# LOGISTIC REGRESSION on K12DH inversion (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_inv_gbk= glm(inversion ~ midpoint, family=binomial, tmp_bd, control = list(maxit=1000))
print(summary(sig_log_reg_inv_gbk))

print("BW25113")
tmp_bd <- gei_dat %>% filter(strain == "BW25113")
head(tmp_bd)
print("###################")
print("# LOGISTIC REGRESSION on BW25113 rev_comp (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_inv_gbk= glm(rev_comp ~ midpoint, family=binomial, tmp_bd, control = list(maxit=1000))
print(summary(sig_log_reg_inv_gbk))
print("###################")
print("# LOGISTIC REGRESSION on BW25113 inversion (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_inv_gbk= glm(inversion ~ midpoint, family=binomial, tmp_bd, control = list(maxit=1000))
print(summary(sig_log_reg_inv_gbk))


print("ATCC")
tmp_bd <- gei_dat %>% filter(strain == "ATCC")
head(tmp_bd)
print("###################")
print("# LOGISTIC REGRESSION on ATCC rev_comp (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_inv_gbk= glm(rev_comp ~ midpoint, family=binomial, tmp_bd, control = list(maxit=1000))
print(summary(sig_log_reg_inv_gbk))
print("###################")
print("# LOGISTIC REGRESSION on ATCC inversion (midpoint)")
print("###################")
#glm with logit link
sig_log_reg_inv_gbk= glm(inversion ~ midpoint, family=binomial, tmp_bd, control = list(maxit=1000))
print(summary(sig_log_reg_inv_gbk))



print("#############                 ")
print("## DESeq expression comparison")
print("#############                 ")
print("getting total combos of inversions")
#head(blocks_new)
#inver_combos <- blocks_new %>% select(block,rev_comp,strain)
#inver_combos <- unique(inver_combos)
#print("test less taxa blocks")
#print("K12MG")
#inver_combos %>% filter(strain == "K12MG") %>% filter(rev_comp == 1)
#print("K12DH")
#inver_combos %>% filter(strain == "K12MG") %>% filter(rev_comp == 1)
#print("BW25113")
#inver_combos %>% filter(strain == "BW25113") %>% filter(rev_comp == 1)
#print("ATCC")
#inver_combos %>% filter(strain == "ATCC") %>% filter(rev_comp == 1)
#taxa_order_combo <- as.vector(unlist(unique(inver_combos$strain)))
#print("order of taxa for inversion combos")
#taxa_order_combo
#rev1 <- inver_combos %>% filter(rev_comp == 1)
#rev1
#r2 <- inver_combos %>%
#         group_by(block) %>% 
#         summarise(combo = toString(rev_comp))
##         summarise(combo = list(rev_comp))
#r2 <- r2 %>%
#    ungroup %>%
#     unnest()
#r2
#
#r2[which(r2$block == "Block62"),]
#uniq_combos <- as.data.frame(unique(r2$combo))
#colnames(uniq_combos) <- "pattern"
#head(uniq_combos)
#print("inver_combos")
#uniq_combos <- uniq_combos %>% separate(pattern,
#                taxa_order_combo,", ")
#inver_combos_df <- as.data.frame(t(uniq_combos))
##colnames(inver_combos_df) <- c("A","B")
#inver_combos_df <- tibble::rownames_to_column(inver_combos_df, "strain")
#inver_combos_df

#read in the rev_comp data file:
rev_comp_inf <- read.csv("../new_dataframes/raw_rev_comp.csv", header = TRUE)
rev_comp_inf <- subset(rev_comp_inf, select = -c(1:4))
inver_combos_df <- unique(rev_comp_inf)
inver_combos_df
inver_combo <- rownames(inver_combos_df)
inver_combos_df <- cbind(inver_combo=inver_combo, inver_combos_df)
inver_combos_df <- gather(inver_combos_df, sample, treatment, ATCC_GSE94978_1:K12MG_GSE60522_3, factor_key=TRUE)
#get strain info from sample ids
inver_combos_df$strain <- gsub("_.*","",inver_combos_df$sample)
#get experiment info from sample ids
inver_combos_df$tmp_expID <- gsub("^.*_G","G",inver_combos_df$sample)
inver_combos_df$expID <- gsub("_.*","",inver_combos_df$tmp_expID)
inver_combos_df <- inver_combos_df %>% select(inver_combo,treatment,sample,strain,expID)
inver_combos_df$expID <- as.factor(inver_combos_df$expID)
head(inver_combos_df)

##set up raw exp df
##this will need to be changed once the real data comes in
##********************
print("read in raw data file of all combined experiements")
raw_file <- as.character(args[7])
raw_dat <- read.csv(raw_file, header = TRUE)
#ensure that all raw values are integers
raw_dat <- raw_dat %>% 
                mutate_at(seq(5,23,1), as.integer)
head(raw_dat)
#dealing with one K12 gene that matches to multiple ATCC genes: b3894
raw_dat <- raw_dat[-44,]
print("duplicates")
n_occur <- data.frame(table(raw_dat$MG_names))
raw_dat[raw_dat$MG_names %in% n_occur$Var1[n_occur$Freq > 1],]
#make df with JUST expression and experiment values
raw_deseq <- raw_dat %>% select(MG_names,ATCC_GSE94978_1:K12MG_GSE60522_3)
head(raw_deseq)
raw_deseqW <- raw_deseq
rownames(raw_deseqW) <- raw_deseq$MG_names
raw_deseqW <- raw_deseqW[,-1]
head(raw_deseqW)



#print(unique(raw_dat$replicates))
##make df with JUST expression and experiment values
##raw_deseq <- subset(raw_dat, select = c("Locus_tag","gene_id","inversion","replicates","raw_exp","strain","exp"))
#raw_deseq <- raw_dat %>% select(gene_id,replicates,raw_exp)
#head(raw_deseq)
##raw_deseq[c(564,2602,189,2264),]
##raw_deseqW <- spread(raw_deseq,replicates, raw_ex)
##head(raw_deseqW)
#
#print("#toy df with just K12MG experiments")
##raw_deseq <- subset(raw_dat, select = c("gene_id","inversion","replicates","raw_exp"))
#raw_deseq <- raw_dat %>% select(gene_id,replicates,raw_exp,strain)
#raw_deseq$raw_exp <- as.integer(raw_deseq$raw_exp)
#raw_deseqSubset <- raw_deseq[grep("K12MG", raw_deseq$strain), ]
#head(raw_deseqSubset)
#raw_deseqSubset <- raw_deseqSubset[c(-2602,-222,-2264,-5660,-3280,-5322),]
#raw_deseqW <- spread(raw_deseqSubset,replicates, raw_exp)
#raw_deseqW <- raw_deseqW[complete.cases(raw_deseqW),]
#head(raw_deseqW)
#
#print("get file with sample info: ALL TAXA")
#exp_info <- raw_dat %>% select(replicates,strain,exp)
#exp_info <- unique(exp_info)
#treatment_A <- c()
#for(i in 1:length(exp_info$replicate)) {
#    st <- exp_info$strain[i]
#    inver_p <- inver_combos_df$A[which(inver_combos_df$strain == st)]
#    treatment_A <- c(treatment_A,as.numeric(as.character(inver_p)))
#}
#treatment_B <- c()
#for(i in 1:length(exp_info$replicate)) {
#    st <- exp_info$strain[i]
#    inver_p <- inver_combos_df$B[which(inver_combos_df$strain == st)]
#    treatment_B <- c(treatment_B,as.numeric(as.character(inver_p)))
#}
#exp_info <- cbind(exp_info,treatment_A,treatment_B)
#exp_info
#
#print("toy sample info: only DH")
#sample_inf <- exp_info %>% filter(strain == "K12MG")
#sample_inf$treatment_A[2] <- 1
#class(sample_inf$starin)
##sample_inf$strain[4] <- "FAKE"
##sample_inf$strain <- replace(sample_inf$strain, 4, "FAKE")
#sample_inf
#
print("#############################################################")
print("DESeq on inversion combo 1")
print("#############################################################")
print("re-format sample info info format DESeq can recognize")
#sampleData <- inver_combos_df %>% filter(inver_combo == 1)
sampleData <- inver_combos_df %>% filter(inver_combo == 2179)
sampleData <- sampleData %>% select(sample, strain, expID, treatment)
colnames(sampleData) <- c("replicates","strain","expID","treatment")
rownames(sampleData) <- sampleData$replicates
sampleData$expID <- factor(sampleData$expID)
sampleData$treatment <- factor(sampleData$treatment)
sampleData$replicates <- factor(sampleData$replicates)
sampleData$strain <- factor(sampleData$strain)
sampleData <- sampleData[,-1]
ind.n <- factor(c(1,2,1,2,3,1,1,2,1,2,3,4,1,2,1,2,1,2,3))
sampleData$ind.n <- ind.n
#sampleData <- sampleData %>% select(expID, treatment)
sampleData <- tibble::rownames_to_column(sampleData, "sample")
rownames(sampleData) <- sampleData$sample
sampleData$sample<-as.factor(gsub("^.*\\_G","",sampleData$sample))
sampleData
write.table(sampleData, 'deseqExpDesign.csv', sep = ",")
m1 <- model.matrix(~ expID + expID:ind.n + expID:treatment, sampleData)
all.zero <- apply(m1, 2, function(x) all(x==0))
idx <- which(all.zero)
m1 <- m1[,-idx]
m1
m1 <- as.data.frame(m1)
print("Put the columns of the count data in the same order as rows names of the sample info, then make sure it worked")
raw_deseqW <- raw_deseqW[,unique(rownames(sampleData))]
head(raw_deseqW)
all(colnames(raw_deseqW) == rownames(sampleData))
print("Order the treatments so that it is sensible: non-inversion (control) -> inversion")
sampleData$treatment <- factor(sampleData$treatment, levels=c("0", "1"))
print("#############                 ")
print("## DESeq analysis accounting for experiment effect")
print("#############                 ")
print("Create the DEseq2DataSet object")
#sampleData <- sampleData[c(1,3,4,7),]
#raw_deseqW <- raw_deseqW %>% 
#              select(ATCC_GSE94978_1,BW25113_GSE73673_6,BW25113_GSE73673_7,K12DH_GSE98890_1)
head(raw_deseqW)
head(sampleData)
print("DE btwn strains")
deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~ strain)
deseq_results <- DESeq(deseq2Data)
print("done deseq")
deseq_results2 <- results(deseq_results, alpha=0.05)
# alpha = 0.05 is the  "cut-off" for significance (not really - I will
# discuss).
print("summary")
summary(deseq_results2)
print("DE btwn experiments")
deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~  expID)
deseq_results <- DESeq(deseq2Data)
print("done deseq")
deseq_results2 <- results(deseq_results,alpha=0.05)
# alpha = 0.05 is the  "cut-off" for significance (not really - I will
# discuss).
print("summary")
summary(deseq_results2)
print("DE btwn treatments")
deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~  treatment)
deseq_results <- DESeq(deseq2Data)
print("done deseq")
str(deseq_results)
deseq_results2 <- results(deseq_results, contrast = c("treatment", 1,0),  alpha=0.05)
# alpha = 0.05 is the  "cut-off" for significance (not really - I will
# discuss).
deseq_results2 <- deseq_results2[order(deseq_results2$pvalue),]
print("summary")
summary(deseq_results2)
print("check PCA")
for_pca <- rlog(deseq_results, 
                blind = TRUE)
pdf("DESeq_all_tax_2Nov20.pdf")
plotPCA(for_pca, 
        intgroup=c("expID"),
        ntop = 10000) 
dev.off()
print("check pvalue dist")
pdf("DESeqtreatment_pval_dist.pdf")
ggplot(data = as.data.frame(deseq_results2), aes(x = pvalue)) + 
  geom_histogram(bins = 100)
dev.off()
##deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~ strain + treatment)
##deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~ expID + expID:ind.n + expID:treatment)
##deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=m1, design= ~ expIDGSE94978:treatment1 + expIDGSE98890:treatment1)
##deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=m1)
##deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~ strain)
#deseq_results <- DESeq(deseq2Data)
#print("done deseq")
#deseq_results2 <- results(deseq_results, alpha=0.05)
## alpha = 0.05 is the  "cut-off" for significance (not really - I will
## discuss).
#print("summary")
#summary(deseq_results2)
##deseq2Data <- DESeqDataSetFromTximport(countData=raw_deseqW, colData=sampleData, design= ~ treatment)
## Coerce to a data frame
#deseq2ResDF <- as.data.frame(deseq_results2)
#
## Examine this data frame
#head(deseq2ResDF)
#
## Set a boolean column for significance
#deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)
#
## Plot the results similar to DEseq2
#p <- (ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant))
#+ geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish)
#+ scale_x_log10() 
#+ geom_hline(yintercept = 0, colour="tomato1",size=2) 
#+ labs(x="Mean of Normalized Counts", y="Log Fold Change") 
#+ scale_colour_manual(name="q-value", values=("Significant"="red"),na.value="grey50") 
#)
#
#pdf("DESeq_all_tax_2Nov20.pdf")
##plotMA(deseq_results2)
#p
#dev.off()
#
#
#
#
#
#
##sampleData <- sample_inf
##rownames(sampleData) <- sampleData$replicate
###keep <- c("treatment_A","exp","strain")
###sampleData <- sampleData[,keep]
##colnames(sampleData) <- c("replicates","strain","expID","treatment")
##sampleData$expID <- factor(sampleData$expID)
##sampleData$treatment <- factor(sampleData$treatment)
##sampleData$replicates <- factor(sampleData$replicates)
##sampleData$strain <- factor(sampleData$strain)
##sampleData
##print("Put the columns of the count data in the same order as rows names of the sample info, then make sure it worked")
##raw_deseqW <- raw_deseqW[,unique(rownames(sampleData))]
##head(raw_deseqW)
##all(colnames(raw_deseqW) == rownames(sampleData))
##print("Order the treatments so that it is sensible: non-inversion (basically control) -> inversion")
##sampleData$treatment <- factor(sampleData$treatment, levels=c("0", "1"))
##print("#############                 ")
##print("## DESeq looking at experiment effect")
##print("#############                 ")
##deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~ expID)
##deseq_results <- DESeq(deseq2Data)
##print("done deseq")
##deseq_results2 <- results(deseq_results, alpha=0.05)
### alpha = 0.05 is the  "cut-off" for significance (not really - I will
### discuss).
##print("summary")
##summary(deseq_results2)
##print("we see that a lot of genes are differentially expressed (up or down reg), therefore we need to account for experiment effects")
##
##print("#############                 ")
##print("## DESeq analysis accounting for experiment effect")
##print("#############                 ")
##print("Create the DEseq2DataSet object")
##deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~ expID + treatment)
##deseq_results <- DESeq(deseq2Data)
##print("done deseq")
##deseq_results2 <- results(deseq_results, alpha=0.05)
### alpha = 0.05 is the  "cut-off" for significance (not really - I will
### discuss).
##print("summary")
##summary(deseq_results2)
###deseq2Data <- DESeqDataSetFromTximport(countData=raw_deseqW, colData=sampleData, design= ~ treatment)
### Coerce to a data frame
##deseq2ResDF <- as.data.frame(deseq_results2)
##
### Examine this data frame
##head(deseq2ResDF)
##
### Set a boolean column for significance
##deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1,
##"Significant", NA)
##
### Plot the results similar to DEseq2
##p <- (ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant))
##+ geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish)
##+ scale_x_log10() 
##+ geom_hline(yintercept = 0, colour="tomato1",size=2) 
##+ labs(x="Mean of Normalized Counts", y="Log Fold Change") 
##+ scale_colour_manual(name="q-value", values=("Significant"="red"),na.value="grey50") 
##)
##
##pdf("DESeq_K12MG_DE_genes_8Oct2020.pdf")
###plotMA(deseq_results2)
##p
##dev.off()
##print("#############                 ")
##print("## DESeq without accounting for experiment effect")
##print("#############                 ")
##deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~ treatment)
##deseq_results <- DESeq(deseq2Data)
##print("done deseq")
##deseq_results2 <- results(deseq_results, alpha=0.05)
### alpha = 0.05 is the  "cut-off" for significance (not really - I will
### discuss).
##print("summary")
##summary(deseq_results2)
##
##
##
##
###deseq2Data <- DESeqDataSetFromTximport(countData=raw_deseqW, colData=sampleData, design= ~ treatment)
####load.model <- formula(~ expID+ strain + treatment)
###test_rep_effects <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design=load.model)
###test_rep_effects2 <- DESeq(test_rep_effects)
###print("done deseq")
###test_rep_effects2_results <- results(test_rep_effects2, alpha=0.05)
#### alpha = 0.05 is the  "cut-off" for significance (not really - I will
#### discuss).
###print("summary")
###summary(test_rep_effects2_results)
#### 2 genes which may show  evidence of lane effects, but this is a bit
#### incomplete for the full data set.
###print("plot results")
###pdf("exp_effect_deseq.pdf")
###plotMA(test_rep_effects2_results)
###dev.off()
###
####ians test for lane effects
###for_pca <- rlog(test_rep_effects2, 
###                blind = TRUE)
###
###pdf("PCA_effect_deseq.pdf")
###plotPCA(for_pca, 
###        intgroup=c("expID"),
###        ntop = 500) 
###dev.off()
####head(test_rep_effects2_results)
##
### We now fit the simple model
###rownames(raw_deseqW) <- raw_deseqW$gene_id
###raw_deseqW <- raw_deseqW[,c(-1,-2)]
###inversion <- raw_deseqW$inversion
###head(raw_deseqW)
###tail(raw_deseqW)
#### set up experimental design
###experimental_design = data.frame(
###  sample_names = colnames(raw_deseqW),  # sample name
####  individual = factor(colnames(exp_df)), # each individual strain
###  treatment = factor(inversion)  # inversion = 1 or no inversion = 0
####  lane = factor(parse_names[,3])      # Which lane on the Illumina flowcell.
###)
###DESeq_data <- DESeqDataSetFromMatrix(raw_deseqW, experimental_design, design = formula(~ treatment))
###
###
###
####read in gene mapping file to tell us which genes are homologous
###gmap_file <- as.character(args[8])
###gmap <- read.table(gmap_file, sep = "\t", header = FALSE)
####reshape
###gmaps <- subset(gmap, select = c("V2","V4","V6","V8"))
###colnames(gmaps) <- c("CP009072","CP009273","NC_010473","U00096000") 
###
####tmp_df <- gei_dat[,c("norm_exp","gene_id","strain","inversion")]
####DH <- tmp_df[which(tmp_df$strain == "K12DH"),]
####MG <- unique(tmp_df[which(tmp_df$strain == "K12MG"),])
####DH <- data.frame(DH[!duplicated(DH[ , c("gene_id")]),])
####levels(DH$strain) <- c(levels(DH$strain), "K12DH_i")
####DH$strain[DH$inversion == 1] <- "K12DH_i"
####MG <- data.frame(MG[!duplicated(MG[ , c("gene_id")]),])
####levels(MG$strain) <- c(levels(MG$strain), "K12MG_i")
####MG$strain[MG$inversion == 1] <- "K12MG_i"
####gene_names <- unique(MG$gene_id)
####DH <- unique(DH[which(DH$gene_id %in% gene_names),])
####exp_df <- rbind(DH,MG)
####exp_df <- spread(exp_df, strain, norm_exp)
####rownames(exp_df) <- exp_df$gene_id
####exp_df <- exp_df[,-1]
#####inversions_col <- exp_df$inversion
####exp_df <- exp_df[,-1]
#####********************
####t_count = 1
####t <- vector(mode='character', length = length(colnames(exp_df)))
####for (l in colnames(exp_df)) {
####  tmp_vec <- strsplit(l, "_")
####  tmp_l <- length(unlist(tmp_vec))
####  if (tmp_l > 1){
####    t[t_count] <- "inversion"
####  } else {
####    t[t_count] <- "none"
####  }
####  t_count <- t_count +1
####}
#####deal with NA values
####exp_df <- exp_df %>% replace(is.na(.), 0)
####
####print(t)
##### set up experimental design
####experimental_design = data.frame(
####  sample_names = colnames(exp_df),  # sample name
#####  individual = factor(colnames(exp_df)), # each individual strain
####  treatment = factor(t)  # inversion = 1 or no inversion = 0
#####  lane = factor(parse_names[,3])      # Which lane on the Illumina flowcell.
####)
####
####DESeq_data <- DESeqDataSetFromMatrix(exp_df, experimental_design, design = formula(~ treatment))
####
###
###
###############
#### ALL inversions
###############
###
###
###
####t.test(gei_dat$norm_exp[gei_dat$block == "Block800",], gei_dat$norm_exp[gei_dat$block == "Block800",])
####
####test_inver <- gei_dat %>% filter(inversion == 1)
####test_revcomp <- gei_dat %>% filter(rev_comp == 1)
####
####gei_dat %>% filter(block == "Block800") %>%
####  filter(rev_comp == 1)
####gei_dat %>% filter(block == "Block800") %>%
####  filter(rev_comp == 0)
####
####gei_dat %>% select(block, inversion, rev_comp, norm_exp) %>%
####  group_by(block, rev_comp)
####
####gei_dat %>%
####  select(filter(rev_comp == 1),block, norm_exp, inversion)
####
####
####gei_dat %>%
####  filter(block == "Block800" && block == "Block799") %>%
####  group_by(block) %>%
####  do(tidy(t.test(norm_exp ~ rev_comp, data = .)))
####
#####inverted rows we want
####gei_dat[gei_dat$block %in% gei_dat$block[gei_dat$rev_comp == 1],]
####gei_dat %>%
####  filter(rev_comp == 1)
####
##
##
##################################################
print("HNS BINDING")
print("##################################################")


#read in the data
#and format

# K12 MG1655
print("Grainger 2006 cutoff")
graniger_dat <- read.csv("../HNS_protein/raw_data_files/Grainger_2006_HNS_binding_cutoff_coding_and_noncoding.csv", header = TRUE)
sub_g_dat <- graniger_dat %>%
             select(Gene.name,Coordinates)
colnames(sub_g_dat) <- c("gene_name","Coords")
coords <- str_split_fixed(sub_g_dat$Coords, "-", 2)
colnames(coords) <- c("start","end")
sub_g_dat <- cbind(sub_g_dat,coords)
sub_g_dat <- sub_g_dat %>%
             select(gene_name,start,end) %>%
             filter(gene_name != "Between convergent ORFs") %>%
             filter(start != "N/A")
gene_inf <- read.table("../Genomes/Ecoli_K12_MG1655_chrom_U00096_gene_info.txt", header = TRUE)
gene_inf <- gene_inf %>% select(gbk_strand,gbk_gene_id)
colnames(gene_inf) <- c("gbk_strand","gene_name")
head(gene_inf)
sub_g_dat <- merge(sub_g_dat,gene_inf, by = "gene_name")
head(sub_g_dat)

#below has to be combined with gene info file for K12MG
# W3
print("Higashi 2016 coding")
file <- "../HNS_protein/raw_data_files/Higashi_2016_HNS_binding_sites_coding.csv"
higashi_dat <- read.csv(file, header = TRUE)
h_hgt_dat <- higashi_dat[,c(2,5)]
colnames(h_hgt_dat) <- c("gene_name","HGT")
gene_inf <- read.table("../Genomes/Ecoli_K12_MG1655_chrom_U00096_gene_info.txt", header = TRUE)
gene_inf <- gene_inf %>% select(gbk_start,gbk_end,gbk_midpoint,gbk_gene_id)
colnames(gene_inf) <- c("start","end","gbk_midpoint","gene_name")
h_hgt_dat  <- merge(h_hgt_dat,gene_inf, by = "gene_name")
sub_h_dat <- higashi_dat[,c(2,6,7,8)]
colnames(sub_h_dat) <- c("gene_name","HNS_binding","HNS_cutoff","HNS_transcript")
gene_inf <- read.table("../Genomes/Ecoli_K12_MG1655_chrom_U00096_gene_info.txt", header = TRUE)
gene_inf <- gene_inf %>% select(gbk_start,gbk_end,gbk_midpoint,gbk_gene_id)
colnames(gene_inf) <- c("start","end","gbk_midpoint","gene_name")
sub_h_df <- merge(sub_h_dat,gene_inf, by = "gene_name")
hns_viz_h1 <- sub_h_df
head(sub_h_df)
#non-coding Higashi
print("Higashi 2016 non-coding")
file <- "../HNS_protein/raw_data_files/higashi_nc_dat.csv"
higashi_nc_dat <- read.csv(file, header = TRUE)
sub_hnc_dat <- higashi_nc_dat[,c(6,7,8,10,11)]
head(sub_hnc_dat)

# W3
print("Ueda 2013")
ueda_dat <- read.csv("../HNS_protein/raw_data_files/Ueda_2013_HNS_binding_sites_W3110.csv", header = TRUE)
sub_u_dat <- ueda_dat[,c(1,2)]
colnames(sub_u_dat) <- c("start","end")
sub_u_dat <- na.omit(sub_u_dat)
head(sub_u_dat)


print("Lang 2007")
print("all data is an HNS binding site")
lang_dat <- read.csv("../HNS_protein/raw_data_files/Lang_2007_HNS_binding.csv", header = TRUE)
sub_l_df <- na.omit(lang_dat)
head(sub_l_df)

print("Oshima 2006")
oshima_dat <- read.csv("../HNS_protein/raw_data_files/Oshima_2006_HGT_HNS.csv", header = TRUE)
oshima_dat$Nakamura[is.na(oshima_dat$Nakamura)] <- 0
oshima_dat$Lawrence[is.na(oshima_dat$Lawrence)] <- 0
oshima_dat$hns[is.na(oshima_dat$hns)] <- 0
oshima_dat <- oshima_dat %>% drop_na(start)
sub_o_df <- filter(oshima_dat, hns == 1)
head(sub_o_df)

print("combine HNS binary info to inversion df")
print("THIS IS NOT BIDIRECTIONAL BC IT IS THE START AND ENDS AND NOT THE MIDPOINT!")
head(inver_dat_bidir)
inver_cor_d <- inver_dat_bidir %>% filter(strain == "K12MG") %>%
            select(gbk_start,gbk_end,inversion,block)
#            select(start,end)
#colnames(inver_cor_d) <- c("start1","end1")
colnames(inver_cor_d) <- c("start1","end1","Inversion","block")
inver_cor_d$Inversion <- as.character(inver_cor_d$Inversion)
inver_cor_d <- unique(inver_cor_d)
head(inver_cor_d)
write.table(inver_cor_d, 'inver_cor_data.csv', sep = "\t")

print("################################################################################")
print("CORRELATION TESTS FOR GRAINGER 2006 cutoff DATA")
print("################################################################################")
print("hns")
head(sub_g_dat)
sub_g_dat$start <- as.numeric(as.character(sub_g_dat$start))
sub_g_dat$end <- as.numeric(as.character(sub_g_dat$end))
hns_dat <- within(sub_g_dat, midpoint <- (end + start) /2)
colnames(hns_dat)[colnames(hns_dat) == "midpoint"] <- "midpoint"
class2 <- rep("G_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding")
inver_cor_d$G_HNS_binding <- as.integer(inver_cor_d$G_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start1 == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$G_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
#filter by K12 strain and ensure that we are only looking at the wilcoxon test results
cor_dat <- cor_dat %>% filter(strain == "K12MG") %>% filter(test == "block_w_pvalue")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
summary(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "G_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$G_HNS_binding <- as.integer(cor_dat$G_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
cor_dat$inversion <- as.integer(cor_dat$inversion)
head(cor_dat)
tail(cor_dat)
summary(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$G_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$G_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR UEDA 2013 cutoff DATA")
print("################################################################################")
print("hns")
head(sub_u_dat)
sub_u_dat$start <- as.numeric(as.character(sub_u_dat$start))
sub_u_dat$end <- as.numeric(as.character(sub_u_dat$end))
hns_dat <- within(sub_u_dat, midpoint <- (end + start) /2)
colnames(hns_dat)[colnames(hns_dat) == "midpoint"] <- "midpoint"
class2 <- rep("U_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding")
inver_cor_d$U_HNS_binding <- as.integer(inver_cor_d$U_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$U_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "U_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$U_HNS_binding <- as.integer(cor_dat$U_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$U_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$U_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR Higashi 2016 criteria 1: HNS binding DATA")
print("################################################################################")
print("hns")
head(sub_h_df)
sub_h_df <- sub_h_df %>%
            filter(HNS_binding == TRUE)
sub_h_df$start <- as.numeric(as.character(sub_h_df$start))
sub_h_df$end <- as.numeric(as.character(sub_h_df$end))
hns_dat <- sub_h_df
class2 <- rep("H1_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding","H1_HNS_binding")
inver_cor_d$H1_HNS_binding <- as.integer(inver_cor_d$H1_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$H1_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "H1_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$H1_HNS_binding <- as.integer(cor_dat$H1_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$H1_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$H1_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR Higashi 2016 criteria 1: HNS binding WITH NON-CODING DATA")
print("################################################################################")
print("hns")
sub_hnc_dat <- sub_hnc_dat %>%
            filter(HNS == TRUE) %>%
            select(start,end)
sub_hnc_dat$start <- as.numeric(as.character(sub_hnc_dat$start))
sub_hnc_dat$end <- as.numeric(as.character(sub_hnc_dat$end))
#combine coding higashi with non-cod criteria 1
sub_h_df <- sub_h_df %>%
            filter(HNS_transcript == TRUE)
tmp_df <- sub_h_df %>% select(start,end)
tmp_df$start <- as.numeric(tmp_df$start)
tmp_df$end <- as.numeric(tmp_df$end)
tmp_df <- rbind(tmp_df,sub_hnc_dat)
hns_dat <- tmp_df
class2 <- rep("H1nc1_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding","H1_HNS_binding","H1nc1_HNS_binding")
inver_cor_d$H1nc1_HNS_binding <- as.integer(inver_cor_d$H1nc1_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$H1nc1_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "H1nc1_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$H1nc1_HNS_binding <- as.integer(cor_dat$H1nc1_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$H1nc1_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$H1nc1_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR Higashi 2016 criteria 1: HNS binding WITH NON-CODING 2 DATA")
print("################################################################################")
print("hns")
sub_hnc_dat <- sub_hnc_dat %>%
            filter(TtoT == TRUE) %>%
            select(start,end)
sub_hnc_dat$start <- as.numeric(as.character(sub_hnc_dat$start))
sub_hnc_dat$end <- as.numeric(as.character(sub_hnc_dat$end))
#combine coding higashi with non-cod criteria 1
sub_h_df <- sub_h_df %>%
            filter(HNS_transcript == TRUE) 
tmp_df <- sub_h_df %>% select(start,end)
tmp_df$start <- as.numeric(tmp_df$start)
tmp_df$end <- as.numeric(tmp_df$end)
tmp_df <- rbind(tmp_df,sub_hnc_dat)
hns_dat <- tmp_df
class2 <- rep("H1nc2_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding","H1_HNS_binding","H1nc1_HNS_binding","H1nc2_HNS_binding")
inver_cor_d$H1nc2_HNS_binding <- as.integer(inver_cor_d$H1nc2_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$H1nc2_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "H1nc2_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$H1nc2_HNS_binding <- as.integer(cor_dat$H1nc2_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$H1nc2_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$H1nc2_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR Higashi 2016 criteria 1: HNS binding WITH NON-CODING 3 DATA")
print("################################################################################")
print("hns")
sub_hnc_dat <- sub_hnc_dat %>%
            filter(known_promoter == TRUE) %>%
            select(start,end)
sub_hnc_dat$start <- as.numeric(as.character(sub_hnc_dat$start))
sub_hnc_dat$end <- as.numeric(as.character(sub_hnc_dat$end))
#combine coding higashi with non-cod criteria 1
sub_h_df <- sub_h_df %>%
            filter(HNS_transcript == TRUE) 
tmp_df <- sub_h_df %>% select(start,end)
tmp_df$start <- as.numeric(tmp_df$start)
tmp_df$end <- as.numeric(tmp_df$end)
tmp_df <- rbind(tmp_df,sub_hnc_dat)
hns_dat <- tmp_df
class2 <- rep("H1nc3_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding","H1_HNS_binding","H1nc1_HNS_binding","H1nc2_HNS_binding","H1nc3_HNS_binding")
inver_cor_d$H1nc3_HNS_binding <- as.integer(inver_cor_d$H1nc3_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$H1nc3_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "H1nc3_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$H1nc3_HNS_binding <- as.integer(cor_dat$H1nc3_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$H1nc3_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$H1nc3_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR Higashi 2016 criteria 2: HNS cutoff DATA")
print("################################################################################")
print("hns")
sub_h_df <- sub_h_df %>%
            filter(HNS_cutoff == TRUE)
sub_h_df$start <- as.numeric(as.character(sub_h_df$start))
sub_h_df$end <- as.numeric(as.character(sub_h_df$end))
hns_dat <- sub_h_df
class2 <- rep("H2_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding","H1_HNS_binding","H1nc1_HNS_binding","H1nc2_HNS_binding","H1nc3_HNS_binding","H2_HNS_binding")
inver_cor_d$H2_HNS_binding <- as.integer(inver_cor_d$H2_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$H2_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "H2_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$H2_HNS_binding <- as.integer(cor_dat$H2_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$H2_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$H2_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR Higashi 2016 criteria 3: HNS transcript DATA")
print("################################################################################")
print("hns")
sub_h_df <- sub_h_df %>%
            filter(HNS_transcript == TRUE)
sub_h_df$start <- as.numeric(as.character(sub_h_df$start))
sub_h_df$end <- as.numeric(as.character(sub_h_df$end))
hns_dat <- sub_h_df
class2 <- rep("H3_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding","H1_HNS_binding","H1nc1_HNS_binding","H1nc2_HNS_binding","H1nc3_HNS_binding","H2_HNS_binding","H3_HNS_binding")
inver_cor_d$H3_HNS_binding <- as.integer(inver_cor_d$H3_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$H3_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "H3_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$H3_HNS_binding <- as.integer(cor_dat$H3_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$H3_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$H3_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR Lang 2007: composite data")
print("################################################################################")
print("hns")
sub_l_df$start <- as.numeric(as.character(sub_l_df$start))
sub_l_df$end <- as.numeric(as.character(sub_l_df$end))
hns_dat <- sub_l_df
class2 <- rep("L_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <-
c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding","H1_HNS_binding","H1nc1_HNS_binding","H1nc2_HNS_binding","H1nc3_HNS_binding","H2_HNS_binding","H3_HNS_binding","L_HNS_binding")
inver_cor_d$L_HNS_binding <- as.integer(inver_cor_d$L_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$L_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "L_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$L_HNS_binding <- as.integer(cor_dat$L_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$L_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$L_HNS_binding,
         method = "pearson")

print("################################################################################")
print("CORRELATION TESTS FOR Oshima 2006: composite data")
print("################################################################################")
print("hns")
sub_o_df$start <- as.numeric(as.character(sub_o_df$start))
sub_o_df$end <- as.numeric(as.character(sub_o_df$end))
hns_dat <- sub_o_df
class2 <- rep("O_HNS_Binding",length(hns_dat$start))
hns_dat <- cbind(hns_dat,class2)
hns_cor_d <- hns_dat
head(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% select(start,end,class2)
hns_cor_d <- hns_cor_d %>% select(start,end)
summary(hns_cor_d)
#hns_cor_d <- hns_cor_d %>% 
#    mutate(class2 = recode(class2, 
#                      "No_HNS_Binding" = "0", 
#                      "HNS_Binding" = "1"))
hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
head(hns_cor_d)
write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")


ir1 = with(inver_cor_d, IRanges(start1, end1))
ir2 = with(hns_cor_d, IRanges(start, end))
inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
inver_cor_d[which(inver_cor_d$start1 == 383921),]
colnames(inver_cor_d) <- c("start1","end1","Inversion","block","G_HNS_binding","U_HNS_binding","H1_HNS_binding","H1nc1_HNS_binding","H1nc2_HNS_binding","H1nc3_HNS_binding","H2_HNS_binding","H3_HNS_binding","L_HNS_binding", "O_HNS_binding")
inver_cor_d$O_HNS_binding <- as.integer(inver_cor_d$O_HNS_binding)
inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
head(inver_cor_d)
inver_cor_d[which(inver_cor_d$start == 383921),]

print("################################################################################")
print("correlation test btwn inversion/non-inversion and hns binding")
print("################################################################################")
cor.test(inver_cor_d$Inversion, inver_cor_d$O_HNS_binding,
         method = "pearson")

print("get new df with hns binding + inversions + sig block info")
cor_dat <- cor_dat %>% filter(strain == "K12MG")
cor_dat <- cor_dat %>%
    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
head(cor_dat)
ir1 = with(cor_dat, IRanges(gbk_start, gbk_end))
cor_dat$overlap = countOverlaps(ir1, ir2) != 0
names(cor_dat)[names(cor_dat)=="overlap"] <- "O_HNS_binding"
cor_dat[which(cor_dat$start == 383921),]
cor_dat$O_HNS_binding <- as.integer(cor_dat$O_HNS_binding)
cor_dat$sig <- as.integer(cor_dat$sig)
head(cor_dat)


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
print("################################################################################")
cor.test(cor_dat$sig, cor_dat$O_HNS_binding,
         method = "pearson")


print("################################################################################")
print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
print("################################################################################")
#cor_dat <- cor_dat %>% filter(inversion == 1)
cor.test(cor_dat$inversion, cor_dat$O_HNS_binding,
         method = "pearson")







print("checking overlap between HNS binding sites in all datasets")
print("total HNS binding sites for each dataset")
inver_cor_d %>%
gather(x, value, G_HNS_binding:O_HNS_binding)%>%
group_by(x)%>%
tally(value == 1)

tmp_d <- inver_cor_d%>% select(G_HNS_binding:O_HNS_binding) %>% 
         filter( G_HNS_binding ==1 | U_HNS_binding == 1 | H1_HNS_binding == 1 | H1nc1_HNS_binding == 1| H1nc2_HNS_binding == 1 | H1nc3_HNS_binding == 1| H2_HNS_binding == 1 | H3_HNS_binding == 1 | L_HNS_binding == 1 |O_HNS_binding) 
print("total number of genes/blocks")
length(tmp_d$G_HNS_binding)
print("total number of rows where HNS binding is the same for all datasets")
sum(apply(tmp_d, 1, function(x) length(unique(x))==1))


print("HNS and expression in Oshima, Lang and Higashi (H1)")
#get blocks where Higashi (H1) or Lang had HNS bound
head(inver_cor_d)
l_h_blocks <- inver_cor_d %>%
         filter(H1_HNS_binding == 1 | L_HNS_binding == 1 | O_HNS_binding ==1) %>%
         select(block)
l_h_blocks <- unique(l_h_blocks)
l_h_blocks <- as.character(sort(l_h_blocks$block))
#l_h_blocks <- l_h_blocks[-184]
class(l_h_blocks)
print("#get gei_dat that is just ^ blocks")
#gei_dat %>% filter(block == "Block789")
#gei_dat %>% filter(block == "Block129")
#gei_dat %>% filter(block == "Block141")
l_h_dat <- gei_dat %>%
           filter(block %in% l_h_blocks)
g_l <- as.character(sort(unique(l_h_dat$block)))
class(g_l)
print("check if all HNS blocks were grabbed")
identical(g_l, l_h_blocks)


print("wilcoxon sign ranked test to see if there is a diff in exp btwn HNS bound inversions and HNS bound non-inversions")
wilcox.test(l_h_dat$norm_exp[l_h_dat$inversion == 1], l_h_dat$norm_exp[l_h_dat$inversion == 0])

print("Avg inver exp with HNS binding")
mean(l_h_dat$norm_exp[l_h_dat$inversion == 1])
print("Avg NON-inver exp with HNS binding")
mean(l_h_dat$norm_exp[l_h_dat$inversion == 0])

print("wilcoxon sign ranked test to see if there is a diff in exp btwn sig HNS bound inversions and non-sig HNS bound non-inversions")
#make df with sig diff in gene exp btwn inversions
sig_gei <- l_h_dat %>%
    mutate(sig = if_else(block_w_pvalue <= 0.05, 1, 0))
head(sig_gei)
wilcox.test(sig_gei$norm_exp[sig_gei$sig == 1], sig_gei$norm_exp[sig_gei$sig == 0])


print("create column with HNS binding and non-binding")
gei_dat_hns <- gei_dat %>%
               mutate(HNS = ifelse(block %in% l_h_blocks, 1, 0))
print("wilcoxon sign ranked test to see if there is a diff in exp btwn HNS bound and HNS un-bound blocks")
wilcox.test(gei_dat_hns$norm_exp[gei_dat_hns$HNS == 1], gei_dat_hns$norm_exp[gei_dat_hns$HNS == 0])
print("Avg exp with HNS binding")
mean(gei_dat_hns$norm_exp[gei_dat_hns$HNS == 1])
print("Avg exp WITHOUT HNS binding")
mean(gei_dat_hns$norm_exp[gei_dat_hns$HNS == 0])

print("overlap btwn Higashi (H1) and Lang and Oshima")
l_h_overlap <- inver_cor_d %>%
         filter(H1_HNS_binding == 1 & L_HNS_binding == 1 & O_HNS_binding == 1)
length(l_h_overlap$Inversion)

print("##################")
print("##### HGT and inversions (in the context of K12MG) ######")
print("##################")
print("*****filtering HGT data to include all categories of HGT (including nearby genes that are HGT")
print("Higashi HGT")
h_hgt_dat <- h_hgt_dat %>% filter(HGT != "x")
hgt_inver <- gei_dat %>% filter(strain == "K12MG")
ir1 = with(hgt_inver, IRanges(gbk_start, gbk_end))
ir2 = with(h_hgt_dat, IRanges(start, end))
hgt_inver$HGT = countOverlaps(ir1, ir2) != 0
hgt_inver$HGT <- as.numeric(hgt_inver$HGT)
hgt_inver <- hgt_inver %>%
               mutate(sig = ifelse(block_w_pvalue <= 0.05, 1, 0))
head(hgt_inver)
summary(hgt_inver)
print("################################################################################")
print("correlation test btwn inversion/non-inversion and HGT")
print("################################################################################")
cor.test(hgt_inver$inversion, hgt_inver$HGT,
         method = "pearson")
print("################################################################################")
print("correlation test btwn sig inversion/non-sig inversion and HGT")
print("################################################################################")
tmp_hgt <- filter(hgt_inver, inversion ==1)
cor.test(tmp_hgt$sig, tmp_hgt$HGT, method = "pearson")

print("*****filtering HGT data to include all categories of HGT ")
print("Oshima HGT")
o_hgt_dat <- sub_o_df %>% filter(Nakamura == 1 | Lawrence == 1)
head(o_hgt_dat)
hgt_inver <- gei_dat %>% filter(strain == "K12MG")
ir1 = with(hgt_inver, IRanges(gbk_start, gbk_end))
ir2 = with(o_hgt_dat, IRanges(start, end))
hgt_inver$HGT = countOverlaps(ir1, ir2) != 0
hgt_inver$HGT <- as.numeric(hgt_inver$HGT)
hgt_inver <- hgt_inver %>%
               mutate(sig = ifelse(block_w_pvalue <= 0.05, 1, 0))
head(hgt_inver)
summary(hgt_inver)
print("################################################################################")
print("correlation test btwn inversion/non-inversion and HGT")
print("################################################################################")
cor.test(hgt_inver$inversion, hgt_inver$HGT,
         method = "pearson")
print("################################################################################")
print("correlation test btwn sig inversion/non-sig inversion and HGT")
print("################################################################################")
tmp_hgt <- filter(hgt_inver, inversion ==1)
cor.test(tmp_hgt$sig, tmp_hgt$HGT, method = "pearson")





print("HNS binding and inversions viz")
print("Using Higashi criteria 1 for HNS binding")
inver_cor_d <- within(inver_cor_d, midpoint <- (end1 + start1) /2)
inver_cor_d$midpoint = inver_cor_d$midpoint / 1000000
fake_val <- rep(10,length(inver_cor_d$midpoint))
inver_cor_d <- cbind(inver_cor_d, fake_val)
inver_cor_d <- filter(inver_cor_d, H1_HNS_binding == 1)
print("inver_cor_d")
head(inver_cor_d)
summary(inver_cor_d)
fake_val <- rep(10,length(cor_dat$block))
fake_val2 <- rep(9.95,length(cor_dat$block))
hns_inver <- cbind(cor_dat,fake_val,fake_val2)
hns_inver <- hns_inver %>% filter(inversion == 1)
hns_inver$non_bidir_midpoint = hns_inver$non_bidir_midpoint / 1000000
print("hns_inver")
head(hns_inver)
summary(hns_inver)
hns_sig <- hns_inver %>% filter(sig == 1)
print("hns_sig")
head(hns_sig)
#hns_bind <- hns_inver %>% filter(H1_HNS_binding == 1)
hns_bind <- hns_viz_h1 %>% filter(HNS_binding == TRUE)
hns_bind$gbk_midpoint = hns_bind$gbk_midpoint / 1000000
fake_val <- rep(10,length(hns_bind$gbk_midpoint))
hns_bind <- cbind(hns_bind, fake_val)
head(hns_bind)


override.shape <- c(1,2,16)
override.size <- c(5,5,10)
p <- (ggplot(hns_inver, aes(x=non_bidir_midpoint, y=fake_val, colour = strain))
   #inversions in grey
   + geom_point(size = 10, aes(colour = "#BEBEBE"))
   + geom_point(data = inver_cor_d, aes(x=midpoint, y=fake_val,  color = "#2E294E"), size = 3,shape = 1)
   # sig inversions
   + geom_point(data = hns_sig, aes(x=non_bidir_midpoint, y=fake_val2,  color = "#AA5042"), size = 3,shape = 2)
   + labs(title= "H-NS Binding and Inversions",x = "Genomic Location (Mbp)", y = "") 
   + scale_color_manual(values = c("#2E294E","#AA5042","#BEBEBE"), labels = c("H-NS Binding","Significant Inversion","Inversion"))
   + guides(colour = guide_legend(override.aes = list(shape = override.shape, size = override.size, fill=NA)))
   + scale_y_continuous(limits = c(9.5, 10.5))
   + theme(axis.text.y = element_blank(),
          legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
          )
)
pdf("hns_all_inversions.pdf")
p
dev.off()


###############NEW HNS GRAPH Density plot
fake2 <- rep(0.025,length(hns_inver$non_bidir_midpoint)) 
hns_inver <- cbind(hns_inver,fake2)
p <- (ggplot(inver_cor_d, aes(x=midpoint))
  + geom_density(alpha=0.5, aes(fill = "#2E294E"))
  + geom_density(data = hns_sig, aes(x=non_bidir_midpoint, fill = "#AA5042"),alpha=0.4)
#  + geom_density(data = hns_inver, aes(x=non_bidir_midpoint),alpha=0.4, fill = "#BEBEBE")
  + geom_point(data = hns_inver, aes(x = non_bidir_midpoint,y=fake2,colour = "#42474C"), size=10)
   + labs(title= "H-NS Binding and Inversions",x = "Genomic Location (Mbp)", y = "Density") 
   + scale_color_manual(values = c("#42474C"), labels = c("Inversion"))
   + scale_fill_manual(values = c("#2E294E","#AA5042","#42474C" ), labels = c("H-NS Binding","Significant Inversion"))
#  + scale_fill_manual(name = "", labels = c("H-NS Binding", "Significant Inversions"))
  + theme(legend.position = "top")
)
pdf("hns_all_inversions_hist.pdf")
p
dev.off()

##################################################################################
##################################################################################
print("### HNS GRAPH HISTOGRAM")
##################################################################################
##################################################################################

#all HNS binding data is the hns_bind df
head(hns_bind)
hns_bind <- hns_bind %>% select(gbk_midpoint,HNS_binding)
hns_bind$HNS_binding <- as.numeric(hns_bind$HNS_binding)
hns_bind$gbk_midpoint <- hns_bind$gbk_midpoint *1000000
head(hns_bind)
#new min and max for the scaled positions (rounded to chunk len
#specified in args)
chunklen <- as.numeric(args[10])
print("chunklen")
chunklen
#nmin_pos <- round_any(min(hns_sig$non_bidir_midpoint), chunklen, f=floor)
nmin_pos <- 0
print("min (rounded)")
nmin_pos
#nmax_pos <- max(hns_sig$non_bidir_midpoint)
nmax_pos <- as.numeric(args[2])
print("max (rounded)")
nmax_pos

#expty vec to hold the rows that split the dat into X kb chunks
subs_rows_to_split_dat <- vector()
#create vector of chunks
chunk_len_of_genome <- round_any(nmax_pos, chunklen, f=ceiling)
chunk_len_of_genome
chunks <- seq(nmin_pos, chunk_len_of_genome, chunklen)
chunks
##checking if the lengths are equal for the later bind
#if (length(unlist(tot_gene_subs_10kb)) == length(chunks)) {
#    chunks_pos_NOzero <- chunks
#} else {
    chunks_pos_NOzero <- chunks[which(chunks != 0)]
#}
chunks_pos_NOzero

#add in fake data for chunks so it will split properly
#len_of_chunks <- length(chunks_pos_NOzero)
len_of_chunks <- length(chunks)
fhns <- rep(as.numeric(1), len_of_chunks)
#fake_dat <- cbind(chunks_pos_NOzero, fhns)
fake_dat <- cbind(chunks, fhns)
colnames(fake_dat) <- c("gbk_midpoint","HNS_binding")
head(fake_dat)
head(hns_bind)
#add fake data to orig data
hns_bind <- rbind(hns_bind,fake_dat)
head(hns_bind)
tail(hns_bind)
hns_bind$gbk_midpoint <- as.numeric(hns_bind$gbk_midpoint)
hns_bind$HNS_binding <- as.numeric(hns_bind$HNS_binding)
head(hns_bind)


#order the data by positions, only if its NOT strep
if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name == "S.meliloti"| bac_name == "Streptomyces") {
        ord_pos <- order(hns_bind$gbk_midpoint)
        hns_bind$HNS_binding <- hns_bind$HNS_binding[ord_pos]
        hns_bind$gbk_midpoint <- hns_bind$gbk_midpoint[ord_pos]
}
print("SAVED BIDIRECTIONAL DATA TO FILE")
write.table(hns_bind, 'hns_bind.csv', sep = ",")

print("checking order")
head(hns_bind)
tail(hns_bind)

sig_rows_to_split_dat <- vector()
for (i in chunks) {
  subs_rows <-
which(abs(hns_bind$gbk_midpoint-i)==min(abs(hns_bind$gbk_midpoint-i)))
# finding the closest number to each 10kb without going over it
  subs_max_row <- max(subs_rows)
  subs_rows_to_split_dat <- c(subs_rows_to_split_dat, subs_max_row)
}#for

#split data by the above specified rows
#sometimes the closest number was technically in the next 10kb chunk of seq
#but its fine.
data_just_subs <- hns_bind$HNS_binding
list_dat_sets_just_subs <- split(data_just_subs, findInterval(1:nrow(hns_bind), subs_rows_to_split_dat))
list_dat_sets_just_subs
print("split subs into 10kb chunks")


##########
#add up all subs in each 10kb section
##########
list_tot_add_subs_10kb <- lapply(list_dat_sets_just_subs, sum)
tot_gene_subs_10kb <- as.data.frame(matrix(unlist(list_tot_add_subs_10kb), byrow = F))
print("head/tail of total_subs")
head(tot_gene_subs_10kb)
tail(tot_gene_subs_10kb)
nrow(tot_gene_subs_10kb)
#length(chunks_pos_NOzero)
length(chunks)
class(tot_gene_subs_10kb)
#remove last element if lengths do not match (for cbind later)
#if (nrow(tot_gene_subs_10kb) == length(chunks_pos_NOzero)) {
if (nrow(tot_gene_subs_10kb) == length(chunks)) {
    tot_gene_subs_10kb <- tot_gene_subs_10kb
} else {
    tot_gene_subs_10kb <- tot_gene_subs_10kb[-nrow(tot_gene_subs_10kb),]
    tot_gene_subs_10kb <- tot_gene_subs_10kb[-1]
}

#tot_gene_subs_10kb <- as.data.frame(cbind(tot_gene_subs_10kb, chunks_pos_NOzero))
tot_gene_subs_10kb <- as.data.frame(cbind(tot_gene_subs_10kb, chunks))
head(tot_gene_subs_10kb)
tail(tot_gene_subs_10kb)
#subtract 1 from HNS binding column (bc of the fake data that was added)
tot_gene_subs_10kb$V1 <- tot_gene_subs_10kb$V1 -1
#add 10000 to each bar so that it represents all pts within that chunk of the genome
tot_gene_subs_10kb$chunks <- tot_gene_subs_10kb$chunks +10000
print("tot_gene_subs_10kb")
write.csv(tot_gene_subs_10kb)

#####################################################
print("sig inversions hist data")
#####################################################
#all sig inversions data is the inver_dat_bidir df
inver_dat_bidir$block_w_pvalue <- as.numeric(inver_dat_bidir$block_w_pvalue)
sig_id <- inver_dat_bidir %>%
          select(start,end,strain,block,inversion,block_w_pvalue) %>%
          filter(strain == "K12MG") %>%
          filter(block_w_pvalue <= 0.05)
sig_id <- unique(sig_id)
sig_id <- within(sig_id, midpoint <- (end + start) /2)
colnames(sig_id)[colnames(sig_id) == "midpoint"] <- "midpoint"
head(sig_id)
sig_id <- sig_id %>% select(inversion,midpoint)
head(sig_id)

#add in fake data for chunks so it will split properly
#len_of_chunks <- length(chunks_pos_NOzero)
len_of_chunks <- length(chunks)
finv <- rep(as.numeric(1), len_of_chunks)
#fake_dat <- cbind(chunks_pos_NOzero, fhns)
fake_dat <- cbind(chunks, finv)
colnames(fake_dat) <- c("midpoint","inversion")
head(fake_dat)
#add fake data to orig data
sig_bind <- rbind(sig_id,fake_dat)
head(sig_bind)
tail(sig_bind)
sig_bind$gbk_midpoint <- as.numeric(sig_bind$midpoint)
sig_bind$HNS_binding <- as.numeric(sig_bind$inversion)


#order the data by positions, only if its NOT strep
if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name == "S.meliloti"| bac_name == "Streptomyces") {
        ord_pos <- order(sig_bind$gbk_midpoint)
        sig_bind$HNS_binding <- sig_bind$HNS_binding[ord_pos]
        sig_bind$gbk_midpoint <- sig_bind$gbk_midpoint[ord_pos]
}
print("SAVED BIDIRECTIONAL DATA TO FILE")
write.table(sig_bind, 'sig_bind.csv', sep = ",")

print("checking order")
head(sig_bind)
tail(sig_bind)

for (i in chunks) {
  sig_rows <-
which(abs(sig_bind$gbk_midpoint-i)==min(abs(sig_bind$gbk_midpoint-i)))
# finding the closest number to each 10kb without going over it
  sig_max_row <- max(sig_rows)
  sig_rows_to_split_dat <- c(sig_rows_to_split_dat, sig_max_row)
}#for

#split data by the above specified rows
#sometimes the closest number was technically in the next 10kb chunk of seq
#but its fine.
data_just_sig <- sig_bind$inversion
list_dat_sets_just_sig <- split(data_just_sig, findInterval(1:nrow(sig_bind), sig_rows_to_split_dat))
list_dat_sets_just_sig
print("split sig into 10kb chunks")


##########
#add up all sig in each 10kb section
##########
list_tot_add_sig_10kb <- lapply(list_dat_sets_just_sig, sum)
tot_gene_sig_10kb <- as.data.frame(matrix(unlist(list_tot_add_sig_10kb), byrow = F))
print("head/tail of total_sig")
head(tot_gene_sig_10kb)
tail(tot_gene_sig_10kb)
nrow(tot_gene_sig_10kb)
#length(chunks_pos_NOzero)
length(chunks)
class(tot_gene_sig_10kb)
#remove last element if lengths do not match (for cbind later)
#if (nrow(tot_gene_subs_10kb) == length(chunks_pos_NOzero)) {
if (nrow(tot_gene_sig_10kb) == length(chunks)) {
    tot_gene_sig_10kb <- tot_gene_sig_10kb
} else {
    tot_gene_sig_10kb <- tot_gene_sig_10kb[-nrow(tot_gene_sig_10kb),]
    tot_gene_sig_10kb <- tot_gene_sig_10kb[-1]
}

#tot_gene_sig_10kb <- as.data.frame(cbind(tot_gene_sig_10kb, chunks_pos_NOzero))
tot_gene_sig_10kb <- as.data.frame(cbind(tot_gene_sig_10kb, chunks))
head(tot_gene_sig_10kb)
tail(tot_gene_sig_10kb)
#subtract 1 from HNS binding column (bc of the fake data that was added)
tot_gene_sig_10kb$V1 <- tot_gene_sig_10kb$V1 -1
#add 10000 to each bar so that it represents all pts within that chunk of the genome
tot_gene_sig_10kb$chunks <- tot_gene_sig_10kb$chunks +10000
print("tot_gene_sig_10kb")
write.csv(tot_gene_sig_10kb)

######################################################
# PLOT
######################################################
options(scipen=3)
tot_gene_subs_10kb$chunks <- tot_gene_subs_10kb$chunks / 1000000
tot_gene_sig_10kb$chunks <- tot_gene_sig_10kb$chunks / 1000000
subs_hist <- (ggplot(tot_gene_subs_10kb, aes(x = chunks, y = V1))
  + geom_histogram(stat = "identity", fill= "#2E294E", alpha = 0.9)
  + geom_histogram(data = tot_gene_sig_10kb, aes(x=chunks, y=V1),stat = "identity", fill= "#AA5042")
  + geom_point(data = hns_inver, aes(x = non_bidir_midpoint,y=fake2,colour = "#42474C"), size=10)
   + scale_fill_manual(values = c("#2E294E","#AA5042","#42474C" ), labels = c("H-NS Binding","Significant Inversion"))
  + theme(legend.position = "top")
#  geom_histogram(stat = "identity", fill= "#FE5F55")
  + labs(x = "Distance from the Origin of Replication (Mbp)", y = "Number of H-NS Binding Sites per 10Kbp")
##FOR PRESENTATION ONLY: NEXT LINE
#  + labs(x = "", y = "")
  + ggtitle("H-NS binding")
#   + geom_vline(xintercept = 0, colour = "red")
 + scale_x_continuous(labels = function(x) ifelse(x ==0, "0", x))
 + scale_y_continuous(labels = function(x) ifelse(x ==0,"0", x))
)
###########################
pdf("hns_hist.pdf")
subs_hist
#grid.newpage()
#grid.draw(rbind(ggplotGrob(gene_num_hist), ggplotGrob(exp_bar_top), size = "last"))
## next line adds border
#grid.rect(width = 0.99, height = 0.99, gp = gpar(lwd = 2, col = "black", fill = NA))
dev.off()
###########################



#################
###########
##get total number of sites in each 10kb section
###########
#list_freq <- lapply(list_dat_sets_just_subs, table)
#list_freq
#list_freq_sum <- lapply(list_freq, sum)
#list_freq_sum
#tot_num_sites <- as.data.frame(matrix(unlist(list_freq_sum), byrow=F))
#tot_num_sites
#length(tot_num_sites)
##remove last element if lengths do not match (for cbind later)
#if (nrow(tot_num_sites) == length(chunks_pos_NOzero)) {
#    tot_num_sites <- tot_num_sites
#} else {
#    tot_num_sites <- tot_num_sites[-nrow(tot_num_sites),]
#    tot_num_sites <- tot_num_sites[-1]
#}
#freq_dat <- as.data.frame(cbind(tot_num_sites,chunks_pos_NOzero))
#print("freq_dat")
#head(freq_dat)
#tail(freq_dat)
#write.csv(freq_dat,"hns_freq_dat.csv")
############
###weighted dataframe
############
#weighted_data <- cbind(tot_gene_subs_10kb,freq_dat)
#colnames(weighted_data) <- c("tot_gene_subs_10kb", "chunks_pos_NOzero", "tot_num_sites", "chunks_pos_NOzero")
#head(weighted_data)
#tail(weighted_data)
#length(weighted_data$tot_gene_subs_10kb)
#length(weighted_data$tot_num_sites)
#weighted_data$weighted_subs <- weighted_data$tot_gene_subs_10kb / weighted_data$tot_num_sites
#head(weighted_data)
#tail(weighted_data)
#
#weighted_data <- weighted_data[-2]
#head(weighted_data)
#tail(weighted_data)
#write.csv(weighted_data,"weighted_dat.csv")


#gather(x, value, G_HNS_binding:H3_HNS_binding)%>%
#group_by(x)%>%
#tally(value == 1)

##print("combine all the hns binding info into one df")
###except W3 bc it has different gene names
##hns_bind <- merge(sub_g_dat,sub_h_df, by = "gene_name")
##head(hns_bind)
#
#print ("merge HNS with K12 MG data")
#k12MG_df <- inver_dat_bidir %>% filter(strain == "K12MG") %>%
#            select(midpoint,rev_comp)
#colnames(k12MG_df) <- c("midpoint","class2")
#k12MG_df$class2 <- as.character(k12MG_df$class2)
#head(k12MG_df)
#summary(k12MG_df)
##inver_k12mg <- k12MG_df %>% filter(class == "block_avg_exp_invert")
##class2 <- rep("Inversion",length(inver_k12mg$class))
##inver_k12mg <- cbind(inver_k12mg,class2)
##inver_k12mg <- inver_k12mg %>% select(midpoint, class2)
##inver_k12mg <- k12MG_df %>% select(midpoint, class)
##change names in class col
#cor_inver_df <- k12MG_df
#inver_k12mg <- k12MG_df %>% 
#    mutate(class2 = recode(class2, 
#                      "0" = "No Inversion", 
#                      "1" = "Inversion"))
##lst <- c(1 = "Inversion", 0 = "No Inversion")
##inver_k12mg <- k12MG_df
##inver_k12mg$class2 <- as.character(lst[inver_k12mg$rev_comp])
##inver_k12mg <- inver_k12mg %>% select(midpoint,class2)
#print("inver_k12mg head")
#inver_k12mg <- unique(inver_k12mg)
#head(inver_k12mg)
#
#head(sub_g_dat)
#sub_g_dat$end <- as.numeric(as.character(sub_g_dat$end))
#sub_g_dat$start <- as.numeric(as.character(sub_g_dat$start))
#hns_dat <- within(sub_g_dat, midpoint <- (end + start) /2)
#colnames(hns_dat)[colnames(hns_dat) == "midpoint"] <- "midpoint"
#class2 <- rep("HNS_Binding",length(hns_dat$start))
#hns_dat <- cbind(hns_dat,class2)
#hns_cor_d <- hns_dat
#hns_dat <- hns_dat %>% select(midpoint, class2,gbk_strand)
#head(hns_dat)
#
##print("check higashi non coding overlap with data")
##head(sub_hnc_dat)
##sub_hnc_dat <- sub_hnc_dat %>% filter(HNS == "True")
##head(sub_hnc_dat)
##sub_hnc_dat$end <- as.numeric(as.character(sub_hnc_dat$end))
##sub_hnc_dat$start <- as.numeric(as.character(sub_hnc_dat$start))
##hns_dat <- within(sub_hnc_dat, midpoint <- (end + start) /2)
##colnames(hns_dat)[colnames(hns_dat) == "midpoint"] <- "midpoint"
##class2 <- rep("HNS_Binding",length(hns_dat$start))
##hns_dat <- cbind(hns_dat,class2)
##hns_cor_d <- hns_dat
##hns_dat <- hns_dat %>% select(start,end, class2)
##hns_dat <- hns_dat[complete.cases(hns_dat),]
##head(hns_dat)
##hns_cor_d <- hns_dat
##
##
##print("get new df with non-cod hns binding + inversions + sig block info")
##inver_cor_d <- inver_tmp_k12 %>% 
##            select(start,end,inversion)
##colnames(inver_cor_d) <- c("start1","end1","Inversion")
##inver_cor_d$Inversion <- as.character(inver_cor_d$Inversion)
##inver_cor_d <- unique(inver_cor_d)
##head(inver_cor_d)
##
##ir1 = with(inver_cor_d, IRanges(start1, end1))
##ir2 = with(hns_cor_d, IRanges(start, end))
##print("did ir1 and ir2")
##inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
##inver_cor_d[which(inver_cor_d$start1 == 383921),]
##colnames(inver_cor_d) <- c("start","end","Inversion","HNS_binding")
##inver_cor_d$HNS_binding <- as.integer(inver_cor_d$HNS_binding)
##inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
##head(inver_cor_d)
###inver_cor_d[which(inver_cor_d$start == 383921),]
#
#print("################################################################################")
#print("#ORIGIN SCALING AND BIDIRECTIONALITY HNS                                            ")
#print("################################################################################")
##first scaling things to the origin (if necessary)
#max_pos <- as.numeric(args[2])
#print("max_pos")
#max_pos
#oriC_pos <- as.numeric(args[3])
#print("oriC")
#oriC_pos
#terminus <- as.numeric(args[4])
#print("ter")
#terminus
#new_pos <- hns_dat$midpoint
#tmp_pos <- hns_dat$midpoint
#print("MIN POS")
#min(hns_dat$midpoint)
#
#if (bac_name == "E.coli" | replicon == "pSymA") {
#  to_shift_ter <- max_pos - oriC_pos
#  shifted_ter <-terminus + to_shift_ter
#  terminus <- shifted_ter
#}
#print("shifted ter")
#terminus
#
#if (replicon == "pSymB") {
#  shifted_ter <- terminus - oriC_pos
#  terminus <- shifted_ter
#}
#
#if (bac_name == "E.coli" | replicon == "pSymA" | replicon == "pSymB")
#{
#  for(i in 1:length(tmp_pos)) {
#    if (tmp_pos[i] >= oriC_pos) {
#      new_pos[i] <- tmp_pos[i] - oriC_pos
#    } else {
#      tmp_end <- max_pos - oriC_pos
#      new_pos[i] <- tmp_pos[i] + tmp_end
#    }
#  }
#  tmp_pos <- new_pos
#}
#
#
##now accounting for the bidirectionality. if things are between the start pos and
##the terminus then they will stay as the same position. If not, then they will be
##changed to a new position starting at 1 and going to the terminus
#new_pos2 <- tmp_pos
##also have to account for bidirectional replicaion in the strand, in
##the left replichore a complemented gene (1) is actually on the leading
##strand. so for the left replichore, all 0 -> 1, and 1 -> 0
#new_strand <- hns_dat$gbk_strand
#if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name ==
#"S.meliloti") {
#  for(i in 1:length(tmp_pos)) {
#    #left replichore
#    if (tmp_pos[i] > terminus) {
#      new_pos2[i] <- max_pos - tmp_pos[i]
#      # making sure the strand column accounts for bidirectional rep
#      if (hns_dat$gbk_strand[i] == 0) {
#        new_strand[i] <- 1
#      } else {
#        new_strand[i] <- 0
#      }
#    } else {
#    }
#  }
#  tmp_pos <- new_pos2
#
#  print("max tmp_pos")
#  max(tmp_pos)
#}
#
#
#if (bac_name == "Streptomyces") {
#  for(i in 1:length(tmp_pos)) {
#    # right replichore
#    if (tmp_pos[i] >= oriC_pos) {
#      new_pos[i] <- tmp_pos[i] - oriC_pos
#    }#if btwn origin and end of genome
#    # left replichore
#    if (tmp_pos[i] <= oriC_pos) {
#      new_pos[i] <- -1 * (oriC_pos - tmp_pos[i])
#      # making sure the strand column accounts for bidirectional rep
#      if (hns_dat$gbk_strand[i] == 0) {
#        new_strand[i] <- 1
#      } else {
#        new_strand[i] <- 0
#      }
#    }#if btwn origin and beginning of genome
#    if (tmp_pos[i] == oriC_pos) {
#      new_pos[i] <- 0
#    }#if equal to origin
#  }
#  tmp_pos <- new_pos
#}
#
#
#hns_dat$gbk_strand <- new_strand
#
#hns_dat$midpoint <- tmp_pos
#head(hns_dat)
#max(hns_dat$midpoint)
#min(hns_dat$midpoint)
#hns_dat <- hns_dat %>% select(midpoint,class2)
#
#hns_inver <- rbind(inver_k12mg,hns_dat)
#head(hns_inver)
#tail(hns_inver)
#
#print("combine HNS binary info to inversion df")
#print("THIS IS NOT BIDIRECTIONAL BC IT IS THE START AND ENDS AND NOT THE MIDPOINT!")
#head(inver_dat_bidir)
#inver_cor_d <- inver_dat_bidir %>% filter(strain == "K12MG") %>%
#            select(start,end,inversion)
##            select(start,end)
##colnames(inver_cor_d) <- c("start1","end1")
#colnames(inver_cor_d) <- c("start1","end1","Inversion")
#inver_cor_d$Inversion <- as.character(inver_cor_d$Inversion)
#inver_cor_d <- unique(inver_cor_d)
#head(inver_cor_d)
#write.table(inver_cor_d, 'inver_cor_data.csv', sep = "\t")
#
#print("hns")
#head(hns_cor_d)
##hns_cor_d <- hns_cor_d %>% select(start,end,class2)
#hns_cor_d <- hns_cor_d %>% select(start,end)
#summary(hns_cor_d)
##hns_cor_d <- hns_cor_d %>% 
##    mutate(class2 = recode(class2, 
##                      "No_HNS_Binding" = "0", 
##                      "HNS_Binding" = "1"))
#hns_cor_d <- hns_cor_d[order(hns_cor_d$start),]
#head(hns_cor_d)
#write.table(hns_cor_d, 'hns_cor_data.csv', sep = "\t")
#
#
#ir1 = with(inver_cor_d, IRanges(start1, end1))
#ir2 = with(hns_cor_d, IRanges(start, end))
#inver_cor_d$overlap = countOverlaps(ir1, ir2) != 0
#inver_cor_d[which(inver_cor_d$start1 == 383921),]
#colnames(inver_cor_d) <- c("start","end","Inversion","HNS_binding")
#inver_cor_d$HNS_binding <- as.integer(inver_cor_d$HNS_binding)
#inver_cor_d$Inversion <- as.integer(inver_cor_d$Inversion)
#head(inver_cor_d)
#inver_cor_d[which(inver_cor_d$start == 383921),]
#
#print("################################################################################")
#print("correlation test btwn inversion/non-inversion and hns binding")
#print("################################################################################")
#cor.test(inver_cor_d$Inversion, inver_cor_d$HNS_binding,
#         method = "pearson")
#
#print("get new df with hns binding + inversions + sig block info")
#cor_dat <- cor_dat %>% filter(strain == "K12MG")
#cor_dat <- cor_dat %>%
#    mutate(sig = if_else(pvalue <= 0.05, 1, 0))
#head(cor_dat)
#ir1 = with(cor_dat, IRanges(start, end))
#cor_dat$overlap = countOverlaps(ir1, ir2) != 0
#names(cor_dat)[names(cor_dat)=="overlap"] <- "HNS_binding"
#cor_dat[which(cor_dat$start == 383921),]
#cor_dat$HNS_binding <- as.integer(cor_dat$HNS_binding)
#cor_dat$sig <- as.integer(cor_dat$sig)
#head(cor_dat)
#
#
#print("################################################################################")
#print("correlation test btwn sig inversion/non-inversion and hns binding ALL BLOCKS")
#print("################################################################################")
#cor.test(cor_dat$sig, cor_dat$HNS_binding,
#         method = "pearson")
#
#
#print("################################################################################")
#print("correlation test btwn sig inversion/non-inversion and hns binding ONLY INVERSION BLOCKS")
#print("################################################################################")
##cor_dat <- cor_dat %>% filter(inversion == 1)
#cor.test(cor_dat$inversion, cor_dat$HNS_binding,
#         method = "pearson")
#
#
#print("test plot of HNS binding and inversions")
#pdf("hns_inver_plot_test.pdf")
#ggplot(hns_inver, aes(x=midpoint, y=class2, color=class2))+
##  geom_jitter(aes(tt, val), data = df, colour = I("red"), 
##               position = position_jitter(width = 0.05)) +
##  geom_point(size = 3) +
#  geom_point(size = 2, alpha=0.4)
##  geom_errorbar(aes(ymin=val-sd, ymax=val+sd), width = 0.01, size = 1)
#dev.off()
#
#
#
#
##print("###############################################################################")
##print("INVERSION VISUALIZATION PARALLEL SETS")
##print("###############################################################################")
##print("read in block info file")
##file_name <- as.character(args[9])
##block_inf <- read.table(file_name,sep = "\t", header = FALSE)
##colnames(block_inf) <- c("block","strain","start","end","rev_comp","inversion")
##print("make column of midpoint of each block")
##block_inf <- within(block_inf, midpoint <- (end + start) /2)
##colnames(block_inf)[colnames(block_inf) == "midpoint"] <- "midpoint"
##head(block_inf)
##bi_dat <- block_inf %>% select(block,strain,midpoint)
##bi_dat <-  spread(bi_dat, strain, midpoint)
##head(bi_dat)
##
##print("test parallel sets plot")
##ps <- bi_dat %>%
##  gather_set_data(2:3) %>%
##head(ps)
###  ggplot(aes(x, id = id, split = y, value = 1))  +
###  geom_parallel_sets(aes(fill = engine)) 
###pdf("parallel_sets.pdf")
###ps
###dev.off()
##
##
