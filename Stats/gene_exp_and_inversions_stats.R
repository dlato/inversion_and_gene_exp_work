#test script for inversions and gene expression analysis
#library(tidyverse)
library(dplyr)
#library(tidyr)
#library(broom)
#library(rstatix)
#library(ggpubr)
#library(DESeq2)
library(reshape2)
library("RColorBrewer")
library("gplots")
library("ggplot2")
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
#                  legend.title = element_blank(),
                  legend.text = element_text(size = 18),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facest
                  panel.spacing.x=unit(0, "lines"),
#                  legend.key = element_blank(),
#                  legend.background=element_blank(),
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
gei_dat <- within(gei_dat, midpoint <- (end - start) /2)
colnames(gei_dat)[colnames(gei_dat) == "midpoint"] <- "midpoint"
head(gei_dat)
bac_name <- as.character(args[5])
replicon <- as.character(args[6])

options("scipen"=100, "digits"=10)
print("################################################################################")
print("#ORIGIN SCALING AND BIDIRECTIONALITY                                            ")
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


##now accounting for the bidirectionality. if things are between the start pos and
##the terminus then they will stay as the same position. If not, then they will be
##changed to a new position starting at 1 and going to the terminus
#new_pos2 <- tmp_pos
##also have to account for bidirectional replicaion in the strand, in
##the left replichore a complemented gene (1) is actually on the leading
##strand. so for the left replichore, all 0 -> 1, and 1 -> 0
#new_strand <- gei_dat$strand
#if (bac_name == "E.coli" | bac_name == "B.subtilis" | bac_name ==
#"S.meliloti") {
#  for(i in 1:length(tmp_pos)) {
#    #left replichore
#    if (tmp_pos[i] > terminus) {
#      new_pos2[i] <- max_pos - tmp_pos[i]
#      # making sure the strand column accounts for bidirectional rep
#      if (gei_dat$strand[i] == 0) {
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
#      if (gei_dat$strand[i] == 0) {
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
#gei_dat$strand <- new_strand
#
gei_dat$gbk_midpoint <- tmp_pos
#gei_dat <- as.data.frame(cbind(gei_dat$block, gei_dat$gene, gei_dat$sec, tmp_pos, gei_dat$dS, gei_dat$dN, gei_dat$omega, gei_dat$sec_len))
#colnames(gei_dat) <- c("block","end","gbk_end","gbk_gene_id","gbk_locus_tag","gbk_midpoint","gbk_old_locus_tag","gbk_start","gbk_strand","gene_id","inversion","locus_tag","norm_exp","rev_comp","start","strain","taxa","gene_name")
head(gei_dat)
max(gei_dat$gbk_midpoint)
min(gei_dat$gbk_midpoint)

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
print("means")
print(i_len)
print(ni_len)
    uniq_avg_exp_invert <- c(uniq_avg_exp_invert,i_avg)
    uniq_avg_exp_noninvert <- c(uniq_avg_exp_noninvert,ni_avg)
    uniq_avg_len_invert <- c(uniq_avg_len_invert,i_len)
    uniq_avg_len_noninvert <- c(uniq_avg_len_noninvert,ni_len)
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

print("SAVED DATA TO FILE")
write.table(gei_dat, 'inversions_gene_exp_wtest_data.csv', sep = "\t")

print("make df with just block info")
block_df <- subset(gei_dat,select = c("block","start","end","midpoint","gbk_midpoint","block_w_pvalue","block_t_pvalue","block_t_stat","block_t_confinf1","block_t_confinf2", "block_avg_exp_invert", "block_avg_exp_noninvert","block_avg_len_invert", "block_avg_len_noninvert","strain"))
#block_df_uniq <- unique(block_df)
block_df_uniq <- block_df

#flip so ggplot can use
block_df_w <- melt(block_df_uniq,
        # ID variables - all the variables to keep but not split apart
        # on
    id.vars=c("block", "start", "end","midpoint","gbk_midpoint","block_t_stat","block_t_confinf1","block_t_confinf2","block_avg_exp_invert", "block_avg_exp_noninvert","block_avg_len_invert", "block_avg_len_noninvert", "strain"),
        # The source columns
    measure.vars=c("block_w_pvalue", "block_t_pvalue"),
        # Name of the destination column that will identify the
        # original
        # column that the measurement came from
    variable.name="test",
    value.name="pvalue"
)
print("non zero pvals") 
complete_block_df <- block_df_w[which(block_df_w$pvalue != "NA"),]
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


print("HEAD")
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
tmp2 <- subset(tmp,select = c("block","pvalue","block_avg_exp_invert", "block_avg_exp_noninvert","block_avg_len_invert", "block_avg_len_noninvert"))
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
up <- df[which(df$len_reg == "up"),]
length(unique(up$block))
print("percent of SIG blocks with inversions len > noninversion len")
(length(unique(up$block)) / length(unique(df$block))) *100
print("number of SIG blocks with inversion gene len < noninversions gene len")
up <- df[which(df$len_reg == "down"),]
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




#############
## DESeq expression comparison
#############
##set up df
##this will need to be changed once the real data comes in
##********************
#tmp_df <- gei_dat[,c("norm_exp","gene_id","strain","inversion")]
#DH <- tmp_df[which(tmp_df$strain == "K12DH"),]
#MG <- unique(tmp_df[which(tmp_df$strain == "K12MG"),])
#DH <- data.frame(DH[!duplicated(DH[ , c("gene_id")]),])
#levels(DH$strain) <- c(levels(DH$strain), "K12DH_i")
#DH$strain[DH$inversion == 1] <- "K12DH_i"
#MG <- data.frame(MG[!duplicated(MG[ , c("gene_id")]),])
#levels(MG$strain) <- c(levels(MG$strain), "K12MG_i")
#MG$strain[MG$inversion == 1] <- "K12MG_i"
#gene_names <- unique(MG$gene_id)
#DH <- unique(DH[which(DH$gene_id %in% gene_names),])
#exp_df <- rbind(DH,MG)
#exp_df <- spread(exp_df, strain, norm_exp)
#rownames(exp_df) <- exp_df$gene_id
#exp_df <- exp_df[,-1]
##inversions_col <- exp_df$inversion
#exp_df <- exp_df[,-1]
##********************
#t_count = 1
#t <- vector(mode='character', length = length(colnames(exp_df)))
#for (l in colnames(exp_df)) {
#  tmp_vec <- strsplit(l, "_")
#  tmp_l <- length(unlist(tmp_vec))
#  if (tmp_l > 1){
#    t[t_count] <- "inversion"
#  } else {
#    t[t_count] <- "none"
#  }
#  t_count <- t_count +1
#}
##deal with NA values
#exp_df <- exp_df %>% replace(is.na(.), 0)
#
#print(t)
## set up experimental design
#experimental_design = data.frame(
#  sample_names = colnames(exp_df),  # sample name
##  individual = factor(colnames(exp_df)), # each individual strain
#  treatment = factor(t)  # inversion = 1 or no inversion = 0
##  lane = factor(parse_names[,3])      # Which lane on the Illumina flowcell.
#)
#
#DESeq_data <- DESeqDataSetFromMatrix(exp_df, experimental_design, design = formula(~ treatment))
#


############
# ALL inversions
############



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
