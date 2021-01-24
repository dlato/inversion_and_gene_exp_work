#library(tidyverse)
library(dplyr)
library(tidyr)
library("stringr")
#library(broom)
#library(rstatix)
#library(ggpubr)
library(DESeq2)
#library(scales) # needed for oob parameter
#library(viridis)
##library(tximport)
#library(readr)
#library(reshape2)
#library(data.table)
#library("RColorBrewer")
#library("gplots")
#library("ggplot2")
#library("GenomicRanges")
print("#############                 ")
print("## DESeq expression comparison")
print("#############                 ")
raw_file <- "raw_expression_DESeq.csv"
raw_dat <- read.csv(raw_file, header = TRUE)
#ensure that all raw values are integers
raw_deseqW <- raw_dat %>%
                mutate_at(seq(2,19,1), as.integer)
head(raw_deseqW)
sample_file <- "deseq_ExpDesign2.csv"
sampleData <- read.csv(sample_file, header = TRUE)
sampleData
s_cols <- row.names(sampleData)
df.sub <- raw_deseqW[, s_cols]
df.sub
print("Put the columns of the count data in the same order as rows
names of the sample info, then make sure it worked")
raw_deseqW <- raw_deseqW[,unique(rownames(sampleData))]
all(colnames(raw_deseqW) == rownames(sampleData))
print("Order the treatments so that it is sensible: non-inversion (control) -> inversion")
sampleData$treatment <- factor(sampleData$treatment, levels=c("0", "1"))
print("#############                 ")
print("## DESeq analysis accounting for experiment effect")
print("#############                 ")
print("Create the DEseq2DataSet object")
print("DE btwn strains")
deseq2Data <- DESeqDataSetFromMatrix(countData=raw_deseqW, colData=sampleData, design= ~ expID + treatment)
deseq_results <- DESeq(deseq2Data)
print("done deseq")
deseq_results2 <- results(deseq_results, alpha=0.05)
print("summary")
summary(deseq_results2)


