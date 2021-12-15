# Testing gene expression difference between different media

library(DESeq2)
library(dplyr)
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
##############
# SET SEED
##############
set.seed(925)

rawexp <- read.csv("new_dataframes/raw_exp.csv", header = TRUE)
#select only BW25113 cols
bw_exp <- rawexp %>% select(BW_names, BW25113_GSE73673_6,BW25113_GSE73673_7,BW25113_GSE73673_8, BW25113_GSE85914)
#deal with non-unique gene names
bw_exp <- bw_exp[!duplicated(bw_exp[ , c("BW_names")]), ]
#make gene names row names
rownames(bw_exp) <- bw_exp$BW_names
bw_exp <- bw_exp %>% select(-BW_names)
head(bw_exp)
#create sample info
samples <- data.frame (sample  = c("BW25113_GSE73673_6","BW25113_GSE73673_7","BW25113_GSE73673_8", "BW25113_GSE85914"),
                       experiment  = c("GSE73673","GSE73673","GSE73673", "GSE85914"),
                       media = c("LB","LB","LB","M9")
)
rownames(samples) <- samples$sample
head(samples)

#DESeq analysis
deseq2Data <- DESeqDataSetFromMatrix(countData=bw_exp,
                                     colData=samples,
                                     design= ~ media)
deseq_results <- DESeq(deseq2Data)
print("done deseq")
#alpha cutoff of 0.05
deseq_results2 <- results(deseq_results, alpha=0.05)
print("summary")
summary(deseq_results2)

pdf("DESeq_BWA25113_media.pdf")
plotMA(deseq_results2)
dev.off()
