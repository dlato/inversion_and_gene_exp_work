# COG analysis for gene expression and inversions:
# COGs in significant blocks (blocks with sig diff in exp btwn homologous inverted and non-inverted genes via permutation test)

library(ggplot2)
library(dplyr)

dat <- read.table("C:/Users/synch/Documents/PhD/inversion_and_gene_exp_work/Stats/inversions_gene_exp_sig_K12MG_COGs.csv")
colnames(dat) <- c("geneID", "COGID")
head(dat)
dat$COGID <- as.factor(dat$COGID)
cog_def <- read.delim("C:/Users/synch/Documents/cog.def.tab", sep = "\t")
colnames(cog_def) <- c("COGID", "COGcat", "disc", "geneName","disc2", "disc3", "disc4")
head(cog_def)
summary(cog_def)
#add col that has the COG broad category as a col
cog_def$COGID <- as.factor(cog_def$COGID)
cog_def$COGcat <- as.factor(cog_def$COGcat)
cog_def %>% filter(COGID == "COG0156")
cog_def %>% filter(COGID == "COG0001")
dat %>% filter(COGID == "COG0156")
summary(cog_def)
dat$COGcat <- cog_def$COGcat[match(dat$COGID,cog_def$COGID)]
head(dat)
g <- (ggplot(dat, aes(x = COGcat)) # Plot with values on top
  + geom_bar(stat = "count") 
  #geom_text(aes(label = Freq), vjust = 0)
)
g
