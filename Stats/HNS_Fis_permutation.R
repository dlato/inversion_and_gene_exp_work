#HNS and Fis permutation test

library(regioneR)
library(dplyr)
library(stringr)
set.seed(369)
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

#read in inversion data
inver <- read.table("./Stats/inver_cor_data.csv")
sig_inver <- read.table("./Stats/sig_inversion_dat.csv")
#make proper bed file format
chrom <- rep("chr1", length(nrow(inver)))
inver <- cbind(chrom,inver)
head(inver)
chrom <- rep("chr1", length(nrow(sig_inver)))
sig_inver <- cbind(chrom,sig_inver)
head(sig_inver)

#create universe (all of my alignment data)
uni <- inver %>%
       select(c("chrom","start1","end1"))
head(uni)
summary(uni)

#create df with inversion regions only
inver_a <- inver %>%
       filter(Inversion == 1) %>%
       select(c("chrom","start1","end1"))
head(inver_a)
summary(inver_a)

#create df with sig inversion regions only
inver_s <- sig_inver %>%
       filter(sig == 1) %>%
       select(c("chrom","gbk_start","gbk_end"))
head(inver_s)
summary(inver_s)

#convert to proper object
hnsR_ALL <- toGRanges(hns_ALL)
hnsR_g <- toGRanges(hns_g)
hnsR_h1 <- toGRanges(hns_h1)
hnsR_h2 <- toGRanges(hns_h2)
hnsR_h3 <- toGRanges(hns_h3)
hnsR_h4 <- toGRanges(hns_h4)
hnsR_h5 <- toGRanges(hns_h5)
hnsR_h6 <- toGRanges(hns_h6)
hnsR_l <- toGRanges(hns_g)
hnsR_o <- toGRanges(hns_l)
hnsR_u <- toGRanges(hns_u)
uniR <- toGRanges(uni)
inverR_a <- toGRanges(inver_a)
inverR_s <- toGRanges(inver_s)
head(hnsR_g)
head(uniR)
head(inverR_a)
head(inverR_s)

#permutation tests
#hns binding with inversions
iterations <- 1000

print("HNS_ALL")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_ALL, verbose=FALSE)
summary(pt)

print("HNS_G")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_g, verbose=FALSE)
summary(pt)

print("HNS_H1")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h1, verbose=FALSE)
summary(pt)
print("HNS_H2")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h2, verbose=FALSE)
summary(pt)
print("HNS_H3")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h3, verbose=FALSE)
summary(pt)
print("HNS_H4")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h4, verbose=FALSE)
summary(pt)
print("HNS_H5")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h5, verbose=FALSE)
summary(pt)
print("HNS_H6")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h6, verbose=FALSE)
summary(pt)

print("HNS_L")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_l, verbose=FALSE)
summary(pt)

print("HNS_O")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_o, verbose=FALSE)
summary(pt)

print("HNS_U")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_u, verbose=FALSE)
summary(pt)

#permutation tests with sig inversions
print("HNS_ALL")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_ALL, verbose=FALSE)
summary(pt)

print("HNS_G")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_g, verbose=FALSE)
summary(pt)

print("HNS_H1")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h1, verbose=FALSE)
summary(pt)
print("HNS_H2")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h2, verbose=FALSE)
summary(pt)
print("HNS_H3")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h3, verbose=FALSE)
summary(pt)
print("HNS_H4")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h4, verbose=FALSE)
summary(pt)
print("HNS_H5")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h5, verbose=FALSE)
summary(pt)
print("HNS_H6")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_h6, verbose=FALSE)
summary(pt)

print("HNS_L")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_l, verbose=FALSE)
summary(pt)

print("HNS_O")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_o, verbose=FALSE)
summary(pt)

print("HNS_U")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=hnsR_u, verbose=FALSE)
summary(pt)


#############
# Fis
#############


#read in data
g_dat <- read.csv("./fis_binding/Grainger_2006_fis_binding.csv")
sub_g_dat <- g_dat %>%
        select(Gene_name,Coordinates)
colnames(sub_g_dat) <- c("gene_name","Coords")
coords <- str_split_fixed(sub_g_dat$Coords, "-", 2)
colnames(coords) <- c("start","end")
sub_g_dat <- cbind(sub_g_dat,coords)
sub_g_dat <- sub_g_dat %>%
        select(gene_name,start,end) %>%
        filter(gene_name != "Between convergent ORFs") %>%
        filter(start != "N/A")

head(sub_g_dat)

#make proper bed file format
chrom <- rep("chr1", length(nrow(sub_g_dat)))
fis_d <- cbind(chrom,sub_g_dat)
head(fis_d)
summary(fis_d)

fis_d <- fis_d %>%
         select(c("chrom","start","end"))
fis_d$start <- as.numeric(as.character(fis_d$start))
fis_d$end <- as.numeric(as.character(fis_d$end))
head(fis_d)
summary(fis_d)

#convert to proper object
fisR <- toGRanges(fis_d)
head(fisR)

#permutation tests
#fis and inversions (all)
print("Fis_G")
pt <- permTest(A=inverR_a, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=fisR, verbose=FALSE)
summary(pt)


#fis and sig inversions
print("Fis_G")
pt <- permTest(A=inverR_s, ntimes=iterations, randomize.function=resampleRegions, universe=uniR, evaluate.function=numOverlaps, B=fisR, verbose=FALSE)
summary(pt)
