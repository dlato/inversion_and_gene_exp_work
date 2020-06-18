
#load library
library(data.table)
library(edgeR)
library(magrittr)
library(dplyr)
library(ggplot2)

options(scipen = 999)

setwd("/home/dlato/Queenie_2019_2020/Upload052420/Gene_expression_data")

getwd()

##just housekeeping
files <- list.files()

#ATCC_GSE94978 simpliying df-------------------
#no duplicates; different gene id than the rest of the data files
ATCC_GSE94978 <- read.delim(files[1])

exp <- DGEList(counts=ATCC_GSE94978[,9:10], group = rep(1,2))
norm <- cpm(calcNormFactors(exp), normalized.lib.sizes=FALSE) 
norm <- as.data.frame(norm) 
colnames(norm) <- c("norm1", "norm2")

ATCC_GSE94978 %<>% cbind(norm) %>% group_by(gene_id) %>%
  mutate(count = median(c(norm1, norm2)),
         ATCC_GSE94978 = 1) %>% select(-c(V13))


#BW25113 simplying df----------------
#BW25113_GSE73673
BW25113_GSE73673_06 <- read.delim(files[2])
BW25113_GSE73673_07 <- read.delim(files[3])
BW25113_GSE73673_08 <- read.delim(files[4])
BW25113_GSE73673 <- list(BW25113_GSE73673_06,
                         BW25113_GSE73673_07,
                         BW25113_GSE73673_08) %>% 
  reduce(left_join, by="V1")

remove(BW25113_GSE73673_06, BW25113_GSE73673_07, BW25113_GSE73673_08)

exp <- DGEList(count=BW25113_GSE73673[,2:4], group=rep(1, 3))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2', 'norm3')


BW25113_GSE73673 %<>%
  cbind(norm) %>% 
  group_by(gene_id = V1) %>%
  mutate(norm_count = median(c(norm1, norm2, norm3)),
          raw_count = median(c(V2.x, V2.y, V2)),
            BW25113_GSE73673 = 1) %>% ungroup() %>%
  select(gene_id, norm_count, raw_count, BW25113_GSE73673) %>% arrange()





#BW25113_GSE85914

BW25113_GSE85914 <- read.delim(files[5])
BW25113_GSE85914 %<>% filter(!(substring(Name, 1, 1) == "-")) 
#remove predicted DNA (uncertain variants)

exp <- DGEList(count=BW25113_GSE85914[,3], group=rep(1, 1))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm')

BW25113_GSE85914 %<>% 
  cbind(norm) %>% 
  group_by(gene_id = Synonym) %>%
  mutate(norm_count = norm, 
         raw_count = `Expression WT`, BW25113_GSE85914 = 1) %>%
  select(Name, gene_id, norm_count, raw_count, BW25113_GSE85914) %>% arrange()


#merge 2 files
BW25113 <- BW25113_GSE73673 %>% merge(BW25113_GSE85914, by= "gene_id", all = TRUE) 




##This below block of code was to compare another normalizing method  

# #normalizing common (because DGEList cant handle NA) raw counts for comparison 
# BW25113_raw <- BW25113[complete.cases(BW25113),]
# exp <- DGEList(count=BW25113_raw[,c(3,7)], group=rep(1, 2))
# norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
# names(norm) <- c('raw_norm1', 'raw_norm2')
# BW25113_raw %<>% cbind(norm) %>% 
#   group_by(gene_id, Name) %>% 
#   summarise(raw_norm = median(c(raw_norm1, raw_norm2)),
#          norm = median(c(norm_count.x, norm_count.y)))
# 
# #plot to compare 2 normalizing method
# plot(BW25113_raw$norm, BW25113_raw$raw_norm)
# abline(coef=c(0,1),col='red')


BW25113 %<>%
  group_by(Locus_tag=gene_id) %>%
  mutate(norm_count = median(c(norm_count.x, norm_count.y), na.rm = TRUE),
         raw_count = median(c(raw_count.x, raw_count.y), na.rm = TRUE)) %>%
  ungroup() %>%
  select(gene_id=Name, Locus_tag, norm_count, raw_count, 
                BW25113_GSE73673, BW25113_GSE85914) %>%
  arrange()
#All unique Name and gene_id 





#K12_DH10B  simplying df-------------
#K12_DH10B_GSE98890
#4362-4237 = 125 dups nrow(K12_DH10B_GSE98890)-length(unique(K12_DH10B_GSE98890$gene_name))

K12_DH10B_GSE98890 %<>% group_by(gene_name) %>%
  mutate(temp=length(unlist(strsplit(gene_name,"_"))),
         locus_tag= paste(unlist(strsplit(gene_name,"_"))[temp-1], 
                          unlist(strsplit(gene_name, "_"))[temp], sep="_"),
         gene_id= ifelse(unlist(strsplit(gene_name,"_"))[2] == "locusTag", 
                          unlist(strsplit(gene_name, "_"))[1],locus_tag),
         gene_id = ifelse(is.na(gene_id), gsub("_", "", locus_tag), gene_id)) %>%
  ungroup() %>%
  select(gene_id, locus_tag, Control_rep1, Control_rep2)
#all unique now 

#normalize
exp <- DGEList(count=K12_DH10B_GSE98890[,3:4], group=rep(1, 2))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2')

K12_DH10B_GSE98890 %<>% cbind(norm) %>%
  group_by(gene_id, locus_tag) %>%
  summarise(norm_count = median(c(norm1, norm2)),
            raw_count = median(c(Control_rep1, Control_rep2)), 
            K12_DH10B_GSE98890 = 1) 



##This below block of code was to compare another normalizing method 
# #plot to compare normalizing methods
# plot(K12_DH10B_GSE98890$norm_count, K12_DH10B_GSE98890$raw_count)
# abline(coef=c(0,1),col='red')


#K12_MG1655  simplying df----------
#K12_MG1655_GSE114917 
K12_MG1655_GSE114917 <- read.delim(files[7])
#normalize
exp <- DGEList(count=K12_MG1655_GSE114917[,2:5], group=rep(1, 4))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2', 'norm3', 'norm4')

K12_MG1655_GSE114917 %<>% cbind(norm) %>% 
  group_by(gene_id = GeneName) %>% #Alphabet
  summarise(norm_count = median(c(norm1, norm2, norm3, norm4)),
            raw_count = median(c(LP1, LP2, LP3, LP4)),
            K12_MG1655_GSE114917 = 1) 


#K12_MG1655_GSE40313
K12_MG1655_GSE40313_1 <- read.delim(files[8])
K12_MG1655_GSE40313_2 <- read.delim(files[9])
K12_MG1655_GSE40313 <- list(K12_MG1655_GSE40313_1, 
                            K12_MG1655_GSE40313_2) %>%
  reduce(left_join, by="V1") 

exp <- DGEList(count=K12_MG1655_GSE40313[,2:3], group=rep(1, 2))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2')

K12_MG1655_GSE40313 %<>% cbind(norm) %>%
  group_by(Locus_tag = V1) %>% #bxxxx
  summarise(norm_count = median(c(norm1, norm2)),
            raw_count = median(c(V2.x, V2.y)),
            K12_MG1655_GSE40313 = 1)


##K12_MG1655_GSE54199 
#the column names are the same for both DNA and RNA files

#DNA
K12_MG1655_GSE54199_DNA <- read.delim(files[10]) 
colnames(K12_MG1655_GSE54199_DNA) <- as.character(K12_MG1655_GSE54199_DNA[1,])
K12_MG1655_GSE54199_DNA <- K12_MG1655_GSE54199_DNA[-1,] %>% 
  mutate(K12_MG1655_GSE54199_DNA = 1, 
         Raw_counts_replicate_1=as.numeric(Raw_counts_replicate_1),
         Raw_counts_replicate_2=as.numeric(Raw_counts_replicate_2)) 

exp <- DGEList(count=K12_MG1655_GSE54199_DNA[,c(6,7)], group=rep(1, 2))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2')

K12_MG1655_GSE54199_DNA %<>% cbind(norm) %>%
  group_by(Locus_tag) %>% #bxxxx
  mutate(norm_count = median(c(norm1, norm2)),
            raw_count = median(c(Raw_counts_replicate_1, Raw_counts_replicate_2)),
            K12_MG1655_GSE40313_DNA = 1) %>% ungroup() %>%
  select(Locus_tag, Translation_start, Translation_stop, Strand, 
         gene_id=Gene, RPKM, norm_count, raw_count, K12_MG1655_GSE54199_DNA)

#RNA
K12_MG1655_GSE54199_RNA <- read.delim(files[11]) %>% 
  mutate(K12_MG1655_GSE54199_RNA = 1) 

exp <- DGEList(count=K12_MG1655_GSE54199_RNA[,c(6,7)], group=rep(1, 2))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2')

K12_MG1655_GSE54199_RNA %<>% cbind(norm) %>%
  group_by(Locus_tag) %>% #bxxxx
  mutate(norm_count = median(c(norm1, norm2)),
            raw_count = median(c(Raw_counts_replicate_1, Raw_counts_replicate_2)),
            K12_MG1655_GSE40313_RNA = 1) %>% ungroup() %>%
  select(Locus_tag, Translation_start, Translation_stop, Strand, 
         gene_id=Gene, RPKM, norm_count, raw_count, K12_MG1655_GSE54199_RNA) 

#merging
K12_MG1655_GSE54199 <- merge(K12_MG1655_GSE54199_DNA, 
                                 K12_MG1655_GSE54199_RNA, 
                             by="Locus_tag", all.x = TRUE) 

K12_MG1655_GSE54199 %<>% group_by(gene_id.x, Locus_tag) %>%
  mutate(norm = median(c(norm_count.x, norm_count.y)),
         raw_norm = median(c(raw_count.x, raw_count.y))) %>% ungroup() %>%
  select(Locus_tag, gene_id = gene_id.x, Translation_start = Translation_start.x,
                Translation_stop = Translation_stop.x, Strand = Strand.x, RPKM = RPKM.x,
                norm, raw_norm, K12_MG1655_GSE54199_DNA, K12_MG1655_GSE54199_RNA) %>%
  arrange()



#K12_MG1655_GSE60522
K12_MG1655_GSE60522 <- read.delim(files[12])
K12_MG1655_GSE60522 %<>% filter(!(gene == ""))

exp <- DGEList(count=K12_MG1655_GSE60522[,2:4], group=rep(1, 3))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2', 'norm3')

K12_MG1655_GSE60522 %<>% cbind(norm) %>%
  group_by(gene_id = gene) %>% #Alphabet
  mutate(norm_count = median(c(norm1, norm2, norm3)),
         raw_count = median(c(`WT untr 1`, `WT untr 2`, `WT untr 3`)),
         K12_MG1655_GSE60522 = 1) %>%
  select(gene_id, norm_count, raw_count, K12_MG1655_GSE60522)




#Merging Order as label for count.x
#1. K12_MG1655_GSE60522_count (Alp, gene_id)
#2. K12_MG1655_GSE114917_count (Alp, gene_id)
#3. K12_MG1655_GSE54199 (Alp, num, gene_id, Name)
#4. K12_MG1655_GSE40313 

K12MG1655 <- merge(K12_MG1655_GSE60522, 
                       K12_MG1655_GSE114917, by="gene_id", all = TRUE) %>%
  group_by(gene_id) %>% #by alpabet gene id 
  mutate(norm1 = norm_count.x, norm2 = norm_count.y,
         raw1 = raw_count.x, raw2 = raw_count.y) %>% 
  select(-c(raw_count.x, raw_count.y, norm_count.x, norm_count.y)) %>% #4573
  merge(K12_MG1655_GSE54199, by="gene_id", all = TRUE) %>% #4589 - 16 from sum(!(unique(K12_MG1655_GSE54199$gene_id) %in% unique(K12MG1655$gene_id)))
  mutate(norm3 = norm, raw3 = raw_norm) %>% select(-c(norm, raw_norm)) %>%
  merge(K12_MG1655_GSE40313, by = "Locus_tag", all = TRUE) %>% 
  mutate(norm4 = norm_count, raw4 = raw_count) %>% 
  select(-c(norm_count, raw_count)) %>% #4814
  mutate(unique_id = paste0(gene_id, Locus_tag))


##This below block of code was to compare another normalizing method 
# K12MG1655_raw <- K12MG1655[complete.cases(K12MG1655),]
# exp <- DGEList(count=K12MG1655_raw[,c(7,8,16,19)], group=rep(1, 4))
# norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
# names(norm) <- c('norm_1', 'norm_2', 'norm_3', 'norm_4')
# 
# K12MG1655_raw %<>% cbind(norm) %>% group_by(Locus_tag, gene_id) %>%
#   summarise(raw_norm = median(c(norm_1, norm_2, norm_3, norm_4)),
#             norm = median(c(norm1, norm2, norm3, norm4)))
# 
# plot(K12MG1655_raw$norm, K12MG1655_raw$raw_norm)
# abline(coef=c(0,1),col='red')


K12MG1655 %<>% group_by(unique_id) %>% 
  mutate(norm = median(c(norm1, norm2, norm3, norm4), na.rm = TRUE),
         raw = median(c(raw1, raw2, raw3, raw4))) %>%
  ungroup() %>%
  select(Locus_tag, gene_id, Strand, RPKM, Translation_start, Translation_stop,
                raw, norm,K12_MG1655_GSE60522, K12_MG1655_GSE114917, 
                K12_MG1655_GSE54199_DNA, K12_MG1655_GSE54199_RNA, 
                K12_MG1655_GSE40313) %>% arrange()







#O157_H7_GSE46120  simplying df--------------
H7_GSE46120_36 <- read.delim(files[13], header = FALSE) %>% 
  mutate(V9 = as.character(V9), gene_id = substring(V9, 9, nchar(V9))) %>% 
  group_by(gene_id, V1, V2, V3, V7, V8) %>%
  summarise(V4=min(V4), V5= max(V5), raw = sum(V6)) 

exp <- DGEList(count=H7_GSE46120_36[,9], group=rep(1, 1))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1')

H7_GSE46120_36 %<>% as.data.frame() %>% cbind(norm) %>% 
  mutate(H7_GSE46120_36 = 1)


H7_GSE46120_37 <- read.delim(files[14], header = FALSE) %>% 
  mutate(V9 = as.character(V9), gene_id = substring(V9, 9, nchar(V9))) %>%
  group_by(gene_id, V1, V2, V3, V7, V8) %>%
  summarise(V4=min(V4), V5= max(V5), raw = sum(V6)) 

exp <- DGEList(count=H7_GSE46120_37[,9], group=rep(1, 1))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm2')

H7_GSE46120_37 %<>% as.data.frame() %>% cbind(norm) %>% 
  mutate(H7_GSE46120_37 = 1)


#merge
H7_GSE46120 <- merge(H7_GSE46120_36, H7_GSE46120_37, by="gene_id", all = TRUE) 


##This below block of code was to compare another normalizing method 
# H7_GSE46120_raw <- H7_GSE46120[complete.cases(H7_GSE46120),] 
# exp <- DGEList(count=H7_GSE46120_raw[,c(9, 18)], group=rep(1, 2))
# norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
# names(norm) <- c('raw_norm1','raw_norm2')
# H7_GSE46120_raw %<>% cbind(norm) %>%
#   group_by(gene_id) %>%
#   mutate(raw_norm=median(c(raw_norm1, raw_norm2)), 
#          norm=median(c(norm1, norm2)))
# 
# plot(H7_GSE46120_raw$norm, H7_GSE46120_raw$raw_norm)
# abline(coef=c(0,1),col='red')

H7_GSE46120 %<>% group_by(gene_id) %>%
  mutate(norm_count = median(c(norm1, norm2), na.rm = TRUE),
         raw_count = median(c(raw.x, raw.y), na.rm = TRUE)) %>% ungroup() %>%
  mutate(V1=coalesce(V1.x, V1.y),
        V2=coalesce(V2.x, V2.y),
        V3=coalesce(V3.x, V3.y),
        V4=coalesce(V4.x, V4.y),
        V5=coalesce(V5.x, V5.y),
        V7=coalesce(V7.x, V7.y)) %>% 
  select(gene_id, V1, V2, V3, V4, V5, V7, raw_count, norm_count,H7_GSE46120_36,H7_GSE46120_37)
  






# %>%
#   filter(is.na(str_extract(gene_id, "ECs"))) %>%
#   filter(is.na(str_extract(gene_id, "rrn"))) %>%
#   filter(is.na(str_extract(gene_id, "INT"))) %>%
#   filter(is.na(str_extract(gene_id, "no matches"))) #33 rows


## adding info using the gene_info files-------
ATCC_info <- read.delim("../Genomes/Ecoli_ATCC_25922_NZ_CP009072_gene_info.txt")
#ATCC_info_old <- read.delim("../Genomes/Ecoli_ATCC_25922_NZ_CP009072_13May20_gene_info.txt")
BW_info <- read.delim("../Genomes/Ecoli_BW25113_NZ_CP009273_gene_info.txt")
#BW_info_old <- read.delim("../Genomes/Ecoli_BW25113_NZ_CP009273_13May20_gene_info.txt")
H7_info <- read.delim("../Genomes/Ecoli_0157H7_chrom_Sakai_BA000007_gene_info.txt")
#H7_info_old <- read.delim("../Genomes/Ecoli_0157H7_chrom_Sakai_NC_002695_3Mar20_gene_info.txt")
K12_DH_info <- read.delim("../Genomes/Ecoli_K12_DH10B_NC_010473_gene_info.txt")
#K12_DH_info_old <- read.delim("../Genomes/Ecoli_K12_DH10B_NC_010473_13May20_gene_info.txt")
K12_MG_info <- read.delim("../Genomes/Ecoli_K12_MG1655_chrom_U00096_gene_info.txt")
#K12_MG_info_old <- read.delim("../Genomes/Ecoli_K12_MG1655_chrom_U00096_13May20_gene_info.txt")


##ATCC info--------------------------
# nrow(ATCC_GSE94978) #5180
# nrow(ATCC_info) #4939 old #5272 new 
ATCC_GSE94978_final <- ATCC_GSE94978 %>% merge(ATCC_info, by.x="gene_id",
                                               by.y="gbk_locus_tag") %>%
  select(gene_id, gbk_start, gbk_end, gbk_midpoint, gbk_gene_id, gbk_old_locus_tag,
                norm_exp=count)


#nrow(ATCC_GSE94978_final[is.na(ATCC_GSE94978_final$gbk_start)==1,]) 
#355 not found in gnk

##BW info -----------------

# nrow(BW25113) #4507
# nrow(BW_info) #4464 #4775
# sum(is.na(BW_info$gbk_old_locus_tag)) #403old #572
# sum(!is.na(BW_info$gbk_gene_id)) #174old #188 there're only 188 gene id 
# sum(is.na(BW25113$gene_id)) #189
# sum(is.na(BW25113$Locus_tag)) #0 



BW25113 %<>% mutate(key=substring(Locus_tag,2,5))
BW_info %<>% mutate(key=substring(gbk_old_locus_tag,9,12))

##*** Checking Inconsistencies between gene names between BW and K12_MG********

# BW25113_test <- BW25113 %>% merge(BW_info, by="key", all.x=T)
# View(BW25113_test[is.na(BW25113_test$gbk_start),]) #307 inconsistenices
# 
# BW25113_test <- BW25113_test[,-c(8:14)]
# BW25113_test %<>% merge(BW_info, by.x = "gene_id", by.y = "gbk_gene_id", all.x = T)
# 
# 
# #Adding new bxxxx gene id info from K12 
# BW25113_names <- BW25113 %>% select(gene_id, Locus_tag) #4507
# K12MG1655_names <- K12MG1655 %>% select(gene_id, Locus_tag) #4898
# check <- merge(BW25113_names, K12MG1655_names, by="Locus_tag")#4409
# check2 <- check[complete.cases(check),] #4203 
# check3 <- anti_join(check, check2) #206, K12 doesnt have the locus tag 
# diff <- check2[!check2$gene_id.x==check2$gene_id.y,] #54 diff gene id with same locus tag
# 
# #check if BW_info 
# BW_info_names <- BW_info %>% select(key, gbk_gene_id) %>%
#   mutate(Locus_tag=paste0("b",key)) %>% filter(!is.na(gbk_gene_id)) %>% 
#   filter(!is.na(key)) %>% 
#   merge(BW25113_names, by="Locus_tag")
# View(BW_info_names[BW_info_names$gbk_gene_id!=BW_info_names$gene_id,]) #31 inconsistencies 



#*****************************
#check if K12_info same as BW df 
# check4 <- merge(BW25113_names, K12_MG_info, by.x="Locus_tag", by.y="gbk_locus_tag") %>%
#   select(Locus_tag, gene_id, gbk_gene_id) 
# check4 <- check4[complete.cases(check4),]
# View(check4[check4$gene_id!=check4$gbk_gene_id,]) #only 19 inconsistencies o/f 4318

  
BW25113_comp <- BW25113[complete.cases(BW25113[, c(5:6)]),]
# nrow(BW25113_comp) #4110
# 
# sum(is.na(BW25113_comp$gene_id)) #0
# sum(is.na(BW25113_comp$Locus_tag)) #0 

BW25113_comp %<>% merge(BW_info, by="key", all.x = T) #merged by 
#sum(is.na(BW25113_comp$gbk_start)) #63 NAs
BW25113_final <- BW25113_comp[is.na(BW25113_comp$gbk_start),] %>% select(-c(gbk_start:gbk_strand)) %>%
  merge(BW_info, by.x = "gene_id", by.y = "gbk_gene_id", all.x = T) %>% mutate(key=key.x,
                                                                               gbk_gene_id=gene_id) %>%
  select(-c(key.y,key.x)) 
BW25113_final <- BW25113_final[, c(13, 1:9, 14, 10:12)]
#sum(is.na(BW25113_final2$gbk_start))#61 



BW25113_final %<>% rbind(BW25113_comp[!is.na(BW25113_comp$gbk_start),])

#****
BW25113_final <- BW25113_final[complete.cases(BW25113_final$gbk_strand),] 
BW25113_final %<>% select(gene_id, Locus_tag, gbk_start, gbk_end, gbk_midpoint, gbk_gene_id, gbk_locus_tag,
                norm_exp=norm_count)
#nrow(BW25113_final)#4049

##K12_DH_info-----------------

##Checks 
# nrow(K12_DH_info) #4257old #4537
# sum(is.na(K12_DH_info$gbk_gene_id)) #alphas 1630NAs | 1875NAs
# sum(is.na(K12_DH_info$gbk_old_locus_tag)) #200 NAs | 285NAs
# 
# nrow(K12_DH10B_GSE98890) #4362
# sum(is.na(K12_DH10B_GSE98890$gene_id)) #alpha no NAs
# sum(is.na(K12_DH10B_GSE98890$locus_tag)) #ECDH10B no NAs


K12_DH10B_final <- K12_DH10B_GSE98890 %>% merge(K12_DH_info, by.x="locus_tag", 
                                               by.y="gbk_old_locus_tag") %>%
                                               by.y="gbk_old_locus_tag") %>%
  select(locus_tag, gene_id, gbk_start, gbk_end, gbk_midpoint, gbk_gene_id, gbk_locus_tag,
                norm_exp=norm_count)
#115 not found in gbk
#4247


##K12_MG_info-----------
##Checks 
# nrow(K12_MG_info) #4141 old #4498 locus_tag->bxxxx, gene_id->alpha 
# sum(!is.na(K12_MG_info$gbk_locus_tag)) #no NAs  
# sum(!is.na(K12_MG_info$gbk_gene_id)) #no NAs  
# 
# 
# nrow(K12MG1655) #4898 locus_tag->bxxxx, gene_id->alpha
# sum(is.na(K12MG1655$Locus_tag)) #262 NAs
# sum(is.na(K12MG1655$gene_id)) #263NAs

K12MG_final <- K12MG1655[complete.cases(K12MG1655[,c(9:13)]),] #3891

K12MG_final %<>% merge(K12_MG_info, by.x="Locus_tag", by.y="gbk_locus_tag", all.x=T) %>%
  select(Locus_tag, gene_id, gbk_start, gbk_end, gbk_midpoint, gbk_gene_id, gbk_old_locus_tag,
                norm_exp=norm)
#3891


##Write to dataframe 
fwrite(ATCC_GSE94978_final, "../Final_expression_051820/Ecoli_ATCC_25922_final_exp.csv", row.names = F)
fwrite(BW25113_final, "../Final_expression_051820/Ecoli_BW25113_final_exp.csv", row.names = F)
fwrite(K12_DH10B_final, "../Final_expression_051820/Ecoli_K12_DH10B_final_exp.csv", row.names = F)
fwrite(K12MG_final, "../Final_expression_051820/Ecoli_K12_MG1655_final_exp.csv", row.names = F)

##compare 
BW25113_list <- BW25113_final$gene_id #4049
K12MG_list <- K12MG_final$gene_id #3891
sum(BW25113_list %in% K12MG_list) #3838

##H7_info--------------------
##Checks 
# nrow(H7_info) #5491
# nrow(H7_info_old) #5330
# nrow(H7_GSE46120) #1039

H7_GSE46120_comp <- H7_GSE46120[complete.cases(H7_GSE46120[,10:11]),] #39 

H7_GSE46120 %<>% mutate(gene_id_new=if_else(substring(gene_id, nchar(gene_id),nchar(gene_id))==";",
                                        substring(gene_id, 1, nchar(gene_id)-1),
                                        gene_id))

print("code finished running")


# colnames(K12_MG_info) <- c("Start", "Stop", "Midpoint", "gene_id", "Locus_tag", "Strand")
# K12MG1655_final <- K12MG1655 %>% merge(K12_MG_info, by="Locus_tag", all.x = T) 
# temp <- K12MG1655_final %>% filter(is.na(Start)==1) %>% 
#   select(-c(Start, Stop, Midpoint, gene_id.y, Strand.y)) %>%
#   merge(K12_MG_info, by.x="gene_id.x", by.y="gene_id", all.x=T)
# temp2 <- temp %>% filter(is.na(Start)==1)   



## Compare gene ids for different strains-------------
# the final dataframes are:
# ATCC_GSE94978
# BW25113 #4507
# K12_DH10B_GSE98890 #4362
# K12MG1655 #4898
# H7_GSE46120 #1039

# length(intersect(BW25113[["gene_id"]], K12_DH10B_GSE98890[["gene_id"]]))
# #3584 from 4507, 4362
# View(K12_DH10B_GSE98890[BW25113[["gene_id"]] %in% K12_DH10B_GSE98890[["gene_id"]],])
# 
# 
# length(intersect(K12MG1655[["gene_id"]], K12_DH10B_GSE98890[["gene_id"]]))
# #3735 from 4898, 4362
# View(K12_DH10B_GSE98890[K12MG1655$gene_id %in% K12_DH10B_GSE98890$gene_id,])
# 
# length(intersect(BW25113[["gene_id"]], K12MG1655[["gene_id"]]))
# #4277 from 4507, 4898
# View(BW25113[K12MG1655[["gene_id"]] %in% BW25113[["gene_id"]],])




