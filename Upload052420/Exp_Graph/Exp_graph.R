pacman::p_load(tidyverse, data.table, edgeR, magrittr, dplyr, ggplot2)
options(scipen = 999)

setwd("/Users/quee/Backup/Gene_expression_data/")
files <- list.files()


##BW_25113------------------------

#BW25113_GSE85913
BW25113_GSE73673_06 <- fread(files[2])
BW25113_GSE73673_07 <- fread(files[3])
BW25113_GSE73673_08 <- fread(files[4])
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
  mutate(BW25113_GSE73673 = median(c(norm1, norm2, norm3))) %>% ungroup() %>%
  dplyr::select(gene_id, BW25113_GSE73673) %>% arrange()



#BW25113_GSE85914

BW25113_GSE85914 <- fread(files[5])
BW25113_GSE85914 %<>% filter(!(substring(Name, 1, 1) == "-")) 
#remove predicted DNA (uncertain variants)

exp <- DGEList(count=BW25113_GSE85914[,3], group=rep(1, 1))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm')

BW25113_GSE85914 %<>%  
  cbind(norm) %>% 
  group_by(gene_id = Synonym) %>%
  mutate(BW25113_GSE85914 = norm) %>%
  dplyr::select(gene_id, BW25113_GSE85914) %>% arrange()


#merge 2 files
BW25113 <- BW25113_GSE73673 %>% merge(BW25113_GSE85914, by= "gene_id") 
BW25113 %<>% group_by(gene_id) %>% mutate(norm = median(c(BW25113_GSE73673, BW25113_GSE85914), na.rm = T),
         GSE73673 = log10(abs(BW25113_GSE73673-norm)),
         GSE85914 = log10(abs(BW25113_GSE85914-norm)))

BW25113 %<>% gather(GSE73673, GSE85914, key="Dataset", value = "exp")


##Graph BW 
# ggplot(BW25113, aes(x=gene_id, y=exp, color=Dataset)) + geom_point() +
#   labs(title=expression(~italic(Ecoli)~' BW25113'),
#        x="Gene id", y = "log10(|expression - median expression|)")


##K12MG----------------------

K12_MG1655_GSE114917 <- fread(files[7])
#normalize
exp <- DGEList(count=K12_MG1655_GSE114917[,2:5], group=rep(1, 4))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2', 'norm3', 'norm4')

K12_MG1655_GSE114917 %<>% cbind(norm) %>% 
  group_by(gene_id = GeneName) %>% #Alphabet
  summarise(K12_MG1655_GSE114917 = median(c(norm1, norm2, norm3, norm4))) 


#K12_MG1655_GSE40313
K12_MG1655_GSE40313_1 <- fread(files[8])
K12_MG1655_GSE40313_2 <- fread(files[9])
K12_MG1655_GSE40313 <- list(K12_MG1655_GSE40313_1, 
                            K12_MG1655_GSE40313_2) %>%
  reduce(left_join, by="V1") 

exp <- DGEList(count=K12_MG1655_GSE40313[,2:3], group=rep(1, 2))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2')

K12_MG1655_GSE40313 %<>% cbind(norm) %>%
  group_by(Locus_tag = V1) %>% #bxxxx
  summarise(K12_MG1655_GSE40313 = median(c(norm1, norm2)))

rm(K12_MG1655_GSE40313_1, K12_MG1655_GSE40313_2)

##K12_MG1655_GSE54199 
#the column names are the same for both DNA and RNA files

#DNA
K12_MG1655_GSE54199_DNA <- fread(files[10], fill = TRUE) 
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
  dplyr::select(Locus_tag, Translation_start, Translation_stop, Strand, 
                gene_id=Gene, RPKM, norm_count, raw_count, K12_MG1655_GSE54199_DNA)

#RNA
K12_MG1655_GSE54199_RNA <- fread(files[11]) %>% 
  mutate(K12_MG1655_GSE54199_RNA = 1) 

exp <- DGEList(count=K12_MG1655_GSE54199_RNA[,c(6,7)], group=rep(1, 2))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2')

K12_MG1655_GSE54199_RNA %<>% cbind(norm) %>%
  group_by(Locus_tag) %>% #bxxxx
  mutate(norm_count = median(c(norm1, norm2)),
         raw_count = median(c(Raw_counts_replicate_1, Raw_counts_replicate_2)),
         K12_MG1655_GSE40313_RNA = 1) %>% ungroup() %>%
  dplyr::select(Locus_tag, Translation_start, Translation_stop, Strand, 
                gene_id=Gene, RPKM, norm_count, raw_count, K12_MG1655_GSE54199_RNA) 

#merging
K12_MG1655_GSE54199 <- merge(K12_MG1655_GSE54199_DNA, 
                             K12_MG1655_GSE54199_RNA, 
                             by="Locus_tag") 

K12_MG1655_GSE54199 %<>% group_by(gene_id.x, Locus_tag) %>%
  summarise(gene_id = gene_id.x,
            K12_MG1655_GSE54199 = median(c(norm_count.x, norm_count.y))) %>% 
  ungroup() %>% dplyr::select(-gene_id.x)



#K12_MG1655_GSE60522
K12_MG1655_GSE60522 <- fread(files[12])
K12_MG1655_GSE60522 %<>% filter(!(gene == ""))

exp <- DGEList(count=K12_MG1655_GSE60522[,2:4], group=rep(1, 3))
norm <- as.data.frame(cpm(calcNormFactors(exp), normalized.lib.sizes = FALSE))
names(norm) <- c('norm1', 'norm2', 'norm3')

K12_MG1655_GSE60522 %<>% cbind(norm) %>%
  group_by(gene_id = gene) %>% #Alphabet
  summarise(K12_MG1655_GSE60522 = median(c(norm1, norm2, norm3))) 




#Merging Order as label for count.x
#1. K12_MG1655_GSE60522_count (Alp, gene_id)
#2. K12_MG1655_GSE114917_count (Alp, gene_id)
#3. K12_MG1655_GSE54199 (Alp, num, gene_id, Name)
#4. K12_MG1655_GSE40313 

K12MG1655 <- merge(K12_MG1655_GSE60522, 
                   K12_MG1655_GSE114917, by="gene_id", all = TRUE) %>%#by alpabet gene id #4573
  merge(K12_MG1655_GSE54199, by="gene_id", all = TRUE) %>% #4589 - 16 from sum(!(unique(K12_MG1655_GSE54199$gene_id) %in% unique(K12MG1655$gene_id)))
  merge(K12_MG1655_GSE40313, by = "Locus_tag", all = TRUE) %>% #4814
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


K12MG1655_graph <- K12MG1655[complete.cases(K12MG1655),]
K12MG1655_graph %<>% group_by(unique_id) %>% 
  mutate(norm = median(c(K12_MG1655_GSE60522,
                         K12_MG1655_GSE114917,
                         K12_MG1655_GSE54199,
                         K12_MG1655_GSE40313)),
         GSE60522 = log10(abs(K12_MG1655_GSE60522-norm)),
         GSE114917 = log10(abs(K12_MG1655_GSE114917-norm)),
         GSE54199 = log10(abs(K12_MG1655_GSE54199-norm)),
         GSE40313 = log10(abs(K12_MG1655_GSE40313-norm))) %>% 
  gather(GSE60522,
         GSE114917,
         GSE54199,
         GSE40313, key="Dataset", value = "exp")

K12plot <- ggplot(data = K12MG1655_graph, 
       mapping = aes(x=unique_id, y=exp,color=Dataset)) + 
  geom_point(shape=1) +
  labs(title=expression(~italic(Ecoli)~' K12_MG1655'),
       x="Gene id", y = "log10(|expression - median expression|)") 
ggsave("../Exp_Graph/K12MG_graph.pdf", K12plot)



