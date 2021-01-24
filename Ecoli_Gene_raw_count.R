pacman::p_load(tidyverse, data.table, ggplot2, magrittr)
options(scipen = 999)

setwd("./Gene_expression_data/") 

files <- list.files()

#ATCC_GSE94978 simpliying df-------------------
#no duplicates; different gene id than the rest of the data files
ATCC_GSE94978 <- fread(files[1], fill=TRUE) 

ATCC_GSE94978 <- ATCC_GSE94978 %>% select(c(gene_id, gene_name, ATCC2_count, ATCC3_count)) %>%
  mutate(strain ='ATCC') %>%
  select(gene_id, Locus_tag=gene_name, GSE94978_1=ATCC2_count,
         GSE94978_2=ATCC3_count, strain)



#BW25113 simplying df----------------
#BW25113_GSE73673
BW25113_GSE73673_06 <- fread(files[2]) #4299
BW25113_GSE73673_07 <- fread(files[3]) #4299
BW25113_GSE73673_08 <- fread(files[4]) #4299

BW25113_GSE73673 <- list(BW25113_GSE73673_06,
                         BW25113_GSE73673_07,
                         BW25113_GSE73673_08) %>%
  reduce(left_join, by="V1")

remove(BW25113_GSE73673_06, BW25113_GSE73673_07, BW25113_GSE73673_08) #12897

BW25113_GSE73673 <- BW25113_GSE73673 %>%
  select(Locus_tag=V1, GSE73673_6=V2.x, GSE73673_7=V2.y, GSE73673_8=V2) #4299


#BW25113_GSE85914
BW25113_GSE85914 <- fread(files[5])
BW25113_GSE85914 <- BW25113_GSE85914 %>% filter(!(substring(Name, 1, 1) == "-")) %>% #remove predicted DNA (uncertain variants)
  dplyr::select(gene_id=Name,Locus_tag=Synonym, GSE85914=`Expression WT`) #4318

BW25113 <- inner_join(BW25113_GSE73673, BW25113_GSE85914, by=c('Locus_tag')) #4299




#K12_DH10B  simplying df-------------
#K12_DH10B_GSE98890
#4362-4237 = 125 dups nrow(K12_DH10B_GSE98890)-length(unique(K12_DH10B_GSE98890$gene_name))
K12_DH10B_GSE98890 <- fread(files[6], fill = TRUE)

K12_DH10B_GSE98890 <- K12_DH10B_GSE98890 %>% group_by(gene_name) %>%
  mutate(temp=length(unlist(strsplit(gene_name,"_"))),
         locus_tag= paste(unlist(strsplit(gene_name,"_"))[temp-1], 
                          unlist(strsplit(gene_name, "_"))[temp], sep="_"),
         gene_id= ifelse(unlist(strsplit(gene_name,"_"))[2] == "locusTag", 
                          unlist(strsplit(gene_name, "_"))[1],locus_tag),
         gene_id = ifelse(is.na(gene_id), gsub("_", "", locus_tag), gene_id)) %>%
  ungroup() %>% #all unique locus tag
  select(gene_id, Locus_tag=locus_tag, GSE98890_1=Control_rep1,
         GSE98890_2=Control_rep2)


#K12_MG1655  simplying df----------
#K12_MG1655_GSE114917 
K12_MG1655_GSE114917 <- fread(files[7]) %>%
  mutate(exp='GSE114917') %>%
  select(gene_id=GeneName, GSE114917_LP1= LP1, GSE114917_LP2=LP2,
         GSE114917_LP3=LP3, GSE114917_LP4=LP4) #4452
 

#K12_MG1655_GSE40313
K12_MG1655_GSE40313_1 <- fread(files[8]) 
K12_MG1655_GSE40313_2 <- fread(files[9]) 

K12_MG1655_GSE40313 <- left_join(K12_MG1655_GSE40313_1, 
                                 K12_MG1655_GSE40313_2, by='V1') %>%
  select(Locus_tag=V1, GSE40313_1=V2.x, GSE40313_2=V2.y) #4245

remove(K12_MG1655_GSE40313_1, K12_MG1655_GSE40313_2)


##K12_MG1655_GSE54199 
#the column names are the same for both DNA and RNA files

#DNA
K12_MG1655_GSE54199_DNA <- fread(files[10], fill = TRUE) 
colnames(K12_MG1655_GSE54199_DNA) <- as.character(K12_MG1655_GSE54199_DNA[1,])

K12_MG1655_GSE54199_DNA <- K12_MG1655_GSE54199_DNA[-1,] %>% 
  group_by(Locus_tag) %>% 
  mutate(GSE54199_DNA=mean(c(as.numeric(Raw_counts_replicate_1), 
                      as.numeric(Raw_counts_replicate_2)))) %>% ungroup() %>% 
  select(gene_id=Gene, Locus_tag, GSE54199_DNA)
#4373


#RNA
K12_MG1655_GSE54199_RNA <- fread(files[11]) %>% 
  group_by(Locus_tag) %>% 
  mutate(GSE54199_RNA=mean(c(Raw_counts_replicate_1, 
                      Raw_counts_replicate_2))) %>% ungroup() %>% 
  select(gene_id=Gene, Locus_tag, GSE54199_RNA)
#4368

K12_MG1655_GSE54199 <- left_join(K12_MG1655_GSE54199_RNA, K12_MG1655_GSE54199_DNA, 
                                 by=c('Locus_tag', 'gene_id')) 
#4368
remove(K12_MG1655_GSE54199_DNA, K12_MG1655_GSE54199_RNA)


#K12_MG1655_GSE60522
K12_MG1655_GSE60522 <- fread(files[12]) %>% filter(!(gene == "")) %>%
  select(gene_id=gene, GSE60522_1=`WT untr 1`, GSE60522_2=`WT untr 2`,
         GSE60522_3=`WT untr 3`)
#4497



#Merging Order as label for count.x
#1. K12_MG1655_GSE60522_count (Alp, gene_id)
#2. K12_MG1655_GSE114917_count (Alp, gene_id)
#3. K12_MG1655_GSE54199 (Alp, num, gene_id, Name)
#4. K12_MG1655_GSE40313 

K12MG1655 <- merge(K12_MG1655_GSE60522, 
                       K12_MG1655_GSE114917, by="gene_id") %>% #4376
  merge(K12_MG1655_GSE54199, by="gene_id") %>% #4235
  merge(K12_MG1655_GSE40313, by = "Locus_tag") #3891



## adding info using the gene_info files-------
ATCC_info <- fread("../Genomes/Ecoli_ATCC_25922_NZ_CP009072_gene_info.txt")
BW_info <- fread("../Genomes/Ecoli_BW25113_NZ_CP009273_gene_info.txt")
H7_info <- fread("../Genomes/Ecoli_0157H7_chrom_Sakai_BA000007_gene_info.txt")
K12_DH_info <- fread("../Genomes/Ecoli_K12_DH10B_NC_010473_gene_info.txt")
K12_MG_info <- fread("../Genomes/Ecoli_K12_MG1655_chrom_U00096_gene_info.txt")


##ATCC info--------------------------
ATCC_GSE94978_final <- ATCC_GSE94978 %>% merge(ATCC_info, by.x="gene_id", #5180
                                               by.y="gbk_locus_tag") %>%
  select(gene_id, gbk_start, gbk_end, gbk_midpoint, gbk_gene_id, gbk_old_locus_tag,
         gbk_strand, GSE94978_1,GSE94978_2) #4825
  
#nrow(ATCC_GSE94978_final[is.na(ATCC_GSE94978_final$gbk_start)==1,]) 
#355 not found in gnk -> (5180-355)*2==9680


#Mapping gene names to BW from K12 --------------
gene_map <- fread("../BW25113_and_K12_MG1655_gene_map.txt") #4201
gene_map %<>% dplyr::select(-c(V1, V3,V5))
colnames(gene_map) <- c('gbk_locus_tag', 'mapped_name')

BW_info_mapped <- merge(BW_info, gene_map, by='gbk_locus_tag', all.x=T)

#4110 -> 4077 (lose 33 genes)
BW_final <- merge(BW25113, BW_info_mapped, by.x='Locus_tag', by.y='mapped_name')



##K12_DH_info-----------------

K12_DH10B_final <- merge(K12_DH10B_GSE98890, K12_DH_info, by.x="Locus_tag", 
                                               by.y="gbk_old_locus_tag") 
#115 not found in gbk
#4247


##K12_MG_info-----------
K12MG_final <- merge(K12MG1655, K12_MG_info, by.x="Locus_tag", by.y="gbk_locus_tag", all.x=T) 
#3891



##Write to dataframe 
fwrite(ATCC_GSE94978_final, "../Final_expression_051820/Ecoli_ATCC_25922_final_RawCounts.csv", row.names = F)
fwrite(BW_final, "../Final_expression_051820/Ecoli_BW25113_final_RawCounts.csv", row.names = F) 
fwrite(K12_DH10B_final, "../Final_expression_051820/Ecoli_K12_DH10B_final_RawCounts.csv", row.names = F)
fwrite(K12MG_final, "../Final_expression_051820/Ecoli_K12_MG1655_final_RawCounts.csv", row.names = F)


