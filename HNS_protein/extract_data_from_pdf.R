#this will extract the data in tables from Oshima 2006 and Lang 2007 PDFs with H-NS binding data
library(tabulizer)
library(dplyr)
library(tidyr)
library(readr) ## parse_number()


#the number at the end of this next line indicates the single page to extract the information from
#tt4 <- tabulizer::extract_areas("HNS_Protein/raw_data_files/Oshima_2006_HNS_binding_sites_and_HGT_supptab2.pdf")[[2]]
#below will extract information from ALL pages
#tt4 <- tabulizer::extract_areas("HNS_Protein/raw_data_files/Oshima_2006_HNS_binding_sites_and_HGT_supptab2.pdf")


o_tab <- extract_tables("HNS_Protein/raw_data_files/Oshima_2006_HNS_binding_sites_and_HGT_supptab2.pdf")
head(o_tab[[1]])
head(o_tab[[2]])
length(o_tab[[2]])
el1 <- o_tab[[2]]
ncol(el1)
length(o_tab)

#column names for final table
col_names_tab <- c("locus_tag", "gene", "lengthnt","start","end", "strand", "product","Nakamura","Lawrence","hns","conservation")
#getting total number of columns for each pdf page
cols_pages <- unlist(lapply(o_tab, function(x) ncol(x)))
#only pages with the correct total 11 columns
full_pages <- which(cols_pages == 11)
#add full pages to final table
final_o_tab <- as.data.frame(do.call(rbind, o_tab[full_pages]))
colnames(final_o_tab) <- col_names_tab

#dealing with pages with 7 columns
t7 <- as.data.frame(o_tab[[67]])
t7 <- t7 %>%
       separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- cbind(t7[,1:7], Nakamura="", test=t7[,8],hns="", t2=t7[,ncol(t7)])
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

t7 <- as.data.frame(o_tab[[93]])
t7 <- t7 %>% select(-V2) %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- t7 %>%
  separate(V1, c("lt", "gene"), " ")
t7 <- cbind(t7[,1:7], test=t7[,8], t3="",hns="", t2=t7[,ncol(t7)])
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

t7 <- as.data.frame(o_tab[[95]])
t7 <- t7 %>% select(-V2) %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- t7 %>%
  separate(V1, c("lt", "gene"), " ")
t7 <- cbind(t7[,1:7], test=t7[,8], t3="",hns="", t2=t7[,ncol(t7)])
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

#dealing with pages with 8 columns
which(cols_pages == 8)
n_h <- c(55,57,61,82,83,92,97)
t7 <- as.data.frame(do.call(rbind,o_tab[n_h]))
t7 <- t7 %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- cbind(t7[,1:7], test=t7[,8], t3="", t2=t7[,9:ncol(t7)])
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

n_h_s1 <- c(90,103)
t7 <- as.data.frame(do.call(rbind,o_tab[n_h_s1]))
t7 <- t7 %>% select(-V2) %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- t7 %>%
  separate(V1, c("lt", "gene"), " ")
t7 <- cbind(t7[,1:7], test=t7[,8], t3="", t2=t7[,9:ncol(t7)])
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)


n_l <- c(60)
t7 <- as.data.frame(do.call(rbind,o_tab[n_l]))
t7 <- t7 %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- cbind(t7[,1:7], test=t7[,8], t2=t7[,9],t3="", t=t7[,ncol(t7)] )
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)


n_l_s1 <- c(99)
t7 <- as.data.frame(do.call(rbind,o_tab[n_l_s1]))
t7 <- t7 %>% select(-V2) %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- t7 %>%
  separate(V1, c("lt", "gene"), " ")
t7 <- cbind(t7[,1:7], test=t7[,8], t2=t7[,9],t3="", t=t7[,ncol(t7)] )
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

#dealing with pages with 9 columns
which(cols_pages == 9)

s3 <- c(27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,43,44,45,46,47,48,50,51,52,53,56,62,63,64,65,68,69,70,72,74,75,
        76,79,80,81,84,85,86,88,89,91,96,100,101,104,105,107)
t7 <- as.data.frame(do.call(rbind,o_tab[s3]))
t7 <- t7 %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

s1_s3 <- c(25,26,42,49,54,58,59,66,71,73,77,78,87,94,98,102,106,108)
t7 <- as.data.frame(do.call(rbind,o_tab[s1_s3]))
t7 <- t7 %>% select(-V2) %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- t7 %>%
  separate(V1, c("lt", "gene"), " ")
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

n <- c(11)
t7 <- as.data.frame(do.call(rbind,o_tab[n]))
t7 <- cbind(t7[,1:8],t3="",t4="", t=t7[,ncol(t7)] )
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

n_s1 <- c(6)
t7 <- as.data.frame(do.call(rbind,o_tab[n_s1]))
t7 <- t7 %>% select(-V2) %>%
  separate(V1, c("lt", "gene"), " ")
t7 <- cbind(t7[,1:8],t3="",t4="", t=t7[,ncol(t7)] )
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

n_h_s3 <- c(55)
t7 <- as.data.frame(do.call(rbind,o_tab[n_h_s3]))
t7 <- t7 %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- cbind(t7[,1:8],t3="",t=t7[,9:ncol(t7)] )
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

#dealing with pages with 10 columns
which(cols_pages == 10)
p <- c(1)
t7 <- as.data.frame(do.call(rbind,o_tab[p]))
t7 <- t7 %>%
  separate(V4, c("start", "end"), " ")
colnames(t7) <- col_names_tab
t7 <- t7[-1,]
final_o_tab <- rbind(final_o_tab, t7)

n_h_s1 <- c(19)
t7 <- as.data.frame(do.call(rbind,o_tab[n_h_s1]))
t7 <- t7 %>% select(-V2) %>%
  separate(V1, c("lt", "gene"), " ")
t7 <- cbind(t7[,1:8],t3="",t=t7[,9:ncol(t7)] )
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

n_h <- c(3)
t7 <- as.data.frame(do.call(rbind,o_tab[n_h]))
t7 <- cbind(t7[,1:8],t3="",t=t7[,9:ncol(t7)] )
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)

length(final_o_tab$locus_tag)

write.csv(final_o_tab, "HNS_Protein/raw_data_files/Oshima_2006_HGT_HNS.csv")

