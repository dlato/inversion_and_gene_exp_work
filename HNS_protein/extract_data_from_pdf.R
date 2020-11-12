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
t7 <- as.data.frame(o_tab[[55]])
t7 <- t7 %>% select(-V2) %>%
  separate(V3, c("lengthnt", "start", "end"), " ")
t7 <- t7 %>%
  separate(V1, c("lt", "gene"), " ")
t7 <- cbind(t7[,1:7], test=t7[,8], t3="",hns="", t2=t7[,ncol(t7)])
colnames(t7) <- col_names_tab
final_o_tab <- rbind(final_o_tab, t7)



final_o_tab <- do.call(rbind, tt4[-length(tt4)])



m <- matrix(parse_number(c(tt4)), ncol=ncol(tt4))
