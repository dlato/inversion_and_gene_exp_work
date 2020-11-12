#this will extract the data in tables from Oshima 2006 and Lang 2007 PDFs with H-NS binding data
library(tabulizer)
library(dplyr)
library(readr) ## parse_number()


#the number at the end of this next line indicates the single page to extract the information from
#tt4 <- tabulizer::extract_areas("HNS_Protein/raw_data_files/Oshima_2006_HNS_binding_sites_and_HGT_supptab2.pdf")[[2]]
o_tab <- extract_tables("HNS_Protein/raw_data_files/Oshima_2006_HNS_binding_sites_and_HGT_supptab2.pdf")
final_o_tab <- do.call(rbind, o_tab[-length(o_tab)])

m <- matrix(parse_number(c(tt4)), ncol=ncol(tt4))
