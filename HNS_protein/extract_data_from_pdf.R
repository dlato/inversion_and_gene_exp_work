#this will extract the data in tables from Oshima 2006 and Lang 2007 PDFs with H-NS binding data
library(tabulizer)
library(dplyr)


tt4 <- tabulizer::extract_areas("raw_data_files/Oshima_2006_HNS_binding_sites_and_HGT.pdf")[[1]]
