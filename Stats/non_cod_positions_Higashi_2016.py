#!/usr/local/bin/python3
#
# Program will go through the inverted mafft blocks and grab the block
# number, taxa names, start position, end position, and if it is
# reversed.
#
# to run: ~/CODE/inversions_block_info.py *.mafft
#
# file columns:
# block_name	taxa_name	block_start	block_end	rev_comp	inversion	
#
# NOTE: 1 = inversion, and 1 = reversed complement
#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import Bio
import numpy as np
from Bio import SeqIO
#import pdb # FOR DEBUGGING ONLY import pdb
#import pandas as pd #importing pandas
#import numpy as np #importing Numpy
#import csv
import os #importing operating system commands
import glob #importing glob which helps read files
from Bio import AlignIO
import pandas as pd
#from functools import partial
#from collections import defaultdict
# FOR DEBUGGING ONLY pdb.set_trace()



#read in Higashi non-coding file
h_nc = sys.argv[1]
hnc_dat = pd.read_csv(h_nc)
hnc_dat.columns =['index','L_gene','R_gene','avg_dist','HNS','TtoT','known_promoter','TF_name']

#read in gene info file
inff = sys.argv[2]
ginfo = pd.read_csv(inff,sep="\t")


starts = []
ends = []
for index, row in hnc_dat.iterrows():
    print(hnc_dat.loc[index,"L_gene"])
    print(hnc_dat.iloc[index]["L_gene"])
    start = ginfo.gbk_start[ginfo['gbk_gene_id'] == hnc_dat.loc[index,"L_gene"]]
    if (len(start) > 0):
        start += 1
        print(ginfo.loc[ginfo['gbk_gene_id'] == hnc_dat.loc[index,"L_gene"],]) 
        print(start.to_string(index=False))
        start = start.to_string(index=False)
        end = ginfo.loc[ginfo['gbk_gene_id'] == hnc_dat.loc[index,"L_gene"],"gbk_end"]
        end -= 1
        end = end.to_string(index=False)
        print(end)
        print("--------------------------")
        starts.append(start)
        ends.append(end)
    else:
        start = "NA"
        end = "NA"
        starts.append(start)
        ends.append(end)
#    hnc_dat.loc[ginfo['gbk_gene_id'] == hnc_df.loc[index,"L_gene"],'start'] = ginfo.loc[ginfo['gbk_gene_id'] == hnc_df.loc[index,"L_gene"],"gbk_start"]
#ginfo

hnc_dat['start'] = starts
hnc_dat['end'] = ends

print("HEAD")
print(hnc_dat.head())

# saving the DataFrame as a CSV file 
hnc_dat_csv = hnc_dat.to_csv('../HNS_protein/raw_data_files/higashi_nc_dat.csv', index = True) 
print('\nCSV String:\n', hnc_dat_csv) 


