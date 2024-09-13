import pysam
import pandas as pd
import re
import subprocess
import os
import sys
import time
import argparse

# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Split bam files according to contigs for each species     .\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-s", "--sample", help="The name of the mammuth sample")
parser.add_argument("-c", "--contig_file",help="The contig file with all contigs for all species and all mammuth samples")
parser.add_argument("-o", "--output",help="The tab output")
parser.add_argument("-t", "--type",help="If all takes all the samples")
args = parser.parse_args()

# Usage exemple python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts/Create_split_table.py -s Mammuthus-P033  -c /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output/Contig_species_df.csv -o TOP2000_taxonomy_to_split.txt


contig_file= args.contig_file
#contig_file="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_P033/TOP2000_indiv_filtered/Contig_species_df2.csv"

if args.type =="all":
 Contig_species_df=pd.read_csv(contig_file, sep="\t",low_memory=False)
 Contig_species_df=Contig_species_df.loc[Contig_species_df["Name"].str.contains("Mammuthus")]
else:
 Name=args.sample
 Name=re.sub("_ALL","",Name)
 Output=args.output
 Contig_species_df=pd.read_csv(contig_file, sep="\t",low_memory=False)
 Contig_species_df=Contig_species_df.loc[Contig_species_df["Name"].eq(Name)]

print(Contig_species_df)
print("############################")
print("Create species table for further extraction for sample : ",Name)
print("############################")
print(Contig_species_df)
#Contig_species_df=Contig_species_df.loc[Contig_species_df["species"].str.contains("Pseudomonas_fluorescens")]

## Remove duplicated
Contig_species_df2=Contig_species_df.drop_duplicates(subset=['Taxonomy'], keep='first')
Contig_species_df2=Contig_species_df2["Taxonomy"]
#Contig_species_df2.to_csv(Output,index=False)

if not os.path.exists("tmp"):
	os.makedirs("tmp")
for i in Contig_species_df["Taxonomy"].unique():
	with open("tmp/"+i+".tmp", 'w') as f:
    		pass

print("Split file table written to : ", Output) 
