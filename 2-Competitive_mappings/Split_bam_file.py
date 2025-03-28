import pysam
import pandas as pd 
import re 
import subprocess
import os
import sys 
import time
import argparse
import math 

# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Split bam files according to contigs for each species     .\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-c", "--contig_file",help="The contig file with all contigs for all species and all mammuth samples")
parser.add_argument("-b", "--bam_file",help="The bam file with all reads for all species and all mammuth samples")
parser.add_argument("-t", "--taxonomy",help="The name of the taxonomy to extract")
parser.add_argument("-d", "--dir_output",help="Desired output dire (default the dir where the script is run)")

args = parser.parse_args()

#Get the uniprot IDS

# Usage exemple python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts/Split_bam_file.py  -t Glacieibacterium_sp018982925 -c /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output/Contig_species_df.csv -b /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts/P033_test/work_Bowtie_mapping_nextflow_P033_concat/34/193bdfcc3d5a2db0e797014cfd3dce/Mammuthus-P033_ALL.bam

directory_output=args.dir_output
contig_file= args.contig_file
taxonomy=args.taxonomy
Bam_file=args.bam_file
Name=re.sub("_ALL.bam","",os.path.basename(Bam_file))
Name=re.sub(".bam","",Name)
Name=re.sub("_sed","",Name)
Name=re.sub("_lab","",Name)
print("Name : ", Name)
print("Bamfile: ", Bam_file)
Contig_file=args.contig_file

#Contig_file="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output/Contig_species_df.csv"
#Contig_file="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_P033/TOP2000_indiv_filtered/Contig_species_df2.csv"
#Name="Mammuthus-MD075"
#Bam_file="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts/work/66/dde8efcf28362d5b90310e540b7848/Mammuthus-MD075_ALL.bam"

# Paths for input and output BAM files
# remove duplicated row in column "Fasta_name" for Kraken_report_ALL_S_copy2

Contig_species_df=pd.read_csv(Contig_file, sep="\t")

print(Contig_species_df)
print( "Taxo : ", taxonomy)
print(Contig_species_df.loc[Contig_species_df["Taxonomy"].str.contains(taxonomy.split("_")[0])])
print(Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1])])
print(len(Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1])]))

if len(Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1])]) > 0:
 print("OOOOKKKK")
Contig_species_df=Contig_species_df.loc[Contig_species_df["Name"].eq(Name)]
if len(Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy)]) > 0:
 print("1")
 Contig_species_df=Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy)]
if len(Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1])]) > 0:
  print("2")
  Contig_species_df=Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1] )]
if len(taxonomy.split("_")) >2:
 if len (Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1]+"_"+taxonomy.split("_")[2] )]) > 0 :
   print("3")
   Contig_species_df=Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1]+"_"+taxonomy.split("_")[2] )]
if len(taxonomy.split("_")) >3:
 if len(Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1]+"_"+taxonomy.split("_")[2]+"_"+taxonomy.split("_")[3])]) >0 :
   Contig_species_df=Contig_species_df.loc[Contig_species_df["Taxonomy"].eq(taxonomy.split("_")[0]+"_"+taxonomy.split("_")[1]+"_"+taxonomy.split("_")[2]+"_"+taxonomy.split("_")[3])]
   print("Nothin found")
  
print(Contig_species_df)
print("############################")
print("Extract reads for sample : ",Name)
print("Species : ", taxonomy)
print("############################")


print(Contig_species_df)

for index, row in Contig_species_df.drop_duplicates(subset=['Fasta_name']).iterrows():
  if directory_output:
    contig_names = Contig_species_df.loc[Contig_species_df["Fasta_name"].str.contains(row["Fasta_name"])]["record.id"].unique() 
    modified_contig_names = [item.replace("|", r"\|") for item in contig_names]
    if len(modified_contig_names) > 500:
        print("To much scaffold, creating chuncks...")
        num_chunks = int(math.ceil(len(modified_contig_names) / 500.0))
        for i in range(num_chunks):
            chunk = modified_contig_names[i * 500: (i + 1) * 500]
            if "_sed.bam" in os.path.basename(Bam_file):
                extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", f"_extracted_part_sed_{i}.bam")
                print("extracted_bam_file :", extracted_bam_file)
            elif "_lab.bam" in os.path.basename(Bam_file):
                extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", f"_extracted_part_lab_{i}.bam")
            else:
                extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", f"_extracted_part_{i}.bam")
            extracted_bam_file_PMDscores = re.sub(r".bam", r"_PMDscores.txt", extracted_bam_file)
            extracted_bam_file = re.sub(r"\(", r"_", extracted_bam_file)
            extracted_bam_file = re.sub(r"\)", r"_", extracted_bam_file)
            extracted_bam_file = re.sub(r"\[", r"_", extracted_bam_file)
            extracted_bam_file = re.sub(r"\/", r"_", extracted_bam_file)
            extracted_bam_file = re.sub(r"\]", r"_", extracted_bam_file)
            print("Extract reads for species: ", row["Taxonomy"])
            command1 = ["samtools", "view -h -b", Bam_file, " ".join(list(set(chunk))), ">", directory_output + extracted_bam_file]
            subprocess.run(" ".join(command1), shell=True, check=True)
        merged_bam_file = extracted_bam_file.replace(".bam", "_merged.bam")
        extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_extracted.bam")
        extracted_bam_file_PMDscores = re.sub(r".bam", r"_PMDscores.txt", extracted_bam_file)
        extracted_bam_file = re.sub(r"\(", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\)", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\[", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\/", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\]", r"_", extracted_bam_file)
        if "_sed.bam" in os.path.basename(Bam_file):
            extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_sed_extracted.bam")
            print("extracted_bam_file2 ",extracted_bam_file)
            print(re.sub("sed_extracted.bam","",extracted_bam_file)) 
            merge_command = "samtools merge "+ directory_output + extracted_bam_file  + " " +  directory_output + re.sub("sed_extracted.bam","",extracted_bam_file) + "*_part_sed_*.bam"
        elif "_lab.bam" in os.path.basename(Bam_file):
            extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_lab_extracted.bam")
            merge_command = "samtools merge "+ directory_output + extracted_bam_file + " " + directory_output+re.sub("lab_extracted.bam","",extracted_bam_file) +"*_part_lab_*.bam"
        else:
            extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_extracted.bam")
            merge_command = "samtools merge "+ directory_output + extracted_bam_file + " *_part_*.bam"
        print(merge_command)
        subprocess.run(merge_command,shell=True)
        print("All reads extracted and saved here: ", extracted_bam_file)
    else:
        if "_sed.bam" in os.path.basename(Bam_file):
            extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_sed_extracted.bam")
        elif "_lab.bam" in os.path.basename(Bam_file):
            extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_lab_extracted.bam")
        else:
            extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_extracted.bam")
        extracted_bam_file_PMDscores = re.sub(r".bam", r"_PMDscores.txt", extracted_bam_file)
        extracted_bam_file = re.sub(r"\(", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\)", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\[", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\/", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\]", r"_", extracted_bam_file)
        print("Extract reads for species: ", row["Taxonomy"])
        command1 = ["samtools", "view -h -b", Bam_file, " ".join(list(set(modified_contig_names))), ">",directory_output + extracted_bam_file]
        print(command1)
        subprocess.run(" ".join(command1), shell=True, check=True)
        print("2-All reads extracted and saved here: ", extracted_bam_file)
  else:
    contig_names = Contig_species_df.loc[Contig_species_df["Fasta_name"].str.contains(row["Fasta_name"])]["record.id"].unique()
    modified_contig_names = [item.replace("|", r"\|") for item in contig_names]
    if len(modified_contig_names) > 500:
        num_chunks = int(math.ceil(len(modified_contig_names) / 500.0))
        for i in range(num_chunks):
            chunk = modified_contig_names[i * 500: (i + 1) * 500]
            extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", f"_extracted_part_{i}.bam")
            extracted_bam_file_PMDscores = re.sub(r".bam", r"_PMDscores.txt", extracted_bam_file)
            extracted_bam_file = re.sub(r"\(", r"_", extracted_bam_file)
            extracted_bam_file = re.sub(r"\)", r"_", extracted_bam_file)
            extracted_bam_file = re.sub(r"\[", r"_", extracted_bam_file)
            extracted_bam_file = re.sub(r"\/", r"_", extracted_bam_file)
            extracted_bam_file = re.sub(r"\]", r"_", extracted_bam_file)
            print("Extract reads for species: ", row["Taxonomy"])
            command1 = ["samtools", "view -h -b", Bam_file, " ".join(list(set(chunk))), ">", extracted_bam_file]
            subprocess.run(" ".join(command1), shell=True, check=True)
        merged_bam_file = extracted_bam_file.replace(".bam", "_merged.bam")
        extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_extracted.bam")
        extracted_bam_file_PMDscores = re.sub(r".bam", r"_PMDscores.txt", extracted_bam_file)
        extracted_bam_file = re.sub(r"\(", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\)", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\[", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\/", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\]", r"_", extracted_bam_file)
        merge_command = "samtools merge "+extracted_bam_file + " *_part_*.bam"
        print(merge_command)
        #for i in range(num_chunks):
        #    sub_bam_file = extracted_bam_file.replace(".bam", f"_part_{i}.bam")
        #    merge_command.append(sub_bam_file)
        subprocess.run(merge_command,shell=True)
        print("All reads extracted and saved here: ", merged_bam_file)
    else:
        extracted_bam_file = Name + "_" + row["Taxonomy"] + "_" + row["Fasta_name"].replace(".fna", "_extracted.bam")
        extracted_bam_file_PMDscores = re.sub(r".bam", r"_PMDscores.txt", extracted_bam_file)
        extracted_bam_file = re.sub(r"\(", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\)", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\[", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\/", r"_", extracted_bam_file)
        extracted_bam_file = re.sub(r"\]", r"_", extracted_bam_file)
        print("Extract reads for species: ", row["Taxonomy"])
        command1 = ["samtools", "view -h -b", Bam_file, " ".join(list(set(modified_contig_names))), ">", extracted_bam_file]
        print(command1)
        subprocess.run(" ".join(command1), shell=True, check=True)
        print("All reads extracted and saved here: ", extracted_bam_file)
