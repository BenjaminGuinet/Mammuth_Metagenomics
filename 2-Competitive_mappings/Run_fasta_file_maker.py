import pandas as pd
import os
from Bio import SeqIO
import re
import warnings
from sys import argv
import numpy as np

"""
# Eg usage : python3 Run_fasta_file_maker.py /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/\
 /crex/proj/snic2022-6-144/nobackup/NIKOLAY/GTDB/GTDB_KRAKEN2_DB/genomes_from_gtdbtk/GTDB_KRAKEN2_DB/library/added/\
 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/Genome_name_files_GTDB.txt\
 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/Kraken_report_ALL_S_copy.csv\
 Mammuthus-FK008
"""

Output_dir = argv[1]
Genome_path = argv[2]
GTDB_genome_files_df = argv[3]
Kraken_report_ALL_S_table = argv[4]
Name = argv[5]

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


for type in ["ALL_concat"]:
 print("Processing : ", type)
 print("\n")
 #Output_dir = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_10_best/"
 #Genome_path = "/crex/proj/snic2022-6-144/nobackup/NIKOLAY/GTDB/GTDB_KRAKEN2_DB/genomes_from_gtdbtk/GTDB_KRAKEN2_DB/library/added/"
 #GTDB_genome_files_df = pd.read_csv(os.path.join(Output_dir, "Genome_name_files_GTDB.txt"), sep="\t", header=None)
 Contig_species_data = []
 Kraken_report_ALL_S_copy = pd.read_csv(Kraken_report_ALL_S_table, sep="\t",low_memory=False)
 Kraken_report_ALL_S_copy = Kraken_report_ALL_S_copy.loc[Kraken_report_ALL_S_copy["Name"].eq(Name)]
 Kraken_report_ALL_S_copy.loc[Kraken_report_ALL_S_copy["species"].isna(),"species"]=Kraken_report_ALL_S_copy["Taxonomy"]
 Kraken_report_ALL_S_copy['species2'] = Kraken_report_ALL_S_copy['species'].str.replace("uncultured_","")
 Kraken_report_ALL_S_copy['genus'] = Kraken_report_ALL_S_copy['species2'].str.split('_').str[0]
 Kraken_report_ALL_S_copy['genus2'] = Kraken_report_ALL_S_copy['Taxonomy'].str.split('_').str[0]
 Kraken_report_ALL_S_copy = Kraken_report_ALL_S_copy.sort_values(by=['ReadsDirect'], ascending=False)
 Kraken_report_ALL_S_copy["genus"]=Kraken_report_ALL_S_copy["genus"].str.replace("]","")
 Kraken_report_ALL_S_copy["genus"]=Kraken_report_ALL_S_copy["genus"].str.replace("[","")
 if type =="ALL_concat":
  Kraken_report_ALL_S_copy = Kraken_report_ALL_S_copy.drop_duplicates(subset=['genus',"Name"], keep="first")
  Concatenated="YES"
 else:
  Concatenated="NO"
 print(len(Kraken_report_ALL_S_copy["Name"].unique()))
 i=1
 for Name in Kraken_report_ALL_S_copy["Name"].unique():
    print("############################")
    print("Copy files for sample: ", Name)
    print("############################")
    working_directory = os.path.join(Output_dir, "Concat_refs")
    os.makedirs(working_directory,exist_ok=True)

    Kraken_report_ALL_S_copy2 = Kraken_report_ALL_S_copy.loc[Kraken_report_ALL_S_copy["Name"].str.contains(Name)]
    #Kraken_report_ALL_S_copy2 = Kraken_report_ALL_S_copy2.sort_values(by=["ReadsDirect"], ascending=False).head(Number_top_hits)
    print(Kraken_report_ALL_S_copy2[["Name","species","ReadsDirect"]])
    print(Kraken_report_ALL_S_copy2)

    fasta_path = Name+"_ALL.fna"
    if Concatenated=="YES":
     with open(fasta_path, 'w') as output_fasta_all:
        print("Writting a concatenation of all fasta files...")
        for index, row in Kraken_report_ALL_S_copy2.iterrows():
            file_path = os.path.join(Genome_path, row["Fasta_name"])
            record_dict = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
            sequences_without_acgt = [name for name, record in record_dict.items() if all(base  not in record.seq for base in ['A', 'C', 'T', 'G', 'a', 'c', 't', 'g'])]

            for name in sequences_without_acgt:
                del record_dict[name]

            for record_id, record in record_dict.items():
                sequence = record.seq
                print(f">{record.id}", file=output_fasta_all)
                print(sequence, file=output_fasta_all)
                new_row = {"Name": Name, "species": row["species"],"Taxonomy": row["Taxonomy"], "Kraken_taxid" : row["Kraken_taxid"], "record.id": record_id, "Fasta_name": row["Fasta_name"]}
                Contig_species_data.append(new_row)
    elif Concatenated=="NO":
        print("No concatenation, writting each fasta individually...")
        for index, row in Kraken_report_ALL_S_copy2.iterrows():
          species_name=row["Taxonomy"]
          species_name=re.sub(r"\(", r"_", species_name)
          species_name=re.sub(r"\)", r"_", species_name)
          species_name=re.sub(r"\[", r"_", species_name)
          species_name=re.sub(r"\/", r"_", species_name)
          species_name=re.sub(r"\]", r"_", species_name)
          fasta_path = working_directory+"/"+Name+"_"+species_name+"_ALL.fna2"
          with open(fasta_path, 'w') as output_fasta_all:
            file_path = os.path.join(Genome_path, row["Fasta_name"])
            record_dict = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
            #print(file_path)
            sequences_without_acgt = [name for name, record in record_dict.items() if all(base not in record.seq for base in ['A', 'C', 'T', 'G', 'a', 'c', 't', 'g'])]
            for name in sequences_without_acgt:
                del record_dict[name]

            for record_id, record in record_dict.items():
                sequence = record.seq
                print(f">{record.id}", file=output_fasta_all)
                print(sequence, file=output_fasta_all)
                new_row = {"Name": Name, "species": row["species"], "Taxonomy": row["Taxonomy"], "Kraken_taxid" : row["Kraken_taxid"],"record.id": record_id, "Fasta_name": row["Fasta_name"]}
                Contig_species_data.append(new_row)

    i+=1
    print(i,"/",len(Kraken_report_ALL_S_copy["Name"].unique()))
    print("All done.\n")

 Contig_species_df = pd.DataFrame(Contig_species_data, columns=["Name", "species","Taxonomy","Kraken_taxid", "record.id", "Fasta_name"])
 Contig_species_df.to_csv(Output_dir+"/Contig_species_df.csv", sep="\t", index=False, mode='a') 
 #if os.path.exists(Output_dir+"/Contig_species_df.csv"):
 # Contig_species_df_main = pd.read_csv(Output_dir+"/Contig_species_df.csv", sep="\t")
 # Contig_species_df_main=pd.concat([Contig_species_df_main, Contig_species_df])
 # Contig_species_df_main.to_csv(Output_dir+"/Contig_species_df.csv", sep="\t", index=False)
 #else:
 # Contig_species_df.to_csv(Output_dir+"/Contig_species_df.csv", sep="\t", index=False)
 print("Type : ", type, " done.")
 print("\n")
