import zipfile 
import glob
from Bio import Entrez
import json
import os, re 
import urllib
import gzip
import pandas as pd 
Entrez.email = "Benjamin.guinet95@gmail.com"
#import pysam
from Bio import SeqIO 
import argparse
import math 

ALL_assembly="NO"

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(term, download=True, path='assemblies'):
    links=0
    if links==0:
        print('No links found, trying with datasets')
        import subprocess
        term = re.sub(" ","",term)
        # Run this subprocess "datasets download genome accession term" to get the ftp link
        command = f"/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/datasets download genome accession {term}"
        subprocess.run(command, shell=True, cwd=Output_directory)
        print(command)
        # Rename the file created "ncbi_dataset.zip" to the desired filename
        # Unzip the downloaded file
        with zipfile.ZipFile(Output_directory+ "/ncbi_dataset.zip", 'r') as zip_ref:
            zip_ref.extractall(Output_directory)
        # Move the file in ncbi_dataset.zip to the desired folder
        command3 = f"mv {Output_directory}/ncbi_dataset/data/{term}/*.fna {Output_directory}/"
        subprocess.run(command3, shell=True, cwd=Output_directory)
        print(command3)
        # Remove the ncbi_dataset folder
        command4 = f"rm -r {Output_directory}/ncbi_dataset"
        subprocess.run(command4, shell=True, cwd=Output_directory)
        # Remove the ncbi_dataset.zip file
        command5 = f"rm {Output_directory}/ncbi_dataset.zip"
        subprocess.run(command5, shell=True, cwd=Output_directory)
        # Remove README files
        command6 = f"rm {Output_directory}/*.md"
        subprocess.run(command6, shell=True, cwd=Output_directory)
        links = glob.glob(os.path.join(Output_directory, "*"+term+"*.fna"))
        links = links[0]
    return links


#Finds the ids associated with the assembly
def get_ids(term):
    ids = []
    handle = Entrez.esearch(db="assembly", term=term)
    record = Entrez.read(handle)
    ids.append(record["IdList"])
    return ids


#Fetch raw output
def get_raw_assembly_summary(id):
    handle = Entrez.esummary(db="assembly",id=id,report="full")
    record = Entrez.read(handle)
    #Return individual fields
    #XML output: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=79781&report=%22full%22
    #return(record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']) #This will return the Assembly name
    return(record)


#JSON formatted output
def get_assembly_summary_json(id):
    handle = Entrez.esummary(db="assembly",id=id,report="full")
    record = Entrez.read(handle)
    #Convert raw output to json
    return(json.dumps(record, sort_keys=True,indent=4, separators=(',', ': ')))


# Finds the IDs associated with the BIOSAMPLE
def get_ids_biosample(term):
    ids = []
    handle = Entrez.esearch(db="biosample", term=term)
    record = Entrez.read(handle)
    ids.extend(record["IdList"])
    return ids


# Fetch raw output
def get_raw_biosample_summary(id):
    handle = Entrez.esummary(db="biosample", id=id, report="full")
    record = Entrez.read(handle)
    return record


# JSON formatted output
def get_biosample_summary_json(id):
    handle = Entrez.esummary(db="biosample", id=id)
    record = Entrez.read(handle)
    return record


# Print out a message when the program is initiated.
print('-------------------------------------------------------------------------------------------------------\n')
print('                        DOWNLOADING ASSEMBLIES AND METADATA FROM NCBI                                \n')
print('-------------------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-o", "--out_dir",help="The output directory were the files will be created", required=True)
# Create another arg to ask if the user wants to download the assemblies
parser.add_argument("-d", "--download_assemblies",help="Download the assemblies from NCBI", required=False, default=True)
# Create another arg to ask if the user wants to download the metadata from NCBI
parser.add_argument("-m", "--download_metadata",help="Download the metadata from NCBI", required=False, default=True)
# Ask for the table with one unique column containing the accessions
parser.add_argument("-t", "--table",help="The table with the accessions (Header should be named Accessions )", required=True)


# PLEASE RUN module load NCBI-datasets/15.29.0 before running this script

# Example usage : python3 Download_NCBI_assemblies.py -t Accession_table.txt -d -o /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Test_code

args = parser.parse_args()

Output_directory = args.out_dir
Table_with_accessions = pd.read_csv(args.table,sep="\t")

"""
# Toy examples 
Output_directory = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Test"
Table_with_accessions = pd.read_csv("/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/All_species_accessions.csv",sep="\t")
"""

# Remove space in the column Accession
Table_with_accessions['Accessions'] = Table_with_accessions['Accessions'].str.replace(" ", "")

#Output_directory = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Test"

# Specify the directory where you want to save the downloaded files
if not os.path.exists(Output_directory):
        os.makedirs(Output_directory)

#Create an empty df with two columns Accession and Path
df = pd.DataFrame(columns=['Accessions', 'Path'])

accession_numbers = Table_with_accessions['Accessions'].tolist()


if args.download_assemblies is True:
        for accession in accession_numbers:
            if ALL_assembly=="YES":
               assembly="yes"
            else:
               assembly="no"
            if "GCA" in accession :
                assembly="yes"
            if "GCF" in accession :
                assembly="yes"
            if assembly =="yes":
                links = get_assemblies(accession, download=True, path=Output_directory)
                try:
                    file_name_with_extension = os.path.basename(links)
                    print("Accession : ", accession, " downloaded.")
                    # Add the values in the df usinf concat
                    df = pd.concat([df, pd.DataFrame({'Accessions': [accession], 'Path': [Output_directory+"/"+file_name_with_extension]})], ignore_index=True)
                except: 
                    print("Accession : ", accession, " not downloaded WARNING.")
            else: 
                # Fetch the sequence from NCBI
                handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(handle, "fasta")
                handle.close()
                # Save the sequence to a FASTA file
                output_file = os.path.join(Output_directory, f"{accession}.fna")
                SeqIO.write(seq_record, output_file, "fasta")
                print(f"Sequence saved to {output_file}")
                df = pd.concat([df, pd.DataFrame({'Accessions': [accession], 'Path': output_file})], ignore_index=True)
else:
        df["Accessions"] = accession_numbers
        df["Path"] = "Not downloaded"
# Get assembly informations 
print(df)
df = df.dropna(subset=['Path'])
df_no_assembly =df[~df['Accessions'].str.contains("GCA|GCF")]
final_df = pd.DataFrame()
for term in df['Accessions'].unique():
        if ALL_assembly=="YES":
          assembly="yes"
        else:
          assembly="no"
        if "GCA" in term :
                assembly="yes"
        if "GCF" in term :
                assembly="yes"
        if assembly =="yes":
            for id in get_ids(term):
                    #print(get_raw_assembly_summary(id)) #For raw output
                    try:
                        data = json.loads(get_assembly_summary_json(id))
                    except:
                        data = json.loads(get_assembly_summary_json(id[0]))
                    biosample = data['DocumentSummarySet']['DocumentSummary'][0]['BioSampleAccn']
                    print (biosample)
                    accession = term
                    Supressed="NO"
                    if "suppressed" in str(data):
                        Supressed="YES"
                    # Concat the data to the df
                    for id2 in get_ids_biosample(biosample):
                                # print(get_raw_biosample_summary(id)) # For raw output
                                data2=get_biosample_summary_json(id2) # JSON Formatted
                                pattern = r'<Attribute.*?display_name="(.*?)">(.*?)</Attribute>'
                                matches = re.findall(pattern, data2['DocumentSummarySet']['DocumentSummary'][0]['SampleData'])
                                result_list = [(display_name, value) for display_name, value in matches]
                                sub_df=pd.DataFrame(result_list, columns=['Attribute', 'Value'])
                                # put Attibute as columns and value as rows
                                # Only keep columns host   sample type    strain and remove NaN
                                # Pivot the DataFrame
                                df_pivoted = sub_df.pivot(columns='Attribute', values='Value')
                                # Drop NaN columns (like 'collection date' and 'geographic location')
                                df_pivoted = df_pivoted.dropna(axis=1)
                                #
                                sub_df= pd.DataFrame(result_list, columns=['Attribute', 'Value'])
                                # Pivot the DataFrame
                                df_pivoted = pd.pivot_table(sub_df, index=sub_df.index, columns='Attribute', values='Value', aggfunc=lambda x: x.iloc[0])
                                # Reorder columns based on the desired order
                                # Remove every column that are not in that column names [['host', 'sample type', 'strain',"geographic location"]]
                                desired_columns = ['host', 'sample type', 'strain', 'geographic location',"isolation source"]
                                # Remove columns not in the desired column names
                                columns_to_remove = [col for col in df_pivoted.columns if col not in desired_columns]
                                df_pivoted = df_pivoted.drop(columns=columns_to_remove)
                                # Drop rows where all values are NaN
                                df_pivoted = df_pivoted.dropna(how='all')
                                df_pivoted = df_pivoted.fillna(method='ffill')
                                # Group by the columns with common values and concatenate them
                                df_pivoted  = df_pivoted.groupby(df_pivoted .columns.tolist()).first().reset_index()
                                print(df_pivoted)
                                accession = data2['DocumentSummarySet']['DocumentSummary'][0]['Accession']
                                taxonomy = data2['DocumentSummarySet']['DocumentSummary'][0]['Taxonomy']
                                organism = data2['DocumentSummarySet']['DocumentSummary'][0]['Organism']
                                # Add the accession, taxonomy and organism to the DataFrame
                                df_pivoted['BioSample'] = accession
                                df_pivoted['Taxonomy'] = taxonomy
                                df_pivoted['Organism'] = organism
                                df_pivoted['Accessions'] = term
                                df_pivoted['Supressed'] = Supressed
                                final_df = pd.concat([final_df, df_pivoted], ignore_index=True)
                                print(final_df)
                                # Creating DataFrame
# Merge the two DataFrame
print(final_df)
print(df)
df = pd.merge(df, final_df, on='Accessions', how='left')
# Remove every column were Accession is duplicated and keep the row were the value of the column BioSample is not NaN
# Step 2: Identify duplicates in the Accession column
# Replace nan value in "host" by "Missing
df['host'] = df['host'].fillna("missing_host")
# Replace 'not available', 'missing', or 'Not Applicable' in 'host' column with values from 'isolation_source'
df['host'] =df.apply(lambda row: row['isolation source'] if row['host'] in ['not available',"missing_host", 'missing', 'Not Applicable',"nan",'not available: to be reported later'] else row['host'], axis=1)
# If 'host' column still contains 'not available', 'missing', or 'Not Applicable', add 'Missing' to the 'host' column
df['host'] = df.apply(lambda row: 'Missing' if row['host'] in ['not available', 'missing',"missing_host", 'Not Applicable',"nan",'not available: to be reported later'] else row['host'], axis=1)
   # Replace "pig" and pigs by "Pig" in host
df['host'] = df['host'].replace("Missing", 'missing_host')
df['geographic location'] = df['geographic location'].replace(['not determined','Not Applicable'], 'missing_location')
# Add a maj to only the first letter of the host
df['host'] = df['host'].str.capitalize()
df['geographic location'] = df['geographic location'].str.capitalize()
df['geographic location'] = df['geographic location'].str.replace(' ', '')
# replace space by underscore
df['host'] = df['host'].str.replace(' ', '_')
# Add the column host and geographic location to the full path separated by a "_" and add .fna at the end
# Replace any nan in the table by "Missing"
df = df.fillna("Missing")
# Add df_no_assembly rows to df 
df = pd.concat([df, df_no_assembly], ignore_index=True)
print(df)
df.to_csv(Output_directory+"/Genome_metadata.csv", index=False)
print("Metadata file saved to :", Output_directory, "/Genome_metadata.csv")
