from sys import argv
import os
import subprocess
import pandas as pd 
import re 

# USAGE EXAMPLE 
# python3 Gather_coverage_from_bam.py /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts/P033_test/work_Bowtie_mapping_nextflow_P033_indiv_filtered/d7/ecb2bcf8dcf7c1a8d334753ed224b0/Mammuthus-P033_Homoserinimonas_aerilata_ALL_mapping0.bam
# Define the paths
#input_directory = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_P033/TOP2000_indiv_filtered/All_extracted_bam"
#output_file = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_P033/TOP2000_indiv_filtered/Bam_files.bamcov"

# List all files ending with '.bam' in the specified directory
#bam_files = [file for file in os.listdir(input_directory) if file.endswith(".bam")]

bam_file=argv[1]

# Create en empty tab with the columns name Scaffold, start, end, N_reads, N_covered_bases, Percent_covered, Avg_cov, Avg_baseq, Avg_mapq
Cov_tab = pd.DataFrame(columns=["Scaffold", "start", "end", "N_reads", "N_covered_bases", "Percent_covered", "Avg_cov", "Avg_baseq", "Avg_mapq"])

# Execute bamcov command for each bam file and save the output to the output file

print(f"Processing file {bam_file}")
#bam_file_path = os.path.join(input_directory, bam_file)

filename = os.path.basename(bam_file)
command = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/TOOLS/bamcov/bamcov -H "+bam_file+" > "+ re.sub("\\.bam",".cov",filename)
subprocess.run(command, shell=True)
print("Process done, adding columns")
# read the output
tab=pd.read_csv(re.sub("\\.bam",".cov",filename), sep='\t', header=None)
# If a genome present zero cov everywhere, remove this : 
tab.columns = ["Scaffold", "start", "end", "N_reads", "N_covered_bases", "Percent_covered", "Avg_cov", "Avg_baseq", "Avg_mapq"]
tab['Genome'] = tab['Scaffold'].str.split('_', n=2).str[:2].str.join('_')
tab = tab.groupby('Genome').filter(lambda x: not (x['N_reads'] == 0).all())
tab.to_csv(re.sub("\\.bam",".cov",filename),sep=";",index=False)
