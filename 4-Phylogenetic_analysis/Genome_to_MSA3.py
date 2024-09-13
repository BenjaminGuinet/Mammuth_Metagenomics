import argparse
from Bio import AlignIO
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os,re 
import pandas as pd
import subprocess
import numpy as np 
import time
import random

print("\n")

print("""
██████╗ ███████╗███╗   ██╗ ██████╗ ███╗   ███╗███████╗███████╗    ████████╗ ██████╗     ███╗   ███╗███████╗ █████╗ 
██╔════╝ ██╔════╝████╗  ██║██╔═══██╗████╗ ████║██╔════╝██╔════╝    ╚══██╔══╝██╔═══██╗    ████╗ ████║██╔════╝██╔══██╗
██║  ███╗█████╗  ██╔██╗ ██║██║   ██║██╔████╔██║█████╗  ███████╗       ██║   ██║   ██║    ██╔████╔██║███████╗███████║
██║   ██║██╔══╝  ██║╚██╗██║██║   ██║██║╚██╔╝██║██╔══╝  ╚════██║       ██║   ██║   ██║    ██║╚██╔╝██║╚════██║██╔══██║
╚██████╔╝███████╗██║ ╚████║╚██████╔╝██║ ╚═╝ ██║███████╗███████║       ██║   ╚██████╔╝    ██║ ╚═╝ ██║███████║██║  ██║
 ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝╚══════╝       ╚═╝    ╚═════╝     ╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝
""")

print("\n")

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow to generate a MSA from a list of genomes using Sibeliaz.')
parser.add_argument("-g","--genome_path", help="A txt file with one column with the full path of the genomes to align")
parser.add_argument("-o","--output_dir", help="The output directory where will be saved files including the alignement file")
parser.add_argument("-nsites","--nsites", help="The number of desired sites in the final MSA (taken randomly)")
parser.add_argument("-b","--b", help="Length of the  chains of bubbles default is 200")
parser.add_argument("-K","--K", help="For small datasets, like bacteria, we recommend k=15, and for mammalian-sized genomes k=25. The default is 25.")
parser.add_argument("-t","--threads", help="Number of threads to use.")
#parser.add_argument("-Taxon","--Taxon", help="comma separated list of taxa for wich we want to force to include all sites in the MSA, in that case number of final sites >= nsites. Name of Taxa should be the name of the Genome file without the full path and remove everythin after the '.'")
parser.add_argument("-NJ","--Run_NJ_phylogeny", help="Set to yes if you want a quick phylogeny made with rapidnj.")
parser.add_argument("-read_map_tab","--read_map_tab", help="The 2 column table with in the first column the full path to the fastq file and in the second the name of the refence genome (as to be in -g)")
parser.add_argument("-map_prog", "--Mapping_soft", help="Write bowtie2 or bwa")

args = parser.parse_args()

# usage exemple : python3 Genome_to_MSA.py -g Candidates_Pasteurella.txt -o Genome_to_MSA_dir -nsites 10000 -K 15 -t 16 -Taxon Mammuthus-FK033_ALL_new_Bisgaard_extracted_angsd_mapq30_mindepth1,Mammuthus-MD024_ALL_new_Bisgaard_extracted_angsd_mapq30_mindepth1 --Run_NJ_phylogeny yes


print ("Prior running this script, you will need to load the conda environment with the sibeliaz software")
print ("Can be installed using : conda install sibeliaz (but you will need to first have you own conda environment)")
print("\n")
print("Then load your env :")
print("module load conda")
print("source conda_init.sh")
print("conda activate My_conda")

print("\n")

print("Sibeliaz is a whole-genome alignment and locally-coliinear blocks construction pipeline and was made by :")
print("\n")
print("- Ilia Minkin")
print("- Paul Medvedev")
print("\n")

print("Sibeliaz method targets closely related sequences, such as strains of the same species.")

sibeliaz_output_dir=args.output_dir
Input_genome_tab=args.genome_path
nsites=args.nsites
K=args.K
bubble=args.b
#Taxon=args.Taxon
#if args.Taxon:
# Taxon_selected_list = Taxon.split(",")
threads=args.threads
Run_NJ_phylogeny=args.Run_NJ_phylogeny
Read_map_tab=args.read_map_tab

#Read_map_tab="Erysipelothrix_read_to_map.txt"

"""
# For test purpose 
Input_genome_tab="Candidates_Pasteurella.txt"
sibeliaz_output_dir="Genome_to_MSA_dir"
nsites=10000
K=15
threads=16
Taxon="Mammuthus-FK033_ALL_new_Bisgaard_extracted_angsd_mapq30_mindepth1,Mammuthus-MD024_ALL_new_Bisgaard_extracted_angsd_mapq30_mindepth1"
Taxon_selected_list = Taxon.split(",")
"""


# get the time 
start_time = time.time()

Candidate_genome_tab=pd.read_csv(Input_genome_tab,sep="\t",header=None)

# Rename columns
Candidate_genome_tab.columns = ["Genome_path"]
print("\t")
print("#########################################")
print("Number of genomes : ", len(Candidate_genome_tab["Genome_path"].unique()))
print("#########################################")
# Function to read contig names from a FASTA file
def get_contig_names(fasta_path):
    contig_names = [record.id for record in SeqIO.parse(fasta_path, "fasta")]
    return contig_names

# Apply the function to each row in the DataFrame and store the results in a new column
Candidate_genome_tab['contig_names'] = Candidate_genome_tab["Genome_path"].apply(get_contig_names)
# Add a third column called Species_name which is the name of the file without the full path and remove everything after the "." 
Candidate_genome_tab['Species_name'] = Candidate_genome_tab["Genome_path"].apply(lambda x: re.sub(r'\..*$', '', os.path.basename(x)))


# if sibeliaz_output_dir+"/output_sequences.fasta" already exists

if os.path.exists(sibeliaz_output_dir+"/output_sequences.fasta"):
    print("\n")
    print("The file ",sibeliaz_output_dir+"/output_sequences.fasta"," already exists. The script will not run again.")
    print("\n")
    print("If you want to run the script again, please remove the file ",sibeliaz_output_dir+"/output_sequences.fasta"," and run the script again.")
    print("\n")
else:
    print("\n")
    print("Running sibeliaz aligner ... ")
    print("\n")
    #
    # Run the aligner
    # Joining the paths from the DataFrame's column with spaces
    paths = " ".join(Candidate_genome_tab['Genome_path'])
    # Command prefix
    command_prefix = "sibeliaz"
    # Full command
    command = f"{command_prefix} -o {sibeliaz_output_dir} -k {K} -b {bubble} -t {threads} {paths}"
    # Execute the command
    print(command)
    subprocess.run(command,shell=True, universal_newlines=True)
    print(command)
for i in ["A"]:
    # Parse the MAF file geenrated by Sibeliaz
    def parse_maf(file_path):
        alignment_blocks = []
        current_block = None  # Initialize current_block outside the loop
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('a'):  # Start of a new alignment block
                    if current_block is not None and len(current_block['sequences']) > 1 :
                        alignment_blocks.append(current_block)  # Add only if it has exactly 3 sequences
                    current_block = {'sequences': []}
                elif line.startswith('s') and current_block is not None:  # Ensure current_block is initialized
                    parts = line.strip().split()
                    # Extracting sequence details: src, start, size, strand, srcSize, sequence
                    sequence_info = {
                        'src': parts[1],
                        'start': int(parts[2]),
                        'size': int(parts[3]),
                        'strand': parts[4],
                        'srcSize': int(parts[5]),
                        'sequence': parts[6]
                    }
                    current_block['sequences'].append(sequence_info)
                # You might want to handle lines that do not start with 'a' or 's' here
            if current_block is not None and len(current_block['sequences']) > 1:  # Check for the last block
                alignment_blocks.append(current_block)
        return alignment_blocks
    #
    print("\n")
    print("Sibeliaz done.")
    print("\n")
    print("\n")
    print("Processing maf results ...")
    print("\n")
    file_path = sibeliaz_output_dir+'/alignment.maf'
    alignment_blocks = parse_maf(file_path)
    # Make a DataFrame with the sequences with on each row the blocks and on each column the species
    df_columns = Candidate_genome_tab["Species_name"].unique()
    df2 = pd.DataFrame(columns=df_columns)
    rows_list = []  # Initialize an empty list to store rows before concatenation
    for data in alignment_blocks:
        row = {file: "" for file in df_columns}  # Initialize a row with empty strings for each file
        for seq in data['sequences']:
            # Try to find the corresponding file name using the mapping DataFrame
            try: 
                file_name = Candidate_genome_tab[Candidate_genome_tab['contig_names'].apply(lambda x: seq['src'] in x)]["Species_name"].iloc[0]
                row[file_name] = seq['sequence']  # Update the row with the sequence
            except IndexError:
                # This handles cases where the src is not found in the mapping DataFrame
                continue  # Skip to the next sequence if src not found
        rows_list.append(row)
    # Use pd.concat to append all rows at once for better performance
    df2 = pd.concat([df2, pd.DataFrame(rows_list)], ignore_index=True)
    # Replace empty cells wich are block not presenting any homology with the other species with "-" repeated of the sizes of the block alignment
    def replace_empty_with_avg(row):
        # Filter non-empty cells and calculate average length of non-empty cells
        non_empty_cells = [str(cell) for cell in row if cell]
        avg_length = np.mean([len(cell) for cell in non_empty_cells]) if non_empty_cells else 0
        # Round the average length to nearest whole number
        avg_length_rounded = int(round(avg_length))
        # Replace empty cells with '-' repeated 'avg_length_rounded' times
        return [cell if cell else '-' * avg_length_rounded for cell in row]
    # Applying the function to each row
    df3 = df2.apply(replace_empty_with_avg, axis=1, result_type='expand')
    df3.columns = df_columns  
    # Remove sites only composed of gaps
    def remove_common_gaps_df(df):
        # Create a copy of the DataFrame to avoid modifying the original data
        modified_df = df.copy() 
        # Iterate over each row by index
        for index, row in df.iterrows():
            # Convert the row to a list of strings
            row_list = row.tolist()
            row_length = len(row_list[0])  # Assuming all cells in a row are of the same length
            common_dash_positions = []  # To store positions where dash is common across all cells
            # Identify common dash positions
            for i in range(row_length):
                if all(cell[i] == '-' for cell in row_list):
                    common_dash_positions.append(i)      
            # Remove dashes from identified positions for each cell in the row
            for col in df.columns:
                modified_df.at[index, col] = ''.join([char for idx, char in enumerate(df.at[index, col]) if idx not in common_dash_positions])
        return modified_df
    #
    df4= remove_common_gaps_df(df3)
    #
    # Write the result into a final MSA file
    fasta_content = ""
    for col in df4.columns:
        # Concatenate all cell values in this column, excluding "--"
        sequence = ''.join(df4[col][df4[col] != "--"])
        # Add this column's data to the FASTA content
        fasta_content += f">{col}\n{sequence}\n"
    #
    # Save the FASTA content to a file
    with open(sibeliaz_output_dir+"/output_sequences.fasta", "w") as fasta_file:
        fasta_file.write(fasta_content)

print("\n")
print ("Mapping reads back to aligned best references in the alignment ...") 

Candidate_read_to_map_tab=pd.read_csv(Read_map_tab,sep="\t",header=None)

# Remove everything after the first "." in the column 1

Candidate_read_to_map_tab[1]=Candidate_read_to_map_tab[1].apply(lambda x: re.sub(r'\..*$', '', x))

#sibeliaz_output_dir="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/ZOE_project/Erysipelothrix_genomes_ali_phylo/Genome_to_MSA_2"

multifasta = SeqIO.to_dict(SeqIO.parse(sibeliaz_output_dir+"/output_sequences.fasta", "fasta"))

# Read dataframe
# Create a directory to store the bowtie2 index
os.makedirs(sibeliaz_output_dir+"/Aligned_contig_to_map/", exist_ok=True)

# Loop through each row
for index, row in Candidate_read_to_map_tab.iterrows():
    SeqIO.write(multifasta[row[1]], sibeliaz_output_dir+"/Aligned_contig_to_map/"+row[1]+".fasta", "fasta")
    # Build bowtie2 index for the reference genome
    ref_genome = sibeliaz_output_dir+"/Aligned_contig_to_map/"+row[1]+".fasta"
    fastq_file = row[0]
    # Change dir to the output dir
    os.chdir(sibeliaz_output_dir+"/Aligned_contig_to_map/")
    # Constructing the bowtie2-build command
    if args.Mapping_soft=="bowtie2":
     command1 = ["bash /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/ZOE_project/Scripts/Map_to_ref_bowtie2.sh ", ref_genome, " ", fastq_file," ", sibeliaz_output_dir+"/Aligned_contig_to_map/"]
     subprocess.run("".join(command1), shell=True)
    else:
     print("Default mapping with  BWA")
     command1 = ["bash /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/ZOE_project/Scripts/Map_to_ref_bwa2.sh ", ref_genome, " ", fastq_file," ", sibeliaz_output_dir+"/Aligned_contig_to_map/ ", "clip_no"]
     print("COMMAND TEST ")
     print("".join(command1))
     subprocess.run("".join(command1), shell=True)
# Now within the directory sibeliaz_output_dir+"/Aligned_contig_to_map/", cat all the files ending with "_for_phylogeny.fa" within the file sibeliaz_output_dir+"/output_sequences.fasta" using subprocess.run
# Constructing the cat command
command2 = ["cat *_for_phylogeny.fa >> ", sibeliaz_output_dir+"/output_sequences.fasta"]
subprocess.run("".join(command2), shell=True)


print("\n")
print("MSA file containing all sites saved to : ",sibeliaz_output_dir+"/output_sequences.fasta")
print("\n")


print("\n")
print("Trimming sites ... ")
print("\n")

# Make taxon_selected_list being the name of the contigs from all files  ending with "_for_phylogeny.fa" in the directory sibeliaz_output_dir+"/Aligned_contig_to_map/"
# Get the list of files ending with "_for_phylogeny.fa"
files = [f for f in os.listdir(sibeliaz_output_dir+"/Aligned_contig_to_map/") if f.endswith("_for_phylogeny.fa")]
# Get the list of contigs from all files ending with "_for_phylogeny.fa"
Taxon_selected_list = []
for file in files:
    contigs = list(SeqIO.parse(sibeliaz_output_dir+"/Aligned_contig_to_map/"+file, "fasta"))
    for contig in contigs:
        Taxon_selected_list.append(contig.id)

def extract_sites(input_fasta, output_fasta, taxon_selected_list, nsites):
    nsites = int(nsites)
    sequences = list(SeqIO.parse(input_fasta, "fasta"))  # Read the FASTA file
    total_sites = len(sequences[0].seq)
    coverage = np.zeros(total_sites)  # Array to hold coverage count for each site
    # Update coverage for each site based on all sequences
    for record in sequences:
        for i, base in enumerate(str(record.seq)):
            if base not in ['N', '-']:  # Count only valid bases
                coverage[i] += 1
    selected_sites = set()
    # If taxa are specified, prioritize their sites
    if Taxon_selected_list is not None:
        for record in sequences:
            if record.id in Taxon_selected_list:
                for i, base in enumerate(str(record.seq)):
                    if base not in ['N', '-']:
                        selected_sites.add(i)
                        coverage[i] = -1  # Mark as selected to avoid re-selection
    # Sites not yet selected
    remaining_sites = np.argsort(coverage)[::-1]  # Sort sites by coverage, descending
    # Select additional sites based on coverage, excluding already selected
    additional_sites_needed = nsites - len(selected_sites)
    for i in range(total_sites):
        if additional_sites_needed <= 0:
            break
        if coverage[remaining_sites[i]] > 0:  # If site has positive coverage and not already selected
            selected_sites.add(remaining_sites[i])
            additional_sites_needed -= 1
    # Create trimmed sequences
    trimmed_sequences = []
    for record in sequences:
        trimmed_seq = ''.join([str(record.seq[i]) for i in sorted(selected_sites)])
        trimmed_sequences.append(SeqIO.SeqRecord(Seq(trimmed_seq), id=record.id, description=""))
    # Write the output FASTA file
    SeqIO.write(trimmed_sequences, output_fasta, "fasta")

# if Taxon is specified : 
if 'Taxon_selected_list' in locals():
    extract_sites(sibeliaz_output_dir+"/output_sequences.fasta", sibeliaz_output_dir+"/output_sequences_trimmed.fasta", Taxon_selected_list,nsites)
else:
    extract_sites(sibeliaz_output_dir+"/output_sequences.fasta", sibeliaz_output_dir+"/output_sequences_trimmed.fasta", None,nsites)

print("\n")
print("\n")

def count_non_dash_n(sequence):
    return sum(1 for base in sequence if base not in ['-', 'N'])

def create_dataframe_from_fasta(file_path):
    data = {'Contig': [], 'Non-Dash-N_Sites': []}
    for record in SeqIO.parse(file_path, 'fasta'):
        contig_name = record.id
        non_dash_n_count = count_non_dash_n(record.seq)
        data['Contig'].append(contig_name)
        data['Non-Dash-N_Sites'].append(non_dash_n_count)
    df = pd.DataFrame(data)
    return df

# Example usage:
df = create_dataframe_from_fasta(sibeliaz_output_dir+"/output_sequences.fasta")
Sub_trimmed_df = create_dataframe_from_fasta(sibeliaz_output_dir+"/output_sequences_trimmed.fasta")

if Run_NJ_phylogeny == "yes":
                print ("Running NJ phylogeny method")
                rapidnj_path = "/home/benjguin/TOOLS/rapidNJ/bin/rapidnj"
                command_phylo = [rapidnj_path, sibeliaz_output_dir+"/output_sequences_trimmed.fasta", "-i", "fa", "-c", str(threads), "-b", "100", "-t", "d", "-x", sibeliaz_output_dir+"/output_sequences_trimmed.nwk"]
                subprocess.run(command_phylo)
                print("\n")
                print("Final phylogeny written to : ",sibeliaz_output_dir+"/output_sequences_trimmed.nwk")
                phylo = Phylo.read(sibeliaz_output_dir+"/output_sequences_trimmed.nwk", "newick")
                Phylo.draw_ascii(phylo)

print("\n")
print("\n")

# end time
end_time = time.time()

print("MSA trimmed file saved to : ",sibeliaz_output_dir+"/output_sequences_trimmed.fasta")
print("Time to run the script : ",end_time - start_time)
print("-----------------------------------------------------------------------------------\n")


