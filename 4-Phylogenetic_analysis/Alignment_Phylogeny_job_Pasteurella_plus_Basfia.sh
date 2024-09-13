#!/bin/sh
#SBATCH -A naiss2023-5-251
#SBATCH -p node
#SBATCH -M snowy
#SBATCH -t 48:00:00
#SBATCH -J Pasteurella_ali
#SBATCH -o /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Script_ali_and_phylo_animals/Basfia_plus_Pasteurella_ali.out
#SBATCH -e /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Script_ali_and_phylo_animals/Basfia_plus_Pasteurella.error

module load conda
source conda_init.sh
conda activate My_conda
module load iqtree/2.2.2.6-omp-mpi

cd /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia_plus_Pasteurella/

python3  /crex/proj/snic2022-6-144/nobackup/BENJAMIN/TOOLS/Genome_to_MSA3.py\
 -g /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia_plus_Pasteurella/Genome_assemblies.txt\
 -o /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia_plus_Pasteurella/Genome_to_MSA_corrected/\
 -nsites 5000000\
 -K 15 -b 600 -t 16\
 -NJ no\
 -map_prog bowtie2\
 -read_map_tab /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia_plus_Pasteurella/Read_to_map.txt

iqtree2 -s /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia_plus_Pasteurella/Genome_to_MSA_corrected/output_sequences.fasta -m TVM+F+I+R9 --abayes --lbp 1000  -B 1000 -alrt 1000 -nt 16 -redo -bnni

#cp /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Pasteurella/Genome_to_MSA/*.treefile /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/All_trees/
#cp /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Pasteurella/*.csv crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/All_trees/
