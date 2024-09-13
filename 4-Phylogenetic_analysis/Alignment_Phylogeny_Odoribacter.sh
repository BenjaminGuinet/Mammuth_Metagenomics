#!/bin/sh
#SBATCH -A naiss2024-5-131
#SBATCH -p node
#SBATCH -t 48:00:00
#SBATCH -J Odoribacter_ali
#SBATCH -o /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Script_ali_and_phylo_animals/Odoribacter_ali.out
#SBATCH -e /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Script_ali_and_phylo_animals/Odoribacter_ali.error

module load conda
source conda_init.sh
conda activate My_conda
module load iqtree/2.2.2.6-omp-mpi

cd /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Odoribacter/

python3  /crex/proj/snic2022-6-144/nobackup/BENJAMIN/TOOLS/Genome_to_MSA3.py\
 -g /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Odoribacter/Genome_assemblies.txt\
 -o /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Odoribacter/Genome_to_MSA/\
 -nsites 5000000\
 -K 15 -b 600 -t 16\
 -NJ no\
 -map_prog bowtie2\
 -read_map_tab /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Odoribacter/Read_to_map.txt

iqtree2 -s /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Odoribacter/Genome_to_MSA/output_sequences.fasta -m K3Pu+F+I+R10 --abayes --lbp 1000  -B 1000 -alrt 1000 -nt 16 -redo -bnni

#cp /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Odoribacter/Genome_to_MSA/*.treefile /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/All_trees/
#cp /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Odoribacter/*.csv crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/All_trees/
