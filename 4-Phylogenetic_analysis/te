

conda activate My_conda
module load singularity
module load bowtie2
module load angsd

Sp_name="Pasteurella_mapq30_ALL"
output_main_dir="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/${Sp_name}"

mkdir -p ${output_main_dir}
cd ${output_main_dir}

# first download the assemblies python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Download_NCBI_assemblies.py -t List_assemblies.txt -o .

list_modern_genomes="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Script_pannacota/${Sp_name}_assemblies.txt"
cp ${list_modern_genomes} .


panacota="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/panacota.img"

# Annotate the genomes 
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA annotate -l ${Sp_name}_assemblies.txt -d . -r 2-res-prokka -n GENO  --threads 16

cat 2-res-prokka/Proteins/* >> 2-res-prokka/Proteins/GENO3.All.prt

# Create the pangenome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt /${panacota} PanACoTA pangenome -l ${Sp_name}_assemblies.txt -d 2-res-prokka/Proteins  -i 0.7  -o 3-pangenome -n GENO3 --threads 16

# Create the core and persistent genome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA corepers -p 3-pangenome/PanGenome-GENO3.All.prt-clust-0.7-mode1-th16.lst -o 4-corepers  -t 0.3
#The persistent genome contains 1931 families, each one having exactly 1 member from at least 30.0% of the 23 different genomes (that is 7 genomes). The other genomes are absent from the family
# Align the genomes to the core and persistent genome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA align -c 4-corepers/PersGenome_PanGenome-GENO3.All.prt-clust-0.7-mode1-th16.lst-all_0.3.lst -l 2-res-prokka/LSTINFO-${Sp_name}_assemblies.lst  -n GENO3_1 -d 2-res-prokka -o 5-align --threads 16

# mapp the reads against the closest reference genome 

# Extract only the part of the closest reference genome to mapp the ancient reads
# Define the input file pattern
input_pattern="5-align/Align-GENO3_1/GENO3_1-mafft-prt2nuc.*.aln"

# Loop through each matching file
for file in $input_pattern; do
    # Extract the filename without the path
    base_name=$(basename "$file")
    # Define the output file name
    output_file="5-align/Align-GENO3_1/Isolated_${base_name}"
    # Extract sequences where the header starts with GEN1
    awk '/^>/ {printit = ($0 ~ /^>GEN1/)} printit' "$file" > "$output_file"
    echo "Extracted sequences saved to: $output_file"
done

rm -rf 5-align/All_mapped_reads
mkdir 5-align/All_mapped_reads

fastq_files=(
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Pasteurella_mapq30/Mammuthus-FK033_ALL_new_Bisgaard_extracted_mapped.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Pasteurella_mapq30/Mammuthus-MD024_ALL_new_Bisgaard_extracted_mapped.fastq.gz"
)


output_dir="5-align/All_mapped_reads"
output_dir2="5-align/Align-GENO3_1"


for file in ${output_main_dir}/5-align/Align-GENO3_1/Isolated_GENO3_1-mafft-prt2nuc.*.aln; do
    echo "Processing file: $file"

    ref_genome=$file
    ref_genome_built=$(basename "$ref_genome" .aln)
    
    cd "$output_dir"

    # Create symbolic links for each fastq file
    for fastq_file in "${fastq_files[@]}"; do
        ln -sf "$fastq_file" .
    done

    threads_MAP_READS=4

    # Build Bowtie2 index
    bowtie2-build --quiet --threads $threads_MAP_READS -f "$ref_genome" "$ref_genome_built"
    echo "Mapping..."

    for fastq_file in "${fastq_files[@]}"; do
        fastq_basename=$(basename "$fastq_file")
        bam_filename="${fastq_basename%.fastq.gz}_vs_${ref_genome_built}"

        bowtie2 --very-sensitive --threads $threads_MAP_READS -x "$ref_genome_built" -U "$fastq_basename" | \
            samtools view --verbosity 0 -b -F 4 -q 30 -@ $threads_MAP_READS | \
            samtools sort --verbosity 0 -@ $threads_MAP_READS -O bam -o "${bam_filename}.bam"

        python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Filter_edit_reads_bam.py "${bam_filename}.bam" "${bam_filename}_edit0-2.bam"

        echo "Sorting..."
        samtools sort -@ $threads_MAP_READS -O bam -o "${bam_filename}_sorted.bam" "${bam_filename}.bam"
        samtools sort -@ $threads_MAP_READS -O bam -o "${bam_filename}_sorted_edit0-2.bam" "${bam_filename}_edit0-2.bam"
        
        echo "Indexing..."
        samtools index -@ $threads_MAP_READS "${bam_filename}_sorted.bam"
        samtools index -@ $threads_MAP_READS "${bam_filename}_sorted_edit0-2.bam"
    
        angsd -doFasta 2 -doCounts 1 -minInd 1 -minMapQ 30 -i "${bam_filename}_sorted.bam" -out "${bam_filename}_sorted_for_phylogeny"
        angsd -doFasta 2 -doCounts 1 -minMapQ 30 -minInd 1 -setMinDepth 2 -i "${bam_filename}_sorted_edit0-2.bam" -out "${bam_filename}_sorted_for_phylogeny_edit0-2_minDepth2"
        gunzip "${bam_filename}_sorted_for_phylogeny.fa.gz" --force
        gunzip "${bam_filename}_sorted_for_phylogeny_edit0-2_minDepth2.fa.gz" --force
        

        new_name=$(basename "$fastq_basename" .fastq.gz)
        sed -i "1s/.*/>${new_name}/" "${bam_filename}_sorted_for_phylogeny.fa"
        sed -i "1s/.*/>${new_name}/" "${bam_filename}_sorted_for_phylogeny_edit0-2_minDepth2.fa"
    done
    
    new_file="${file/Isolated_/}"
    cat_file="${file/Isolated_/ALL_cat_}"
    cat_file_edit0_2_minDepth2="${file/Isolated_/ALL_cat_edit0-2_minDepth2_}"
    # remove if file already exists
    rm -f "$cat_file"
    rm -f "$cat_file_edit0_2_minDepth2"
    rm {$output_dir2}*Isolated*
    rm $output_dir2/ALL_cat_*

    passed="NO"

    for fastq_file in "${fastq_files[@]}"; do
        fastq_basename=$(basename "$fastq_file")
        bam_filename="${fastq_basename%.fastq.gz}_vs_${ref_genome_built}"
        
        if [  -s "${bam_filename}_sorted_for_phylogeny.fa" ]; then
            passed="YES"
            sed -i "s/GEN1/${fastq_basename%.fastq.gz}/g" "${bam_filename}_sorted_for_phylogeny.fa"
            sed -i "s/GEN1/${fastq_basename%.fastq.gz}/g" "${bam_filename}_sorted_for_phylogeny_edit0-2_minDepth2.fa"
            cat echo "${bam_filename}_sorted_for_phylogeny.fa" >> "${cat_file}"
            cat "${bam_filename}_sorted_for_phylogeny_edit0-2_minDepth2.fa" >> "${cat_file_edit0_2_minDepth2}"
            echo "Done."
        else
            echo "Skipped: ${bam_filename}_sorted_for_phylogeny.fa is empty."
        fi
    done

    if [ "$passed" == "YES" ]; then
        cat "$new_file" >> "$cat_file"
        cat "$new_file" >> "$cat_file_edit0_2_minDepth2"
    fi

    echo "Done processing file: $file"
done



rm *.bt2
find .r -type f -name "*.bam" ! -name "*_sorted.bam" -exec rm {} +


# Create partitions with files 
cd ../../
mkdir 6-phylogeny

cd 6-phylogeny



python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py ${output_main_dir}/5-align/Align-GENO3_1/ ${output_main_dir}/6-phylogeny/All_partitions_merged.aln ${output_main_dir}/6-phylogeny/All_partitions_merged.part --tag ALL_cat_GENO3
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py ${output_main_dir}/5-align/Align-GENO3_1/ ${output_main_dir}/6-phylogeny/All_partitions_merged_edit0-2_minDepth2.aln ${output_main_dir}/6-phylogeny/All_partitions_merged_edit0-2_minDepth2.part --tag ALL_cat_edit0-2_minDepth2_GENO3

# Update the names 
# Read the file and format the information into key-value pairs
awk -F' :: ' '{print $2"="$1}' "$list_modern_genomes" | sed 's/_genomic\.fna//g' > gene_names.txt
# Loop over each line in the mapping file
while IFS="=" read -r gen value; do
  # Use sed to replace GENx with the corresponding value in the original file
  sed -i "s/$gen/$value/g" All_partitions*
done < gene_names.txt




#!/bin/bash
#SBATCH -A naiss2024-22-84
#SBATCH -p shared
#SBATCH -n 7
#SBATCH -t 10:00:00
#SBATCH -e /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Pasteurella_mapq30_ALL/PHYLOGENY.err
#SBATCH -o /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Pasteurella_mapq30_ALL/PHYLOGENY.out
#SBATCH --job-name PHYLOGENY_Pasteurella

# Run iqtree 
module load iqtree/2.3.5-cpeGNU-23.12

cd cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Pasteurella_mapq30_ALL/6-phylogeny

echo "Done."
rm All_partitions_merged_Sub*
for file in All_partitions*.aln; do
    iqtree2 -s $file -p "${file%.aln}.part"  -m MFP --threads 7 -B 1000 -alrt 1000 -redo
done

# Count the shared sites
for file in All_partitions_merged*.aln; do
    python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Count_shared_sites_MSA.py $file "${file%.aln}.csv"  
done





###########################
# Phylogenetic placement analysis 
##########################


python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Remove_seq_from_MSA.py All_partitions_merged_edit0-2_minDepth2.aln -delete  Mammuthus-FK033_ALL_new_Bisgaard_extracted_mapped Mammuthus-MD024_ALL_new_Bisgaard_extracted_mapped -o No_mammoths_All_partitions_merged_edit0-2_minDepth2.aln
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Remove_seq_from_MSA.py All_partitions_merged_edit0-2_minDepth2.aln -delete  GCA_900454705.1_50465_G01 GCA_020810595.1_ASM2081059v1 GCA_000262245.1_ASM26224v1 GCA_002850605.1_ASM285060v1 GCA_003261335.1_ASM326133v1 GCA_030783905.1_ASM3078390v1 GCA_900638445.1_57675_E01 GCF_900187275.1_51765_F01 GCA_002073255.2_ASM207325v2 GCA_013377295.1_ASM1337729v1 GCA_900454845.1_50465_C01 GCA_000298675.1_P1059v1 GCA_018139065.1_ASM1813906v1 GCA_014338465.1_ASM1433846v1 GCA_003261515.1_ASM326151v1 GCA_900636625.1_43295_B02 GCF_000973525.1_ASM97352v1 GCA_020810675.1_ASM2081067v1 GCA_900186835.1_49950_F01 GCA_003261435.1_ASM326143v1 GCA_011390865.1_ASM1139086v1 GCA_963693435.1_Pasteurella_atlantica_F6K1 GCA_003096995.1_ASM309699v1 -o No_modern_All_partitions_merged_edit0-2_minDepth2.aln

iqtree2 -s No_mammoths_All_partitions_merged_edit0-2_minDepth2.aln -p All_partitions_merged_edit0-2_minDepth2.part -m MFP --threads 14 -B 1000 -alrt 1000 --prefix No_mammoths_All_partitions_merged_edit0-2_minDepth2
rm -f ./epa_info.log
epa-ng -t No_mammoths_All_partitions_merged_edit0-2_minDepth2.treefile  --ref-msa No_mammoths_All_partitions_merged_edit0-2_minDepth2.aln  --query No_modern_All_partitions_merged_edit0-2_minDepth2.aln -m GTR+G
/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/gappa/bin/gappa examine assign --jplace-path epa_result.jplace --taxon-file taxon_file.txt --per-query-results --allow-file-overwriting





# PHYLETPATH ANALYSIS 

# First rename the Tips 
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Rename_phylogeny_and_MSA.py ${output_main_dir}/Genome_metadata.csv ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2.aln ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2.treefile

## create VCF with snp-sites from multiple sequence alignment (MSA) file ###
snp-sites ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.aln -c -v -o ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2.vcf

module load PDC/23.12 R/4.4.1-cpeGNU-23.12
## Bianca’s script to change missing data coding in vcf ###
Rscript  /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/fix_vcf.R ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2.vcf ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_fixed.vcf

Ref_name="Bisgaard_Taxon_45"

## fix naming problem in vcf file (replace 1’s with the name of the closest ref) ###
sed "s/^1\t/${Ref_name}\t/" "${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_fixed.vcf" > "${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_fixed_new_name.vcf"


## Mapping sample reads to closest ref  ###

# first run bowtie2 to map the reads against the reference genome

fastq_files=(
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Pasteurella_mapq30/Mammuthus-FK033_ALL_new_Bisgaard_extracted_mapped.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Pasteurella_mapq30/Mammuthus-MD024_ALL_new_Bisgaard_extracted_mapped.fastq.gz"
)

#extract the Bisgaard Taxon 45 aligned part 
awk -v Ref_name="$Ref_name" '/^>/ {printit = ($0 ~ "^>"Ref_name)} printit'  No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.aln > Target_REF_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.aln

ref_genome_built="Target_REF_No_mammoths_All_partitions_merged_edit0-2_minDepth2"
# create bowtie2 index
bowtie2-build --quiet --threads 8 -f Target_REF_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.aln ${ref_genome_built}

# Mapp each mammoth fastq file to the Bisgaard taxon 45 reference genome
for fastq_file in "${fastq_files[@]}"; do
        fastq_basename=$(basename "$fastq_file")
        bam_filename="${fastq_basename%.fastq.gz}_vs_${ref_genome_built}"
        bowtie2 --very-sensitive --threads 6 -x Target_REF_No_mammoths_All_partitions_merged_edit0-2_minDepth2 -U "$fastq_file" | \
            samtools view --verbosity 0 -b -F 4 -q 30 -@ 6 | \
            samtools sort --verbosity 0 -@ 6 -O bam -o "${bam_filename}.bam"
        # echo the full paht and name of the bam file to bam.list
        echo "$(pwd)/${bam_filename}.bam" >> bam.list
done 


### Run pathPhynder ###

#Edit bootstrap values 

python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Edit_newick_bootstraps.py -input ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile -output ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile
Rscript /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Root_tree.R ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile  ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted 


# 1) Assign informative SNPs to tree branches
/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/phynder -B -o branches.snp ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_fixed_new_name.vcf

# 2) Run pathPhynder to call those SNPs in a given dataset of ancient samples  and find the best path and branch where these can be mapped in the tree
#Prepare data - this will output a bed file for calling variants and tables for pylogenetic placement
#The -G parameter is optional and in this case adds ISOGG haplogroup information to each variant.
module unload R/4.4.1-cpeGNU-23.12
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data:/data /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/pathphynder_b8532c.sif pathPhynder -s prepare -i ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted -p Pasteurella_test -f branches.snp 

# 3) Run pathPhynder best path, call variants, place samples, plot results (the -G can be used to identify haplogroups and it is optional)
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data:/data /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/pathphynder_b8532c.sif pathPhynder  -i ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted   -p tree_data/Pasteurella_test -l bam.list -s all -t 100 -r Target_REF_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.aln

#######################
# MAximums likelihood  #

#convert calls to vcf
Rscript /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/R/make_vcf.R intree_folder/ $Ref_name ancient_calls.vcf

#place samples with phynder
/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/phynder -q ancient_calls.vcf -p 0.01 -o query.phy ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_fixed_new_name.vcf

#plot results
module load R/4.4.1-cpeGNU-23.12
Rscript /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/R/plot_likes.R ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted query.phy results_folder






import pysam
import pandas as pd

# Load the BAM file

df = pd.read_csv("coverage.txt", sep="\t", header=None)
# keep only column 2 > 0
df = df[df[2] > 0]
print(df)

# sort by Position 

df.sort_values(by=['Position'], inplace=True)

tab = pd.read_csv("tree_data/Pasteurella_test.sites.bed", sep="\t", header=None)
tab[0] = "CA_030783905.1_ASM3078390v1"
# Display the DataFrame
print(df)
# save the file tab
tab.to_csv("tree_data/New_Pasteurella_test.sites.bed", sep="\t", header=False, index=False)

# output positions in df are also in tab[1]
df.loc[df[1].isin(tab[1])]

tab2 = pd.read_csv("Output_mpileup.txt", sep="\t", header=None)


tab2.loc[tab2[1].isin(df[1])]









# 1) Assign SNPs to branches of the tree, as above. Skip if you have done this before.
/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/phynder -B -o branches.snp /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data/BigTree_Y/bigtree_annotated_V1.nwk /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data/BigTree_Y/BigTree.Y.201219.vcf.gz

# 2) Prepare sites (writes bed files for variant calling and other files for phylogenetic placement). Skip if you have done this before.
#The -G parameter is optional and in this case adds ISOGG haplogroup information to each variant.
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data:/data /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/pathphynder_b8532c.sif pathPhynder -s prepare -i /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data/BigTree_Y/bigtree_annotated_V1.nwk -p BigTree_Y_data -f branches.snp 



# 3) Run pathPhynder best path, call variants, place samples, plot results (the -G can be used to identify haplogroups and it is optional)
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data:/data /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/pathphynder_b8532c.sif pathPhynder  -i /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data/BigTree_Y/bigtree_annotated_V1.nwk   -p tree_data/BigTree_Y_data -l bam.list -s all -t 100 





#convert calls to vcf
Rscript /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/R/make_vcf.R tree_data chrY ancient_calls.vcf

#place samples with phynder
/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/phynder -q ancient_calls.vcf -p 0.01 -o query.phy ../data/BigTree_Y/bigtree_annotated_V1.nwk ../data/BigTree_Y/BigTree.Y.201219.vcf.gz

#plot results
Rscript /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/R/plot_likes.R /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data/BigTree_Y/bigtree_annotated_V1.nwk query.phy results_folder



/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/phynder -q ${output_main_dir}/6-phylogeny/All_partitions_merged_edit0-2_minDepth2.vcf -p 0.01 -o query.phy ${output_main_dir}/6-phylogeny/All_partitions_merged_edit0-2_minDepth2.part_bootstrap_edited.treefile ${output_main_dir}/6-phylogeny/All_partitions_merged_edit0-2_minDepth2.vcf






# compute individual phylogenies

# first create the sub files in order to keep the partitions and remove the N sites in the targeted species
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Extract_sequences_from_MSA2.py -o .\
    --species Mammuthus-FK033_ALL_new_Bisgaard_extracted_mapped Mammuthus-MD024_ALL_new_Bisgaard_extracted_mapped\
     --alignment_dir /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Pasteurella_mapq30_ALL/5-align/Align-GENO3_1/


# Now we will merge all these files 
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_Mammuthus-MD024.aln All_partitions_merged_Sub_Mammuthus-MD024.part --tag  Sub_Mammuthus-MD024
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_Mammuthus-FK033.aln All_partitions_merged_Sub_Mammuthus-FK033.part  --tag  Sub_Mammuthus-FK033


while IFS="=" read -r gen value; do
  # Use sed to replace GENx with the corresponding value in the original file
  sed -i "s/$gen/$value/g" All_partitions_merged_Sub*
done < gene_names.txt

# remove useless files
rm Sub_*

for file in All_partitions_merged_Sub*.aln; do
    iqtree2 -s $file -p "${file%.aln}.part"  -m MFP --threads 7 -B 1000 -alrt 1000 
done



# Plot edit distance and read length distribution
mkdir 7-plots
cd 7-plots


samtools merge -@ 8 Pasteurella_Mammuthus-FK033.bam ../5-align/All_mapped_reads/Mammuthus-FK033*sorted.bam
samtools merge -@ 8 Pasteurella_Mammuthus-MD024.bam ../5-align/All_mapped_reads/Mammuthus-MD024*sorted.bam
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Plot_edit_distance_and_percidentity3.py --bam_files Pasteurella_Mammuthus-FK033.bam Pasteurella_Mammuthus-MD024.bam  --out Pasteurella_read_stats.pdf

