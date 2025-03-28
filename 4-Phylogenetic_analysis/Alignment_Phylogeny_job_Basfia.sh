conda activate My_conda
module load singularity
module load bowtie2
module load angsd

Sp_name="Basfia_mapq30_ALL"
Sp_name2="Basfia"
output_main_dir="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/${Sp_name}"

mkdir -p ${output_main_dir}
cd ${output_main_dir}

# first download the assemblies python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Download_NCBI_assemblies.py -t List_assemblies.txt -o .

list_modern_genomes="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Script_pannacota/${Sp_name}_assemblies.txt"
cp ${list_modern_genomes} .

panacota="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/panacota.img"

# Annotate the genomes 
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA annotate -l ${Sp_name}_assemblies.txt -d . -r 2-res-prokka -n GENO  --threads 30

cat 2-res-prokka/Proteins/* >> 2-res-prokka/Proteins/GENO3.All.prt

# Create the pangenome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt /${panacota} PanACoTA pangenome -l ${Sp_name}_assemblies.txt -d 2-res-prokka/Proteins  -i 0.7  -o 3-pangenome -n GENO3 --threads 30

# Create the core and persistent genome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA corepers -p 3-pangenome/PanGenome-GENO3.All.prt-clust-0.7-mode1-th30.lst -o 4-corepers  -t 0.3
# The persistent genome contains 2063 families, each one having exactly 1 member from at least 30.0% of the 15 different genomes (that is 5 genomes). The other genomes are absent from the family.
# Align the genomes to the core and persistent genome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA align -c 4-corepers/PersGenome_PanGenome-GENO3.All.prt-clust-0.7-mode1-th30.lst-all_0.3.lst -l 2-res-prokka/LSTINFO-${Sp_name}_assemblies.lst  -n GENO3_1 -d 2-res-prokka -o 5-align --threads 20

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

# remove if file already exists
rm -f "5-align/All_mapped_reads/"
mkdir 5-align/All_mapped_reads

fastq_files=(
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-MD228_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L414_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L423_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L389_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L426_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L418_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-FK012_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-M13_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-M40_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
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
#SBATCH -e /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Basfia_mapq30_ALL/PHYLOGENY.err
#SBATCH -o /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Basfia_mapq30_ALL/PHYLOGENY.out
#SBATCH --job-name PHYLOGENY_Basfia

# Run iqtree 
module load iqtree/2.3.5-cpeGNU-23.12

cd cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Basfia_mapq30_ALL/6-phylogeny
# 
rm All_partitions_merged_Sub*
for file in *.aln; do
    iqtree2 -s $file -p "${file%.aln}.part"  -m MFP --threads 7 -B 1000 -alrt 1000 
done


# Count the shared sites
for file in All_partitions_merged*.aln; do
    python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Count_shared_sites_MSA.py $file "${file%.aln}.csv"  
done


# Phylogenetic placement analysis 
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Remove_seq_from_MSA.py All_partitions_merged_edit0-2_minDepth2.aln -delete  sorted_Mammuthus-M40_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-FK012_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L418_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-MD228_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L389_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L423_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L414_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L426_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-M13_Basfia_A_sp014843355_g1GiUvRh2j_extracted -o No_mammoths_All_partitions_merged_edit0-2_minDepth2.aln
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Remove_seq_from_MSA.py All_partitions_merged_edit0-2_minDepth2.aln -delete  GCA_003101015.1_ASM310101v1 GCA_901764995.1_5970_8_10_v2 GCA_900110905.1_IMG-taxon_2593339215_annotated_assembly GCA_900445015.1_55147_H01 GCA_000739435.1_ASM73943v1 GCA_000175195.1_ASM17519v1 GCA_900638445.1_57675_E01 GCA_900460625.1_49642_E01 GCA_004519515.1_ASM451951v1 GCA_000374285.1_ASM37428v1 GCA_000801145.1_ASM80114v1 GCA_000931575.1_ASM93157v1 GCF_011455875.1_ASM1145587v1 GCF_014843355.1_ASM1484335v1 GCA_900101275.1_IMG-taxon_2593339214_annotated_assembly -o No_modern_All_partitions_merged_edit0-2_minDepth2.aln

iqtree2 -s No_mammoths_All_partitions_merged_edit0-2_minDepth2.aln -p All_partitions_merged_edit0-2_minDepth2.part -m MFP --threads 20 -B 1000 -alrt 1000 --prefix No_mammoths_All_partitions_merged_edit0-2_minDepth2
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

Ref_name="Erysipelothrix_tonsillarum_DSM_14972"

## fix naming problem in vcf file (replace 1’s with the name of the closest ref) ###
sed "s/^1\t/${Ref_name}\t/" "${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_fixed.vcf" > "${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_fixed_new_name.vcf"


## Mapping sample reads to closest ref  ###

# first run bowtie2 to map the reads against the reference genome

fastq_files=(
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-MD228_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L414_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L423_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L389_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L426_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-L418_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-FK012_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-M13_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Basfia/sorted_Mammuthus-M40_Basfia_A_sp014843355_g1GiUvRh2j_extracted.fastq.gz"
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
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data:/data /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/pathphynder_b8532c.sif pathPhynder -s prepare -i ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted  -p Pathfinder_dir -f branches.snp 

# 3) Run pathPhynder best path, call variants, place samples, plot results (the -G can be used to identify haplogroups and it is optional)
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/data:/data /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/pathphynder_b8532c.sif pathPhynder  -i ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted  -p tree_data/Pathfinder_dir -l bam.list -s all -t 100 -r Target_REF_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.aln

#######################
# MAximums likelihood  #

#convert calls to vcf
Rscript /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/R/make_vcf.R intree_folder/ ${Ref_name} ancient_calls.vcf

#place samples with phynder
/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/phynder -q ancient_calls.vcf -p 0.01 -o query.phy ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted ${output_main_dir}/6-phylogeny/No_mammoths_All_partitions_merged_edit0-2_minDepth2_fixed_new_name.vcf

#plot results
module load R/4.4.1-cpeGNU-23.12
Rscript /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/phynder/pathPhynder/R/plot_likes.R ${output_main_dir}/6-phylogeny/Bootstrap_edited_No_mammoths_All_partitions_merged_edit0-2_minDepth2_NewTipName.treefile_rooted query.phy results_folder












# compute individual phylogenies

# first create the sub files in order to keep the partitions and  remove the N sites in the targeted species
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Extract_sequences_from_MSA2.py -o .\
    --species sorted_Mammuthus-L423_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L426_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L418_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-FK012_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-M13_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-M40_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-MD228_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L389_Basfia_A_sp014843355_g1GiUvRh2j_extracted sorted_Mammuthus-L414_Basfia_A_sp014843355_g1GiUvRh2j_extracted\
     --alignment_dir /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Basfia_mapq30_ALL/5-align/Align-GENO3_1/


# Now we will merge all these files 
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-L414.aln All_partitions_merged_Sub_sorted_Mammuthus-L414.part --tag  Sub_sorted_Mammuthus-L414
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-L423.aln All_partitions_merged_Sub_sorted_Mammuthus-L423.part --tag  Sub_sorted_Mammuthus-L423
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-L426.aln All_partitions_merged_Sub_sorted_Mammuthus-L426.part --tag  Sub_sorted_Mammuthus-L426
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-L418.aln All_partitions_merged_Sub_sorted_Mammuthus-L418.part --tag  Sub_sorted_Mammuthus-L418
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-FK012.aln All_partitions_merged_Sub_sorted_Mammuthus-FK012.part --tag  Sub_sorted_Mammuthus-FK012
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-M13.aln All_partitions_merged_Sub_sorted_Mammuthus-M13.part --tag  Sub_sorted_Mammuthus-M13
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-M40.aln All_partitions_merged_Sub_sorted_Mammuthus-M40.part --tag  Sub_sorted_Mammuthus-M40
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-MD228.aln All_partitions_merged_Sub_sorted_Mammuthus-MD228.part --tag  Sub_sorted_Mammuthus-MD228
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-L389.aln All_partitions_merged_Sub_sorted_Mammuthus-L389.part --tag  Sub_sorted_Mammuthus-L389


while IFS="=" read -r gen value; do
  # Use sed to replace GENx with the corresponding value in the original file
  sed -i "s/$gen/$value/g" All_partitions_merged_Sub_sorted*
done < gene_names.txt

# remove useless files
rm Sub*

for file in All_partitions_merged_Sub*.aln; do
    iqtree2 -s $file -p "${file%.aln}.part"  -m MFP --threads 7 -B 1000 -alrt 1000 
done




mkdir 7-plots
cd 7-plots

samtools merge -@ 8 Basfia_Mammuthus-L414.bam ../5-align/All_mapped_reads/sorted_Mammuthus-L414*sorted.bam -f
samtools merge -@ 8 Basfia_Mammuthus-L423.bam ../5-align/All_mapped_reads/sorted_Mammuthus-L423*sorted.bam -f 
samtools merge -@ 8 Basfia_Mammuthus-L426.bam ../5-align/All_mapped_reads/sorted_Mammuthus-L426*sorted.bam -f 
samtools merge -@ 8 Basfia_Mammuthus-L418.bam ../5-align/All_mapped_reads/sorted_Mammuthus-L418*sorted.bam -f 
samtools merge -@ 8 Basfia_Mammuthus-FK012.bam ../5-align/All_mapped_reads/sorted_Mammuthus-FK012*sorted.bam -f 
samtools merge -@ 8 Basfia_Mammuthus-M13.bam ../5-align/All_mapped_reads/sorted_Mammuthus-M13*sorted.bam -f 
samtools merge -@ 8 Basfia_Mammuthus-M40.bam ../5-align/All_mapped_reads/sorted_Mammuthus-M40*sorted.bam -f 
samtools merge -@ 8 Basfia_Mammuthus-MD228.bam ../5-align/All_mapped_reads/sorted_Mammuthus-MD228*sorted.bam -f 
samtools merge -@ 8 Basfia_Mammuthus-L389.bam ../5-align/All_mapped_reads/sorted_Mammuthus-L389*sorted.bam -f 


python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Plot_edit_distance_and_percidentity3.py --bam_files Basfia_Mammuthus-L414.bam Basfia_Mammuthus-L423.bam Basfia_Mammuthus-L426.bam Basfia_Mammuthus-L418.bam Basfia_Mammuthus-FK012.bam Basfia_Mammuthus-M13.bam Basfia_Mammuthus-M40.bam Basfia_Mammuthus-MD228.bam Basfia_Mammuthus-L389.bam --out Basfia_read_stats.pdf



