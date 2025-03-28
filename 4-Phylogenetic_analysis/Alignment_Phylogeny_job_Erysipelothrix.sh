


conda activate My_conda
module load singularity
module load bowtie2
module load angsd


Sp_name="Erysipelothrix"
output_main_dir="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/${Sp_name}"

mkdir -p ${output_main_dir}
cd ${output_main_dir}

list_modern_genomes="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Script_pannacota/${Sp_name}_assemblies.txt"
cp  ${list_modern_genomes} .

panacota="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/panacota.img"

# Annotate the genomes 
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA annotate -l ${Sp_name}_assemblies.txt -d . -r 2-res-prokka -n GENO  --threads 30

cat 2-res-prokka/Proteins/* >> 2-res-prokka/Proteins/GENO3.All.prt

# Create the pangenome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt /${panacota} PanACoTA pangenome -l ${Sp_name}_assemblies.txt -d 2-res-prokka/Proteins  -i 0.8  -o 3-pangenome -n GENO3 --threads 30

# Create the core and persistent genome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA corepers -p 3-pangenome/PanGenome-GENO3.All.prt-clust-0.8-mode1-th30.lst -o 4-corepers  -t 0.5

# Align the genomes to the core and persistent genome
singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt ${panacota} PanACoTA align -c 4-corepers/PersGenome_PanGenome-GENO3.All.prt-clust-0.8-mode1-th30.lst-all_0.5.lst -l 2-res-prokka/LSTINFO-${Sp_name}_assemblies.lst  -n GENO3_1 -d 2-res-prokka -o 5-align --threads 20

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

mkdir 5-align/All_mapped_reads

fastq_files=(
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Erysipelothrix/sorted_Mammuthus-FK001_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Erysipelothrix/sorted_Mammuthus-MD024_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Erysipelothrix/sorted_Mammuthus-P033_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Erysipelothrix/sorted_Mammuthus-G.ERR2260501_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Erysipelothrix/sorted_Mammuthus-MD213_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted.fastq.gz"	
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Adycha_sample/Outputs/to_keep/sortedP033_nonUDG_treated_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted_clip.fastq.gz"
    "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_ali_and_phylo_animals/Erysipelothrix/sorted_Mammuthus-P033_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted_merged_P033_nonUDG.fastq.gz"
)

output_dir="5-align/All_mapped_reads"

for file in ${output_main_dir}/5-align/Align-GENO3_1/Isolated_GENO3_1-mafft-prt2nuc.*.aln; do
    echo "Processing file: $file"

    ref_genome=$file
    ref_genome_built=$(basename "$ref_genome" .aln)
    
    cd "$output_dir"

    # Create symbolic links for each fastq file
    for fastq_file in "${fastq_files[@]}"; do
        ln -s "$fastq_file" .
    done

    threads_MAP_READS=16

    # Build Bowtie2 index
    bowtie2-build --quiet --threads $threads_MAP_READS -f "$ref_genome" "$ref_genome_built"
    echo "Mapping..."

    for fastq_file in "${fastq_files[@]}"; do
        fastq_basename=$(basename "$fastq_file")
        bam_filename="${fastq_basename%.fastq.gz}_vs_${ref_genome_built}"

        bowtie2 --very-sensitive --threads $threads_MAP_READS -x "$ref_genome_built" -U "$fastq_basename" | \
            samtools view --verbosity 0 -b -F 4 -q 25 -@ $threads_MAP_READS | \
            samtools sort --verbosity 0 -@ $threads_MAP_READS -O bam -o "${bam_filename}.bam"

        echo "Sorting..."
        samtools sort -@ $threads_MAP_READS -O bam -o "${bam_filename}_sorted.bam" "${bam_filename}.bam"
        
        echo "Indexing..."
        samtools index -@ $threads_MAP_READS "${bam_filename}_sorted.bam"
    
        angsd -doFasta 2 -doCounts 1 -minInd 1 -minMapQ 25 -i "${bam_filename}_sorted.bam" -out "${bam_filename}_sorted_for_phylogeny"
        gunzip "${bam_filename}_sorted_for_phylogeny.fa.gz" --force
        
        new_name=$(basename "$fastq_basename" .fastq.gz)
        sed -i "1s/.*/>${new_name}/" "${bam_filename}_sorted_for_phylogeny.fa"
    done
    
    new_file="${file/Isolated_/}"
    cat_file="${file/Isolated_/ALL_cat_}"

    passed="NO"

    for fastq_file in "${fastq_files[@]}"; do
        fastq_basename=$(basename "$fastq_file")
        bam_filename="${fastq_basename%.fastq.gz}_vs_${ref_genome_built}"
        
        if [  -s "${bam_filename}_sorted_for_phylogeny.fa" ]; then
            passed="YES"
            sed -i "s/GEN1/${fastq_basename%.fastq.gz}/g" "${bam_filename}_sorted_for_phylogeny.fa"
            cat "${bam_filename}_sorted_for_phylogeny.fa" >> "${cat_file}"
            echo "Done."
        else
            echo "Skipped: ${bam_filename}_sorted_for_phylogeny.fa is empty."
        fi
    done

    if [ "$passed" == "YES" ]; then
        cat "$new_file" >> "$cat_file"
    fi

    echo "Done processing file: $file"
done

rm *.bt2
find .r -type f -name "*.bam" ! -name "*_sorted.bam" -exec rm {} +


# Create partitions with files 
cd ../../
mkdir 6-phylogeny

cd 6-phylogeny

python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py ${output_main_dir}/5-align/Align-GENO3_1/ ${output_main_dir}/6-phylogeny/All_partitions_merged_filtred.aln\
 ${output_main_dir}/6-phylogeny/All_partitions_merged_filtred.part\
  --species_pair  sorted_Mammuthus-P033_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted sorted_Mammuthus-MD024_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted --tag ALL_cat_
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py ${output_main_dir}/5-align/Align-GENO3_1/ ${output_main_dir}/6-phylogeny/All_partitions_merged.aln ${output_main_dir}/6-phylogeny/All_partitions_merged.part --tag ALL_cat_

python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py ${output_main_dir}/5-align/Align-GENO3_1/ ${output_main_dir}/6-phylogeny/All_partitions_merged_filtred_full.aln\
 ${output_main_dir}/6-phylogeny/All_partitions_merged_filtred_full.part\
  --species_pair  sorted_Mammuthus-MD024_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted sorted_Mammuthus-P033_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted_merged_P033_nonUDG sorted_Mammuthus-FK001_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted --tag ALL_cat_
# in this case we could not use Mammuthus-G.ERR2260501 and Mammuthus-MD213 as filter as they shared to less sites with the other mammuths

# Update the names 
# Read the file and format the information into key-value pairs
awk -F' :: ' '{print $2"="$1}' "$list_modern_genomes" | sed 's/_genomic\.fna//g' > gene_names.txt
# Loop over each line in the mapping file
while IFS="=" read -r gen value; do
  # Use sed to replace GENx with the corresponding value in the original file
  sed -i "s/$gen/$value/g" All_partitions_merged*
done < gene_names.txt

# Run iqtree 
module load iqtree/2.3.5-cpeGNU-23.12

for file in All_partitions_merged*.aln; do
    iqtree2 -s $file -p "${file%.aln}.part"  -m MFP --threads 7 -B 1000 -alrt 1000 
done

# Count the shared sites
for file in All_partitions_merged*.aln; do
    python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Count_shared_sites_MSA.py $file "${file%.aln}.csv"  
done


# compute individual phylogenies

# first create the sub files in order to keep the partitions and  remove the N sites in the targeted species
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Extract_sequences_from_MSA2.py -o .\
    --species sorted_Mammuthus-FK001_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted sorted_Mammuthus-G sorted_Mammuthus-MD213_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted sortedP033_nonUDG_treated_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted_clip sorted_Mammuthus-MD024_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted sorted_Mammuthus-P033_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted\
     --alignment_dir /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Meta_mammuth_project/Test_folder/Output_pannacota/Erysipelothrix/5-align/Align-GENO3_1/

sorted_Mammuthus-P033_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted
sortedP033_nonUDG_treated_Erysipelothrix_tonsillarum_ZKselpDrZ__extracted_clip
# Now we will merge all these files 
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-FK001.aln All_partitions_merged_Sub_sorted_Mammuthus-FK001.part --tag  Sub_sorted_Mammuthus-FK001
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-G.aln All_partitions_merged_Sub_sorted_Mammuthus-G.part  --tag  Sub_sorted_Mammuthus-G
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-MD213.aln All_partitions_merged_Sub_sorted_Mammuthus-MD213.part  --tag  Sub_sorted_Mammuthus-MD213
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_sortedP033_nonUDG_treated.aln All_partitions_merged_Sub_sorted_sortedP033_nonUDG_treated.part  --tag  Sub_sorted_sortedP033_nonUDG_treated
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-MD024.aln All_partitions_merged_Sub_sorted_Mammuthus-MD024.part  --tag  Sub_sorted_Mammuthus-MD024
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py . All_partitions_merged_Sub_sorted_Mammuthus-P033.aln All_partitions_merged_Sub_sorted_Mammuthus-P033.part  --tag  Sub_sorted_Mammuthus-P033

while IFS="=" read -r gen value; do
  # Use sed to replace GENx with the corresponding value in the original file
  sed -i "s/$gen/$value/g" All_partitions_merged_Sub_sorted*
done < gene_names.txt

# remove useless files
rm Sub_sorted*

for file in All_partitions_merged_Sub*.aln; do
    iqtree2 -s $file -p "${file%.aln}.part"  -m MFP --threads 7 -B 1000 -alrt 1000 
done
