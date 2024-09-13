
nextflow.enable.dsl=2
params.fastq_output = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output/"
params.inputFile = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/Name_list.txt"
params.output_concat = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/Concat_refs"
params.main_dir = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping"
params.output_dir_extracted_bam = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/All_extracted_bam"


process CONCATENATE {
    // Set the maximum number of threads to use for this process
    errorStrategy { task.exitStatus == 143 ? 'ignore' : 'terminate' }

    publishDir params.output_concat, mode: 'symlink'

    // 3.1 Define the input file
    input:
    each (samples)

    // 3.2 Define the output file
    output:
    path "*_ALL.fna", emit: concatenate_ref

    // 3.3 Define the script
    script:
    """
    python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts_mapping/Run_fasta_file_maker.py ${params.main_dir}\
     /crex/proj/snic2022-6-144/nobackup/NIKOLAY/GTDB/GTDB_KRAKEN2_DB/genomes_from_gtdbtk/GTDB_KRAKEN2_DB/library/added/\
     /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/Genome_name_files_GTDB.txt\
     /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/Kraken_report_ALL_S_copy.csv\
     ${samples} 
    """
}


process MAP_READS {
    // Set the maximum number of threads to use for this process
    errorStrategy { task.exitStatus == 143 ? 'ignore' : 'terminate' }

    publishDir params.output_dir_extracted_bam, mode: 'symlink'    

    // 3.1 Define the input file
    input:
    each (input_file)

    // 3.2 Define the output file
    output:
    path "tmp/*.tmp", emit: species_to_split

    // 3.3 Define the script
    script:
    """
    filename=\$(basename ${input_file} _ALL.fna)
    bam_filename="${params.fastq_output}\${filename}.unmapped.fastq.gz"
    echo indexing
    bowtie2-build --quiet --threads ${params.threads_MAP_READS}  -f ${input_file} \${filename}
    echo mapping
    bowtie2 --quiet --very-sensitive --threads ${params.threads_MAP_READS} -x \${filename} -U \${bam_filename} | samtools view --verbosity 0 -b -F 4 -@ ${params.threads_MAP_READS} | samtools sort --verbosity 0 -@ ${params.threads_MAP_READS} -O bam -o \${filename}.bam
    echo indexing 
    samtools index -@ ${params.threads_MAP_READS} \${filename}.bam
    echo splitting bam
    python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts_mapping/Create_split_table.py -s ${input_file.baseName} -c ${params.main_dir}/Contig_species_df.csv -o ${input_file.baseName}_taxonomy_to_split.txt
    """
}

process ANCIEN_DAMAGE {

    publishDir params.output_dir_extracted_bam, mode: 'symlink'
   
    errorStrategy { task.exitStatus == 140 ? 'ignore' : 'terminate' }        

    // 3.1 Define the input file
    input:
    each (bam_to_split)

    output:
    path("*_extracted.bam"), emit: splited_bam_files
    path("*_dna_damage.txt"), emit: damage_files
    path("*_dna_damage_plot.pdf"), emit: damage_plots
    path("*.cov"), emit: cov_files
    path("*_extracted.bam.bai"), emit: index_files
    // 3.2 Define the script

    script:
    """ 
    directory=\$(dirname "${bam_to_split}")
    trimmed_directory=\$(echo "\${directory}" | sed "s|tmp.*||")
    ln -s \${trimmed_directory}*.bam .
    ln -s \${trimmed_directory}*.bai .
    Taxonomy_name=\$(basename ${bam_to_split} .tmp)
    python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts_mapping/Split_bam_file.py -t \${Taxonomy_name} -c ${params.main_dir}/Contig_species_df.csv -b *.bam

    Bam_name=\$(find . -type f -name '*_extracted.bam' -exec basename {} \\;)
    
    for file in \$Bam_name; do
       base=\$(basename "\$file" "_extracted.bam")
       samtools sort -@ ${params.threads_ANCIEN_DAMAGE} --verbosity 0 "\$file" -o "sorted_\${base}_extracted.bam" 2> Samtools_sort_log_file.txt
    done
    samtools index -@ ${params.threads_ANCIEN_DAMAGE} sorted_*extracted.bam
    python /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts_mapping/Calculate_dna_damage.py sorted_*extracted.bam
    python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts_mapping/Gather_coverage_from_bam.py sorted_*extracted.bam
    """
}

workflow {
   samples = file(params.inputFile).readLines() 
   
   CONCATENATE(samples)
 
   fasta_files=CONCATENATE.out.concatenate_ref.collect()

   MAP_READS(fasta_files)

   bam_to_split = MAP_READS.out.species_to_split.collect()

   ANCIEN_DAMAGE(bam_to_split)
}
