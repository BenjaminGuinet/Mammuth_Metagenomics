
params.fastq_output = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output/"
params.inputFile = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/Name_list.txt"
params.output_concat = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/Concat_refs"
params.main_dir = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping"
params.output_dir_extracted_bam = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/All_extracted_bam"

process MAP_READS {
    // Set the maximum number of threads to use for this process
    errorStrategy { task.exitStatus == 143 ? 'ignore' : 'terminate' }

    publishDir params.output_dir_extracted_bam, mode: 'symlink'    

    // 3.1 Define the input file
    input:
    each (input_file)

    // 3.2 Define the output file
    output:
    path "*.bam", emit: bam_files

    // 3.3 Define the script
    script:
    """
    filename=\$(basename ${input_file} _ALL.fna)
    bam_filename_sedi="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Blank_contamination_reads/Sediment_metagenomes/ALL_merged.sediment.fastq.gz"
    bam_filename_lab="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Blank_contamination_reads/Laboratory_blanks_reads/ALL_merged.lab.fastq.gz"
    echo LAB BLANK ANALYSIS
    echo indexing
    bowtie2-build --quiet --threads ${params.threads_MAP_READS}  -f ${input_file} \${filename}
    echo mapping
    bowtie2 --quiet --very-sensitive --threads ${params.threads_MAP_READS} -x \${filename} -U \${bam_filename_lab} | samtools view --verbosity 0 -b -F 4 -@ ${params.threads_MAP_READS} | samtools sort --verbosity 0 -@ ${params.threads_MAP_READS} -O bam -o \${filename}_lab.bam 2> Samtools_sort_log_file.txt
    echo indexing 
    samtools index -@ ${params.threads_MAP_READS} \${filename}_lab.bam
    echo SEDIMENT ANALYSIS
    echo indexing
    bowtie2-build --quiet --threads ${params.threads_MAP_READS}  -f ${input_file} \${filename}
    echo mapping
    bowtie2 --quiet --very-sensitive --threads ${params.threads_MAP_READS} -x \${filename} -U \${bam_filename_sedi} | samtools view --verbosity 0 -b -F 4 -@ ${params.threads_MAP_READS} | samtools sort --verbosity 0 -@ ${params.threads_MAP_READS} -O bam -o \${filename}_sed.bam 2> Samtools_sort_log_file.txt
    echo indexing 
    samtools index -@ ${params.threads_MAP_READS} \${filename}_sed.bam
    echo splitting bam
    #python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Scripts_mapping/Create_split_table.py -s ${input_file.baseName} -c ${params.main_dir}/Contig_species_df.csv -o ${input_file.baseName}_taxonomy_to_split.txt
    """
}

workflow {

   fasta_files=reads = Channel.fromPath("${params.output_concat}/*_ALL.fna")

   MAP_READS(fasta_files)

}

