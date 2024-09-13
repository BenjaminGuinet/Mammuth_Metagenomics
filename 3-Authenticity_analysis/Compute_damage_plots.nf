params.output_dir_extracted_bam = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/All_extracted_bam"
params.output_dir_extracted_damage = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/All_extracted_damage_plots"
// params.inputFile="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/All_extracted_bam/All_sorted_bam_files.txt"
// params.inputFile="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/All_extracted_bam/All_sorted_bam_files_sub.txt"
params.inputFile="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_mapping/All_extracted_bam/missing_files.csv"
process DAMAGE_PLOTS {
    // Set the maximum number of threads to use for this process
    errorStrategy 'ignore'

    publishDir params.output_dir_extracted_damage, mode: 'symlink'    

    // 3.1 Define the input file
    input:
    each (input_file)

    // 3.2 Define the output file
    output:
    path "*dna_damage.csv", emit: damage_table
    path "*dna_damage_plots.pdf", emit: damage_plot
    // 3.3 Define the script
    script:
    """
    cp ${input_file} .
    filename=\$(basename ${input_file})
    python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/TOOLS/Make_multiple_plotV2.py -b \${filename} &> Log_file.txt
    rm \${filename}
    """
}

workflow {
   bam_files=reads = file(params.inputFile).readLines()
   DAMAGE_PLOTS(bam_files)
}

