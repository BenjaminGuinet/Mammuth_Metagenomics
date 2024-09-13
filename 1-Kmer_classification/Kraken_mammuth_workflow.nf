// 1. Define the input directory
// params.input_dir = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Data2/"
params.input_dir = "/crex/proj/sllstore2017093/mammoth/final_bams/ancient_dna/"

// 2. Define the output directory
params.output_dir = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output"

params.threads="16"

// 3. Define the process using nextflow.enable.dsl=2

process EXTRACT_UNMAPPED {

    // Set the maximum number of threads to use for this process
    cpus params.threads

    publishDir params.output_dir, mode: 'symlink'

    // Set the maximum concurrent jobs at 100
    maxForks 190

    // 3.1 Define the input file
    input:
    each (input_file)

    // 3.2 Define the output file
    output:
    path("${input_file.baseName}.unmapped.fastq.gz")

    // 3.3 Define the script
    script:
    """
    echo ""
    echo "############################################################################################################"
    echo "###    STEP1"
    echo "###    Extracting unmapped reads of : ${input_file.baseName}                                            "
    echo "###    All output will be written in :  ${params.output_dir}                                             "
    echo "############################################################################################################"
    echo ""
    time_start=\$(date +%s)
    mkdir -p ${params.output_dir}
    echo ""
    samtools fastq -@ ${params.threads} -f 4 ${input_file} | pigz -p ${params.threads} -c > ${input_file.baseName}.unmapped.fastq.gz
    echo "Unmapped reads extracted succesully"
    echo ""
    time_end=\$(date +%s)
    time_diff=\$((time_end - time_start))
    echo "Time spent for extracting unmapped reads of ${input_file.baseName} : \$time_diff seconds"
    echo ""
    """

}

// 4. Run Kraken2 classification in the .unmapped.fastq files
process KRAKEN2 {

    // Set the maximum concurrency to 1
    maxForks 400

    // Set the maximum number of threads to use for this process
    cpus params.threads

    publishDir params.output_dir, mode: 'symlink'

    // 3.1 Define the input file
    input:
    each (input_file) 

    // 3.2 Define the output file
    output:
    tuple path("${input_file.baseName}.kraken2.report"),
    path("${input_file.baseName}_unclass.fq.gz"),
    path("${input_file.baseName}.kraken2.out.gz")


    // 3.3 Define the script
    script:
    """
    echo ""
    echo "############################################################################################################"
    echo "###    STEP2"
    echo "###    Running Kakraken2 on unmapped reads of : ${input_file.baseName}   	                           "
    echo "############################################################################################################"
    echo ""
    time_start=\$(date +%s)
    kraken2 --db /crex/proj/snic2022-6-144/nobackup/BENJAMIN/GTDB_KRAKEN2_DB/ --unclassified-out ${input_file.baseName}_unclass.fq --report-minimizer-data --threads ${params.threads} --report ${input_file.baseName}.kraken2.report --output ${input_file.baseName}.kraken2.out ${input_file}
    pigz  -p ${params.threads} ${input_file.baseName}.kraken2.out
    pigz  -p ${params.threads} ${input_file.baseName}_unclass.fq
    # Remove the uncompressed files
    echo "Kraken2 succesully ran on unmapped reads..."
    echo ""
    time_end=\$(date +%s)
    time_diff=\$((time_end - time_start))
    echo "Time spent for running Kraken2 on unmapped reads of ${input_file.baseName} : \$time_diff seconds"
    """
}


// 4. Run the workflow
workflow {
    unmapped = Channel.fromPath("${params.input_dir}/*.bam")
    EXTRACT_UNMAPPED(unmapped)
    KRAKEN2(EXTRACT_UNMAPPED.output)
}
