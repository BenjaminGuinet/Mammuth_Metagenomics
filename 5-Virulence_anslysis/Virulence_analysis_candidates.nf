nextflow.enable.dsl=2

params.output_dir = "/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_virulence/"
params.input_fasta_dir ="/crex/proj/snic2022-6-144/nobackup/BENJAMIN/Meta_mammuth_project/Test_folder/Output_virulence/Mappings/"

process RUN_MMSEQS {
    // Set the maximum number of threads to use for this process

    // errorStrategy 'ignore'

    publishDir params.output_dir, mode: 'symlink'

    input:
    each (input_fasta) 

    // 3.2 Define the output file
    output:
    path "*.m8", emit: mmseqs_search_result
    // path "*.csv", emit: filtred_mmseqs_search_result

    // 3.3 Define the script
    script:
    """
    input_fasta_basename=\$(basename ${input_fasta} .fa)
    mmseqs createdb ${input_fasta} \${input_fasta_basename}_db
    ln -s ${params.output_dir}VFDB_setB_nt_db* .
    mmseqs search \${input_fasta_basename}_db VFDB_setB_nt_db \${input_fasta_basename}_vs_VFDB_setB_nt_result \${input_fasta_basename}_vs_VFDB_setB_nt_tmp -s 7 --search-type 3 --threads ${params.threads_RUN_MMSEQS} -a 
    mmseqs convertalis \${input_fasta_basename}_db VFDB_setB_nt_db \${input_fasta_basename}_vs_VFDB_setB_nt_result \${input_fasta_basename}_vs_VFDB_setB_nt_result.m8 --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,qcov,tcov,evalue,bits,qlen,tlen,qaln,taln' --format-mode 4 --search-type 3
    #python3 /crex/proj/snic2022-6-144/nobackup/BENJAMIN/TOOLS/Filter_blast.py -b \${input_fasta_basename}_vs_VFDB_setB_nt_result.m8
    """
}

workflow {
    input_fasta_ch = Channel
        .fromPath("${params.input_fasta_dir}/*.fa") # this correspondes to the fasta genome conscensus generated with ANGSD for each clade
        .toList()
    RUN_MMSEQS(input_fasta_ch)
}
