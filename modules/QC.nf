/* =========================================
   Quality Control Processes
   ========================================= */

nextflow.enable.dsl = 2

/* =========================================
   PROCESS: FASTQC - Quality Assessment
   ========================================= */
process FASTQC {

    tag { sample_id }

    publishDir "${params.output_dir}/${params.qc_dir}", mode: 'rellink'
    
    errorStrategy 'finish'

    container "community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29"
    
    input:
    tuple val(sample_id), path(read1), path(read2), val(forward_primer), val(reverse_primer)
    
    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip


    script:
      """
      fastqc ${read1} ${read2}
      """

}

/* =========================================
   PROCESS: MULTIQC - Comprehensive Quality Report
   ========================================= */
process MULTIQC {
    tag { "MultiQC_report" }
    publishDir "${params.output_dir}/${params.multiqc_dir}", mode: 'rellink'
    
    container 'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b'
    
    input:
    path qc_files
    
    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data
    
    script:
    """
    multiqc .
    """
}
