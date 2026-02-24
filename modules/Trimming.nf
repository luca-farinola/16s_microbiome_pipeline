/* =========================================
   Trimming with Cutadapt
   ========================================= */

nextflow.enable.dsl = 2

/* =========================================
   PROCESS: CUTADAPT - Primer Removal and Quality Filtering
   ========================================= */
process CUTADAPT {
    tag { sample_id }
    publishDir "${params.output_dir}/${params.CUTADAPTtrim_dir}", mode: 'rellink'
    
    container 'community.wave.seqera.io/library/cutadapt:5.2--6b9a64599468f0c4'
    
    input:
    tuple val(sample_id), path(read1), path(read2), val(forward_primer), val(reverse_primer)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: fastq
    path "${sample_id}.trim.log", emit: log
    
    // e = maximum fraction of mismatches allowed Primer Length × Error Rate (e)
    // O = minimum overlap length between the read and the primer sequence
    // --rc: also search for the reverse complement of the primer sequences
    // --minimum-length: discard reads shorter than this length after trimming

    script:
    """
    cutadapt -g ${forward_primer} -G ${reverse_primer} \
        --cores ${task.cpus} \
        --rc \
        -e 0.1 \
        -O 3 \
        --minimum-length 1 \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -p ${sample_id}_R2_trimmed.fastq.gz \
        ${read1} ${read2} > ${sample_id}.trim.log
    """
}

process DADA2_FILTNTRIM {
    tag { sample_id }
    label 'process_low'
    publishDir "${params.output_dir}/${params.DADA2trim_dir}", mode: 'rellink'

    conda "bioconda::bioconductor-dada2=1.34.0 bioconda::bioconductor-biostrings conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_1.filt.fastq.gz"), path("${sample_id}_2.filt.fastq.gz"), emit: fastq
    path "${sample_id}.filter_stats.tsv", emit: stats

    // truncLen: truncate reads to this length (default: no truncation). 
    // Reads shorter than this after truncation will be discarded. 
    // This should be set to the position where quality drops off, typically determined by inspecting quality profiles

    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    out <- filterAndTrim(
        "${read1}", "${sample_id}_1.filt.fastq.gz",
        "${read2}", "${sample_id}_2.filt.fastq.gz",
        compress = TRUE,
        multithread = ${task.cpus},
        verbose = TRUE
    )

    # Add sample ID to output
    out <- cbind(out, sample_id = "${sample_id}")

    # Write statistics
    write.table(out, file = "${sample_id}.filter_stats.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
    """
}