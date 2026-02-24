#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import only needed modules
include { FASTQC } from './modules/QC.nf'
include { MULTIQC } from './modules/QC.nf'
include { CUTADAPT } from './modules/Trimming.nf'
include { DADA2_FILTNTRIM } from './modules/Trimming.nf'
include { DADA2_ERR } from './modules/Infer_ASVs.nf'
include { DADA2_DENOISING } from './modules/Infer_ASVs.nf'
include { DADA2_MERGE } from './modules/Infer_ASVs.nf'
include { DADA2_SEQTABLE } from './modules/Infer_ASVs.nf'
include { KRAKEN2_TAXONOMY; KRAKEN2_TABLE } from './modules/Taxonomy_Kraken_SILVA.nf'
include { DADA2_ASSIGNTAXONOMY; DADA2_TAX_TABLE } from './modules/Taxonomy_Dada2_SILVA.nf'

// Create input channel from CSV
Channel
    .fromPath(params.sample_csv)
    .splitCsv(header: true)
    .map { row -> 
        [row.sample_id, file(row.read1_path), file(row.read2_path), row.forward_primer, row.reverse_primer]
    }
    .set { raw_reads_ch }

workflow {
    
    // QC and Trimming
    FASTQC(raw_reads_ch)
    CUTADAPT(raw_reads_ch)
    MULTIQC( FASTQC.out.zip.collect())
    
    // DADA2 workflow - ASVs
    DADA2_FILTNTRIM(CUTADAPT.out.fastq)
    DADA2_ERR(DADA2_FILTNTRIM.out.fastq)
    filtered_reads_ch = DADA2_FILTNTRIM.out.fastq
    errF_ch = DADA2_ERR.out.errF
    errR_ch = DADA2_ERR.out.errR

    filtered_reads_ch
        .join(errF_ch)
        .join(errR_ch)
        .set { denoising_input_ch }

    DADA2_DENOISING(denoising_input_ch)

    dadaF_ch = DADA2_DENOISING.out.dadaF
    dadaR_ch = DADA2_DENOISING.out.dadaR

    filtered_reads_ch
        .join(dadaF_ch)
        .join(dadaR_ch)
        .set { merge_input_ch }

    DADA2_MERGE(merge_input_ch)

    // After DADA2_MERGE
    merged_ch = DADA2_MERGE.out.merged

    // Extract only the merged file paths and collect into a list
    merged_ch.map { it[1] }.collect().set { merged_files_ch }

    DADA2_SEQTABLE(merged_files_ch)

    KRAKEN2_TAXONOMY(
        DADA2_SEQTABLE.out.asv_fasta,
        params.kraken2_db
    )

    KRAKEN2_TABLE(
        KRAKEN2_TAXONOMY.out.kraken_output,
        KRAKEN2_TAXONOMY.out.kraken_report,
        DADA2_SEQTABLE.out.asvtable
    )

    DADA2_ASSIGNTAXONOMY(
        DADA2_SEQTABLE.out.asv_fasta,
        params.silva_train_set,
        params.silva_species_set
    )
}
                                                          