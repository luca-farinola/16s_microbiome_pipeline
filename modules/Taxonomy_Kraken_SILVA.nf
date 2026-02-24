/* =========================================
   Taxonomy Assignment with Kraken2 + SILVA 138
   ========================================= */

nextflow.enable.dsl = 2

/* =========================================
   PROCESS: KRAKEN2_TAXONOMY
   Assign taxonomy to ASV sequences using Kraken2 and SILVA 138 database
   ========================================= */

process KRAKEN2_TAXONOMY {

    tag "Kraken2 Taxonomy Assignment"
    
    publishDir "${params.output_dir}/${params.taxonomy_dir_kraken}", mode: 'copy'
    
    container "community.wave.seqera.io/library/kraken2:2.17.1--11b4c2ed8e92f26c"

    input:
    path(asv_fasta)
    path(kraken2_db)

    output:
    path("kraken2_output.txt"), emit: kraken_output
    path("kraken2_report.txt"), emit: kraken_report

    script:
    """
    # Run Kraken2 on ASV sequences
    kraken2 --db ${kraken2_db} \
            --output kraken2_output.txt \
            --report kraken2_report.txt \
            --threads ${task.cpus} \
            ${asv_fasta}
    """
}

process KRAKEN2_TABLE {

    tag "Create Kraken2 Abundance Table"
    
    publishDir "${params.output_dir}/${params.taxonomy_dir_kraken}", mode: 'copy'
    
    container "community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186"

    input:
    path(kraken_output)
    path(kraken_report)
    path(asvtable_csv)

    output:
    path("kraken2_taxonomy.txt"), emit: asv_tax

    script:
    """
    Rscript - <<'EOF'
    library(dplyr)
    library(tidyr)
    library(tibble)
    
    # Read Kraken2 output (format: C/U, seqID, taxid, seqlen, kmers)
    kraken_output <- read.table("${kraken_output}", sep="\\t", header=FALSE, stringsAsFactors=FALSE)
    colnames(kraken_output) <- c("classified", "asv_id", "taxid", "seqlen", "kmers")
    
    # Read Kraken2 report (format: percent, clade_reads, direct_reads, rank, taxid, name)
    kraken_report <- read.table("${kraken_report}", sep="\\t", header=FALSE, stringsAsFactors=FALSE)
    colnames(kraken_report) <- c("percent", "clade_reads", "direct_reads", "rank", "taxid", "name")
    
    # Read ASV abundance table (now ASVs are column names: ASV1, ASV2, etc.)
    asv_table <- read.csv("${asvtable_csv}", row.names=1, check.names=FALSE)
    
    # Clean up ASV IDs in kraken output
    kraken_output\$asv_id <- gsub("^>", "", kraken_output\$asv_id)
    
    # Transpose ASV table to get ASVs as rows, samples as columns
    asv_table_t <- t(asv_table) %>% as.data.frame() %>%
      tibble::rownames_to_column("asv_id")
    
    # Merge with taxonomy (now direct match on asv_id)
    asv_with_tax <- asv_table_t %>%
      left_join(kraken_output %>% select(asv_id, taxid, classified), by = "asv_id") %>%
      left_join(kraken_report %>% select(taxid, name, rank) %>% distinct(), by = "taxid") %>%
      select(asv_id, taxid, name, rank, classified, everything())
    
    # Save outputs
    write.table(asv_with_tax, "kraken2_taxonomy.txt", sep="\\t", quote=FALSE, row.names=FALSE)
    """
}
