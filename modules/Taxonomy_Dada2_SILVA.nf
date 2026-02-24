/* =========================================
   Taxonomy Assignment with DADA2 + SILVA
   ========================================= */

nextflow.enable.dsl = 2

/* =========================================
   PROCESS: DADA2_ASSIGNTAXONOMY
   Assign taxonomy using DADA2's naive Bayes classifier with SILVA reference
   ========================================= */

process DADA2_ASSIGNTAXONOMY {

    tag "DADA2 Taxonomy Assignment"
    
    publishDir "${params.output_dir}/${params.taxonomy_dir_dada2}/dada2_taxonomy", mode: 'copy'
    
   conda "bioconda::bioconductor-dada2=1.34.0 bioconda::bioconductor-biostrings conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"
    input:
    path(asv_fasta)
    path(silva_train_set)
    path(silva_species_set)

    output:
    path("tax_assignment.rds"), emit: tax_rds
    path("tax_assignment.csv"), emit: tax_csv

    script:
    """
    Rscript - <<'EOF'
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(Biostrings))

    # Read ASV sequences from FASTA
    asv_seqs_dna <- readDNAStringSet("${asv_fasta}")
    
    # Assign taxonomy using SILVA
    taxa <- assignTaxonomy(asv_seqs_dna, "${silva_train_set}", multithread=${task.cpus}, verbose=TRUE, minBoot=50)
    
    # Add species assignment if available
    taxa_species <- addSpecies(taxa, "${silva_species_set}", verbose=TRUE, allowMultiple=FALSE)
    
    # Save as RDS
    saveRDS(taxa_species, "tax_assignment.rds")
    
    # Convert to data frame and save as CSV
    taxa_df <- as.data.frame(taxa_species)
    taxa_df\$ASV <- names(asv_seqs_dna)
    taxa_df <- taxa_df[, c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
    write.csv(taxa_df, "tax_assignment.csv", row.names=FALSE, quote=FALSE)
    EOF
    """
}

/* =========================================
   PROCESS: DADA2_TAX_TABLE
   Create abundance table aggregated by taxonomy
   ========================================= */

process DADA2_TAX_TABLE {

    tag "Create DADA2 Taxonomy Abundance Table"
    
    publishDir "${params.output_dir}/${params.Final_dir_dada2}/dada2_taxonomy", mode: 'copy'
    
    container "community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186"

    input:
    path(asvtable_csv)
    path(tax_assignment_csv)

    output:
    path("dada2_abundance_by_rank.txt"), emit: by_rank
    path("dada2_abundance_aggregated.txt"), emit: aggregated
    path("asv_with_dada2_taxonomy.csv"), emit: asv_tax
    path("dada2_tax_summary.txt"), emit: summary

    script:
    """
    Rscript - <<'EOF'

    library(dplyr)
    library(tidyr)
    library(tibble)

    # Read ASV abundance table
    asv_table <- read.csv("${asvtable_csv}", row.names=1, stringsAsFactors=FALSE)
    
    # Read DADA2 taxonomy assignments
    tax_assign <- read.csv("${tax_assignment_csv}", stringsAsFactors=FALSE)
    
    # Transpose ASV table (ASVs as rows, samples as columns)
    asv_table_t <- t(asv_table) %>% as.data.frame() %>%
      tibble::rownames_to_column("sample")


    # Create ASV-abundance long format
    asv_long <- asv_table_t %>%
      pivot_longer(cols = -sample, names_to = "asv_id", values_to = "abundance") %>%
      filter(abundance > 0)

    asv_long$ASV <- asv_long$sample
    asv_long$sample <- NULL
    
    # Merge with taxonomy
    asv_with_tax <- asv_long %>%
      left_join(tax_assign, by = "ASV")
    
    # Save full table with taxonomy
    write.csv(asv_with_tax, "asv_with_dada2_taxonomy.csv", row.names=FALSE, quote=FALSE)
    
    """
}
