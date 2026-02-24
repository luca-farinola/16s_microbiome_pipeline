/* =========================================
   DADA2 ASV Inference
   ========================================= */

nextflow.enable.dsl = 2

/* =========================================
   PROCESS: DADA2 - Learn Error Models and Infer ASVs
   ========================================= */

process DADA2_ERR {

      tag { sample_id }

      publishDir "${params.output_dir}/${params.DADA2ERR_dir}", mode: 'rellink'
   
      conda "bioconda::bioconductor-dada2=1.34.0 bioconda::bioconductor-biostrings conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"
   
      input:
      tuple val(sample_id), path(read1), path(read2)
   
      output:
         tuple val(sample_id), path("*.errF.rds"), emit: errF
         tuple val(sample_id), path("*.errR.rds"), emit: errR
         tuple val(sample_id), path("*.errF.pdf"), emit: pdfF
         tuple val(sample_id), path("*.errR.pdf"), emit: pdfR
         tuple val(sample_id), path("*.errF.svg"), emit: svgF
         tuple val(sample_id), path("*.errR.svg"), emit: svgR
         tuple val(sample_id), path("*.errF.log"), emit: logF
         tuple val(sample_id), path("*.errR.log"), emit: logR

         script:
         """
         #!/usr/bin/env Rscript
            suppressPackageStartupMessages(library(dada2))

            sample_id <- "${sample_id}"
            set.seed(1234)

            # Find filtered forward and reverse reads
            fnFs <- sort(list.files(".", pattern = paste0(sample_id, "_1.filt.fastq.gz"), full.names = TRUE))
            fnRs <- sort(list.files(".", pattern = paste0(sample_id, "_2.filt.fastq.gz"), full.names = TRUE))

            # Learn error models
            sink(file = paste0(sample_id, ".errF.log"))
            errF <- learnErrors(fnFs, multithread = ${task.cpus}, verbose = TRUE)
            saveRDS(errF, paste0(sample_id, ".errF.rds"))
            sink(file = NULL)

            sink(file = paste0(sample_id, ".errR.log"))
            errR <- learnErrors(fnRs, multithread = ${task.cpus}, verbose = TRUE)
            saveRDS(errR, paste0(sample_id, ".errR.rds"))
            sink(file = NULL)

            # Diagnostic plots
            pdf(paste0(sample_id, ".errF.pdf"))
            plotErrors(errF, nominalQ = TRUE)
            dev.off()
            pdf(paste0(sample_id, ".errR.pdf"))
            plotErrors(errR, nominalQ = TRUE)
            dev.off()
            svg(paste0(sample_id, ".errF.svg"))
            plotErrors(errF, nominalQ = TRUE)
            dev.off()
            svg(paste0(sample_id, ".errR.svg"))
            plotErrors(errR, nominalQ = TRUE)
            dev.off()

            # Convergence check
            sink(file = paste0(sample_id, ".errF.convergence.txt"))
            try(dada2:::checkConvergence(errF), silent = FALSE)
            sink(file = NULL)
            sink(file = paste0(sample_id, ".errR.convergence.txt"))
            try(dada2:::checkConvergence(errR), silent = FALSE)
            sink(file = NULL)

            # Record invocation args
            writeLines(paste0("learnErrorsF\tcpus=", ${task.cpus}), con = "learnErrorsF.args.txt")
            writeLines(paste0("learnErrorsR\tcpus=", ${task.cpus}), con = "learnErrorsR.args.txt")

            # Versions provenance
            writeLines(c(paste0('"', "${task.process}", '":'),
                               paste0("  R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
                               paste0("  dada2: ", packageVersion("dada2")) ), con = "versions.yml")

            """}


process DADA2_DENOISING {
   tag { sample_id }

   publishDir "${params.output_dir}/${params.DADA2denoise_dir}", mode: 'rellink'
   
   conda "bioconda::bioconductor-dada2=1.34.0 bioconda::bioconductor-biostrings conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"

   input:
      tuple val(sample_id), path(read1), path(read2), path(errF), path(errR)

   output:
      tuple val(sample_id), path("*.dadaF.rds"), emit: dadaF
      tuple val(sample_id), path("*.dadaR.rds"), emit: dadaR
      tuple val(sample_id), path("*.dada.log"), emit: log

   script:
   """
   #!/usr/bin/env Rscript
   suppressPackageStartupMessages(library(dada2))

   sample_id <- "${sample_id}"
   read1_file <- "${read1}"
   read2_file <- "${read2}"
   errF_file <- "${errF}"
   errR_file <- "${errR}"

   set.seed(1234)

   # Load error models
   errF <- readRDS(errF_file)
   errR <- readRDS(errR_file)

   # Denoise forward and reverse reads
   sink(file = paste0(sample_id, ".dada.log"))
   dadaF <- dada(read1_file, err = errF, multithread = ${task.cpus}, verbose = TRUE)
   dadaR <- dada(read2_file, err = errR, multithread = ${task.cpus}, verbose = TRUE)
   saveRDS(dadaF, paste0(sample_id, ".dadaF.rds"))
   saveRDS(dadaR, paste0(sample_id, ".dadaR.rds"))
   sink(file = NULL)
   """
}

process DADA2_MERGE {

   tag { sample_id }

   publishDir "${params.output_dir}/${params.DADA2merge_dir}", mode: 'rellink'

   conda "bioconda::bioconductor-dada2=1.34.0 bioconda::bioconductor-biostrings conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"

   input:
      tuple val(sample_id), path(read1), path(read2), path(dadaF), path(dadaR)

   output:
      tuple val(sample_id), path("*.merged.rds"), emit: merged

   script:
   """
   #!/usr/bin/env Rscript
   suppressPackageStartupMessages(library(dada2))

   sample_id <- "${sample_id}"
   read1_file <- "${read1}"
   read2_file <- "${read2}"
   dadaF_file <- "${dadaF}"
   dadaR_file <- "${dadaR}"

   # Load denoised objects
   dadaF <- readRDS(dadaF_file)
   dadaR <- readRDS(dadaR_file)

   # Load filtered reads
   filtFs <- read1_file
   filtRs <- read2_file

   # Merge pairs
   mergers <- mergePairs(dadaF, filtFs, dadaR, filtRs, verbose=TRUE)
   saveRDS(mergers, paste0(sample_id, ".merged.rds"))
   """
}


process DADA2_SEQTABLE {

   publishDir "${params.output_dir}/${params.DADA2asvtable_dir}", mode: 'rellink'

   conda "bioconda::bioconductor-dada2=1.34.0 bioconda::bioconductor-biostrings conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"

   input:
      path(merged_files)

   output:
      path("asvtable.csv"), emit: asvtable
      path("asvtable.rds"), emit: asvtable_rds
      path("*.fasta"), emit: asv_fasta

script:
"""
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(Biostrings))

merger_files <- list.files(pattern = "*.merged.rds")
sample_names <- gsub(".merged.rds", "", merger_files)
mergers <- lapply(merger_files, readRDS)
names(mergers) <- sample_names
seqtab <- makeSequenceTable(mergers)
rownames(seqtab) <- sample_names
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Create ASV IDs and rename columns
asv_seqs <- colnames(seqtab_nochim)
asv_ids <- paste0("ASV", seq_along(asv_seqs))
colnames(seqtab_nochim) <- asv_ids

# Write ASV table with ASV IDs as column names
write.csv(seqtab_nochim, file = "asvtable.csv", quote = FALSE, row.names = TRUE)
saveRDS(seqtab_nochim, "asvtable.rds")

# Write FASTA with sequences
dna_seqs <- DNAStringSet(asv_seqs)
names(dna_seqs) <- asv_ids
writeXStringSet(dna_seqs, "ASVs.fasta", format="fasta", width=200)
"""
}