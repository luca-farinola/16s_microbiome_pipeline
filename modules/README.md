# Nextflow Modules Documentation

## Overview

This directory contains modularized Nextflow processes for the 16S Microbiome Pipeline. Each module is self-contained and can be reused in other pipelines.

## Module Structure

### Directory Layout
```
modules/
├── quality_control.nf     # Quality assessment and filtering
├── assembly.nf            # Read merging and format conversion
├── taxonomy.nf            # Taxonomic classification
├── quantification.nf      # OTU abundance quantification
├── analysis.nf            # R-based analysis (phyloseq)
└── reporting.nf           # Quality reporting
```

## Modules Description

### 1. quality_control.nf
**Purpose**: Initial quality assessment and trimming of raw reads

**Processes**:
- `FASTQC` - Quality control metrics on raw reads
  - Input: Raw paired-end FASTQ files
  - Output: HTML reports + ZIP files
  
- `TRIMMING` - Quality filtering with fastp
  - Input: Raw paired-end FASTQ files
  - Output: Quality-trimmed FASTQ files + statistics

**Key Features**:
- Auto-detects adapters
- Quality threshold: Q20
- Generates HTML and JSON reports

---

### 2. assembly.nf
**Purpose**: Merge paired-end reads and prepare for taxonomy assignment

**Processes**:
- `MERGE_READS` - Assemble overlapping paired reads with vsearch
  - Input: Quality-trimmed paired FASTQ files
  - Output: Merged FASTQ files + merge statistics
  
- `FASTQ_TO_FASTA` - Convert FASTQ to FASTA format
  - Input: Merged FASTQ files
  - Output: FASTA files

**Key Features**:
- Uses vsearch for efficient merging
- Relabels reads with sample prefix
- Prepares for pathoscope input

---

### 3. taxonomy.nf
**Purpose**: Assign taxonomy using Pathoscope against SILVA database

**Processes**:
- `PATHOSCOPE_CLASSIFY` - Taxonomic classification
  - Input: FASTA sequences
  - Output: Pathoscope reports + SAM alignment files

**Key Features**:
- Uses SILVA 16S reference database
- Generates statistical confidence scores
- Produces bowtie2 alignment files

---

### 4. quantification.nf
**Purpose**: Quantify OTU abundances and combine samples

**Processes**:
- `GENERATE_OTU_TABLE` - Parse pathoscope output to OTU counts
  - Input: Pathoscope reports
  - Output: Per-sample OTU count tables (CSV format)
  
- `COMPILE_OTU_TABLE` - Combine all samples into single OTU matrix
  - Input: All per-sample OTU tables
  - Output: Combined OTU table (TSV format) **⭐ RAW OTU COUNTS**

**Key Features**:
- Preserves raw (non-normalized) counts
- Handles missing values across samples
- Ready for downstream analysis

---

### 5. analysis.nf
**Purpose**: Create R-ready phyloseq object from OTU data

**Processes**:
- `CREATE_PHYLOSEQ` - Generate phyloseq object in R
  - Input: Combined OTU table
  - Output: phyloseq object (RDS format) + summary **⭐ FOR R ANALYSIS**

**Key Features**:
- Installs R packages automatically
- Includes OTU table, taxonomy, and sample metadata
- Generates object summary statistics

---

### 6. reporting.nf
**Purpose**: Generate comprehensive quality reports

**Processes**:
- `MULTIQC` - Aggregate quality metrics
  - Input: FastQC and fastp reports
  - Output: Interactive HTML report + data folder **⭐ QUALITY REPORT**

**Key Features**:
- Aggregates reports from all samples
- Interactive HTML visualization
- Downloadable data tables

---

## Using Modules in main_pathscope_16s.nf

### Import Syntax
```nextflow
include { FASTQC; TRIMMING } from './modules/quality_control.nf'
include { MERGE_READS; FASTQ_TO_FASTA } from './modules/assembly.nf'
include { PATHOSCOPE_CLASSIFY } from './modules/taxonomy.nf'
include { GENERATE_OTU_TABLE; COMPILE_OTU_TABLE } from './modules/quantification.nf'
include { CREATE_PHYLOSEQ } from './modules/analysis.nf'
include { MULTIQC } from './modules/reporting.nf'
```

### Invoking Processes
```nextflow
workflow {
    // Use named emit outputs for clarity
    FASTQC(raw_reads_ch)
    TRIMMING(raw_reads_ch)
    MERGE_READS(TRIMMING.out.fastq)
    FASTQ_TO_FASTA(MERGE_READS.out.fastq)
    PATHOSCOPE_CLASSIFY(FASTQ_TO_FASTA.out.fasta)
    GENERATE_OTU_TABLE(PATHOSCOPE_CLASSIFY.out.report)
    COMPILE_OTU_TABLE(GENERATE_OTU_TABLE.out.counts.collect())
    CREATE_PHYLOSEQ(COMPILE_OTU_TABLE.out.combined)
    MULTIQC(FASTQC.out.zip.collect(), TRIMMING.out.json.collect())
}
```

---

## Module Output Naming

Each module uses **named emits** for clear, documented outputs:

| Module | Process | Output Name | Type |
|--------|---------|------------|------|
| quality_control | FASTQC | html, zip | HTML reports, ZIP files |
| quality_control | TRIMMING | fastq, json, html | Trimmed reads, statistics |
| assembly | MERGE_READS | fastq, stats | Merged reads, merge stats |
| assembly | FASTQ_TO_FASTA | fasta | FASTA sequences |
| taxonomy | PATHOSCOPE_CLASSIFY | report, sam | Taxonomy report, alignments |
| quantification | GENERATE_OTU_TABLE | table, counts | OTU matrix, CSV counts |
| quantification | COMPILE_OTU_TABLE | combined | Combined OTU table (TSV) |
| analysis | CREATE_PHYLOSEQ | object, summary | RDS object, summary text |
| reporting | MULTIQC | html, data | HTML report, QC data |

---

## Extending the Pipeline

### Adding a New Process to Existing Module

1. Edit the appropriate module file (e.g., `quantification.nf`)
2. Add a new `process` block with proper inputs/outputs
3. Use named `emit` directives for clear output naming
4. Import in main file with `include { NEW_PROCESS }`

### Creating a New Module

1. Create new file: `modules/new_feature.nf`
2. Add DSL2 header: `nextflow.enable.dsl = 2`
3. Define processes with named emits
4. Import in main file

### Example Adding a New Process:

```nextflow
// In modules/quantification.nf

process NORMALIZE_OTU_TABLE {
    publishDir "${params.output_dir}/08.Normalized_OTU", mode: 'copy'
    
    container 'docker://biocontainers/r-base:4.3.1_cv2'
    
    input:
    path otu_table
    
    output:
    path "normalized_otu_table.tsv", emit: normalized
    
    script:
    """
    R --vanilla <<'EOF'
    # Normalization code here
    EOF
    """
}
```

Then in `main_pathscope_16s.nf`:

```nextflow
include { GENERATE_OTU_TABLE; COMPILE_OTU_TABLE; NORMALIZE_OTU_TABLE } from './modules/quantification.nf'

workflow {
    // ... existing code ...
    COMPILE_OTU_TABLE(GENERATE_OTU_TABLE.out.counts.collect())
    NORMALIZE_OTU_TABLE(COMPILE_OTU_TABLE.out.combined)
}
```

---

## Module Parameters

All modules use global `params` defined in main file:

```nextflow
params.output_dir = "../02.Analysis"
params.silva_db = "/path/to/silva/database"
params.pathscope_ref = "/path/to/pathoscope/reference"
```

These are accessible in all modules without explicit passing.

---

## Container Images Used

Each module specifies its own container image:

| Module | Container | Source |
|--------|-----------|--------|
| quality_control | biocontainers/fastqc:v0.11.9_cv8 | Docker |
| quality_control | biocontainers/fastp:v0.23.2_cv2 | Docker |
| assembly | biocontainers/vsearch:v2.24.0_cv1 | Docker |
| assembly | biocontainers/seqtk:v1.3_cv3 | Docker |
| taxonomy | pathoscope:latest | Custom |
| quantification | biocontainers/qiime2:2023.9_cv1 | Docker |
| quantification | biocontainers/r-base:4.3.1_cv2 | Docker |
| analysis | biocontainers/r-base:4.3.1_cv2 | Docker |
| reporting | biocontainers/multiqc:v1.14_cv7 | Docker |

---

## Module Dependencies

**Data Flow** (→ indicates dependency):

```
raw_reads
    ↓
[FASTQC + TRIMMING] → MERGE_READS → FASTQ_TO_FASTA → PATHOSCOPE_CLASSIFY 
                                                           ↓
                                                  GENERATE_OTU_TABLE 
                                                           ↓
                                                  COMPILE_OTU_TABLE 
                                                           ↓
                                                  CREATE_PHYLOSEQ
                                                  
FASTQC output ──┐
                ├─→ MULTIQC
TRIMMING output ┘
```

---

## Best Practices

1. **Use Named Emits**: Makes output referencing clear and maintainable
2. **Consistent Naming**: Use consistent process naming conventions
3. **Documentation**: Add comments explaining process purpose
4. **Error Handling**: Containers should handle errors gracefully
5. **Container Versioning**: Specify exact version tags (not `latest`)
6. **publishDir**: Always specify output directory for public results

---

## Troubleshooting

### Module Not Found
- Check import statement matches module file name
- Verify relative path is correct (use `./modules/`)
- Ensure DSL2 is enabled: `nextflow.enable.dsl = 2`

### Output Not Found
- Verify `emit` name in process matches usage in workflow
- Check `publishDir` is set correctly
- Review `.nextflow.log` for errors

### Container Issues
- Verify internet connectivity for pulling images
- Check Singularity cache permissions
- Pre-pull images: `singularity pull docker://image_name`

---

## Changes from Monolith to Modular

| Aspect | Before | After |
|--------|--------|-------|
| File Size | ~420 lines | ~115 lines (main) |
| Process Management | All inline | Organized in 6 modules |
| Reusability | Low | High |
| Maintainability | Difficult | Easy |
| Scalability | Limited | Extensible |
| Output Organization | Unclear | Named emits |

---

## Version

**Created**: February 2025  
**Module Structure**: DSL2 (Nextflow ≥21.10.0)  
**Status**: Production Ready ✅

---

For questions or improvements, see main pipeline documentation files.
