# 16S rRNA Microbiome Analysis Pipeline

A complete **Nextflow pipeline** for 16S rRNA gene amplicon sequencing analysis. Processes raw paired-end reads through quality control, primer trimming, ASV inference (DADA2), and dual taxonomic classification (Kraken2 + SILVA and DADA2 + SILVA).

## Pipeline Overview

```
Raw Reads (R1 & R2)
      ↓
  [FastQC] → Quality Assessment
      ↓
  [Cutadapt] → Primer Removal
      ↓
  [DADA2 FilterAndTrim] → Quality-based Filtering
      ↓
  [DADA2 Learn Error Models] → Error Rate Profiling
      ↓
  [DADA2 Denoising] → ASV Inference (R1 & R2)
      ↓
  [DADA2 Merge Pairs] → Combine R1 + R2 Reads
      ↓
  [DADA2 Sequence Table] → ASV Abundance Table (Chimera Removed)
      ↓
  ┌─────────────────────────────────────┐
  │ Parallel Taxonomy Classification    │
  ├─────────────────────────────────────┤
  │ [Kraken2 + SILVA 138]  (Standard)   │
  │ [DADA2 + SILVA 138]    (Naive Bayes)│
  └─────────────────────────────────────┘
      ↓
  [MultiQC] → Aggregate QC Report
      ↓
  [Output: Abundance Tables + Taxonomy]
```

## Features

✅ **Dual Taxonomy Methods**
  - **Kraken2 + SILVA 138**: Fast k-mer based classification
  - **DADA2 + SILVA 138**: Sensitive naive Bayes classifier with species assignment

✅ **Complete Quality Pipeline**
  - FastQC on raw reads
  - Cutadapt primer trimming
  - DADA2 quality filtering
  - Error model learning (per-sample)
  - Chimera detection & removal

✅ **Robust ASV Inference**
  - Per-sample error model learning
  - Forward and reverse read denoising
  - Paired-end read merging
  - Duplicate sequence merging
  - Bimera (chimera) removal (consensus method)

✅ **Comprehensive Output**
  - ASV abundance tables (CSV format)
  - Representative sequences (FASTA)
  - Taxonomic assignments (two methods)
  - Error profile visualizations (PDF/SVG)
  - MultiQC aggregated report

✅ **HPC Ready**
  - SLURM integration (tested on PALMA cluster)
  - Singularity containerization
  - Resource-aware task scheduling

## Requirements

### Software
- **Nextflow** ≥ 22.10.0
- **Singularity** (HPC) or **Docker** (local)

## Setup

### 3. Prepare Sample Sheet
Create `samples.csv`:
```csv
sample_id,read1_path,read2_path,forward_primer,reverse_primer
sample001,01.RawData/sample001/sample001_R1.fastq.gz,01.RawData/sample001/sample001_R2.fastq.gz,AGAGTTTGATCCTGGCTCAG,ATTACCGCGGCTGCTGG
sample002,01.RawData/sample002/sample002_R1.fastq.gz,01.RawData/sample002/sample002_R2.fastq.gz,AGAGTTTGATCCTGGCTCAG,ATTACCGCGGCTGCTGG
sample003,01.RawData/sample003/sample003_R1.fastq.gz,01.RawData/sample003/sample003_R2.fastq.gz,AGAGTTTGATCCTGGCTCAG,ATTACCGCGGCTGCTGG
```

## Running the Pipeline

```bash

nextflow run main_16s.nf \
    --output_dir output \
    -profile palma \
    --sample_csv  samples.csv\ 
```

## Output Structure

```
02.Analysis/
├── 01.FastQC/
│   ├── sample001_R1_fastqc.html
│   ├── sample001_R2_fastqc.html
│   └── sample001_*_fastqc.zip
│
├── 02.CUTADAPTTrimm/
│   ├── sample001_R1_trimmed.fastq.gz
│   ├── sample001_R2_trimmed.fastq.gz
│   └── sample001.trim.log
│
├── 03.DADA2Trimm/
│   ├── sample001_1.filt.fastq.gz
│   ├── sample001_2.filt.fastq.gz
│   └── sample001.filter_stats.tsv
│
├── 04.DADA2ErrorModel/
│   ├── sample001.errF.rds           ← Error model (R1)
│   ├── sample001.errR.rds           ← Error model (R2)
│   ├── sample001.errF.pdf           ← Error profile plot
│   ├── sample001.errR.pdf
│   ├── sample001.errF.svg
│   └── sample001.errR.svg
│
├── 05.DADA2Denoise/
│   ├── sample001.dadaF.rds          ← Denoised reads (R1)
│   ├── sample001.dadaR.rds          ← Denoised reads (R2)
│   └── sample001.dada.log
│
├── 06.DADA2Merge/
│   └── sample001.merged.rds         ← Merged paired-end reads
│
├── 07.ASVTable/
│   ├── asvtable.csv                 ← **ASV abundance (rows=samples, cols=ASVs)**
│   ├── asvtable.rds                 ← R format
│   └── ASVs.fasta                   ← Representative sequences
│
├── 08.MultiQC_Report/
│   ├── multiqc_report.html          ← **Aggregated QC report**
│   └── multiqc_data/
│
├── 09.Taxonomy_kraken/
│   ├── kraken2_output.txt           ← Per-ASV Kraken2 classifications
│   └── kraken2_report.txt           ← Kraken2 summary report
│
├── 09.Taxonomy_dada2/
│   ├── dada2_taxonomy/
│   │   ├── tax_assignment.rds       ← R format taxonomy object
│   │   └── tax_assignment.csv       ← Taxonomy table (Kingdom→Species)
│
├── 10.Aboundance_table_kraken/
│   └── kraken2_taxonomy.txt         ← **ASVs + Kraken2 taxonomy + counts**
│
├── 10.Aboundance_table_dada2/
│   ├── asv_with_dada2_taxonomy.csv  ← **ASVs + DADA2 taxonomy + counts**
│   ├── dada2_abundance_by_rank.txt  ← Aggregated by taxonomic rank
│   └── dada2_tax_summary.txt        ← Summary statistics
│
├── pipeline_report.html             ← Nextflow execution report
├── timeline.html                    ← Task execution timeline
└── trace.tsv                        ← Detailed task metrics
```

## Key Output Files

### 1. ASV Abundance Table (`07.ASVTable/asvtable.csv`)
Rows = samples, Columns = ASV IDs (ASV1, ASV2, ..., ASVn)
```
,ASV1,ASV2,ASV3,...
sample001,1234,567,89,...
sample002,2345,678,12,...
```

### 2. Kraken2 Taxonomy (`10.Aboundance_table_kraken/kraken2_taxonomy.txt`)
Combined ASV counts + Kraken2 classifications:
```
asv_id  taxid  name                              rank       classified  sample001  sample002
ASV1    573   Bacteria                          domain     C           1234       2345
ASV2    562   Bacteroides genus                 genus      C           567        678
```

### 3. DADA2 Taxonomy (`10.Aboundance_table_dada2/asv_with_dada2_taxonomy.csv`)
Full lineage assignments with bootstrap confidence:
```
ASV,Kingdom,Phylum,Class,Order,Family,Genus,Species,sample001,sample002
ASV1,Bacteria,Firmicutes,Bacilli,Bacillales,Bacillaceae,Bacillus,subtilis,1234,2345
ASV2,Bacteria,Bacteroidota,Bacteroidia,Bacteroidales,Bacteroidaceae,Bacteroides,vulgatus,567,678
```

### 4. Representative Sequences (`07.ASVTable/ASVs.fasta`)
FASTA with unique ASV sequences:
```
>ASV1
AGTTTGATCCTGGCTCAGTACACGTGCACAAGGCGATTTGAGACTGATGTTGTGTGT...
>ASV2
AGTTTGATCCTGGCTCAGTACACGTGCACAAGGCGATCCGAAACTGATGCTGTGTGT...
```

### 5. Error Profiles (`04.DADA2ErrorModel/*.pdf`)
Visual quality control: Observed vs. Expected error rates per quality score

### 6. MultiQC Report (`08.MultiQC_Report/multiqc_report.html`)
Interactive HTML dashboard:
- FastQC per-base quality scores
- Cutadapt trimming statistics
- DADA2 ASV inference summary
- Per-sample QC metrics



## things to fix 

these lines expose internal server paths:

 - kraken2_db = "/cloud/wwu1/e_expat/expatstorage/reference_annotations/Microbiome/Kraken2_SILVA/16S_SILVA138_k2db"
 - silva_train_set = "/cloud/wwu1/e_expat/expatstorage/reference_annotations/Microbiome/DADA2/silva_nr99_v138.2_toGenus_trainset.fa.gz"
 - silva_species_set = "/cloud/wwu1/e_expat/expatstorage/reference_annotations/Microbiome/DADA2/silva_v138.2_assignSpecies.fa.gz"

create docker immage that people can pull instead of sing seqera 
make a flag so that you could choose wehter to use kraken or DADA2 for taxonomy 
