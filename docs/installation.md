---
layout: page
title: "Installation"
nav_order: 2
---

# Installation

## Quick-start

1. Install Miniconda3 (macOS or Linux).
2. Create and activate the conda environment:
   ```bash
   conda create -n microbiota16s -c conda-forge -c bioconda \
     r-base=4.2 "cutadapt>=4.0" "vsearch>=2.23" "mafft>=7.5" "snp-dists>=0.8" "blast>=2.14"
   conda activate microbiota16s
   ```
3. Install R packages (see [step 3](#3-install-r-packages) below).
4. Download SILVA 138.2 and RefSeq RNA databases (see [step 4](#4-download-reference-databases) below).
5. Clone the pipeline and copy the config template to your project directory:
   ```bash
   git clone https://github.com/BioinfoLab-UIMEIP/microbiota-pipeline.git
   cp microbiota-pipeline/config/project_template.R /path/to/my_project/my_study.R
   ```
6. Edit `my_study.R` with paths and labels for your study. Minimal required variables:
   ```r
   RAW_DIR       <- "/path/to/raw_fastq"
   RESULTS_BASE  <- "/path/to/results"
   KIT           <- "QIAseq"
   REGION        <- "V4V5"
   METADATA_FILE <- "/path/to/metadata.csv"
   SAMPLE_ID_COL <- "Sample_ID"
   GROUP_COL     <- "Group"
   CASE_LABEL    <- "Case"
   CTRL_LABEL    <- "Control"
   SILVA_TRAINSET <- "/path/to/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
   SILVA_SPECIES  <- "/path/to/silva_v138.2_assignSpecies.fa.gz"
   REFSEQ_BLASTDB <- "/path/to/16S_ribosomal_RNA"
   ```
7. Run the full pipeline:
   ```bash
   bash microbiota-pipeline/run_pipeline.sh --config /path/to/my_project/my_study.R
   ```

---

## System requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| Operating system | macOS 12+ or Linux | macOS 14+ (Darwin 24.x) or Ubuntu 22.04 |
| RAM | 8 GB | 16 GB or more |
| Disk space (pipeline + databases) | 30 GB | 50 GB |
| CPU cores | 4 | 8 or more |
| R version | 4.2 | 4.2 |
| conda | Any recent version | Mambaforge or Miniconda |

The pipeline has been tested on macOS (Darwin 24.x) with R 4.2. Linux is expected to work without modification.

---

## 1. Install conda

If conda is not already present, install Mambaforge or Miniconda:

```bash
# macOS (Apple Silicon or Intel)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh

# Verify
conda --version
```

---

## 2. Create the conda environment

```bash
conda create -n microbiota16s \
  -c conda-forge -c bioconda -c defaults \
  r-base=4.2 \
  "cutadapt>=4.0" \
  "vsearch>=2.23" \
  "mafft>=7.5" \
  "snp-dists>=0.8" \
  "blast>=2.14"

conda activate microbiota16s
```

Confirm that all external binaries are on `PATH`:

```bash
cutadapt --version
vsearch --version
mafft --version
snp-dists --version
blastn -version
```

---

## 3. Install R packages

Open an R session **inside the activated conda environment** (`conda activate microbiota16s` first):

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("dada2", "ShortRead", "Biostrings"))

# CRAN packages
install.packages(c(
  "vegan",
  "ape",
  "randomForest",
  "MASS",
  "ggplot2",
  "ggrepel",
  "gridExtra",
  "dplyr",
  "tidyr",
  "rmarkdown",
  "flexdashboard",
  "plotly",
  "DT",
  "randomcoloR",
  "reshape2",
  "tibble"
))
```

> **CRITICAL:** Install all R packages inside the conda environment that will run the pipeline. Installing them in a separate R library causes import errors at runtime.

---

## 4. Download reference databases

### SILVA 138.2 (mandatory)

Download both files from Zenodo:

```
https://zenodo.org/record/8392695
```

You need:
- `silva_nr99_v138.2_train_set.fa.gz` — training set for `assignTaxonomy()`
- `silva_species_assignment_v138.2.fa.gz` — species FASTA for `addSpecies()`

Set these paths in your project config:

```r
SILVA_TRAINSET <- "/path/to/databases/silva_nr99_v138.2_train_set.fa.gz"
SILVA_SPECIES  <- "/path/to/databases/silva_species_assignment_v138.2.fa.gz"
```

### Specialized database (optional — second taxonomy tier)

A site-specific reference database can improve species-level resolution for taxa not well covered by SILVA. For oral and gastric studies, HOMD is a common choice; for other environments, substitute an appropriate DADA2-compatible trainset.

```r
HOMD_TRAINSET <- "/path/to/databases/specialized_trainset.fa"   # e.g., HOMD_16S_rRNA_RefSeq_DADA2.fa
HOMD_SPECIES  <- "/path/to/databases/specialized_species.fa"    # optional species-level FASTA
```

Leave `HOMD_TRAINSET <- ""` (or omit the variable) to skip this tier and go directly from SILVA to RefSeq BLAST.

### NCBI RefSeq RNA BLAST database

The BLAST step handles residual unresolved features. The database download takes approximately 2–10 GB.

```bash
# Inside the conda environment:
cd /path/to/databases/
update_blastdb.pl --decompress 16S_ribosomal_RNA
```

If `update_blastdb.pl` is not available, download manually:

```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
tar -xzf 16S_ribosomal_RNA.tar.gz
```

Set in config:

```r
REFSEQ_BLASTDB <- "/path/to/databases/16S_ribosomal_RNA"
```

---

## 5. Clone the pipeline

```bash
git clone https://github.com/BioinfoLab-UIMEIP/microbiota-pipeline.git
cd microbiota-pipeline
```

Confirm the `scripts/primers.tsv` file is present — it is required for step 2 (primer trimming).

---

## 6. Verification

Run the following inside the activated environment to confirm all components are accessible:

```bash
# External binaries
cutadapt --version
vsearch --version
mafft --version
snp-dists --version
blastn -version

# R libraries
Rscript -e "
  for (pkg in c('dada2','ShortRead','Biostrings','vegan','randomForest',
                'MASS','ggplot2','rmarkdown','flexdashboard','plotly','DT')) {
    loaded <- requireNamespace(pkg, quietly=TRUE)
    cat(sprintf('%-20s %s\n', pkg, if (loaded) 'OK' else 'MISSING'))
  }
"
```

All packages should report `OK` before proceeding to configuration.

---

## Environment file

A pre-configured environment specification is available at the repository root as `environment.yml`:

```bash
conda env create -f environment.yml
conda activate microbiota16s
```

After activating, install the R packages as described in step 3 above.

---

Next: [Configuration reference](configuration.md)
