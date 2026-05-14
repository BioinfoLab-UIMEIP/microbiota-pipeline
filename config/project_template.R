# =============================================================================
# Config:      project_template.R
# Pipeline:    Microbiota 16S Amplicon Analysis Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, UIMEIP, CMN SXXI, IMSS, Mexico City
# Description: Template configuration — copy to your project directory and
#              fill in the paths and parameters for your study.
#              Save as: /path/to/my_project/my_study.R
# =============================================================================

# ── REQUIRED ──────────────────────────────────────────────────────────────────
# Directory containing raw paired-end FASTQ files (R1/R2 pairs).
RAW_DIR      <- "/path/to/raw_fastq/"

# Base directory for pipeline outputs. A timestamped subdirectory
# (results_YYYYMMDD_HHMMSS) will be created here on each run.
RESULTS_BASE <- "/path/to/my_project/"

# Sequencing kit (must match an entry in scripts/primers.tsv).
KIT    <- "QIAseq"   # e.g. "QIAseq", "EMP", "NEBNext"

# 16S hypervariable region (must match an entry in scripts/primers.tsv).
REGION <- "V4V5"     # e.g. "V4V5", "V3V4", "V1V2"

# ── SAMPLE SELECTION ──────────────────────────────────────────────────────────
# Filter the sample manifest to a specific tissue/sample type before any
# processing. Set to NULL to include all samples.
# The manifest is joined with METADATA_FILE on SAMPLE_ID_COL so that only
# samples matching the filter value(s) are retained from step 1 onward.
SAMPLE_TYPE_COL    <- "sample_type"        # column name in metadata CSV
SAMPLE_TYPE_FILTER <- "Biopsy"             # "Biopsy" | "Saliva" | c("Biopsy","Saliva") | NULL

# ── COMPARISON GROUPS ─────────────────────────────────────────────────────────
# Path to sample metadata CSV. Required columns: SAMPLE_ID_COL, GROUP_COL,
# and SAMPLE_TYPE_COL (if SAMPLE_TYPE_FILTER is used).
METADATA_FILE <- "/path/to/my_project/metadata.csv"

# Column in metadata that matches the FASTQ-derived sample identifier.
SAMPLE_ID_COL <- "Sample_ID"

# Column and labels that define the case/control contrast.
GROUP_COL  <- "Group"
CASE_LABEL <- "Case"      # value in GROUP_COL for the case group
CTRL_LABEL <- "Control"   # value in GROUP_COL for the control group

# ── TAXONOMY FILTER ───────────────────────────────────────────────────────────
# Kingdoms to retain after assignment. Set KEEP_UNASSIGNED_KINGDOM = TRUE to
# keep ASVs whose kingdom could not be assigned (avoids silent data loss).
KEEP_KINGDOMS           <- c("Bacteria", "Archaea")
KEEP_UNASSIGNED_KINGDOM <- TRUE

# ── REFERENCE DATABASES ───────────────────────────────────────────────────────
# SILVA 138.2 — mandatory (download from https://zenodo.org/record/8392695)
SILVA_TRAINSET <- "/path/to/db/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
SILVA_SPECIES  <- "/path/to/db/silva_v138.2_assignSpecies.fa.gz"

# HOMD — optional; recommended for oral/gastric samples
# (download from https://www.homd.org/ftp/)
# Set to "" or comment out to skip the HOMD tier.
HOMD_TRAINSET  <- "/path/to/db/HOMD_trainset_toGenus.fasta"
HOMD_SPECIES   <- "/path/to/db/HOMD_assignSpecies_dada2.fa"

# NCBI RefSeq RNA BLAST database — for residual unresolved features
# (download with: update_blastdb.pl --decompress 16S_ribosomal_RNA)
REFSEQ_BLASTDB <- "/path/to/db/blastdb/refseq_rna"

# ── OPTIONAL OVERRIDES ────────────────────────────────────────────────────────
# Leave commented out to use the conda-environment defaults.

# External binaries (defaults assume they are on PATH via conda env)
# CUTADAPT_BIN    <- "cutadapt"
# VSEARCH_BIN     <- "vsearch"
# MAFFT_SMART_BIN <- "mafft"   # or path to mafft_smart.sh wrapper
# SNP_DISTS_BIN   <- "snp-dists"
# BLASTN_BIN      <- "blastn"
# BSDTAR_BIN      <- "bsdtar"

# DADA2 truncation lengths [R1, R2] in bp.
# If not set, values are read from qc/trunc_suggestion.txt produced by 01_qc.R.
# TRUNC_LEN <- c(260L, 220L)

# ASV collapse parameters
ENABLE_ASV_COLLAPSE <- TRUE    # run the post-denoising collapse step
USE_COLLAPSED_ASVS  <- TRUE    # use collapsed ASVs in downstream steps
VSEARCH_IDENTITY    <- 0.999   # clustering identity threshold
HC_SNP_THRESHOLD    <- 1L      # SNP distance cut-off for hierarchical clustering
HC_METHOD           <- "complete"
MAFFT_MODE          <- "smart"
MAFFT_THREADS       <- 4L      # adjust to available CPU cores

# Taxonomy assignment
TAX_CHUNK_SIZE <- 250L   # features per BLAST chunk
TAX_N_JOBS     <- 4L     # parallel BLAST chunks (set to available cores)
KEEP_TAX_INTERMEDIATES <- FALSE

# Feature filtering (applied in 06_filter.R)
FEATURE_LEVEL              <- "species"
MIN_PREVALENCE             <- 2L
MIN_TOTAL_ABUNDANCE_UNIQUE <- 10L

# Beta diversity
PERMANOVA_NPERM <- 999L

# Differential abundance
DA_MIN_PREV    <- 2L
DA_PSEUDOCOUNT <- 1e-5

# Machine learning
RF_NTREE <- 500L
RF_NFOLD <- 5L
TOP_N_RF <- 20L
