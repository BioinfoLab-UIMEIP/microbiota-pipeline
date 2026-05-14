---
layout: page
title: "Configuration"
nav_order: 3
---

# Configuration

The pipeline is driven by a single R configuration file (a plain `.R` script that is `source()`d at the start of each step). Every parameter has a documented default inside the step scripts, so only the parameters relevant to your study need to be set.

---

## Creating a project configuration file

```bash
# Copy the template to your project directory
cp /path/to/pipeline/config/project_template.R /path/to/my_project/my_study.R

# Edit the copy — do NOT edit the original template
```

> **CRITICAL:** Keep the configuration file in your project directory alongside your data, not inside the pipeline source tree. Pass it explicitly every time a script or `run_pipeline.sh` is called.

---

## Core parameters

### File paths

| Parameter | Type | Description |
|-----------|------|-------------|
| `RAW_DIR` | character | Directory containing raw paired-end FASTQ files (Illumina naming: `*_L001_R1_001.fastq.gz`) |
| `RESULTS_BASE` | character | Base directory for output. The pipeline creates `results_YYYYMMDD_HHMMSS/` inside this directory on the first run |
| `METADATA_FILE` | character | Full path to the sample metadata CSV or TSV file |
| `MANIFEST_FILE` | character | Optional. Path to a pre-built manifest TSV with columns `sample_id`, `read1`, `read2`. If empty, the manifest is built automatically from `RAW_DIR` |

### Metadata columns

| Parameter | Default | Description |
|-----------|---------|-------------|
| `SAMPLE_ID_COL` | `"Sample_ID"` | Column in the metadata CSV that contains sample identifiers. Must match the prefix of the FASTQ filenames after removing `_L001_R1_001.fastq.gz` and converting non-alphanumeric characters to underscores |
| `GROUP_COL` | `"Group"` | Column containing case/control labels |
| `SAMPLE_TYPE_COL` | `"sample_type"` | Column containing tissue-type or sample-type labels |
| `CASE_LABEL` | — | Exact string value for the case group in `GROUP_COL` (e.g., `"Case"`) |
| `CTRL_LABEL` | — | Exact string value for the control group in `GROUP_COL` (e.g., `"Control"`) |
| `SAMPLE_TYPE_FILTER` | `NULL` | Character vector of tissue types to retain (e.g., `c("Biopsy")`). Set to `NULL` to include all sample types |

> **CRITICAL:** `SAMPLE_ID_COL` values must match FASTQ filename prefixes exactly (after the alphanumeric normalization applied in step 1). Mismatches silently exclude samples from downstream joins.

### Sequencing library parameters

| Parameter | Description |
|-----------|-------------|
| `KIT` | Library preparation kit. Supported values: `"canonical"`, `"EMP"`, `"Illumina"`, `"QIAseq"`, `"Zymo"`. Must exist in `scripts/primers.tsv` |
| `REGION` | Hypervariable region. Supported values depend on kit. Common values: `"V4"`, `"V4V5"`, `"V3V4"`, `"V1V2"`, `"V2V3"`, `"V5V7"`, `"V7V9"` |
| `PRIMER_FWD` | Optional. Override the forward primer from `primers.tsv` with this sequence |
| `PRIMER_REV` | Optional. Override the reverse primer from `primers.tsv` with this sequence |

### DADA2 parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TRUNC_LEN` | `c(0, 0)` | Forward and reverse truncation lengths in cycles. If both are 0, the pipeline reads the suggestion from `qc/trunc_suggestion.txt` automatically |
| `MAX_EE` | `c(2, 2)` | Maximum expected errors for forward and reverse reads |
| `MAX_N` | `0` | Maximum ambiguous base calls allowed |
| `RM_PHIX` | `TRUE` | Remove PhiX reads |
| `MIN_OVERLAP` | `12` | Minimum overlap (bp) for paired-end merging |
| `MAX_MISMATCH` | `0` | Maximum mismatches allowed in the merge overlap |
| `POOL_METHOD` | `"pseudo"` | DADA2 pooling strategy. `"pseudo"` increases sensitivity; `"independent"` is faster |

### ASV collapse parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ENABLE_ASV_COLLAPSE` | `FALSE` | Set to `TRUE` to run step 4 (`04_collapse.R`). Strongly recommended for species-level analyses |
| `USE_COLLAPSED_ASVS` | Same as `ENABLE_ASV_COLLAPSE` | Whether downstream steps (taxonomy, filter, statistics) use the collapsed feature table |
| `VSEARCH_IDENTITY` | `0.999` | Minimum identity threshold for vsearch pre-clustering (99.9%) |
| `HC_SNP_THRESHOLD` | `1` | SNP distance cutoff for complete-linkage hierarchical clustering |
| `HC_METHOD` | `"complete"` | Linkage method for `hclust()`. `"complete"` is conservative and prevents chaining |
| `MAFFT_THREADS` | `8` | Threads for MAFFT alignment |

### Taxonomy parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `SILVA_TRAINSET` | `""` | Full path to the SILVA 138.2 training set FASTA (gzipped accepted) |
| `SILVA_SPECIES` | `""` | Full path to the SILVA species-assignment FASTA |
| `HOMD_TRAINSET` | `""` | Full path to the HOMD genus trainset. Leave empty to skip the HOMD tier |
| `REFSEQ_BLASTDB` | `""` | Path prefix for the NCBI RefSeq RNA BLAST database (e.g., `/path/to/16S_ribosomal_RNA`) |
| `TAX_N_JOBS` | `10` | Number of parallel jobs for taxonomy and BLAST chunks |
| `TAX_CHUNK_SIZE` | `250` | Sequences per chunk for parallel processing |
| `KEEP_TAX_INTERMEDIATES` | `FALSE` | Retain per-chunk intermediate files for debugging |

### Filtering parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FEATURE_LEVEL` | `"asv"` | Feature level for the final table. Use `"species"` to aggregate by species label before statistical tests |
| `KEEP_KINGDOMS` | `c("Bacteria", "Archaea")` | Kingdoms retained after the kingdom filter. Eukaryota and Viruses are removed by default |
| `KEEP_UNASSIGNED_KINGDOM` | `TRUE` | Retain features with no kingdom assignment (unresolved but potentially prokaryotic) |
| `MIN_PREVALENCE` | `2` | Minimum number of samples in which a feature must appear |
| `MIN_TOTAL_ABUNDANCE_UNIQUE` | `10` | Minimum total count across all samples (used with `FEATURE_LEVEL="species"`) |
| `ASV_LENGTH_MIN` | Auto from `primers.tsv` | Minimum acceptable ASV length in bp. Automatically set from the primer catalog if not specified |
| `ASV_LENGTH_MAX` | Auto from `primers.tsv` | Maximum acceptable ASV length in bp |

### Statistical analysis parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DA_PSEUDOCOUNT` | `1e-5` | Pseudocount added before log2FC calculation on TSS-normalized abundances |
| `DA_MIN_PREV` | `2` | Minimum prevalence for features included in differential abundance testing |
| `PERMANOVA_NPERM` | `999` | Number of permutations for `vegan::adonis2()` |
| `RF_NTREE` | `500` | Number of trees in the Random Forest |
| `RF_NFOLD` | `5` | Folds for stratified cross-validation |
| `TOP_N_RF` | `20` | Top features from Random Forest passed to LDA |

---

## Metadata CSV format

The metadata file must contain at minimum:

| Column | Requirement | Example |
|--------|-------------|---------|
| Sample identifier (`SAMPLE_ID_COL`) | Unique per sample; must match FASTQ filename prefix | `Sample_ID` |
| Group label (`GROUP_COL`) | Case or control label; values must match `CASE_LABEL`/`CTRL_LABEL` | `Group` |
| Sample type (`SAMPLE_TYPE_COL`) | Optional but required if `SAMPLE_TYPE_FILTER` is used | `sample_type` |

**Example metadata CSV:**

```
Sample_ID,Group,Disease,sample_type,Sex,Age
P001,Case,Gastric cancer,Biopsy,M,62
P002,Control,Non-atrophic gastritis,Biopsy,F,55
P003,Case,Gastric cancer,Saliva,M,68
```

CSV files with `.csv` extension are read with `sep=","`. All other extensions are treated as TSV (`sep="\t"`).

> **Note:** Rows with `Sample_ID == "Unk"` or any identifier that does not match a FASTQ filename are excluded silently. Verify that identifiers are consistent before running the pipeline.

---

## Supported primer kits and regions

The pipeline reads primer sequences from `scripts/primers.tsv`. Supported combinations:

| KIT | REGION |
|-----|--------|
| `canonical` | V1V2, V2V3, V3V4, V4, V4V5, V5V7, V7V9 |
| `EMP` | V4, V4V5 |
| `Illumina` | V3V4 |
| `QIAseq` | V1V2, V2V3, V3V4, V4V5, V5V7, V7V9, ITS1 |
| `Zymo` | V1V2, V1V3, V4 |

Set `KIT` and `REGION` exactly as shown above. These values are case-sensitive.

---

Next: [Usage and output structure](usage.md)
