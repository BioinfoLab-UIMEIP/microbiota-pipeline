# 16S Microbiota Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A54.2-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Linux-lightgrey.svg)](#)
[![conda](https://img.shields.io/badge/install%20with-conda-green.svg)](#requirements-r-environment-conda)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20185069.svg)](https://doi.org/10.5281/zenodo.20185069)

> Systematic identification of disease-associated microbial species from 16S rRNA amplicon data in case-control studies.

[Installation](docs/installation.md) · [Configuration](docs/configuration.md) · [Usage](docs/usage.md) · [Methods](docs/methods.md) · [Troubleshooting](docs/troubleshooting.md)

---

## About

A modular R pipeline for case-control 16S amplicon microbiota analysis. Starting from raw paired-end FASTQ files, it covers quality profiling, primer trimming, DADA2 denoising, optional near-identical ASV collapsing, and taxonomy assignment through a SILVA → specialized-database → RefSeq cascade, then produces alpha and beta diversity analyses, compositional summaries, differential abundance testing, Random Forest and LDA-based informative-feature analysis.

---

## Quick-start

1. Install Miniconda3 on macOS 12+ or Linux.
2. Create and activate the conda environment (see [Requirements](#requirements-r-environment-conda)).
3. Install R packages from `conda-forge` and `bioconda`.
4. Download SILVA 138.2 and RefSeq RNA databases (see [Reference databases](#reference-databases)).
5. Copy the config template to your project directory and fill in paths and labels:

```bash
cp /path/to/pipeline/config/project_template.R /path/to/my_project/my_study.R
```

Minimal required variables in `my_study.R`:

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

6. Run the pipeline:

```bash
# Full run from step 1
./run_pipeline.sh --config /path/to/my_project/my_study.R

# Resume from a specific step
./run_pipeline.sh --config /path/to/my_project/my_study.R --from 6 --outdir /path/to/results_YYYYMMDD_HHMMSS

# Run a single step manually
Rscript scripts/07_alpha.R /path/to/my_project/my_study.R --outdir /path/to/results_YYYYMMDD_HHMMSS
```

7. Review outputs: `flags/`, `alpha/`, `beta/`, `differential/`, `informative/`, `report.html`.

---

## Requirements: System

- macOS 12+ (Monterey or later) recommended
- Conda / Miniconda3 for R and external tools
- At least 16 GB RAM; 32 GB recommended for large DADA2 runs
- Approximately 50 GB of free disk space for reference databases

## Requirements: R environment (conda)

Create and activate the conda environment:

```bash
conda create -n microbiota_pipeline r-base=4.3 -c conda-forge -y
conda activate microbiota_pipeline
```

Install Bioconductor packages:

```bash
conda install -c conda-forge -c bioconda \
  bioconductor-dada2 \
  bioconductor-biostrings \
  bioconductor-shortread \
  -y
```

Install CRAN packages from conda-forge:

```bash
conda install -c conda-forge \
  r-ggplot2 \
  r-dplyr \
  r-tidyr \
  r-vegan \
  r-ggrepel \
  r-gridextra \
  r-randomforest \
  r-mass \
  r-scales \
  r-patchwork \
  r-rmarkdown \
  r-knitr \
  r-plotly \
  r-dt \
  r-htmltools \
  r-flexdashboard \
  r-ggvenn \
  -y
```

`parallel`, `grDevices`, and `stats` are part of base R and do not need separate installation.

## Requirements: External binaries

All external tools can be installed into the same `microbiota_pipeline` conda environment used for R.

```bash
conda install -c bioconda -c conda-forge \
  cutadapt vsearch mafft snp-dists blast \
  -y
```

| Tool | Used in | Install | Verify |
|------|---------|---------|--------|
| `cutadapt` >= 4.0 | `02_demux.R` for primer trimming | `conda install -c bioconda cutadapt` | `cutadapt --version` |
| `vsearch` >= 2.22 | `04_collapse.R` for ASV clustering | `conda install -c bioconda vsearch` | `vsearch --version` |
| `MAFFT` >= 7.5 | `04_collapse.R` for multiple sequence alignment | `conda install -c bioconda mafft` | `mafft --version` |
| `snp-dists` >= 0.8 | `04_collapse.R` for SNP distance matrices | `conda install -c bioconda snp-dists` | `snp-dists -v` |
| `BLAST+` / `blastn` >= 2.13 | `05_tax.R` for RefSeq taxonomy refinement | `conda install -c bioconda blast` | `blastn -version` |
| `bsdtar` | `05_tax.R` for archive extraction | built into macOS via `libarchive` | `bsdtar --version` |

## Reference databases

The taxonomy cascade is configured in the project config file. Each database path must be assigned explicitly.

| Database | Used in | Config variable | Download |
|----------|---------|-----------------|----------|
| SILVA 138.2 NR99 trainset + species | `05_tax.R` | `SILVA_TRAINSET`, `SILVA_SPECIES` | <https://zenodo.org/record/8392695> |
| Specialized database (optional; e.g., HOMD for oral/gastric studies) | `05_tax.R` | `HOMD_TRAINSET`, `HOMD_SPECIES` | site-specific (e.g., <https://www.homd.org/ftp/>) |
| RefSeq RNA BLAST database | `05_tax.R` | `REFSEQ_BLASTDB` | NCBI `update_blastdb.pl 16S_ribosomal_RNA` |

Example SILVA configuration:

```r
SILVA_TRAINSET <- "/path/to/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
SILVA_SPECIES  <- "/path/to/silva_v138.2_assignSpecies.fa.gz"
```

For RefSeq RNA, download and format the BLAST database in advance, for example:

```bash
update_blastdb.pl refseq_rna
```

Then point `REFSEQ_BLASTDB` to the local `refseq_rna` database prefix in your config.

## Configuration

The pipeline is driven by an R config file sourced by every step. The intended template path is `config/project_template.R`; if that file is not present in your checkout, use the same variable pattern shown below or copy one of the existing example configs under `config/`. The standard entrypoint is:

```bash
./run_pipeline.sh --config <config.R> [--from N] [--outdir <dir>]
```

Required variables are:

- `RAW_DIR`
- `RESULTS_BASE`
- `KIT`
- `REGION`
- `METADATA_FILE`
- `SAMPLE_ID_COL`
- `GROUP_COL`
- `CASE_LABEL`
- `CTRL_LABEL`

Minimal config example:

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

CUTADAPT_BIN  <- "cutadapt"
VSEARCH_BIN   <- "vsearch"
MAFFT_SMART_BIN <- "mafft"
SNP_DISTS_BIN <- "snp-dists"
BLASTN_BIN    <- "blastn"
BSDTAR_BIN    <- "bsdtar"

SILVA_TRAINSET <- "/path/to/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
SILVA_SPECIES  <- "/path/to/silva_v138.2_assignSpecies.fa.gz"
HOMD_TRAINSET  <- "/path/to/HOMD_trainset_toGenus.fasta"
HOMD_SPECIES   <- "/path/to/HOMD_assignSpecies_dada2.fa"
REFSEQ_BLASTDB <- "/path/to/refseq_rna"
```

## Running the pipeline

Typical commands:

```bash
./run_pipeline.sh --config config/my_config.R
./run_pipeline.sh --config config/my_config.R --from 6
./run_pipeline.sh --config config/my_config.R --from 7 --outdir /path/to/results
```

Pipeline steps:

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `01_qc.R` | Quality profile visualisation and truncation-length suggestion |
| 2 | `02_demux.R` | Primer trimming and paired-end demultiplexing with `cutadapt` |
| 3 | `03_dada2.R` | DADA2 filtering, denoising, merging, and chimera removal |
| 4 | `04_collapse.R` | Optional clustering/collapse of near-identical ASVs |
| 5 | `05_tax.R` | Taxonomy assignment using SILVA -> HOMD -> RefSeq |
| 6 | `06_filter.R` | Length, kingdom, prevalence, and feature-level filtering |
| 7 | `07_alpha.R` | Alpha diversity metrics and group comparisons |
| 8 | `08_beta.R` | Bray-Curtis and Jaccard beta diversity, PCoA, PERMANOVA, betadisper |
| 9 | `09_composition.R` | Taxonomic composition, prevalence, and shared/unique feature summaries |
| 10 | `10_differential.R` | Differential abundance analysis between case and control groups |
| 11 | `11_informative.R` | Random Forest feature ranking and LDA visualisation |
| 12 | `12_report.R` | Rendering of the self-contained HTML dashboard |

## Output structure

Each run writes into a timestamped directory under `RESULTS_BASE`. The following outputs are created:

- `flags/flags_NN_<step>.tsv`: per-step status flags and run diagnostics
- `qc/`: quality profile PDFs and truncation suggestions
- `manifest.tsv`: auto-generated sample manifest when `MANIFEST_FILE` is not provided
- `demux/trim_summary.tsv`: cutadapt retention summary
- `dada2/asv_table.tsv`: raw ASV abundance table
- `dada2/rep_seqs.fasta`: representative ASV sequences
- `taxonomy/`: taxonomy assignments and optional redundancy-collapse outputs
- `filter/asv_table_filtered.tsv`: filtered feature table
- `filter/asv_table_tss.tsv`: total-sum-scaled feature table
- `alpha/alpha_diversity.tsv` and `alpha/alpha_diversity.pdf`: alpha diversity tables and plots
- `beta/permanova.tsv`, `beta/betadisper.tsv`, and `beta/pcoa.pdf`: beta diversity statistics and ordinations
- `composition/composition.pdf`: stacked composition plots, plus prevalence and shared-feature tables
- `differential/da_results.tsv`, `differential/volcano.pdf`, and `differential/ma.pdf`: differential abundance results and visual summaries
- `informative/rf_importance.tsv`, `informative/rf_importance_all.pdf`, `informative/rf_importance_filtered.pdf`, `informative/lda_scores.tsv`, and `informative/lda_scores.pdf`: informative-feature outputs
- `report.html`: self-contained interactive dashboard


---
## Citation

If you use jdist, please cite:

Torres, R.C., Meléndez-Sánchez D., Payró-González M., Almaguer-Molina D., Torres J. Protocol for systematic identification of disease-associated microbial species from 16S rRNA amplicon data in case-control microbiome studies. (in preparation).

---
## License

Released under the MIT License.

---
## Contact

Roberto C. Torres
Medical Research Unit on Infectious and Parasitic Diseases (UIMEIP)
Instituto Mexicano del Seguro Social (IMSS) Mexico City, Mexico

