---
layout: page
title: "16S Microbiota Pipeline"
nav_order: 1
---

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A54.2-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Linux-lightgrey.svg)](#)
[![conda](https://img.shields.io/badge/install%20with-conda-green.svg)](installation.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20185069.svg)](https://doi.org/10.5281/zenodo.20185069)

**Systematic identification of disease-associated microbial species from 16S rRNA amplicon data in case-control studies**

*Roberto C. Torres, PhD. — Bioinformatics Lab, Infectious Diseases Research Unit (UIMEIP), CMN SXXI, IMSS, Mexico City*

---

## What this pipeline does

This pipeline converts raw paired-end 16S rRNA amplicon FASTQ files into a ranked, directional list of candidate disease-associated microbial species. It is designed for case-control study designs and produces an interactive HTML report alongside TSV tables and PDF plots.

Twelve sequential R scripts orchestrated by `run_pipeline.sh` cover every stage from quality profiling through machine-learning prioritization:

| Step | Script | Function |
|------|--------|----------|
| 1 | `01_qc.R` | Quality profiling, TRUNC_LEN suggestion, manifest building |
| 2 | `02_demux.R` | Primer trimming via cutadapt |
| 3 | `03_dada2.R` | DADA2 denoising with pseudo-pooling and paired-end merging |
| 4 | `04_collapse.R` | Post-denoising ASV collapse (vsearch + MAFFT + snp-dists + HC) |
| 5 | `05_tax.R` | Two-level taxonomic cascade: SILVA → HOMD/site-specific DB → BLAST |
| 6 | `06_filter.R` | Kingdom filter, species-level collapse, prevalence/abundance filter |
| 7 | `07_alpha.R` | Alpha diversity (richness, Pielou, Shannon, Simpson), Wilcoxon + BH |
| 8 | `08_beta.R` | Beta diversity PCoA (Bray-Curtis/Jaccard), PERMANOVA, betadisper |
| 9 | `09_composition.R` | Prevalence and relative-abundance profiles by group |
| 10 | `10_differential.R` | Wilcoxon DA testing, log2FC, volcano plot |
| 11 | `11_informative.R` | Random Forest (stratified 5-fold CV) + LDA prioritization |
| 12 | `12_report.R` | Interactive HTML dashboard (flexdashboard + plotly) |

---

## Key innovations

### 1. Structured discovery path

Each stage answers a specific biological question — whether reads are technically adequate, which features are credible, whether communities differ globally, which species differ between groups, and which taxa drive multifeature discrimination. The result is a ranked candidate list rather than only a community profile.

### 2. Post-denoising ASV collapse

Near-identical ASVs can fragment abundance across multiple table rows and weaken species-level statistical power. This pipeline reduces that redundancy after DADA2 by combining:
- `vsearch` clustering at >=99.9% identity
- `MAFFT` multiple alignment of cluster representatives
- `snp-dists` pairwise SNP-distance calculation
- Complete-linkage hierarchical clustering at SNP threshold = 1
- Abundance summation onto the most abundant representative

### 3. Two-level taxonomic cascade

SILVA 138.2 is the mandatory general-purpose tier. A context-dependent specialized database (HOMD for oral/gastric studies) is applied second. NCBI BLAST against RefSeq RNA resolves residual unclassified features. This staged approach recovers clinically relevant taxa that a single database would leave unresolved.

---

## Quick install

```bash
# 1. Create conda environment
conda create -n microbiota16s -c conda-forge -c bioconda \
  r-base=4.2 cutadapt vsearch mafft snp-dists blast
conda activate microbiota16s

# 2. Install R packages (inside the conda environment)
Rscript -e "
  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')
  BiocManager::install(c('dada2','ShortRead','Biostrings'))
  install.packages(c('vegan','ape','randomForest','MASS','ggplot2','ggrepel',
    'gridExtra','dplyr','tidyr','rmarkdown','flexdashboard','plotly','DT',
    'randomcoloR','reshape2','tibble'))
"
```

See the [full installation guide](installation.md) for reference database setup and verification.

---

## Quick start

```bash
# Copy and edit the project configuration template
cp /path/to/pipeline/config/project_template.R /path/to/my_project/my_study.R
# Edit my_study.R: set RESULTS_BASE, RAW_DIR, METADATA_FILE, database paths,
#                  CASE_LABEL, CTRL_LABEL, KIT, REGION

# Activate the environment
conda activate microbiota16s

# Run the full pipeline
bash /path/to/pipeline/run_pipeline.sh \
  --config /path/to/my_project/my_study.R

# Resume from a specific step (example: restart from step 7)
bash /path/to/pipeline/run_pipeline.sh \
  --config /path/to/my_project/my_study.R \
  --from 7 \
  --outdir /path/to/results_20260505_120000
```

Outputs are written to a timestamped subdirectory: `results_YYYYMMDD_HHMMSS/`

---

## Documentation

- [Installation guide](installation.md)
- [Configuration reference](configuration.md)
- [Usage and output structure](usage.md)
- [Analytical methods](methods.md)
- [Troubleshooting](troubleshooting.md)

---

## Citation

If you use this pipeline, please cite:

> Torres, R.C. (2026). Protocol for systematic identification of disease-associated microbial species from 16S rRNA amplicon data in case-control microbiome studies. *STAR Protocols* (in preparation).

---

## Contact

**Roberto C. Torres, PhD.**  
Bioinformatics Lab, Infectious Diseases Research Unit (UIMEIP)  
Centro Médico Nacional Siglo XXI (CMN SXXI), IMSS  
Mexico City, Mexico  
[torres.roberto.c@gmail.com](mailto:torres.roberto.c@gmail.com)
