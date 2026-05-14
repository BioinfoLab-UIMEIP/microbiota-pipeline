---
layout: page
title: "Usage"
nav_order: 4
---

# Usage

## Prerequisites

- conda environment activated: `conda activate microbiota16s`
- Project configuration file created and edited (see [Configuration](configuration.md))
- FASTQ files accessible at `RAW_DIR` or via `MANIFEST_FILE`
- Reference databases downloaded and paths set in the config

---

## Running the full pipeline

```bash
bash /path/to/pipeline/run_pipeline.sh \
  --config /path/to/my_project/my_study.R
```

`run_pipeline.sh` will:
1. Read `RESULTS_BASE` from the config
2. Create a timestamped run directory: `<RESULTS_BASE>/results_YYYYMMDD_HHMMSS/`
3. Execute steps 1–12 in sequence, stopping at the first failure

---

## Resuming from a specific step

If a step fails or you need to re-run from a specific point, use `--from` together with `--outdir` pointing to the existing run directory:

```bash
bash /path/to/pipeline/run_pipeline.sh \
  --config /path/to/my_project/my_study.R \
  --from 7 \
  --outdir /path/to/results_20260505_120000
```

`--outdir` is required for `--from > 1`. The pipeline will skip steps 1–6 and continue from step 7 inside the existing run directory.

> **Note:** `--outdir` is ignored when `--from 1`. In that case, a fresh timestamped directory is always created.

---

## Running individual scripts

Each R script can be called directly outside of `run_pipeline.sh`. This is useful for debugging or interactive exploration:

```bash
Rscript scripts/01_qc.R /path/to/my_project/my_study.R \
  --outdir /path/to/results_20260505_120000

Rscript scripts/05_tax.R /path/to/my_project/my_study.R \
  --outdir /path/to/results_20260505_120000
```

Both `<config.R>` and `--outdir <dir>` are required for all scripts.

---

## Command-line reference

```
Usage: run_pipeline.sh [OPTIONS]

Options:
  --config <path>   Path to project config R file (required)
  --outdir <path>   Existing results_* run directory (required when --from > 1)
  --from <N>        Start from step N (1-12, default: 1)
  --help            Show this help message
```

---

## Output directory structure

For each run, the pipeline creates:

```
<RESULTS_BASE>/
  results_YYYYMMDD_HHMMSS/
    manifest.tsv               # Sample manifest (sample_id, read1, read2)
    report.html                # Final interactive HTML dashboard
    qc/
      qc_report.pdf            # Per-cycle quality profiles (R1 and R2)
      trunc_suggestion.txt     # Suggested TRUNC_LEN_R1 and TRUNC_LEN_R2
    demux/
      <sample_id>_R1.fastq.gz  # Primer-trimmed FASTQ pairs
      <sample_id>_R2.fastq.gz
      trim_summary.tsv         # Per-sample cutadapt retention percentages
    dada2/
      asv_table.tsv            # Raw ASV count table (samples x ASVs)
      rep_seqs.fasta           # Representative sequences for all ASVs
      track.tsv                # Per-sample read counts at each DADA2 stage
    taxonomy/
      taxonomy_ASVs.tsv        # Full taxonomy table (or taxonomy_nrASVs.tsv if collapsed)
      taxonomy_assigned_collapsed.tsv  # Taxonomy after filtering (used by steps 7-12)
      taxonomy_performance.pdf         # Stacked bar: assignment proportion by database tier
      redundancy/              # ASV collapse intermediates (present only if ENABLE_ASV_COLLAPSE=TRUE)
        asv_table_collapsed.tsv
        rep_seqs_collapsed.fasta
        cluster_map.tsv
        cluster_representatives.tsv
        centroids_hc_cutoff.pdf
    filter/
      asv_table_filtered.tsv   # Count table after length, kingdom, and prevalence filters
      asv_table_tss.tsv        # Total-sum-scaled (TSS) relative-abundance table
      taxonomy_filtered.tsv    # Taxonomy table matching filtered features
    alpha/
      alpha_diversity.tsv      # Per-sample richness, evenness, Shannon, Simpson
      alpha_boxplots.pdf       # Group comparison boxplots with Wilcoxon p-values
    beta/
      pcoa_braycurtis.pdf      # Bray-Curtis PCoA plot
      pcoa_jaccard.pdf         # Jaccard PCoA plot
      permanova_results.tsv    # adonis2 R2 and p-values
      betadisper_results.tsv   # Dispersion homogeneity test results
    composition/
      phylum_barplot.pdf       # Mean relative abundance by phylum
      genus_barplot.pdf        # Top genera by group
      species_barplot.pdf      # Top species by group
      prevalence_summary.tsv   # Per-feature prevalence by group
    da/
      da_results.tsv           # Differential abundance results (all features)
      da_significant.tsv       # Significant features only (BH p < 0.05)
      volcano_plot.pdf         # Volcano plot (log2FC vs -log10 BH p-value)
    informative/
      rf_importance.tsv        # Random Forest mean-decrease-accuracy per feature
      rf_importance_plot.pdf   # Importance bar plot
      lda_coefficients.tsv     # LDA LD1 coefficients
      lda_plot.pdf             # LD1 directionality plot
      rf_cv_accuracy.tsv       # Cross-validated accuracy per fold
    flags/
      flags_01_qc.tsv
      flags_02_demux.tsv
      flags_03_dada2.tsv
      flags_04_collapse.tsv
      flags_05_tax.tsv
      flags_06_filter.tsv
      flags_07_alpha.tsv
      flags_08_beta.tsv
      flags_09_composition.tsv
      flags_10_differential.tsv
      flags_11_informative.tsv
      flags_12_report.tsv
```

---

## Flags TSV format

Every script writes a flags file in `flags/` when it completes. The format is:

```
step    key                    value        status
01_qc   n_samples_raw          120          ok
01_qc   n_samples_after_filter 60           ok
01_qc   sample_type_filter     Biopsy       ok
01_qc   trunc_r1               270          ok
01_qc   trunc_r2               210          warn
01_qc   out_dir                /path/...    ok
```

Status values:
- `ok` — value is within expected range
- `warn` — value is outside a soft threshold (e.g., truncation length = 0, retention < 70%, merge < 30%)
- `error` — a failure condition that will cause the pipeline to abort

To check the overall status of a completed run:

```bash
grep -h "warn\|error" /path/to/results_*/flags/flags_*.tsv
```

---

Next: [Analytical methods](methods.md)
