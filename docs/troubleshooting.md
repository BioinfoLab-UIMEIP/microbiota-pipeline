---
layout: page
title: "Troubleshooting"
nav_order: 6
---

# Troubleshooting

## Quick diagnostics

Before investigating individual problems, check the flags files for the step that failed:

```bash
# Show all warn/error entries across all steps
grep -h "warn\|error" /path/to/results_*/flags/flags_*.tsv

# Show flags for a specific step
cat /path/to/results_20260505_120000/flags/flags_03_dada2.tsv
```

---

## Common problems and solutions

### 1. `primers.tsv` not found at runtime

**Symptom:** `02_demux.R` stops with `No entry for kit=... region=... in primers.tsv` or cannot locate the file.

**Possible cause:** The script resolves `primers.tsv` relative to its own location (`scripts/`). Running from a different working directory, or referencing the scripts by an absolute path that differs from the pipeline installation path, can change the expected lookup location.

**Solution:**
- Confirm that `scripts/primers.tsv` exists in the `scripts/` directory alongside the R scripts.
- If running from a non-standard directory, set the environment variable `SCRIPTS_DIR` to the absolute path of the `scripts/` directory before running:
  ```bash
  export SCRIPTS_DIR=/path/to/pipeline/scripts
  Rscript /path/to/pipeline/scripts/02_demux.R /path/to/my_study.R --outdir /path/to/outdir
  ```

---

### 2. `manifest.tsv` not found by steps 03–12 when running `run_pipeline.sh`

**Symptom:** Steps 3 through 12 fail immediately with a message that `manifest.tsv` cannot be found.

**Possible cause:** When running scripts manually (not via `run_pipeline.sh`), `--outdir` must point to the full timestamped run directory (`results_YYYYMMDD_HHMMSS/`), not to `RESULTS_BASE`. `run_pipeline.sh` handles this automatically, but manual runs require the full path.

**Solution:**
- Pass the full path to the timestamped subdirectory:
  ```bash
  Rscript scripts/03_dada2.R /path/to/my_study.R \
    --outdir /path/to/results_20260505_120000
  ```
- Do not pass `RESULTS_BASE` as `--outdir` for steps 2–12.

---

### 3. Very low merge retention (< 5%) in most samples

**Symptom:** `flags_03_dada2.tsv` reports a high number of `low_merge` samples. The `track.tsv` shows the majority of reads are lost at the merging step.

**Possible causes:**
1. Truncation lengths are too short to preserve a 12 bp overlap between forward and reverse reads.
2. The library is dominated by non-specific amplification or host DNA (common in low-biomass biopsy samples).

**Solutions:**
1. Inspect `qc/qc_report.pdf` and verify that the truncation settings preserve at least 12 bp of overlap. For V4-V5 with a ~380 bp insert, `TRUNC_LEN = c(270, 210)` gives 270 + 210 - 380 = 100 bp of overlap — well above the minimum. If quality drops earlier, try `c(250, 200)`.
2. If overlap is adequate and merge rates remain low, the problem is biological: low microbial biomass and non-specific amplification. Consider excluding samples with fewer than 1,000 merged reads. Do not increase `MAX_MISMATCH` or reduce `MIN_OVERLAP` below 10 bp to compensate.

---

### 4. BLAST runs extremely slowly or appears to hang

**Symptom:** Step 5 (`05_tax.R`) is running for many hours with no apparent progress on the BLAST tier. System CPU is idle or not peaking.

**Possible cause:** BLAST is running serially over many feature chunks, the database is very large, or `blastn` is not found on `PATH` inside the conda environment.

**Solutions:**
- Increase `TAX_N_JOBS` in the config to parallelize BLAST over more CPU cores (e.g., `TAX_N_JOBS <- 16`).
- Confirm that `blastn` is on `PATH` inside the active conda environment: `which blastn`.
- Refresh the database if it was downloaded long ago: `update_blastdb.pl --decompress 16S_ribosomal_RNA`.
- If the BLAST database is on a network file system, move it to local storage to avoid I/O latency.
- Set `KEEP_TAX_INTERMEDIATES = TRUE` in the config; if step 5 is re-run, completed chunks are read from cache and only incomplete chunks are repeated.

---

### 5. The HTML report fails to render (plotly or geom_violin error)

**Symptom:** Step 12 (`12_report.R`) fails with an error referencing `geom_violin`, a plotly conversion error, or a missing object.

**Possible cause:**
- The `report_template.Rmd` was edited and a `geom_violin` call was reintroduced into an interactive (plotly) panel. Plotly cannot render `geom_violin` directly.
- A group has too few observations for violin rendering.
- An upstream step failed and a required input file is absent.

**Solutions:**
- Use the shipped `report_template.Rmd` without modification. Interactive panels should use `geom_boxplot + geom_jitter`, not `geom_violin`.
- Check the flags files for all prior steps: `grep -h "error" flags/flags_*.tsv`. An empty or corrupted upstream output (e.g., `da_results.tsv` missing) will cause the report to fail on that panel.
- If only one panel fails, set `eval=FALSE` for that code chunk in the template temporarily to identify the offending section.

---

### 6. Sample count drops unexpectedly between steps 1 and 7

**Symptom:** The manifest built in step 1 contains N samples, but by step 7 (alpha diversity) the merged data frame has substantially fewer samples.

**Possible causes:**
- `SAMPLE_TYPE_FILTER` excludes samples after the metadata join in steps 7–11.
- Rows with `Sample_ID == "Unk"` (or any ID that does not match a FASTQ filename) are excluded automatically.
- Case-sensitive mismatch between `SAMPLE_TYPE_FILTER` values and the actual values in `SAMPLE_TYPE_COL`.

**Solutions:**
- Check the `n_samples` key in each flags file to track where sample count drops.
- Verify that `SAMPLE_TYPE_FILTER` uses exactly the strings that appear in the metadata column (e.g., `c("Biopsy")` not `c("biopsy")`).
- Inspect sample identifiers: `cut -f1 manifest.tsv | sort > manifest_ids.txt` and compare with the Sample_ID column in the metadata CSV.

---

### 7. LDA fails with "variables are collinear" or "fewer observations than variables"

**Symptom:** Step 11 (`11_informative.R`) stops with an error from `MASS::lda()` about collinear variables or too few observations.

**Possible cause:** The default `TOP_N_RF = 20` features passed to LDA exceeds the effective degrees of freedom for a small dataset. `lda()` requires the number of features to be strictly less than the number of samples minus one per class.

**Solution:** Reduce `TOP_N_RF` in the config so that `TOP_N_RF <= n_samples - 3`. Practical guidelines:
- Fewer than 20 total samples: `TOP_N_RF <- 5`
- 20–30 total samples: `TOP_N_RF <- 10`
- 30–55 total samples: `TOP_N_RF <- 15`
- More than 55 total samples: the default of 20 is safe

---

### 8. SILVA `assignTaxonomy()` is very slow or the R process runs out of memory

**Symptom:** Step 5 stalls for hours on the SILVA tier, or the R session is killed by the OS due to memory exhaustion.

**Possible cause:** The SILVA 138.2 training set is large (approximately 2–4 GB when uncompressed). Systems with less than 8 GB RAM can experience severe slowdown or crash during multithreaded assignment.

**Solutions:**
- Ensure at least 8 GB RAM is available and close other memory-intensive processes before running step 5.
- Reduce `TAX_N_JOBS` to 2–4 to lower peak memory pressure (fewer concurrent `assignTaxonomy` workers).
- Use a gzip-compressed trainset file (`silva_nr99_v138.2_train_set.fa.gz`); DADA2 reads it directly and does not require full decompression to disk.
- Consider running step 5 on a machine with more RAM, or splitting the query into smaller chunks by temporarily reducing `TAX_CHUNK_SIZE`.

---

## Additional diagnostics

### Verify the BLAST database index

```bash
blastdbcheck -db /path/to/databases/16S_ribosomal_RNA -verbosity 1
```

### Inspect DADA2 read tracking

```bash
column -t /path/to/results_*/dada2/track.tsv | head -20
```

### Check cutadapt retention

```bash
column -t /path/to/results_*/demux/trim_summary.tsv
```

### List samples below 30% merge retention

```bash
awk -F'\t' 'NR>1 && $4 < 30' /path/to/results_*/dada2/track.tsv
```

---

If a problem is not covered here, consult the flags files for the failed step and review the R script source code in `scripts/`. Each script prints diagnostic `message()` calls that are written to stderr and captured in the terminal output when running `run_pipeline.sh`.
