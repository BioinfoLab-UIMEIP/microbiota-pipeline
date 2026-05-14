---
layout: page
title: "Methods"
nav_order: 5
---

# Analytical Methods

This page describes the biological question answered at each pipeline step, the methodological choice, and the key parameters.

---

## Step 1: Quality profiling and manifest building (`01_qc.R`)

**Biological question:** Are the raw reads of sufficient quality to support DADA2 error modelling and paired-end merging?

**Method:** `ShortRead::FastqSampler` draws 5,000 reads from each of up to 8 randomly selected libraries (controlled by `QC_N_SAMPLES` and `QC_READS_PER_FILE`). Per-cycle median Phred scores and interquartile ranges are computed. Truncation lengths are suggested as the last cycle where median Phred remains at or above `Q_THRESHOLD` (default Q25), clipped to the 5th percentile of observed read lengths. A length histogram and quality profile PDF are written to `qc/`.

If `SAMPLE_TYPE_FILTER` is set, only samples belonging to the specified tissue types are included in the manifest and all downstream steps. This is the single point where a tissue-type subset is selected.

---

## Step 2: Primer trimming (`02_demux.R`)

**Biological question:** Have the primer sequences been removed so that DADA2 error models are not distorted by primer heterogeneity?

**Method:** `cutadapt` is called with `--discard-untrimmed`, error rate 0.1, and minimum post-trim length 100 bp. Primer sequences are read from `scripts/primers.tsv` by matching `KIT` and `REGION`. Read-through adapters (present in some kits such as Illumina Nextera XT or QIAseq) are trimmed in the same call using the `-a`/`-A` adapter flags. Retention percentages are recorded per sample.

A `warn` flag is set for any sample retaining fewer than 70% of reads, and for the overall median when it falls below 70%.

---

## Step 3: DADA2 denoising (`03_dada2.R`)

**Biological question:** What are the exact amplicon sequence variants (ASVs) present in this dataset, and how abundant is each one in each sample?

**Method:** DADA2 v1.26+ implements an error-correcting denoising model. The pipeline applies:
- `filterAndTrim()`: `maxEE = c(2, 2)`, `maxN = 0`, `rm.phix = TRUE`
- Error models learned from a random subsample of reads using `learnErrors()`
- `dada()` with `pool = "pseudo"` — a two-pass approach that increases sensitivity for reproducible low-abundance variants while controlling runtime
- Paired-end merging with `MIN_OVERLAP = 12` bp and `MAX_MISMATCH = 0`
- Chimera removal with `removeBimeraDenovo(method = "consensus")`

Pseudo-pooling was chosen over full pooling because it provides comparable sensitivity for rare ASVs with substantially lower memory requirements for moderate-sized datasets. For very large studies (>200 samples), independent sample processing may be preferred.

Truncation lengths are read from `qc/trunc_suggestion.txt` if `TRUNC_LEN` is not set explicitly in the config.

> **Truncation guidance for V4-V5 (520 bp amplicon):** Typical starting values are `TRUNC_LEN = c(270, 210)`. The merged length must remain within the expected amplicon range (~340-420 bp for V4-V5). At minimum, the sum of truncated lengths minus the amplicon length must be at least 12 bp (the `MIN_OVERLAP` requirement).

---

## Step 4: ASV collapse (`04_collapse.R`)

**Biological question:** Are there near-identical ASVs that fragment the same biological entity across multiple table rows, reducing species-level statistical power?

**Rationale:** DADA2 inference is exact-sequence, which maximises resolution. However, single-nucleotide sequencing errors that escape the DADA2 model, genuine intraspecific variation within a 16S hypervariable region, and slight truncation differences can produce multiple ASVs for the same organism. This step consolidates such near-identical redundancy before taxonomy and statistics.

**Method (four-stage collapse):**

1. **vsearch pre-clustering at 99.9% identity:** Sequences sorted by decreasing total abundance are clustered with `--cluster_fast`. This groups the vast majority of near-identical ASVs efficiently.

2. **MAFFT multiple alignment:** Vsearch centroid sequences (cluster representatives) are aligned with MAFFT using the configured mode (set via `MAFFT_MODE` in the config; default: `"smart"` for auto-selection).

3. **SNP-distance matrix:** `snp-dists` computes pairwise SNP distances from the multiple alignment.

4. **Complete-linkage hierarchical clustering at SNP = 1:** `hclust(method = "complete")` is applied to the distance matrix. `cutree()` at `HC_SNP_THRESHOLD = 1` defines the final collapse groups. Within each group, the most abundant representative sequence is retained. All counts are summed onto that representative.

The `ENABLE_ASV_COLLAPSE` flag must be set to `TRUE` in the config. When `USE_COLLAPSED_ASVS = TRUE`, steps 5–12 operate on the collapsed feature table (`asv_table_collapsed.tsv` and `rep_seqs_collapsed.fasta`).

> **CRITICAL:** Do not skip ASV collapse for species-level analyses. Studies in which collapse was applied showed 5–30% reduction in feature count without any loss of total abundance, resulting in stronger per-species statistical signal.

---

## Step 5: Taxonomic assignment cascade (`05_tax.R`)

**Biological question:** What organisms do the representative sequences belong to, resolved to the best taxonomic depth available from current reference databases?

**Rationale:** No single reference database covers all environments equally. SILVA is comprehensive but may lack species-level resolution for certain ecological niches. A site-specific specialized database (e.g., HOMD for oral/gastric studies) provides higher resolution for those taxa. BLAST against RefSeq RNA captures taxa absent from curated trainsets.

**Method (three-tier cascade):**

1. **SILVA 138.2 (mandatory):** `dada2::assignTaxonomy()` with `minBoot = 50` and `tryRC = TRUE` assigns Kingdom through Genus. `dada2::addSpecies()` with the SILVA species FASTA adds species when an exact match exists. Sequences fully resolved to species at this stage are not passed to subsequent tiers.

2. **Site-specific database (optional):** Features not resolved to species by SILVA are queried against `HOMD_TRAINSET` (set to a specialized trainset appropriate for your sample type) using the same DADA2 classifier. Assignments are grafted onto the SILVA lineage when SILVA already resolved to genus, extending resolution one level deeper.

3. **NCBI RefSeq RNA BLAST (residuals):** Features still unresolved after tiers 1–2 are queried with `blastn -task megablast -max_target_seqs 5 -max_hsps 1 -evalue 1e-20`. The top hit scientific name is used to construct a species label grafted onto the existing SILVA lineage when possible. BLAST runs in parallel across sequence chunks (`TAX_N_JOBS` chunks of `TAX_CHUNK_SIZE` sequences each).

4. **UNITE (optional, for fungi):** Completely unresolved features can optionally be queried against a UNITE ITS trainset if `UNITE_TRAINSET` is set. Primarily relevant for mixed-kingdom datasets.

A cascade logic (`choose_better`) always retains the assignment with the deepest resolved lineage, so no information is discarded.

**Outputs:**
- `taxonomy_ASVs.tsv` (or `taxonomy_nrASVs.tsv` for collapsed mode): full taxonomy table with columns `ASV_ID`, `nr_group`, Kingdom–Species, `dada_species`, `refseq_species`, `unite_species`, `Species`, `species_source`, `assignment_level`, `source_db`
- `taxonomy_performance.pdf`: stacked bar chart showing proportion of features resolved per database tier

---

## Step 6: Filtering (`06_filter.R`)

**Biological question:** Which features are of appropriate length, bacterial origin, and sufficient prevalence to be included in ecological and statistical analyses?

**Method:**

1. **Length filter:** ASV lengths are extracted from the taxonomy table's `sequence` column. Features outside `[ASV_LENGTH_MIN, ASV_LENGTH_MAX]` are removed. Default values are read from `scripts/primers.tsv` based on the expected amplicon range for the configured KIT/REGION combination.

2. **Kingdom filter:** Features assigned to kingdoms not in `KEEP_KINGDOMS` (default: `c("Bacteria", "Archaea")`) are removed. Eukaryota, Viruses, and host sequences are excluded. Features with no kingdom assignment are retained by default (`KEEP_UNASSIGNED_KINGDOM = TRUE`).

3. **Feature-level collapse (optional):** When `FEATURE_LEVEL = "species"`, all features sharing the same species label are aggregated by summing their counts. The most abundant representative sequence is retained as the row identity.

4. **Prevalence filter:** Features present in fewer than `MIN_PREVALENCE` samples are removed. At `FEATURE_LEVEL = "species"`, a minimum total abundance of `MIN_TOTAL_ABUNDANCE_UNIQUE` counts is also applied.

5. **TSS normalization:** A total-sum-scaled (TSS) relative-abundance table is computed by dividing each sample's counts by its total.

Outputs:
- `filter/asv_table_filtered.tsv` — count table after all filters
- `filter/asv_table_tss.tsv` — TSS relative-abundance table (input for beta diversity and differential abundance)
- `filter/taxonomy_filtered.tsv` — filtered taxonomy table

---

## Step 7: Alpha diversity (`07_alpha.R`)

**Biological question:** Do case and control groups differ in within-sample species richness or evenness?

**Method:** Four metrics are computed per sample from the filtered count table:
- **Observed richness:** count of features with abundance > 0
- **Pielou evenness:** Shannon H' / ln(richness); 0 = one dominant feature, 1 = perfectly even
- **Shannon H':** `-sum(p * ln(p))` across non-zero features
- **Simpson 1-D:** `1 - sum(p^2)`

Groups are compared with Wilcoxon rank-sum tests (two-sided). P-values are adjusted with the Benjamini-Hochberg procedure. Boxplots with jittered points and annotated Wilcoxon p-values are written to `alpha/`.

---

## Step 8: Beta diversity (`08_beta.R`)

**Biological question:** Do the overall community compositions differ between case and control groups?

**Method:**

1. **Dissimilarity matrices:** Bray-Curtis (`vegdist(method="bray")`) and binary Jaccard (`vegdist(method="jaccard", binary=TRUE)`) dissimilarities are computed from the TSS table.

2. **PCoA:** `cmdscale()` with `k=2` and `eig=TRUE`. Percent variance explained by each axis is reported as `eigenvalue / sum(eigenvalues)`.

3. **PERMANOVA:** `vegan::adonis2()` with `PERMANOVA_NPERM` permutations (default 999) tests whether group centroids differ. The R2 statistic is the effect size.

4. **Dispersion test:** `vegan::betadisper()` followed by `permutest()` tests whether the within-group spread (dispersion) differs between groups. A significant dispersion difference should be noted when interpreting PERMANOVA results.

> When betadisper is significant, PERMANOVA significance may partly reflect variance differences rather than centroid separation.

---

## Step 9: Compositional profiling (`09_composition.R`)

**Biological question:** What is the taxonomic composition of each group, and which taxa are uniquely or differentially prevalent?

**Method:** Mean relative abundances by group are computed from the TSS table at phylum, genus, and species levels. The top `TOP_N_PHYLA`, `TOP_N_GENERA`, and `TOP_N_SPECIES` are plotted as stacked or side-by-side bar charts. Prevalence (fraction of samples with count > 0) is summarized per group.

---

## Step 10: Differential abundance (`10_differential.R`)

**Biological question:** Which individual species are significantly more or less abundant in cases versus controls?

**Method:**

- **Test:** Two-sided Wilcoxon rank-sum test on raw counts, applied to all features passing `DA_MIN_PREV`.
- **Multiple testing:** Benjamini-Hochberg FDR correction.
- **Effect size:** log2 fold change computed from TSS-normalized abundances after adding `DA_PSEUDOCOUNT = 1e-5` to zero values.

Wilcoxon was chosen over parametric alternatives because it is distribution-free, robust to the zero-inflation characteristic of microbiome count data, and easy to audit. For confirmatory analyses, results can be compared against ANCOM-BC applied to the same feature table.

**Output:** `da/da_results.tsv` contains Species, log2FC, Wilcoxon W statistic, p-value, BH-adjusted p-value, and mean TSS abundance in each group. The volcano plot visualizes all features with significant hits colored by enrichment direction.

---

## Step 11: Machine-learning prioritization (`11_informative.R`)

**Biological question:** Which of the differentially abundant species contribute most to case-control discrimination, and in which direction?

**Method (two-stage):**

1. **Random Forest (stratified 5-fold cross-validation):**
   - Input: CLR-transformed species counts (zero imputation by half-minimum before log-transform)
   - `randomForest::randomForest()` with `ntree = RF_NTREE` (default 500)
   - Stratified folds preserve case:control ratio across folds
   - Feature importance: mean decrease in accuracy (MDA) averaged over `RF_NFOLD` folds
   - Cross-validated accuracy is reported in the flags file

2. **LDA (`MASS::lda()`):**
   - Input: the top `TOP_N_RF` features by mean MDA on CLR-transformed data
   - LD1 coefficients are normalized to a [-100, +100] scale for visualization
   - Positive LD1 aligns with the case group, negative LD1 with control

The combined RF + LDA output ranks species by discriminative importance and assigns directionality. This is the primary output used to nominate candidate taxa for downstream experimental validation.

> **CRITICAL:** LDA requires that the number of features not exceed the effective degrees of freedom (`n_samples - 2`). For studies with fewer than 55 samples, reduce `TOP_N_RF` so that `TOP_N_RF <= n_samples - 3`.

---

## Step 12: Interactive HTML report (`12_report.R`)

**Purpose:** Compile all pipeline outputs into a single, self-contained HTML dashboard for integrated review, interpretation, and sharing.

**Method:** `rmarkdown::render()` processes `scripts/report_template.Rmd` using `flexdashboard` layout with `plotly` interactive panels and `DT` searchable tables. The report integrates outputs from all previous steps. It is self-contained HTML and can be shared without a web server.

---

## Statistical choices and their rationale

| Choice | Rationale |
|--------|-----------|
| Wilcoxon DA testing | Distribution-free; robust to zero-inflation and non-normality common in 16S count data |
| TSS normalization | Simple and widely used; avoids assumptions about library-size composition |
| CLR transformation for RF/LDA | Compositionally aware; centers log-ratios to remove spurious negative correlations |
| Complete-linkage HC for ASV collapse | Conservative; prevents chaining of distant sequences into the same cluster |
| Pseudo-pooling in DADA2 | Better sensitivity for low-abundance variants than per-sample mode, with lower memory than full pooling |
| BH FDR correction | Standard for exploratory discovery; controls expected proportion of false positives in the candidate set |

---

Next: [Troubleshooting](troubleshooting.md)
