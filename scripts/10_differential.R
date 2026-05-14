# =============================================================================
# Script:      10_differential.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Differential abundance analysis between case and control groups
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 10_differential.R <config.R> --outdir <dir>")
source(cfg)
if (!exists("METADATA_FILE") || !nzchar(METADATA_FILE)) stop("METADATA_FILE must be set in config")
if (!exists("GROUP_COL")) GROUP_COL <- "Group"
if (!exists("SAMPLE_ID_COL")) SAMPLE_ID_COL <- "Sample_ID"
if (!exists("CASE_LABEL")) CASE_LABEL <- NULL
if (!exists("CTRL_LABEL")) CTRL_LABEL <- NULL
if (!exists("DA_PSEUDOCOUNT")) DA_PSEUDOCOUNT <- 1e-5
if (!exists("DA_MIN_PREV")) DA_MIN_PREV <- 2L

od_idx  <- which(args == "--outdir")
out_dir <- if (length(od_idx)) path.expand(args[od_idx[1] + 1]) else stop("--outdir <dir> is required")
run_seed <- suppressWarnings(as.integer(gsub("[^0-9]", "", basename(out_dir))))
set.seed(if (!is.na(run_seed) && run_seed > 0L) run_seed else 42L)
pastel_n <- function(n) {
  hues <- seq(15, 375, length.out = n + 1L)[seq_len(n)]
  grDevices::hcl(hues, c = 35, l = 82)
}
group_levels <- c(CTRL_LABEL, CASE_LABEL)
group_colors <- setNames(pastel_n(2L), group_levels)
taxa_palette <- function(n) pastel_n(n)
asv_path <- file.path(out_dir, "filter", "asv_table_filtered.tsv"); tax_path <- file.path(out_dir, "taxonomy", "taxonomy_assigned_collapsed.tsv")
if (!file.exists(asv_path)) stop("asv_table_filtered.tsv not found: ", asv_path)
if (!file.exists(tax_path)) stop("taxonomy_filtered.tsv not found: ", tax_path)
if (!file.exists(path.expand(METADATA_FILE))) stop("Metadata file not found: ", METADATA_FILE)

asv <- read.table(asv_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tss_path <- file.path(out_dir, "filter", "asv_table_tss.tsv")
tss_for_mean <- if (file.exists(tss_path)) {
  t <- read.table(tss_path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE)
  setNames(colMeans(t[, setdiff(colnames(t), "sample_id"), drop = FALSE]),
    setdiff(colnames(t), "sample_id"))
} else NULL
tax <- read.table(tax_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
sep_m <- if (grepl("\\.csv$", METADATA_FILE, ignore.case = TRUE)) "," else "\t"
meta <- read.table(path.expand(METADATA_FILE), sep = sep_m, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample_id" %in% colnames(asv)) stop("asv_table_filtered.tsv must contain sample_id column")
if (!exists("SAMPLE_TYPE_COL")) SAMPLE_TYPE_COL <- "sample_type"
if (GROUP_COL %in% colnames(meta) && !is.null(CASE_LABEL) && !is.null(CTRL_LABEL)) {
  meta <- meta[meta[[GROUP_COL]] %in% c(CASE_LABEL, CTRL_LABEL), , drop = FALSE]
  if (exists("SAMPLE_TYPE_FILTER") && !is.null(SAMPLE_TYPE_FILTER) &&
      length(SAMPLE_TYPE_FILTER) > 0 && nzchar(SAMPLE_TYPE_FILTER[1]) &&
      SAMPLE_TYPE_COL %in% colnames(meta)) {
    meta <- meta[meta[[SAMPLE_TYPE_COL]] %in% SAMPLE_TYPE_FILTER, , drop = FALSE]
    message("Sample type filter applied: ",
            paste(SAMPLE_TYPE_FILTER, collapse = ", "),
            " (", nrow(meta), " samples retained)")
  }
}
merged <- merge(asv, meta, by.x = "sample_id", by.y = SAMPLE_ID_COL)
if (!GROUP_COL %in% colnames(merged)) stop("GROUP_COL '", GROUP_COL, "' not in metadata")
grps <- sort(unique(merged[[GROUP_COL]])); if (length(grps) < 2) stop("Need at least two groups in metadata")
if (is.null(CASE_LABEL)) CASE_LABEL <- grps[1]
if (is.null(CTRL_LABEL)) CTRL_LABEL <- grps[2]
group_levels <- c(CTRL_LABEL, CASE_LABEL)
group_colors <- setNames(pastel_n(2L), group_levels)
case_rows <- merged[[GROUP_COL]] == CASE_LABEL; ctrl_rows <- merged[[GROUP_COL]] == CTRL_LABEL
if (!any(case_rows)) stop("No samples found for CASE_LABEL: ", CASE_LABEL)
if (!any(ctrl_rows)) stop("No samples found for CTRL_LABEL: ", CTRL_LABEL)
seq_cols <- setdiff(colnames(asv), "sample_id")
case_mat <- data.matrix(merged[case_rows, seq_cols, drop = FALSE]); ctrl_mat <- data.matrix(merged[ctrl_rows, seq_cols, drop = FALSE])

tss_ps <- function(m) { rs <- rowSums(m); out <- sweep(m, 1, ifelse(rs == 0, 1, rs), "/") + DA_PSEUDOCOUNT; out[rs == 0, ] <- DA_PSEUDOCOUNT; out }
case_tss <- tss_ps(case_mat); ctrl_tss <- tss_ps(ctrl_mat)
log2fc <- colMeans(log2(case_tss)) - colMeans(log2(ctrl_tss))
pvals <- vapply(seq_along(seq_cols), function(i) {
  x <- case_mat[, i]; y <- ctrl_mat[, i]
  if (sum(x > 0) + sum(y > 0) < DA_MIN_PREV) return(NA_real_)
  suppressWarnings(wilcox.test(x, y)$p.value)
}, numeric(1))
padj <- p.adjust(pvals, method = "BH")

genus_col <- intersect(c("Genus", "genus"), colnames(tax))[1]; phylum_col <- intersect(c("Phylum", "phylum"), colnames(tax))[1]
da_df <- data.frame(feature_id = seq_cols, log2fc = round(log2fc, 4), pval = pvals, padj = padj, sig = !is.na(padj) & padj < 0.05, stringsAsFactors = FALSE)
if ("feature_id" %in% colnames(tax) && all(da_df$feature_id %in% tax$feature_id)) {
  idx <- match(da_df$feature_id, tax$feature_id)
} else if ("nr_group" %in% colnames(tax) && all(da_df$feature_id %in% tax$nr_group)) {
  idx <- match(da_df$feature_id, tax$nr_group)
} else if (all(da_df$feature_id %in% tax$ASV_ID)) {
  idx <- match(da_df$feature_id, tax$ASV_ID)
} else if (all(da_df$feature_id %in% tax$sequence)) {
  idx <- match(da_df$feature_id, tax$sequence)
} else {
  stop("Filtered feature columns do not match taxonomy feature_id, nr_group, ASV_ID or sequence columns")
}
da_df$sequence <- tax$sequence[idx]
da_df$ASV_ID <- tax$ASV_ID[idx]
if ("feature_id" %in% colnames(tax)) da_df$taxonomy_feature_id <- tax$feature_id[idx]
if ("nr_group" %in% colnames(tax)) da_df$nr_group <- tax$nr_group[idx]
if (!is.na(genus_col)) da_df$Genus <- tax[[genus_col]][idx]
if (!is.na(phylum_col)) da_df$Phylum <- tax[[phylum_col]][idx]
da_df$baseMean_rel <- if (!is.null(tss_for_mean)) {
  tss_for_mean[match(da_df$feature_id, names(tss_for_mean))]
} else 0
da_df$baseMean_rel[is.na(da_df$baseMean_rel)] <- 0
da_df <- da_df[order(da_df$padj, na.last = TRUE), ]

dir.create(file.path(out_dir, "differential"), recursive = TRUE, showWarnings = FALSE)
write.table(da_df, file.path(out_dir, "differential", "da_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
})
da_df$category <- ifelse(da_df$sig & da_df$log2fc > 0, CASE_LABEL,
  ifelse(da_df$sig & da_df$log2fc < 0, CTRL_LABEL, "NS"))
sig_label <- !is.na(da_df$padj) & da_df$padj < 0.05 & abs(da_df$log2fc) >= 1
label_col <- if ("Genus" %in% colnames(da_df)) "Genus" else "feature_id"
cat_colors <- c(setNames(group_colors[CASE_LABEL], CASE_LABEL),
  setNames(group_colors[CTRL_LABEL], CTRL_LABEL),
  NS = "lightgray")
p_volc <- ggplot(da_df, aes(log2fc, -log10(pval),
  color = category, size = baseMean_rel)) +
  geom_point(alpha = 0.9) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = cat_colors, name = "Enriched in") +
  scale_size_continuous(range = c(1.5, 8), name = "Mean relative abundance") +
  theme_minimal(base_size = 12) +
  labs(title = paste0(CASE_LABEL, " vs ", CTRL_LABEL),
    x = "log2(FC)", y = "-log10(p-value)")
if (requireNamespace("ggrepel", quietly = TRUE) && any(sig_label)) {
  p_volc <- p_volc +
    ggrepel::geom_text_repel(data = da_df[sig_label, , drop = FALSE],
      aes(label = .data[[label_col]]),
      size = 2.5, max.overlaps = 20, show.legend = FALSE)
}
pdf(file.path(out_dir, "differential", "volcano.pdf"), width = 12, height = 9)
print(p_volc)
dev.off()

p_ma <- ggplot(da_df, aes(log10(baseMean_rel + 1e-6), log2fc, color = category)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = cat_colors, name = "Enriched in") +
  theme_minimal(base_size = 12) +
  labs(title = paste0("MA plot: ", CASE_LABEL, " vs ", CTRL_LABEL),
    x = "log10(mean relative abundance + 1e-6)", y = "log2(FC)")
pdf(file.path(out_dir, "differential", "ma.pdf"), width = 8, height = 5)
print(p_ma)
dev.off()

dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
n_sig <- sum(da_df$sig, na.rm = TRUE)
flags <- data.frame(step = "10_differential", key = c("case_label", "ctrl_label", "n_tested", "n_sig_asvs"),
  value = c(CASE_LABEL, CTRL_LABEL, sum(!is.na(da_df$pval)), n_sig), status = c("ok", "ok", "ok", if (n_sig == 0) "warn" else "ok"), stringsAsFactors = FALSE)
write.table(flags, file.path(out_dir, "flags", "flags_10_differential.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("DA: case=", CASE_LABEL, " ctrl=", CTRL_LABEL, " | significant features: ", n_sig)
