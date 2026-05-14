# =============================================================================
# Script:      06_filter.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Length filter, non-bacterial removal, feature-level collapse, prevalence filter
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 06_filter.R <config.R> --outdir <dir>")
source(cfg)

if (!exists("ASV_LENGTH_MIN")) ASV_LENGTH_MIN <- NA_integer_
if (!exists("ASV_LENGTH_MAX")) ASV_LENGTH_MAX <- NA_integer_
if (!exists("MIN_PREVALENCE")) MIN_PREVALENCE <- 2L
if (!exists("USE_COLLAPSED_ASVS")) USE_COLLAPSED_ASVS <- FALSE
if (!exists("FEATURE_LEVEL")) FEATURE_LEVEL <- "asv"
if (!exists("MIN_TOTAL_ABUNDANCE_UNIQUE")) MIN_TOTAL_ABUNDANCE_UNIQUE <- 10L
if (!exists("KEEP_KINGDOMS")) KEEP_KINGDOMS <- c("Bacteria", "Archaea")
if (!exists("KEEP_UNASSIGNED_KINGDOM")) KEEP_UNASSIGNED_KINGDOM <- TRUE

od_idx  <- which(args == "--outdir")
out_dir <- if (length(od_idx)) path.expand(args[od_idx[1] + 1]) else stop("--outdir <dir> is required")
asv_path <- if (USE_COLLAPSED_ASVS) file.path(out_dir, "taxonomy", "redundancy", "asv_table_collapsed.tsv") else file.path(out_dir, "dada2", "asv_table.tsv")
tax_path <- if (USE_COLLAPSED_ASVS) file.path(out_dir, "taxonomy", "taxonomy_nrASVs.tsv") else file.path(out_dir, "taxonomy", "taxonomy_ASVs.tsv")
if (!file.exists(asv_path)) stop("ASV table not found: ", asv_path)
if (!file.exists(tax_path)) stop("taxonomy table not found: ", tax_path)
asv <- read.table(asv_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tax <- read.table(tax_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample_id") %in% colnames(asv))) stop("asv_table.tsv must contain sample_id column")
if (!all(c("ASV_ID", "sequence") %in% colnames(tax))) stop("taxonomy.tsv must contain ASV_ID and sequence columns")

input_cols <- setdiff(colnames(asv), "sample_id")
if (!exists("KEEP_KINGDOMS")) KEEP_KINGDOMS <- c("Bacteria", "Archaea")
if (!exists("KEEP_UNASSIGNED_KINGDOM")) KEEP_UNASSIGNED_KINGDOM <- TRUE
kingdom_col <- intersect(c("Kingdom", "kingdom"), colnames(tax))[1]
if (!is.na(kingdom_col)) {
  blank_kingdom <- is.na(tax[[kingdom_col]]) | !nzchar(trimws(tax[[kingdom_col]])) |
                   trimws(tax[[kingdom_col]]) %in% c("NA", "N/A", "NULL", "null", "Unassigned")
  keep_rows <- (trimws(tax[[kingdom_col]]) %in% KEEP_KINGDOMS) |
               (KEEP_UNASSIGNED_KINGDOM & blank_kingdom)
  n_removed <- sum(!keep_rows)
  if (n_removed > 0) {
    removed_kingdoms <- table(tax[[kingdom_col]][!keep_rows])
    message("Kingdom filter: removed ", n_removed, " ASVs — ",
            paste(names(removed_kingdoms), removed_kingdoms, sep = "=", collapse = ", "))
  }
  tax <- tax[keep_rows, , drop = FALSE]
  keep_cols <- c("sample_id",
    intersect(input_cols,
      unique(c(tax$nr_group, tax$ASV_ID, tax$sequence,
               if ("feature_id" %in% colnames(tax)) tax$feature_id))))
  asv <- asv[, intersect(colnames(asv), keep_cols), drop = FALSE]
  input_cols <- setdiff(colnames(asv), "sample_id")
}
if ("nr_group" %in% colnames(tax) && all(input_cols %in% tax$nr_group)) {
  tax_idx <- match(input_cols, tax$nr_group)
} else if (all(input_cols %in% tax$ASV_ID)) {
  tax_idx <- match(input_cols, tax$ASV_ID)
} else if (all(input_cols %in% tax$sequence)) {
  tax_idx <- match(input_cols, tax$sequence)
} else {
  stop("ASV table columns do not match taxonomy ASV_ID or sequence columns")
}
asv_lens <- nchar(tax$sequence[tax_idx])
db_path <- file.path(dirname(dirname(normalizePath(cfg))), "scripts", "primers.tsv")
if (file.exists(db_path)) {
  db <- read.table(db_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  row <- db[db$kit == KIT & db$region == REGION, , drop = FALSE]
  if (nrow(row) == 1) {
    if (is.na(ASV_LENGTH_MIN)) ASV_LENGTH_MIN <- as.integer(row$amplicon_min_bp)
    if (is.na(ASV_LENGTH_MAX)) ASV_LENGTH_MAX <- as.integer(row$amplicon_max_bp)
  }
}
if (is.na(ASV_LENGTH_MIN)) ASV_LENGTH_MIN <- 0L
if (is.na(ASV_LENGTH_MAX)) ASV_LENGTH_MAX <- .Machine$integer.max

keep_seq_cols <- input_cols[asv_lens >= ASV_LENGTH_MIN & asv_lens <= ASV_LENGTH_MAX]
asv_len <- asv[, c("sample_id", keep_seq_cols), drop = FALSE]
tax_len <- tax[tax_idx[match(keep_seq_cols, input_cols)], , drop = FALSE]

if (identical(tolower(FEATURE_LEVEL), "species")) {
  if (!"Species" %in% colnames(tax_len)) stop("taxonomy table must contain Species column for FEATURE_LEVEL='species'")
  species_ids <- trimws(as.character(tax_len$Species))
  species_ids[!nzchar(species_ids) | species_ids %in% c("NA", "Na", "N/A", "NULL")] <- paste0("nr_", seq_len(sum(!nzchar(species_ids) | species_ids %in% c("NA", "Na", "N/A", "NULL"))))
  mat_len <- data.matrix(asv_len[, -1, drop = FALSE])
  colnames(mat_len) <- species_ids
  agg_mat <- t(rowsum(t(mat_len), group = species_ids, reorder = FALSE))
  prevalence <- if (ncol(agg_mat)) colSums(agg_mat > 0) else numeric(0)
  total_abundance <- if (ncol(agg_mat)) colSums(agg_mat) else numeric(0)
  keep_prev <- names(prevalence)[prevalence >= MIN_PREVALENCE | total_abundance >= MIN_TOTAL_ABUNDANCE_UNIQUE]
  asv_filt <- data.frame(sample_id = asv_len$sample_id, agg_mat[, keep_prev, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)

  tax_len$feature_id <- species_ids
  if (!"centroid_total_abundance" %in% colnames(tax_len)) tax_len$centroid_total_abundance <- 0
  split_tax <- split(tax_len, tax_len$feature_id)
  tax_filt <- do.call(rbind, lapply(split_tax, function(df) {
    df[which.max(df$centroid_total_abundance), , drop = FALSE]
  }))
  tax_filt <- tax_filt[match(keep_prev, tax_filt$feature_id), , drop = FALSE]
  tax_filt$member_n <- as.integer(table(species_ids)[tax_filt$feature_id])
  rownames(tax_filt) <- NULL
} else {
  prevalence <- if (ncol(asv_len) > 1) colSums(asv_len[, -1, drop = FALSE] > 0) else numeric(0)
  keep_prev <- names(prevalence)[prevalence >= MIN_PREVALENCE]
  asv_filt <- asv_len[, c("sample_id", keep_prev), drop = FALSE]
  if ("nr_group" %in% colnames(tax) && all(keep_prev %in% tax$nr_group)) {
    tax_filt <- tax[tax$nr_group %in% keep_prev, , drop = FALSE]
    tax_filt <- tax_filt[match(keep_prev, tax_filt$nr_group), , drop = FALSE]
  } else if (all(keep_prev %in% tax$ASV_ID)) {
    tax_filt <- tax[tax$ASV_ID %in% keep_prev, , drop = FALSE]
    tax_filt <- tax_filt[match(keep_prev, tax_filt$ASV_ID), , drop = FALSE]
  } else {
    tax_filt <- tax[tax$sequence %in% keep_prev, , drop = FALSE]
    tax_filt <- tax_filt[match(keep_prev, tax_filt$sequence), , drop = FALSE]
  }
  tax_filt$feature_id <- keep_prev
}

mat <- data.matrix(asv_filt[, -1, drop = FALSE])
row_sums <- rowSums(mat)
tss_mat <- if (ncol(mat)) sweep(mat, 1, row_sums, "/") else mat
if (nrow(tss_mat)) tss_mat[row_sums == 0, ] <- 0
asv_tss <- data.frame(sample_id = asv_filt$sample_id, tss_mat, check.names = FALSE)

dir.create(file.path(out_dir, "filter"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "taxonomy"), recursive = TRUE, showWarnings = FALSE)
write.table(asv_filt, file.path(out_dir, "filter", "asv_table_filtered.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(asv_tss, file.path(out_dir, "filter", "asv_table_tss.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tax_filt, file.path(out_dir, "filter", "taxonomy_filtered.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tax_filt, file.path(out_dir, "taxonomy", "taxonomy_assigned_collapsed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

n_asvs_raw <- length(input_cols); n_asvs_len <- length(keep_seq_cols); n_asvs_prev <- ncol(asv_filt) - 1L
flags <- data.frame(
  step = "06_filter",
  key = c("n_asvs_raw", "n_asvs_len_pass", "n_features_pass", "length_min", "length_max", "min_prevalence", "min_total_abundance_unique", "n_samples", "input_mode", "feature_level"),
  value = c(n_asvs_raw, n_asvs_len, n_asvs_prev, ASV_LENGTH_MIN, ASV_LENGTH_MAX, MIN_PREVALENCE, MIN_TOTAL_ABUNDANCE_UNIQUE, nrow(asv_filt), if (USE_COLLAPSED_ASVS) "collapsed" else "original", FEATURE_LEVEL),
  status = c("ok", if (n_asvs_len < n_asvs_raw * 0.5) "warn" else "ok", if (n_asvs_prev < n_asvs_raw * 0.3) "warn" else "ok", "ok", "ok", "ok", "ok", "ok", "ok", "ok"),
  stringsAsFactors = FALSE
)
dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
write.table(flags, file.path(out_dir, "flags", "flags_06_filter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Filter: ", n_asvs_raw, " ASVs -> ", n_asvs_len, " (length) -> ", n_asvs_prev, " ", FEATURE_LEVEL, " features")
