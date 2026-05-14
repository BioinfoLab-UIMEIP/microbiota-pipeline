# =============================================================================
# Script:      07_alpha.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Alpha diversity metrics (richness, evenness, Shannon, Simpson) and group comparison
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 07_alpha.R <config.R> --outdir <dir>")
source(cfg)
if (!exists("GROUP_COL")) GROUP_COL <- "Group"
if (!exists("SAMPLE_ID_COL")) SAMPLE_ID_COL <- "Sample_ID"
if (!exists("CASE_LABEL")) CASE_LABEL <- NULL
if (!exists("CTRL_LABEL")) CTRL_LABEL <- NULL
if (!exists("METADATA_FILE") || !nzchar(METADATA_FILE)) stop("METADATA_FILE must be set in config")

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
asv_path <- file.path(out_dir, "filter", "asv_table_filtered.tsv"); meta_path <- path.expand(METADATA_FILE)
if (!file.exists(asv_path)) stop("asv_table_filtered.tsv not found: ", asv_path)
if (!file.exists(meta_path)) stop("Metadata file not found: ", meta_path)

message("Loading filtered ASV table...")
asv <- read.table(asv_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample_id" %in% colnames(asv)) stop("asv_table_filtered.tsv must contain sample_id column")
message("Loading metadata...")
sep <- if (grepl("\\.csv$", meta_path, ignore.case = TRUE)) "," else "\t"
meta <- read.table(meta_path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!SAMPLE_ID_COL %in% colnames(meta)) stop("Metadata sample ID column not found: ", SAMPLE_ID_COL)
if (!exists("SAMPLE_TYPE_COL")) SAMPLE_TYPE_COL <- "sample_type"
if (exists("SAMPLE_TYPE_FILTER") && !is.null(SAMPLE_TYPE_FILTER) &&
    length(SAMPLE_TYPE_FILTER) > 0 && nzchar(SAMPLE_TYPE_FILTER[1]) &&
    SAMPLE_TYPE_COL %in% colnames(meta)) {
  meta <- meta[meta[[SAMPLE_TYPE_COL]] %in% SAMPLE_TYPE_FILTER, , drop = FALSE]
  message("Sample type filter applied: ",
          paste(SAMPLE_TYPE_FILTER, collapse = ", "),
          " (", nrow(meta), " samples retained)")
}

mat <- as.matrix(asv[, -1, drop = FALSE]); rownames(mat) <- asv$sample_id
obs <- rowSums(mat > 0)
shannon <- apply(mat, 1, function(x) { x <- x[x > 0]; if (!length(x)) return(0); p <- x / sum(x); -sum(p * log(p)) })
simpson <- apply(mat, 1, function(x) { x <- x[x > 0]; if (!length(x)) return(0); p <- x / sum(x); 1 - sum(p^2) })
evenness <- ifelse(obs <= 1L, 0, shannon / log(obs))
alpha_df <- data.frame(sample_id = rownames(mat),
  observed = obs, evenness = round(evenness, 4),
  shannon = round(shannon, 4), simpson = round(simpson, 4),
  stringsAsFactors = FALSE)

message("Merging alpha diversity with metadata...")
merged <- merge(alpha_df, meta, by.x = "sample_id", by.y = SAMPLE_ID_COL, all.x = TRUE)
group_ok <- GROUP_COL %in% colnames(merged)
if (!group_ok) { message("WARNING: GROUP_COL '", GROUP_COL, "' not found in metadata"); GROUP_COL <- NULL }

dir.create(file.path(out_dir, "alpha"), recursive = TRUE, showWarnings = FALSE)
write.table(merged, file.path(out_dir, "alpha", "alpha_diversity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(GROUP_COL)) {
  suppressPackageStartupMessages(library(ggplot2))
  has_gridextra <- requireNamespace("gridExtra", quietly = TRUE)
  if (has_gridextra) suppressPackageStartupMessages(library(gridExtra))
  merged <- merged[!is.na(merged[[GROUP_COL]]) & nzchar(trimws(as.character(merged[[GROUP_COL]]))), , drop = FALSE]
  if (!is.null(CASE_LABEL) && !is.null(CTRL_LABEL)) {
    merged <- merged[merged[[GROUP_COL]] %in% c(CASE_LABEL, CTRL_LABEL), , drop = FALSE]
    merged[[GROUP_COL]] <- factor(merged[[GROUP_COL]], levels = c(CTRL_LABEL, CASE_LABEL))
  }
  metrics <- c("observed", "evenness", "shannon", "simpson")
  labels <- c("Richness", "Evenness", "Shannon", "Simpson")
  metric_p <- setNames(rep(NA_real_, length(metrics)), metrics)
  for (m in metrics) {
    keep <- !is.na(merged[[m]]) & !is.na(merged[[GROUP_COL]])
    metric_p[m] <- if (sum(keep) > 1L && length(unique(merged[[GROUP_COL]][keep])) == 2L) {
      suppressWarnings(wilcox.test(merged[[m]][keep] ~ merged[[GROUP_COL]][keep])$p.value)
    } else {
      NA_real_
    }
  }
  metric_padj <- p.adjust(metric_p, method = "BH")
  metric_padj[is.na(metric_p)] <- NA_real_
  alpha_plots <- vector("list", length(metrics))
  for (i in seq_along(metrics)) {
    m <- metrics[i]
    keep <- !is.na(merged[[m]]) & !is.na(merged[[GROUP_COL]])
    plot_df <- merged[keep, , drop = FALSE]
    alpha_plots[[i]] <- ggplot(plot_df, aes(x = .data[[GROUP_COL]], y = .data[[m]], fill = .data[[GROUP_COL]], color = .data[[GROUP_COL]])) +
      geom_violin(trim = FALSE, alpha = 0.45) +
      geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8, color = "black") +
      geom_jitter(width = 0.12, size = 1, alpha = 0.5, show.legend = FALSE) +
      scale_fill_manual(values = group_colors) +
      scale_color_manual(values = group_colors) +
      labs(
        title = labels[i],
        subtitle = paste0("p adj = ", signif(metric_padj[m], 3)),
        x = "",
        y = labels[i]
      ) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")
  }
  if (has_gridextra) {
    pdf(file.path(out_dir, "alpha", "alpha_diversity.pdf"), width = 15, height = 5)
    do.call(gridExtra::grid.arrange, c(alpha_plots, nrow = 1))
    dev.off()
  } else {
    pdf(file.path(out_dir, "alpha", "alpha_diversity.pdf"), width = 5, height = 4)
    for (p in alpha_plots) print(p)
    dev.off()
  }
  exclude_cols <- unique(c("sample_id", GROUP_COL, SAMPLE_ID_COL))
  numeric_covars <- character(0)
  for (col in setdiff(colnames(merged), exclude_cols)) {
    vals <- suppressWarnings(as.numeric(as.character(merged[[col]])))
    if (is.numeric(vals) && any(!is.na(vals)) && mean(is.na(vals)) < 0.5) numeric_covars <- c(numeric_covars, col)
  }
  if (!length(numeric_covars)) {
    message("Figure alpha_covariates skipped: no numeric covariates")
  } else {
    cov_plots <- list()
    cov_metrics <- c("shannon", "simpson")
    for (m in cov_metrics) {
      for (covar in numeric_covars) {
        plot_df <- merged[, c(GROUP_COL, m, covar), drop = FALSE]
        plot_df[[covar]] <- suppressWarnings(as.numeric(as.character(plot_df[[covar]])))
        plot_df <- plot_df[!is.na(plot_df[[m]]) & !is.na(plot_df[[covar]]) & !is.na(plot_df[[GROUP_COL]]), , drop = FALSE]
        cov_plots[[paste(m, covar, sep = "_")]] <- ggplot(plot_df, aes(x = .data[[covar]], y = .data[[m]])) +
          geom_point(aes(color = .data[[GROUP_COL]]), alpha = 0.7, size = 1.5) +
          geom_smooth(aes(color = .data[[GROUP_COL]]), method = "lm", se = TRUE, formula = y ~ x, alpha = 0.2) +
          scale_color_manual(values = group_colors) +
          labs(title = paste(labels[match(m, metrics)], "vs", covar), x = covar, y = labels[match(m, metrics)], color = GROUP_COL) +
          theme_minimal(base_size = 10)
      }
    }
    if (has_gridextra) {
      pdf(file.path(out_dir, "alpha", "alpha_covariates.pdf"), width = max(18, 4 * length(numeric_covars)), height = 10)
      do.call(gridExtra::grid.arrange, c(cov_plots, nrow = length(cov_metrics), ncol = length(numeric_covars)))
      dev.off()
    } else {
      pdf(file.path(out_dir, "alpha", "alpha_covariates.pdf"), width = 6, height = 4)
      for (p in cov_plots) print(p)
      dev.off()
    }
  }
  comparison_label <- if (!is.null(CASE_LABEL) && !is.null(CTRL_LABEL)) {
    paste0(CASE_LABEL, "_vs_", CTRL_LABEL)
  } else {
    grp_levels <- unique(as.character(merged[[GROUP_COL]]))
    if (length(grp_levels) >= 2L) paste0(grp_levels[1], "_vs_", grp_levels[2]) else NA_character_
  }
  wres_df <- data.frame(
    metric = metrics,
    comparison = comparison_label,
    p_adj = unname(metric_padj[metrics]),
    stringsAsFactors = FALSE
  )
  wres_df <- wres_df[!is.na(wres_df$p_adj), , drop = FALSE]
  write.table(wres_df, file.path(out_dir, "alpha", "alpha_wilcoxon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

flags <- data.frame(
  step = "07_alpha",
  key = c("n_samples", "group_col", "median_shannon", "median_observed"),
  value = c(nrow(merged), if (is.null(GROUP_COL)) "missing" else GROUP_COL, round(stats::median(merged$shannon, na.rm = TRUE), 4), stats::median(merged$observed, na.rm = TRUE)),
  status = c("ok", if (is.null(GROUP_COL)) "warn" else "ok", "ok", "ok"),
  stringsAsFactors = FALSE
)
dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
write.table(flags, file.path(out_dir, "flags", "flags_07_alpha.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Alpha diversity done: ", nrow(merged), " samples")
