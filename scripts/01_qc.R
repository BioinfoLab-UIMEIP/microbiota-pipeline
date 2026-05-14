# =============================================================================
# Script:      01_qc.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Quality profile visualization and TRUNC_LEN suggestion from raw FASTQs
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- if (length(args) >= 1) args[1] else stop("Usage: Rscript scripts/01_qc.R <config.R> --outdir <dir>")
source(cfg)

if (!exists("MANIFEST_FILE"))     MANIFEST_FILE     <- ""
if (!exists("QC_N_SAMPLES"))      QC_N_SAMPLES      <- 8L
if (!exists("Q_THRESHOLD"))       Q_THRESHOLD       <- 25L
if (!exists("QC_READS_PER_FILE")) QC_READS_PER_FILE <- 5000L
if (!exists("SAMPLE_TYPE_FILTER")) SAMPLE_TYPE_FILTER <- NULL
if (!exists("SAMPLE_TYPE_COL"))    SAMPLE_TYPE_COL    <- "sample_type"
if (!exists("SAMPLE_ID_COL"))      SAMPLE_ID_COL      <- "Sample_ID"

suppressPackageStartupMessages({
  library(ShortRead)
  library(ggplot2)
  library(gridExtra)
})

od_idx  <- which(args == "--outdir")
out_dir <- if (length(od_idx)) path.expand(args[od_idx[1] + 1]) else stop("--outdir <dir> is required")
dir.create(file.path(out_dir, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
message("Preparing manifest")
if (nzchar(MANIFEST_FILE)) {
  manifest <- read.delim(MANIFEST_FILE, sep = "\t", stringsAsFactors = FALSE)
} else {
  r1 <- list.files(RAW_DIR, pattern = "_L001_R1_001\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE)
  r1 <- r1[file.info(r1)$size >= 1e6]
  r2 <- sub("_R1_", "_R2_", r1)
  keep <- file.exists(r2)
  manifest <- data.frame(
    sample_id = gsub("[^A-Za-z0-9]+", "_", sub("_L001_R1_001\\.fastq\\.gz$", "", basename(r1[keep]))),
    read1 = r1[keep], read2 = r2[keep], stringsAsFactors = FALSE
  )
}
if (!nrow(manifest)) stop("No valid FASTQ pairs found")

# Apply SAMPLE_TYPE_FILTER: join with metadata and keep only matching sample types.
# This filters the manifest before any processing so all downstream steps see only
# the requested sample type subset.
n_total <- nrow(manifest)
if (!is.null(SAMPLE_TYPE_FILTER) && length(SAMPLE_TYPE_FILTER) > 0 &&
    nzchar(SAMPLE_TYPE_FILTER[1])) {
  if (!nzchar(METADATA_FILE) || !file.exists(METADATA_FILE))
    stop("METADATA_FILE must be set and exist to apply SAMPLE_TYPE_FILTER")
  meta     <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)
  if (!SAMPLE_ID_COL   %in% colnames(meta)) stop("SAMPLE_ID_COL '",   SAMPLE_ID_COL,   "' not found in metadata")
  if (!SAMPLE_TYPE_COL %in% colnames(meta)) stop("SAMPLE_TYPE_COL '", SAMPLE_TYPE_COL, "' not found in metadata")
  # Normalise metadata IDs the same way manifest sample_id is derived from filenames
  keep_ids <- gsub("[^A-Za-z0-9]+", "_", meta[[SAMPLE_ID_COL]][meta[[SAMPLE_TYPE_COL]] %in% SAMPLE_TYPE_FILTER])
  manifest <- manifest[manifest$sample_id %in% keep_ids, , drop = FALSE]
  if (!nrow(manifest))
    stop("No samples remain after SAMPLE_TYPE_FILTER='",
         paste(SAMPLE_TYPE_FILTER, collapse = ", "),
         "'. Check that SAMPLE_ID_COL and SAMPLE_TYPE_COL are correct.")
  message("SAMPLE_TYPE_FILTER='", paste(SAMPLE_TYPE_FILTER, collapse = ", "),
          "': ", nrow(manifest), " / ", n_total, " samples retained")
}

write.table(manifest, file.path(out_dir, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

profile_fastq <- function(path) {
  sampler <- FastqSampler(path, n = QC_READS_PER_FILE)
  on.exit(close(sampler), add = TRUE)
  fq <- yield(sampler)
  qmat <- as(quality(fq), "matrix")
  list(
    qual = data.frame(
      cycle    = seq_len(ncol(qmat)),
      median_q = apply(qmat, 2, median, na.rm = TRUE),
      q25      = apply(qmat, 2, quantile, probs = 0.25, na.rm = TRUE),
      q75      = apply(qmat, 2, quantile, probs = 0.75, na.rm = TRUE)
    ),
    lengths = as.integer(width(sread(fq)))
  )
}

aggregate_profiles <- function(paths) {
  lst <- lapply(paths, profile_fastq)
  max_cycle <- max(vapply(lst, function(x) nrow(x$qual), integer(1)))
  agg <- do.call(rbind, lapply(seq_len(max_cycle), function(i) {
    rows <- do.call(rbind, lapply(lst, function(x) {
      if (nrow(x$qual) >= i) x$qual[i, c("median_q", "q25", "q75")] else NULL
    }))
    data.frame(cycle = i, median_q = median(rows$median_q), q25 = median(rows$q25), q75 = median(rows$q75))
  }))
  list(profile = agg, lengths = unlist(lapply(lst, `[[`, "lengths"), use.names = FALSE))
}

set.seed(42)
sel <- sample(seq_len(nrow(manifest)), min(QC_N_SAMPLES, nrow(manifest)))
message("Profiling ", length(sel), " samples")
r1_res <- aggregate_profiles(manifest$read1[sel])
r2_res <- aggregate_profiles(manifest$read2[sel])
get_trunc <- function(profile, lengths) {
  ok <- which(profile$median_q >= Q_THRESHOLD)
  last_ok_cycle <- if (length(ok)) max(ok) else 0L
  p05_len <- floor(as.numeric(quantile(lengths, 0.05, na.rm = TRUE)))
  trunc <- min(last_ok_cycle, p05_len)
  if (trunc < 100L) 0L else as.integer(trunc)
}
trunc_r1 <- get_trunc(r1_res$profile, r1_res$lengths)
trunc_r2 <- get_trunc(r2_res$profile, r2_res$lengths)

plot_profile <- function(df, trunc, label) {
  max_cycle <- max(df$cycle)
  p <- ggplot(df, aes(cycle, median_q)) +
    # quality zone background
    annotate("rect", xmin = -Inf, xmax = Inf, ymin =  0, ymax = 20, fill = "#f8d7da", alpha = 0.5) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20, ymax = 28, fill = "#fff3cd", alpha = 0.5) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 28, ymax = 40, fill = "#d4edda", alpha = 0.5) +
    # zone labels (right margin)
    annotate("text", x = max_cycle, y = 10, label = "Low (<Q20)",  hjust = 1, size = 2.8, color = "#721c24") +
    annotate("text", x = max_cycle, y = 24, label = "Medium",       hjust = 1, size = 2.8, color = "#856404") +
    annotate("text", x = max_cycle, y = 38, label = "High (≥Q28)", hjust = 1, size = 2.8, color = "#155724") +
    # IQR ribbon + median line
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey50", alpha = 0.35) +
    geom_line(color = "steelblue", linewidth = 0.9) +
    # threshold line + label
    geom_hline(yintercept = Q_THRESHOLD, linetype = "dashed", color = "red", linewidth = 0.7) +
    annotate("text", x = 1, y = Q_THRESHOLD + 1, label = paste0("Threshold Q", Q_THRESHOLD),
             hjust = 0, size = 2.8, color = "red") +
    scale_y_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 28, 30, 40)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(title = label, x = "Ciclo", y = "Calidad (Phred)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  # truncation line + label
  if (trunc > 0)
    p <- p +
      geom_vline(xintercept = trunc, linetype = "dashed", color = "darkgreen", linewidth = 0.9) +
      annotate("text", x = trunc + 1, y = 2, label = paste0("Trunc=", trunc),
               hjust = 0, size = 2.8, color = "darkgreen")
  p
}

plot_hist <- function(lengths, label) {
  p05 <- floor(quantile(lengths, 0.05))
  ggplot(data.frame(length = lengths), aes(length)) +
    geom_histogram(binwidth = 2, fill = "steelblue", color = "white") +
    geom_vline(xintercept = p05, linetype = "dashed", color = "darkorange", linewidth = 0.8) +
    annotate("text", x = p05 + 1, y = Inf, label = paste0("p05=", p05),
             hjust = 0, vjust = 1.5, size = 2.8, color = "darkorange") +
    labs(title = label, x = "Longitud (bp)", y = "Lecturas") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

n_label <- paste0("(", length(sel), " muestras, ", QC_READS_PER_FILE, " lecturas c/u)")
pdf(file.path(out_dir, "qc", "qc_report.pdf"), width = 13, height = 9)
grid.arrange(
  plot_profile(r1_res$profile, trunc_r1, paste("R1 — perfil de calidad", n_label)),
  plot_hist(r1_res$lengths,              paste("R1 — distribución de longitudes", n_label)),
  plot_profile(r2_res$profile, trunc_r2, paste("R2 — perfil de calidad", n_label)),
  plot_hist(r2_res$lengths,              paste("R2 — distribución de longitudes", n_label)),
  ncol = 2, nrow = 2
)
dev.off()

writeLines(c(sprintf("TRUNC_LEN_R1=%d", trunc_r1), sprintf("TRUNC_LEN_R2=%d", trunc_r2)),
           file.path(out_dir, "qc", "trunc_suggestion.txt"))
sample_type_label <- if (!is.null(SAMPLE_TYPE_FILTER) && nzchar(SAMPLE_TYPE_FILTER[1]))
  paste(SAMPLE_TYPE_FILTER, collapse = ", ") else "all"
flags <- data.frame(
  step = "01_qc",
  key = c("n_samples_raw", "n_samples_after_filter", "sample_type_filter",
          "n_qc", "trunc_r1", "trunc_r2", "out_dir", "region", "kit"),
  value = c(n_total, nrow(manifest), sample_type_label,
            length(sel), trunc_r1, trunc_r2, out_dir, REGION, KIT),
  status = c("ok", "ok", "ok",
             "ok", if (trunc_r1 > 0) "ok" else "warn", if (trunc_r2 > 0) "ok" else "warn",
             "ok", "ok", "ok"),
  stringsAsFactors = FALSE
)
write.table(flags, file.path(out_dir, "flags", "flags_01_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("QC done \u2014 ", out_dir, " | TRUNC R1=", trunc_r1, " R2=", trunc_r2)
