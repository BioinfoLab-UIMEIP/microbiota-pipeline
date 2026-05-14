# =============================================================================
# Script:      04_collapse.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: ASV clustering by sequence identity; collapses near-identical ASVs
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg  <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 04_collapse.R <config.R> --outdir <dir>")
source(cfg)

if (!exists("ENABLE_ASV_COLLAPSE")) ENABLE_ASV_COLLAPSE <- FALSE
if (!exists("USE_COLLAPSED_ASVS"))  USE_COLLAPSED_ASVS  <- ENABLE_ASV_COLLAPSE
if (!exists("VSEARCH_BIN"))         VSEARCH_BIN         <- "vsearch"
if (!exists("MAFFT_SMART_BIN"))     MAFFT_SMART_BIN     <- "mafft"
if (!exists("SNP_DISTS_BIN"))       SNP_DISTS_BIN       <- "snp-dists"
if (!exists("VSEARCH_IDENTITY"))    VSEARCH_IDENTITY    <- 0.999
if (!exists("HC_SNP_THRESHOLD"))    HC_SNP_THRESHOLD    <- 1
if (!exists("HC_METHOD"))           HC_METHOD           <- "complete"
if (!exists("MAFFT_MODE"))          MAFFT_MODE          <- "smart"
if (!exists("MAFFT_THREADS"))       MAFFT_THREADS       <- 8L

od_idx  <- which(args == "--outdir")
out_dir <- if (length(od_idx)) path.expand(args[od_idx[1] + 1]) else stop("--outdir <dir> is required")

if (!ENABLE_ASV_COLLAPSE) {
  message("ASV collapse disabled via ENABLE_ASV_COLLAPSE=FALSE; skipping 04_collapse.R")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({ library(Biostrings) })

require_bin <- function(bin_path) {
  resolved <- Sys.which(bin_path)
  if (!nzchar(resolved)) stop("Required binary not found in PATH or config: ", bin_path)
  resolved
}

vsearch_exec <- require_bin(VSEARCH_BIN)
mafft_exec <- require_bin(MAFFT_SMART_BIN)
snp_exec <- require_bin(SNP_DISTS_BIN)

asv_path   <- file.path(out_dir, "dada2", "asv_table.tsv")
fasta_path <- file.path(out_dir, "dada2", "rep_seqs.fasta")
if (!file.exists(asv_path)) stop("asv_table.tsv not found: ", asv_path)
if (!file.exists(fasta_path)) stop("rep_seqs.fasta not found: ", fasta_path)

taxonomy_dir <- file.path(out_dir, "taxonomy")
redundancy_dir <- file.path(taxonomy_dir, "redundancy")
dir.create(redundancy_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)

asv <- read.table(asv_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample_id" %in% colnames(asv)) stop("asv_table.tsv must contain sample_id column")
seqs_dna <- readDNAStringSet(fasta_path)
seqs <- as.character(seqs_dna)
ids <- names(seqs_dna)
seq_cols <- setdiff(colnames(asv), "sample_id")

if (length(seq_cols) != length(ids)) stop("asv_table.tsv column count does not match rep_seqs.fasta sequence count")
if (!identical(unname(seq_cols), unname(seqs))) stop("rep_seqs.fasta sequences do not match asv_table.tsv column names in order")

abund <- asv[, c("sample_id", seq_cols), drop = FALSE]
colnames(abund)[-1] <- ids
mat <- as.matrix(abund[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
total_abund <- colSums(mat)
seq_by_id <- setNames(seqs, ids)

aggregate_by_group <- function(group_ids, mat_in) {
  grp_levels <- unique(group_ids)
  agg <- vapply(grp_levels, function(g) rowSums(mat_in[, group_ids == g, drop = FALSE]), numeric(nrow(mat_in)))
  if (is.null(dim(agg))) agg <- matrix(agg, ncol = 1, dimnames = list(NULL, grp_levels))
  colnames(agg) <- grp_levels
  agg
}

order_ids <- names(sort(total_abund, decreasing = TRUE))
vsearch_in <- file.path(redundancy_dir, "vsearch_input.fasta")
writeXStringSet(DNAStringSet(seq_by_id[order_ids], use.names = TRUE), vsearch_in)

uc_file <- file.path(redundancy_dir, "vsearch_0999.uc")
centroids_fa <- file.path(redundancy_dir, "vsearch_centroids_0999.fasta")
system2(vsearch_exec, c(
  "--cluster_fast", vsearch_in,
  "--id", format(VSEARCH_IDENTITY, scientific = FALSE, trim = TRUE),
  "--strand", "plus",
  "--centroids", centroids_fa,
  "--uc", uc_file
))

uc <- read.delim(uc_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "")
if (ncol(uc) < 10) stop("Unexpected UC format in ", uc_file)
colnames(uc)[1:10] <- c("record_type", "cluster_num", "seq_length", "pct_id", "strand", "unused1", "unused2", "cigar", "query", "target")

v_map <- do.call(rbind, lapply(seq_len(nrow(uc)), function(i) {
  row <- uc[i, ]
  if (row$record_type == "S") {
    data.frame(ASV_ID = row$query, vsearch_centroid_id = row$query, stringsAsFactors = FALSE)
  } else if (row$record_type == "H") {
    data.frame(ASV_ID = row$query, vsearch_centroid_id = row$target, stringsAsFactors = FALSE)
  } else {
    NULL
  }
}))
v_map <- unique(v_map)
if (!nrow(v_map)) stop("No S/H records parsed from ", uc_file)
if (!all(ids %in% v_map$ASV_ID)) stop("Some ASV IDs missing from vsearch mapping")
v_map <- v_map[match(ids, v_map$ASV_ID), , drop = FALSE]

centroid_ids <- unique(v_map$vsearch_centroid_id)
centroid_total <- tapply(total_abund[v_map$ASV_ID], v_map$vsearch_centroid_id, sum)

if (length(centroid_ids) > 1) {
  centroid_seqs <- seq_by_id[centroid_ids]
  hc_input <- file.path(redundancy_dir, "centroids_for_hc.fasta")
  aligned_fa <- file.path(redundancy_dir, "centroids_for_hc.aligned.fasta")
  snp_matrix_file <- file.path(redundancy_dir, "centroids_for_hc.snpdist.tsv")

  writeXStringSet(DNAStringSet(centroid_seqs, use.names = TRUE), hc_input)
  system2(mafft_exec, c("-i", hc_input, "-o", aligned_fa, "-m", MAFFT_MODE, "-t", as.character(MAFFT_THREADS)))
  system2(snp_exec, c("-b", aligned_fa), stdout = snp_matrix_file)

  dist_matrix <- read.table(snp_matrix_file, header = TRUE, row.names = 1, check.names = FALSE)
  hc <- hclust(as.dist(dist_matrix), method = HC_METHOD)
  hc_groups <- cutree(hc, h = HC_SNP_THRESHOLD)

  group_df <- data.frame(vsearch_centroid_id = names(hc_groups), hc_group = as.integer(hc_groups), stringsAsFactors = FALSE)
  write.csv(group_df, file.path(redundancy_dir, "centroids_hc_groups.csv"), row.names = FALSE)
  pdf(file.path(redundancy_dir, "centroids_hc_cutoff.pdf"), width = 13, height = 6)
  plot(hc, labels = FALSE, xlab = "Centroids", ylab = "SNP distance")
  abline(h = HC_SNP_THRESHOLD, col = "blue", lty = 2, lwd = 2)
  dev.off()
} else {
  hc_groups <- setNames(1L, centroid_ids)
}

centroid_group_df <- data.frame(
  vsearch_centroid_id = names(hc_groups),
  hc_group = as.integer(hc_groups),
  centroid_total_abundance = as.numeric(centroid_total[names(hc_groups)]),
  stringsAsFactors = FALSE
)

rep_df <- do.call(rbind, lapply(split(centroid_group_df, centroid_group_df$hc_group), function(df) df[which.max(df$centroid_total_abundance), , drop = FALSE]))
rep_df <- rep_df[order(rep_df$hc_group), , drop = FALSE]
rep_df$final_cluster_id <- sprintf("NR_%04d", rep_df$hc_group)

cluster_key <- rep_df[, c("hc_group", "final_cluster_id", "vsearch_centroid_id"), drop = FALSE]
colnames(cluster_key)[3] <- "final_representative_asv_id"
centroid_group_df <- merge(centroid_group_df, cluster_key, by = "hc_group", all.x = TRUE, sort = FALSE)

map_df <- merge(v_map, centroid_group_df, by = "vsearch_centroid_id", all.x = TRUE, sort = FALSE)
map_df <- map_df[match(ids, map_df$ASV_ID), , drop = FALSE]
map_df$original_sequence <- seq_by_id[map_df$ASV_ID]
map_df$original_total_abundance <- total_abund[map_df$ASV_ID]
map_df$final_representative_sequence <- seq_by_id[map_df$final_representative_asv_id]

final_ids <- rep_df$final_cluster_id
final_group_ids <- map_df$final_cluster_id[match(colnames(mat), map_df$ASV_ID)]
collapsed_mat <- aggregate_by_group(final_group_ids, mat)
collapsed_mat <- collapsed_mat[, final_ids, drop = FALSE]

collapsed_tbl <- data.frame(sample_id = abund$sample_id, collapsed_mat, check.names = FALSE, stringsAsFactors = FALSE)
write.table(collapsed_tbl, file.path(redundancy_dir, "asv_table_collapsed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

final_seq_map <- setNames(rep_df$vsearch_centroid_id, rep_df$final_cluster_id)
final_seqs <- seq_by_id[final_seq_map]
names(final_seqs) <- rep_df$final_cluster_id
writeXStringSet(DNAStringSet(final_seqs, use.names = TRUE), file.path(redundancy_dir, "rep_seqs_collapsed.fasta"))

rep_summary <- data.frame(
  final_cluster_id = rep_df$final_cluster_id,
  hc_group = rep_df$hc_group,
  final_representative_asv_id = rep_df$vsearch_centroid_id,
  final_representative_sequence = seq_by_id[rep_df$vsearch_centroid_id],
  centroid_total_abundance = rep_df$centroid_total_abundance,
  stringsAsFactors = FALSE
)

write.table(map_df, file.path(redundancy_dir, "cluster_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rep_summary, file.path(redundancy_dir, "cluster_representatives.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

stopifnot(sum(mat) == sum(collapsed_mat))

flags <- data.frame(
  step = "04_collapse",
  key = c("n_asvs_original", "n_vsearch_centroids", "n_final_clusters", "vsearch_identity", "hc_snp_threshold", "hc_method", "use_collapsed_asvs"),
  value = c(length(ids), length(centroid_ids), length(final_ids), VSEARCH_IDENTITY, HC_SNP_THRESHOLD, HC_METHOD, USE_COLLAPSED_ASVS),
  status = c("ok", "ok", if (length(final_ids) >= length(ids)) "warn" else "ok", "ok", "ok", "ok", if (USE_COLLAPSED_ASVS) "ok" else "warn"),
  stringsAsFactors = FALSE
)
write.table(flags, file.path(out_dir, "flags", "flags_04_collapse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("ASV collapse done — ", length(ids), " original -> ", length(centroid_ids), " vsearch centroids -> ", length(final_ids), " final clusters")
