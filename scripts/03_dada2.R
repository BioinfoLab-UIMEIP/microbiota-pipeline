# =============================================================================
# Script:      03_dada2.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: DADA2 denoising, merging and chimera removal; produces ASV table
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg  <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 03_dada2.R <config.R> --outdir <dir>")
source(cfg)

if (!exists("TRUNC_LEN"))    TRUNC_LEN    <- c(0L, 0L)
if (!exists("TRIM_LEFT"))    TRIM_LEFT    <- c(0L, 0L)
if (!exists("TRUNC_Q"))      TRUNC_Q      <- 2L
if (!exists("MAX_EE"))       MAX_EE       <- c(2, 2)
if (!exists("MAX_N"))        MAX_N        <- 0L
if (!exists("RM_PHIX"))      RM_PHIX      <- TRUE
if (!exists("MIN_OVERLAP"))  MIN_OVERLAP  <- 12L
if (!exists("MAX_MISMATCH")) MAX_MISMATCH <- 0L
if (!exists("POOL_METHOD"))  POOL_METHOD  <- "pseudo"

suppressPackageStartupMessages({ library(dada2) })

od_idx  <- which(args == "--outdir")
out_dir <- if (length(od_idx)) path.expand(args[od_idx[1] + 1]) else stop("--outdir <dir> is required")
if (all(TRUNC_LEN == 0L)) {
  trunc_file <- file.path(out_dir, "qc", "trunc_suggestion.txt")
  if (file.exists(trunc_file)) {
    lines <- readLines(trunc_file)
    TRUNC_LEN[1] <- as.integer(sub("TRUNC_LEN_R1=", "", lines[grep("TRUNC_LEN_R1", lines)]))
    TRUNC_LEN[2] <- as.integer(sub("TRUNC_LEN_R2=", "", lines[grep("TRUNC_LEN_R2", lines)]))
    message("TRUNC_LEN from QC: R1=", TRUNC_LEN[1], " R2=", TRUNC_LEN[2])
  } else message("WARNING: trunc_suggestion.txt not found, using TRUNC_LEN=c(0,0) (no truncation)")
}

message("Scanning demuxed FASTQs")
demux_dir <- file.path(out_dir, "demux")
fwd <- sort(list.files(demux_dir, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE))
rev <- sub("_R1\\.fastq\\.gz$", "_R2.fastq.gz", fwd)
sid <- sub("_R1\\.fastq\\.gz$", "", basename(fwd))
ok_files <- file.exists(rev) & file.info(fwd)$size > 0 & file.info(rev)$size > 0
if (any(!ok_files)) message("Skipping ", sum(!ok_files), " samples with missing/empty mates")
fwd <- fwd[ok_files]; rev <- rev[ok_files]; sid <- sid[ok_files]
if (!length(fwd)) stop("No valid demuxed paired FASTQs found in ", demux_dir)
names(fwd) <- names(rev) <- sid

dada2_dir <- file.path(out_dir, "dada2")
filt_dir <- file.path(dada2_dir, "filtered")
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
fwd_filt <- setNames(file.path(filt_dir, paste0(sid, "_R1.fastq.gz")), sid)
rev_filt <- setNames(file.path(filt_dir, paste0(sid, "_R2.fastq.gz")), sid)

message("Running filterAndTrim on ", length(sid), " samples")
filt <- filterAndTrim(fwd, fwd_filt, rev, rev_filt, truncLen = TRUNC_LEN, trimLeft = TRIM_LEFT,
                      truncQ = TRUNC_Q, maxEE = MAX_EE, maxN = MAX_N, rm.phix = RM_PHIX,
                      compress = TRUE, multithread = TRUE)
keep <- filt[, "reads.out"] > 0
message("Samples passing filterAndTrim: ", sum(keep), "/", nrow(filt))
if (!any(keep)) stop("No samples passed filterAndTrim")
fwd_filt_ok <- fwd_filt[keep]; rev_filt_ok <- rev_filt[keep]

message("Learning error models")
err_fwd <- learnErrors(fwd_filt_ok, multithread = TRUE)
err_rev <- learnErrors(rev_filt_ok, multithread = TRUE)
pdf(file.path(dada2_dir, "error_model_R1.pdf")); plotErrors(err_fwd, nominalQ = TRUE); dev.off()
pdf(file.path(dada2_dir, "error_model_R2.pdf")); plotErrors(err_rev, nominalQ = TRUE); dev.off()

message("Running dada")
dada_fwd <- dada(fwd_filt_ok, err = err_fwd, pool = POOL_METHOD, multithread = TRUE)
dada_rev <- dada(rev_filt_ok, err = err_rev, pool = POOL_METHOD, multithread = TRUE)
message("Merging pairs")
mergers <- mergePairs(dada_fwd, fwd_filt_ok, dada_rev, rev_filt_ok, minOverlap = MIN_OVERLAP,
                      maxMismatch = MAX_MISMATCH)
if (is.null(names(mergers)) && length(mergers) == length(fwd_filt_ok)) {
  names(mergers) <- names(fwd_filt_ok)
}
seqtab <- makeSequenceTable(mergers)
if (is.null(rownames(seqtab)) && nrow(seqtab) == length(mergers)) {
  rownames(seqtab) <- names(mergers)
}
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
if (is.null(rownames(seqtab_nochim)) && nrow(seqtab_nochim) == nrow(seqtab)) {
  rownames(seqtab_nochim) <- rownames(seqtab)
}

all_sid <- rownames(filt)
track <- data.frame(sample_id = all_sid, input = filt[, "reads.in"], filtered = filt[, "reads.out"],
                    denoised_fwd = 0L, denoised_rev = 0L, merged = 0L, nonchim = 0L, stringsAsFactors = FALSE)
getN <- function(x) sum(getUniques(x))
track$denoised_fwd[keep] <- vapply(dada_fwd, getN, integer(1))
track$denoised_rev[keep] <- vapply(dada_rev, getN, integer(1))
track$merged[keep] <- vapply(mergers, getN, integer(1))
nc <- rowSums(seqtab_nochim)
if (is.null(names(nc)) && length(nc) == nrow(seqtab_nochim)) {
  names(nc) <- rownames(seqtab_nochim)
}
nc_match <- match(names(nc), track$sample_id)
track$nonchim[nc_match[!is.na(nc_match)]] <- nc[!is.na(nc_match)]
track$percent_merged <- ifelse(track$input > 0, round(track$nonchim / track$input * 100, 1), NA_real_)
write.table(track, file.path(dada2_dir, "track.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

asv_df <- data.frame(sample_id = rownames(seqtab_nochim), seqtab_nochim, check.names = FALSE)
write.table(asv_df, file.path(dada2_dir, "asv_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
asv_seqs <- colnames(seqtab_nochim); pad <- max(4L, nchar(length(asv_seqs)))
hdr <- sprintf(">ASV_%0*d", pad, seq_along(asv_seqs))
writeLines(as.vector(rbind(hdr, asv_seqs)), file.path(dada2_dir, "rep_seqs.fasta"))

n_input <- length(sid); n_filtered <- sum(keep); n_asvs <- ncol(seqtab); n_asvs_nochim <- ncol(seqtab_nochim)
med_pct <- median(track$percent_merged, na.rm = TRUE)
flags <- data.frame(
  step = "03_dada2",
  key = c("n_input_samples", "n_filtered_samples", "n_asvs", "n_asvs_nochim", "median_pct_merged", "trunc_r1", "trunc_r2"),
  value = c(n_input, n_filtered, n_asvs, n_asvs_nochim, round(med_pct, 1), TRUNC_LEN[1], TRUNC_LEN[2]),
  status = c("ok", if (n_filtered < n_input * 0.8) "warn" else "ok", "ok",
             if (n_asvs_nochim < n_asvs * 0.5) "warn" else "ok", if (is.na(med_pct) || med_pct < 30) "warn" else "ok", "ok", "ok"),
  stringsAsFactors = FALSE
)
write.table(flags, file.path(out_dir, "flags", "flags_03_dada2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("DADA2 done — ", n_filtered, "/", n_input, " filtered | ", n_asvs_nochim, " non-chimeric ASVs")
