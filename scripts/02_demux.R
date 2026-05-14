# =============================================================================
# Script:      02_demux.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Primer trimming and demultiplexing via cutadapt
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg  <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 02_demux.R <config.R> --outdir <dir>")
source(cfg)

if (!exists("CUTADAPT_BIN"))        CUTADAPT_BIN        <- "cutadapt"
if (!exists("CUTADAPT_MIN_LENGTH")) CUTADAPT_MIN_LENGTH <- 100L
if (!exists("CUTADAPT_ERROR_RATE")) CUTADAPT_ERROR_RATE <- 0.1
if (!exists("PRIMER_FWD"))          PRIMER_FWD          <- ""
if (!exists("PRIMER_REV"))          PRIMER_REV          <- ""

od_idx  <- which(args == "--outdir")
out_dir <- if (length(od_idx)) path.expand(args[od_idx[1] + 1]) else stop("--outdir <dir> is required")
message("Results dir: ", out_dir)

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b
# Load primers from catalog — resolve relative to this script, not the config
.script_self <- tryCatch(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  mustWork = FALSE), error = function(e) "")
.scripts_dir <- if (nzchar(.script_self) && !is.na(.script_self) && file.exists(.script_self)) {
  dirname(.script_self)
} else {
  Sys.getenv("SCRIPTS_DIR", dirname(normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE)))
}
db_file <- file.path(.scripts_dir, "primers.tsv")
db      <- read.delim(db_file, stringsAsFactors = FALSE)
pr      <- db[db$kit == KIT & db$region == REGION, ]
if (!nrow(pr)) stop("No entry for kit='", KIT, "' region='", REGION, "' in primers.tsv")

fwd      <- if (nzchar(PRIMER_FWD)) PRIMER_FWD else pr$fwd_seq
rev_p    <- if (nzchar(PRIMER_REV)) PRIMER_REV else pr$rev_seq
adapt_r1 <- if (!is.na(pr$adapter_r1) && nzchar(trimws(pr$adapter_r1))) trimws(pr$adapter_r1) else ""
adapt_r2 <- if (!is.na(pr$adapter_r2) && nzchar(trimws(pr$adapter_r2))) trimws(pr$adapter_r2) else ""

message("Kit: ", KIT, " | Region: ", REGION,
        " | Primers: ", pr$fwd_name, " / ", pr$rev_name)
if (nzchar(adapt_r1)) message("Read-through adapters will also be trimmed")

# Manifest and output dirs
manifest  <- read.delim(file.path(out_dir, "manifest.tsv"), stringsAsFactors = FALSE)
demux_dir <- file.path(out_dir, "demux")
dir.create(demux_dir, recursive = TRUE, showWarnings = FALSE)

# Parse pairs-retained % from cutadapt log
parse_pct <- function(log_lines) {
  target <- grep("Pairs written \\(passing filters\\)", log_lines, value = TRUE)
  if (!length(target)) return(NA_real_)
  m <- regmatches(target[1], regexpr("[0-9]+\\.[0-9]+", target[1]))
  if (!length(m)) return(NA_real_)
  as.numeric(m)
}

# Run cutadapt for one sample, return summary row
trim_sample <- function(i, row) {
  message(sprintf("  [%d/%d] %s", i, nrow(manifest), row$sample_id))
  o1  <- file.path(demux_dir, paste0(row$sample_id, "_R1.fastq.gz"))
  o2  <- file.path(demux_dir, paste0(row$sample_id, "_R2.fastq.gz"))

  cut_args <- c(
    "-g", fwd, "-G", rev_p,
    if (nzchar(adapt_r1)) c("-a", adapt_r1),
    if (nzchar(adapt_r2)) c("-A", adapt_r2),
    "--discard-untrimmed",
    "--minimum-length", CUTADAPT_MIN_LENGTH,
    "-e", CUTADAPT_ERROR_RATE,
    "-o", o1, "-p", o2,
    row$read1, row$read2
  )

  log  <- system2(CUTADAPT_BIN, cut_args, stdout = TRUE, stderr = TRUE)
  ret  <- attr(log, "status"); if (is.null(ret)) ret <- 0L
  pct  <- parse_pct(log)
  # cutadapt 5.x exits with 1 for low-retention warnings (not fatal); 2+ = real error
  ok   <- ret <= 1L && file.exists(o1) && file.size(o1) > 0L

  if (!ok) message("    ERROR: cutadapt failed (code ", ret, ")")

  data.frame(sample_id    = row$sample_id,
             exit_code    = ret,
             ok           = ok,
             pct_retained = pct,
             stringsAsFactors = FALSE)
}

message("Trimming ", nrow(manifest), " samples...")
summary_df <- do.call(rbind, lapply(seq_len(nrow(manifest)),
                                    function(i) trim_sample(i, manifest[i, ])))

write.table(summary_df, file.path(demux_dir, "trim_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Flags
n_ok     <- sum(summary_df$ok)
n_fail   <- sum(!summary_df$ok)
med_pct  <- round(median(summary_df$pct_retained, na.rm = TRUE), 1)
low_ret  <- sum(summary_df$pct_retained < 70, na.rm = TRUE)

flags <- data.frame(
  step   = "02_demux",
  key    = c("n_samples", "n_ok", "n_failed", "median_pct_retained", "n_low_retention",
             "fwd_primer", "rev_primer"),
  value  = c(nrow(manifest), n_ok, n_fail, med_pct, low_ret, fwd, rev_p),
  status = c("ok", "ok",
             if (n_fail   > 0)  "warn" else "ok",
             if (med_pct  < 70) "warn" else "ok",
             if (low_ret  > 0)  "warn" else "ok",
             "ok", "ok"),
  stringsAsFactors = FALSE
)
write.table(flags, file.path(out_dir, "flags", "flags_02_demux.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("Demux done — %d/%d OK | median retained %.1f%% | %d samples <70%%",
                n_ok, nrow(manifest), med_pct, low_ret))
