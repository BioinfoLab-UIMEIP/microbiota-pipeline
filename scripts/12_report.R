# =============================================================================
# Script:      12_report.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Renders interactive HTML dashboard report from pipeline outputs
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg  <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 12_report.R <config.R> --outdir <dir>")
source(cfg)
if (!exists("SAMPLE_TYPE_FILTER")) SAMPLE_TYPE_FILTER <- NULL
if (!exists("SAMPLE_TYPE_COL"))    SAMPLE_TYPE_COL    <- "sample_type"
od_idx  <- which(args == "--outdir")
out_dir <- if (length(od_idx)) path.expand(args[od_idx[1] + 1]) else stop("--outdir <dir> is required")

# Resolve scripts directory: argv[0] when called via Rscript, else SCRIPTS_DIR env
script_self <- tryCatch(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))][1],
  mustWork = FALSE), error = function(e) "")
script_self <- sub("--file=", "", script_self)
scripts_dir <- if (nzchar(script_self) && !is.na(script_self) && file.exists(script_self)) {
  dirname(script_self)
} else {
  Sys.getenv("SCRIPTS_DIR", dirname(normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE)))
}
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b
template <- file.path(scripts_dir, "report_template.Rmd")
if (!file.exists(template))
  stop("report_template.Rmd not found in: ", scripts_dir,
       "\nSet SCRIPTS_DIR env variable to the pipeline scripts/ directory.")

out_html <- file.path(out_dir, "report.html")
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("rmarkdown required. Install with: conda install -c conda-forge r-rmarkdown")
}

rmarkdown::render(
  input = template,
  output_file = out_html,
  params = list(
    out_dir = out_dir,
    case_label = CASE_LABEL,
    ctrl_label = CTRL_LABEL,
    group_col = GROUP_COL,
    meta_file = METADATA_FILE,
    sample_id_col = SAMPLE_ID_COL,
    sample_type_filter = SAMPLE_TYPE_FILTER,
    sample_type_col = SAMPLE_TYPE_COL
  ),
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)

dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
flags <- data.frame(
  step = "12_report",
  key = "report_html",
  value = "report.html",
  status = "ok",
  stringsAsFactors = FALSE
)
write.table(
  flags,
  file.path(out_dir, "flags", "flags_12_report.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Report saved: ", out_html)
