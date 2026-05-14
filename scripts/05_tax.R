# =============================================================================
# Script:      05_tax.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Taxonomic assignment cascade: SILVA → HOMD → RefSeq BLAST
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg  <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 05_tax.R <config.R> --outdir <dir>")
source(cfg)

if (!exists("USE_COLLAPSED_ASVS")) USE_COLLAPSED_ASVS <- FALSE
if (!exists("TAX_CHUNK_SIZE")) TAX_CHUNK_SIZE <- 250L
if (!exists("TAX_N_JOBS")) TAX_N_JOBS <- 10L
if (!exists("BLASTN_BIN")) BLASTN_BIN <- "blastn"
if (!exists("BSDTAR_BIN")) BSDTAR_BIN <- "bsdtar"
if (!exists("KEEP_TAX_INTERMEDIATES")) KEEP_TAX_INTERMEDIATES <- FALSE
if (!exists("SILVA_TRAINSET")) SILVA_TRAINSET <- ""
if (!exists("SILVA_SPECIES")) SILVA_SPECIES <- ""
if (!exists("HOMD_TRAINSET")) HOMD_TRAINSET <- ""
if (!exists("UNITE_TRAINSET")) UNITE_TRAINSET <- ""
if (!exists("REFSEQ_BLASTDB")) REFSEQ_BLASTDB <- if (exists("REFSEQ_BLAST_DB")) REFSEQ_BLAST_DB else ""
if (!exists("TAX_DB")) TAX_DB <- ""
if (!exists("SPECIES_DB")) SPECIES_DB <- ""

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(parallel)
})

od_idx  <- which(args == "--outdir")
out_dir <- if (length(od_idx)) path.expand(args[od_idx[1] + 1]) else stop("--outdir <dir> is required")

rank_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "Na", "N/A", "NULL")] <- NA_character_
  x
}

blank <- function(x) {
  y <- trim_na(x)
  is.na(y)
}

require_file <- function(path, label) {
  if (!nzchar(path)) stop(label, " must be set in config for 05_tax.R")
  if (!file.exists(path)) stop(label, " not found: ", path)
  path
}

require_bin <- function(bin_path) {
  resolved <- Sys.which(bin_path)
  if (!nzchar(resolved)) stop("Required binary not found in PATH or config: ", bin_path)
  resolved
}

require_blast_db <- function(db_prefix, label) {
  if (!nzchar(db_prefix)) stop(label, " must be set in config for 05_tax.R")
  single_ok <- all(file.exists(paste0(db_prefix, c(".nhr", ".nin", ".nsq"))))
  multi_ok <- any(file.exists(Sys.glob(paste0(db_prefix, ".*.nhr"))))
  if (!single_ok && !multi_ok) {
    stop(label, " BLAST index files not found for prefix: ", db_prefix)
  }
  db_prefix
}

materialize_trainset <- function(path, tax_dir, label, archive_bin = BSDTAR_BIN) {
  require_file(path, label)
  if (!grepl("\\.(7z|zip)$", path, ignore.case = TRUE)) return(path)

  archive_exec <- require_bin(archive_bin)
  cache_dir <- file.path(tax_dir, "reference_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(cache_dir, paste0(tolower(label), "_trainset.fasta"))

  if (file.exists(out_path) && file.info(out_path)$size > 0 && file.info(out_path)$mtime >= file.info(path)$mtime) {
    return(out_path)
  }

  listing <- suppressWarnings(system2(archive_exec, c("-tf", path), stdout = TRUE, stderr = TRUE))
  fasta_name <- listing[grepl("\\.(fa|fasta|fna|fas)(\\.gz)?$", listing, ignore.case = TRUE)][1]
  if (!length(fasta_name) || is.na(fasta_name) || !nzchar(fasta_name)) {
    stop("No FASTA file found inside archive for ", label, ": ", path)
  }

  cmd <- paste(
    shQuote(archive_exec),
    paste(shQuote(c("-xOf", path, fasta_name)), collapse = " "),
    ">",
    shQuote(out_path)
  )
  status <- system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(out_path) || file.info(out_path)$size <= 0) {
    stop("Failed to extract ", label, " archive: ", path)
  }
  out_path
}

split_ids <- function(ids, chunk_size) {
  if (!length(ids)) return(list())
  split(ids, ceiling(seq_along(ids) / chunk_size))
}

empty_rank_df <- function(ids) {
  out <- as.data.frame(matrix(NA_character_, nrow = length(ids), ncol = length(rank_cols)))
  colnames(out) <- rank_cols
  rownames(out) <- ids
  out
}

normalize_kingdom <- function(x) {
  x <- trim_na(x)
  ifelse(
    grepl("bacteria", x, ignore.case = TRUE), "Bacteria",
    ifelse(
      grepl("archaea", x, ignore.case = TRUE), "Archaea",
      ifelse(
        grepl("eukary", x, ignore.case = TRUE), "Eukaryota",
        ifelse(grepl("virus|viruses", x, ignore.case = TRUE), "Viruses", x)
      )
    )
  )
}

normalize_tax_df <- function(x, mode = c("generic", "homd")) {
  mode <- match.arg(mode)
  if (is.null(x) || !nrow(x)) return(empty_rank_df(character()))
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in rank_cols) {
    if (!nm %in% colnames(df)) df[[nm]] <- NA_character_
    df[[nm]] <- trim_na(df[[nm]])
  }
  if (mode == "homd") {
    df$Kingdom <- normalize_kingdom(df$Kingdom)
  }
  df[, rank_cols, drop = FALSE]
}

sanitize_lineage_vec <- function(v) {
  v <- trim_na(v)
  first_blank <- match(TRUE, is.na(v), nomatch = 0L)
  if (first_blank > 0L) v[first_blank:length(v)] <- NA_character_
  v
}

sanitize_lineage_df <- function(df) {
  if (!nrow(df)) return(df)
  out <- t(apply(df[, rank_cols, drop = FALSE], 1, sanitize_lineage_vec))
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(out) <- rank_cols
  rownames(out) <- rownames(df)
  out
}

lineage_depth <- function(v) {
  v <- sanitize_lineage_vec(v)
  sum(!is.na(v))
}

assignment_level <- function(v) {
  depth <- lineage_depth(v)
  if (depth <= 0L) return("Unresolved")
  rank_cols[depth]
}

choose_better <- function(current, candidate, current_source, candidate_source) {
  current_depth <- lineage_depth(current)
  candidate_depth <- lineage_depth(candidate)
  if (candidate_depth > current_depth) return(list(tax = candidate, source = candidate_source))
  list(tax = current, source = current_source)
}

candidate_from_homd <- function(silva_tax, homd_raw) {
  homd_raw <- trim_na(homd_raw)
  standalone <- sanitize_lineage_vec(homd_raw)
  best <- standalone

  graft <- silva_tax
  genus <- homd_raw["Genus"]
  if (!is.na(genus) && lineage_depth(silva_tax) >= 5L) {
    if (is.na(graft["Genus"]) || identical(graft["Genus"], genus)) {
      graft["Genus"] <- genus
    }
  }
  graft <- sanitize_lineage_vec(graft)
  if (lineage_depth(graft) > lineage_depth(best)) best <- graft
  best
}

candidate_from_refseq <- function(silva_tax, sci_name) {
  sci_name <- trim_na(sci_name)
  if (is.na(sci_name)) return(rep(NA_character_, length(rank_cols)))

  genus <- trim_na(sub(" .*", "", sci_name))
  standalone <- setNames(rep(NA_character_, length(rank_cols)), rank_cols)
  standalone["Genus"] <- genus
  standalone["Species"] <- sci_name
  standalone <- sanitize_lineage_vec(standalone)

  best <- standalone
  graft <- silva_tax
  if (lineage_depth(silva_tax) >= 5L && !is.na(genus)) {
    if (is.na(graft["Genus"]) || identical(graft["Genus"], genus)) {
      graft["Genus"] <- genus
      graft["Species"] <- sci_name
    }
  }
  graft <- sanitize_lineage_vec(graft)
  if (lineage_depth(graft) > lineage_depth(best)) best <- graft
  best
}

write_chunk_fastas <- function(chunk_ids, seqs, chunk_dir, prefix) {
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  out <- character(length(chunk_ids))
  for (i in seq_along(chunk_ids)) {
    fa <- file.path(chunk_dir, sprintf("%s_chunk_%03d.fasta", prefix, i))
    writeXStringSet(DNAStringSet(seqs[chunk_ids[[i]]], use.names = TRUE), fa)
    out[i] <- fa
  }
  names(out) <- sprintf("%s_chunk_%03d", prefix, seq_along(chunk_ids))
  out
}

read_tax_chunk_outputs <- function(chunk_ids, chunk_dir, stage_name) {
  if (!length(chunk_ids)) return(empty_rank_df(character()))
  ids <- unlist(chunk_ids, use.names = FALSE)
  out_files <- file.path(chunk_dir, sprintf("%s_chunk_%03d.taxonomy.tsv", stage_name, seq_along(chunk_ids)))
  if (!all(file.exists(out_files))) return(NULL)
  res <- lapply(out_files, function(path) {
    df <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(df) <- df$ASV_ID
    df[, rank_cols, drop = FALSE]
  })
  combined <- do.call(rbind, res)
  combined <- combined[ids, rank_cols, drop = FALSE]
  rownames(combined) <- ids
  combined
}

cleanup_tax_intermediates <- function(tax_dir) {
  if (isTRUE(KEEP_TAX_INTERMEDIATES)) return(invisible(FALSE))
  for (dir_path in c(file.path(tax_dir, "chunks"), file.path(tax_dir, "reference_cache"))) {
    if (dir.exists(dir_path)) {
      unlink(dir_path, recursive = TRUE, force = TRUE)
    }
  }
  invisible(TRUE)
}

write_assignment_outputs <- function(out_df, taxonomy_dir, stage_sources) {
  level_vals <- if ("assignment_level" %in% colnames(out_df)) out_df$assignment_level else vapply(seq_len(nrow(out_df)), function(i) assignment_level(unlist(out_df[i, rank_cols], use.names = TRUE)), character(1))
  level_source <- data.frame(
    assignment_level = level_vals,
    source_db = ifelse(blank(out_df$source_db), "Unresolved", out_df$source_db),
    stringsAsFactors = FALSE
  )
  level_order <- c("Unresolved", rank_cols)
  source_order <- c("SILVA", "HOMD", "RefSeq", "UNITE", "Unresolved")
  level_source$assignment_level <- factor(level_source$assignment_level, levels = level_order)
  level_source$source_db <- factor(level_source$source_db, levels = source_order)

  count_tbl <- as.data.frame(table(level_source$assignment_level, level_source$source_db), stringsAsFactors = FALSE)
  colnames(count_tbl) <- c("assignment_level", "source_db", "n")
  count_tbl <- count_tbl[count_tbl$n > 0, , drop = FALSE]
  count_tbl$prop <- ave(count_tbl$n, count_tbl$assignment_level, FUN = function(x) x / sum(x))
  write.table(count_tbl, file.path(taxonomy_dir, "work", "taxonomy_assignment_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  stage_df <- do.call(rbind, lapply(names(stage_sources), function(stage_name) {
    src <- trim_na(stage_sources[[stage_name]])
    src[is.na(src)] <- "Unresolved"
    tab <- as.data.frame(table(factor(src, levels = c("SILVA", "HOMD", "RefSeq", "UNITE", "Unresolved"))), stringsAsFactors = FALSE)
    colnames(tab) <- c("source_db", "n")
    tab$stage <- stage_name
    tab$prop <- tab$n / sum(tab$n)
    tab
  }))
  stage_df <- stage_df[stage_df$n > 0, , drop = FALSE]
  write.table(stage_df, file.path(taxonomy_dir, "work", "taxonomy_performance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  stage_labels <- c(
    "raw" = "Raw",
    "silva" = "SILVA",
    "silva_homd" = "SILVA + HOMD",
    "silva_homd_refseq" = "SILVA + HOMD + RefSeq",
    "silva_homd_refseq_unite" = "SILVA + HOMD + RefSeq + UNITE",
    "single_db" = "Single DB"
  )
  stage_df$stage_label <- unname(stage_labels[stage_df$stage])
  palette <- c("SILVA" = "#1b9e77", "HOMD" = "#7570b3", "RefSeq" = "#d95f02", "UNITE" = "#66a61e", "Unresolved" = "#bdbdbd")
  mat <- xtabs(prop ~ stage_label + source_db, data = stage_df)
  pdf(file.path(taxonomy_dir, "taxonomy_performance.pdf"), width = 10, height = 6)
  barplot(
    t(mat),
    beside = FALSE,
    horiz = FALSE,
    col = palette[colnames(t(mat))],
    ylim = c(0, 1),
    las = 2,
    ylab = "Proportion of nrASVs",
    xlab = "",
    main = "Taxonomic assignment performance by source"
  )
  legend("topright", legend = colnames(t(mat)), fill = palette[colnames(t(mat))], bty = "n", cex = 0.9)
  dev.off()
}

run_taxonomy_chunks <- function(seqs_subset, trainset, species_db, stage_name, tax_dir, n_jobs, chunk_size, mode = c("generic", "homd")) {
  mode <- match.arg(mode)
  ids <- names(seqs_subset)
  if (!length(ids)) return(empty_rank_df(ids))
  chunk_ids <- split_ids(ids, chunk_size)
  chunk_dir <- file.path(tax_dir, "chunks", stage_name)
  write_chunk_fastas(chunk_ids, seqs_subset, chunk_dir, stage_name)
  cached <- read_tax_chunk_outputs(chunk_ids, chunk_dir, stage_name)
  if (!is.null(cached)) return(cached)

  worker <- function(i) {
    chunk_seq <- seqs_subset[chunk_ids[[i]]]
    tax <- assignTaxonomy(chunk_seq, trainset, multithread = FALSE, minBoot = 50, tryRC = TRUE, verbose = FALSE)
    if (nzchar(species_db) && file.exists(species_db)) {
      tax <- tryCatch(
        addSpecies(tax, species_db, allowMultiple = FALSE, tryRC = TRUE, verbose = FALSE),
        error = function(e) tax
      )
    }
    tax_df <- normalize_tax_df(tax, mode = mode)
    rownames(tax_df) <- names(chunk_seq)
    out_tsv <- file.path(chunk_dir, sprintf("%s_chunk_%03d.taxonomy.tsv", stage_name, i))
    write.table(
      data.frame(ASV_ID = rownames(tax_df), tax_df, stringsAsFactors = FALSE, check.names = FALSE),
      out_tsv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    tax_df
  }

  jobs <- max(1L, min(as.integer(n_jobs), length(chunk_ids)))
  res <- if (jobs > 1L && length(chunk_ids) > 1L) {
    mclapply(seq_along(chunk_ids), worker, mc.cores = jobs)
  } else {
    lapply(seq_along(chunk_ids), worker)
  }

  combined <- do.call(rbind, res)
  combined <- combined[ids, rank_cols, drop = FALSE]
  rownames(combined) <- ids
  combined
}

run_blast_chunks <- function(seqs_subset, blast_db, blast_bin, tax_dir, n_jobs, chunk_size) {
  ids <- names(seqs_subset)
  if (!length(ids)) return(NULL)
  chunk_ids <- split_ids(ids, chunk_size)
  chunk_dir <- file.path(tax_dir, "chunks", "refseq")
  query_fastas <- write_chunk_fastas(chunk_ids, seqs_subset, chunk_dir, "refseq")

  worker <- function(i) {
    out_tsv <- file.path(chunk_dir, sprintf("refseq_chunk_%03d.blast.tsv", i))
    out_log <- file.path(chunk_dir, sprintf("refseq_chunk_%03d.log", i))
    out_empty <- file.path(chunk_dir, sprintf("refseq_chunk_%03d.empty", i))
    out_started <- file.path(chunk_dir, sprintf("refseq_chunk_%03d.started", i))
    out_done <- file.path(chunk_dir, sprintf("refseq_chunk_%03d.done", i))
    if (file.exists(out_done) && file.exists(out_tsv) && file.info(out_tsv)$size > 0) {
      return(read.delim(
        out_tsv,
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        col.names = c("ASV_ID", "sseqid", "pident", "align_len", "evalue", "bitscore", "sscinames", "stitle")
      ))
    }
    if (file.exists(out_done) && file.exists(out_empty)) return(NULL)
    cat(format(Sys.time(), "%F %T"), "start\n", file = out_started)
    blast_args <- c(
      "-task", "megablast",
      "-num_threads", "1",
      "-db", blast_db,
      "-query", query_fastas[[i]],
      "-max_target_seqs", "5",
      "-max_hsps", "1",
      "-evalue", "1e-20",
      "-outfmt", "6 qseqid sseqid pident length evalue bitscore sscinames stitle",
      "-out", out_tsv
    )
    cmd <- paste(
      shQuote(blast_bin),
      paste(shQuote(blast_args), collapse = " "),
      ">",
      shQuote(out_log),
      "2>&1"
    )
    writeLines(cmd, con = file.path(chunk_dir, sprintf("refseq_chunk_%03d.cmd", i)))
    status <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (status != 0) stop("blastn failed for chunk ", i, " with status ", status)
    cat(format(Sys.time(), "%F %T"), "done\n", file = out_done)
    if (!file.exists(out_tsv) || file.info(out_tsv)$size <= 0) {
      file.create(out_empty)
      return(NULL)
    }
    read.delim(
      out_tsv,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      col.names = c("ASV_ID", "sseqid", "pident", "align_len", "evalue", "bitscore", "sscinames", "stitle")
    )
  }

  jobs <- max(1L, min(as.integer(n_jobs), length(chunk_ids)))
  res <- if (jobs > 1L && length(chunk_ids) > 1L) {
    mclapply(seq_along(chunk_ids), worker, mc.cores = jobs)
  } else {
    lapply(seq_along(chunk_ids), worker)
  }

  blast_df <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(blast_df) || !nrow(blast_df)) return(NULL)
  blast_df <- blast_df[!duplicated(blast_df$ASV_ID), , drop = FALSE]
  rownames(blast_df) <- blast_df$ASV_ID
  blast_df
}

taxonomy_dir <- file.path(out_dir, "taxonomy")
tax_dir <- file.path(taxonomy_dir, "work")
fasta_path <- if (USE_COLLAPSED_ASVS) file.path(taxonomy_dir, "redundancy", "rep_seqs_collapsed.fasta") else file.path(out_dir, "dada2", "rep_seqs.fasta")
tax_out <- if (USE_COLLAPSED_ASVS) file.path(taxonomy_dir, "taxonomy_nrASVs.tsv") else file.path(taxonomy_dir, "taxonomy_ASVs.tsv")
dir.create(taxonomy_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tax_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(fasta_path)) stop("Input FASTA not found: ", fasta_path)

seqs_dna <- readDNAStringSet(fasta_path)
seqs <- as.character(seqs_dna)
names(seqs) <- names(seqs_dna)
n_jobs <- max(1L, as.integer(TAX_N_JOBS))
chunk_size <- max(1L, as.integer(TAX_CHUNK_SIZE))

final_tax <- empty_rank_df(names(seqs))
final_tax$sequence <- seqs
final_tax$length <- nchar(seqs)
final_tax$dada_species <- NA_character_
final_tax$refseq_species <- NA_character_
final_tax$unite_species <- NA_character_
final_tax$source_db <- NA_character_

stage_summary <- list()
stage_sources <- list(raw = setNames(rep("Unresolved", length(seqs)), names(seqs)))

if (nzchar(SILVA_TRAINSET)) {
  require_file(SILVA_TRAINSET, "SILVA_TRAINSET")
  if (nzchar(SILVA_SPECIES) && !file.exists(SILVA_SPECIES)) stop("SILVA_SPECIES not found: ", SILVA_SPECIES)

  message("Stage 1/3: SILVA on ", length(seqs), " sequences in ", ceiling(length(seqs) / chunk_size), " chunks with ", n_jobs, " jobs")
  silva_raw <- run_taxonomy_chunks(seqs, SILVA_TRAINSET, SILVA_SPECIES, "silva", tax_dir, n_jobs, chunk_size, mode = "generic")
  silva_clean <- sanitize_lineage_df(silva_raw)
  silva_depth <- vapply(seq_len(nrow(silva_clean)), function(i) lineage_depth(unlist(silva_clean[i, rank_cols])), integer(1))
  names(silva_depth) <- rownames(silva_clean)
  stage_sources[["silva"]] <- ifelse(silva_depth > 0L, "SILVA", "Unresolved")
  stage_summary[["silva_candidates"]] <- length(seqs)
  stage_summary[["silva_chunks"]] <- ceiling(length(seqs) / chunk_size)

  nonspecies_ids <- names(silva_depth)[silva_depth < length(rank_cols)]
  homd_ids <- nonspecies_ids[blank(silva_clean[nonspecies_ids, "Kingdom"]) | silva_clean[nonspecies_ids, "Kingdom"] != "Eukaryota"]
  stage_summary[["homd_candidates"]] <- length(homd_ids)

  homd_raw <- empty_rank_df(character())
  if (length(homd_ids) && nzchar(HOMD_TRAINSET)) {
    require_file(HOMD_TRAINSET, "HOMD_TRAINSET")
    message("Stage 2/3: HOMD on ", length(homd_ids), " sequences in ", ceiling(length(homd_ids) / chunk_size), " chunks with ", n_jobs, " jobs")
    homd_raw <- run_taxonomy_chunks(seqs[homd_ids], HOMD_TRAINSET, "", "homd", tax_dir, n_jobs, chunk_size, mode = "homd")
    stage_summary[["homd_chunks"]] <- ceiling(length(homd_ids) / chunk_size)
  } else {
    stage_summary[["homd_chunks"]] <- 0L
  }

  refseq_ids <- nonspecies_ids
  stage_summary[["refseq_candidates"]] <- length(refseq_ids)
  refseq_blast <- NULL
  if (length(refseq_ids) && nzchar(REFSEQ_BLASTDB)) {
    require_blast_db(REFSEQ_BLASTDB, "REFSEQ_BLASTDB")
    blast_exec <- require_bin(BLASTN_BIN)
    message("Stage 3/3: RefSeq RNA BLAST on ", length(refseq_ids), " sequences in ", ceiling(length(refseq_ids) / chunk_size), " chunks with ", n_jobs, " jobs")
    refseq_blast <- run_blast_chunks(seqs[refseq_ids], REFSEQ_BLASTDB, blast_exec, tax_dir, n_jobs, chunk_size)
    stage_summary[["refseq_chunks"]] <- ceiling(length(refseq_ids) / chunk_size)
  } else {
    stage_summary[["refseq_chunks"]] <- 0L
  }

  ids <- names(seqs)
  source_after_homd <- stage_sources[["silva"]]
  source_after_refseq <- stage_sources[["silva"]]
  source_after_unite <- stage_sources[["silva"]]
  for (id in ids) {
    silva_tax <- unlist(silva_clean[id, rank_cols], use.names = TRUE)
    best_tax <- silva_tax
    best_source <- if (lineage_depth(best_tax) > 0L) "SILVA" else NA_character_

    if (id %in% rownames(homd_raw)) {
      homd_candidate <- candidate_from_homd(silva_tax, unlist(homd_raw[id, rank_cols], use.names = TRUE))
      chosen <- choose_better(best_tax, homd_candidate, best_source, "HOMD")
      best_tax <- chosen$tax
      best_source <- chosen$source
    }
    source_after_homd[id] <- ifelse(blank(best_source), "Unresolved", best_source)

    if (!is.null(refseq_blast) && id %in% rownames(refseq_blast)) {
      sci_name <- trim_na(refseq_blast[id, "sscinames"])
      if (is.na(sci_name)) sci_name <- trim_na(refseq_blast[id, "stitle"])
      refseq_candidate <- candidate_from_refseq(silva_tax, sci_name)
      chosen <- choose_better(best_tax, refseq_candidate, best_source, "RefSeq")
      best_tax <- chosen$tax
      best_source <- chosen$source
      if (!is.na(sci_name)) final_tax[id, "refseq_species"] <- sci_name
    }
    source_after_refseq[id] <- ifelse(blank(best_source), "Unresolved", best_source)

    final_tax[id, rank_cols] <- best_tax
    final_tax[id, "source_db"] <- best_source
    if (!is.na(silva_tax["Species"])) final_tax[id, "dada_species"] <- silva_tax["Species"]
  }
  stage_sources[["silva_homd"]] <- source_after_homd
  stage_sources[["silva_homd_refseq"]] <- source_after_refseq

  unresolved_ids <- ids[vapply(ids, function(id) lineage_depth(unlist(final_tax[id, rank_cols], use.names = TRUE)) == 0L, logical(1))]
  stage_summary[["unite_candidates"]] <- length(unresolved_ids)
  if (length(unresolved_ids) && nzchar(UNITE_TRAINSET)) {
    unite_trainset <- materialize_trainset(UNITE_TRAINSET, tax_dir, "UNITE")
    message("Stage 4/4: UNITE on ", length(unresolved_ids), " unresolved sequences in ", ceiling(length(unresolved_ids) / chunk_size), " chunks with ", n_jobs, " jobs")
    unite_raw <- run_taxonomy_chunks(seqs[unresolved_ids], unite_trainset, "", "unite", tax_dir, n_jobs, chunk_size, mode = "generic")
    unite_clean <- sanitize_lineage_df(unite_raw)
    stage_summary[["unite_chunks"]] <- ceiling(length(unresolved_ids) / chunk_size)
    for (id in unresolved_ids) {
      unite_tax <- unlist(unite_clean[id, rank_cols], use.names = TRUE)
      chosen <- choose_better(unlist(final_tax[id, rank_cols], use.names = TRUE), unite_tax, final_tax[id, "source_db"], "UNITE")
      final_tax[id, rank_cols] <- chosen$tax
      final_tax[id, "source_db"] <- chosen$source
      source_after_unite[id] <- ifelse(blank(chosen$source), "Unresolved", chosen$source)
      if (!is.na(unite_tax["Species"])) final_tax[id, "unite_species"] <- unite_tax["Species"]
    }
    stage_sources[["silva_homd_refseq_unite"]] <- source_after_unite
  } else {
    stage_summary[["unite_chunks"]] <- 0L
    stage_sources[["silva_homd_refseq_unite"]] <- source_after_refseq
  }
} else {
  require_file(TAX_DB, "TAX_DB")
  if (nzchar(SPECIES_DB) && !file.exists(SPECIES_DB)) stop("SPECIES_DB not found: ", SPECIES_DB)
  message("Single-DB mode: ", length(seqs), " sequences in ", ceiling(length(seqs) / chunk_size), " chunks with ", n_jobs, " jobs")
  simple_raw <- run_taxonomy_chunks(seqs, TAX_DB, SPECIES_DB, "simple", tax_dir, n_jobs, chunk_size, mode = "generic")
  simple_clean <- sanitize_lineage_df(simple_raw)
  final_tax[, rank_cols] <- simple_clean[rownames(final_tax), rank_cols, drop = FALSE]
  final_tax$source_db <- ifelse(vapply(seq_len(nrow(simple_clean)), function(i) lineage_depth(unlist(simple_clean[i, rank_cols])) > 0L, logical(1)), basename(TAX_DB), NA_character_)
  final_tax$dada_species <- simple_clean$Species
  stage_summary[["silva_candidates"]] <- 0L
  stage_summary[["silva_chunks"]] <- 0L
  stage_summary[["homd_candidates"]] <- 0L
  stage_summary[["homd_chunks"]] <- 0L
  stage_summary[["refseq_candidates"]] <- 0L
  stage_summary[["refseq_chunks"]] <- 0L
  stage_summary[["unite_candidates"]] <- 0L
  stage_summary[["unite_chunks"]] <- 0L
  stage_sources[["single_db"]] <- ifelse(vapply(seq_len(nrow(simple_clean)), function(i) lineage_depth(unlist(simple_clean[i, rank_cols])) > 0L, logical(1)), basename(TAX_DB), "Unresolved")
}

species_fallback_label <- function(lineage_vec, nr_group) {
  suffix <- sub("^.*?(\\d+)$", "\\1", nr_group)
  if (identical(suffix, nr_group)) suffix <- gsub("[^A-Za-z0-9]+", "_", nr_group)
  clean_tax <- sanitize_lineage_vec(lineage_vec)
  depth <- lineage_depth(clean_tax)
  if (depth <= 0L) return(paste0("nr_", suffix))
  label <- clean_tax[depth]
  label <- gsub("[^A-Za-z0-9]+", "_", label)
  paste0(label, "_", suffix)
}

normalize_binomial_species <- function(genus, lineage_species, refseq_species, unite_species, dada_species, nr_group, lineage_vec, assignment_level, source_db) {
  genus <- trim_na(genus)
  if (!identical(assignment_level, "Species")) {
    return(list(
      label = species_fallback_label(lineage_vec, nr_group),
      source = if (lineage_depth(lineage_vec) > 0L) "fallback_taxon" else "fallback_nr"
    ))
  }

  source_db <- trim_na(source_db)
  source_candidates <- switch(
    source_db,
    "RefSeq" = list(
      list(value = trim_na(lineage_species), source = "lineage"),
      list(value = trim_na(refseq_species), source = "refseq_species")
    ),
    "UNITE" = list(
      list(value = trim_na(lineage_species), source = "lineage"),
      list(value = trim_na(unite_species), source = "unite_species")
    ),
    "HOMD" = list(
      list(value = trim_na(lineage_species), source = "lineage")
    ),
    "SILVA" = list(
      list(value = trim_na(lineage_species), source = "lineage"),
      list(value = trim_na(dada_species), source = "dada_species")
    ),
    list(
      list(value = trim_na(lineage_species), source = "lineage"),
      list(value = trim_na(refseq_species), source = "refseq_species"),
      list(value = trim_na(unite_species), source = "unite_species"),
      list(value = trim_na(dada_species), source = "dada_species")
    )
  )

  candidates <- list(
    source_candidates[[1]]
  )
  if (length(source_candidates) > 1L) {
    candidates <- c(candidates, source_candidates[2:length(source_candidates)])
  }

  clean_tokens <- function(x) {
    x <- trim_na(x)
    if (is.na(x)) return(character())
    x <- gsub("[\\[\\]\\(\\),;:]+", " ", x)
    toks <- strsplit(x, "\\s+")[[1]]
    toks <- toks[nzchar(toks)]
    toks
  }
  valid_species_tokens <- function(x) length(x) >= 2L && all(grepl("[A-Za-z]", x[1:2]))
  valid_epithet <- function(x) length(x) == 1L && grepl("[A-Za-z]", x)

  for (cand in candidates) {
    toks <- clean_tokens(cand$value)
    if (!length(toks)) next
    if (!is.na(genus)) {
      if (valid_species_tokens(toks) && identical(toks[1], genus)) {
        return(list(label = paste(toks[1:2], collapse = " "), source = cand$source))
      }
      if (valid_epithet(toks) && !identical(toks[1], genus)) {
        return(list(label = paste(genus, toks[1]), source = cand$source))
      }
    }
    if (valid_species_tokens(toks)) {
      return(list(label = paste(toks[1:2], collapse = " "), source = cand$source))
    }
  }

  list(label = species_fallback_label(lineage_vec, nr_group), source = if (lineage_depth(lineage_vec) > 0L) "fallback_taxon" else "fallback_nr")
}

raw_species <- final_tax$Species
assignment_level_raw <- vapply(seq_len(nrow(final_tax)), function(i) assignment_level(unlist(final_tax[i, rank_cols], use.names = TRUE)), character(1))

out_df <- final_tax[, c("sequence", "length", rank_cols, "dada_species", "refseq_species", "unite_species", "source_db"), drop = FALSE]
out_df <- data.frame(nr_group = rownames(final_tax), out_df, stringsAsFactors = FALSE, check.names = FALSE)
map_df_expanded <- NULL

if (USE_COLLAPSED_ASVS) {
  rep_path <- file.path(taxonomy_dir, "redundancy", "cluster_representatives.tsv")
  map_path <- file.path(taxonomy_dir, "redundancy", "cluster_map.tsv")
  if (file.exists(rep_path)) {
    rep_df <- read.table(rep_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    rep_df <- rep_df[, c("final_cluster_id", "hc_group", "final_representative_asv_id", "centroid_total_abundance"), drop = FALSE]
    colnames(rep_df)[1] <- "nr_group"
    out_df <- merge(out_df, rep_df, by = "nr_group", all.x = TRUE, sort = FALSE)
    out_df <- out_df[match(rownames(final_tax), out_df$nr_group), , drop = FALSE]
  }
  out_df$ASV_ID <- out_df$final_representative_asv_id
  out_df$representative_asv_id <- out_df$final_representative_asv_id
  if (file.exists(map_path)) {
    map_df_expanded <- read.table(map_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
} else {
  out_df$ASV_ID <- rownames(final_tax)
  out_df$representative_asv_id <- rownames(final_tax)
}

out_df$lineage_species <- raw_species
out_df$assignment_level <- assignment_level_raw
species_norm <- lapply(seq_len(nrow(out_df)), function(i) {
  normalize_binomial_species(
    genus = out_df$Genus[i],
    lineage_species = out_df$lineage_species[i],
    refseq_species = out_df$refseq_species[i],
    unite_species = out_df$unite_species[i],
    dada_species = out_df$dada_species[i],
    nr_group = out_df$nr_group[i],
    lineage_vec = unlist(out_df[i, rank_cols], use.names = TRUE),
    assignment_level = out_df$assignment_level[i],
    source_db = out_df$source_db[i]
  )
})
out_df$Species <- vapply(species_norm, `[[`, character(1), "label")
out_df$species_source <- vapply(species_norm, `[[`, character(1), "source")

rank_cols_no_species <- setdiff(rank_cols, "Species")
out_df <- out_df[, c(
  "ASV_ID", "nr_group", "representative_asv_id", "sequence", "length", rank_cols_no_species,
  "lineage_species", "dada_species", "refseq_species", "unite_species", "Species",
  "species_source", "assignment_level", "source_db",
  intersect(c("hc_group", "centroid_total_abundance", "final_representative_asv_id"), colnames(out_df))
), drop = FALSE]

write.table(out_df, tax_out, sep = "\t", quote = FALSE, row.names = FALSE)

if (USE_COLLAPSED_ASVS && !is.null(map_df_expanded)) {
  expanded_df <- merge(
    map_df_expanded[, c("ASV_ID", "final_cluster_id", "vsearch_centroid_id", "hc_group", "final_representative_asv_id", "original_total_abundance"), drop = FALSE],
    out_df,
    by.x = "final_cluster_id",
    by.y = "nr_group",
    all.x = TRUE,
    sort = FALSE
  )
  colnames(expanded_df)[colnames(expanded_df) == "ASV_ID.x"] <- "original_asv_id"
  colnames(expanded_df)[colnames(expanded_df) == "ASV_ID.y"] <- "ASV_ID"
  colnames(expanded_df)[colnames(expanded_df) == "final_cluster_id"] <- "nr_group"
  keep_cols <- intersect(c(
    "original_asv_id", "nr_group", "ASV_ID", "representative_asv_id", "vsearch_centroid_id", "hc_group",
    "final_representative_asv_id", "original_total_abundance", "sequence", "length", rank_cols_no_species,
    "lineage_species", "dada_species", "refseq_species", "unite_species", "Species",
    "species_source", "assignment_level", "source_db", "centroid_total_abundance"
  ), colnames(expanded_df))
  expanded_df <- expanded_df[, keep_cols, drop = FALSE]
  write.table(expanded_df, file.path(tax_dir, "taxonomy_nrASVs_expanded.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

n_kingdom <- sum(!blank(out_df$Kingdom))
n_genus <- sum(!blank(out_df$Genus))
n_species <- sum(!blank(out_df$lineage_species))
n_valid_species <- sum(out_df$assignment_level == "Species")
n_source_silva <- sum(out_df$source_db == "SILVA", na.rm = TRUE)
n_source_homd <- sum(out_df$source_db == "HOMD", na.rm = TRUE)
n_source_refseq <- sum(out_df$source_db == "RefSeq", na.rm = TRUE)
n_source_unite <- sum(out_df$source_db == "UNITE", na.rm = TRUE)
n_unresolved <- sum(vapply(seq_len(nrow(out_df)), function(i) lineage_depth(unlist(out_df[i, rank_cols])) == 0L, logical(1)))

summary_df <- data.frame(
  metric = c(
    "n_asvs_total",
    "n_with_kingdom",
    "n_with_genus",
    "n_with_species",
    "n_valid_species",
    "n_source_silva",
    "n_source_homd",
    "n_source_refseq",
    "n_source_unite",
    "n_still_unresolved",
    "tax_n_jobs",
    "tax_chunk_size",
    names(stage_summary)
  ),
  value = c(
    nrow(out_df),
    n_kingdom,
    n_genus,
    n_species,
    n_valid_species,
    n_source_silva,
    n_source_homd,
    n_source_refseq,
    n_source_unite,
    n_unresolved,
    n_jobs,
    chunk_size,
    unlist(stage_summary, use.names = FALSE)
  ),
  stringsAsFactors = FALSE
)
write.table(summary_df, file.path(tax_dir, "taxonomy_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write_assignment_outputs(out_df, taxonomy_dir, stage_sources)

flags <- data.frame(
  step = "05_tax",
  key = c("n_asvs", "n_classified_kingdom", "n_classified_genus", "n_valid_species", "tax_n_jobs", "tax_chunk_size", "input_mode"),
  value = c(nrow(out_df), n_kingdom, n_genus, n_valid_species, n_jobs, chunk_size, if (USE_COLLAPSED_ASVS) "collapsed" else "original"),
  status = c(
    "ok",
    if (n_kingdom < nrow(out_df) * 0.8) "warn" else "ok",
    if (n_genus < nrow(out_df) * 0.5) "warn" else "ok",
    if (n_valid_species < nrow(out_df) * 0.05) "warn" else "ok",
    "ok",
    "ok",
    "ok"
  ),
  stringsAsFactors = FALSE
)
write.table(flags, file.path(out_dir, "flags", "flags_05_tax.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cleanup_tax_intermediates(tax_dir)
message("Taxonomy done — ", nrow(out_df), " ASVs | kingdom classified: ", n_kingdom, " | genus classified: ", n_genus, " | valid species: ", n_valid_species)
