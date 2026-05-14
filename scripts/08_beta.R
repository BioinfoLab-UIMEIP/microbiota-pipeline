# =============================================================================
# Script:      08_beta.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Beta diversity PCoA and PERMANOVA with Bray-Curtis and Jaccard distances
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 08_beta.R <config.R> --outdir <dir>")
source(cfg)
if (!exists("GROUP_COL")) GROUP_COL <- "Group"
if (!exists("SAMPLE_ID_COL")) SAMPLE_ID_COL <- "Sample_ID"
if (!exists("CASE_LABEL")) CASE_LABEL <- NULL
if (!exists("CTRL_LABEL")) CTRL_LABEL <- NULL
if (!exists("PERMANOVA_NPERM")) PERMANOVA_NPERM <- 999L
if (!exists("MAFFT_BIN")) MAFFT_BIN <- "mafft"
if (!exists("MAFFT_THREADS")) MAFFT_THREADS <- 4L
if (!exists("UNIFRAC_MAX_FEAT")) UNIFRAC_MAX_FEAT <- 500L
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
tss_path <- file.path(out_dir, "filter", "asv_table_tss.tsv"); meta_path <- path.expand(METADATA_FILE)
if (!file.exists(tss_path)) stop("asv_table_tss.tsv not found: ", tss_path)
if (!file.exists(meta_path)) stop("Metadata file not found: ", meta_path)

message("Loading TSS ASV table...")
tss <- read.table(tss_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample_id" %in% colnames(tss)) stop("asv_table_tss.tsv must contain sample_id column")
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
merged <- merge(tss, meta, by.x = "sample_id", by.y = SAMPLE_ID_COL)
if (!nrow(merged)) stop("No overlapping samples between asv_table_tss.tsv and metadata")
mat <- as.matrix(merged[, colnames(tss)[-1], drop = FALSE]); rownames(mat) <- merged$sample_id
message("Computing distances...")
suppressPackageStartupMessages(library(vegan))
bc_d <- vegdist(mat, method = "bray")
jac_d <- vegdist(mat > 0, method = "jaccard", binary = TRUE)
message("Running PCoA...")
pcoa_bc <- cmdscale(bc_d, k = 2, eig = TRUE); pcoa_jac <- cmdscale(jac_d, k = 2, eig = TRUE)

dir.create(file.path(out_dir, "beta"), recursive = TRUE, showWarnings = FALSE)
grp_ok <- GROUP_COL %in% colnames(merged) && length(unique(merged[[GROUP_COL]])) > 1
perm_bc_r2 <- perm_bc_p <- perm_jac_r2 <- perm_jac_p <- NA_real_
disp_bc_p <- disp_jac_p <- NA_real_
perm_wuf_r2 <- perm_wuf_p <- perm_uwuf_r2 <- perm_uwuf_p <- NA_real_
disp_wuf_p <- disp_uwuf_p <- NA_real_
if (grp_ok) {
  grp <- merged[[GROUP_COL]]; set.seed(42)
  perm_bc <- adonis2(bc_d ~ grp, permutations = PERMANOVA_NPERM)
  perm_jac <- adonis2(jac_d ~ grp, permutations = PERMANOVA_NPERM)
  model_row_bc <- setdiff(rownames(perm_bc), c("Residual", "Total"))[1]
  model_row_jac <- setdiff(rownames(perm_jac), c("Residual", "Total"))[1]
  perm_bc_r2 <- unname(perm_bc[model_row_bc, "R2"]); perm_bc_p <- unname(perm_bc[model_row_bc, "Pr(>F)"])
  perm_jac_r2 <- unname(perm_jac[model_row_jac, "R2"]); perm_jac_p <- unname(perm_jac[model_row_jac, "Pr(>F)"])
  bd_bc <- betadisper(bc_d, grp)
  bd_jac <- betadisper(jac_d, grp)
  disp_bc_p <- permutest(bd_bc, permutations = PERMANOVA_NPERM)$tab[1, "Pr(>F)"]
  disp_jac_p <- permutest(bd_jac, permutations = PERMANOVA_NPERM)$tab[1, "Pr(>F)"]
  write.table(data.frame(metric = c("bray_curtis", "jaccard"), R2 = c(perm_bc_r2, perm_jac_r2), p_value = c(perm_bc_p, perm_jac_p), stringsAsFactors = FALSE), file.path(out_dir, "beta", "permanova.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(metric = c("bray_curtis", "jaccard"), p_value = c(disp_bc_p, disp_jac_p), stringsAsFactors = FALSE), file.path(out_dir, "beta", "betadisper.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("PERMANOVA Bray-Curtis R2=", round(perm_bc_r2, 3), " p=", perm_bc_p)
  suppressPackageStartupMessages(library(ggplot2))
  has_gridextra <- requireNamespace("gridExtra", quietly = TRUE)
  if (has_gridextra) suppressPackageStartupMessages(library(gridExtra))
  mk <- function(p, g, s) {
    eig <- p$eig; eig[eig < 0] <- 0
    list(
      df = data.frame(PC1 = p$points[, 1], PC2 = p$points[, 2], group = g, sample_id = s, stringsAsFactors = FALSE),
      pct = round(eig / sum(eig) * 100, 1)
    )
  }
  plt <- function(res, title, r2 = NA, pval = NA) {
    ann <- if (!is.na(r2) && !is.na(pval))
      paste0("PERMANOVA R^2=", round(r2, 3), "  p=", round(pval, 3))
    else ""
    ggplot(res$df, aes(PC1, PC2, color = group)) +
      stat_ellipse(level = 0.95, linetype = 2) +
      geom_point(size = 2, alpha = 0.8) +
      annotate("text", x = -Inf, y = Inf, label = ann,
        hjust = -0.05, vjust = 1.3, size = 3) +
      scale_color_manual(values = group_colors) +
      labs(title = title,
        x = paste0("PCoA1 (", res$pct[1], "%)"),
        y = paste0("PCoA2 (", res$pct[2], "%)"),
        color = GROUP_COL) +
      theme_minimal(base_size = 11)
  }
  p_bc <- plt(mk(pcoa_bc, grp, merged$sample_id), "PCoA Bray-Curtis", perm_bc_r2, perm_bc_p)
  p_jac <- plt(mk(pcoa_jac, grp, merged$sample_id), "PCoA Jaccard", perm_jac_r2, perm_jac_p)
  if (has_gridextra) {
    pdf(file.path(out_dir, "beta", "pcoa.pdf"), width = 14, height = 5)
    gridExtra::grid.arrange(p_bc, p_jac, nrow = 1)
    dev.off()
  } else {
    pdf(file.path(out_dir, "beta", "pcoa.pdf"), width = 7, height = 5)
    print(p_bc)
    print(p_jac)
    dev.off()
  }

  wuf_d <- uwuf_d <- NULL
  pcoa_wuf <- pcoa_uwuf <- NULL

  if (requireNamespace("ape", quietly = TRUE) && requireNamespace("GUniFrac", quietly = TRUE)) {
    feat_seqs <- colnames(mat)
    n_feat <- length(feat_seqs)
    is_dna_cols <- mean(grepl("^[ACGTacgt]+$", feat_seqs)) > 0.9

    build_unifrac_tree <- function(seqs_vec, ids_vec) {
      n <- length(seqs_vec)
      if (n > UNIFRAC_MAX_FEAT) stop(n, " features > UNIFRAC_MAX_FEAT (", UNIFRAC_MAX_FEAT, ").")
      short_ids <- paste0("F", seq_len(n))
      id_map    <- setNames(ids_vec, short_ids)
      tmp_fa    <- tempfile(fileext = ".fasta")
      tmp_aln   <- tempfile(fileext = ".aln.fasta")
      writeLines(paste0(">", short_ids, "\n", seqs_vec), tmp_fa)
      mafft_cmd <- paste(shQuote(MAFFT_BIN), "--auto --thread", MAFFT_THREADS,
                         "--quiet", shQuote(tmp_fa), ">", shQuote(tmp_aln))
      ret <- system(mafft_cmd, intern = FALSE)
      if (ret != 0L) stop("MAFFT returned exit code ", ret)
      aln  <- ape::read.dna(tmp_aln, format = "fasta")
      dmat <- ape::dist.dna(aln, model = "raw", pairwise.deletion = TRUE)
      if (any(is.na(dmat))) dmat[is.na(dmat)] <- max(dmat, na.rm = TRUE)
      tree <- ape::nj(dmat)
      tree <- ape::root(tree, outgroup = short_ids[1], resolve.root = TRUE)
      tree$tip.label <- id_map[tree$tip.label]
      tree
    }

    if (!is_dna_cols) {
      tax_path <- file.path(out_dir, "filter", "taxonomy_filtered.tsv")
      if (!file.exists(tax_path)) {
        message("UniFrac skipped: taxonomy_filtered.tsv not found (expected at ", tax_path, ").")
      } else {
        tryCatch({
          message("Feature columns are taxonomy labels — loading representative sequences from taxonomy_filtered.tsv...")
          tax <- read.table(tax_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                            quote = "", comment.char = "")
          req_tax_cols <- c("feature_id", "sequence", "centroid_total_abundance")
          if (!all(req_tax_cols %in% colnames(tax)))
            stop("taxonomy_filtered.tsv missing columns: ",
                 paste(setdiff(req_tax_cols, colnames(tax)), collapse = ", "))
          tax <- tax[tax$feature_id %in% feat_seqs & nchar(trimws(tax$sequence)) > 0, ]
          # Per species: pick ASV with highest centroid_total_abundance as representative
          rep_seq <- tapply(seq_len(nrow(tax)), tax$feature_id, function(idx) {
            sub <- tax[idx, , drop = FALSE]
            sub$sequence[which.max(sub$centroid_total_abundance)]
          })
          valid_ids  <- intersect(feat_seqs, names(rep_seq))
          valid_seqs <- unname(rep_seq[valid_ids])
          if (length(valid_ids) < 4L)
            stop("Fewer than 4 species have representative sequences in taxonomy_filtered.tsv.")
          message("Building phylogenetic tree from ", length(valid_ids),
                  " representative sequences (species-level, most abundant ASV per species)...")
          tree <- build_unifrac_tree(valid_seqs, valid_ids)
          otu_mat <- mat[, valid_ids, drop = FALSE]
          uf_res  <- GUniFrac::GUniFrac(otu_mat, tree, alpha = c(0, 1))
          wuf_d   <- as.dist(uf_res$unifracs[, , "d_1"])
          uwuf_d  <- as.dist(uf_res$unifracs[, , "d_UW"])
          message("UniFrac distances computed successfully (", length(valid_ids),
                  " species, representative sequences).")
        }, error = function(e) {
          message("WARNING: UniFrac skipped — ", conditionMessage(e))
        })
      }
    } else if (n_feat >= 4L && n_feat <= UNIFRAC_MAX_FEAT) {
      tryCatch({
        message("Building phylogenetic tree for UniFrac (", n_feat, " ASV sequences)...")
        tree   <- build_unifrac_tree(feat_seqs, feat_seqs)
        otu_mat <- mat
        uf_res  <- GUniFrac::GUniFrac(otu_mat, tree, alpha = c(0, 1))
        wuf_d   <- as.dist(uf_res$unifracs[, , "d_1"])
        uwuf_d  <- as.dist(uf_res$unifracs[, , "d_UW"])
        message("UniFrac distances computed successfully.")
      }, error = function(e) {
        message("WARNING: UniFrac skipped — ", conditionMessage(e))
      })
    } else {
      if (n_feat < 4L) message("UniFrac skipped: fewer than 4 features.")
      if (n_feat > UNIFRAC_MAX_FEAT) message("UniFrac skipped: ", n_feat,
                                             " features > UNIFRAC_MAX_FEAT (", UNIFRAC_MAX_FEAT, ").")
    }
  } else {
    missing_pkgs <- c(
      if (!requireNamespace("ape", quietly = TRUE)) "ape",
      if (!requireNamespace("GUniFrac", quietly = TRUE)) "GUniFrac"
    )
    message("UniFrac skipped: missing packages: ", paste(missing_pkgs, collapse = ", "))
  }

  if (!is.null(wuf_d) && !is.null(uwuf_d)) {
    tryCatch({
      pcoa_wuf <- cmdscale(wuf_d, k = 2, eig = TRUE)
      pcoa_uwuf <- cmdscale(uwuf_d, k = 2, eig = TRUE)

      run_perm_bd <- function(d, grp) {
        pm <- adonis2(d ~ grp, permutations = PERMANOVA_NPERM)
        mrow <- setdiff(rownames(pm), c("Residual", "Total"))[1]
        r2 <- unname(pm[mrow, "R2"])
        pv <- unname(pm[mrow, "Pr(>F)"])
        bd <- betadisper(d, grp)
        bp <- permutest(bd, permutations = PERMANOVA_NPERM)$tab[1, "Pr(>F)"]
        list(r2 = r2, p = pv, disp_p = bp)
      }
      set.seed(42)
      res_wuf <- run_perm_bd(wuf_d, grp)
      res_uwuf <- run_perm_bd(uwuf_d, grp)
      perm_wuf_r2 <- res_wuf$r2; perm_wuf_p <- res_wuf$p; disp_wuf_p <- res_wuf$disp_p
      perm_uwuf_r2 <- res_uwuf$r2; perm_uwuf_p <- res_uwuf$p; disp_uwuf_p <- res_uwuf$disp_p
      message("UniFrac PERMANOVA done. Weighted R2=", round(perm_wuf_r2, 3), " Unweighted R2=", round(perm_uwuf_r2, 3))

      perm_path <- file.path(out_dir, "beta", "permanova.tsv")
      perm_tbl <- read.table(perm_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      perm_new <- data.frame(
        metric = c("weighted_unifrac", "unweighted_unifrac"),
        R2 = c(perm_wuf_r2, perm_uwuf_r2),
        p_value = c(perm_wuf_p, perm_uwuf_p),
        stringsAsFactors = FALSE
      )
      write.table(rbind(perm_tbl, perm_new), perm_path, sep = "\t", quote = FALSE, row.names = FALSE)

      disp_path <- file.path(out_dir, "beta", "betadisper.tsv")
      disp_tbl <- read.table(disp_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      disp_new <- data.frame(
        metric = c("weighted_unifrac", "unweighted_unifrac"),
        p_value = c(disp_wuf_p, disp_uwuf_p),
        stringsAsFactors = FALSE
      )
      write.table(rbind(disp_tbl, disp_new), disp_path, sep = "\t", quote = FALSE, row.names = FALSE)

      save_pcoa_tsv <- function(pcoa_obj, grp_vec, sid_vec, metric_name) {
        eig <- pcoa_obj$eig; eig[eig < 0] <- 0
        pct <- round(100 * eig / sum(eig), 1)
        df <- data.frame(
          sample_id = sid_vec,
          PC1 = pcoa_obj$points[, 1],
          PC2 = pcoa_obj$points[, 2],
          group = as.character(grp_vec),
          pct1 = pct[1],
          pct2 = pct[2],
          stringsAsFactors = FALSE
        )
        write.table(df, file.path(out_dir, "beta", paste0("pcoa_", metric_name, ".tsv")),
          sep = "\t", quote = FALSE, row.names = FALSE)
        invisible(df)
      }
      df_wuf <- save_pcoa_tsv(pcoa_wuf, grp, merged$sample_id, "wunifrac")
      df_uwuf <- save_pcoa_tsv(pcoa_uwuf, grp, merged$sample_id, "uwunifrac")

      mk_wuf <- function(df, pct, title, r2, pval) {
        ann <- paste0("PERMANOVA R^2=", round(r2, 3), "  p=", round(pval, 3))
        ggplot(df, aes(PC1, PC2, color = group)) +
          stat_ellipse(level = 0.95, linetype = 2) +
          geom_point(size = 2, alpha = 0.8) +
          annotate("text", x = -Inf, y = Inf, label = ann, hjust = -0.05, vjust = 1.3, size = 3) +
          scale_color_manual(values = group_colors) +
          labs(title = title,
            x = paste0("PCoA1 (", pct[1], "%)"),
            y = paste0("PCoA2 (", pct[2], "%)"),
            color = GROUP_COL) +
          theme_minimal(base_size = 11)
      }
      eig_wuf <- pcoa_wuf$eig; eig_wuf[eig_wuf < 0] <- 0
      eig_uwuf <- pcoa_uwuf$eig; eig_uwuf[eig_uwuf < 0] <- 0
      pct_wuf <- round(100 * eig_wuf / sum(eig_wuf), 1)
      pct_uwuf <- round(100 * eig_uwuf / sum(eig_uwuf), 1)
      p_wuf <- mk_wuf(df_wuf, pct_wuf, "PCoA Weighted UniFrac", perm_wuf_r2, perm_wuf_p)
      p_uwuf <- mk_wuf(df_uwuf, pct_uwuf, "PCoA Unweighted UniFrac", perm_uwuf_r2, perm_uwuf_p)

      if (has_gridextra) {
        pdf(file.path(out_dir, "beta", "pcoa.pdf"), width = 14, height = 10)
        gridExtra::grid.arrange(p_bc, p_jac, p_wuf, p_uwuf, nrow = 2)
      } else {
        pdf(file.path(out_dir, "beta", "pcoa.pdf"), width = 7, height = 5)
        print(p_bc); print(p_jac); print(p_wuf); print(p_uwuf)
      }
      dev.off()
      message("UniFrac PCoA saved.")
    }, error = function(e) {
      message("WARNING: UniFrac PCoA/PERMANOVA failed — ", conditionMessage(e))
    })
  }
} else {
  message("Skipping PERMANOVA and plots: GROUP_COL missing or only one group")
}

dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
flags <- data.frame(
  step = "08_beta",
  key = c("n_samples",
          "permanova_bc_R2", "permanova_bc_p", "betadisper_bc_p",
          "permanova_jac_R2", "permanova_jac_p", "betadisper_jac_p",
          "permanova_wuf_R2", "permanova_wuf_p", "betadisper_wuf_p",
          "permanova_uwuf_R2", "permanova_uwuf_p", "betadisper_uwuf_p"),
  value = c(nrow(merged),
            perm_bc_r2, perm_bc_p, disp_bc_p,
            perm_jac_r2, perm_jac_p, disp_jac_p,
            perm_wuf_r2, perm_wuf_p, disp_wuf_p,
            perm_uwuf_r2, perm_uwuf_p, disp_uwuf_p),
  status = c(
    "ok",
    if (is.na(perm_bc_r2) || perm_bc_r2 < 0.05) "warn" else "ok",
    if (is.na(perm_bc_p) || perm_bc_p > 0.05) "warn" else "ok",
    if (is.na(disp_bc_p) || disp_bc_p < 0.05) "warn" else "ok",
    if (is.na(perm_jac_r2) || perm_jac_r2 < 0.05) "warn" else "ok",
    if (is.na(perm_jac_p) || perm_jac_p > 0.05) "warn" else "ok",
    if (is.na(disp_jac_p) || disp_jac_p < 0.05) "warn" else "ok",
    if (is.na(perm_wuf_r2)) "skip" else if (perm_wuf_r2 < 0.05) "warn" else "ok",
    if (is.na(perm_wuf_p)) "skip" else if (perm_wuf_p > 0.05) "warn" else "ok",
    if (is.na(disp_wuf_p)) "skip" else if (disp_wuf_p < 0.05) "warn" else "ok",
    if (is.na(perm_uwuf_r2)) "skip" else if (perm_uwuf_r2 < 0.05) "warn" else "ok",
    if (is.na(perm_uwuf_p)) "skip" else if (perm_uwuf_p > 0.05) "warn" else "ok",
    if (is.na(disp_uwuf_p)) "skip" else if (disp_uwuf_p < 0.05) "warn" else "ok"
  ),
  stringsAsFactors = FALSE
)
write.table(flags, file.path(out_dir, "flags", "flags_08_beta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Beta diversity done: ", nrow(merged), " samples")
