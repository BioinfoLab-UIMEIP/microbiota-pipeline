# =============================================================================
# Script:      11_informative.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Random Forest feature selection and LDA discriminant analysis
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 11_informative.R <config.R> --outdir <dir>")
source(cfg)
if (!exists("GROUP_COL")) GROUP_COL <- "Group"
if (!exists("SAMPLE_ID_COL")) SAMPLE_ID_COL <- "Sample_ID"
if (!exists("CASE_LABEL")) CASE_LABEL <- NULL
if (!exists("CTRL_LABEL")) CTRL_LABEL <- NULL
if (!exists("RF_NTREE")) RF_NTREE <- 500L
if (!exists("RF_NFOLD")) RF_NFOLD <- 5L
if (!exists("TOP_N_RF")) TOP_N_RF <- 20L
if (!exists("METADATA_FILE") || !nzchar(METADATA_FILE)) stop("METADATA_FILE must be set in config")
if (is.null(CASE_LABEL) || is.null(CTRL_LABEL)) stop("CASE_LABEL and CTRL_LABEL must be set in config")
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
asv_path <- file.path(out_dir, "filter", "asv_table_filtered.tsv"); tax_path <- file.path(out_dir, "taxonomy", "taxonomy_assigned_collapsed.tsv"); meta_path <- path.expand(METADATA_FILE)
if (!file.exists(asv_path)) stop("asv_table_filtered.tsv not found: ", asv_path)
if (!file.exists(tax_path)) stop("taxonomy_filtered.tsv not found: ", tax_path)
if (!file.exists(meta_path)) stop("Metadata file not found: ", meta_path)

message("Loading informative-feature inputs...")
asv <- read.table(asv_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tax <- read.table(tax_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
sep_m <- if (grepl("\\.csv$", meta_path, ignore.case = TRUE)) "," else "\t"
meta <- read.table(meta_path, sep = sep_m, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample_id" %in% colnames(asv)) stop("asv_table_filtered.tsv must contain sample_id column")
if (!all(c("ASV_ID", "sequence") %in% colnames(tax))) stop("taxonomy_filtered.tsv must contain ASV_ID and sequence columns")
if (!SAMPLE_ID_COL %in% colnames(meta)) stop("Metadata sample ID column not found: ", SAMPLE_ID_COL)
if (!GROUP_COL %in% colnames(meta)) stop("Metadata group column not found: ", GROUP_COL)
meta <- meta[meta[[GROUP_COL]] %in% c(CASE_LABEL, CTRL_LABEL), , drop = FALSE]
if (!exists("SAMPLE_TYPE_COL")) SAMPLE_TYPE_COL <- "sample_type"
if (exists("SAMPLE_TYPE_FILTER") && !is.null(SAMPLE_TYPE_FILTER) &&
    length(SAMPLE_TYPE_FILTER) > 0 && nzchar(SAMPLE_TYPE_FILTER[1]) &&
    SAMPLE_TYPE_COL %in% colnames(meta)) {
  meta <- meta[meta[[SAMPLE_TYPE_COL]] %in% SAMPLE_TYPE_FILTER, , drop = FALSE]
  message("Sample type filter applied: ",
          paste(SAMPLE_TYPE_FILTER, collapse = ", "),
          " (", nrow(meta), " samples retained)")
}
merged <- merge(asv, meta, by.x = "sample_id", by.y = SAMPLE_ID_COL, all.x = FALSE, all.y = FALSE)
if (nrow(merged) < 4L) stop("Need at least 4 samples across case/control after filtering")
grp_factor <- factor(merged[[GROUP_COL]], levels = c(CASE_LABEL, CTRL_LABEL))
if (any(table(grp_factor) < 2L)) stop("Each group must have at least 2 samples for CV")
seq_cols <- setdiff(colnames(asv), "sample_id")
mat <- as.matrix(merged[, seq_cols, drop = FALSE]); storage.mode(mat) <- "numeric"
keep_var <- apply(mat, 2, var) > 0
mat <- mat[, keep_var, drop = FALSE]
if (!ncol(mat)) stop("No non-zero-variance features available after filtering")

clr_transform <- function(m) {
  m_ps <- m
  for (i in seq_len(nrow(m_ps))) {
    nz <- m_ps[i, ] > 0
    if (!any(nz)) stop("Cannot CLR-transform sample with all-zero counts: ", rownames(m_ps)[i] %||% i)
    m_ps[i, !nz] <- min(m_ps[i, nz]) / 2
  }
  log_m <- log(m_ps)
  sweep(log_m, 1, rowMeans(log_m), "-")
}
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x
mat_clr <- clr_transform(mat)

suppressPackageStartupMessages(library(randomForest))
message("Running Random Forest ", RF_NFOLD, "-fold stratified CV...")
n <- nrow(mat); folds <- vector("list", RF_NFOLD)
for (g in levels(grp_factor)) {
  idx_shuf <- sample(which(grp_factor == g))
  for (k in seq_len(RF_NFOLD)) folds[[k]] <- c(folds[[k]], idx_shuf[seq(k, length(idx_shuf), RF_NFOLD)])
}
oof_pred <- character(n); imp_list <- vector("list", RF_NFOLD)
for (k in seq_len(RF_NFOLD)) {
  test_i <- folds[[k]]; train_i <- setdiff(seq_len(n), test_i)
  rf_k <- randomForest(x = mat[train_i, , drop = FALSE], y = grp_factor[train_i], ntree = RF_NTREE, importance = TRUE)
  oof_pred[test_i] <- as.character(predict(rf_k, mat[test_i, , drop = FALSE]))
  imp_list[[k]] <- importance(rf_k, type = 1)[, 1]
}
imp_mat <- do.call(cbind, lapply(imp_list, function(x) setNames(x, names(x))))
imp_avg <- rowMeans(imp_mat, na.rm = TRUE)
imp_df <- data.frame(sequence = names(imp_avg), importance = imp_avg, stringsAsFactors = FALSE)
imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), , drop = FALSE]
cv_acc <- mean(oof_pred == as.character(grp_factor), na.rm = TRUE)
message("RF CV accuracy: ", round(cv_acc * 100, 1), "%")

lda_df <- NULL; lda_scores <- NULL
top_seqs <- head(imp_df$sequence, min(50L, nrow(imp_df)))
if (requireNamespace("MASS", quietly = TRUE)) {
  message("Running LDA on CLR-transformed top RF features...")
  lda_fit <- MASS::lda(mat_clr[, top_seqs, drop = FALSE], grouping = grp_factor)
  lda_pred <- predict(lda_fit, mat_clr[, top_seqs, drop = FALSE])
  lda_df <- data.frame(LD1 = lda_pred$x[, 1], group = as.character(grp_factor), stringsAsFactors = FALSE)
  lda_scores <- data.frame(sequence = top_seqs, LD1 = lda_fit$scaling[, 1], stringsAsFactors = FALSE)
  lda_scores <- lda_scores[order(abs(lda_scores$LD1), decreasing = TRUE), , drop = FALSE]
} else message("MASS not available - skipping LDA")

genus_col <- intersect(c("Genus", "genus"), colnames(tax))[1]; phylum_col <- intersect(c("Phylum", "phylum"), colnames(tax))[1]
annotate_taxa <- function(df) {
  if ("feature_id" %in% colnames(tax) && all(df$sequence %in% tax$feature_id)) {
    idx <- match(df$sequence, tax$feature_id)
  } else if ("nr_group" %in% colnames(tax) && all(df$sequence %in% tax$nr_group)) {
    idx <- match(df$sequence, tax$nr_group)
  } else if (all(df$sequence %in% tax$ASV_ID)) {
    idx <- match(df$sequence, tax$ASV_ID)
  } else if (all(df$sequence %in% tax$sequence)) {
    idx <- match(df$sequence, tax$sequence)
  } else {
    stop("Feature IDs do not match taxonomy feature_id, nr_group, ASV_ID or sequence columns")
  }
  df$feature_id <- df$sequence
  df$sequence <- tax$sequence[idx]
  df$ASV_ID <- tax$ASV_ID[idx]
  if ("nr_group" %in% colnames(tax)) df$nr_group <- tax$nr_group[idx]
  if (!is.na(genus_col)) df$Genus <- tax[[genus_col]][idx]
  if (!is.na(phylum_col)) df$Phylum <- tax[[phylum_col]][idx]
  df
}
imp_df <- annotate_taxa(imp_df)
dir.create(file.path(out_dir, "informative"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)
write.table(imp_df, file.path(out_dir, "informative", "rf_importance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(lda_scores)) write.table(annotate_taxa(lda_scores), file.path(out_dir, "informative", "lda_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

suppressPackageStartupMessages(library(ggplot2))
top_imp <- head(imp_df, TOP_N_RF); top_imp$label <- ifelse(!is.na(top_imp$Genus) & nzchar(top_imp$Genus), top_imp$Genus, top_imp$ASV_ID)
top_feature_ids <- top_imp$feature_id
rel_mat <- sweep(mat, 1, rowSums(mat), "/")
rel_mat[!is.finite(rel_mat)] <- 0
rf_comp <- data.frame(sample_id = merged$sample_id, group = grp_factor, rel_mat[, top_feature_ids, drop = FALSE], check.names = FALSE)
rf_comp_long <- reshape(
  rf_comp,
  varying = top_feature_ids,
  v.names = "rel_abundance",
  timevar = "feature_id",
  times = top_feature_ids,
  direction = "long"
)
rf_comp_long <- rf_comp_long[, c("sample_id", "group", "feature_id", "rel_abundance")]
rf_comp_long$feature_id <- as.character(rf_comp_long$feature_id)
rf_comp_long$label <- top_imp$label[match(rf_comp_long$feature_id, top_imp$feature_id)]
rf_comp_mean <- aggregate(rel_abundance ~ group + label, data = rf_comp_long, FUN = mean)
write.table(rf_comp_mean, file.path(out_dir, "informative", "rf_top_features_group_mean.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
save_pdf_plot <- function(path, plot_obj, width, height) {
  tryCatch({
    pdf(path, width = width, height = height)
    print(plot_obj)
    dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("WARNING: figure failed: ", conditionMessage(e))
  })
}

imp_all <- imp_df
imp_all$label <- ifelse(!is.na(imp_all$Genus) & nzchar(imp_all$Genus),
  imp_all$Genus, imp_all$ASV_ID)
imp_all$label <- factor(make.unique(imp_all$label),
  levels = rev(make.unique(imp_all$label)))
p_all <- ggplot(imp_all, aes(importance, label)) +
  geom_point(size = 3, color = "black") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 2.5),
    axis.text.x = element_text(size = 8)) +
  labs(title = "Species importance - whole model",
    x = "Mean Decrease Accuracy", y = "")
save_pdf_plot(file.path(out_dir, "informative", "rf_importance_all.pdf"), p_all, width = 5, height = 7)

q_breaks <- seq(0.50, 0.95, by = 0.05)
cutoffs <- quantile(imp_df$importance, probs = q_breaks, na.rm = TRUE)
n_feats <- sapply(cutoffs, function(q) sum(imp_df$importance >= q, na.rm = TRUE))
k <- -log(1 - max(0.5, cv_acc - 0.5) / 0.5 + 1e-9) /
  max(1, n_feats[which.min(abs(n_feats - TOP_N_RF))])
k <- max(k, 0.01)
acc_sim <- 0.5 + 0.5 * (1 - exp(-k * n_feats))
gain <- c(0, diff(acc_sim))
opt_idx <- which.max(gain)
cutoff_df <- data.frame(q = q_breaks, n_features = n_feats, accuracy = acc_sim,
  stringsAsFactors = FALSE)
p_cut <- ggplot(cutoff_df, aes(n_features, accuracy)) +
  geom_line() + geom_point(size = 2) +
  geom_point(data = cutoff_df[opt_idx, , drop = FALSE],
    color = "red", size = 4) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_text(size = 8)) +
  labs(title = "RF feature cutoff curve",
    x = "Number of features", y = "Simulated accuracy")
save_pdf_plot(file.path(out_dir, "informative", "rf_cutoff.pdf"), p_cut, width = 5, height = 7)

top_imp$label <- factor(make.unique(top_imp$label),
  levels = rev(make.unique(top_imp$label)))
p_filt <- ggplot(top_imp, aes(importance, label)) +
  geom_point(size = 3, color = "black") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 8)) +
  labs(title = "Most informative species",
    x = "Mean Decrease Accuracy", y = "")
save_pdf_plot(file.path(out_dir, "informative", "rf_importance_filtered.pdf"), p_filt, width = 5, height = 7)

rf_labels <- make.unique(as.character(top_imp$label))
rf_comp_mean$label <- factor(rf_comp_mean$label, levels = rf_labels)
p_comp <- ggplot(rf_comp_mean, aes(group, rel_abundance, fill = label)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = taxa_palette(length(rf_labels))) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(title = "Top RF features: mean relative abundance by group",
    x = "", y = "Fraction", fill = NULL) +
  theme_minimal(base_size = 10)
save_pdf_plot(file.path(out_dir, "informative", "informative_composition.pdf"), p_comp, width = 8, height = 6)

if (!is.null(lda_scores) && !is.null(mat_clr)) {
  top_lda2 <- head(annotate_taxa(lda_scores), TOP_N_RF)
  # Use species feature_id as label (already the collapsed species name)
  top_lda2$label <- make.unique(top_lda2$feature_id)
  top_seqs_lda <- top_lda2$feature_id[top_lda2$feature_id %in% colnames(mat_clr)]

  if (length(top_seqs_lda) >= 2L) {
    case_idx <- which(grp_factor == CASE_LABEL)
    ctrl_idx <- which(grp_factor == CTRL_LABEL)
    # Use TSS-based means so case bars always go right (+) and ctrl bars always go left (-)
    rs <- rowSums(mat); rs[rs == 0] <- 1
    mat_tss_local <- sweep(mat, 1, rs, "/")
    tss_sub       <- mat_tss_local[, top_seqs_lda, drop = FALSE]
    mean_tss_case <- colMeans(tss_sub[case_idx, , drop = FALSE], na.rm = TRUE)
    mean_tss_ctrl <- colMeans(tss_sub[ctrl_idx, , drop = FALSE], na.rm = TRUE)
    total_tss     <- mean_tss_case + mean_tss_ctrl
    total_tss[total_tss == 0] <- 1
    adj_case <-  (mean_tss_case / total_tss) * 100   # always >= 0
    adj_ctrl <- -(mean_tss_ctrl / total_tss) * 100   # always <= 0

    bar_df <- rbind(
      data.frame(species = top_seqs_lda, side = CASE_LABEL,
                 mean_norm_abundance = adj_case, stringsAsFactors = FALSE),
      data.frame(species = top_seqs_lda, side = CTRL_LABEL,
                 mean_norm_abundance = adj_ctrl, stringsAsFactors = FALSE)
    )
    bar_df$LD1 <- top_lda2$LD1[match(bar_df$species, top_lda2$feature_id)]
    bar_df$label <- top_lda2$label[match(bar_df$species, top_lda2$feature_id)]
    bar_df$group <- ifelse(bar_df$LD1 >= 0, CASE_LABEL, CTRL_LABEL)
    # Scale LD1 to ±100 so dots are visible on the same axis range as bars
    ld1_max <- max(abs(bar_df$LD1), na.rm = TRUE)
    if (ld1_max > 0) bar_df$LD1_plot <- (bar_df$LD1 / ld1_max) * 100 else bar_df$LD1_plot <- 0
    case_present_lda <- colSums(mat[which(grp_factor == CASE_LABEL),
                                    top_seqs_lda, drop = FALSE] > 0) > 0
    ctrl_present_lda <- colSums(mat[which(grp_factor == CTRL_LABEL),
                                    top_seqs_lda, drop = FALSE] > 0) > 0
    gs_flag <- as.integer(!(case_present_lda & ctrl_present_lda))
    bar_df$group_specific <- gs_flag[match(bar_df$species, top_seqs_lda)]
    bar_df$group_specific[is.na(bar_df$group_specific)] <- 0L

    n_f <- length(top_seqs_lda)
    p_lda2 <- ggplot(bar_df) +
      geom_point(aes(x = reorder(label, LD1), y = LD1_plot, color = group),
                 stat = "identity", size = 2, show.legend = FALSE) +
      geom_bar(aes(x = reorder(label, LD1), y = mean_norm_abundance, fill = side),
               alpha = 0.5,
               stat = "identity",
               linewidth = ifelse(bar_df$group_specific == 1L, 0.3, 0),
               color = ifelse(bar_df$group_specific == 1L, "darkred", "black")) +
      geom_hline(yintercept = 0, color = "black") +
      coord_flip() +
      scale_fill_manual(values = group_colors, name = "Group") +
      scale_color_manual(values = group_colors, name = "Group") +
      scale_y_continuous(
        labels = abs,
        breaks = seq(-100, 100, by = 25),
        sec.axis = sec_axis(
          transform = ~ . * 1,
          labels = abs,
          breaks = seq(-100, 100, by = 25),
          name = "Adjusted mean normalized abundance (%)"
        )
      ) +
      labs(title = "",
           x = paste0(n_f, " most informative species"),
           y = "Absolute LDA coefficient") +
      theme_minimal(base_size = 9) +
      theme(axis.text.y = element_text(size = 6))

    tryCatch({
      pdf(file.path(out_dir, "informative", "lda_scores.pdf"), width = 5, height = 7)
      print(p_lda2)
      dev.off()
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      message("WARNING: lda_scores.pdf failed: ", conditionMessage(e))
    })
  }
}

if (requireNamespace("vegan", quietly = TRUE)) {
  top_seqs_mat <- mat[, top_seqs[top_seqs %in% colnames(mat)], drop = FALSE]
  if (ncol(top_seqs_mat) >= 2L && nrow(top_seqs_mat) >= 4L) {
    rel_info <- sweep(top_seqs_mat, 1, rowSums(top_seqs_mat) + 1e-9, "/")
    d_info <- vegan::vegdist(rel_info, method = "bray")
    pc_info <- cmdscale(d_info, k = 2, eig = TRUE)
    eig_info <- pc_info$eig; eig_info[eig_info < 0] <- 0
    pct_info <- round(eig_info / sum(eig_info) * 100, 1)
    df_info <- data.frame(PC1 = pc_info$points[, 1], PC2 = pc_info$points[, 2],
      group = as.character(grp_factor), stringsAsFactors = FALSE)
    perm_info <- vegan::adonis2(d_info ~ grp_factor, permutations = 999)
    mrow <- setdiff(rownames(perm_info), c("Residual", "Total"))[1]
    r2_i <- round(perm_info[mrow, "R2"], 3)
    pv_i <- round(perm_info[mrow, "Pr(>F)"], 3)
    ann_info <- paste0("PERMANOVA R^2=", r2_i, "  p=", pv_i)
    p_ib <- ggplot(df_info, aes(PC1, PC2, color = group)) +
      stat_ellipse(level = 0.95, linetype = 2) +
      geom_point(size = 2, alpha = 0.8) +
      annotate("text", x = -Inf, y = Inf, label = ann_info,
        hjust = -0.05, vjust = 1.3, size = 3) +
      scale_color_manual(values = group_colors) +
      theme_minimal(base_size = 11) +
      labs(title = "PCoA Bray-Curtis - informative species",
        x = paste0("PCoA1 (", pct_info[1], "%)"),
        y = paste0("PCoA2 (", pct_info[2], "%)"),
        color = GROUP_COL)
    save_pdf_plot(file.path(out_dir, "informative", "informative_beta.pdf"), p_ib, width = 7, height = 5)
    write.table(
      data.frame(
        sample_id = rownames(df_info),
        PC1 = df_info$PC1,
        PC2 = df_info$PC2,
        group = df_info$group,
        pct1 = pct_info[1],
        pct2 = pct_info[2],
        permanova_R2 = r2_i,
        permanova_p = pv_i,
        stringsAsFactors = FALSE
      ),
      file.path(out_dir, "informative", "informative_pcoa.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

flags <- data.frame(step = "11_informative", key = c("n_features_rf", "rf_cv_accuracy", "lda_available"), value = c(nrow(imp_df), sprintf("%.2f", cv_acc), if (is.null(lda_scores)) "no" else "yes"), status = c("ok", if (cv_acc < 0.6) "warn" else "ok", "ok"), stringsAsFactors = FALSE)
write.table(flags, file.path(out_dir, "flags", "flags_11_informative.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Informative taxa done: ", nrow(imp_df), " features retained")
