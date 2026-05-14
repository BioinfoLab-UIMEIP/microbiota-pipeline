# =============================================================================
# Script:      09_composition.R
# Pipeline:    16S Microbiota Pipeline
# Author:      Roberto C. Torres, PhD. <torres.roberto.c@gmail.com>
# Institution: Bioinformatics Lab, Infectious Diseases Research Unit,
#              CMN SXXI, IMSS, Mexico City
# Description: Taxonomic composition plots and group-specific species identification
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 09_composition.R <config.R> --outdir <dir>")
source(cfg)
if (!exists("GROUP_COL")) GROUP_COL <- "Group"
if (!exists("SAMPLE_ID_COL")) SAMPLE_ID_COL <- "Sample_ID"
if (!exists("CASE_LABEL")) CASE_LABEL <- "Case"
if (!exists("CTRL_LABEL")) CTRL_LABEL <- "Control"
if (!exists("TOP_N_PHYLA")) TOP_N_PHYLA <- 10L
if (!exists("TOP_N_GENERA")) TOP_N_GENERA <- 15L
if (!exists("TOP_N_SPECIES")) TOP_N_SPECIES <- 20L
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
tss_path <- file.path(out_dir, "filter", "asv_table_tss.tsv")
tax_path <- file.path(out_dir, "taxonomy", "taxonomy_assigned_collapsed.tsv")
meta_path <- path.expand(METADATA_FILE)
if (!file.exists(tss_path)) stop("asv_table_tss.tsv not found: ", tss_path)
if (!file.exists(tax_path)) stop("taxonomy_assigned_collapsed.tsv not found: ", tax_path)
if (!file.exists(meta_path)) stop("Metadata file not found: ", meta_path)

message("Loading composition inputs...")
tss <- read.table(tss_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tax <- read.table(tax_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
sep_m <- if (grepl("\\.csv$", meta_path, ignore.case = TRUE)) "," else "\t"
meta <- read.table(meta_path, sep = sep_m, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample_id" %in% colnames(tss)) stop("asv_table_tss.tsv must contain sample_id column")
if (!SAMPLE_ID_COL %in% colnames(meta)) stop("Metadata sample ID column not found: ", SAMPLE_ID_COL)
if (!GROUP_COL %in% colnames(meta)) stop("Metadata group column not found: ", GROUP_COL)

meta <- meta[!is.na(meta[[GROUP_COL]]) & meta[[GROUP_COL]] %in% c(CTRL_LABEL, CASE_LABEL), , drop = FALSE]
if (!exists("SAMPLE_TYPE_COL")) SAMPLE_TYPE_COL <- "sample_type"
if (exists("SAMPLE_TYPE_FILTER") && !is.null(SAMPLE_TYPE_FILTER) &&
    length(SAMPLE_TYPE_FILTER) > 0 && nzchar(SAMPLE_TYPE_FILTER[1]) &&
    SAMPLE_TYPE_COL %in% colnames(meta)) {
  meta <- meta[meta[[SAMPLE_TYPE_COL]] %in% SAMPLE_TYPE_FILTER, , drop = FALSE]
  message("Sample type filter applied: ",
          paste(SAMPLE_TYPE_FILTER, collapse = ", "),
          " (", nrow(meta), " samples retained)")
}
meta[[GROUP_COL]] <- factor(meta[[GROUP_COL]], levels = c(CTRL_LABEL, CASE_LABEL))
meta <- meta[!duplicated(meta[[SAMPLE_ID_COL]]), , drop = FALSE]
tss <- tss[tss$sample_id %in% meta[[SAMPLE_ID_COL]], , drop = FALSE]
meta <- meta[match(tss$sample_id, meta[[SAMPLE_ID_COL]]), , drop = FALSE]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

feature_cols <- setdiff(colnames(tss), "sample_id")
match_feature_tax <- function(feature_ids, tax_df) {
  if ("feature_id" %in% colnames(tax_df) && all(feature_ids %in% tax_df$feature_id)) return(match(feature_ids, tax_df$feature_id))
  if ("nr_group" %in% colnames(tax_df) && all(feature_ids %in% tax_df$nr_group)) return(match(feature_ids, tax_df$nr_group))
  if (all(feature_ids %in% tax_df$ASV_ID)) return(match(feature_ids, tax_df$ASV_ID))
  if (all(feature_ids %in% tax_df$sequence)) return(match(feature_ids, tax_df$sequence))
  stop("Feature columns do not match taxonomy feature_id, nr_group, ASV_ID or sequence columns")
}
tax_idx <- match_feature_tax(feature_cols, tax)

aggregate_rank <- function(rank_col) {
  rank_vals <- tax[[rank_col]][tax_idx]
  rank_vals[is.na(rank_vals) | !nzchar(rank_vals)] <- "Unclassified"
  mat <- as.matrix(tss[, feature_cols, drop = FALSE])
  agg <- rowsum(t(mat), group = rank_vals, reorder = FALSE)
  agg <- t(agg)
  data.frame(sample_id = tss$sample_id, agg, check.names = FALSE, stringsAsFactors = FALSE)
}

collapse_top_n <- function(df, top_n) {
  taxa_cols <- setdiff(colnames(df), "sample_id")
  means <- colMeans(df[, taxa_cols, drop = FALSE], na.rm = TRUE)
  keep <- names(sort(means, decreasing = TRUE))[seq_len(min(top_n, length(means)))]
  out <- df[, c("sample_id", keep), drop = FALSE]
  other <- setdiff(taxa_cols, keep)
  if (length(other)) out$Other <- rowSums(df[, other, drop = FALSE])
  out
}

to_long <- function(df) {
  long <- pivot_longer(df, cols = -sample_id, names_to = "taxon", values_to = "abundance")
  long[[GROUP_COL]] <- meta[[GROUP_COL]][match(long$sample_id, meta[[SAMPLE_ID_COL]])]
  long
}

plot_sample_facets <- function(df, title) {
  long <- to_long(df)
  n_taxa <- length(unique(long$taxon))
  ord <- meta$sample_id[order(meta[[GROUP_COL]], meta[[SAMPLE_ID_COL]])]
  long$sample_id <- factor(long$sample_id, levels = ord)
  ggplot(long, aes(sample_id, abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "fill", width = 0.95) +
    scale_fill_manual(values = taxa_palette(n_taxa)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    facet_grid(. ~ .data[[GROUP_COL]], scales = "free_x", space = "free_x") +
    labs(title = title, x = "", y = "Relative abundance", fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.spacing.x = unit(0.3, "lines"))
}

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-3) return("p < 0.001")
  paste0("p = ", signif(p, 3))
}

species_df <- data.frame(sample_id = tss$sample_id, tss[, feature_cols, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
phylum_df <- collapse_top_n(aggregate_rank("Phylum"), TOP_N_PHYLA)
genus_df <- collapse_top_n(aggregate_rank("Genus"), TOP_N_GENERA)
species_plot_df <- collapse_top_n(species_df, TOP_N_SPECIES)

group_mean_species <- to_long(species_df) %>%
  group_by(.data[[GROUP_COL]], taxon) %>%
  summarise(mean_rel_abundance = mean(abundance, na.rm = TRUE), .groups = "drop")
top20_group_species <- group_mean_species %>%
  group_by(taxon) %>%
  summarise(score = mean(mean_rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  slice_max(score, n = TOP_N_SPECIES, with_ties = FALSE) %>%
  pull(taxon)
group_mean_species <- group_mean_species %>% filter(taxon %in% top20_group_species)
top_sp_levels <- top20_group_species
comp_df <- group_mean_species[group_mean_species$taxon %in% top_sp_levels, , drop = FALSE]
comp_df$taxon <- factor(comp_df$taxon, levels = sort(unique(comp_df$taxon)))

top20_plot <- ggplot(comp_df,
    aes(x = .data[[GROUP_COL]], y = mean_rel_abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8, colour = NA) +
  scale_fill_manual(values = taxa_palette(length(levels(comp_df$taxon))),
                    name = paste0(TOP_N_SPECIES, " most abundant species")) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(x = "", y = "Fraction (group mean)",
       title = paste0("Top ", TOP_N_SPECIES, " species by mean relative abundance")) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title   = element_text(size = 14),
        legend.text  = element_text(size = 9),
        legend.title = element_text(size = 11))

top20_long <- to_long(species_df) %>%
  filter(taxon %in% top20_group_species) %>%
  mutate(taxon = factor(taxon, levels = top20_group_species))
top20_pvals <- top20_long %>%
  group_by(taxon) %>%
  summarise(
    p_value = suppressWarnings(
      if (length(unique(.data[[GROUP_COL]])) == 2L) wilcox.test(abundance ~ .data[[GROUP_COL]])$p.value else NA_real_
    ),
    ymax = max(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = vapply(p_value, fmt_p, character(1)))
top20_boxplot <- ggplot(top20_long, aes(.data[[GROUP_COL]], abundance, fill = .data[[GROUP_COL]])) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.25, color = "black") +
  facet_wrap(~ taxon, scales = "free_y", ncol = 5) +
  geom_text(data = top20_pvals, aes(x = 1.5, y = ymax, label = label), inherit.aes = FALSE, vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Top 20 species: sample-level relative abundance", x = "", y = "Relative abundance") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none")

species_presence <- (species_df[, -1, drop = FALSE] > 0)
case_ids <- meta[[SAMPLE_ID_COL]][meta[[GROUP_COL]] == CASE_LABEL]
ctrl_ids <- meta[[SAMPLE_ID_COL]][meta[[GROUP_COL]] == CTRL_LABEL]
case_mat <- species_presence[species_df$sample_id %in% case_ids, , drop = FALSE]
ctrl_mat <- species_presence[species_df$sample_id %in% ctrl_ids, , drop = FALSE]
species_in_case <- colnames(case_mat)[colSums(case_mat) > 0]
species_in_ctrl <- colnames(ctrl_mat)[colSums(ctrl_mat) > 0]
shared_species <- intersect(species_in_case, species_in_ctrl)
case_unique <- setdiff(species_in_case, species_in_ctrl)
ctrl_unique <- setdiff(species_in_ctrl, species_in_case)

venn_plot <- NULL
if (requireNamespace("ggvenn", quietly = TRUE)) {
  venn_plot <- ggvenn::ggvenn(list(Cases = species_in_case, Controls = species_in_ctrl), fill_color = unname(group_colors[c(CASE_LABEL, CTRL_LABEL)]), stroke_color = NA)
}

prev_tbl <- to_long(species_df) %>%
  group_by(taxon, .data[[GROUP_COL]]) %>%
  summarise(prevalence = sum(abundance > 0), mean_rel_abundance = mean(abundance, na.rm = TRUE), .groups = "drop")
present2 <- prev_tbl %>%
  group_by(taxon) %>%
  summarise(total_prev = sum(prevalence), .groups = "drop") %>%
  filter(total_prev >= 2) %>%
  pull(taxon)
prev_tbl <- prev_tbl %>% filter(taxon %in% present2)
spec_class <- data.frame(
  taxon = c(shared_species, case_unique, ctrl_unique),
  class = c(rep("Shared", length(shared_species)), rep(CASE_LABEL, length(case_unique)), rep(CTRL_LABEL, length(ctrl_unique))),
  stringsAsFactors = FALSE
)
prev_tbl <- left_join(prev_tbl, spec_class, by = "taxon")
prev_tbl$class[is.na(prev_tbl$class)] <- "Shared"
prev_scatter <- ggplot(prev_tbl, aes(prevalence, mean_rel_abundance, color = .data[[GROUP_COL]])) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Prevalence and mean relative abundance of species present in at least two individuals", x = "Prevalence", y = "Mean relative abundance") +
  theme_minimal(base_size = 10)

make_prev_abund_plot <- function(species_list, axis_label) {
  sp_data <- prev_tbl[prev_tbl$taxon %in% species_list, , drop = FALSE]
  if (!nrow(sp_data)) return(NULL)
  max_prev <- tapply(sp_data$prevalence, sp_data$taxon, max)
  keep_sp  <- names(max_prev)[max_prev >= 3L]
  if (!length(keep_sp)) keep_sp <- names(sort(max_prev, decreasing = TRUE))[seq_len(min(20L, length(max_prev)))]
  sp_data  <- sp_data[sp_data$taxon %in% keep_sp, , drop = FALSE]
  case_prev_ord <- sp_data[sp_data[[GROUP_COL]] == CASE_LABEL, , drop = FALSE]
  ord <- case_prev_ord$taxon[order(case_prev_ord$prevalence, decreasing = FALSE)]
  remaining <- setdiff(unique(sp_data$taxon), ord)
  ord <- c(remaining, ord)
  sp_data$taxon <- factor(sp_data$taxon, levels = ord)
  ggplot(sp_data, aes(x = taxon, fill = .data[[GROUP_COL]])) +
    geom_bar(aes(y = prevalence), stat = "identity",
             position = position_dodge(width = 0), alpha = 0.5, width = 0.7) +
    geom_point(aes(y = mean_rel_abundance * 100, color = .data[[GROUP_COL]]),
               position = position_dodge(width = 0.7), size = 3,
               show.legend = FALSE) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    scale_y_continuous(
      name = "Prevalence (bars)",
      sec.axis = sec_axis(~ ., name = "Mean relative abundance (%) (dots)")
    ) +
    labs(x = axis_label, fill = NULL) +
    coord_flip() +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 9),
          axis.title   = element_text(size = 11),
          legend.text  = element_text(size = 10),
          legend.title = element_text(size = 11))
}

case_prev_plot <- make_prev_abund_plot(case_unique, paste0(CASE_LABEL, "-specific species"))
ctrl_prev_plot <- make_prev_abund_plot(ctrl_unique, paste0(CTRL_LABEL, "-specific species"))
group_specific_plot <- case_prev_plot

dir.create(file.path(out_dir, "composition"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "flags"), recursive = TRUE, showWarnings = FALSE)

write.table(data.frame(species = shared_species), file.path(out_dir, "composition", "species_shared.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(species = case_unique), file.path(out_dir, "composition", "species_unique_cases.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(species = ctrl_unique), file.path(out_dir, "composition", "species_unique_controls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prev_tbl, file.path(out_dir, "composition", "species_prevalence_mean_relative_abundance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(group_mean_species, file.path(out_dir, "composition", "top20_species_group_mean.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top20_pvals, file.path(out_dir, "composition", "top20_species_wilcoxon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Generating composition plots...")
pdf(file.path(out_dir, "composition", "composition.pdf"), width = 12, height = 6)
print(plot_sample_facets(phylum_df, "Composition by Phylum"))
print(plot_sample_facets(genus_df, "Composition by Genus"))
print(plot_sample_facets(species_plot_df, "Composition by Species"))
print(top20_plot)
print(top20_boxplot)
if (!is.null(venn_plot)) print(venn_plot)
print(prev_scatter)
if (!is.null(group_specific_plot)) print(group_specific_plot)
if (!is.null(ctrl_prev_plot)) print(ctrl_prev_plot)
dev.off()

flags <- data.frame(
  step = "09_composition",
  key = c("n_samples", "n_phyla", "n_genera", "n_species", "n_shared_species", "n_unique_case_species", "n_unique_control_species"),
  value = c(nrow(tss), ncol(phylum_df) - 1L, ncol(genus_df) - 1L, ncol(species_df) - 1L, length(shared_species), length(case_unique), length(ctrl_unique)),
  status = "ok",
  stringsAsFactors = FALSE
)
write.table(flags, file.path(out_dir, "flags", "flags_09_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Composition done: ", nrow(tss), " samples")
