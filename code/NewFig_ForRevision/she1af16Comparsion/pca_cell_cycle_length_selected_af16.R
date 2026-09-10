#!/usr/bin/env Rscript

# PCA of cell-cycle length using ceCompare.R CE/SHE1/C. nigoni settings and
# the six AF16 embryos specified for this exploratory package.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2); library(ggrepel); library(tibble); library(stringr)
})

base_dir <- path.expand("~/Downloads/NCRevision/she1_af16_normalized_concordance/R_exploratory_z_slope_0.95_1.00_six/")
output_dir <- file.path(base_dir, "pca_cell_cycle_length")
figure_dir <- file.path(output_dir, "figures")
source_dir <- file.path(output_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

datasets_config <- tribble(
  ~group, ~id, ~file_path, ~time_limit, ~minutes_per_tp,
  "ce", 1, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD191108plc1p1.csv", 205, 1.43,
  "ce", 2, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200109plc1p1.csv", 205, 1.43,
  "ce", 3, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200113plc1p3.csv", 195, 1.44,
  "ce", 4, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200113plc1p2.csv", 205, 1.44,
  "ce", 5, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200322plc1p2.csv", 195, 1.44,
  "ce", 6, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200323plc1p1.csv", 185, 1.44,
  "ce", 7, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200326plc1p3.csv", 220, 1.44,
  "ce", 8, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200326plc1p4.csv", 195, 1.44,
  "new_ce", 1, "~/Downloads/NCRevision/ceStrain/ce/CD260707his72p1.csv", 210, 1.38,
  "new_ce", 2, "~/Downloads/NCRevision/ceStrain/ce/CD260707his72p2.csv", 220, 1.38,
  "new_ce", 3, "~/Downloads/NCRevision/ceStrain/ce/CD260707his72p3.csv", 215, 1.38,
  "she1", 1, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD240731cbhis72p1.csv", 165, 1.57,
  "she1", 2, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD240731cbhis72p2.csv", 175, 1.57,
  "she1", 3, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD240731cbhis72p3.csv", 170, 1.57,
  "she1", 4, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD241202cbhis72p1.csv", 160, 1.58,
  "she1", 5, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD241202cbhis72p2.csv", 165, 1.58,
  "she1", 6, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD241202cbhis72p4.csv", 180, 1.58,
  "cni", 1, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD241202cnhis72p1.csv", 235, 1.58,
  "cni", 2, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD240712cnhis72p1.csv", 235, 1.60,
  "cni", 3, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD240712cnhis72p2.csv", 235, 1.60,
  "cni", 4, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD240712cnhis72p3.csv", 235, 1.60,
  "cni", 5, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD241207cnhis72p1.csv", 230, 1.65,
  "cni", 6, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD241207cnhis72p3.csv", 230, 1.65,
  "cni", 7, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD241202cnhis72p2.csv", 200, 1.58,
  "af16", 1, "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260412cbhis72p1.csv", 110, 1.640,
  "af16", 2, "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260412cbhis72p3.csv", 115, 1.640,
  "af16", 3, "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260413cbhis72p5.csv", 105, 1.800,
  "af16", 4, "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260415cbhis72p1.csv", 160, 1.240,
  "af16", 5, "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260417cbhis72p2.csv", 140, 1.300,
  "af16", 6, "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260505cbhis72p2.csv", 125, 1.547
) %>% mutate(embryo = paste0(group, id))
stopifnot(all(file.exists(path.expand(datasets_config$file_path))))

process_one <- function(group, id, file_path, time_limit, minutes_per_tp, embryo) {
  read_csv(path.expand(file_path), show_col_types = FALSE) %>%
    transmute(cell = as.character(cell), time = as.numeric(time)) %>%
    filter(is.finite(time), time <= time_limit) %>%
    group_by(cell) %>%
    summarise(cell_cycle_length_min = (max(time) - min(time)) * minutes_per_tp, .groups = "drop") %>%
    mutate(group, embryo)
}

length_long <- pmap_dfr(datasets_config, process_one)
lineage <- read_tsv(file.path(base_dir, "input_metadata", "AllLineage.tsv"), show_col_types = FALSE) %>%
  transmute(cell = as.character(CellName), lineage_id = as.character(ID), fate = as.character(CellFate)) %>%
  filter(!is.na(cell), !is.na(lineage_id)) %>% distinct(cell, .keep_all = TRUE)
observed_ids <- lineage %>% filter(cell %in% unique(length_long$cell)) %>% pull(lineage_id)
dividing_cells <- lineage %>%
  rowwise() %>%
  mutate(has_observed_daughter = any(nchar(observed_ids) == nchar(lineage_id) + 1L & startsWith(observed_ids, lineage_id))) %>%
  ungroup() %>%
  filter(fate != "Death", !cell %in% c("Zygote", "AB", "P1", "ABa", "ABp", "EMS", "P2"), has_observed_daughter) %>%
  pull(cell)

length_filtered <- length_long %>% filter(cell %in% dividing_cells)
wide <- length_filtered %>%
  select(cell, embryo, cell_cycle_length_min) %>%
  pivot_wider(
    names_from = embryo,
    values_from = cell_cycle_length_min
  ) %>%
  drop_na()
complete <- wide %>% drop_na() %>% column_to_rownames("cell") %>% as.matrix()
if (nrow(complete) < 2) stop("Fewer than two shared dividing cells remain for PCA.")

pca <- prcomp(t(complete), center = TRUE, scale. = TRUE)
variance <- summary(pca)$importance[2, 1:2] * 100
scores <- as.data.frame(pca$x[, 1:2]) %>% rownames_to_column("embryo") %>%
  left_join(datasets_config %>% select(embryo, group), by = "embryo")
centroids <- scores %>% group_by(group) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
ellipse_scores <- scores %>% add_count(group, name = "n_group") %>% filter(n_group >= 4L)

group_labels <- c(ce = "C. elegans", new_ce = "New C. elegans", she1 = "C. briggsae SHE1", af16 = "C. briggsae AF16", cni = "C. nigoni")
group_colours <- c(ce = "#C95745", new_ce = "#923B4A", she1 = "#3A7CA5", af16 = "#2E8B57", cni = "#C69214")
p <- ggplot(scores, aes(PC1, PC2, colour = group)) +
  stat_ellipse(data = ellipse_scores, aes(group = group), type = "norm", level = .68, linewidth = .45, linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 2.6, alpha = .88) +
  geom_point(data = centroids, shape = 4, size = 4.2, stroke = 1.1, show.legend = FALSE) +
  geom_text_repel(data = centroids, aes(PC1, PC2, label = unname(group_labels[group])), colour = "black", inherit.aes = FALSE, size = 3, fontface = "bold", box.padding = .28, segment.colour = "grey45") +
  scale_colour_manual(values = group_colours, labels = group_labels, name = NULL) +
  labs(title = "PCA of cell-cycle-length profiles", subtitle = paste0(nrow(complete), " shared non-death dividing cells; selected AF16 six-embryo set"), x = sprintf("PC1 (%.1f%%)", variance[[1]]), y = sprintf("PC2 (%.1f%%)", variance[[2]])) +
  theme_classic(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8), legend.position = "bottom")
print(p)

write_csv(datasets_config, file.path(source_dir, "PCA_embryo_metadata.csv"))
write_csv(tibble(cell = rownames(complete)), file.path(source_dir, "PCA_shared_dividing_cells.csv"))
write_csv(as.data.frame(complete) %>% rownames_to_column("cell"), file.path(source_dir, "PCA_cell_cycle_length_matrix.csv"))
write_csv(scores, file.path(source_dir, "PCA_scores_cell_cycle_length.csv"))
write_csv(centroids, file.path(source_dir, "PCA_group_centroids.csv"))
write_csv(tibble(component = c("PC1", "PC2"), variance_explained_percent = variance), file.path(source_dir, "PCA_variance_explained.csv"))
ggsave(file.path(figure_dir, "Figure_PCA_cell_cycle_length_selected_AF16.pdf"), p, width = 6.7, height = 5.1)
ggsave(file.path(figure_dir, "Figure_PCA_cell_cycle_length_selected_AF16.png"), p, width = 6.7, height = 5.1, dpi = 600)
message("Completed PCA: ", output_dir)

