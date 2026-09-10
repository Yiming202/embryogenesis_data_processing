#!/usr/bin/env Rscript

# SHE1 versus AF16 3D position comparison.
# Temporal correction follows she1af16.R. The primary endpoint follows
# con_position.R; a companion coordinate plot uses centroid- and
# embryo-length-normalized coordinates in the common imaging-axis frame.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
})

base_dir <- path.expand("~/Downloads/NCRevision/she1_af16_normalized_concordance/R_exploratory_z_slope_0.95_1.00_six//")
source_dir <- file.path(base_dir, "source_data")
figure_dir <- file.path(base_dir, "figures")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Time limits and minutes per time point exactly follow she1af16.R.
# Both SHE1 and AF16 are treated with the C. briggsae z calibration used in
# the positional-normalization workflow. Confirm this is the correct voxel
# calibration for the AF16 acquisitions before manuscript submission.
embryos <- tribble(
  ~group, ~embryo, ~file, ~time_limit, ~minutes_per_tp, ~z_factor,
  "she1", "she1p1", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD240731cbhis72p1.csv", 165, 1.57, 4.78,
  "she1", "she1p2", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD240731cbhis72p2.csv", 175, 1.57, 4.78,
  "she1", "she1p3", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD240731cbhis72p3.csv", 170, 1.57, 4.78,
  "she1", "she1p4", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD241202cbhis72p1.csv", 160, 1.58, 4.78,
  "she1", "she1p5", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD241202cbhis72p2.csv", 165, 1.58, 4.78,
  "she1", "she1p6", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD241202cbhis72p4.csv", 180, 1.58, 4.78,
  "af16", "af16p1", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260412cbhis72p1.csv", 110, 1.640, 4.78,
  "af16", "af16p2", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260412cbhis72p3.csv", 115, 1.640, 4.78,
  "af16", "af16p3", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260413cbhis72p5.csv", 105, 1.800, 4.78,
  "af16", "af16p4", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260415cbhis72p1.csv", 160, 1.240, 4.78,
  "af16", "af16p5", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260417cbhis72p2.csv", 140, 1.300, 4.78,
  "af16", "af16p6", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260505cbhis72p2.csv", 125, 1.547, 4.78
)

lineage_file <- file.path(base_dir, "input_metadata", "AllLineage.tsv")
coefficient_file <- file.path(source_dir, "she1_af16_time_normalization_coefficients.csv")
stage_file <- file.path(base_dir, "input_metadata", "Length_Time_Absolute.csv")
founder_cells <- c("Zygote", "AB", "P1", "ABa", "ABp", "EMS", "P2")
quantiles <- c(.20, .40, .60, .80, 1.00)
reference_embryo <- "she1p4"
required_n <- c(she1 = 6L, af16 = 6L)

stopifnot(file.exists(lineage_file), file.exists(coefficient_file), file.exists(stage_file), all(file.exists(path.expand(embryos$file))))

read_positions <- function(group, embryo, file, time_limit, minutes_per_tp, z_factor) {
  d <- read_csv(path.expand(file), show_col_types = FALSE)
  needed <- c("cell", "time", "x", "y", "z")
  if (!all(needed %in% names(d))) stop("Missing position columns in ", file)
  d %>%
    transmute(
      group, embryo, cell = as.character(cell), time = as.numeric(time),
      time_min = as.numeric(time) * minutes_per_tp,
      x = as.numeric(x), y = as.numeric(y), z = as.numeric(z) * z_factor
    ) %>%
    filter(is.finite(time), is.finite(time_min), is.finite(x), is.finite(y), is.finite(z), time <= time_limit)
}

tracks_all <- pmap_dfr(embryos, read_positions)
ab2_cells <- c("ABa", "ABp", "EMS", "P2")
ab2_offsets <- tracks_all %>%
  filter(cell %in% ab2_cells) %>%
  group_by(group, embryo) %>%
  summarise(n_ab2_cells = n_distinct(cell), ab2_offset_min = max(time_min), .groups = "drop")
if (nrow(ab2_offsets) != nrow(embryos) || any(ab2_offsets$n_ab2_cells < length(ab2_cells))) {
  stop("All four AB2 cells are required to apply the she1af16.R temporal normalization.")
}
time_coefficients <- read_csv(coefficient_file, show_col_types = FALSE) %>%
  transmute(time_column = as.character(embryo), intercept = as.numeric(intercept), slope = as.numeric(slope), r_squared = as.numeric(r_squared)) %>%
  filter(time_column %in% paste0("Tim_", embryos$embryo))
if (nrow(time_coefficients) != nrow(embryos) || any(!is.finite(time_coefficients$slope)) || any(time_coefficients$slope <= 0)) {
  stop("The SHE1-AF16 time-normalization coefficient table is incomplete or invalid.")
}
tracks_all <- tracks_all %>%
  mutate(time_column = paste0("Tim_", embryo)) %>%
  left_join(ab2_offsets %>% select(group, embryo, ab2_offset_min), by = c("group", "embryo")) %>%
  left_join(time_coefficients, by = "time_column") %>%
  mutate(
    ab2_aligned_time_min = time_min - ab2_offset_min,
    normalized_time = (ab2_aligned_time_min - intercept) / slope
  ) %>%
  filter(is.finite(normalized_time)) %>%
  select(-time_column)
lineage <- read_tsv(lineage_file, show_col_types = FALSE) %>%
  transmute(cell = as.character(CellName), lineage_id = as.character(ID), fate = as.character(CellFate)) %>%
  filter(!is.na(cell), !is.na(lineage_id)) %>% distinct(cell, .keep_all = TRUE)
stage_map <- read_csv(stage_file, show_col_types = FALSE) %>%
  transmute(cell = as.character(cell), stage = as.character(stage)) %>% distinct(cell, .keep_all = TRUE)

# A retained parent must be non-death, present in every embryo, and have both
# direct daughters observed in every embryo after the configured editing limit.
daughter_map <- lineage %>%
  mutate(parent_id = substr(lineage_id, 1, nchar(lineage_id) - 1)) %>%
  filter(nchar(lineage_id) > 1) %>%
  group_by(parent_id) %>% summarise(daughters = list(cell), .groups = "drop")

observed <- tracks_all %>% distinct(embryo, cell)
cell_coverage <- observed %>% count(cell, name = "n_embryos")
candidate_cells <- lineage %>%
  filter(fate != "Death", !cell %in% founder_cells) %>%
  inner_join(cell_coverage, by = "cell") %>%
  filter(n_embryos == nrow(embryos)) %>%
  left_join(daughter_map, by = c("lineage_id" = "parent_id")) %>%
  mutate(daughters = map(daughters, ~ if (is.null(.x)) character() else .x)) %>%
  filter(map_int(daughters, length) == 2L)

has_two_daughters_in_all_embryos <- function(parent, daughters) {
  all(map_lgl(unique(embryos$embryo), function(e) {
    present <- observed %>% filter(embryo == e) %>% pull(cell)
    all(daughters %in% present)
  }))
}

retained_cells <- candidate_cells %>%
  mutate(daughters_in_all_embryos = map2_lgl(cell, daughters, has_two_daughters_in_all_embryos)) %>%
  filter(daughters_in_all_embryos) %>%
  inner_join(stage_map, by = "cell")
if (nrow(retained_cells) < 10L) stop("Fewer than ten shared, non-death cells with both daughters observed in all twelve embryos.")

tracks <- tracks_all %>%
  semi_join(retained_cells %>% select(cell), by = "cell") %>%
  inner_join(retained_cells %>% select(cell, stage), by = "cell") %>%
  group_by(group, embryo, cell) %>%
  mutate(time_percent = (normalized_time - min(normalized_time)) / (max(normalized_time) - min(normalized_time))) %>%
  ungroup() %>% filter(is.finite(time_percent))

# Same embryo-length estimator used in the prior revised spatial package.
embryo_length <- tracks_all %>%
  group_by(group, embryo, cell) %>% slice_max(time, n = 2, with_ties = FALSE) %>% ungroup() %>%
  group_by(group, embryo) %>% group_modify(~ {
    distances <- sort(as.numeric(dist(as.matrix(.x[, c("x", "y", "z")]))), decreasing = TRUE)
    tibble(embryo_length = median(distances[seq_len(min(10L, length(distances)))], na.rm = TRUE))
  }) %>% ungroup()

selected <- crossing(embryos %>% select(group, embryo), q = quantiles) %>%
  mutate(q_label = factor(paste0("q", q * 100), levels = paste0("q", quantiles * 100))) %>%
  left_join(tracks, by = c("group", "embryo"), relationship = "many-to-many") %>%
  group_by(group, embryo, q, q_label, cell) %>%
  slice_min(abs(time_percent - q), n = 1, with_ties = FALSE) %>%
  ungroup() %>% left_join(embryo_length, by = c("group", "embryo"))

# Mean-coordinate comparison requires a common positional origin. At each
# embryo/stage/time-quantile combination, use the centroid of retained cells
# as the origin and divide all three axes by that embryo's length. This keeps
# the imaging x/y/z axes while removing translation and uniform size effects.
selected_coordinates <- selected %>%
  group_by(group, embryo, stage, q, q_label) %>%
  mutate(
    x_normalized = (x - mean(x)) / embryo_length,
    y_normalized = (y - mean(y)) / embryo_length,
    z_normalized = (z - mean(z)) / embryo_length
  ) %>%
  ungroup()

# Imaging x/y orientation can differ across embryos. Rigidly register each
# centered, length-normalized configuration to the SHE1 reference embryo; no
# reflection or additional scaling is permitted. Pairwise distances are thus
# unchanged, while coordinate axes become directly comparable.
coordinate_reference <- selected_coordinates %>%
  filter(embryo == reference_embryo) %>%
  select(stage, q, cell, x_reference = x_normalized, y_reference = y_normalized, z_reference = z_normalized)

rigid_register <- function(d, ref) {
  matched <- d %>%
    select(cell, x_normalized, y_normalized, z_normalized) %>%
    inner_join(ref, by = "cell") %>%
    arrange(cell)
  if (nrow(matched) < 3L) {
    return(d %>% mutate(
      x_registered = NA_real_,
      y_registered = NA_real_,
      z_registered = NA_real_
    ))
  }
  x_mat <- as.matrix(matched[, c("x_normalized", "y_normalized", "z_normalized")])
  y_mat <- as.matrix(matched[, c("x_reference", "y_reference", "z_reference")])
  s <- svd(crossprod(x_mat, y_mat))
  rotation <- s$u %*% t(s$v)
  if (det(rotation) < 0) {
    s$u[, 3] <- -s$u[, 3]
    rotation <- s$u %*% t(s$v)
  }
  registered <- as.matrix(d[, c("x_normalized", "y_normalized", "z_normalized")]) %*% rotation
  d %>% mutate(
    x_registered = registered[, 1],
    y_registered = registered[, 2],
    z_registered = registered[, 3]
  )
}

selected_registered_coordinates <- selected_coordinates %>%
  group_by(group, embryo, stage, q, q_label) %>%
  group_modify(~ {
    ref <- coordinate_reference %>% filter(stage == .y$stage[[1]], q == .y$q[[1]])
    rigid_register(.x, ref)
  }) %>%
  ungroup()

coordinate_group_mean <- selected_registered_coordinates %>%
  select(group, embryo, stage, q, q_label, cell, x_registered, y_registered, z_registered) %>%
  pivot_longer(c(x_registered, y_registered, z_registered), names_to = "axis", values_to = "normalized_coordinate") %>%
  mutate(axis = sub("_registered$", "", axis)) %>%
  group_by(group, stage, q, q_label, cell, axis) %>%
  summarise(n_embryos = n_distinct(embryo), mean_normalized_coordinate = mean(normalized_coordinate), .groups = "drop") %>%
  filter(n_embryos == required_n[group])

coordinate_pairs <- coordinate_group_mean %>%
  select(group, stage, q, q_label, cell, axis, mean_normalized_coordinate) %>%
  pivot_wider(names_from = group, values_from = mean_normalized_coordinate) %>%
  filter(is.finite(she1), is.finite(af16)) %>%
  mutate(difference = af16 - she1, absolute_difference = abs(difference))

make_distance_table <- function(d) {
  cells <- sort(unique(d$cell))
  if (length(cells) < 3L) return(tibble())
  xyz <- d %>% arrange(cell) %>% select(x, y, z) %>% as.matrix()
  distance_matrix <- as.matrix(dist(xyz)) / d$embryo_length[1]
  pairs <- which(upper.tri(distance_matrix), arr.ind = TRUE)
  tibble(
    cell_1 = cells[pairs[, 1]],
    cell_2 = cells[pairs[, 2]],
    normalized_distance = distance_matrix[pairs]
  )
}

# This is the con_position.R endpoint: within-stage cell-cell distances divided
# by embryo length. It requires neither coordinate centering nor rotation.
distance_long <- selected %>%
  group_by(group, embryo, stage, q, q_label) %>%
  group_modify(~ make_distance_table(.x)) %>%
  ungroup()

distance_group_mean <- distance_long %>%
  group_by(group, stage, q, q_label, cell_1, cell_2) %>%
  summarise(n_embryos = n_distinct(embryo), normalized_distance_mean = mean(normalized_distance), .groups = "drop") %>%
  filter(n_embryos == required_n[group])

distance_pairs <- distance_group_mean %>%
  select(group, stage, q, q_label, cell_1, cell_2, normalized_distance_mean) %>%
  pivot_wider(names_from = group, values_from = normalized_distance_mean) %>%
  filter(is.finite(she1), is.finite(af16)) %>%
  mutate(difference = af16 - she1, absolute_difference = abs(difference))

safe_cor <- function(x, y, method = "pearson") if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) NA_real_ else cor(x, y, method = method)
agreement <- distance_pairs %>%
  summarise(
    n_pairs = n(), pearson_r = safe_cor(she1, af16), spearman_rho = safe_cor(she1, af16, "spearman"),
    slope = unname(coef(lm(af16 ~ she1))[2]), intercept = unname(coef(lm(af16 ~ she1))[1]),
    r_squared = summary(lm(af16 ~ she1))$r.squared, rmse = sqrt(mean(difference^2)), mae = mean(absolute_difference), mean_difference = mean(difference)
  )
agreement_by_stage_time <- distance_pairs %>%
  group_by(stage, q, q_label) %>%
  summarise(n_pairs = n(), pearson_r = safe_cor(she1, af16), slope = unname(coef(lm(af16 ~ she1))[2]), r_squared = summary(lm(af16 ~ she1))$r.squared, rmse = sqrt(mean(difference^2)), mae = mean(absolute_difference), .groups = "drop")

coordinate_agreement <- coordinate_pairs %>%
  group_by(axis) %>%
  summarise(
    n_coordinates = n(),
    pearson_r = safe_cor(she1, af16),
    pearson_p = cor.test(she1, af16)$p.value,
    spearman_rho = safe_cor(she1, af16, "spearman"),
    slope = unname(coef(lm(af16 ~ she1))[2]),
    intercept = unname(coef(lm(af16 ~ she1))[1]),
    r_squared = summary(lm(af16 ~ she1))$r.squared,
    rmse = sqrt(mean(difference^2)),
    mae = mean(absolute_difference),
    mean_difference = mean(difference),
    .groups = "drop"
  )

stage_order <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32", "Stage_AB64", "Stage_AB128", "Stage_AB256")
theme_set(theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 10), strip.text = element_text(face = "bold")))
label <- agreement %>% transmute(label = paste0("Pearson r = ", sprintf("%.3f", pearson_r), "\nR2 = ", sprintf("%.3f", r_squared), "\nSlope = ", sprintf("%.3f", slope), "\nRMSE = ", sprintf("%.4f", rmse), "\nn = ", n_pairs))

coordinate_labels <- coordinate_agreement %>%
  transmute(
    axis,
    label = paste0("r = ", sprintf("%.3f", pearson_r),"\nP = ", format(pearson_p, scientific = FALSE, digits = 3), "\nR2 = ", sprintf("%.3f", r_squared), "\nSlope = ", sprintf("%.3f", slope), "\nRMSE = ", sprintf("%.4f", rmse), "\nn = ", n_coordinates)
  )

p_coordinate_agreement <- ggplot(coordinate_pairs, aes(she1, af16)) +
  geom_abline(slope = 1,intercept = 0,colour = "black",linewidth = 1) +
  geom_point(size = 0.75,alpha = 0.55,colour = "#A8C6E4") +
 # geom_smooth(method = "lm", se = FALSE, colour = "#C95745", linewidth = .55) +
  geom_label(data = coordinate_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.04, size = 2.45, label.size = 0) +
  facet_wrap(~axis, nrow = 1) +
  coord_equal() +
  #facet_wrap(~axis, nrow = 1, scales = "free")+
  theme(
    panel.border = element_rect(colour = "black",fill = NA,linewidth = 0.45),axis.line = element_blank(),axis.ticks = element_line(linewidth = 0.35,colour = "black"),
    axis.text = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title = element_text(
      size = 9,
      colour = "black"
    ),
    strip.text = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5
    ),
    panel.grid = element_blank()
  ) +
  labs(title = "SHE1 versus AF16 mean normalized coordinates", subtitle = "Stage-centroided, embryo-length-normalized and rigidly registered; x, y, and z shown separately", x = "SHE1 mean normalized coordinate", y = "AF16 mean normalized coordinate")
print(p_coordinate_agreement)

# write_csv(embryos, file.path(source_dir, "embryo_metadata.csv"))
# write_csv(ab2_offsets, file.path(source_dir, "AB2_time_alignment_offsets.csv"))
# write_csv(time_coefficients, file.path(source_dir, "CE_template_time_normalization_coefficients.csv"))
# write_csv(retained_cells %>% select(cell, stage, fate, lineage_id), file.path(source_dir, "retained_shared_dividing_cells.csv"))
# write_csv(embryo_length, file.path(source_dir, "embryo_length.csv"))
# write_csv(selected, file.path(source_dir, "selected_positions_by_cell_quantile.csv"))
# write_csv(selected_coordinates, file.path(source_dir, "per_embryo_centroided_normalized_coordinates.csv"))
# write_csv(selected_registered_coordinates, file.path(source_dir, "per_embryo_registered_normalized_coordinates.csv"))
# write_csv(coordinate_pairs, file.path(source_dir, "SHE1_AF16_mean_normalized_coordinate_pairs.csv"))
# write_csv(coordinate_agreement, file.path(source_dir, "SHE1_AF16_normalized_coordinate_agreement_summary.csv"))
# write_csv(distance_long, file.path(source_dir, "per_embryo_normalized_pairwise_distances.csv"))
# write_csv(distance_pairs, file.path(source_dir, "SHE1_AF16_mean_normalized_distance_pairs.csv"))
# write_csv(agreement, file.path(source_dir, "SHE1_AF16_normalized_distance_agreement_summary.csv"))
# write_csv(agreement_by_stage_time, file.path(source_dir, "SHE1_AF16_normalized_distance_agreement_by_stage_quantile.csv"))

# ggsave(file.path(figure_dir, "Figure2_SHE1_AF16_normalized_pairwise_distance_agreement.pdf"), p_agreement, width = 6.5, height = 4.3)
# ggsave(file.path(figure_dir, "Figure2_SHE1_AF16_normalized_pairwise_distance_agreement.png"), p_agreement, width = 6.5, height = 4.3, dpi = 600)
# ggsave(file.path(figure_dir, "FigureS1_SHE1_AF16_normalized_pairwise_distance_differences.pdf"), p_difference, width = 10, height = 4.2)
# ggsave(file.path(figure_dir, "FigureS1_SHE1_AF16_normalized_pairwise_distance_differences.png"), p_difference, width = 10, height = 4.2, dpi = 600)
ggsave(file.path(figure_dir, "Figure3_SHE1_AF16_mean_normalized_xyz_agreement.pdf"), p_coordinate_agreement, width = 8.6, height = 3.5)
ggsave(file.path(figure_dir, "Figure3_SHE1_AF16_mean_normalized_xyz_agreement.png"), p_coordinate_agreement, width = 8.6, height = 3.5, dpi = 600)

writeLines(c(
  "# SHE1 versus AF16 normalized relative 3D arrangement", "",
  "## Input and cell eligibility", "The analysis uses the six SHE1 and six AF16 CD files and the editing limits/time resolutions specified in `she1af16.R`. Cells annotated as Death, founder cells, and parents without both direct daughters observed in every one of the 12 embryos are excluded.", "",
  "## Temporal and spatial normalization", "Temporal normalization exactly follows `she1af16.R`. For each embryo, the maximum terminal time of ABa, ABp, EMS, and P2 is subtracted from every frame. The CE-template regression correction is then applied: normalized time = (AB2-aligned time in minutes - embryo-specific intercept) / embryo-specific slope. Positions are selected at five within-cell quantiles of this normalized time. Following `con_position.R`, the Euclidean distance between every pair of cells within each stage is divided by that embryo's length. These unitless distances are invariant to microscope translation, embryo orientation, and uniform embryo-size differences.", "",
  "## Mean-coordinate companion analysis", "For the mean-coordinate comparison, the retained-cell centroid is subtracted independently for every embryo, developmental stage, and temporal quantile, and x, y, and z are divided by the embryo length. Each configuration is then rigidly registered to the SHE1 reference embryo (she1p4) using the shared cells at the same stage and temporal quantile. Registration permits rotation only; it does not permit reflection or additional scaling. Group means are calculated only when all six embryos contribute. This yields relative coordinates in a shared embryo frame; values can be negative because the stage-specific centroid is set to zero.", "",
  "## Group comparison", "Mean normalized pairwise distances and mean normalized coordinates are calculated only when all six SHE1 and all six AF16 embryos contribute. Pearson/Spearman correlation, slope, R2, RMSE, MAE, and mean difference are descriptive agreement measures. They are not treated as independent-cell hypothesis tests.", "",
  "## Calibration note", "The z correction of 4.78 follows the C. briggsae position-normalization workflow. Confirm that the AF16 data were acquired with the same z-step calibration before using the result in a manuscript."
), file.path(base_dir, "METHODS.md"))

writeLines(c("# Contents", "", "- `code/revised_she1_af16_normalized_xyz.R`: reproducible R analysis.", "- `figures/`: normalized pairwise-distance agreement, difference distribution, and mean x/y/z-coordinate agreement.", "- `source_data/`: every retained cell, per-embryo normalized position/distance, group mean, and agreement statistic.", "", "Run with:", "", "```bash", "Rscript code/revised_she1_af16_normalized_xyz.R", "```"), file.path(base_dir, "README.md"))
message("Completed: ", base_dir)

