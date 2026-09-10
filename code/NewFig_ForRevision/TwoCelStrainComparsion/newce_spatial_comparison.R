#!/usr/bin/env Rscript

# Spatial comparison of old CE and new CE embryos.
# Primary endpoint: embryo-length-scaled within-stage pairwise 3D distances.
# Secondary endpoint: Procrustes-registered, embryo-length-scaled mean coordinates.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- path.expand("~/Downloads/NCRevision/ceStrain/revised")
source_dir <- file.path(base_dir, "source_data")
figure_dir <- file.path(base_dir, "figures")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Values for the old CE embryos follow normalize_pos.R. The CE4 cutoff is 200,
# matching that normalization script. Verify the three new-CE temporal and
# z-calibration values against the original acquisition record before submission.
embryos <- tribble(
  ~group, ~embryo, ~file, ~time_limit, ~intercept, ~slope, ~z_factor,
  "ce", "ce1", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD191108plc1p1.csv", 205, 7.10, 1.00, 4.67,
  "ce", "ce2", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200109plc1p1.csv", 205, -9.57, 1.05, 4.67,
  "ce", "ce3", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200113plc1p3.csv", 195, -3.95, 1.00, 4.67,
  "ce", "ce4", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200113plc1p2.csv", 200, 7.04, 1.04, 4.67,
  "ce", "ce5", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200322plc1p2.csv", 195, 1.79, 0.95, 4.67,
  "ce", "ce6", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200323plc1p1.csv", 185, -7.77, 0.95, 4.67,
  "ce", "ce7", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200326plc1p3.csv", 220, 4.46, 1.04, 4.67,
  "ce", "ce8", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200326plc1p4.csv", 195, 0.89, 0.98, 4.67,
  "new_ce", "new_ce1", "~/Downloads/NCRevision/ceStrain/ce/CD260707his72p1.csv", 210, 2.63, 1.03, 4.78,
  "new_ce", "new_ce2", "~/Downloads/NCRevision/ceStrain/ce/CD260707his72p2.csv", 220, 20.00, 1.07, 4.78,
  "new_ce", "new_ce3", "~/Downloads/NCRevision/ceStrain/ce/CD260707his72p3.csv", 215, 7.01, 1.01, 4.78
)

stage_file <- path.expand("~/Desktop/HybridAnalysis/cbcn/Length_Time_Absolute.csv")
founder_cells <- c("Zygote", "AB", "P1", "ABa", "ABp", "EMS", "P2")
quantiles <- c(0.20, 0.40, 0.60, 0.80, 1.00)
reference_embryo <- "ce4"
required_group_n <- c(ce = 8L, new_ce = 3L)

if (!file.exists(stage_file)) stop("Stage annotation file not found: ", stage_file)
if (!all(file.exists(path.expand(embryos$file)))) stop("At least one configured CD file does not exist.")

standardize_cd <- function(group, embryo, file, time_limit, intercept, slope, z_factor) {
  d <- read_csv(path.expand(file), show_col_types = FALSE)
  needed <- c("cell", "time", "x", "y", "z")
  if (!all(needed %in% names(d))) stop("Missing cell/time/x/y/z columns in ", file)
  d %>%
    transmute(cell = as.character(cell), time = as.numeric(time), x = as.numeric(x), y = as.numeric(y), z = as.numeric(z) * z_factor) %>%
    filter(is.finite(time), is.finite(x), is.finite(y), is.finite(z), time <= time_limit) %>%
    mutate(group = group, embryo = embryo, nortime = (time - intercept) / slope) %>%
    select(group, embryo, cell, time, nortime, x, y, z)
}

tracks_all <- pmap_dfr(embryos, standardize_cd)
# cell_sets <- tracks_all %>% distinct(embryo, cell) %>% group_by(cell) %>% summarise(n_embryos = n_distinct(embryo), .groups = "drop")
# complete_cells <- cell_sets %>% filter(n_embryos == nrow(embryos), !cell %in% founder_cells) %>% pull(cell)
# if (length(complete_cells) < 10L) stop("Fewer than ten cells are shared by all 11 embryos after founder-cell exclusion.")


############################################################
## Identify cells with observed daughter cells
############################################################

lineage <- read_tsv(
  "AllLineage.tsv",
  show_col_types = FALSE
)

required_lineage_cols <- c("CellName", "ID", "CellFate")

if (!all(required_lineage_cols %in% names(lineage))) {
  stop(
    "AllLineage.tsv must contain: ",
    paste(required_lineage_cols, collapse = ", ")
  )
}

lineage_id <- lineage %>%
  select(CellName, ID, CellFate) %>%
  filter(
    !is.na(CellName),
    !is.na(ID)
  ) %>%
  distinct(CellName, .keep_all = TRUE)

death_cells <- lineage_id %>%
  filter(CellFate == "Death") %>%
  pull(CellName) %>%
  unique()

observed_cells <- tracks_all %>%
  pull(cell) %>%
  unique()

observed_lineage_id <- lineage_id %>%
  filter(CellName %in% observed_cells)

observed_ids <- observed_lineage_id$ID

dividing_cell_table <- observed_lineage_id %>%
  rowwise() %>%
  mutate(
    has_observed_daughter = any(
      nchar(observed_ids) == nchar(ID) + 1 &
        startsWith(observed_ids, ID)
    )
  ) %>%
  ungroup()

dividing_cells <- dividing_cell_table %>%
  filter(has_observed_daughter) %>%
  pull(CellName) %>%
  unique()

no_daughter_cells <- dividing_cell_table %>%
  filter(!has_observed_daughter) %>%
  pull(CellName) %>%
  unique()

cat("\n========== Dividing-cell filtering summary ==========\n")
cat("Observed cells:", length(observed_cells), "\n")
cat("Observed cells with lineage ID:", nrow(observed_lineage_id), "\n")
cat("Death cells:", length(death_cells), "\n")
cat("Founder cells:", length(founder_cells), "\n")
cat("Cells with observed daughter cells:", length(dividing_cells), "\n")
cat("Cells without observed daughter cells:", length(no_daughter_cells), "\n")
cat("====================================================\n\n")

############################################################
## Shared complete dividing cells
############################################################

cell_sets <- tracks_all %>%
  distinct(embryo, cell) %>%
  group_by(cell) %>%
  summarise(
    n_embryos = n_distinct(embryo),
    .groups = "drop"
  )

complete_cells <- cell_sets %>%
  filter(
    n_embryos == nrow(embryos),
    !cell %in% founder_cells,
    !cell %in% death_cells,
    cell %in% dividing_cells
  ) %>%
  pull(cell)

cat("Complete cells before dividing/death/founder filter:", sum(cell_sets$n_embryos == nrow(embryos)), "\n")
cat("Complete dividing cells retained:", length(complete_cells), "\n\n")

if (length(complete_cells) < 10L) {
  stop("Fewer than ten cells are shared by all 11 embryos after founder/death/dividing-cell filtering.")
}


stage_map <- read_csv(stage_file, show_col_types = FALSE) %>%
  transmute(cell = as.character(cell), stage = as.character(stage)) %>% distinct(cell, .keep_all = TRUE)

tracks <- tracks_all %>%
  filter(cell %in% complete_cells) %>%
  inner_join(stage_map, by = "cell") %>%
  group_by(group, embryo, cell) %>%
  mutate(time_percent = (nortime - min(nortime)) / (max(nortime) - min(nortime))) %>%
  ungroup() %>%
  filter(is.finite(time_percent))

# Same length estimator as con_position.R: median of the ten greatest pairwise
# distances among the final two observations of each tracked cell.
embryo_length <- tracks_all %>%
  group_by(group, embryo, cell) %>%
  slice_max(time, n = 2, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(group, embryo) %>%
  group_modify(~ {
    xyz <- as.matrix(.x[, c("x", "y", "z")])
    distances <- sort(as.numeric(dist(xyz)), decreasing = TRUE)
    tibble(embryo_length = median(distances[seq_len(min(10L, length(distances)))], na.rm = TRUE))
  }) %>%
  ungroup()
if (any(!is.finite(embryo_length$embryo_length)) || any(embryo_length$embryo_length <= 0)) stop("Invalid embryo-length estimate.")

# Select one position per cell at each within-cell temporal quantile, following
# the time-quantile logic in con_position.R.
selected_positions <- crossing(embryos %>% select(group, embryo), q = quantiles) %>%
  mutate(q_label = paste0("q", q * 100)) %>%
  left_join(tracks, by = c("group", "embryo"), relationship = "many-to-many") %>%
  group_by(group, embryo, q, q_label, cell) %>%
  slice_min(abs(time_percent - q), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(embryo_length, by = c("group", "embryo"))

coverage <- selected_positions %>% count(group, embryo, stage, q_label, name = "n_cells")
write_csv(embryos, file.path(source_dir, "embryo_metadata.csv"))
write_csv(cell_sets, file.path(source_dir, "all_embryo_cell_coverage.csv"))
write_csv(coverage, file.path(source_dir, "selected_position_coverage.csv"))
write_csv(embryo_length, file.path(source_dir, "embryo_length.csv"))

make_distance_table <- function(d) {
  cells <- sort(unique(d$cell))
  if (length(cells) < 3L) return(tibble())
  xyz <- d %>% arrange(cell) %>% select(x, y, z) %>% as.matrix()
  dm <- as.matrix(dist(xyz)) / d$embryo_length[1]
  pairs <- which(upper.tri(dm), arr.ind = TRUE)
  tibble(cell_1 = cells[pairs[, 1]], cell_2 = cells[pairs[, 2]], normalized_distance = dm[pairs])
}

distance_long <- selected_positions %>%
  group_by(group, embryo, stage, q, q_label) %>%
  group_modify(~ make_distance_table(.x)) %>%
  ungroup()

# Mean distance is calculated only when every embryo in the relevant group
# contributed the same cell pair.
distance_group_mean <- distance_long %>%
  group_by(group, stage, q, q_label, cell_1, cell_2) %>%
  summarise(n_embryos = n_distinct(embryo), normalized_distance_mean = mean(normalized_distance), .groups = "drop") %>%
  filter(n_embryos == required_group_n[group])

distance_pair <- distance_group_mean %>%
  select(group, stage, q, q_label, cell_1, cell_2, normalized_distance_mean) %>%
  pivot_wider(names_from = group, values_from = normalized_distance_mean) %>%
  filter(is.finite(ce), is.finite(new_ce)) %>%
  mutate(difference = new_ce - ce, absolute_difference = abs(difference))

safe_cor <- function(x, y, method = "pearson") {
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = method)
}

distance_agreement <- distance_pair %>%
  group_by(stage, q, q_label) %>%
  summarise(
    n_pairs = n(),
    pearson_r = safe_cor(ce, new_ce),
    spearman_rho = safe_cor(ce, new_ce, "spearman"),
    slope = unname(coef(lm(new_ce ~ ce))[2]),
    intercept = unname(coef(lm(new_ce ~ ce))[1]),
    r_squared = summary(lm(new_ce ~ ce))$r.squared,
    rmse = sqrt(mean(difference^2)),
    mae = mean(absolute_difference),
    mean_difference = mean(difference),
    .groups = "drop"
  )

write_csv(distance_long, file.path(source_dir, "per_embryo_normalized_pairwise_distances.csv"))
write_csv(distance_pair, file.path(source_dir, "CE_newCE_mean_normalized_distance_pairs.csv"))
write_csv(distance_agreement, file.path(source_dir, "CE_newCE_distance_agreement_by_stage_quantile.csv"))

# Secondary coordinate analysis. Each within-stage configuration is centered,
# scaled by embryo length, then rotated to the CE4 reference configuration by
# an orientation-preserving Procrustes transform. Raw coordinates are never
# compared across embryos.
reference_configs <- selected_positions %>%
  filter(embryo == reference_embryo) %>%
  select(stage, q, q_label, cell, x_ref = x, y_ref = y, z_ref = z, ref_length = embryo_length)

register_config <- function(d, reference) {
  common <- inner_join(d %>% select(cell, x, y, z, embryo_length), reference, by = "cell") %>% arrange(cell)
  if (nrow(common) < 3L) return(tibble())
  target <- as.matrix(common[, c("x", "y", "z")]) / common$embryo_length[1]
  ref <- as.matrix(common[, c("x_ref", "y_ref", "z_ref")]) / common$ref_length[1]
  target <- sweep(target, 2, colMeans(target), "-")
  ref <- sweep(ref, 2, colMeans(ref), "-")
  s <- svd(t(target) %*% ref)
  rotation <- s$u %*% t(s$v)
  if (det(rotation) < 0) {
    s$u[, 3] <- -s$u[, 3]
    rotation <- s$u %*% t(s$v)
  }
  aligned <- target %*% rotation
  tibble(cell = common$cell, x = aligned[, 1], y = aligned[, 2], z = aligned[, 3])
}

registered_coordinates <- selected_positions %>%
  group_by(group, embryo, stage, q, q_label) %>%
  group_modify(~ {
    ref <- reference_configs %>% filter(stage == .y$stage, q == .y$q) %>% select(cell, x_ref, y_ref, z_ref, ref_length)
    register_config(.x, ref)
  }) %>%
  ungroup()

coordinate_group_mean <- registered_coordinates %>%
  group_by(group, stage, q, q_label, cell) %>%
  summarise(n_embryos = n_distinct(embryo), x = mean(x), y = mean(y), z = mean(z), .groups = "drop") %>%
  filter(n_embryos == required_group_n[group])

coordinate_pair <- coordinate_group_mean %>%
  pivot_longer(c(x, y, z), names_to = "axis", values_to = "coordinate") %>%
  select(group, stage, q, q_label, cell, axis, coordinate) %>%
  pivot_wider(names_from = group, values_from = coordinate) %>%
  filter(is.finite(ce), is.finite(new_ce)) %>%
  mutate(difference = new_ce - ce)

coordinate_agreement <- coordinate_pair %>%
  group_by(axis) %>%
  summarise(
    n_points = n(), pearson_r = safe_cor(ce, new_ce), spearman_rho = safe_cor(ce, new_ce, "spearman"),
    slope = unname(coef(lm(new_ce ~ ce))[2]), r_squared = summary(lm(new_ce ~ ce))$r.squared,
    pearson_p = cor.test(ce, new_ce, method = "pearson")$p.value,
    rmse = sqrt(mean(difference^2)), mae = mean(abs(difference)), .groups = "drop"
  )

write_csv(registered_coordinates, file.path(source_dir, "per_embryo_Procrustes_registered_coordinates.csv"))
write_csv(coordinate_pair, file.path(source_dir, "CE_newCE_registered_mean_coordinate_pairs.csv"))
write_csv(coordinate_agreement, file.path(source_dir, "CE_newCE_registered_coordinate_agreement.csv"))

stage_order <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32", "Stage_AB64", "Stage_AB128", "Stage_AB256")
theme_set(theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 10), strip.text = element_text(face = "bold"), legend.position = "top"))

overall_distance <- distance_pair %>% summarise(
  n_pairs = n(), pearson_r = safe_cor(ce, new_ce), spearman_rho = safe_cor(ce, new_ce, "spearman"),
  slope = unname(coef(lm(new_ce ~ ce))[2]), r_squared = summary(lm(new_ce ~ ce))$r.squared,
  rmse = sqrt(mean(difference^2)), mae = mean(absolute_difference)
)
write_csv(overall_distance, file.path(source_dir, "CE_newCE_overall_distance_agreement.csv"))

distance_label <- overall_distance %>% transmute(label = paste0("Pearson r = ", sprintf("%.3f", pearson_r), "\nR2 = ", sprintf("%.3f", r_squared), "\nSlope = ", sprintf("%.3f", slope), "\nRMSE = ", sprintf("%.4f", rmse), "\nn = ", n_pairs))


coord_label <- coordinate_agreement %>% mutate(label = paste0("r = ", sprintf("%.3f", pearson_r), "\nR2 = ", sprintf("%.3f", r_squared), "\nRMSE = ", sprintf("%.4f", rmse)))
coord_label <- coordinate_agreement %>%
  mutate(
    p_label = paste0(
      "P = ",
      format(
        pearson_p,
        scientific = FALSE,
        digits = 4
      )
    ),
    label = paste0(
      "n = ", n_points,
      "\nPearson r = ", sprintf("%.3f", pearson_r),
      "\n", p_label,
      "\nR² = ", sprintf("%.3f", r_squared),
      "\nSlope = ", sprintf("%.3f", slope),
      "\nRMSE = ", sprintf("%.4f", rmse)
    )
  )

p_coordinate <- ggplot(coordinate_pair, aes(ce, new_ce)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    size = 0.75,
    alpha = 0.55,
    colour = "#A8C6E4"
  ) +
 # geom_smooth(method = "lm", se = FALSE, colour = "#C95745", linewidth = .55) +
  geom_label(data = coord_label, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.04, size = 2.5, label.size = 0) +
  facet_wrap(~axis, nrow = 1) +
  coord_equal() +
  theme_classic(
    base_size = 8,
    base_family = "sans"
  ) +
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
  labs(title = "Registered mean-coordinate agreement", subtitle = "Secondary analysis: centroided, embryo-length-scaled, orientation-preserving Procrustes registration", x = "CE registered coordinate", y = "new CE registered coordinate")

print(p_coordinate)

#ggsave(file.path(figure_dir, "FigureS1_registered_mean_coordinate_agreement.pdf"), p_coordinate, width = 7.2, height = 3.6)

writeLines(c(
  "# Revised CE versus new-CE spatial comparison", "",
  "## Primary comparison", "For each embryo and within-stage temporal quantile, coordinates of every cell pair are converted to Euclidean distances and divided by that embryo's length. Group means use only pairs represented by all eight CE embryos and all three new-CE embryos. This endpoint is invariant to image translation, embryo orientation, and overall embryo size.", "",
  "## Temporal selection", "Each cell track is transformed to normalized time using the configured embryo-specific intercept and slope. One row closest to each within-cell temporal quantile (20%, 40%, 60%, 80%, and 100%) is selected, following the reference positional-analysis logic. A single old-CE reference embryo is not needed for the primary distance endpoint because pairwise distances are coordinate-system invariant.", "",
  "## Secondary mean-coordinate analysis", "Mean x/y/z coordinates are reported only after each within-stage configuration has been centered, divided by embryo length, and aligned to the CE4 configuration using an orientation-preserving orthogonal Procrustes transform. These coordinates are unitless registered coordinates, not raw microscope coordinates.", "",
  "## Interpretation", "Pearson/Spearman correlation, slope, RMSE, and MAE are descriptive agreement measures. Correlation is not treated as independent-cell hypothesis testing. The three new-CE embryos are biological replicates; the source tables retain their contribution counts.", "",
  "## Calibration note", "The old-CE values come from the reference positional-normalization workflow. The new-CE timing and z-calibration parameters are taken from the supplied comparison script and must be checked against acquisition metadata before manuscript submission."
), file.path(base_dir, "METHODS.md"))

writeLines(c(
  "# Contents", "", "- `code/revised_ce_newce_spatial_comparison.R`: reproducible analysis.", "- `source_data/`: input configuration, embryo lengths, all pairwise distances, group means, and agreement statistics.", "- `figures/`: primary distance-based comparison and secondary registered-coordinate comparison.", "", "Run with:", "", "```bash", "Rscript code/revised_ce_newce_spatial_comparison.R", "```"
), file.path(base_dir, "README.md"))

message("Completed revised spatial comparison: ", base_dir)

