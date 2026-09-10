#read file
setwd("~/Downloads/NCRevision/enhanced_raw/")
getwd()
library(ggplot2);library(dplyr);library(hrbrthemes);library(viridis);
library(tibble);library(tidyverse);library(ggnewscale);library(tidyr)

raw <- read.csv("~/Downloads/NCRevision/enhanced_raw/raw_CD191108plc1p1.csv", header = TRUE)
enhanced <- read.csv("~/Downloads/NCRevision/enhanced_raw/enhanced_CD191108plc1p1.csv", header = TRUE)
death <- read_tsv("~/Downloads/NCRevision/enhanced_raw/AllLineage.tsv")

apoptotic_cells <- death %>%
  filter(CellFate == "Death") %>%
  pull(CellName)

#specify cell order
raw_time <-raw %>%
  filter(!time>205) %>%
  group_by(cell) %>%
  select(cell,time)%>%
  filter(time == max(time))

enhanced_time <-enhanced %>%
  filter(!time>205) %>%
  group_by(cell) %>%
  select(cell,time)%>%
  filter(time == max(time))

# raw_Len <- raw %>%
#   as.tibble() %>%
#   filter(!time>205) %>%
#   group_by(cell) %>%
#   tally() %>% #统计值
#   as.tibble() 
# enhanced_Len <- enhanced %>%
#   as.tibble() %>%
#   filter(!time>205) %>%
#   group_by(cell) %>%
#   tally() %>% #统计值
#   as.tibble() 

time_pair <- inner_join(raw_time, enhanced_time, by = "cell") %>%
  rename(
    raw_time = time.x,
    enhanced_time = time.y
  ) %>%
  filter(!cell %in% apoptotic_cells)

time_pair %>%
  filter(abs(raw_time - enhanced_time) >= 2)


fit <- lm(enhanced_time ~ raw_time, data = time_pair)

pearson_test <- cor.test(
  time_pair$raw_time,
  time_pair$enhanced_time,
  method = "pearson"
)

pearson_r <- unname(pearson_test$estimate)
pearson_p <- pearson_test$p.value

slope <- unname(coef(fit)[2])
intercept <- unname(coef(fit)[1])
r_squared <- summary(fit)$r.squared

label_text <- sprintf(
  "n = %d\nPearson r = %.10f\nP = %.10f\nSlope = %.6f\nIntercept = %.6f\nR² = %.10f",
  nrow(time_pair),
  pearson_r,
  pearson_p,
  slope,
  intercept,
  r_squared
)


p <- ggplot(time_pair, aes(x = raw_time, y = enhanced_time)) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "gray",
    #linetype = "dashed",
    linewidth = 1
  ) +
  geom_point(
    size = 2.0,
    #alpha = 0.3,
    colour = "#A8C6E4"
  ) +
  # geom_abline(
  #   slope = 1,
  #   intercept = 0,
  #   linetype = "dashed",
  #   colour = "#8DB3D0",
  #   linewidth = 0.6
  # ) +
  annotate(
    "label",
    x = Inf,
    y = -Inf,
    label = label_text,
    hjust = 1.04,
    vjust = -0.1,
    size = 6.0,
    lineheight = 0.94,
    linewidth = 0,
    fill = "white",
    alpha = 0.85
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    panel.border = element_rect(colour = "black",fill = NA,linewidth = 0.45),
    #axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    plot.title = element_text(size = 15, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Raw observed terminal time, frames",
    y = "Enhanced observed terminal time, frames",
    title = "Raw vs enhanced observed terminal time"
  )

print(p)

# ggsave(
#   "raw_vs_enhanced_terminal_time_pearson_lm.pdf",
#   p,
#   width = 6,
#   height = 6
# )






#!/usr/bin/env Rscript

# Compare raw and enhanced xyz coordinates after excluding Death cells.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# -----------------------------
# Paths
# -----------------------------

base_dir <- path.expand("~/Downloads/NCRevision/enhanced_raw")

raw_file <- file.path(base_dir, "raw_CD191108plc1p1.csv")
enhanced_file <- file.path(base_dir, "enhanced_CD191108plc1p1.csv")
lineage_file <- file.path(base_dir, "AllLineage.tsv")

out_dir <- file.path(base_dir, "xyz_raw_vs_enhanced")
figure_dir <- file.path(out_dir, "figures")
source_dir <- file.path(out_dir, "source_data")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

end_frame <- 205L

# -----------------------------
# Read data
# -----------------------------

raw <- read.csv(raw_file, stringsAsFactors = FALSE, check.names = FALSE)
enhanced <- read.csv(enhanced_file, stringsAsFactors = FALSE, check.names = FALSE)
lineage <- read_tsv(lineage_file, show_col_types = FALSE)

required_cols <- c("cellTime", "cell", "time", "x", "y", "z")

if (!all(required_cols %in% names(raw))) {
  stop("Raw file must contain: cellTime, cell, time, x, y, z")
}

if (!all(required_cols %in% names(enhanced))) {
  stop("Enhanced file must contain: cellTime, cell, time, x, y, z")
}

if (!all(c("CellName", "CellFate") %in% names(lineage))) {
  stop("AllLineage.tsv must contain CellName and CellFate.")
}

# -----------------------------
# Death-cell exclusion
# -----------------------------

death_cells <- lineage %>%
  filter(CellFate == "Death") %>%
  pull(CellName)

# -----------------------------
# Match raw and enhanced xyz by cellTime
# -----------------------------

xyz_pair <- inner_join(
  raw %>%
    filter(time <= end_frame) %>%
    select(cellTime, cell, time, x, y, z),
  enhanced %>%
    filter(time <= end_frame) %>%
    select(cellTime, cell, time, x, y, z),
  by = "cellTime",
  suffix = c("_raw", "_enhanced")
) %>%
  filter(cell_raw == cell_enhanced) %>%
  filter(time_raw == time_enhanced) %>%
  transmute(
    cell = cell_raw,
    time = time_raw,
    x_raw = x_raw,
    x_enhanced = x_enhanced,
    y_raw = y_raw,
    y_enhanced = y_enhanced,
    z_raw = z_raw,
    z_enhanced = z_enhanced
  ) %>%
  filter(!cell %in% death_cells) %>%
  filter(if_all(
    c(x_raw, x_enhanced, y_raw, y_enhanced, z_raw, z_enhanced),
    is.finite
  ))

# -----------------------------
# Convert to long format
# -----------------------------

xyz_long <- xyz_pair %>%
  pivot_longer(
    cols = c(
      x_raw, x_enhanced,
      y_raw, y_enhanced,
      z_raw, z_enhanced
    ),
    names_to = c("axis", ".value"),
    names_pattern = "([xyz])_(raw|enhanced)"
  )

# -----------------------------
# Summary statistics
# -----------------------------

agreement_summary <- xyz_long %>%
  group_by(axis) %>%
  summarise(
    n = n(),
    n_cells = n_distinct(cell),
    pearson_r = cor(raw, enhanced, method = "pearson"),
    pearson_p = cor.test(raw, enhanced, method = "pearson")$p.value,
    spearman_rho = cor(raw, enhanced, method = "spearman"),
    slope = unname(coef(lm(enhanced ~ raw))[2]),
    intercept = unname(coef(lm(enhanced ~ raw))[1]),
    r_squared = summary(lm(enhanced ~ raw))$r.squared,
    median_difference = median(enhanced - raw),
    mean_difference = mean(enhanced - raw),
    rmse = sqrt(mean((enhanced - raw)^2)),
    mae = mean(abs(enhanced - raw)),
    exact_percent = 100 * mean(raw == enhanced),
    within_1_unit_percent = 100 * mean(abs(enhanced - raw) <= 1),
    within_2_units_percent = 100 * mean(abs(enhanced - raw) <= 2),
    within_3_units_percent = 100 * mean(abs(enhanced - raw) <= 3),
    .groups = "drop"
  )

print(agreement_summary)

# -----------------------------
# Labels for plot
# -----------------------------

axis_labels <- agreement_summary %>%
  mutate(
    label = sprintf(
      "n = %d\nCells = %d\nPearson r = %.10f\nP = %.10f\nSlope = %.6f\nIntercept = %.6f\nR² = %.10f\nRMSE = %.3f\nMAE = %.3f",
      n,
      n_cells,
      pearson_r,
      pearson_p,
      slope,
      intercept,
      r_squared,
      rmse,
      mae
    )
  )
# -----------------------------
# Plot
# -----------------------------

p_xyz <- ggplot(xyz_long, aes(x = raw, y = enhanced)) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "gray30",
    fill = "gray75",
    linewidth = 0.7
  ) +
  geom_point(
    size = 0.45,
    alpha = 0.25,
    colour = "#A8C6E4"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "#4DBBD5",
    linewidth = 0.5
  ) +
  geom_label(
    data = axis_labels,
    aes(x = Inf, y = -Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.04,
    vjust = -0.1,
    size = 2.4,
    lineheight = 0.92,
    linewidth = 0,
    fill = "white",
    alpha = 0.85
  ) +
  facet_wrap(
    ~axis,
    scales = "free",
    ncol = 3
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    panel.border = element_rect(colour = "black",fill = NA,linewidth = 0.45),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    strip.text = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Raw coordinate",
    y = "Enhanced coordinate",
    title = "Raw vs enhanced xyz coordinate agreement"
  )

print(p_xyz)

# -----------------------------
# Save figure
# -----------------------------

ggsave(
  file.path(figure_dir, "raw_vs_enhanced_xyz_coordinate_agreement.pdf"),
  p_xyz,
  width = 9,
  height = 3.5
)

message("XYZ comparison written to: ", out_dir)




















#!/usr/bin/env Rscript

# =============================================================================
# Compare raw and enhanced xyz coordinates
#
# Filters:
#   1. Only frames <= end_frame
#   2. Exclude Death cells
#   3. Only keep cells with observed daughter cells
#      within the current raw/enhanced dataset
#
# A cell is considered to have observed daughter cells if another observed cell
# has an ID that is exactly one character longer and starts with the current ID.
#
# Example:
#   parent ID   = Zaaa
#   daughter ID = Zaaaa or Zaaap
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# =============================================================================
# Paths
# =============================================================================

base_dir <- path.expand("~/Downloads/NCRevision/enhanced_raw")

raw_file <- file.path(base_dir, "raw_CD191108plc1p1.csv")
enhanced_file <- file.path(base_dir, "enhanced_CD191108plc1p1.csv")
lineage_file <- file.path(base_dir, "AllLineage.tsv")

out_dir <- file.path(base_dir, "xyz_raw_vs_enhanced_dividing_cells_only")
figure_dir <- file.path(out_dir, "figures")
source_dir <- file.path(out_dir, "source_data")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

end_frame <- 205L

# =============================================================================
# Read data
# =============================================================================

raw <- read.csv(
  raw_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

enhanced <- read.csv(
  enhanced_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lineage <- read_tsv(
  lineage_file,
  show_col_types = FALSE
)

# =============================================================================
# Check required columns
# =============================================================================

required_xyz_cols <- c("cellTime", "cell", "time", "x", "y", "z")

if (!all(required_xyz_cols %in% names(raw))) {
  stop(
    "Raw file must contain: ",
    paste(required_xyz_cols, collapse = ", ")
  )
}

if (!all(required_xyz_cols %in% names(enhanced))) {
  stop(
    "Enhanced file must contain: ",
    paste(required_xyz_cols, collapse = ", ")
  )
}

required_lineage_cols <- c("CellName", "ID", "CellFate")

if (!all(required_lineage_cols %in% names(lineage))) {
  stop(
    "AllLineage.tsv must contain: ",
    paste(required_lineage_cols, collapse = ", ")
  )
}

# =============================================================================
# Death-cell exclusion
# =============================================================================

death_cells <- lineage %>%
  filter(CellFate == "Death") %>%
  pull(CellName) %>%
  unique()

# =============================================================================
# Identify cells with observed daughter cells
# =============================================================================
# This uses only cells that are observed in BOTH raw and enhanced data
# within end_frame.
#
# A cell is treated as dividing if:
#   another observed cell has ID length = parent ID length + 1
#   and starts with the parent ID.
#
# This matches the logic used in your lineage-tree plotting code.

lineage_id <- lineage %>%
  select(CellName, ID, CellFate) %>%
  filter(
    !is.na(CellName),
    !is.na(ID)
  ) %>%
  distinct(CellName, .keep_all = TRUE)

observed_cells_raw <- raw %>%
  filter(time <= end_frame) %>%
  pull(cell) %>%
  unique()

observed_cells_enhanced <- enhanced %>%
  filter(time <= end_frame) %>%
  pull(cell) %>%
  unique()

observed_cells <- intersect(
  observed_cells_raw,
  observed_cells_enhanced
)

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

cat("\n========== Filtering summary ==========\n")
cat("Observed cells in raw before end_frame:", length(observed_cells_raw), "\n")
cat("Observed cells in enhanced before end_frame:", length(observed_cells_enhanced), "\n")
cat("Observed cells in both datasets:", length(observed_cells), "\n")
cat("Death cells in lineage file:", length(death_cells), "\n")
cat("Cells with observed daughter cells:", length(dividing_cells), "\n")
cat("Cells without observed daughter cells:", length(no_daughter_cells), "\n")
cat("=======================================\n\n")


# =============================================================================
# Match raw and enhanced xyz by cellTime
# =============================================================================

xyz_pair <- inner_join(
  raw %>%
    filter(time <= end_frame) %>%
    select(cellTime, cell, time, x, y, z),
  enhanced %>%
    filter(time <= end_frame) %>%
    select(cellTime, cell, time, x, y, z),
  by = "cellTime",
  suffix = c("_raw", "_enhanced")
) %>%
  filter(cell_raw == cell_enhanced) %>%
  filter(time_raw == time_enhanced) %>%
  transmute(
    cell = cell_raw,
    time = time_raw,
    x_raw = x_raw,
    x_enhanced = x_enhanced,
    y_raw = y_raw,
    y_enhanced = y_enhanced,
    z_raw = z_raw,
    z_enhanced = z_enhanced
  ) %>%
  filter(!cell %in% death_cells) %>%
  filter(cell %in% dividing_cells) %>%
  filter(if_all(
    c(
      x_raw,
      x_enhanced,
      y_raw,
      y_enhanced,
      z_raw,
      z_enhanced
    ),
    is.finite
  ))

cat("Final matched xyz observations:", nrow(xyz_pair), "\n")
cat("Final matched dividing cells:", n_distinct(xyz_pair$cell), "\n\n")

if (nrow(xyz_pair) == 0) {
  stop("No xyz pairs remain after filtering. Please check lineage IDs and cell names.")
}



# =============================================================================
# Convert to long format
# =============================================================================

xyz_long <- xyz_pair %>%
  pivot_longer(
    cols = c(
      x_raw,
      x_enhanced,
      y_raw,
      y_enhanced,
      z_raw,
      z_enhanced
    ),
    names_to = c("axis", ".value"),
    names_pattern = "([xyz])_(raw|enhanced)"
  )

# Make axis order explicit
xyz_long <- xyz_long %>%
  mutate(
    axis = factor(axis, levels = c("x", "y", "z"))
  )

# =============================================================================
# Summary statistics
# =============================================================================

agreement_summary <- xyz_long %>%
  group_by(axis) %>%
  summarise(
    n = n(),
    n_cells = n_distinct(cell),
    pearson_r = cor(raw, enhanced, method = "pearson"),
    pearson_p = cor.test(raw, enhanced, method = "pearson")$p.value,
    spearman_rho = cor(raw, enhanced, method = "spearman"),
    slope = unname(coef(lm(enhanced ~ raw))[2]),
    intercept = unname(coef(lm(enhanced ~ raw))[1]),
    r_squared = summary(lm(enhanced ~ raw))$r.squared,
    median_difference = median(enhanced - raw),
    mean_difference = mean(enhanced - raw),
    rmse = sqrt(mean((enhanced - raw)^2)),
    mae = mean(abs(enhanced - raw)),
    exact_percent = 100 * mean(raw == enhanced),
    within_1_unit_percent = 100 * mean(abs(enhanced - raw) <= 1),
    within_2_units_percent = 100 * mean(abs(enhanced - raw) <= 2),
    within_3_units_percent = 100 * mean(abs(enhanced - raw) <= 3),
    .groups = "drop"
  )

print(agreement_summary)

# =============================================================================
# Labels for plot
# =============================================================================

axis_labels <- agreement_summary %>%
  mutate(
    label = sprintf(
      paste0(
        "n = %d\n",
        "Cells = %d\n",
        "Pearson r = %.10f\n",
        "P = %.10f\n",
        "Slope = %.6f\n",
        "Intercept = %.6f\n",
        "R² = %.10f\n",
        "RMSE = %.3f\n",
        "MAE = %.3f"
      ),
      n,
      n_cells,
      pearson_r,
      pearson_p,
      slope,
      intercept,
      r_squared,
      rmse,
      mae
    )
  )

# =============================================================================
# Plot
# =============================================================================

p_xyz <- ggplot(xyz_long, aes(x = raw, y = enhanced)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    size = 0.45,
    alpha = 0.25,
    colour = "#A8C6E4"
  ) +
  geom_label(
    data = axis_labels,
    aes(
      x = Inf,
      y = -Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.04,
    vjust = -0.1,
    size = 2.4,
    lineheight = 0.92,
    linewidth = 0,
    fill = "white",
    alpha = 0.85
  ) +
  facet_wrap(
    ~axis,
    scales = "free",
    ncol = 3
  ) +
  theme_classic(
    base_size = 8,
    base_family = "sans"
  ) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.45
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(
      linewidth = 0.35,
      colour = "black"
    ),
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
  labs(
    x = "Raw coordinate",
    y = "Enhanced coordinate",
    title = "Raw vs enhanced xyz coordinate agreement in dividing cells"
  )

print(p_xyz)

# =============================================================================
# Save figure
# =============================================================================

ggsave(
  file.path(
    figure_dir,
    "raw_vs_enhanced_xyz_coordinate_agreement_dividing_cells_only.pdf"
  ),
  p_xyz,
  width = 9,
  height = 3.5
)

ggsave(
  file.path(
    figure_dir,
    "raw_vs_enhanced_xyz_coordinate_agreement_dividing_cells_only.png"
  ),
  p_xyz,
  width = 9,
  height = 3.5,
  dpi = 600
)

message("XYZ comparison written to: ", out_dir)

