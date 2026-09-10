#!/usr/bin/env Rscript

# Selected timing, coordinate, and coverage evidence panels.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- path.expand("~/Downloads/NCRevision/enhanced_raw")
raw_file <- file.path(base_dir, "raw_CD191108plc1p1.csv")
enhanced_file <- file.path(base_dir, "enhanced_CD191108plc1p1.csv")
lineage_file <- file.path(base_dir, "AllLineage.tsv")
out_dir <- file.path(base_dir, "selected_all_nondeath_evidence")
figure_dir <- file.path(out_dir, "figures")
source_dir <- file.path(out_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

end_frame <- 205L
raw <- read.csv(raw_file, stringsAsFactors = FALSE, check.names = FALSE)
enhanced <- read.csv(enhanced_file, stringsAsFactors = FALSE, check.names = FALSE)
lineage <- read.delim(lineage_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("cellTime", "cell", "time", "x", "y", "z") %in% names(raw)) ||
    !all(c("cellTime", "cell", "time", "x", "y", "z") %in% names(enhanced))) {
  stop("Both CD files must contain cellTime, cell, time, x, y, and z.")
}
if (!all(c("CellName", "CellFate") %in% names(lineage))) {
  stop("AllLineage.tsv must contain CellName and CellFate.")
}

agreement <- function(x, y, metric) {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  fit <- lm(y ~ x)
  tibble(
    metric = metric,
    n = length(x),
    exact_percent = 100 * mean(x == y),
    median_difference = median(y - x),
    rmse = sqrt(mean((y - x)^2)),
    mae = mean(abs(y - x)),
    slope = unname(coef(fit)[2]),
    r_squared = summary(fit)$r.squared,
    spearman_rho = cor(x, y, method = "spearman"),
    within_3_frames_percent = 100 * mean(abs(y - x) <= 3)
  )
}

extract_cell_metrics <- function(df, suffix) {
  df %>%
    filter(time <= end_frame) %>%
    group_by(cell) %>%
    summarise(
      birth_time = min(time),
      observed_terminal_time = max(time),
      observed_tracking_span = max(time) - min(time) + 1,
      .groups = "drop"
    ) %>%
    rename_with(~paste0(.x, "_", suffix), -cell)
}

raw_cells <- extract_cell_metrics(raw, "raw")
enhanced_cells <- extract_cell_metrics(enhanced, "enhanced")
shared_cells <- inner_join(raw_cells, enhanced_cells, by = "cell")
death_cells <- lineage %>% filter(CellFate == "Death") %>% pull(CellName)
nondeath_cells <- shared_cells %>% filter(!cell %in% death_cells)

nondeath_timing_summary <- bind_rows(
  agreement(nondeath_cells$observed_terminal_time_raw, nondeath_cells$observed_terminal_time_enhanced,
            "Observed terminal time"),
  agreement(nondeath_cells$observed_tracking_span_raw, nondeath_cells$observed_tracking_span_enhanced,
            "Observed tracking span")
)
all_terminal_summary <- agreement(shared_cells$observed_terminal_time_raw,
                                  shared_cells$observed_terminal_time_enhanced,
                                  "Observed terminal time, all shared cells")

# Panel c: per-cell mean xyz coordinate agreement after death exclusion.
shared_xyz <- inner_join(
  raw %>% filter(time <= end_frame) %>% select(cellTime, cell, time, x, y, z),
  enhanced %>% filter(time <= end_frame) %>% select(cellTime, cell, time, x, y, z),
  by = "cellTime", suffix = c("_raw", "_enhanced")
) %>%
  filter(cell_raw == cell_enhanced, time_raw == time_enhanced) %>%
  transmute(cell = cell_raw, time = time_raw,
            x_raw, y_raw, z_raw, x_enhanced, y_enhanced, z_enhanced) %>%
  filter(if_all(c(x_raw, y_raw, z_raw, x_enhanced, y_enhanced, z_enhanced), is.finite)) %>%
  filter(!cell %in% death_cells)

cell_centroids <- shared_xyz %>%
  group_by(cell) %>%
  summarise(
    x_raw = mean(x_raw), x_enhanced = mean(x_enhanced),
    y_raw = mean(y_raw), y_enhanced = mean(y_enhanced),
    z_raw = mean(z_raw), z_enhanced = mean(z_enhanced),
    .groups = "drop"
  )
centroid_long <- cell_centroids %>%
  pivot_longer(cols = -cell, names_to = c("axis", ".value"),
               names_pattern = "([xyz])_(raw|enhanced)")
centroid_summary <- centroid_long %>%
  group_by(axis) %>%
  summarise(
    n_cells = n(),
    slope = unname(coef(lm(enhanced ~ raw))[2]),
    r_squared = summary(lm(enhanced ~ raw))$r.squared,
    rmse = sqrt(mean((enhanced - raw)^2)),
    .groups = "drop"
  )

# Panel d: frame-wise detection/tracking coverage.
frame_raw <- raw %>% filter(time <= end_frame) %>% count(time, name = "n_raw")
frame_enhanced <- enhanced %>% filter(time <= end_frame) %>% count(time, name = "n_enhanced")
frame_coverage <- full_join(frame_raw, frame_enhanced, by = "time") %>%
  complete(time = seq(1, end_frame)) %>%
  mutate(n_raw = replace_na(n_raw, 0L), n_enhanced = replace_na(n_enhanced, 0L),
         difference = n_enhanced - n_raw)

death_filter_audit <- tibble(
  measure = c("Shared cells through frame 205", "Cells annotated Death", "Non-death cells retained"),
  n_cells = c(nrow(shared_cells), sum(shared_cells$cell %in% death_cells), nrow(nondeath_cells))
)

write.csv(nondeath_cells, file.path(source_dir, "nondeath_observed_timing_cell_metrics.csv"), row.names = FALSE)
write.csv(nondeath_timing_summary, file.path(source_dir, "nondeath_observed_timing_summary.csv"), row.names = FALSE)
write.csv(shared_cells, file.path(source_dir, "all_shared_observed_timing_cell_metrics.csv"), row.names = FALSE)
write.csv(all_terminal_summary, file.path(source_dir, "all_shared_observed_terminal_time_summary.csv"), row.names = FALSE)
write.csv(cell_centroids, file.path(source_dir, "nondeath_per_cell_mean_coordinates.csv"), row.names = FALSE)
write.csv(centroid_summary, file.path(source_dir, "nondeath_per_cell_mean_coordinate_summary.csv"), row.names = FALSE)
write.csv(frame_coverage, file.path(source_dir, "frame_detection_tracking_coverage.csv"), row.names = FALSE)
write.csv(death_filter_audit, file.path(source_dir, "death_filter_audit.csv"), row.names = FALSE)

theme_set(
  theme_classic(base_size = 7, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.text = element_text(size = 6.3, face = "bold"),
      plot.title = element_text(size = 8, face = "bold"),
      panel.grid = element_blank()
    )
)

timing_plot_data <- nondeath_cells %>%
  select(cell, observed_terminal_time_raw, observed_terminal_time_enhanced,
         observed_tracking_span_raw, observed_tracking_span_enhanced) %>%
  pivot_longer(cols = -cell, names_to = c("metric", ".value"),
               names_pattern = "(observed_terminal_time|observed_tracking_span)_(raw|enhanced)") %>%
  mutate(metric = recode(metric, observed_terminal_time = "Observed terminal time",
                          observed_tracking_span = "Observed tracking span"))
timing_labels <- nondeath_timing_summary %>%
  mutate(label = sprintf("n = %d\nExact = %.1f%%\nWithin 3 frames = %.1f%%\nSpearman rho = %.3f\nSlope = %.3f; R^2 = %.3f\nRMSE = %.2f",
                         n, exact_percent, within_3_frames_percent, spearman_rho, slope, r_squared, rmse))
fig_a <- ggplot(timing_plot_data, aes(x = raw, y = enhanced)) +
  geom_point(size = 0.55, alpha = 0.3, colour = "grey20") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#4DBBD5", linewidth = 0.5) +
  geom_label(data = timing_labels, aes(x = Inf, y = -Inf, label = label), inherit.aes = FALSE,
             hjust = 1.04, vjust = -0.1, size = 1.65, lineheight = 0.94,
             linewidth = 0, fill = "white", alpha = 0.85) +
  facet_wrap(~metric, scales = "free", ncol = 2, axes = "all", axis.labels = "all") +
  labs(x = "Raw-image CD value (frames)", y = "Enhanced-image CD value (frames)",
       title = "a  Observed timing after excluding cells annotated as Death")

terminal_label <- sprintf("n = %d\nExact = %.1f%%\nSlope = %.3f\nR^2 = %.3f\nRMSE = %.2f",
                          all_terminal_summary$n, all_terminal_summary$exact_percent,
                          all_terminal_summary$slope, all_terminal_summary$r_squared,
                          all_terminal_summary$rmse)
fig_b <- ggplot(shared_cells, aes(x = observed_terminal_time_raw, y = observed_terminal_time_enhanced)) +
  geom_point(size = 0.55, alpha = 0.3, colour = "grey20") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#4DBBD5", linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = terminal_label, hjust = -0.05, vjust = 1.1,
           size = 2.0, lineheight = 0.94) +
  labs(x = "Raw-image observed terminal time (frames)",
       y = "Enhanced-image observed terminal time (frames)",
       title = "b  Observed terminal time in all shared cells")

centroid_labels <- centroid_summary %>%
  mutate(label = sprintf("n = %d\nR^2 = %.4f\nSlope = %.4f\nRMSE = %.2f",
                         n_cells, r_squared, slope, rmse))
fig_c <- ggplot(centroid_long, aes(x = raw, y = enhanced)) +
  geom_point(size = 0.48, alpha = 0.24, colour = "grey20") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#4DBBD5", linewidth = 0.45) +
  geom_label(data = centroid_labels, aes(x = Inf, y = -Inf, label = label), inherit.aes = FALSE,
             hjust = 1.04, vjust = -0.1, size = 1.4, lineheight = 0.92,
             linewidth = 0, fill = "white", alpha = 0.85) +
  facet_wrap(~axis, ncol = 3, scales = "free", axes = "all", axis.labels = "all") +
  labs(x = "Raw mean coordinate", y = "Enhanced mean coordinate",
       title = "c  Mean-coordinate agreement")

fig_d <- ggplot(frame_coverage, aes(x = time)) +
  geom_line(aes(y = n_raw, colour = "Raw"), linewidth = 0.5) +
  geom_line(aes(y = n_enhanced, colour = "Enhanced"), linewidth = 0.5, linetype = "dashed") +
  scale_colour_manual(values = c(Raw = "grey15", Enhanced = "#4DBBD5")) +
  labs(x = "Imaging frame", y = "Number of nuclei", colour = NULL,
       title = "d  Detection and tracking coverage") +
  theme(legend.position = c(0.2, 0.87), legend.text = element_text(size = 6))

selected_figure <- fig_a / (fig_b | fig_c) / fig_d + plot_layout(heights = c(1, 0.85, 0.65))
ggsave(file.path(figure_dir, "Figure_selected_all_nondeath_evidence.png"), selected_figure, width = 7.0, height = 10.2, dpi = 600)
ggsave(file.path(figure_dir, "Figure_selected_all_nondeath_evidence.pdf"), selected_figure, width = 7.0, height = 10.2)

message("Selected all-nondeath evidence package written to: ", out_dir)
