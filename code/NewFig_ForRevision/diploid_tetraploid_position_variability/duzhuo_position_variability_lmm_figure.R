#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(pbkrtest)
})

input_dir <- "/Users/yiming/Downloads/NCRevision/DuZhuo/draft/diploid_tetraploid_position_variability_group_templates"
out_dir <- "/Users/yiming/Downloads/NCRevision/DuZhuo/diploid_tetraploid_position_variability_REML_LMM_cell_mean"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "code"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

ploidy_levels <- c("Diploid", "Tetraploid")
pal <- c(Diploid = "#CD534C", Tetraploid = "#5b66a1")
eps <- 1e-6

theme_set(
  theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 7.5),
      legend.position = "none"
    )
)

format_q <- function(q) {
  ifelse(is.na(q), "not tested", ifelse(q < 0.0001, "q<0.0001", sprintf("q=%.3f", q)))
}

sig_from_q <- function(q) {
  case_when(
    is.na(q) ~ "not tested",
    q < 0.0001 ~ "****",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

save_plot <- function(p, stem, width = 4.8, height = 4.0) {
  ggsave(file.path(out_dir, "figures", paste0(stem, ".png")), p, width = width, height = height, dpi = 600)
  ggsave(file.path(out_dir, "figures", paste0(stem, ".pdf")), p, width = width, height = height)
  ggsave(file.path(out_dir, "figures", paste0(stem, ".svg")), p, width = width, height = height, device = svglite::svglite)
}

cell_mean_values <- read_csv(file.path(input_dir, "tables", "Table2_cell_mean_position_variability.csv"), show_col_types = FALSE) %>%
  mutate(
    ploidy = droplevels(factor(ploidy, levels = ploidy_levels)),
    cell = factor(cell),
    log_mean_RMSD = log(mean_RMSD + eps)
  ) %>%
  filter(is.finite(mean_RMSD), mean_RMSD > 0)

write_csv(cell_mean_values, file.path(out_dir, "tables", "Table1_cell_mean_position_variability_LMM_input_and_plot_source.csv"))

# Make model data visible to pbkrtest when it reconstructs model calls.
lmm_data_for_pbkrtest <- cell_mean_values
full_model <- lmer(
  log_mean_RMSD ~ ploidy + (1 | cell),
  data = lmm_data_for_pbkrtest,
  REML = TRUE
)
reduced_model <- lmer(
  log_mean_RMSD ~ 1 + (1 | cell),
  data = lmm_data_for_pbkrtest,
  REML = TRUE
)

kr <- KRmodcomp(full_model, reduced_model)
sat <- SATmodcomp(full_model, reduced_model)

kr_p <- as.numeric(kr$stats[["p.value"]])
sat_p <- as.numeric(sat$test$p.value)

stats_table <- tibble(
  comparison = "All cells",
  model = "log_mean_RMSD ~ ploidy + (1 | cell)",
  fit_method = "REML",
  response = "log(mean_position_RMSD + 1e-6)",
  fixed_effect = "ploidy",
  random_effects = "cell",
  diploid_n_cell_values = sum(cell_mean_values$ploidy == "Diploid"),
  tetraploid_n_cell_values = sum(cell_mean_values$ploidy == "Tetraploid"),
  n_cell_identities = n_distinct(cell_mean_values$cell),
  kr_p = kr_p,
  kr_q_BH = kr_p,
  kr_significance = sig_from_q(kr_p),
  satterthwaite_p = sat_p,
  satterthwaite_q_BH = sat_p,
  satterthwaite_significance = sig_from_q(sat_p)
)
write_csv(stats_table, file.path(out_dir, "tables", "Table3_position_variability_REML_LMM_statistics.csv"))

summary_table <- cell_mean_values %>%
  group_by(ploidy) %>%
  summarise(
    n_cells = n(),
    mean_position_RMSD = mean(.data$mean_RMSD),
    median_position_RMSD = median(.data$mean_RMSD),
    sd_position_RMSD = sd(.data$mean_RMSD),
    mean_log_position_RMSD = mean(.data$log_mean_RMSD),
    median_log_position_RMSD = median(.data$log_mean_RMSD),
    .groups = "drop"
  )
write_csv(summary_table, file.path(out_dir, "tables", "Table4_position_variability_summary.csv"))

label_df <- tibble(
  x = 1.5,
  y = quantile(cell_mean_values$mean_RMSD, 0.995, na.rm = TRUE),
  label = paste0(
    "REML LMM on log(position RMSD)\n",
    "KR: ", stats_table$kr_significance, " ", format_q(stats_table$kr_q_BH), "\n",
    "Satt.: ", stats_table$satterthwaite_significance, " ", format_q(stats_table$satterthwaite_q_BH), "\n",
    "random: cell"
  )
)

p <- ggplot(cell_mean_values, aes(ploidy, mean_RMSD, color = ploidy, fill = ploidy)) +
  #stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  # 添加须线端帽（横线
  #geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.35, linewidth = 0.35, colour = "#222222") +
  stat_boxplot(
    geom = "errorbar",
    width = 0.42,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    linewidth = 0.35,
    colour = "black",
    # alpha = 0.65
  ) +
  #geom_point(position = position_jitter(width = 0.16, height = 0), size = 0.55, alpha = 0.35) +
  #stat_summary(fun = median, geom = "point", shape = 23, size = 2.7, fill = "white", colour = "#222222") +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = -0.05,
    size = 2.55
  ) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "All-cell 3D position variability",
    subtitle = "Each point is one cell; significance uses REML LMM on log(mean RMSD).",
    x = NULL,
    y = "Mean position RMSD"
  ) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    ),
    legend.title = element_blank()
  )

print(p)
save_plot(p, "Figure1_all_cell_position_variability_REML_LMM", width = 3, height = 5)

methods <- c(
  "# All-cell 3D position variability REML LMM comparison",
  "",
  "This package reuses the position-variability values generated in the group-template position-variability workflow. The plotted values are the same figure-level unit as the original all-cell position-variability figure: one mean position RMSD value per cell and ploidy group, after averaging across group-template aligned time points.",
  "",
  "For statistical testing, the same cell-level mean RMSD values shown in the figure were used. The response was transformed as log(mean_position_RMSD + 1e-6).",
  "",
  "The full REML-fitted linear mixed-effects model was: log_mean_RMSD ~ ploidy + (1 | cell). Ploidy was treated as the fixed effect of interest. Cell identity was included as a random effect because the same cell identities are compared between diploid and tetraploid groups and some cells have consistently higher or lower positional variability.",
  "",
  "The ploidy effect was tested by comparing the full model with a reduced model lacking ploidy: log_mean_RMSD ~ 1 + (1 | cell). Kenward-Roger and Satterthwaite tests were both reported."
)
writeLines(methods, file.path(out_dir, "methods.md"))

readme <- c(
  "# Diploid/tetraploid position variability REML LMM comparison",
  "",
  "figures/: all-cell position-variability figure with REML LMM significance.",
  "tables/: cell-level LMM input/plot source values, statistics, and summary values.",
  "code/: reproducible R script.",
  "",
  "Main test: REML-fitted LMM on log(mean_position_RMSD + 1e-6), with ploidy as fixed effect and cell identity as a random effect."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Completed: ", out_dir)
