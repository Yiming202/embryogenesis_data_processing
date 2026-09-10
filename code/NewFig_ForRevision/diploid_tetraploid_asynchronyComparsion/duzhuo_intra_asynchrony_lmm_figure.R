#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(pbkrtest)
})

input_dir <- "/Users/yiming/Downloads/NCRevision/DuZhuo/draft/diploid_tetraploid_normalized_asynchrony_CE_template"
out_dir <- "/Users/yiming/Downloads/NCRevision/DuZhuo/diploid_tetraploid_intra_asynchrony_LMM_REML"

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

save_plot <- function(p, stem, width = 4.2, height = 3.8) {
  ggsave(file.path(out_dir, "figures", paste0(stem, ".png")), p, width = width, height = height, dpi = 600)
  ggsave(file.path(out_dir, "figures", paste0(stem, ".pdf")), p, width = width, height = height)
  ggsave(file.path(out_dir, "figures", paste0(stem, ".svg")), p, width = width, height = height, device = svglite::svglite)
}

dip_intra <- read_csv(file.path(input_dir, "tables", "TableS3_diploid_intra_pairwise_RMSD_values.csv"), show_col_types = FALSE) %>%
  transmute(
    metric = "Intra",
    ploidy = "Diploid",
    embryo_pair = paste0("Diploid_pair_", rmsd_index),
    cell,
    stage,
    rmsd = diploid_intra_rmsd
  )

tetra_intra <- read_csv(file.path(input_dir, "tables", "TableS4_tetraploid_intra_pairwise_RMSD_values.csv"), show_col_types = FALSE) %>%
  transmute(
    metric = "Intra",
    ploidy = "Tetraploid",
    embryo_pair = paste0("Tetraploid_pair_", rmsd_index),
    cell,
    stage,
    rmsd = tetraploid_intra_rmsd
  )

intra_values <- bind_rows(dip_intra, tetra_intra) %>%
  mutate(
    ploidy = droplevels(factor(ploidy, levels = ploidy_levels)),
    embryo_pair = factor(embryo_pair),
    cell = factor(cell),
    log_rmsd = log(rmsd + eps)
  ) %>%
  filter(!is.na(stage), is.finite(rmsd), rmsd > 0)

write_csv(intra_values, file.path(out_dir, "tables", "Table1_intra_cell_level_RMSD_values.csv"))

# Make model data available to pbkrtest when it reconstructs model calls.
lmm_data_for_pbkrtest <- intra_values
full_model <- lmer(log_rmsd ~ ploidy + (1 | embryo_pair) + (1 | cell), data = lmm_data_for_pbkrtest, REML = TRUE)
reduced_model <- lmer(log_rmsd ~ 1 + (1 | embryo_pair) + (1 | cell), data = lmm_data_for_pbkrtest, REML = TRUE)

kr <- KRmodcomp(full_model, reduced_model)
sat <- SATmodcomp(full_model, reduced_model)

kr_p <- as.numeric(kr$stats[["p.value"]])
sat_p <- as.numeric(sat$test$p.value)

stats_table <- tibble(
  metric = "Intra",
  model = "log_RMSD ~ ploidy + (1 | embryo_pair) + (1 | cell)",
  fit_method = "REML",
  response = "log(cell_RMSD + 1e-6)",
  fixed_effect = "ploidy",
  random_effects = "embryo_pair; cell",
  diploid_n_cell_values = sum(intra_values$ploidy == "Diploid"),
  tetraploid_n_cell_values = sum(intra_values$ploidy == "Tetraploid"),
  diploid_n_embryo_pairs = n_distinct(intra_values$embryo_pair[intra_values$ploidy == "Diploid"]),
  tetraploid_n_embryo_pairs = n_distinct(intra_values$embryo_pair[intra_values$ploidy == "Tetraploid"]),
  kr_p = kr_p,
  kr_q_BH = kr_p,
  kr_significance = sig_from_q(kr_p),
  satterthwaite_p = sat_p,
  satterthwaite_q_BH = sat_p,
  satterthwaite_significance = sig_from_q(sat_p)
)

write_csv(stats_table, file.path(out_dir, "tables", "Table2_intra_LMM_REML_statistics.csv"))

summary_table <- intra_values %>%
  group_by(ploidy) %>%
  summarise(
    n_cell_values = n(),
    n_cells = n_distinct(cell),
    n_embryo_pairs = n_distinct(embryo_pair),
    mean_rmsd = mean(rmsd),
    median_rmsd = median(rmsd),
    sd_rmsd = sd(rmsd),
    mean_log_rmsd = mean(log_rmsd),
    median_log_rmsd = median(log_rmsd),
    .groups = "drop"
  )

write_csv(summary_table, file.path(out_dir, "tables", "Table3_intra_RMSD_summary.csv"))

label_df <- tibble(
  x = 1.5,
  y = quantile(intra_values$rmsd, 0.98, na.rm = TRUE),
  label = paste0(
    "REML LMM on log(cell RMSD)\n",
    "KR: ", stats_table$kr_significance, " ", format_q(stats_table$kr_q_BH), "\n",
    "Satt.: ", stats_table$satterthwaite_significance, " ", format_q(stats_table$satterthwaite_q_BH), "\n",
    "random: embryo-pair + cell"
  )
)

p <- ggplot(intra_values, aes(ploidy, rmsd, color = ploidy)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  # 添加须线端帽（横线）
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.35) +
  #geom_point(position = position_jitter(width = 0.18, height = 0), size = 0.45, alpha = 0.28) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 8) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = -0.05,
    size = 2.55
  ) +
  scale_color_manual(values = pal) +
  labs(
    title = "Intra-asynchrony comparison",
    subtitle = "Per-cell intra-RMSD values; significance uses REML LMM on log(RMSD + 1e-6).",
    x = NULL,
    y = "Intra-asynchrony RMSD"
  ) +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    )
  )

p <- ggplot(intra_values, aes(ploidy, rmsd, fill = ploidy)) +
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
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    size = 8,
    colour = "black"
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = -0.05,
    size = 2.55,
    colour = "black"
  ) +
  scale_fill_manual(values = pal) +
  labs(
    title = "Intra-asynchrony comparison",
    subtitle = "Per-cell intra-RMSD values; significance uses REML LMM on log(RMSD + 1e-6).",
    x = NULL,
    y = "Intra-asynchrony RMSD"
  ) +
  coord_cartesian(ylim = c(0, 0.2)) +
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

print(p)

save_plot(p, "Figure1_intra_asynchrony_LMM_REML", width = 3, height = 5)

methods <- c(
  "# Intra-asynchrony LMM comparison",
  "",
  "Per-cell intra-asynchrony RMSD values were read from the normalized CE-template asynchrony source tables for diploid and tetraploid embryo-pair comparisons.",
  "",
  "The figure plots the original per-cell RMSD values so that the magnitude of asynchrony remains directly visible. For statistical testing, RMSD values were transformed as log(RMSD + 1e-6).",
  "",
  "The primary model was fitted by restricted maximum likelihood (REML): log_RMSD ~ ploidy + (1 | embryo_pair) + (1 | cell). Ploidy was treated as the fixed effect of interest. Embryo-pair identity was included as a random effect because many cell-level RMSD values come from the same embryo-pair comparison. Cell identity was included as a random effect because some cells are naturally more or less variable than others.",
  "",
  "The ploidy effect was tested by comparing the full model with a reduced model lacking ploidy: log_RMSD ~ 1 + (1 | embryo_pair) + (1 | cell). Kenward-Roger and Satterthwaite tests were both reported."
)
writeLines(methods, file.path(out_dir, "methods.md"))

readme <- c(
  "# Diploid/tetraploid intra-asynchrony LMM comparison",
  "",
  "figures/: manuscript-style intra-asynchrony RMSD comparison figure.",
  "tables/: cell-level source data, LMM statistics, and summary values.",
  "code/: reproducible R script.",
  "",
  "Main test: REML-fitted LMM on log(cell_RMSD + 1e-6), with ploidy as fixed effect and embryo-pair plus cell identity as random effects."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Completed: ", out_dir)
