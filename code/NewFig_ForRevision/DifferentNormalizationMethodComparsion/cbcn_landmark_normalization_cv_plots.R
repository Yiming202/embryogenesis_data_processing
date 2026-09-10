#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(nlme)
})

base_dir <- "/Users/yiming/Downloads/NCRevision/Normalization method/landmark_normalization_methods"
out_fig <- file.path(base_dir, "figures")
out_tab <- file.path(base_dir, "tables")
out_code <- file.path(base_dir, "code")

dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_code, recursive = TRUE, showWarnings = FALSE)

method_files <- tibble(
  method = c("Method 1", "Method 2", "Method 3", "Method 4"),
  method_id = c("M1", "M2", "M3", "RAW"),
  description = c(
    "Global AB4-to-AB256 landmark normalization",
    "Stage-wise landmark normalization",
    "Monotonic landmark time-warp normalization",
    "Raw absolute length without normalization"
  ),
  file = c(
    file.path(out_tab, "Method1_global_AB4_to_AB256_landmark_normalized_Length_Time.csv"),
    file.path(out_tab, "Method2_stagewise_landmark_normalized_Length_Time.csv"),
    file.path(out_tab, "Method3_monotonic_landmark_timewarp_normalized_Length_Time.csv"),
    file.path(out_tab, "Method4_raw_absolute_Length_Time_no_normalization.csv")
  )
)

species_levels <- c("ce", "cb", "cn")
species_labels <- c(ce = "C. elegans", cb = "C. briggsae", cn = "C. nigoni")
species_colors <- c(ce = "#CD534C", cb = "#0073C2", cn = "#EFC000")
comparisons <- list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn"))

theme_set(
  theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5),
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
)

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

format_q <- function(q) {
  ifelse(is.na(q), "not tested", ifelse(q < 0.001, "q<0.001", sprintf("q=%.3f", q)))
}

calc_cv <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < 2 || mean(values) == 0) return(NA_real_)
  sd(values) / mean(values)
}

build_cv_data <- function(method, method_id, description, file) {
  dat <- read_csv(file, show_col_types = FALSE)
  map_dfr(species_levels, function(sp) {
    len_cols <- if (method_id == "RAW") {
      grep(paste0("^Len_", sp, "p\\d+$"), names(dat), value = TRUE)
    } else {
      grep(paste0("^", method_id, "_Len_", sp, "p\\d+$"), names(dat), value = TRUE)
    }
    dat %>%
      rowwise() %>%
      mutate(
        cv = calc_cv(c_across(all_of(len_cols))),
        n_embryos = sum(is.finite(c_across(all_of(len_cols))))
      ) %>%
      ungroup() %>%
      transmute(
        method = method,
        method_id = method_id,
        description = description,
        cell,
        stage,
        species = sp,
        n_embryos = n_embryos,
        cv = cv
      )
  })
}

build_length_long_data <- function(method, method_id, description, file) {
  dat <- read_csv(file, show_col_types = FALSE)
  len_cols <- if (method_id == "RAW") {
    grep("^Len_", names(dat), value = TRUE)
  } else {
    grep(paste0("^", method_id, "_Len_"), names(dat), value = TRUE)
  }
  species_pattern <- if (method_id == "RAW") {
    "^Len_([a-z]+)p[0-9]+$"
  } else {
    paste0("^", method_id, "_Len_([a-z]+)p[0-9]+$")
  }
  embryo_pattern <- if (method_id == "RAW") {
    "^Len_"
  } else {
    paste0("^", method_id, "_Len_")
  }
  dat %>%
    select(cell, stage, all_of(len_cols)) %>%
    pivot_longer(
      cols = all_of(len_cols),
      names_to = "sample",
      values_to = "normalized_length"
    ) %>%
    mutate(
      method = method,
      method_id = method_id,
      description = description,
      species = sub(species_pattern, "\\1", sample),
      embryo = sub(embryo_pattern, "", sample),
      log_length = log(normalized_length + 1e-6)
    ) %>%
    filter(species %in% species_levels, is.finite(log_length), !is.na(stage)) %>%
    mutate(
      species = factor(species, levels = species_levels),
      cell = factor(cell),
      embryo = factor(embryo)
    ) %>%
    select(method, method_id, description, cell, stage, species, embryo, normalized_length, log_length)
}

cv_data <- pmap_dfr(method_files, build_cv_data) %>%
  mutate(
    species = factor(species, levels = species_levels),
    method = factor(method, levels = method_files$method)
  ) %>%
  filter(is.finite(cv), !is.na(stage))

write_csv(cv_data, file.path(out_tab, "Table4_landmark_normalized_cell_cycle_length_CV_source.csv"))

length_long_data <- pmap_dfr(method_files, build_length_long_data)
write_csv(length_long_data, file.path(out_tab, "Table5_landmark_normalized_length_long_variance_LMM_input.csv"))

variance_lmm_pair_test <- function(d_pair) {
  d_pair <- d_pair %>%
    mutate(
      species = droplevels(species),
      cell = factor(cell)
    ) %>%
    filter(is.finite(log_length))

  if (n_distinct(d_pair$species) != 2 || any(table(d_pair$species) < 2)) {
    return(tibble(
      variance_lmm_p = NA_real_,
      residual_variance_ratio_group2_to_group1 = NA_real_,
      n_observations = nrow(d_pair),
      n_cells = n_distinct(d_pair$cell),
      n_embryos = n_distinct(d_pair$embryo),
      model = "not tested"
    ))
  }

  tryCatch({
    equal_var <- lme(
      log_length ~ species + cell,
      random = ~1 | embryo,
      data = d_pair,
      method = "REML",
      control = lmeControl(returnObject = TRUE)
    )
    species_var <- update(equal_var, weights = varIdent(form = ~1 | species))
    var_test <- anova(equal_var, species_var)
    var_p <- as.numeric(var_test$`p-value`[2])
    sd_mult <- coef(species_var$modelStruct$varStruct, unconstrained = FALSE)
    var_ratio <- if (length(sd_mult) == 0) {
      1
    } else {
      as.numeric(sd_mult[[1]])^2
    }
    tibble(
      variance_lmm_p = var_p,
      residual_variance_ratio_group2_to_group1 = var_ratio,
      n_observations = nrow(d_pair),
      n_cells = n_distinct(d_pair$cell),
      n_embryos = n_distinct(d_pair$embryo),
      model = "nlme::lme log(length) ~ species + cell + (1 | embryo), varIdent residual variance by species"
    )
  }, error = function(e) {
    tibble(
      variance_lmm_p = NA_real_,
      residual_variance_ratio_group2_to_group1 = NA_real_,
      n_observations = nrow(d_pair),
      n_cells = n_distinct(d_pair$cell),
      n_embryos = n_distinct(d_pair$embryo),
      model = paste("model failed:", conditionMessage(e))
    )
  })
}

stat_one_method <- function(d) {
  map_dfr(comparisons, function(pair) {
    d_pair <- d %>% filter(species %in% pair) %>% mutate(species = droplevels(species))
    variance_lmm_pair_test(d_pair) %>%
      mutate(group1 = pair[1], group2 = pair[2], .before = 1)
  })
}

stats <- length_long_data %>%
  group_by(method, method_id, description) %>%
  group_modify(~stat_one_method(.x)) %>%
  ungroup() %>%
  group_by(method) %>%
  mutate(
    variance_lmm_q_BH = p.adjust(variance_lmm_p, method = "BH"),
    variance_significance = sig_from_q(variance_lmm_q_BH),
    variance_q_label = format_q(variance_lmm_q_BH)
  ) %>%
  ungroup()

write_csv(stats, file.path(out_tab, "Table6_landmark_normalized_CV_variance_LMM_statistics.csv"))

summary_table <- cv_data %>%
  group_by(method, method_id, description, species) %>%
  summarise(
    n_cells = n(),
    median_cv = median(cv),
    mean_cv = mean(cv),
    sd_cv = sd(cv),
    .groups = "drop"
  )

write_csv(summary_table, file.path(out_tab, "Table7_landmark_normalized_CV_summary.csv"))

annotation_df <- stats %>%
  mutate(
    x = as.numeric(factor(group1, levels = species_levels)),
    xend = as.numeric(factor(group2, levels = species_levels))
  ) %>%
  left_join(
    cv_data %>%
      group_by(method) %>%
      summarise(y_base = quantile(cv, 0.98, na.rm = TRUE), .groups = "drop"),
    by = "method"
  ) %>%
  group_by(method) %>%
  arrange(group1, group2, .by_group = TRUE) %>%
  mutate(
    y = y_base + row_number() * 0.020,
    label = paste0("Variance LMM ", variance_significance, "\n", variance_q_label)
  ) %>%
  ungroup()

plot_df <- cv_data %>%
  group_by(method, species) %>%
  mutate(
    q1 = quantile(cv, 0.25, na.rm = TRUE),
    q3 = quantile(cv, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    plot_keep = cv >= q1 - 1.5 * iqr & cv <= q3 + 1.5 * iqr
  ) %>%
  ungroup() %>%
  filter(plot_keep) %>%
  select(-q1, -q3, -iqr, -plot_keep)

make_cv_plot <- function(df, ann, method_name, desc) {
  ggplot(df, aes(species, cv, fill = species)) +
    stat_boxplot(geom = "errorbar", width = 0.36, linewidth = 0.35) +
    geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35, colour = "black") +
    geom_segment(
      data = ann,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.3
    ) +
    geom_segment(
      data = ann,
      aes(x = x, xend = x, y = y - 0.004, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.3
    ) +
    geom_segment(
      data = ann,
      aes(x = xend, xend = xend, y = y - 0.004, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.3
    ) +
    geom_text(
      data = ann,
      aes(x = (x + xend) / 2, y = y + 0.002, label = label),
      inherit.aes = FALSE,
      size = 2.15,
      lineheight = 0.85
    ) +
    scale_fill_manual(values = species_colors, labels = species_labels) +
    scale_x_discrete(labels = species_labels, expand = expansion(mult = c(0.25, 0.25))) +
    coord_cartesian(ylim = c(0, max(0.17, max(ann$y, na.rm = TRUE) + 0.015))) +
    labs(
      title = paste0(method_name, ": cell-cycle length CV"),
      subtitle = paste0(desc, "; labels show variance LMM q"),
      x = NULL,
      y = "CV of cell-cycle length"
    ) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.line = element_blank())
}

plots <- map(method_files$method, function(m) {
  desc <- method_files$description[method_files$method == m]
  p <- make_cv_plot(
    plot_df %>% filter(method == m),
    annotation_df %>% filter(method == m),
    m,
    desc
  )
  stem <- case_when(
    m == "Method 1" ~ "Figure1_method1_global_landmark_length_CV",
    m == "Method 2" ~ "Figure2_method2_stagewise_landmark_length_CV",
    m == "Method 3" ~ "Figure3_method3_monotonic_landmark_length_CV",
    TRUE ~ "Figure4_method4_raw_absolute_length_CV"
  )
  ggsave(file.path(out_fig, paste0(stem, ".pdf")), p, width = 4.6, height = 4.2)
  p
})

print(plots[[1]])
print(combined_plot)

names(plots) <- method_files$method

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined_plot <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) + patchwork::plot_layout(guides = "collect")
  ggsave(file.path(out_fig, "Figure5_method1_to_method4_length_CV_combined.pdf"), combined_plot, width = 6, height = 8.4)
}

methods_text <- c(
  "# CV plots for landmark-normalized cell-cycle length",
  "",
  "For each landmark-normalization method, per-cell CV was calculated separately for CE, CB, and CN as sd(normalized cell-cycle length across embryos) / mean(normalized cell-cycle length across embryos).",
  "",
  "Method 1 uses the global AB4-to-90th-percentile-AB256 landmark-normalized cell-cycle length table. Method 2 uses the stage-wise 90th-percentile landmark-normalized cell-cycle length table. Method 3 uses the monotonic 90th-percentile landmark time-warp normalized table. Method 4 uses the raw absolute cell-cycle-length table without normalization.",
  "",
  "Pairwise species comparisons used variance LMMs on the underlying normalized cell-cycle-length values, not Wilcoxon tests on CV values. For each species pair and normalization method, an equal-residual-variance REML model was fitted as log(length + 1e-6) ~ species + cell + (1 | embryo), with species and cell identity as fixed effects and embryo identity as a random effect. This was compared with a model allowing species-specific residual variance using nlme::varIdent(form = ~1 | species). Figure labels show Benjamini-Hochberg adjusted q values from this residual-variance comparison. For visual clarity only, boxplots omit points beyond 1.5 IQR within each species/method group; all finite CV values are retained in the CV source table."
)
writeLines(methods_text, file.path(base_dir, "methods_CV_plots.md"))

message("Completed CV plots in: ", base_dir)

