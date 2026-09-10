#!/usr/bin/env Rscript

# Direct-value comparison of diploid and tetraploid timing variability using REML LMM with Kenward-Roger mean-effect tests.

suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tidyr)
  library(purrr); library(stringr); library(ggplot2); library(forcats)
  library(nlme); library(lme4); library(pbkrtest)
})

input_dir <- path.expand("~/Downloads/NCRevision/DuZhuo")
out_dir <- file.path(input_dir, "diploid_tetraploid_CV_CE_template_REML_KR_LMM_v2")
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "code"), recursive = TRUE, showWarnings = FALSE)
workbook <- file.path(input_dir, "mmc2 (3).xlsx")
stage_file <- path.expand("~/Desktop/HybridAnalysis/cbcn/Length_Time_Absolute.csv")
fate_file <- path.expand("~/Desktop/HybridAnalysis/cbcn/linage/AllLineage.tsv")
time_per_tp <- 1.25
ab2_cells <- c("ABa", "ABp", "EMS", "P2")

read_tracks <- function(sheet, group) {
  d <- read_excel(workbook, sheet = sheet, range = cell_cols(1:7), .name_repair = "minimal")
  names(d) <- c("embryo", "cell", "tp", "x_px", "y_px", "z_um", "diameter_px")
  d %>% filter(!is.na(embryo), !is.na(cell), !is.na(tp)) %>%
    transmute(ploidy = group, embryo = as.character(embryo), cell = as.character(cell), tp = as.numeric(tp), embryo_id = paste(group, embryo, sep = "_"))
}

tracks <- bind_rows(read_tracks("1", "Diploid"), read_tracks("2", "Tetraploid"))
stage_map <- read_csv(stage_file, show_col_types = FALSE) %>% transmute(cell = as.character(cell), stage = as.character(stage))
fate_map <- read_tsv(fate_file, show_col_types = FALSE) %>% transmute(cell = as.character(CellName), lineage_id = as.character(ID), cell_fate = as.character(CellFate))

# Recreate the eight-CE division-time template used in she1af16.R.
ce_config <- tribble(
  ~run, ~path, ~limit, ~scale,
  "ce1", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD191108plc1p1.csv", 205, 1.43,
  "ce2", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200109plc1p1.csv", 205, 1.43,
  "ce3", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200113plc1p3.csv", 195, 1.44,
  "ce4", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200113plc1p2.csv", 205, 1.44,
  "ce5", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200322plc1p2.csv", 195, 1.44,
  "ce6", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200323plc1p1.csv", 185, 1.44,
  "ce7", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200326plc1p3.csv", 220, 1.44,
  "ce8", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200326plc1p4.csv", 195, 1.44
)
ce_timelines <- pmap(ce_config, function(run, path, limit, scale) {
  read_csv(path.expand(path), show_col_types = FALSE) %>%
    filter(time <= limit) %>% group_by(cell) %>%
    summarise(!!paste0("Tim_", run) := max(time) * scale, .groups = "drop")
})
ce_times <- reduce(ce_timelines, full_join, by = "cell")
ce_time_columns <- grep("^Tim_ce", names(ce_times), value = TRUE)
for (column in ce_time_columns) {
  offset <- max(ce_times[[column]][ce_times$cell %in% ab2_cells], na.rm = TRUE)
  ce_times[[column]] <- ce_times[[column]] - offset
}
ce_template <- ce_times %>%
  transmute(cell, ce_template_division_time = rowMeans(across(all_of(ce_time_columns)), na.rm = TRUE)) %>%
  filter(is.finite(ce_template_division_time))

shared <- tracks %>% distinct(ploidy, cell) %>% count(cell) %>% filter(n == 2) %>% pull(cell)
annotation <- stage_map %>% inner_join(fate_map, by = "cell") %>% filter(cell %in% shared) %>% distinct(cell, .keep_all = TRUE)
stopifnot(nrow(annotation) == 515L)

starts <- tracks %>% group_by(ploidy, embryo, embryo_id) %>% summarise(start_tp = min(tp), .groups = "drop")
ab2_offsets <- tracks %>% filter(cell %in% ab2_cells) %>%
  group_by(ploidy, embryo, embryo_id) %>%
  summarise(n_ab2_cells = n_distinct(cell), ab2_offset_min = max(tp) * time_per_tp, .groups = "drop")
if (nrow(ab2_offsets) != 10L || any(ab2_offsets$n_ab2_cells < 4L)) stop("All four AB2 anchor cells are required in every embryo.")
cell_tracks <- tracks %>% inner_join(annotation, by = "cell") %>%
  group_by(ploidy, embryo, embryo_id, cell, stage, cell_fate, lineage_id) %>%
  summarise(first_tp = min(tp), last_tp = max(tp), .groups = "drop") %>% left_join(starts, by = c("ploidy", "embryo", "embryo_id"))

daughter_map <- fate_map %>% mutate(parent_id = if_else(nchar(lineage_id) > 1, str_sub(lineage_id, 1, -2), NA_character_)) %>%
  filter(!is.na(parent_id)) %>% group_by(parent_id) %>% summarise(daughters = list(cell), .groups = "drop")
daughter_first <- cell_tracks %>% select(ploidy, embryo, embryo_id, cell, first_tp)

events <- cell_tracks %>% left_join(daughter_map, by = c("lineage_id" = "parent_id")) %>%
  mutate(daughters = map(daughters, ~if (is.null(.x)) character() else .x)) %>%
  mutate(daughter_tp = pmap(list(ploidy, embryo, embryo_id, daughters), function(g, e, id, kids) {
    daughter_first %>% filter(ploidy == g, embryo == e, embryo_id == id, cell %in% kids) %>% pull(first_tp)
  })) %>%
  mutate(
    n_daughters = map_int(daughter_tp, length), daughter_min = map_dbl(daughter_tp, ~if(length(.x)) min(.x) else NA_real_), daughter_max = map_dbl(daughter_tp, ~if(length(.x)) max(.x) else NA_real_),
    division_valid = n_daughters == 2 & daughter_max - daughter_min <= 1 & daughter_min >= last_tp & daughter_min <= last_tp + 2,
    division_time_min = if_else(division_valid, last_tp * time_per_tp, NA_real_),
    cycle_valid = division_valid & first_tp > start_tp,
    cycle_length_min = if_else(cycle_valid, (last_tp - first_tp + 1) * time_per_tp, NA_real_)
  ) %>% left_join(ab2_offsets, by = c("ploidy", "embryo", "embryo_id")) %>%
  mutate(division_time_aligned_min = division_time_min - ab2_offset_min)

# Each diploid and tetraploid embryo is mapped to the same aligned eight-CE template.
regression_coefficients <- events %>% filter(division_valid) %>%
  inner_join(ce_template, by = "cell") %>%
  group_by(ploidy, embryo, embryo_id) %>%
  group_modify(~ {
    fit <- lm(division_time_aligned_min ~ ce_template_division_time, data = .x)
    tibble(
      n_template_cells = nrow(.x), intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]),
      r_squared = summary(fit)$r.squared
    )
  }) %>% ungroup()
if (nrow(regression_coefficients) != 10L || any(!is.finite(regression_coefficients$slope)) || any(regression_coefficients$slope <= 0)) stop("Invalid CE-template regression for one or more embryos.")

values <- events %>% left_join(regression_coefficients, by = c("ploidy", "embryo", "embryo_id")) %>%
  transmute(ploidy, embryo_id, cell, stage, cell_fate,
    cell_cycle_length_raw = cycle_length_min,
    cell_cycle_length_normalized = cycle_length_min / slope,
    division_time_raw = if_else(division_valid, division_time_aligned_min, NA_real_),
    division_time_normalized = if_else(division_valid, (division_time_aligned_min - intercept) / slope, NA_real_)) %>%
  pivot_longer(-c(ploidy, embryo_id, cell, stage, cell_fate), names_to = "feature", values_to = "value") %>%
  separate(feature, into = c("endpoint", "scale"), sep = "_(?=[^_]+$)") %>% filter(is.finite(value), value > 0)

eligible <- values %>% group_by(endpoint, scale, cell) %>%
  summarise(diploid_n = n_distinct(embryo_id[ploidy == "Diploid"]), tetraploid_n = n_distinct(embryo_id[ploidy == "Tetraploid"]), .groups = "drop") %>%
  filter(diploid_n == 5L, tetraploid_n == 5L)
values <- values %>% inner_join(eligible %>% select(endpoint, scale, cell), by = c("endpoint", "scale", "cell"))

cell_cv <- values %>% group_by(endpoint, scale, ploidy, cell, stage, cell_fate) %>%
  summarise(n_embryos = n_distinct(embryo_id), mean_value = mean(value), sd_value = sd(value), cv = sd_value / mean_value, .groups = "drop")

# Mean LMM: Kenward-Roger test for ploidy mean effect.
# Variance LMM: REML comparison of equal residual variance versus ploidy-specific residual variance.
fit_lmm_statistics <- function(d) {
  d <- d %>% mutate(ploidy = factor(ploidy), embryo_id = factor(embryo_id), cell = factor(cell), log_value = log(value))
  tryCatch({
    mean_full <- lmer(log_value ~ ploidy + cell + (1 | embryo_id), data = d, REML = TRUE,
      control = lmerControl(check.conv.singular = "ignore"))
    mean_null <- update(mean_full, . ~ . - ploidy)
    mean_kr_p <- KRmodcomp(mean_full, mean_null)$stats["p.value"]

    variance_equal <- lme(log_value ~ ploidy + cell, random = ~ 1 | embryo_id, data = d, method = "REML",
      control = lmeControl(returnObject = TRUE))
    variance_ploidy <- update(variance_equal, weights = varIdent(form = ~ 1 | ploidy))
    variance_lrt <- anova(variance_equal, variance_ploidy)
    variance_p <- variance_lrt$`p-value`[2]
    sd_multipliers <- coef(variance_ploidy$modelStruct$varStruct, unconstrained = FALSE)
    tetraploid_sd_multiplier <- if ("Tetraploid" %in% names(sd_multipliers)) sd_multipliers[["Tetraploid"]] else sd_multipliers[[1]]
    variance_ratio <- as.numeric(tetraploid_sd_multiplier)^2
    tibble(mean_KR_p = as.numeric(mean_kr_p), variance_LMM_p = variance_p, tetraploid_to_diploid_residual_variance_ratio = variance_ratio, model_note = "REML LMM; cell fixed effect and embryo random intercept")
  }, error = function(e) {
    tibble(mean_KR_p = NA_real_, variance_LMM_p = NA_real_, tetraploid_to_diploid_residual_variance_ratio = NA_real_, model_note = paste("Model failed:", conditionMessage(e)))
  })
}

lmm_summary_for <- function(group_column = NULL, min_cells = 10L) {
  keys <- if (is.null(group_column)) values %>% distinct(endpoint_key = endpoint, scale_key = scale) %>% mutate(subgroup = "All cells") else values %>% distinct(endpoint_key = endpoint, scale_key = scale, subgroup = .data[[group_column]])
  out <- pmap_dfr(keys, function(endpoint_key, scale_key, subgroup) {
    d <- if (is.null(group_column)) values[values$endpoint == endpoint_key & values$scale == scale_key, ] else values[values$endpoint == endpoint_key & values$scale == scale_key & values[[group_column]] == subgroup, ]
    n_cells <- n_distinct(d$cell)
    base <- tibble(endpoint = endpoint_key, scale = scale_key, subgroup = as.character(subgroup), n_cells = n_cells)
    if (n_cells < min_cells) base %>% mutate(mean_KR_p = NA_real_, variance_LMM_p = NA_real_, tetraploid_to_diploid_residual_variance_ratio = NA_real_, model_note = paste0("Not tested: fewer than ", min_cells, " eligible cells")) else bind_cols(base, fit_lmm_statistics(d))
  })
  out %>% mutate(
    mean_KR_q_BH = p.adjust(mean_KR_p, method = "BH"),
    variance_LMM_q_BH = p.adjust(variance_LMM_p, method = "BH"),
    mean_significance = case_when(
      is.na(mean_KR_q_BH) ~ "not tested",
      mean_KR_q_BH < .0001 ~ "****",
      mean_KR_q_BH < .001 ~ "***",
      mean_KR_q_BH < .01 ~ "**",
      mean_KR_q_BH < .05 ~ "*",
      TRUE ~ "ns"
    ),
    variance_significance = case_when(
      is.na(variance_LMM_q_BH) ~ "not tested",
      variance_LMM_q_BH < .0001 ~ "****",
      variance_LMM_q_BH < .001 ~ "***",
      variance_LMM_q_BH < .01 ~ "**",
      variance_LMM_q_BH < .05 ~ "*",
      TRUE ~ "ns"
    )
  )
}

lmm_all <- lmm_summary_for()
lmm_stage <- lmm_summary_for("stage")
lmm_fate <- lmm_summary_for("cell_fate")

# The test uses embryo labels as the independent units; figures show only direct CV values.
exact_permutation_p <- function(d) {
  embryo_info <- d %>% distinct(embryo_id, ploidy) %>% arrange(embryo_id)
  wide <- d %>% select(cell, embryo_id, value) %>% pivot_wider(names_from = embryo_id, values_from = value) %>% arrange(cell)
  x <- as.matrix(wide[, embryo_info$embryo_id, drop = FALSE])
  statistic <- function(tet_indices) {
    dip_indices <- setdiff(seq_len(ncol(x)), tet_indices)
    median(apply(x[, tet_indices, drop = FALSE], 1, sd) / rowMeans(x[, tet_indices, drop = FALSE])) -
      median(apply(x[, dip_indices, drop = FALSE], 1, sd) / rowMeans(x[, dip_indices, drop = FALSE]))
  }
  observed <- statistic(which(embryo_info$ploidy == "Tetraploid"))
  null <- apply(combn(seq_len(ncol(x)), 5L), 2, statistic)
  mean(abs(null) >= abs(observed) - .Machine$double.eps^0.5)
}

summary_for <- function(group_column = NULL, min_cells = 10L) {
  keys <- if (is.null(group_column)) values %>% distinct(endpoint_key = endpoint, scale_key = scale) %>% mutate(subgroup = "All cells") else values %>% distinct(endpoint_key = endpoint, scale_key = scale, subgroup = .data[[group_column]])
  pmap_dfr(keys, function(endpoint_key, scale_key, subgroup) {
    if (is.null(group_column)) {
      d <- values[values$endpoint == endpoint_key & values$scale == scale_key, ]
      c <- cell_cv[cell_cv$endpoint == endpoint_key & cell_cv$scale == scale_key, ]
    } else {
      d <- values[values$endpoint == endpoint_key & values$scale == scale_key & values[[group_column]] == subgroup, ]
      c <- cell_cv[cell_cv$endpoint == endpoint_key & cell_cv$scale == scale_key & cell_cv[[group_column]] == subgroup, ]
    }
    med <- c %>% group_by(ploidy) %>% summarise(median_cv = median(cv), .groups = "drop") %>% pivot_wider(names_from = ploidy, values_from = median_cv)
    n_cells <- n_distinct(c$cell)
    tibble(endpoint = endpoint_key, scale = scale_key, subgroup = as.character(subgroup), n_cells = n_cells,
      diploid_median_cv = med$Diploid, tetraploid_median_cv = med$Tetraploid,
      exact_permutation_p = NA_real_,
      test_note = "Not calculated in this LMM-focused package; see the earlier permutation-analysis package if needed.")
  }) %>% mutate(exact_permutation_q_BH = p.adjust(exact_permutation_p, method = "BH"))
}

all_statistics <- summary_for()
stage_statistics <- summary_for("stage")
fate_statistics <- summary_for("cell_fate")

write_csv(annotation, file.path(out_dir, "tables/cell_stage_fate_annotations.csv"))
write_csv(ce_template, file.path(out_dir, "tables/eight_CE_division_time_template.csv"))
write_csv(regression_coefficients, file.path(out_dir, "tables/diploid_tetraploid_CE_template_regression_coefficients.csv"))
write_csv(eligible, file.path(out_dir, "tables/eligible_cells_by_endpoint.csv"))
write_csv(values, file.path(out_dir, "tables/LMM_input_timing_values.csv"))
write_csv(cell_cv %>% select(endpoint, scale, ploidy, cell, stage, cell_fate, cv), file.path(out_dir, "tables/TableS1_per_cell_CV_values.csv"))
write_csv(all_statistics, file.path(out_dir, "tables/Table1_all_cells_diploid_tetraploid_CV.csv"))
write_csv(stage_statistics, file.path(out_dir, "tables/Table2_stage_diploid_tetraploid_CV.csv"))
write_csv(fate_statistics, file.path(out_dir, "tables/Table3_cell_fate_diploid_tetraploid_CV.csv"))
write_csv(lmm_all, file.path(out_dir, "tables/Table4_all_cells_LMM_statistics.csv"))
write_csv(lmm_stage, file.path(out_dir, "tables/Table5_stage_LMM_statistics.csv"))
write_csv(lmm_fate, file.path(out_dir, "tables/Table6_cell_fate_LMM_statistics.csv"))

palette <- c(Diploid = "#CD534C", Tetraploid = "#5b66a1")
theme_set(theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 10), strip.text = element_text(face = "bold"), legend.position = "top"))
save_plot <- function(p, file, width, height) { ggsave(file.path(out_dir, "figures", paste0(file, ".png")), p, width = width, height = height, dpi = 600); ggsave(file.path(out_dir, "figures", paste0(file, ".pdf")), p, width = width, height = height) }
format_q <- function(q) ifelse(is.na(q), "not tested", paste0("q = ", sprintf("%.3f", q)))
significance_note <- "LMM labels: mean KR and variance model, BH-adjusted q; ns q >= 0.05; * q < 0.05; ** q < 0.01; *** q < 0.001."

make_all_plot <- function(metric) {
  labels <- lmm_all %>% filter(endpoint == metric) %>% mutate(label = paste0("n = 5 embryos/group\nMean KR: ", mean_significance, " ", format_q(mean_KR_q_BH), "\nVariance LMM: ", variance_significance, " ", format_q(variance_LMM_q_BH)))
  ggplot(cell_cv %>% filter(endpoint == metric), aes(ploidy, cv, fill = ploidy)) +
    stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  # 添加须线端帽（横线）
    geom_boxplot(width = .48, outlier.shape = NA, alpha = .65, colour = "#333333") +
    #geom_jitter(aes(colour = ploidy), width = .14, height = 0, size = .65, alpha = .36) +
    #geom_point(stat = "summary", fun = median, shape = 23, size = 2.5, fill = "white", colour = "#222222") +
    geom_text(data = labels, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.25, size = 2.8) +
    facet_wrap(~scale, nrow = 1, scales = "free_y") + scale_fill_manual(values = palette) + scale_colour_manual(values = palette) +
    labs(title = paste0("All cells: ", if(metric == "cell_cycle_length") "cell-cycle length" else "division time", " variability"), subtitle = paste("Each point is one cell's CV across five embryos.", significance_note), x = NULL, y = "Per-cell CV (SD / mean)") +
    theme(plot.subtitle = element_text(size = 7), legend.title = element_blank())
}

make_group_plot <- function(metric, group_column, title, x_label, order = NULL) {
  d <- cell_cv %>% filter(endpoint == metric)
  if (!is.null(order)) d <- d %>% mutate(group = factor(.data[[group_column]], levels = order)) else d <- d %>% mutate(group = .data[[group_column]])
  lmm_table <- if (group_column == "stage") lmm_stage else lmm_fate
  labels <- lmm_table %>% filter(endpoint == metric) %>% transmute(group = subgroup, scale, label = paste0("M:", mean_significance, " V:", variance_significance)) %>%
    mutate(group = if (!is.null(order)) factor(group, levels = order) else group)
  ggplot(d, aes(group, cv, fill = ploidy, colour = ploidy)) +
    geom_boxplot(position = position_dodge(width = .72), width = .62, outlier.shape = NA, alpha = .62, linewidth = .4) +
    geom_point(position = position_jitterdodge(jitter.width = .12, dodge.width = .72), size = .45, alpha = .28) +
    geom_text(data = labels, aes(x = group, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.25, size = 2.4, colour = "#222222") +
    facet_wrap(~scale, nrow = 1, scales = "free_y") + scale_fill_manual(values = palette) + scale_colour_manual(values = palette) +
    labs(title = title, subtitle = paste("Each label shows mean-model Kenward-Roger and variance-model LMM significance for the diploid-versus-tetraploid comparison.", significance_note), x = x_label, y = "Per-cell CV (SD / mean)") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.subtitle = element_text(size = 7), legend.title = element_blank())
}

stage_order <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32", "Stage_AB64", "Stage_AB128", "Stage_AB256")
for (metric in c("cell_cycle_length", "division_time")) {
  label <- if(metric == "cell_cycle_length") "Cell-cycle length" else "Division time"
  save_plot(make_all_plot(metric), paste0("Figure_", metric, "_all_cells"), 6.4, 4.3)
  save_plot(make_group_plot(metric, "stage", paste0("By developmental stage: ", label, " variability"), "Developmental stage", stage_order), paste0("Figure_", metric, "_by_stage"), 8.0, 4.5)
  save_plot(make_group_plot(metric, "cell_fate", paste0("By cell fate: ", label, " variability"), "Cell fate"), paste0("Figure_", metric, "_by_cell_fate"), 8.0, 4.5)
}

p1 <- make_all_plot("cell_cycle_length")
print(p1)

p1 <- make_all_plot("cell_cycle_length") +
  coord_cartesian(ylim = c(0, 0.06)) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    )
  )

print(p1)

ggsave(
  filename = file.path(out_dir, "figures", "Revised_cellLength.pdf"),
  plot = p1,
  width = 6.4,
  height = 4.3
)
writeLines(c(
  "# Direct diploid-versus-tetraploid CV values", "", "Figures show the actual per-cell CV values separately for diploid and tetraploid embryos; no difference-value axis is used.",
  "A cell's CV is SD / mean of that cell's endpoint across the five embryos of the corresponding ploidy group.",
  "The normalization exactly follows the eight-CE template approach. For each CE embryo, division times are converted using its acquisition scale and shifted by the maximum ABa/ABp/EMS/P2 time. The template is the per-cell mean of these eight aligned CE division times.",
  "For every diploid and tetraploid embryo, aligned division time is regressed on the CE template: time = intercept + slope * CE template. Normalized division time is (aligned time - intercept) / slope. Normalized cell-cycle length is raw cell-cycle length / slope.",
  "Primary significance analysis used two mixed-effects models on log-transformed timing values. Both models included ploidy and cell identity as fixed effects and embryo as a random intercept. The mean model was fitted by REML with lme4 and the ploidy effect was tested using a Kenward-Roger comparison to the corresponding model without ploidy. The variance model was fitted by REML with nlme and compared equal residual variance against ploidy-specific residual variance using varIdent; its likelihood-ratio P value tests whether residual timing variability differs by ploidy after accounting for cell identity and embryo.",
  "Figure labels report the Benjamini-Hochberg-adjusted q value and ns/*/**/*** code from the variance mixed-effects model. Tables 4-6 report both the mean-model Kenward-Roger and variance-model tests. The legacy permutation-test columns are retained only for compatibility and are not recalculated in this LMM-focused package.",
  "Stage and fate comparisons with fewer than 10 eligible cells are listed but not tested."
), file.path(out_dir, "methods.md"))
writeLines(c(
  "# Mixed-effects model specification", "",
  "## Purpose", "To assess diploid-versus-tetraploid timing differences without counting measurements from the same embryo as independent biological replicates.", "",
  "## Input", "The analysis uses only lineage cells with an observed, QC-passing division in all five diploid and all five tetraploid embryos. Timing values are converted to minutes and normalized by the common eight-C. elegans-embryo template exactly as specified in `she1af16.R`.", "",
  "## Mean model", "For each endpoint and each subgroup, we fitted `log(timing value) ~ ploidy + cell identity + (1 | embryo)`. Ploidy and cell identity are fixed effects; embryo is a random intercept. The ploidy mean effect was evaluated by a Kenward-Roger comparison of the full model with the corresponding no-ploidy model, both fitted by REML.", "",
  "## Variability model", "Using the same fixed-effect and embryo-random-intercept structure, we fitted two REML `nlme::lme` models: (i) equal residual variance in both ploidy groups and (ii) a `varIdent` residual-variance term for ploidy. Their likelihood-ratio test asks whether residual timing variability differs between tetraploid and diploid embryos after accounting for lineage cell identity and the shared embryo effect.", "",
  "## Multiplicity and figure labels", "P values are Benjamini-Hochberg adjusted within each table. The labels in the CV figures refer to the FDR-adjusted variance-model q value: ns, q >= 0.05; *, q < 0.05; **, q < 0.01; ***, q < 0.001. Mean-model and variance-model results are both reported in Tables 4-6.", "",
  "## Interpretation", "A non-significant variance-model result means these data do not provide evidence that residual timing variability differs by ploidy under this model. It does not establish exact equality; effect estimates and confidence intervals should be considered alongside q values."
), file.path(out_dir, "LMM_method.md"))
writeLines(c("# Contents", "figures/: six direct-value figures (all cells, stage, and fate for both endpoints), labeled by variance-LMM significance.", "tables/: direct CV summaries, LMM statistics, and every LMM input timing value.", "code/: reproducible R script."), file.path(out_dir, "README.md"))
message("Completed: ", out_dir)

