library(tidyverse)

#----------------------------
# 1. 基础设置
#----------------------------
base_dir <- path.expand("~/Downloads/NCRevision/she1_af16_normalized_concordance/R_exploratory_z_slope_0.95_1.00_six/")
source_dir <- file.path(base_dir, "source_data")
figure_dir <- file.path(base_dir, "figures")
stage_file <- file.path(base_dir, "input_metadata", "Length_Time_Absolute.csv")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
setwd(source_dir)

# helper: 读取 CSV
read_cd <- function(path) {
  read.csv(path, header = TRUE, stringsAsFactors = FALSE)
}

# helper: 生成列名
make_col <- function(prefix, tag) paste0(prefix, "_", tag)

# 配置表：每行对应一个数据集
cfg <- tribble(
  ~species, ~run,  ~path,                                                                    ~len_limit, ~len_scale, ~tim_limit, ~tim_scale,
  "ce",     "p1", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD191108plc1p1.csv",             205,        1.43,        205,        1.43,
  "ce",     "p2", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200109plc1p1.csv",             205,        1.43,        205,        1.43,
  "ce",     "p3", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200113plc1p3.csv",             195,        1.44,        195,        1.44,
  "ce",     "p4", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200113plc1p2.csv",             205,        1.44,        205,        1.44,
  "ce",     "p5", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200322plc1p2.csv",             195,        1.44,        195,        1.44,
  "ce",     "p6", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200323plc1p1.csv",             185,        1.44,        185,        1.44,
  "ce",     "p7", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200326plc1p3.csv",             220,        1.44,        220,        1.44,
  "ce",     "p8", "~/Desktop/HybridAnalysis/CDFile/C. elegans/CD200326plc1p4.csv",             195,        1.44,        195,        1.44,
  "she1",     "p1", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD240731cbhis72p1.csv",   165,        1.57,        165,        1.57,
  "she1",     "p2", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD240731cbhis72p2.csv",   175,        1.57,        175,        1.57,
  "she1",     "p3", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD240731cbhis72p3.csv",   170,        1.57,        170,        1.57,
  "she1",     "p4", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD241202cbhis72p1.csv",   160,        1.58,        160,        1.58,
  "she1",     "p5", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD241202cbhis72p2.csv",   165,        1.58,        165,        1.58,
  "she1",     "p6", "~/Desktop/HybridAnalysis/CDFile/C. briggsae(she1)/CD241202cbhis72p4.csv",   180,        1.58,        180,        1.58,
  "af16",     "p1", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260412cbhis72p1.csv",                 110,        1.640,        110,        1.640,
  "af16",     "p2", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260412cbhis72p3.csv",                 115,        1.640,        115,        1.640,
  "af16",     "p3", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260413cbhis72p5.csv",                 105,        1.800,        105,        1.800,
  "af16",     "p4", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260415cbhis72p1.csv",                 160,        1.240,        160,        1.240,
  "af16",     "p5", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260417cbhis72p2.csv",                 140,        1.300,        140,        1.300,
  "af16",     "p6", "~/Downloads/LineagingData/CBData_ZhongShanU/CDFiles/All/CD260505cbhis72p2.csv",                 125,        1.547,        125,        1.547
  
  
 )

#----------------------------
# 2. 处理长度/时间
#----------------------------
process_len <- function(df, limit, scale, label) {
  df %>%
    as_tibble() %>%
    filter(time <= limit) %>%
    count(cell, name = label) %>%
    mutate(!!label := !!sym(label) * scale)
}

process_tim <- function(df, limit, scale, label) {
  df %>%
    as_tibble() %>%
    filter(time <= limit) %>%
    group_by(cell) %>%
    summarise(!!label := max(time) * scale, .groups = "drop")
}

# 读取并计算
datasets <- cfg %>%
  mutate(
    raw     = map(path, read_cd),
    len_col = make_col("Len", paste0(species, run)),
    tim_col = make_col("Tim", paste0(species, run))
  ) %>%
  mutate(
    len_tbl = pmap(list(raw, len_limit, len_scale, len_col),
                   ~process_len(..1, ..2, ..3, ..4)),
    tim_tbl = pmap(list(raw, tim_limit, tim_scale, tim_col),
                   ~process_tim(..1, ..2, ..3, ..4))
  )

len_df <- reduce(datasets$len_tbl, full_join, by = "cell")
tim_df <- reduce(datasets$tim_tbl, full_join, by = "cell")

merge_all <- len_df %>%
  full_join(tim_df, by = "cell") 
#filter(!cell %in% c("ABa", "ABp", "EMS", "P2"))

# The original script later filters and exports `stage`; join its source table
# explicitly so the final package runs in a clean R session.
stage_map <- read.csv(stage_file, header = TRUE, stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  select(cell, stage) %>%
  distinct(cell, .keep_all = TRUE)
merge_all <- merge_all %>% left_join(stage_map, by = "cell")


# -------- Tim 列按 AB2 最大值做平移 --------
tim_cols <- grep("^Tim_", names(merge_all), value = TRUE)
ab2_cells <- c("ABa", "ABp", "EMS", "P2")

ab2_max_map <- sapply(
  tim_cols,
  function(col) {
    vals <- merge_all[[col]][merge_all$cell %in% ab2_cells]
    if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
  },
  simplify = TRUE, USE.NAMES = TRUE
)

for (col in tim_cols) {
  offset <- ab2_max_map[[col]]
  if (!is.na(offset)) {
    merge_all[[col]] <- merge_all[[col]] - offset
  }
}

#----------------------------
# 4. Read lineage and identify dividing cells
#----------------------------

lineage_file <- file.path(base_dir, "input_metadata", "AllLineage.tsv")

lineage <- read_tsv(
  lineage_file,
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

early_founder_cells <- c(
  "AB",
  "P1",
  "ABa",
  "ABp",
  "EMS",
  "P2"
)


#----------------------------
# 5. Keep only annotated, non-death, dividing cells
#----------------------------

merge_all_filtered <- merge_all %>%
  filter(!is.na(stage)) %>%
  filter(!cell %in% early_founder_cells) %>%
  filter(!cell %in% death_cells)

observed_cells <- merge_all_filtered %>%
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
cat("Early founder cells:", length(early_founder_cells), "\n")
cat("Cells with observed daughter cells:", length(dividing_cells), "\n")
cat("Cells without observed daughter cells:", length(no_daughter_cells), "\n")
cat("=====================================================\n\n")

merge_all_filtered <- merge_all_filtered %>%
  filter(cell %in% dividing_cells)

cat("Cells after keeping only dividing cells:", nrow(merge_all_filtered), "\n\n")

# write.csv(
#   merge_all_filtered,
#   "p1_Length_Time_Absolute_filtered.csv",
#   row.names = FALSE
# )


# =============================================================================
# Add-on: compare mean corrected division time of she1 vs af16
# =============================================================================

library(tidyverse)
library(ggplot2)

#----------------------------
# 4. Build CE template mean time
#----------------------------

ce_tim_cols <- paste0("Tim_cep", 1:8)
she1_tim_cols <- paste0("Tim_she1p", 1:6)
af16_tim_cols <- paste0("Tim_af16p", 1:6)

mergeTem <- merge_all_filtered %>%
  mutate(
    temTim_mean = rowMeans(
      select(., all_of(ce_tim_cols)),
      na.rm = TRUE
    )
  ) %>%
  select(cell, temTim_mean)

#----------------------------
# 5. Regression normalization helper
#----------------------------

process_tim_data <- function(data, template_tbl, column_name) {
  
  nor_data <- data %>%
    select(cell, all_of(column_name)) %>%
    left_join(template_tbl, by = "cell") %>%
    filter(
      is.finite(.data[[column_name]]),
      is.finite(temTim_mean)
    )
  
  model <- lm(nor_data[[column_name]] ~ temTim_mean, data = nor_data)
  
  intercept <- unname(coef(model)[1])
  slope <- unname(coef(model)[2])
  
  list(
    nor_data = nor_data,
    intercept = intercept,
    slope = slope,
    r_squared = summary(model)$r.squared,
    n_fit = nrow(nor_data)
  )
}

#----------------------------
# 6. Fit all Tim columns
#----------------------------

tim_columns <- c(
  ce_tim_cols,
  she1_tim_cols,
  af16_tim_cols
)

model_results <- purrr::set_names(tim_columns) %>%
  purrr::map(~ process_tim_data(merge_all_filtered, mergeTem, .x))


coeff_summary <- imap_dfr(
  model_results,
  ~ tibble(
    embryo = .y,
    n_fit = .x$n_fit,
    intercept = .x$intercept,
    slope = .x$slope,
    r_squared = .x$r_squared
  )
)

print(coeff_summary)

write.csv(
  coeff_summary,
  "she1_af16_time_normalization_coefficients.csv",
  row.names = FALSE
)

#----------------------------
# 7. Normalize all division-time columns
#----------------------------

norm_tim_all <- merge_all_filtered

for (col in tim_columns) {
  
  corrected_name <- sub("^Tim_", "corTim_", col)
  
  norm_tim_all[[corrected_name]] <-
    (merge_all_filtered[[col]] - model_results[[col]]$intercept) /
    model_results[[col]]$slope
}

#----------------------------
# 8. Mean corrected division time for she1 and af16
#----------------------------

cor_she1_tim_cols <- paste0("corTim_she1p", 1:6)
cor_af16_tim_cols <- paste0("corTim_af16p", 1:6)

mean_tim_compare <- norm_tim_all %>%
  mutate(
    she1_mean_corTim = rowMeans(
      select(., all_of(cor_she1_tim_cols)),
      na.rm = TRUE
    ),
    af16_mean_corTim = rowMeans(
      select(., all_of(cor_af16_tim_cols)),
      na.rm = TRUE
    ),
    n_she1 = rowSums(
      !is.na(select(., all_of(cor_she1_tim_cols)))
    ),
    n_af16 = rowSums(
      !is.na(select(., all_of(cor_af16_tim_cols)))
    )
  ) %>%
  filter(
    n_she1 > 0,
    n_af16 > 0,
    is.finite(she1_mean_corTim),
    is.finite(af16_mean_corTim)
  ) %>%
  select(
    cell,
    stage,
    she1_mean_corTim,
    af16_mean_corTim,
    n_she1,
    n_af16
  )

write.csv(
  mean_tim_compare,
  "she1_vs_af16_mean_corrected_division_time_source.csv",
  row.names = FALSE
)

#----------------------------
# 9. Summary statistics
#----------------------------

agreement_summary <- mean_tim_compare %>%
  summarise(
    measurement = "Corrected division time",
    n = n(),
    n_cells = n_distinct(cell),
    pearson_r = cor(
      she1_mean_corTim,
      af16_mean_corTim,
      method = "pearson"
    ),
    pearson_p = cor.test(
      she1_mean_corTim,
      af16_mean_corTim,
      method = "pearson"
    )$p.value,
    spearman_rho = cor(
      she1_mean_corTim,
      af16_mean_corTim,
      method = "spearman"
    ),
    slope = unname(
      coef(lm(af16_mean_corTim ~ she1_mean_corTim))[2]
    ),
    intercept = unname(
      coef(lm(af16_mean_corTim ~ she1_mean_corTim))[1]
    ),
    r_squared = summary(
      lm(af16_mean_corTim ~ she1_mean_corTim)
    )$r.squared,
    median_difference = median(
      af16_mean_corTim - she1_mean_corTim
    ),
    mean_difference = mean(
      af16_mean_corTim - she1_mean_corTim
    ),
    rmse = sqrt(
      mean((af16_mean_corTim - she1_mean_corTim)^2)
    ),
    mae = mean(
      abs(af16_mean_corTim - she1_mean_corTim)
    ),
    within_1_min_percent = 100 * mean(
      abs(af16_mean_corTim - she1_mean_corTim) <= 1
    ),
    within_2_min_percent = 100 * mean(
      abs(af16_mean_corTim - she1_mean_corTim) <= 2
    ),
    within_3_min_percent = 100 * mean(
      abs(af16_mean_corTim - she1_mean_corTim) <= 3
    )
  )

print(agreement_summary)

write.csv(
  agreement_summary,
  "she1_vs_af16_mean_corrected_division_time_summary.csv",
  row.names = FALSE
)

#----------------------------
# 10. Plot label
#----------------------------

plot_label <- agreement_summary %>%
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

#----------------------------
# 11. Plot she1 vs af16 mean corrected division time
#----------------------------

p_she1_af16_tim <- ggplot(
  mean_tim_compare,
  aes(
    x = she1_mean_corTim,
    y = af16_mean_corTim
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "black",
    linewidth = 1
  ) +
  # geom_smooth(
  #   method = "lm",
  #   se = FALSE,
  #   colour = "#C95745",
  #   linewidth = 0.6
  # ) +
  geom_point(
    size = 2.5,
    alpha = 0.55,
    colour = "#A8C6E4"
  ) +
  geom_label(
    data = plot_label,
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
  coord_equal() +
  xlim(0, 200) +
  ylim(0, 200) +
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
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5
    ),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Mean corrected she1 division time",
    y = "Mean corrected af16 division time",
    title = "Corrected division time comparison: she1 vs af16"
  )

print(p_she1_af16_tim)

ggsave(
  file.path(figure_dir, "Figure1_SHE1_AF16_corrected_division_time_agreement.pdf"),
  p_she1_af16_tim,
  width = 5,
  height = 5
)

ggsave(
  file.path(figure_dir, "Figure1_SHE1_AF16_corrected_division_time_agreement.png"),
  p_she1_af16_tim,
  width = 5,
  height = 5,
  dpi = 600
)

