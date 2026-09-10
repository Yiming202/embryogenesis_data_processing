#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(tibble)
})

out_dir <- "/Users/yiming/Desktop/HybridAnalysis/cbcn_submit/code/landmark_normalization_methods"
absolute_file <- "/Users/yiming/Desktop/HybridAnalysis/NC_Revision/Length_Time_Absolute.csv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "code"), recursive = TRUE, showWarnings = FALSE)

datasets_config <- tribble(
  ~group, ~id, ~file_path, ~time_limit, ~multiplier,
  "ce", 1, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD191108plc1p1.csv", 205, 1.43,
  "ce", 2, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200109plc1p1.csv", 205, 1.43,
  "ce", 3, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200113plc1p3.csv", 195, 1.44,
  "ce", 4, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200113plc1p2.csv", 205, 1.44,
  "ce", 5, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200322plc1p2.csv", 195, 1.44,
  "ce", 6, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200323plc1p1.csv", 185, 1.44,
  "ce", 7, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200326plc1p3.csv", 220, 1.44,
  "ce", 8, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. elegans/CD200326plc1p4.csv", 195, 1.44,
  "cb", 1, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD240731cbhis72p1.csv", 165, 1.57,
  "cb", 2, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD240731cbhis72p2.csv", 175, 1.57,
  "cb", 3, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD240731cbhis72p3.csv", 170, 1.57,
  "cb", 4, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD241202cbhis72p1.csv", 160, 1.58,
  "cb", 5, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD241202cbhis72p2.csv", 165, 1.58,
  "cb", 6, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. briggsae(she1)/CD241202cbhis72p4.csv", 180, 1.58,
  "cn", 1, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD241202cnhis72p1.csv", 235, 1.58,
  "cn", 2, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD240712cnhis72p1.csv", 235, 1.60,
  "cn", 3, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD240712cnhis72p2.csv", 235, 1.60,
  "cn", 4, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD240712cnhis72p3.csv", 235, 1.60,
  "cn", 5, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD241207cnhis72p1.csv", 230, 1.65,
  "cn", 6, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD241207cnhis72p3.csv", 230, 1.65,
  "cn", 7, "~/Desktop/HybridAnalysis/cbcn/CDFile/C. nigoni/CD241202cnhis72p2.csv", 200, 1.58
) %>%
  mutate(
    embryo = paste0(group, "p", id),
    file_path = path.expand(file_path),
    len_col = paste0("Len_", embryo),
    tim_col = paste0("Tim_", embryo)
  )

stage_levels <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                  "Stage_AB64", "Stage_AB128", "Stage_AB256")
start_cells <- c("ABa", "ABp", "EMS", "P2")
landmark_percentile <- 0.90

q90_time <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  unname(quantile(x, probs = landmark_percentile, na.rm = TRUE, type = 7))
}

absolute_data <- read_csv(absolute_file, show_col_types = FALSE) %>%
  mutate(stage = factor(stage, levels = stage_levels)) %>%
  filter(!is.na(stage))

get_start_time <- function(file_path, time_limit, multiplier) {
  cd <- read_csv(file_path, show_col_types = FALSE) %>%
    filter(time <= time_limit, cell %in% start_cells) %>%
    group_by(cell) %>%
    summarise(first_tp = min(time, na.rm = TRUE), .groups = "drop")
  missing_cells <- setdiff(start_cells, cd$cell)
  if (length(missing_cells) > 0) {
    warning("Missing start cells in ", basename(file_path), ": ", paste(missing_cells, collapse = ", "))
    return(NA_real_)
  }
  max(cd$first_tp, na.rm = TRUE) * multiplier
}

landmarks <- pmap_dfr(datasets_config, function(group, id, file_path, time_limit, multiplier, embryo, len_col, tim_col) {
  if (!tim_col %in% names(absolute_data)) {
    warning("Missing time column: ", tim_col)
    return(NULL)
  }
  stage_end <- absolute_data %>%
    group_by(stage) %>%
    summarise(time_min = q90_time(.data[[tim_col]]), .groups = "drop") %>%
    mutate(landmark = as.character(stage)) %>%
    select(landmark, time_min)

  bind_rows(
    tibble(landmark = "t_start_AB4_detected", time_min = get_start_time(file_path, time_limit, multiplier)),
    stage_end
  ) %>%
    mutate(group = group, id = id, embryo = embryo, .before = 1)
})

ce_reference <- landmarks %>%
  filter(group == "ce") %>%
  group_by(landmark) %>%
  summarise(ref_time_min = mean(time_min, na.rm = TRUE), .groups = "drop")

landmarks_with_ref <- landmarks %>%
  left_join(ce_reference, by = "landmark") %>%
  arrange(group, id, match(landmark, c("t_start_AB4_detected", stage_levels)))

write_csv(landmarks_with_ref, file.path(out_dir, "tables", "Table1_embryo_landmarks_and_CE_reference.csv"))

get_landmark_time <- function(embryo, landmark) {
  landmarks_with_ref$time_min[landmarks_with_ref$embryo == embryo & landmarks_with_ref$landmark == landmark][1]
}

get_ref_time <- function(landmark) {
  ce_reference$ref_time_min[ce_reference$landmark == landmark][1]
}

normalize_global <- function(df, config) {
  embryo <- config$embryo
  len_col <- config$len_col
  tim_col <- config$tim_col
  emb_start <- get_landmark_time(embryo, "t_start_AB4_detected")
  emb_end <- get_landmark_time(embryo, "Stage_AB256")
  ref_start <- get_ref_time("t_start_AB4_detected")
  ref_end <- get_ref_time("Stage_AB256")
  scale_to_ref <- (ref_end - ref_start) / (emb_end - emb_start)

  tibble(
    cell = df$cell,
    !!paste0("M1_Len_", embryo) := df[[len_col]] * scale_to_ref,
    !!paste0("M1_Tim_", embryo) := ref_start + (df[[tim_col]] - emb_start) * scale_to_ref
  )
}

normalize_piecewise_one <- function(times, lengths, stages, embryo) {
  landmark_order <- c("t_start_AB4_detected", stage_levels)
  emb_times <- landmarks_with_ref %>%
    filter(embryo == !!embryo, landmark %in% landmark_order) %>%
    arrange(match(landmark, landmark_order)) %>%
    pull(time_min)
  ref_times <- ce_reference %>%
    filter(landmark %in% landmark_order) %>%
    arrange(match(landmark, landmark_order)) %>%
    pull(ref_time_min)

  out_time <- rep(NA_real_, length(times))
  out_len <- rep(NA_real_, length(lengths))

  for (i in seq_along(times)) {
    raw_time <- times[i]
    if (!is.finite(raw_time)) next
    idx <- findInterval(raw_time, emb_times, rightmost.closed = TRUE)
    idx <- max(1L, min(idx, length(emb_times) - 1L))
    emb_left <- emb_times[idx]
    emb_right <- emb_times[idx + 1L]
    ref_left <- ref_times[idx]
    ref_right <- ref_times[idx + 1L]
    if (!all(is.finite(c(emb_left, emb_right, ref_left, ref_right))) || emb_right == emb_left) next
    scale_to_ref <- (ref_right - ref_left) / (emb_right - emb_left)
    out_time[i] <- ref_left + (raw_time - emb_left) * scale_to_ref
    out_len[i] <- lengths[i] * scale_to_ref
  }

  tibble(norm_length = out_len, norm_time = out_time)
}

normalize_piecewise <- function(df, config) {
  embryo <- config$embryo
  len_col <- config$len_col
  tim_col <- config$tim_col
  norm <- normalize_piecewise_one(df[[tim_col]], df[[len_col]], df$stage, embryo)
  tibble(
    cell = df$cell,
    !!paste0("M2_Len_", embryo) := norm$norm_length,
    !!paste0("M2_Tim_", embryo) := norm$norm_time
  )
}

normalize_monotonic_one <- function(times, lengths, embryo) {
  landmark_order <- c("t_start_AB4_detected", stage_levels)
  emb_times <- landmarks_with_ref %>%
    filter(embryo == !!embryo, landmark %in% landmark_order) %>%
    arrange(match(landmark, landmark_order)) %>%
    pull(time_min)
  ref_times <- ce_reference %>%
    filter(landmark %in% landmark_order) %>%
    arrange(match(landmark, landmark_order)) %>%
    pull(ref_time_min)

  keep <- is.finite(emb_times) & is.finite(ref_times)
  emb_times <- emb_times[keep]
  ref_times <- ref_times[keep]
  if (length(emb_times) < 3 || any(diff(emb_times) <= 0) || any(diff(ref_times) <= 0)) {
    return(normalize_piecewise_one(times, lengths, NULL, embryo))
  }

  warp_fun <- splinefun(x = emb_times, y = ref_times, method = "monoH.FC")
  out_time <- warp_fun(times)
  local_scale <- warp_fun(times, deriv = 1)
  local_scale[!is.finite(local_scale) | local_scale <= 0] <- NA_real_
  tibble(norm_length = lengths * local_scale, norm_time = out_time)
}

normalize_monotonic <- function(df, config) {
  embryo <- config$embryo
  len_col <- config$len_col
  tim_col <- config$tim_col
  norm <- normalize_monotonic_one(df[[tim_col]], df[[len_col]], embryo)
  tibble(
    cell = df$cell,
    !!paste0("M3_Len_", embryo) := norm$norm_length,
    !!paste0("M3_Tim_", embryo) := norm$norm_time
  )
}

method1 <- datasets_config %>%
  split(.$embryo) %>%
  map(~normalize_global(absolute_data, .x)) %>%
  reduce(full_join, by = "cell") %>%
  left_join(absolute_data %>% select(cell, stage), by = "cell")

method2 <- datasets_config %>%
  split(.$embryo) %>%
  map(~normalize_piecewise(absolute_data, .x)) %>%
  reduce(full_join, by = "cell") %>%
  left_join(absolute_data %>% select(cell, stage), by = "cell")

method3 <- datasets_config %>%
  split(.$embryo) %>%
  map(~normalize_monotonic(absolute_data, .x)) %>%
  reduce(full_join, by = "cell") %>%
  left_join(absolute_data %>% select(cell, stage), by = "cell")

method4 <- absolute_data %>%
  select(cell, matches("^Len_"), matches("^Tim_"), stage)

scale_global <- datasets_config %>%
  transmute(group, id, embryo) %>%
  rowwise() %>%
  mutate(
    embryo_start = get_landmark_time(embryo, "t_start_AB4_detected"),
    embryo_end_AB256_90pct = get_landmark_time(embryo, "Stage_AB256"),
    reference_start = get_ref_time("t_start_AB4_detected"),
    reference_end_AB256_90pct = get_ref_time("Stage_AB256"),
    scale_to_CE_reference = (reference_end_AB256_90pct - reference_start) / (embryo_end_AB256_90pct - embryo_start)
  ) %>%
  ungroup()

scale_piecewise <- datasets_config %>%
  select(group, id, embryo) %>%
  crossing(interval_index = seq_len(length(stage_levels))) %>%
  mutate(
    left_landmark = c("t_start_AB4_detected", stage_levels[-length(stage_levels)])[interval_index],
    right_landmark = stage_levels[interval_index]
  ) %>%
  rowwise() %>%
  mutate(
    embryo_left = get_landmark_time(embryo, left_landmark),
    embryo_right = get_landmark_time(embryo, right_landmark),
    reference_left = get_ref_time(left_landmark),
    reference_right = get_ref_time(right_landmark),
    scale_to_CE_reference = (reference_right - reference_left) / (embryo_right - embryo_left)
  ) %>%
  ungroup()

write_csv(method1, file.path(out_dir, "tables", "Method1_global_AB4_to_AB256_landmark_normalized_Length_Time.csv"))
write_csv(method2, file.path(out_dir, "tables", "Method2_stagewise_landmark_normalized_Length_Time.csv"))
write_csv(method3, file.path(out_dir, "tables", "Method3_monotonic_landmark_timewarp_normalized_Length_Time.csv"))
write_csv(method4, file.path(out_dir, "tables", "Method4_raw_absolute_Length_Time_no_normalization.csv"))
write_csv(scale_global, file.path(out_dir, "tables", "Table2_method1_global_scale_factors.csv"))
write_csv(scale_piecewise, file.path(out_dir, "tables", "Table3_method2_stagewise_scale_factors.csv"))

methods <- c(
  "# Landmark normalization methods",
  "",
  "This folder contains four normalization-related outputs for CE, CB, and CN embryos. Methods 1-3 use the existing absolute length/time table from S1.R and the original CD files only to determine the AB4 start landmark. Method 4 keeps the raw absolute values as a normalization-independent control.",
  "",
  "Start landmark: for each embryo, t_start is the earliest time point at which ABa, ABp, EMS, and P2 are all detected. In calculation, this is max(first_time_ABa, first_time_ABp, first_time_EMS, first_time_P2), converted to minutes by the embryo-specific time resolution.",
  "",
  "Stage landmark definition: for each embryo and each AB stage, the landmark time is the 90th-percentile division time among cells assigned to that stage. This avoids making the normalization depend on the single last-dividing cell.",
  "",
  "CE reference: each landmark reference time is the mean of the corresponding embryo-level landmark times across the eight CE embryos.",
  "",
  "Method 1: global AB4-to-AB256 landmark normalization. The embryo-specific interval from t_start to the 90th-percentile division time among Stage_AB256 cells is linearly scaled to the CE reference interval. Normalized time = reference_start + (raw_time - embryo_start) * (reference_end - reference_start) / (embryo_end - embryo_start). Normalized length = raw_length * (reference_end - reference_start) / (embryo_end - embryo_start).",
  "",
  "Method 2: stage-wise landmark normalization. Landmarks are t_start plus the 90th-percentile stage times for Stage_AB4, Stage_AB8, Stage_AB16, Stage_AB32, Stage_AB64, Stage_AB128, and Stage_AB256. For each raw division time, the interval containing that time is identified and scaled to the corresponding CE reference interval. Cell-cycle length is scaled by the same interval-specific scale factor used for that cell's division time.",
  "",
  "Method 3: monotonic landmark time-warp normalization. The same 90th-percentile stage landmarks are used to fit a monotone Hermite spline mapping embryo time to CE-reference time. Normalized division time is the spline-mapped time. Cell-cycle length is multiplied by the local derivative of the monotonic time-warp function at that cell's raw division time.",
  "",
  "Method 4: raw absolute values without normalization. The output table contains the original absolute cell-cycle length and division-time values from S1.R and is intended for normalization-independent statistical comparison. It is plotted beside the 90th-percentile landmark-normalized methods as the raw control.",
  "",
  "No figures are generated by this script."
)
writeLines(methods, file.path(out_dir, "methods.md"))

readme <- c(
  "# CE/CB/CN landmark normalization results",
  "",
  "tables/Method1_global_AB4_to_AB256_landmark_normalized_Length_Time.csv: normalized length/time using one global AB4-to-90th-percentile-AB256 scale factor per embryo.",
  "tables/Method2_stagewise_landmark_normalized_Length_Time.csv: normalized length/time using stage-wise 90th-percentile landmark scale factors per embryo.",
  "tables/Method3_monotonic_landmark_timewarp_normalized_Length_Time.csv: normalized length/time using monotonic 90th-percentile landmark time-warping.",
  "tables/Method4_raw_absolute_Length_Time_no_normalization.csv: raw absolute length/time values for normalization-independent analysis.",
  "tables/Table1_embryo_landmarks_and_CE_reference.csv: embryo 90th-percentile stage landmarks and CE reference landmark means.",
  "tables/Table2_method1_global_scale_factors.csv: global scale factors.",
  "tables/Table3_method2_stagewise_scale_factors.csv: interval-specific scale factors.",
  "code/cbcn_landmark_normalization_methods.R: reproducible R code."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Completed: ", out_dir)
