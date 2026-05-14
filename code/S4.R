# ============================
# cell position normalization template selection
# ============================

# Set working directory and load required packages
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(ggsci)

# Preprocessing function with z_factor scaling
process_dis <- function(data, time_limit, intercept, slope, z_factor) {
  data %>% 
    select(2, 3, 9, 10, 11) %>%
    mutate(z = z * z_factor) %>%
    group_by(cell) %>%
    as_tibble() %>%
    filter(time <= time_limit) %>%
    select(1, 2, 4, 5, 3) %>%
    mutate(nortime = round((time - intercept) / slope, 2))
}

# ============================
# 1. CE data processing
# ============================
ce1 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/ce/CD191108plc1p1.csv", header = TRUE), 205, 7.1, 1, 4.67)
ce2 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/ce/CD200109plc1p1.csv", header = TRUE), 205, -9.57, 1.05, 4.67)
ce3 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/ce/CD200113plc1p3.csv", header = TRUE), 195, -3.95, 1, 4.67)
ce4 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/ce/CD200113plc1p2.csv", header = TRUE), 200, 7.04, 1.04, 4.67)
ce5 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/ce/CD200322plc1p2.csv", header = TRUE), 195, 1.79, 0.95, 4.67)
ce6 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/ce/CD200323plc1p1.csv", header = TRUE), 185, -7.77, 0.95, 4.67)
ce7 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/ce/CD200326plc1p3.csv", header = TRUE), 220, 4.46, 1.04, 4.67)
ce8 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/ce/CD200326plc1p4.csv", header = TRUE), 195, 0.89, 0.98, 4.67)
ce_list <- list(ce1 = ce1, ce2 = ce2, ce3 = ce3, ce4 = ce4, ce5 = ce5, ce6 = ce6, ce7 = ce7, ce8 = ce8)

# ============================
# 2. CB data processing
# ============================
cb1 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p1.csv", header = TRUE), 165, 6.62, 0.91, 4.78)
cb2 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p2.csv", header = TRUE), 175, 8.17, 0.92, 4.78)
cb3 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p3.csv", header = TRUE), 170, 9.41, 0.9, 4.78)
cb4 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p1.csv", header = TRUE), 160, 8.57, 0.87, 4.78)
cb5 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p2.csv", header = TRUE), 165, 22.79, 0.88, 4.78)
cb6 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p4.csv", header = TRUE), 180, 25.75, 0.93, 4.78)
cb_list <- list(cb1 = cb1, cb2 = cb2, cb3 = cb3, cb4 = cb4, cb5 = cb5, cb6 = cb6)

# ============================
# 3. CN data processing
# ============================
cn1 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/cn/CD241202cnhis72p1.csv", header = TRUE), 235, 2.92, 1.35, 4.78)
cn2 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p1.csv", header = TRUE), 235, 16.02, 1.36, 4.78)
cn3 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p2.csv", header = TRUE), 235, 18.43, 1.34, 4.78)
cn4 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p3.csv", header = TRUE), 235, 12.5, 1.38, 4.78)
cn5 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p1.csv", header = TRUE), 230, 15.57, 1.33, 4.78)
cn6 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p3.csv", header = TRUE), 230, 34.95, 1.29, 4.78)
cn7 <- process_dis(read.csv("~/Desktop/cbcn/CDFile/cn/CD241202cnhis72p2.csv", header = TRUE), 200, 6.26, 1.17, 4.78)
cn_list <- list(cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cn5 = cn5, cn6 = cn6, cn7 = cn7)

# ============================
# Generic PCA and plotting utility
# ============================
run_pca <- function(data_list, title_sub = "Group-wise PCA distribution") {
  # Combine all data and keep x, y, z columns
  all_data <- do.call(rbind, lapply(data_list, function(df) df[, c("x", "y", "z")]))
  pca_result <- prcomp(all_data, scale. = FALSE)
  pca_scores <- pca_result$x[, 1:2]
  
  group_sizes <- sapply(data_list, nrow)
  group_labels <- rep(names(data_list), group_sizes)
  
  # Compute group centroids in PCA space
  group_centers <- matrix(NA, nrow = length(data_list), ncol = 2)
  names_indices <- c(0, cumsum(group_sizes))
  for (i in seq_along(data_list)) {
    idx <- (names_indices[i] + 1):names_indices[i + 1]
    group_centers[i, ] <- colMeans(pca_scores[idx, , drop = FALSE])
  }
  
  # Compute global centroid and distances
  global_center <- colMeans(pca_scores)
  euclidean_distance <- function(vec1, vec2) sqrt(sum((vec1 - vec2)^2))
  group_distances <- apply(group_centers, 1, function(center) {
    euclidean_distance(center, global_center)
  })
  
  # Pick the template group as the one closest to the global centroid
  template_group <- names(data_list)[which.min(group_distances)]
  
  # Print summary metrics
  cat("Group centroids:\n")
  print(group_centers)
  cat("Global PCA centroid:", global_center, "\n")
  cat("Distances to global centroid:\n")
  print(group_distances)
  cat("Selected template group:", template_group, "\n")
  
  # Visualization
  pca_df <- data.frame(PC1 = pca_scores[, 1],
                       PC2 = pca_scores[, 2],
                       group = group_labels)
  group_centers_df <- data.frame(PC1 = group_centers[, 1],
                                 PC2 = group_centers[, 2],
                                 group = names(data_list))
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(alpha = 0.6) +
    geom_point(data = group_centers_df, aes(x = PC1, y = PC2, color = group),
               shape = 4, size = 4, stroke = 2) +
    annotate("point", x = global_center[1], y = global_center[2],
             color = "black", size = 5, shape = 8) +
    geom_text(data = group_centers_df, aes(label = group), vjust = -1, size = 5, color = "black") +
    labs(title = "PCA of Combined Raw Data",
         subtitle = paste(title_sub, "- Template:", template_group),
         x = "PC1", y = "PC2") +
    theme_minimal()
}

# ============================
# Run PCA analyses and plot
# ============================
p_ce <- run_pca(ce_list, title_sub = "CE data")
p_cb <- run_pca(cb_list, title_sub = "CB data")
p_cn <- run_pca(cn_list, title_sub = "CN data")

print(p_ce)
print(p_cb)
print(p_cn)






# ---------------------------
# Position Normalization based on the selected template
# ---------------------------
# ---------------------------
# Basic setup and packages
# ---------------------------
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(ggpubr)
library(reshape2)
library(RColorBrewer)

# ---------------------------
# Utility functions
# ---------------------------

# 1. Data preprocessing with z_factor support for custom scaling
process_dis <- function(data, time_limit, intercept, slope, z_factor) {
  data %>%
    select(2, 3, 9, 10, 11) %>%
    mutate(z = z * z_factor) %>%
    group_by(cell) %>%
    as_tibble() %>%
    filter(time <= time_limit) %>%
    select(1, 2, 4, 5, 3) %>%
    mutate(nortime = round((time - intercept) / slope, 2))
}

# 2. Keep only cells common to all datasets and remove unwanted cells
filter_common_cells <- function(datasets, remove_cells) {
  common <- Reduce(intersect, lapply(datasets, function(df) unique(df$cell)))
  lapply(datasets, function(df) {
    df %>% filter(cell %in% common, !(cell %in% remove_cells))
  })
}

# 3. Normalize time within each cell and compute time_percent
normalize_time <- function(df) {
  df %>% 
    group_by(cell) %>%
    mutate(time_percent = (nortime - min(nortime)) / (max(nortime) - min(nortime))) %>%
    ungroup()
}

# 4. Align target datasets to the template by matching closest time_percent and carrying over template nortime as con_time
align_to_template <- function(df, template) {
  do.call(rbind, lapply(unique(template$cell), function(cell_val) {
    tmpl_cell <- filter(template, cell == cell_val)
    df_cell <- filter(df, cell == cell_val)
    do.call(rbind, lapply(1:nrow(tmpl_cell), function(i) {
      tmpl_tp <- tmpl_cell$time_percent[i]
      tmpl_nortime <- tmpl_cell$nortime[i]
      if(nrow(df_cell) == 0) return(NULL)
      idx <- which.min(abs(df_cell$time_percent - tmpl_tp))
      row_selected <- df_cell[idx, ]
      row_selected$con_time <- tmpl_nortime
      row_selected
    }))
  }))
}

# 5. General pipeline: read files, preprocess, normalize, align to template, convert con_time
process_group <- function(filelist, params, z_factor, remove_cells, template_name) {
  datasets <- mapply(function(f, p) {
    process_dis(read.csv(f, header = TRUE), p$time_limit, p$intercept, p$slope, z_factor)
  }, filelist, params, SIMPLIFY = FALSE)
  
  # Optionally enforce common cells (currently disabled)
  #datasets <- filter_common_cells(datasets, remove_cells)
  datasets <- lapply(datasets, normalize_time)
  
  template <- datasets[[template_name]]
  
  aligned <- lapply(names(datasets), function(name) {
    if(name == template_name) {
      datasets[[name]] %>% mutate(con_time = nortime)
    } else {
      align_to_template(datasets[[name]], template)
    }
  })
  names(aligned) <- names(datasets)
  
  # Convert con_time to numeric (keep decimals if needed)
  lapply(aligned, function(df) df %>% mutate(con_time = as.numeric(con_time)))
}

# 6. Check alignment gaps by comparing template con_time values against targets
find_missing_alignment <- function(template_df, target_df) {
  do.call(rbind, lapply(unique(template_df$cell), function(cell_val) {
    tmpl <- filter(template_df, cell == cell_val)
    tgt  <- filter(target_df, cell == cell_val)
    missing <- tmpl %>% filter(!(con_time %in% tgt$con_time))
    if(nrow(missing) > 0) missing %>% mutate(cell = cell_val) else NULL
  }))
}

# ---------------------------
# CE group processing
# ---------------------------
ce_files <- list(
  ce1 = "~/Desktop/cbcn/CDFile/ce/CD191108plc1p1.csv",
  ce2 = "~/Desktop/cbcn/CDFile/ce/CD200109plc1p1.csv",
  ce3 = "~/Desktop/cbcn/CDFile/ce/CD200113plc1p3.csv",
  ce4 = "~/Desktop/cbcn/CDFile/ce/CD200113plc1p2.csv",
  ce5 = "~/Desktop/cbcn/CDFile/ce/CD200322plc1p2.csv",
  ce6 = "~/Desktop/cbcn/CDFile/ce/CD200323plc1p1.csv",
  ce7 = "~/Desktop/cbcn/CDFile/ce/CD200326plc1p3.csv",
  ce8 = "~/Desktop/cbcn/CDFile/ce/CD200326plc1p4.csv"
)
ce_params <- list(
  ce1 = list(time_limit = 205, intercept = 7.1,  slope = 1),
  ce2 = list(time_limit = 205, intercept = -9.57, slope = 1.05),
  ce3 = list(time_limit = 195, intercept = -3.95, slope = 1),
  ce4 = list(time_limit = 200, intercept = 7.04, slope = 1.04),
  ce5 = list(time_limit = 195, intercept = 1.79,  slope = 0.95),
  ce6 = list(time_limit = 185, intercept = -7.77, slope = 0.95),
  ce7 = list(time_limit = 220, intercept = 4.46,  slope = 1.04),
  ce8 = list(time_limit = 195, intercept = 0.89,  slope = 0.98)
)
ce_remove <- c("ABa", "ABp", "EMS", "P2", "Cppappa", "Cppappp", "ABalapaaaap", "ABalapaaaaa")
# Align CE group to template ce4
ce_aligned <- process_group(ce_files, ce_params, 4.67, ce_remove, "ce4")

# ---------------------------
# CB group processing
# ---------------------------
cb_files <- list(
  cb1 = "~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p1.csv",
  cb2 = "~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p2.csv",
  cb3 = "~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p3.csv",
  cb4 = "~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p1.csv",
  cb5 = "~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p2.csv",
  cb6 = "~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p4.csv"
)
cb_params <- list(
  cb1 = list(time_limit = 165, intercept = 6.62, slope = 0.91),
  cb2 = list(time_limit = 175, intercept = 8.17, slope = 0.92),
  cb3 = list(time_limit = 170, intercept = 9.41, slope = 0.9),
  cb4 = list(time_limit = 160, intercept = 8.57, slope = 0.87),
  cb5 = list(time_limit = 165, intercept = 22.79, slope = 0.88),
  cb6 = list(time_limit = 180, intercept = 25.75, slope = 0.93)
)
cb_remove <- c("ABa", "ABp", "EMS", "P2", "ABaraaaapaa", "ABaraaaapap")
# Align CB group to template cb4
cb_aligned <- process_group(cb_files, cb_params, 4.78, cb_remove, "cb4")

# ---------------------------
# CN group processing
# ---------------------------
cn_files <- list(
  cn1 = "~/Desktop/cbcn/CDFile/cn/CD241202cnhis72p1.csv",
  cn2 = "~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p1.csv",
  cn3 = "~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p2.csv",
  cn4 = "~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p3.csv",
  cn5 = "~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p1.csv",
  cn6 = "~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p3.csv",
  cn7 = "~/Desktop/cbcn/CDFile/cn/CD241202cnhis72p2.csv"
)
cn_params <- list(
  cn1 = list(time_limit = 235, intercept = 2.92,  slope = 1.35),
  cn2 = list(time_limit = 235, intercept = 16.02, slope = 1.36),
  cn3 = list(time_limit = 235, intercept = 18.43, slope = 1.34),
  cn4 = list(time_limit = 235, intercept = 12.5,  slope = 1.38),
  cn5 = list(time_limit = 230, intercept = 15.57, slope = 1.33),
  cn6 = list(time_limit = 230, intercept = 34.95, slope = 1.29),
  cn7 = list(time_limit = 200, intercept = 6.26,  slope = 1.17)
)
#cn_remove <- c("ABa", "ABp", "EMS", "P2", "ABalaapaapa", "ABalaapaapp", "Capaaaa", "Capaaap")
cn_remove <- c("ABa", "ABp", "EMS", "P2")

# Align CN group to template cn5
cn_aligned <- process_group(cn_files, cn_params, 4.78, cn_remove, "cn5")

# ---------------------------
# Save intermediate objects
# ---------------------------
save(ce_aligned, cb_aligned, cn_aligned, file = "normalized_position.RData")

# ---------------------------
# Alignment sanity check (CB group example: compare each dataset against template cb4)
# ---------------------------
alignment_issues <- lapply(names(cb_aligned), function(name) {
  if(name != "cb4") {
    list(dataset = name, missing = find_missing_alignment(cb_aligned[["cb4"]], cb_aligned[[name]]))
  } else {
    NULL
  }
})
alignment_issues <- alignment_issues[!sapply(alignment_issues, is.null)]
for(issue in alignment_issues) {
  cat("Dataset", issue$dataset, "missing template alignments:\n")
  print(issue$missing)
  cat("\n-------------------------------------\n")
}






# ============================
# embryo length calculation
# ============================
# ------------------------------------------------------------
# 0. 基本设置
# ------------------------------------------------------------
setwd("~/Desktop/cbcn/draft_code/new/")
library(dplyr)
library(readr)

# ------------------------------------------------------------
# 1. 读取所有数据（保持原有文件路径）
# ------------------------------------------------------------
ce_files <- list(
  ce1 = "~/Desktop/cbcn/CDFile/ce/CD191108plc1p1.csv",
  ce2 = "~/Desktop/cbcn/CDFile/ce/CD200109plc1p1.csv",
  ce3 = "~/Desktop/cbcn/CDFile/ce/CD200113plc1p3.csv",
  ce4 = "~/Desktop/cbcn/CDFile/ce/CD200113plc1p2.csv",
  ce5 = "~/Desktop/cbcn/CDFile/ce/CD200322plc1p2.csv",
  ce6 = "~/Desktop/cbcn/CDFile/ce/CD200323plc1p1.csv",
  ce7 = "~/Desktop/cbcn/CDFile/ce/CD200326plc1p3.csv",
  ce8 = "~/Desktop/cbcn/CDFile/ce/CD200326plc1p4.csv"
)
cb_files <- list(
  cb1 = "~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p1.csv",
  cb2 = "~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p2.csv",
  cb3 = "~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p3.csv",
  cb4 = "~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p1.csv",
  cb5 = "~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p2.csv",
  cb6 = "~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p4.csv"
)
cn_files <- list(
  cn1 = "~/Desktop/cbcn/CDFile/cn/CD241202cnhis72p1.csv",
  cn2 = "~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p1.csv",
  cn3 = "~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p2.csv",
  cn4 = "~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p3.csv",
  cn5 = "~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p1.csv",
  cn6 = "~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p3.csv",
  cn7 = "~/Desktop/cbcn/CDFile/cn/CD241202cnhis72p2.csv"
)

all_files <- c(ce_files, cb_files, cn_files)

# ------------------------------------------------------------
# 2. 统一的预处理函数（与“后一脚本”一致）
# ------------------------------------------------------------
process_dis <- function(data, z_factor = 4.78, time_limit) {
  data %>%
    select(2, 3, 9, 10, 11) %>%    # 保留 cell, time, x, y, z
    mutate(z = z * z_factor) %>%
    group_by(cell) %>%
    as_tibble() %>%
    filter(time <= time_limit) %>%
    select(cell, time, x, y, z)
}

# 针对每个文件给出 time_limit（可沿用原表）
time_limits <- list(
  ce1 = 205, ce2 = 205, ce3 = 195, ce4 = 205, ce5 = 195, ce6 = 185, ce7 = 220, ce8 = 195,
  cb1 = 165, cb2 = 175, cb3 = 170, cb4 = 160, cb5 = 165, cb6 = 180,
  cn1 = 235, cn2 = 235, cn3 = 235, cn4 = 235, cn5 = 230, cn6 = 230, cn7 = 200
)

processed_data <- lapply(names(all_files), function(name) {
  df <- read_csv(all_files[[name]])
  process_dis(df, z_factor = 4.78, time_limit = time_limits[[name]])
})
names(processed_data) <- names(all_files)

# ------------------------------------------------------------
# 3. 胚胎长度计算（完全照搬示例脚本逻辑）
# ------------------------------------------------------------
compute_embryo_length <- function(df, embryo_name) {
  df_filtered <- df %>%
    group_by(cell) %>%
    slice_max(order_by = time, n = 2, with_ties = FALSE) %>%
    ungroup()
  
  xyz_mat <- as.matrix(df_filtered[, c("x", "y", "z")])
  dists <- as.vector(dist(xyz_mat, method = "euclidean"))
  
  n_top <- min(10, length(dists))
  top_vals <- sort(dists, decreasing = TRUE)[seq_len(n_top)]
  embryo_length <- median(top_vals)
  
  tibble(embryo = embryo_name, embryo_length = embryo_length)
}

length_results <- bind_rows(lapply(names(processed_data), function(name) {
  compute_embryo_length(processed_data[[name]], name)
}))

write_tsv(length_results, "embryo_lengths.txt")




# ============================
# cell cell distance calculation
# ============================
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

# 加载所需包
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(reshape2)
library(RColorBrewer)

# 加载数据（其中 ce_aligned, cb_aligned, cn_aligned 均已存在）
load("normalized_position.RData")
length_all_results <- read.table("~/Desktop/cbcn/draft_code/new/embryo_lengths.txt", header = TRUE, sep = "\t")

summary(ce_aligned)
head(ce_aligned[[1]], n = 5)
str(ce_aligned)

# 定义函数：针对每个数据，计算每个 con_time 下 cell 之间的欧氏距离矩阵，
# 并用当前胚胎的胚胎长度进行归一化
compute_distance_by_time <- function(df, embryo_name) {
  # 从 length_all_results 中获取当前胚胎的胚胎长度
  embryo_length_current <- length_all_results$embryo_length[
    length_all_results$embryo == embryo_name
  ]
  
  if (length(embryo_length_current) == 0) {
    warning(paste("没有找到胚胎", embryo_name, "对应的胚胎长度，使用1作为默认值"))
    embryo_length_current <- 1
  }
  
  # 获取所有非重复的 con_time 值
  time_points <- unique(df$con_time)
  
  # 用于存储每个时间点对应的距离矩阵
  distances_list <- list()
  
  # 针对每个时间点分组处理
  for (t in time_points) {
    # 过滤出当前时间点的数据
    sub_df <- df %>% filter(con_time == t)
    
    # 选取空间坐标，构造三维坐标矩阵
    coords <- as.matrix(sub_df %>% select(x, y, z))
    
    # 计算欧氏距离矩阵（行和列均对应每个 cell）
    dmat <- as.matrix(dist(coords, method = "euclidean"))
    
    # 使用当前胚胎的胚胎长度做归一化
    dmat <- dmat / embryo_length_current
    
    # 设置距离矩阵的行名与列名为 cell 名称
    rownames(dmat) <- sub_df$cell
    colnames(dmat) <- sub_df$cell
    
    # 将当前 con_time 的归一化距离矩阵存入列表，key 为当前 time（字符形式）
    distances_list[[as.character(t)]] <- dmat
  }
  
  # 返回按时间分组的距离矩阵列表
  return(distances_list)
}

# --------------------------------------------------
# 对 ce、cb、cn 三类数据，根据各自的胚胎名称（列表内的 key）进行归一化距离计算

# 对 ce 类数据（例如 ce1 到 ce8）
distance_ce <- lapply(names(ce_aligned), function(embryo_name) {
  compute_distance_by_time(ce_aligned[[embryo_name]], embryo_name)
})
names(distance_ce) <- names(ce_aligned)

# 对 cb 类数据（例如 cb1 到 cb6）
distance_cb <- lapply(names(cb_aligned), function(embryo_name) {
  compute_distance_by_time(cb_aligned[[embryo_name]], embryo_name)
})
names(distance_cb) <- names(cb_aligned)

# 对 cn 类数据（例如 cn1 到 cn7）
distance_cn <- lapply(names(cn_aligned), function(embryo_name) {
  compute_distance_by_time(cn_aligned[[embryo_name]], embryo_name)
})
names(distance_cn) <- names(cn_aligned)

# --------------------------------------------------
summary(distance_ce)
names(distance_ce)
names(distance_ce$ce1)
save(distance_ce, distance_cb,distance_cn, 
     file = "distance_results.RData")


# ============================
# distance plot
# ============================
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

load("distance_results.RData")

###28cell, 170,300,509
### ce:42.27 (41.31) (40.35), 90.35 (91.31), 120.15, 165.35
### cb:33.83(30.38)(31.53)(32.68)(34.98), 84.4, 108.54(109.69), 149.92
### cn:25.14(25.89)(26.64)(27.39)(28.14)(28.89)(29.65)(30.4)(31.15)(31.9), 78.52, 104.08, 142.43(143.18)

summary(distance_ce)
names(distance_ce)
names(distance_ce$ce1)

# 假设 distance_ce 已经计算好，并且其结构为：
# distance_ce 是一个列表，键为 "ce1", "ce2", ... "ce8"
# 每个 distance_ce[[embryo]] 又是一个列表，键为时间点字符串

# 1. 定义需要提取的时间点
#desired_times <- c("42.27", "41.31", "40.35", "90.35", "91.31", "120.15", "165.35")
desired_times <- c("42.27", "90.35", "120.15", "165.35")

# 2. 提取出每个胚胎中指定时间点的距离矩阵
extracted_matrices <- list()  # 用于保存每个时间点各胚胎的矩阵
for (tm in desired_times) {
  mats <- list()  # 存储各胚胎在该时间点的距离矩阵
  for (emb in names(distance_ce)) {
    # 如果该胚胎的列表中包含时间点 tm，则提取出来
    if (tm %in% names(distance_ce[[emb]])) {
      mats[[emb]] <- distance_ce[[emb]][[tm]]
    } else {
      warning(paste("胚胎", emb, "缺少时间点", tm))
    }
  }
  if (length(mats) > 0) {
    extracted_matrices[[tm]] <- mats
  }
}

# 3. 对每个时间点计算所有胚胎（ce1到ce8）对应矩阵的逐元素平均值，
#    即保持矩阵的结构，计算所有矩阵元素对应位置的平均结果，
#    并将结果保存到 ce_average 中
ce_average <- list()
for (tm in names(extracted_matrices)) {
  mats_list <- extracted_matrices[[tm]]
  # 确保至少有一个矩阵
  if (length(mats_list) > 0) {
    # 使用 Reduce 对所有矩阵逐元素相加，然后除以胚胎数目得到平均矩阵
    avg_mat <- Reduce("+", mats_list) / length(mats_list)
    ce_average[[tm]] <- avg_mat
  }
}

# 输出结果，查看每个时间点的平均矩阵
for (tm in names(ce_average)) {
  cat("时间点:", tm, "\n")
  print(ce_average[[tm]])
  cat("\n")
}



### cb:33.83(30.38)(31.53)(32.68)(34.98), 84.4, 108.54(109.69), 149.92
# 1. 定义需要提取的时间点
desired_times <- c("33.83", "84.4", "108.54", "149.92")

# 2. 提取出每个胚胎中指定时间点的距离矩阵
extracted_matrices <- list()  # 用于保存每个时间点各胚胎的矩阵
for (tm in desired_times) {
  mats <- list()  # 存储各胚胎在该时间点的距离矩阵
  for (emb in names(distance_cb)) {
    # 如果该胚胎的列表中包含时间点 tm，则提取出来
    if (tm %in% names(distance_cb[[emb]])) {
      mats[[emb]] <- distance_cb[[emb]][[tm]]
    } else {
      warning(paste("胚胎", emb, "缺少时间点", tm))
    }
  }
  if (length(mats) > 0) {
    extracted_matrices[[tm]] <- mats
  }
}

# 3. 对每个时间点计算所有胚胎对应矩阵的逐元素平均值，
#    其他时间点使用 Reduce 直接计算均值，
#    对于时间点 "149.92" 使用特殊处理：
#    如果某个位置有6个重复值，则取6个的平均；如果只有5个，则按5个构成平均值。
cb_average <- list()
for (tm in names(extracted_matrices)) {
  mats_list <- extracted_matrices[[tm]]
  
  # 确保至少有一个矩阵
  if (length(mats_list) > 0) {
    if (tm != "149.92") {
      # 对于其他时间点，直接使用 Reduce 对所有矩阵逐元素求和后平均
      avg_mat <- Reduce("+", mats_list) / length(mats_list)
      cb_average[[tm]] <- avg_mat
    } else {
      # 对时间点 "149.92" 采用特殊处理
      # 先确定最终结果矩阵的尺寸：取所有矩阵中的最大行数和列数
      final_nrow <- max(sapply(mats_list, nrow))
      final_ncol <- max(sapply(mats_list, ncol))
      
      # 初始化叠加矩阵和计数矩阵
      sum_mat <- matrix(0, nrow = final_nrow, ncol = final_ncol)
      count_mat <- matrix(0, nrow = final_nrow, ncol = final_ncol)
      
      # 对每个矩阵进行累加，并统计每个位置有效值的数量
      for (mat in mats_list) {
        nr <- nrow(mat)
        nc <- ncol(mat)
        
        # 假定每个矩阵的数据放在左上角（例如 507x507 子矩阵对应于最终结果的前507行和507列）
        sum_mat[1:nr, 1:nc] <- sum_mat[1:nr, 1:nc] + mat
        count_mat[1:nr, 1:nc] <- count_mat[1:nr, 1:nc] + 1
      }
      
      # 计算平均值：对于每个位置用累加和除以该位置的计数
      avg_mat <- sum_mat
      avg_mat[count_mat > 0] <- sum_mat[count_mat > 0] / count_mat[count_mat > 0]
      avg_mat[count_mat == 0] <- NA
      
      cb_average[[tm]] <- avg_mat
    }
  }
}

# 4. 输出结果，查看每个时间点的平均矩阵
for (tm in names(cb_average)) {
  cat("时间点:", tm, "\n")
  print(cb_average[[tm]])
  cat("\n")
}



### cn:25.14(25.89)(26.64)(27.39)(28.14)(28.89)(29.65)(30.4)(31.15)(31.9), 78.52, 104.08, 142.43(143.18)

# 1. 定义 cn 的候选时间点组
desired_times <- c("25.14","78.52","104.08","143.18")
# 2. 提取出每个胚胎中指定时间点的距离矩阵
extracted_matrices <- list()  # 用于保存每个时间点各胚胎的矩阵
for (tm in desired_times) {
  mats <- list()  # 存储各胚胎在该时间点的距离矩阵
  for (emb in names(distance_cn)) {
    # 如果该胚胎的列表中包含时间点 tm，则提取出来
    if (tm %in% names(distance_cn[[emb]])) {
      mats[[emb]] <- distance_cn[[emb]][[tm]]
    } else {
      warning(paste("胚胎", emb, "缺少时间点", tm))
    }
  }
  if (length(mats) > 0) {
    extracted_matrices[[tm]] <- mats
  }
}

# 3. 对每个时间点计算所有胚胎对应矩阵的逐元素平均值，
#    其他时间点使用 Reduce 直接计算均值，
#    对于时间点 "149.92" 使用特殊处理：
#    如果某个位置有6个重复值，则取6个的平均；如果只有5个，则按5个构成平均值。
cn_average <- list()
for (tm in names(extracted_matrices)) {
  mats_list <- extracted_matrices[[tm]]
  
  # 确保至少有一个矩阵
  if (length(mats_list) > 0) {
    if (tm != "143.18") {
      # 对于其他时间点，直接使用 Reduce 对所有矩阵逐元素求和后平均
      avg_mat <- Reduce("+", mats_list) / length(mats_list)
      cn_average[[tm]] <- avg_mat
    } else {
      # 对时间点 "149.92" 采用特殊处理
      # 先确定最终结果矩阵的尺寸：取所有矩阵中的最大行数和列数
      final_nrow <- max(sapply(mats_list, nrow))
      final_ncol <- max(sapply(mats_list, ncol))
      
      # 初始化叠加矩阵和计数矩阵
      sum_mat <- matrix(0, nrow = final_nrow, ncol = final_ncol)
      count_mat <- matrix(0, nrow = final_nrow, ncol = final_ncol)
      
      # 对每个矩阵进行累加，并统计每个位置有效值的数量
      for (mat in mats_list) {
        nr <- nrow(mat)
        nc <- ncol(mat)
        
        # 假定每个矩阵的数据放在左上角（例如 507x507 子矩阵对应于最终结果的前507行和507列）
        sum_mat[1:nr, 1:nc] <- sum_mat[1:nr, 1:nc] + mat
        count_mat[1:nr, 1:nc] <- count_mat[1:nr, 1:nc] + 1
      }
      
      # 计算平均值：对于每个位置用累加和除以该位置的计数
      avg_mat <- sum_mat
      avg_mat[count_mat > 0] <- sum_mat[count_mat > 0] / count_mat[count_mat > 0]
      avg_mat[count_mat == 0] <- NA
      
      cn_average[[tm]] <- avg_mat
    }
  }
}

# 4. 输出结果，查看每个时间点的平均矩阵
for (tm in names(cn_average)) {
  cat("时间点:", tm, "\n")
  print(cn_average[[tm]])
  cat("\n")
}


# 加载所需的包
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggplotify)  # 用于将 grob 转换为 ggplot 对象
library(cowplot)    # 用于组合 ggplot 对象

# 定义注释颜色（与原代码一致）
ann_colors <- list(
  group = c(A = "#F8766D", M = "#A3A500",
            E = "#00BF7D", C = "#00B0F6",
            D = "#E76BF3")
)

# 辅助函数：根据矩阵生成热图，并转换为 ggplot 对象
# 为了让所有热图总体尺寸相同，我们先设定总体绘图区域的宽高（单位 mm）
# 然后根据矩阵的列数和行数计算每个 cell 的尺寸
plot_heatmap <- function(mat, title) {
  # 如果矩阵有行名，则生成行注释，否则设为空
  if (!is.null(rownames(mat))) {
    ann <- data.frame(row = rownames(mat)) %>% 
      mutate(group = substr(row, 1, 1)) %>%  # 取行名首字符进行分组
      column_to_rownames("row") %>% 
      mutate(group = factor(ifelse(group %in% c("P", "Z"), "D", group),
                            levels = c("A", "M", "E", "C", "D")))
  } else {
    ann <- NULL
  }
  
  # 设定“绘图区域”尺寸（仅指热图主体，不含外围 annotation）
  # 这里我们希望每个热图固定占 100 mm x 100 mm 的大小
  desired_width_mm <- 100
  desired_height_mm <- 100
  
  # 计算矩阵行数与列数
  n_cols <- ncol(mat)
  n_rows <- nrow(mat)
  
  # 动态计算每个 cell 的尺寸（单位：mm）
  # 这样不论矩阵多少列、多少行，热图主体的总尺寸都是 desired_width_mm x desired_height_mm
  cellwidth_val <- desired_width_mm / n_cols
  cellheight_val <- desired_height_mm / n_rows
  
  # 使用 pheatmap 绘制热图
  # 注意：这里设置 cellwidth 和 cellheight 为动态计算的数值
  p <- pheatmap(mat,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                show_colnames = FALSE,
                show_rownames = FALSE,
                annotation_row = ann,
                annotation_col = ann,
                cellwidth = cellwidth_val,
                cellheight = cellheight_val,
                annotation_colors = ann_colors,
                color = colorRampPalette(c("navy", "white", "firebrick3"))(20),
                main = title)
  
  # 将 pheatmap 返回的 gtable 转换为 ggplot 对象
  gg_hm <- as.ggplot(p$gtable)
  return(gg_hm)
}

# === 针对 ce_average 生成热图 ===
# 假设 ce_average 是一个列表，每个元素均为 nxn 数值矩阵
ce_plots <- lapply(names(ce_average), function(tm) {
  plot_heatmap(ce_average[[tm]], title = paste0("ce_average: ", tm))
})

# 使用 cowplot::plot_grid 组合时，指定相同的宽高比例，每个子图大小一致
ce_grid <- plot_grid(plotlist = ce_plots, nrow = 1)

# 保存组合后的图形，固定整个画布的宽高
ggsave("ce_average_heatmaps.pdf", ce_grid, width = 15, height = 5)

# === 针对 cb_average 生成热图 ===
cb_plots <- lapply(names(cb_average), function(tm) {
  plot_heatmap(cb_average[[tm]], title = paste0("cb_average: ", tm))
})
cb_grid <- plot_grid(plotlist = cb_plots, nrow = 1)
ggsave("cb_average_heatmaps.pdf", cb_grid, width = 15, height = 5)

# === 针对 cn_average 生成热图 ===
cn_plots <- lapply(names(cn_average), function(tm) {
  plot_heatmap(cn_average[[tm]], title = paste0("cn_average: ", tm))
})
cn_grid <- plot_grid(plotlist = cn_plots, nrow = 1)
ggsave("cn_average_heatmaps.pdf", cn_grid, width = 15, height = 5)




