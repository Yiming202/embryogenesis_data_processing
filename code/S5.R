
# ================================================================
# 0. 基本设置
# ================================================================
root_dir <- "~/Desktop/cbcn/"
setwd(root_dir)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(reshape2)
  library(RColorBrewer)
  library(ggsci)
  library(purrr)
})

# ================================================================
# 1. 阶段定义
# ================================================================
Stage_AB4 <- c("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")
Stage_AB8  <- c('ABala','ABalp','ABara','ABarp','ABpla','ABplp','ABpra','ABprp','MSa','MSp','Ca','Cp','P3')
Stage_AB16 <- c('ABalaa','ABalap','ABalpa','ABalpp','ABaraa','ABarap','ABarpa','ABarpp',
                'ABplaa','ABplap','ABplpa','ABplpp','ABpraa','ABprap','ABprpa','ABprpp',
                'MSaa','MSap','MSpa','MSpp','Ea','Ep','D')
Stage_AB32 <- c("ABalaaa", "ABalaap", "ABalapa", "ABalapp", "ABalpaa", "ABalpap", "ABalppa", "ABalppp",
                "ABaraaa", "ABaraap", "ABarapa", "ABarapp", "ABarpaa", "ABarpap", "ABarppa", "ABarppp",
                "ABplaaa", "ABplaap", "ABplapa", "ABplapp", "ABplpaa", "ABplpap", "ABplppa", "ABplppp",
                "ABpraaa", "ABpraap", "ABprapa", "ABprapp", "ABprpaa", "ABprpap", "ABprppa", "ABprppp",
                "MSaaa", "MSaap", "MSapa", "MSapp", "MSpaa", "MSpap", "MSppa", "MSppp",
                "Eal", "Ear", "Epl", "Epr", "Caa", "Cap", "Cpa", "Cpp", "Da", "Dp", "P4")
Stage_AB64 <- c("ABalaaaa", "ABalaaap", "ABalaapa", "ABalaapp", "ABalapaa", "ABalapap", "ABalappa", "ABalappp",
                "ABalpaaa", "ABalpaap", "ABalpapa", "ABalpapp", "ABalppaa", "ABalppap", "ABalpppa", "ABalpppp",
                "ABaraaaa", "ABaraaap", "ABaraapa", "ABaraapp", "ABarapaa", "ABarapap", "ABarappa", "ABarappp",
                "ABarpaaa", "ABarpaap", "ABarpapa", "ABarpapp", "ABarppaa", "ABarppap", "ABarpppa", "ABarpppp",
                "ABplaaaa", "ABplaaap", "ABplaapa", "ABplaapp", "ABplapaa", "ABplapap", "ABplappa", "ABplappp",
                "ABplpaaa", "ABplpaap", "ABplpapa", "ABplpapp", "ABplppaa", "ABplppap", "ABplpppa", "ABplpppp",
                "ABpraaaa", "ABpraaap", "ABpraapa", "ABpraapp", "ABprapaa", "ABprapap", "ABprappa", "ABprappp",
                "ABprpaaa", "ABprpaap", "ABprpapa", "ABprpapp", "ABprppaa", "ABprppap", "ABprpppa", "ABprpppp",
                "MSaaaa", "MSaaap", "MSaapa", "MSaapp", "MSapaa", "MSapap", "MSappa", "MSappp",
                "MSpaaa", "MSpaap", "MSpapa", "MSpapp", "MSppaa", "MSppap", "MSpppa", "MSpppp",
                "Caaa", "Caap", "Capa", "Capp", "Cpaa", "Cpap", "Cppa", "Cppp")
Stage_AB128 <- c("ABalaaaal","ABalaaaar","ABalaaapa","ABalaaapp","ABalaapaa","ABalaapap","ABalaappa","ABalaappp",
                 "ABalapaaa","ABalapaap","ABalapapa","ABalapapp","ABalappaa","ABalappap","ABalapppa","ABalapppp",
                 "ABalpaaaa","ABalpaaap","ABalpaapa","ABalpaapp","ABalpapaa","ABalpapap","ABalpappa","ABalpappp",
                 "ABalppaaa","ABalppaap","ABalppapa","ABalppapp","ABalpppaa","ABalpppap","ABalppppa","ABalppppp",
                 "ABaraaaaa","ABaraaaap","ABaraaapa","ABaraaapp","ABaraapaa","ABaraapap","ABaraappa","ABaraappp",
                 "ABarapaaa","ABarapaap","ABarapapa","ABarapapp","ABarappaa","ABarappap","ABarapppa","ABarapppp",
                 "ABarpaaaa","ABarpaaap","ABarpaapa","ABarpaapp","ABarpapaa","ABarpapap","ABarpappa","ABarpappp",
                 "ABarppaaa","ABarppaap","ABarppapa","ABarppapp","ABarpppaa","ABarpppap","ABarppppa","ABarppppp",
                 "ABplaaaaa","ABplaaaap","ABplaaapa","ABplaaapp","ABplaapaa","ABplaapap","ABplaappa","ABplaappp",
                 "ABplapaaa","ABplapaap","ABplapapa","ABplapapp","ABplappaa","ABplappap","ABplapppa","ABplapppp",
                 "ABplpaaaa","ABplpaaap","ABplpaapa","ABplpaapp","ABplpapaa","ABplpapap","ABplpappa","ABplpappp",
                 "ABplppaaa","ABplppaap","ABplppapa","ABplppapp","ABplpppaa","ABplpppap","ABplppppa","ABplppppp",
                 "ABpraaaaa","ABpraaaap","ABpraaapa","ABpraaapp","ABpraapaa","ABpraapap","ABpraappa","ABpraappp",
                 "ABprapaaa","ABprapaap","ABprapapa","ABprapapp","ABprappaa","ABprappap","ABprapppa","ABprapppp",
                 "ABprpaaaa","ABprpaaap","ABprpaapa","ABprpaapp","ABprpapaa","ABprpapap","ABprpappa","ABprpappp",
                 "ABprppaaa","ABprppaap","ABprppapa","ABprppapp","ABprpppaa","ABprpppap","ABprppppa","ABprppppp",
                 "MSaaaaa","MSaaaap","MSaaapa","MSaapaa","MSaapap","MSapaaa","MSapapa","MSapapp",
                 "MSpaaaa","MSpaaap","MSpaapa","MSpapaa","MSpapap","MSppaaa","MSppapa","MSppapp",
                 "Eala","Ealp","Eara","Earp","Epla","Eplp","Epra","Eprp",
                 "Caaaa","Caaap","Caapp","Capaa","Capap","Cappa","Cappp",
                 "Cpaaa","Cpaap","Cpapa","Cpapp","Cppaa","Cppap","Cpppa","Cpppp",
                 "Daa","Dap","Dpa","Dpp")
Stage_AB256 <- c("ABalaaappr","ABalaapppa","ABalaapppp","ABalapaapp","ABalapappa","ABalappapp",
                 "ABalapppap","ABalappppa","ABalappppp","ABalpaapaa","ABalpaapap","ABalpaappa",
                 "ABalpaappp","ABalpapaaa","ABalpapaap","ABalpapapa","ABalpapapp","ABalpappaa",
                 "ABalpappap","ABalpapppp","ABalppapaa","ABalppapap","ABalppappa","ABalppappp",
                 "ABalpppapa","ABalpppapp","ABalppppaa","ABalppppap","ABalpppppa","ABalpppppp",
                 "ABaraaapaa","ABaraaapap","ABaraaappa","ABaraaappp","ABaraapaaa","ABaraapaap",
                 "ABaraapapa","ABaraapapp","ABaraappaa","ABaraappap","ABaraapppa","ABaraapppp",
                 "ABarapaapp","ABarapapap","ABarapappp","ABarappaap","ABarappapa","ABarappapp",
                 "ABarapppaa","ABarapppap","ABarappppa","ABarappppp","ABarpapaap",
                 "ABplaaaaap","ABplaapaap","ABplaapapa","ABplaapapp","ABplapaaaa","ABplapaaap",
                 "ABplapappp","ABplapppap","ABplpaaaaa","ABplpaaaap","ABplpaaapa","ABplpaaapp",
                 "ABplpaapaa","ABplpaapap","ABplpaappa","ABplpaappp","ABplpapaaa","ABplpapaap",
                 "ABplpapapa","ABplpapapp","ABplpappaa","ABplpapppa","ABplppaaaa","ABplppaapa",
                 "ABplppaapp","ABplppapaa","ABplppapap","ABplppappa","ABplppappp","ABplpppaaa",
                 "ABplpppaap","ABplpppapa","ABplppppaa","ABplppppap","ABplpppppa","ABplpppppp",
                 "ABpraaappp","ABpraapaap","ABpraapapp","ABprapaaap","ABprpaaaaa","ABprpaaaap",
                 "ABprpaaapa","ABprpaaapp","ABprpaapaa","ABprpaapap","ABprpaappa","ABprpaappp",
                 "ABprpapaaa","ABprpapaap","ABprpapapa","ABprpapapp","ABprpappaa","ABprpappap",
                 "ABprpapppa","ABprpapppp","ABprppaaaa","ABprppaapp","ABprppapaa","ABprppapap",
                 "ABprppappa","ABprppappp","ABprpppaaa","ABprpppaap","ABprpppapa","ABprppppaa",
                 "ABprppppap","ABprpppppa","ABprpppppp",
                 "Caapa","Capaaa","Capaap","Capapa","Capapp","Cappaa","Cappap","Capppa","Capppp",
                 "Cppaaa","Cppaap","Cppapa","Cppapp","Cpppaa","Cpppap","Cppppa","Cppppp",
                 "Daaa","Daap","Dapa","Dapp","Dpaa","Dpap","Dppa","Dppp",
                 "MSaaapp","MSaappa","MSaappp","MSapaap","MSappaa","MSappap","MSapppa","MSapppp",
                 "MSpappa","MSpappp","MSppaap","MSpppaa","MSpppap","MSppppa","MSppppp")

stage_assign <- function(cell) {
  case_when(
    cell %in% Stage_AB4   ~ "Stage_AB4",
    cell %in% Stage_AB8   ~ "Stage_AB8",
    cell %in% Stage_AB16  ~ "Stage_AB16",
    cell %in% Stage_AB32  ~ "Stage_AB32",
    cell %in% Stage_AB64  ~ "Stage_AB64",
    cell %in% Stage_AB128 ~ "Stage_AB128",
    cell %in% Stage_AB256 ~ "Stage_AB256",
    TRUE ~ NA_character_
  )
}

# ================================================================
# 2. 数据集信息
# ================================================================
make_meta <- function(prefix, files, times, scale) {
  tibble(dataset = paste0(prefix, seq_along(files)),
         file    = files,
         time    = times,
         scale   = scale)
}

ce_meta <- make_meta("ce",
                     c("CDFile/ce/CD191108plc1p1.csv",
                       "CDFile/ce/CD200109plc1p1.csv",
                       "CDFile/ce/CD200113plc1p3.csv",
                       "CDFile/ce/CD200113plc1p2.csv",
                       "CDFile/ce/CD200322plc1p2.csv",
                       "CDFile/ce/CD200323plc1p1.csv",
                       "CDFile/ce/CD200326plc1p3.csv",
                       "CDFile/ce/CD200326plc1p4.csv"),
                     c(205, 205, 195, 205, 195, 185, 220, 195),
                     4.67)

cb_meta <- make_meta("cb",
                     c("CDFile/she1/CD240731cbhis72p1.csv",
                       "CDFile/she1/CD240731cbhis72p2.csv",
                       "CDFile/she1/CD240731cbhis72p3.csv",
                       "CDFile/she1/CD241202cbhis72p1.csv",
                       "CDFile/she1/CD241202cbhis72p2.csv",
                       "CDFile/she1/CD241202cbhis72p4.csv"),
                     c(165, 175, 170, 160, 165, 180),
                     4.78)

cn_meta <- make_meta("cn",
                     c("CDFile/cn/CD241202cnhis72p1.csv",
                       "CDFile/cn/CD240712cnhis72p1.csv",
                       "CDFile/cn/CD240712cnhis72p2.csv",
                       "CDFile/cn/CD240712cnhis72p3.csv",
                       "CDFile/cn/CD241207cnhis72p1.csv",
                       "CDFile/cn/CD241207cnhis72p3.csv",
                       "CDFile/cn/CD241202cnhis72p2.csv"),
                     c(235, 235, 235, 235, 230, 230, 200),
                     4.78)

dataset_meta <- bind_rows(ce_meta, cb_meta, cn_meta)

# ================================================================
# 3. 第一条链：final_result & cell division angle
#    ——完全保持原脚本的处理方式（无 stage）
# ================================================================
process_dis_base <- function(data, time_limit, scale_factor) {
  data %>%
    select(2, 3, 9, 10, 11) %>%
    mutate(z = z * scale_factor) %>%
    group_by(cell) %>%
    as_tibble() %>%
    filter(time <= time_limit) %>%
    mutate(parent = case_when(
      cell %in% c("EMS", "P2") ~ "P1",
      cell %in% c("MS", "E")   ~ "EMS",
      cell %in% c("P3", "C")   ~ "P2",
      cell %in% c("P4", "D")   ~ "P3",
      cell %in% c("Z2", "Z3")  ~ "P4",
      TRUE ~ substr(cell, 1, nchar(cell) - 1)
    )) %>%
    rename(child = cell) %>%
    filter(!(child %in% c("AB", "P1"))) %>%
    drop_na()
}

all_list_base <- dataset_meta %>%
  mutate(data = pmap(list(file, time, scale),
                     ~read.csv(..1, header = TRUE) %>%
                       process_dis_base(time_limit = ..2, scale_factor = ..3))) %>%
  select(dataset, data) %>%
  deframe()

get_max_distance_info <- function(df, df_name = NA) {
  child_repr <- df %>%
    group_by(child) %>%
    slice_max(order_by = time, n = 2, with_ties = FALSE) %>%
    ungroup()
  
  coords <- child_repr %>% select(x, y, z)
  dist_mat <- as.matrix(dist(coords, method = "euclidean"))
  diag(dist_mat) <- NA
  
  upper_idx <- which(upper.tri(dist_mat), arr.ind = TRUE)
  distance_df <- tibble(
    i = upper_idx[, 1],
    j = upper_idx[, 2],
    distance = dist_mat[upper.tri(dist_mat)]
  ) %>% arrange(desc(distance))
  
  top10        <- head(distance_df, 10)
  median_index <- floor((nrow(top10) + 1) / 2)
  median_pair  <- top10[median_index, ]
  
  cell1_info <- child_repr[median_pair$i, ]
  cell2_info <- child_repr[median_pair$j, ]
  ba_vec     <- c(cell2_info$x - cell1_info$x,
                  cell2_info$y - cell1_info$y,
                  cell2_info$z - cell1_info$z)
  
  fetch_vec <- function(df, cell_a, cell_b) {
    tmp <- df %>%
      filter(child %in% c(cell_a, cell_b)) %>%
      group_by(child) %>%
      slice_max(time, n = 1, with_ties = FALSE) %>%
      ungroup()
    if (nrow(tmp) == 2) {
      c(
        tmp$x[tmp$child == cell_a] - tmp$x[tmp$child == cell_b],
        tmp$y[tmp$child == cell_a] - tmp$y[tmp$child == cell_b],
        tmp$z[tmp$child == cell_a] - tmp$z[tmp$child == cell_b]
      )
    } else {
      NA
    }
  }
  
  EMS_ABp <- fetch_vec(df, "EMS", "ABp")
  P2_ABa  <- fetch_vec(df, "P2", "ABa")
  
  dot_ba_P2_ba <- if (is.numeric(P2_ABa) && length(P2_ABa) == 3 &&
                      all(!is.na(P2_ABa)) && all(!is.na(ba_vec))) {
    sum(P2_ABa[1:2] * ba_vec[1:2])
  } else {
    NA_real_
  }
  
  if (!is.na(dot_ba_P2_ba) && dot_ba_P2_ba < 0) {
    ba_vec <- -ba_vec
  }
  
  tibble(
    dataset      = df_name,
    max_distance = median_pair$distance,
    cell1        = as.character(cell1_info$child),
    cell1_x      = cell1_info$x,
    cell1_y      = cell1_info$y,
    cell1_z      = cell1_info$z,
    cell2        = as.character(cell2_info$child),
    cell2_x      = cell2_info$x,
    cell2_y      = cell2_info$y,
    cell2_z      = cell2_info$z,
    ba           = list(ba_vec),
    EMS_ABp      = list(EMS_ABp),
    P2_ABa       = list(P2_ABa),
    dot_ba_P2_ba = dot_ba_P2_ba
  )
}

final_result <- imap_dfr(all_list_base, get_max_distance_info)

calculate_angle <- function(u, v) {
  cos_angle <- sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
  cos_angle <- pmin(pmax(cos_angle, -1), 1)
  acos(cos_angle) * 180 / pi
}

process_ce <- function(data, dataset_label, final_tbl) {
  ce_new <- data %>%
    group_by(child) %>%
    filter(min_rank(time) == 2) %>%
    ungroup() %>%
    group_by(parent) %>%
    reframe(
      child1 = child[grepl("[prv34]$", child)],
      child2 = child[grepl("[adlCD]$", child)],
      dc = list(c(
        x = x[grepl("[prv34]$", child)] - x[grepl("[adlCD]$", child)],
        y = y[grepl("[prv34]$", child)] - y[grepl("[adlCD]$", child)],
        z = z[grepl("[prv34]$", child)] - z[grepl("[adlCD]$", child)]
      ))
    ) %>%
    ungroup()
  
  target  <- final_tbl %>% filter(dataset == dataset_label)
  ba_vec  <- target$ba[[1]]
  EMS_vec <- target$EMS_ABp[[1]]
  
  res <- sapply(ce_new$dc, function(dc_vec) {
    angle_val <- calculate_angle(ba_vec, dc_vec)
    EMS_val   <- sum(EMS_vec[1:2] * dc_vec[1:2])
    if (!is.na(EMS_val) && EMS_val < 0) angle_val <- -angle_val
    c(angle = angle_val, EMS = EMS_val)
  })
  
  ce_new$angle <- res["angle", ]
  ce_new$EMS   <- res["EMS", ]
  ce_new
}

division_list <- imap(all_list_base, ~process_ce(.x, .y, final_result))

extract_df <- function(df, label) {
  df %>% select(1:3, angle) %>% rename(!!paste0("angle_", label) := angle)
}

division_df_list <- imap(division_list, extract_df)
common_cols <- names(division_df_list[[1]])[1:3]
merged_df   <- reduce(division_df_list, inner_join, by = common_cols)

output_dir <- file.path(root_dir, "draft_code/new/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(merged_df, file = file.path(output_dir, "all_divisionAngle.csv"), row.names = FALSE)

# ================================================================
# 4. 第二条链：角度范围 & 斜率
#    ——与原脚本一致，含 stage 并对 NA 行执行 na.omit()
# ================================================================
process_dis_stage <- function(data, time_limit, scale_factor) {
  data %>%
    select(2, 3, 9, 10, 11) %>%
    mutate(z = z * scale_factor) %>%
    group_by(cell) %>%
    as_tibble() %>%
    filter(time <= time_limit) %>%
    mutate(parent = case_when(
      cell %in% c("EMS", "P2") ~ "P1",
      cell %in% c("MS", "E")   ~ "EMS",
      cell %in% c("P3", "C")   ~ "P2",
      cell %in% c("P4", "D")   ~ "P3",
      cell %in% c("Z2", "Z3")  ~ "P4",
      TRUE ~ substr(cell, 1, nchar(cell) - 1)
    )) %>%
    rename(child = cell) %>%
    filter(!(child %in% c("AB", "P1"))) %>%
    mutate(
      stage  = stage_assign(child),
      xorder = substr(parent, 1, 1),
      xorder = factor(ifelse(xorder %in% c("P", "Z"), "D", xorder),
                      levels = c("A", "M", "E", "C", "D"))
    ) %>%
    drop_na()
}

dataset_stage_list <- dataset_meta %>%
  mutate(data = pmap(list(file, time, scale),
                     ~read.csv(..1, header = TRUE) %>%
                       process_dis_stage(time_limit = ..2, scale_factor = ..3))) %>%
  select(dataset, data) %>%
  deframe()

options(scipen = 999)

angle_between <- function(v1, v2) {
  dot_product <- sum(v1 * v2)
  norm1 <- sqrt(sum(v1^2))
  norm2 <- sqrt(sum(v2^2))
  cos_angle <- dot_product / (norm1 * norm2)
  cos_angle <- min(max(cos_angle, -1), 1)
  acos(cos_angle)
}

result_list <- map2(dataset_stage_list, names(dataset_stage_list), function(df, name) {
  df %>%
    mutate(dataset = name) %>%
    group_by(child) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(
      base_x    = x[2],
      base_y    = y[2],
      base_z    = z[2],
      base_time = time[2]
    ) %>%
    filter(time > base_time) %>%
    rowwise() %>%
    mutate(
      angle_rad = angle_between(c(x, y, z), c(base_x, base_y, base_z)),
      angle_deg = angle_rad * 180 / pi
    ) %>%
    ungroup()
})

result <- bind_rows(result_list)

angle_summary_range <- result %>%
  group_by(dataset, child, stage) %>%
  summarise(
    min_angle   = min(angle_deg, na.rm = TRUE),
    max_angle   = max(angle_deg, na.rm = TRUE),
    angle_range = max_angle - min_angle,
    .groups = "drop"
  ) %>%
  select(dataset, child, stage, angle_range)

wide_range <- angle_summary_range %>%
  pivot_wider(
    names_from  = dataset,
    values_from = angle_range,
    names_prefix = "angle_range_"
  )

#write.csv(wide_range, file = file.path(output_dir, "con_angle_range.csv"), row.names = FALSE)

angle_summary_slope <- result %>%
  group_by(dataset, child, stage) %>%
  summarise(
    angle_slope = if (n() > 1) coef(lm(angle_deg ~ time))[2] else NA_real_,
    .groups = "drop"
  )

wide_slope <- angle_summary_slope %>%
  pivot_wider(
    names_from  = dataset,
    values_from = angle_slope,
    names_prefix = "angle_slope_"
  )

#write.csv(wide_slope, file = file.path(output_dir, "con_angle_slope.csv"), row.names = FALSE)





setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
library(tidyr)
library(stringr)
library(ggbeeswarm)
library(beeswarm)

data <- read.csv("~/Desktop/cbcn/draft_code/new/all_divisionAngle.csv")
colnames(data)
colnames(data)[1] <- "cell"
colnames(data)

data <- data %>%
  mutate(
    mean_ce = rowMeans(select(., starts_with("angle_ce")), na.rm = TRUE),
    mean_cb = rowMeans(select(., starts_with("angle_cb")), na.rm = TRUE),
    mean_cn = rowMeans(select(., starts_with("angle_cn")), na.rm = TRUE)
  ) %>%
  select(cell, mean_ce, mean_cb, mean_cn) 

colnames(data)


# 假设你的数据存储在变量 data 中
# 利用 tidyr::pivot_longer 将数据从宽格式转换为长格式
data_long <- pivot_longer(
  data,
  cols = c("mean_ce", "mean_cb", "mean_cn"),
  names_to = "group",
  values_to = "value"
)
# 修改分组名称（去掉前缀 "mean_"），使之和颜色映射对应
data_long$group <- gsub("mean_", "", data_long$group)

# 设置分组顺序为 ce, cb, cn
data_long$group <- factor(data_long$group, levels = c("cn", "cb", "ce"))

# 绘制蜂巢图，并设置四边框和自定义颜色
# p_4f <- ggplot(data_long, aes(x = value, y = group, color = group)) +
#   geom_beeswarm(size = 0.6, shape = 16) +
#   labs(
#     title = "distribution of cell division angle",
#     x = "value",
#     y = ""
#   ) +
#   theme_minimal() +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     plot.background = element_blank(),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#     axis.line = element_blank(),
#     axis.ticks = element_line(color = "black"),            # 显示刻度线
#     axis.text = element_text(size = 12),
#     legend.position = "none"
#   ) +
#   scale_color_manual(values = c("ce" = "#CD534C", "cb" = "#0073C2", "cn" = "#EFC000"))
# 
# print(p_4f)
# # 绘制热图
# ggsave(
#   filename = "divangle_distribution.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 6,                            # 图形宽度（单位默认为英寸）
#   height = 3,                            # 图形高度
#   dpi = 300                              # 分辨率
# )
# 如果还没有安装 ggridges，请先安装
#install.packages("ggridges")
library(ggridges)

p_4f <- ggplot(data_long, aes(x = value, y = group, fill = group)) +
  geom_density_ridges(alpha = 0.7, scale = 1.0, 
                      color = "black", size = 0.3) +
  labs(
    title = "Distribution of Cell Division Angle (Ridge Plot)",
    x = "Value",
    y = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("ce" = "#CD534C", "cb" = "#0073C2", "cn" = "#EFC000"))

print(p_4f)
# # # 绘制热图
# ggsave(
#   filename = "divangle_distribution.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 6,                            # 图形宽度（单位默认为英寸）
#   height = 3,                            # 图形高度
#   dpi = 300                              # 分辨率
# )


# Compute the density for each group and extract the peaks
peak_values <- data_long %>%
  group_by(group) %>%
  summarise(
    density_obj = list(density(value, na.rm = TRUE))
  ) %>%
  mutate(
    # For each group's density object, extract the x value where density is maximum.
    peak_value = sapply(density_obj, function(d) d$x[which.max(d$y)]),
    # Optionally, extract the peak density (i.e., maximum y value).
    peak_density = sapply(density_obj, function(d) max(d$y))
  )

# Print the summary to view the peak positions and corresponding density values for ce, cb, and cn
print(peak_values)


#处理模板数据
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_angle_range.csv", header = TRUE)
colnames(data)
colnames(data)[1] <- "cell"

data <- data %>%
  mutate(
    mean_ce = rowMeans(select(., starts_with("angle_range_ce")), na.rm = TRUE),
    mean_cb = rowMeans(select(., starts_with("angle_range_cb")), na.rm = TRUE),
    mean_cn = rowMeans(select(., starts_with("angle_range_cn")), na.rm = TRUE)
  ) %>%
  select(cell, mean_ce, mean_cb, mean_cn) 
colnames(data)


# 假设你的数据存储在变量 data 中
# 利用 tidyr::pivot_longer 将数据从宽格式转换为长格式
data_long <- pivot_longer(
  data,
  cols = c("mean_ce", "mean_cb", "mean_cn"),
  names_to = "group",
  values_to = "value"
)
# 修改分组名称（去掉前缀 "mean_"），使之和颜色映射对应
data_long$group <- gsub("mean_", "", data_long$group)

# 设置分组顺序为 ce, cb, cn
data_long$group <- factor(data_long$group, levels = c("cn", "cb", "ce"))

library(ggridges)
# 绘制蜂巢图，并设置四边框和自定义颜色
ggplot(data_long, aes(x = value, y = group, fill = group)) +
  geom_density_ridges(alpha = 0.7, scale = 1.0, 
                      color = "black", size = 0.3) +
  labs(
    title = "Ange Range",
    x = "value",
    y = "group"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),            # 显示刻度线
    axis.text = element_text(size = 12),
  ) +
  scale_fill_manual(values = c("ce" = "#CD534C", "cb" = "#0073C2", "cn" = "#EFC000"))

# # 绘制热图
# ggsave(
#   filename = "AngleRange_distribution.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 6,                            # 图形宽度（单位默认为英寸）
#   height = 3,                            # 图形高度
#   dpi = 300                              # 分辨率
# )




#处理模板数据
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_angle_slope.csv", header = TRUE)
colnames(data)
colnames(data)[1] <- "cell"

data <- data %>%
  mutate(
    mean_ce = rowMeans(select(., starts_with("angle_slope_ce")), na.rm = TRUE),
    mean_cb = rowMeans(select(., starts_with("angle_slope_cb")), na.rm = TRUE),
    mean_cn = rowMeans(select(., starts_with("angle_slope_cn")), na.rm = TRUE)
  ) %>%
  select(cell, mean_ce, mean_cb, mean_cn) 
colnames(data)


# 假设你的数据存储在变量 data 中
# 利用 tidyr::pivot_longer 将数据从宽格式转换为长格式
data_long <- pivot_longer(
  data,
  cols = c("mean_ce", "mean_cb", "mean_cn"),
  names_to = "group",
  values_to = "value"
)
# 修改分组名称（去掉前缀 "mean_"），使之和颜色映射对应
data_long$group <- gsub("mean_", "", data_long$group)

# 设置分组顺序为 ce, cb, cn
data_long$group <- factor(data_long$group, levels = c("cn", "cb", "ce"))

# 绘制蜂巢图，并设置四边框和自定义颜色
ggplot(data_long, aes(x = value, y = group, fill = group)) +
  geom_density_ridges(alpha = 0.7, scale = 1.0, 
                      color = "black", size = 0.3) +
  labs(
    title = "Angle Slope",
    x = "value",
    y = "group"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),            # 显示刻度线
    axis.text = element_text(size = 12),
  ) +
  scale_fill_manual(values = c("ce" = "#CD534C", "cb" = "#0073C2", "cn" = "#EFC000"))

# # 绘制热图
# ggsave(
#   filename = "angleSlope_distribution.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 6,                            # 图形宽度（单位默认为英寸）
#   height = 3,                            # 图形高度
#   dpi = 300                              # 分辨率
# )




setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
data <- read.csv("~/Desktop/cbcn/draft_code/new/all_divisionAngle.csv")
colnames(data)
colnames(data)[1] <- "cell"
colnames(data)


stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(data, stage, by = "cell")
colnames(data)

data <- data %>%
  mutate(
    mean_ce = rowMeans(select(., starts_with("angle_ce")), na.rm = TRUE),
    mean_cb = rowMeans(select(., starts_with("angle_cb")), na.rm = TRUE),
    mean_cn = rowMeans(select(., starts_with("angle_cn")), na.rm = TRUE)
  ) %>%
  select(cell, mean_ce, mean_cb, mean_cn,stage) 

colnames(data)

colnames(data)


data_long <- data %>%
  pivot_longer(
    cols = c("mean_ce", "mean_cb", "mean_cn"),
    names_to = "measure",
    values_to = "value"
  )
# 2. 修改 factor 标签，使 y 轴显示 ce, cb, cn
data_long$measure <- factor(data_long$measure,
                            levels = c("mean_cn", "mean_cb", "mean_ce"),
                            labels = c("cn", "cb", "ce"))

order_stage  <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                  "Stage_AB64", "Stage_AB128", "Stage_AB256")
order_letter <- c("A", "M", "E", "C", "D", "P")

# 2. 生成新变量：提取 cell 的首字母和计算 cell 名称长度
data_long <- data_long %>%
  mutate(
    # 强制 stage 按照指定顺序排序
    stage = factor(stage, levels = order_stage),
    # 提取 cell 的首字母并转为因子
    first_letter = factor(substr(cell, 1, 1), levels = order_letter),
    # 计算 cell 名称的长度
    cell_length = nchar(cell)
  ) %>%
  # 3. 按照 stage, 首字母, 及 cell 名称长度排序
  arrange(stage, first_letter, cell_length)

# 4. 更新 cell 为因子并按照排序好的顺序
data_long$cell <- factor(data_long$cell, levels = unique(data_long$cell))
colnames(data_long)
# 绘制热图
ggplot(data_long, aes(x = cell, y = measure, fill = value)) +
  geom_tile(color = "white") +
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  scale_fill_gradientn(
    #colors = c("#fbebf0", "#f7d7de", "#edb2bf", "#e38b9f", "#c95f7a"),  # -10 对应蓝色，7 对应红色，中间过渡为白色
    #colors = c("#c5ebde", "#b1d2c7", "#7d9d90", "#496968", "#39515d"), 
    colors = c("#edf6fb", "#c1d6e9", "#95b6d7", "#628db0", "#3674b1"), 
    limits = c(0, 10),
    oob = scales::squish              # 将超出 [-10,7] 的值压缩到边界
  ) +
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x  = element_blank(),  # 移除 x 轴所有 cell 的标签
    axis.text.y = element_text(size = 10),              # y 轴标签字体大小
    legend.position = "right",
    # 为每个 facet 添加黑色边框
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # 以下项移除所有背景和网格线：
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )


# # 绘制热图
# ggsave(
#   filename = "divangle_value.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 10,                            # 图形宽度（单位默认为英寸）
#   height = 2.5,                            # 图形高度
#   dpi = 300                              # 分辨率
# )





#处理模板数据
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_angle_range.csv", header = TRUE)
colnames(data)
colnames(data)[1] <- "cell"

data <- data %>%
  mutate(
    mean_ce = rowMeans(select(., starts_with("angle_range_ce")), na.rm = TRUE),
    mean_cb = rowMeans(select(., starts_with("angle_range_cb")), na.rm = TRUE),
    mean_cn = rowMeans(select(., starts_with("angle_range_cn")), na.rm = TRUE)
  ) %>%
  select(cell, mean_ce, mean_cb, mean_cn) 
colnames(data)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(data, stage, by = "cell")
colnames(data)

data_long <- data %>%
  pivot_longer(
    cols = c("mean_ce", "mean_cb", "mean_cn"),
    names_to = "measure",
    values_to = "value"
  )
# 2. 修改 factor 标签，使 y 轴显示 ce, cb, cn
data_long$measure <- factor(data_long$measure,
                            levels = c("mean_cn", "mean_cb", "mean_ce"),
                            labels = c("cn", "cb", "ce"))

order_stage  <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                  "Stage_AB64", "Stage_AB128", "Stage_AB256")
order_letter <- c("A", "M", "E", "C", "D", "P")

# 2. 生成新变量：提取 cell 的首字母和计算 cell 名称长度
data_long <- data_long %>%
  mutate(
    # 强制 stage 按照指定顺序排序
    stage = factor(stage, levels = order_stage),
    # 提取 cell 的首字母并转为因子
    first_letter = factor(substr(cell, 1, 1), levels = order_letter),
    # 计算 cell 名称的长度
    cell_length = nchar(cell)
  ) %>%
  # 3. 按照 stage, 首字母, 及 cell 名称长度排序
  arrange(stage, first_letter, cell_length)

# 4. 更新 cell 为因子并按照排序好的顺序
data_long$cell <- factor(data_long$cell, levels = unique(data_long$cell))
colnames(data_long)
# 绘制热图
ggplot(data_long, aes(x = cell, y = measure, fill = value)) +
  geom_tile(color = "white") +
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  scale_fill_gradientn(
    #colors = c("#fbebf0", "#f7d7de", "#edb2bf", "#e38b9f", "#c95f7a"),  # -10 对应蓝色，7 对应红色，中间过渡为白色
    #colors = c("#c5ebde", "#b1d2c7", "#7d9d90", "#496968", "#39515d"), 
    colors = c("#edf6fb", "#c1d6e9", "#95b6d7", "#628db0", "#3674b1"), 
    limits = c(0, 10),
    oob = scales::squish              # 将超出 [-10,7] 的值压缩到边界
  ) +
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x  = element_blank(),  # 移除 x 轴所有 cell 的标签
    axis.text.y = element_text(size = 10),              # y 轴标签字体大小
    legend.position = "right",
    # 为每个 facet 添加黑色边框
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # 以下项移除所有背景和网格线：
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )


# # 绘制热图
# ggsave(
#   filename = "angleRange_value.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 10,                            # 图形宽度（单位默认为英寸）
#   height = 2.5,                            # 图形高度
#   dpi = 300                              # 分辨率
# )




#处理模板数据
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_angle_slope.csv", header = TRUE)
colnames(data)
colnames(data)[1] <- "cell"

data <- data %>%
  mutate(
    mean_ce = rowMeans(select(., starts_with("angle_slope_ce")), na.rm = TRUE),
    mean_cb = rowMeans(select(., starts_with("angle_slope_cb")), na.rm = TRUE),
    mean_cn = rowMeans(select(., starts_with("angle_slope_cn")), na.rm = TRUE)
  ) %>%
  select(cell, mean_ce, mean_cb, mean_cn) 
colnames(data)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(data, stage, by = "cell")
colnames(data)

data_long <- data %>%
  pivot_longer(
    cols = c("mean_ce", "mean_cb", "mean_cn"),
    names_to = "measure",
    values_to = "value"
  )
# 2. 修改 factor 标签，使 y 轴显示 ce, cb, cn
data_long$measure <- factor(data_long$measure,
                            levels = c("mean_cn", "mean_cb", "mean_ce"),
                            labels = c("cn", "cb", "ce"))

order_stage  <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                  "Stage_AB64", "Stage_AB128", "Stage_AB256")
order_letter <- c("A", "M", "E", "C", "D", "P")

# 2. 生成新变量：提取 cell 的首字母和计算 cell 名称长度
data_long <- data_long %>%
  mutate(
    # 强制 stage 按照指定顺序排序
    stage = factor(stage, levels = order_stage),
    # 提取 cell 的首字母并转为因子
    first_letter = factor(substr(cell, 1, 1), levels = order_letter),
    # 计算 cell 名称的长度
    cell_length = nchar(cell)
  ) %>%
  # 3. 按照 stage, 首字母, 及 cell 名称长度排序
  arrange(stage, first_letter, cell_length)

# 4. 更新 cell 为因子并按照排序好的顺序
data_long$cell <- factor(data_long$cell, levels = unique(data_long$cell))
colnames(data_long)
# 绘制热图
ggplot(data_long, aes(x = cell, y = measure, fill = value)) +
  geom_tile(color = "white") +
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  scale_fill_gradientn(
    #colors = c("#fbebf0", "#f7d7de", "#edb2bf", "#e38b9f", "#c95f7a"),  # -10 对应蓝色，7 对应红色，中间过渡为白色
    #colors = c("#c5ebde", "#b1d2c7", "#7d9d90", "#496968", "#39515d"), 
    colors = c("#edf6fb", "#c1d6e9", "#95b6d7", "#628db0", "#3674b1"), 
    limits = c(0, 0.8),
    oob = scales::squish              # 将超出 [-10,7] 的值压缩到边界
  ) +
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x  = element_blank(),  # 移除 x 轴所有 cell 的标签
    axis.text.y = element_text(size = 10),              # y 轴标签字体大小
    legend.position = "right",
    # 为每个 facet 添加黑色边框
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # 以下项移除所有背景和网格线：
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )


# # 绘制热图
# ggsave(
#   filename = "angleSlope_value.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 10,                            # 图形宽度（单位默认为英寸）
#   height = 2.5,                            # 图形高度
#   dpi = 300                              # 分辨率
# )

