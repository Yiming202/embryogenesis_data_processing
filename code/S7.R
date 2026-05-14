setwd("~/Desktop/cbcn/draft_code/new/")
getwd()
library(ggplot2);library(dplyr);library(tidyverse);library(ggsci);library(ggpubr);library(reshape2);library(RColorBrewer)

#处理模板数据
ce1 <- read.csv("~/Desktop/cbcn/CDFile/ce/CD191108plc1p1.csv", header = TRUE)
ce2 <- read.csv("~/Desktop/cbcn/CDFile/ce/CD200109plc1p1.csv", header = TRUE)
ce3 <- read.csv("~/Desktop/cbcn/CDFile/ce/CD200113plc1p3.csv", header = TRUE)
ce4 <- read.csv("~/Desktop/cbcn/CDFile/ce/CD200113plc1p2.csv", header = TRUE)
ce5 <- read.csv("~/Desktop/cbcn/CDFile/ce/CD200322plc1p2.csv", header = TRUE)
ce6 <- read.csv("~/Desktop/cbcn/CDFile/ce/CD200323plc1p1.csv", header = TRUE)
ce7 <- read.csv("~/Desktop/cbcn/CDFile/ce/CD200326plc1p3.csv", header = TRUE)
ce8 <- read.csv("~/Desktop/cbcn/CDFile/ce/CD200326plc1p4.csv", header = TRUE)

cb1 <- read.csv("~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p1.csv", header = TRUE)
cb2 <- read.csv("~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p2.csv", header = TRUE)
cb3 <- read.csv("~/Desktop/cbcn/CDFile/she1/CD240731cbhis72p3.csv", header = TRUE)
cb4 <- read.csv("~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p1.csv", header = TRUE)
cb5 <- read.csv("~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p2.csv", header = TRUE)
cb6 <- read.csv("~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p4.csv", header = TRUE)

cn1 <- read.csv("~/Desktop/cbcn/CDFile/cn/CD241202cnhis72p1.csv", header = TRUE)
cn2 <- read.csv("~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p1.csv", header = TRUE)
cn3 <- read.csv("~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p2.csv", header = TRUE)
cn4 <- read.csv("~/Desktop/cbcn/CDFile/cn/CD240712cnhis72p3.csv", header = TRUE)
cn5 <- read.csv("~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p1.csv", header = TRUE)
cn6 <- read.csv("~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p3.csv", header = TRUE)
cn7 <- read.csv("~/Desktop/cbcn/CDFile/cn/CD241202cnhis72p2.csv", header = TRUE)

# ------------------------------------------------------------
# 3. 主循环：逐个处理每个胚胎文件
# ------------------------------------------------------------
# 定义处理函数
#筛选editing tp, correct z
process_dis <- function(data, time_limit, intercept, slope) {
  data %>%
    select(2,3,9,10,11) %>%
    mutate(z = z * 4.67) %>%
    group_by(cell) %>%
    as.tibble() %>%
    filter(time <= time_limit) %>%
    select(1,2,4,5,3) %>%
    mutate(time = (time - intercept) / slope)
}
ce1 <- process_dis(ce1, 205,7.1,1)
ce2 <- process_dis(ce2, 205,-9.57,1.05)
ce3 <- process_dis(ce3, 195,-3.95,1)
ce4 <- process_dis(ce4, 205,7.04,1.04)
ce5 <- process_dis(ce5, 195,1.79,0.95)
ce6 <- process_dis(ce6, 185,-7.77,0.95)
ce7 <- process_dis(ce7, 220,4.46,1.04)
ce8 <- process_dis(ce8, 195,0.89,0.98)

process_dis <- function(data, time_limit, intercept, slope) {
  data %>%
    select(2,3,9,10,11) %>%
    mutate(z = z * 4.78) %>%
    group_by(cell) %>%
    as.tibble() %>%
    filter(time <= time_limit) %>%
    select(1,2,4,5,3) %>%
    mutate(time = (time - intercept) / slope)
}
# 处理每个数据框
cb1 <- process_dis(cb1, 165,6.62,0.91)
cb2 <- process_dis(cb2, 175,8.17,0.92)
cb3 <- process_dis(cb3, 170,9.41,0.9)
cb4 <- process_dis(cb4, 160,8.57,0.87)
cb5 <- process_dis(cb5, 165,22.79,0.88)
cb6 <- process_dis(cb6, 180,25.75,0.93)

cn1 <- process_dis(cn1, 235,2.92,1.35)
cn2 <- process_dis(cn2, 235,16.02,1.36)
cn3 <- process_dis(cn3, 235,18.43,1.34)
cn4 <- process_dis(cn4, 235,12.5,1.38)
cn5 <- process_dis(cn5, 230,15.57,1.33)
cn6 <- process_dis(cn6, 230,34.95,1.29)
cn7 <- process_dis(cn7, 200,6.26,1.17)

# ------------------------------
# 1. 整理数据、保留共有 cell
# ------------------------------

# 将所有 21 个数据整理到一个 list 中
data_list <- list(
  ce1 = ce1, ce2 = ce2, ce3 = ce3, ce4 = ce4, ce5 = ce5, ce6 = ce6, ce7 = ce7, ce8 = ce8,
  cb1 = cb1, cb2 = cb2, cb3 = cb3, cb4 = cb4, cb5 = cb5, cb6 = cb6,
  cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cn5 = cn5, cn6 = cn6, cn7 = cn7
)

library(dplyr)

# 找出所有数据中共有的 cell
common_cells <- Reduce(intersect, lapply(data_list, function(df) unique(df$cell)))

# 对每个数据，先只保留共有 cell
data_list <- lapply(data_list, function(df) {
  df %>% 
    filter(cell %in% common_cells)
})

# ------------------------------
# 2. 计算每个 cell 的时间百分比
# ------------------------------

# 对每个数据，按 cell 分组，计算 time_percent
data_list <- lapply(data_list, function(df) {
  df %>% group_by(cell) %>% 
    mutate(time_percent = (time - min(time)) / (max(time) - min(time))) %>% 
    ungroup()
})

# ------------------------------
# 3. 对每个 cell 进行对齐处理
# ------------------------------
# 针对每个 cell，在 21 个数据中选择行数最少的作为模板，
# 再对其它数据中的该 cell，根据模板的 time_percent 找到对应最接近的行

# 定义对齐单个 cell 的函数
align_cell <- function(cell_name, data_list) {
  # 从各数据中提取该 cell 的数据，并按 time_percent 排序
  cell_data <- lapply(data_list, function(df) {
    df %>% 
      filter(cell == cell_name) %>% 
      arrange(time_percent)
  })
  
  # 找出各数据中该 cell 的行数
  n_rows <- sapply(cell_data, nrow)
  # 选择行数最少的作为模板
  template_index <- which.min(n_rows)
  template_data <- cell_data[[template_index]]
  template_time <- template_data$time_percent  # 模板的 time_percent 序列
  
  # 对于每个数据，按模板时间点选择最接近的一行
  aligned_data <- lapply(cell_data, function(df) {
    indices <- sapply(template_time, function(t_val) {
      which.min(abs(df$time_percent - t_val))
    })
    # 将 indices 转换为 numeric 向量
    indices <- unlist(indices)
    df[indices, ]
  })
  return(aligned_data)
}

# 对所有共有的 cell 分别对齐，返回一个 list，每个元素对应一个 cell，
# 每个元素是一个 list，长度与 data_list 相同
aligned_by_cell <- lapply(common_cells, function(cell) {
  align_cell(cell, data_list)
})
names(aligned_by_cell) <- common_cells

# ------------------------------
# 4. 合并所有 cell 的结果
# ------------------------------
# 对于每个数据（共 21 个），按 cell 合并对齐后的数据，
# 保证每个数据中每个 cell 都选取了模板时间的行，行数一致

library(purrr)
# 初始化一个 list 存放对齐后的 21 个数据
aligned_data_list <- vector("list", length(data_list))
names(aligned_data_list) <- names(data_list)

# 对每个数据（索引 i），遍历所有 cell，从 aligned_by_cell 中提取第 i 个数据框，然后 rbind
for (i in seq_along(data_list)) {
  aligned_data_list[[i]] <- map_dfr(aligned_by_cell, ~ .x[[i]])
}


# 可以查看某个数据的结构
str(aligned_data_list[[1]])
summary(aligned_data_list)


# 加载必要的包
library(dplyr)
library(purrr)

# 读取胚胎长度文件（请确认文件路径及分隔符是否正确）
emblength <- read.table("~/Desktop/cbcn/embryo_lengths.txt", header = TRUE, sep = "\t")

# 定义计算每个胚胎数据中各 cell 总移动轨迹距离并归一化的函数
calculate_distance <- function(df, embryo_id) {
  # 如果原数据中没有 embryo 列，则添加 embryo 标识
  df <- df %>%
    mutate(embryo = embryo_id)
  
  df_distance <- df %>%
    group_by(cell) %>%                            # 按 cell 分组
    arrange(time, .by_group = TRUE) %>%             # 每个 cell 内按 time 升序排序
    mutate(
      # 计算连续时间点之间的三维欧氏距离
      distance = sqrt((lead(x) - x)^2 + (lead(y) - y)^2 + (lead(z) - z)^2),
      # 将每组第1行的 distance 设为 NA
      distance = if_else(row_number() == 1, NA_real_, distance)
    ) %>%
    summarise(total_distance = sum(distance, na.rm = TRUE), .groups = "drop")
  
  # 从 emblength 中查找对应胚胎的长度
  embryo_length_value <- emblength$embryo_length[which(emblength$embryo == embryo_id)]
  
  # 增加归一化后的列，归一化公式：total_distance / embryo_length
  df_distance <- df_distance %>%
    mutate(normalized_distance = total_distance / embryo_length_value)
  
  return(df_distance)
}

# 注意：aligned_data_list 是一个有名称的列表，每个元素的名称如 "ce1", "ce2", …, "cn7"
# 对每个胚胎数据依次计算归一化移动距离
results_list <- mapply(function(embryo_id, df) {
  calculate_distance(df, embryo_id)
}, names(aligned_data_list), aligned_data_list, SIMPLIFY = FALSE)

# 接下来对每个胚胎数据提取 cell 与 normalized_distance 列，
# 并将 normalized_distance 列重命名为对应胚胎标识（例如 ce1_norm, ce2_norm, etc.）
norm_list <- lapply(names(results_list), function(embryo_id) {
  results_list[[embryo_id]] %>% 
    select(cell, normalized_distance) %>% 
    rename(!!paste0(embryo_id, "_norm") := normalized_distance)
})

# 利用 purrr::reduce 进行多次 full_join，保证所有胚胎数据按 cell 对齐合并
combined_norm <- reduce(norm_list, full_join, by = "cell")

# 查看合并后的列名称
print(colnames(combined_norm))

# 将最终结果写入 CSV 文件
write.csv(combined_norm, file = "con_cellMovement.csv", row.names = FALSE)













# 设置工作目录，并加载必要的包
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

# 读取数据
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_cellMovement.csv")
colnames(data)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(data, stage, by = "cell")
colnames(data)

data <- data %>%
  rowwise() %>%
  mutate(
    # 对 ce 系列列（如 ce1_norm ~ ce8_norm）计算均值和标准差，再计算 cv
    ce_mean = mean(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_sd   = sd(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_cv   = ce_sd / ce_mean,         # 若要以百分比表示则乘以100
    
    # 对 cb 系列列（cb1_norm ~ cb6_norm）计算均值和标准差，再计算 cv
    cb_mean = mean(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_sd   = sd(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_cv   = cb_sd / cb_mean,         # 如果需要百分比：cb_sd / cb_mean * 100
    
    # 对 cn 系列列（cn1_norm ~ cn7_norm）计算均值和标准差，再计算 cv
    cn_mean = mean(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_sd   = sd(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_cv   = cn_sd / cn_mean          # 或乘以100得到百分比
  ) %>%
  ungroup()%>%
  select(1, 26, 29, 32, 23) %>%
  na.omit()

colnames(data)
tissue <- read.table("~/Desktop/cbcn/linage/AllLineage.tsv", header = TRUE, sep = "\t")
colnames(tissue)
tissue <- tissue %>% select(1, 3)
names(tissue)[1] <- "cell"
colnames(tissue)

data <- merge(data, tissue, by = "cell")
colnames(data)
head(data)

data_long <- data %>%
  pivot_longer(
    cols = c("ce_cv", "cb_cv", "cn_cv"),
    names_to = "measure",
    values_to = "RMSD"
  )
# 2. 修改 factor 标签，使 y 轴显示 ce, cb, cn
data_long$measure <- factor(data_long$measure,
                            levels = c("ce_cv", "cb_cv", "cn_cv"),
                            labels = c("ce", "cb", "cn"))

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
library(ggpubr)  # 用于添加比较显著性标注

# 绘制热图
ggplot(data_long, aes(x = measure, y = RMSD, fill = measure)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  # 添加须线端帽（横线）
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      # 隐藏异常值，避免图形杂乱
               color = "black",         # 箱线边框为黑色
               size = 0.5) +
  #facet_grid(. ~ first_letter, scales = "free_x", space = "free") +
  #facet_grid(. ~ CellFate, scales = "free_x", space = "free") +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  #scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(title = "",
       x = "movement_cv",
       y = "Value") +
  scale_fill_manual(values = c("ce" = "#CD534C",   # 自定义 ce_mean 颜色
                               "cb" = "#0073C2",   # 自定义 cb_mean 颜色
                               "cn" = "#EFC000")) +# 自定义 cn_mean 颜色
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),            # 显示刻度线
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  # 添加差异显著性比较，这里比较 "ce"、"cb" 和 "cn" 三组之间的差异
  stat_compare_means(comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
                     label = "p.signif",   # 使用星号显示显著性结果
                     method = "wilcox.test",
                     label.y = c(0.28, 0.28, 0.28)) 

# # 绘制热图
# ggsave(
#   filename = "movement_cv_box_all.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 6,                            # 图形宽度（单位默认为英寸）
#   height = 6,                            # 图形高度
#   dpi = 300                              # 分辨率
# )






# 设置工作目录，并加载必要的包
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

# 读取数据
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_cellMovement.csv")
colnames(data)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1, 44)
data <- merge(data, stage, by = "cell")
colnames(data)

# 生成 ce_movement_mean, cb_movement_mean, cn_movement_mean 并计算 diff
data <- data %>%
  mutate(ce_movement_mean = round(rowMeans(select(., ce1_norm:ce8_norm), na.rm = TRUE), 2)) %>%
  mutate(cb_movement_mean = round(rowMeans(select(., cb1_norm:cb6_norm), na.rm = TRUE), 2)) %>%
  mutate(cn_movement_mean = round(rowMeans(select(., cn1_norm:cn7_norm), na.rm = TRUE), 2)) %>%
  select(1, 24:26, 23) %>%
  na.omit()
colnames(data)

# 将数据从宽格式转换为长格式
data_long <- data %>%
  pivot_longer(
    cols = c("ce_movement_mean", "cb_movement_mean", "cn_movement_mean"),
    names_to = "experiment",
    values_to = "movement_mean"
  ) %>%
  # 将 experiment 列中的 "_movement_mean" 部分去掉，保留 ce, cb, cn
  mutate(
    experiment = gsub("_movement_mean", "", experiment),
    experiment = factor(experiment, levels = c("ce", "cb", "cn"))
  )

# 绘制箱线图
p_4d <-ggplot(data_long, aes(x = experiment, y = movement_mean, fill = experiment)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  # 添加须线端帽（横线）
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      # 隐藏异常值，避免图形杂乱
               color = "black",         # 箱线边框为黑色
               size = 0.5) +
  scale_fill_manual(values = c("ce" = "#CD534C", 
                               "cb" = "#0073C2", 
                               "cn" = "#EFC000")) +
  labs(
    x = "Experiment",
    y = "Movement Mean",
    title = "Boxplot of Movement_Mean"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),            # 显示刻度线
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0.2, 1.0), expand = c(0.01, 0.01)) +
  scale_x_discrete(expand = c(0.2, 0.2))         # 缩小 x 轴两端的空白区域

print(p_4d)



setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

data <- read.csv("~/Desktop/cbcn/draft_code/new/con_cellMovement.csv")
colnames(data)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(data, stage, by = "cell")
colnames(data)
data <- data %>%
  mutate(ce_movement_mean = round(rowMeans(select(., ce1_norm:ce8_norm), na.rm = TRUE), 2)) %>%
  mutate(cb_movement_mean = round(rowMeans(select(., cb1_norm:cb6_norm), na.rm = TRUE), 2)) %>%
  mutate(cn_movement_mean = round(rowMeans(select(., cn1_norm:cn7_norm), na.rm = TRUE), 2)) %>%
  select(1, 24:26,23) %>%
  na.omit()
colnames(data)

# 加载必要的包
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 定义第一次分组用的 gastrulating cells 名单
gas1 <- c("Ea", "Ep", "Z2", "Z3", 
          "Daa", "Dap", "Dpa", "Dpp", 
          "Capaa", "Capap", "Cappa", "Cappp", 
          "Cppaa", "Cppap", "Cpppa", "Cpppp", 
          "MSaaaa", "MSaaap", "MSapaa", "MSapap", "MSappa", "MSappp", 
          "MSpaaa", "MSpaap", "MSppaa", "MSppap", "MSpppa", "MSpppp",
          "MSaapaa", "MSaapap", "MSaappa", "MSaappp", 
          "MSpapaa", "MSpapap", "MSpappa", "MSpappp",
          "ABalaap", "ABalapp", "ABalpppa", "ABalpppp", 
          "ABarappp", "ABplppap", "ABprpaap", "ABalpaaap", 
          "ABalpaapp", "ABaraapaa", "ABaraapap", "ABaraappa", "ABaraappp", 
          "ABarapaap", "ABalpaaaaa", "ABalpaaaap", "ABalpaapaa", 
          "ABalpaapap", "ABalpappaa", "ABalpapppa", "ABalpapppp", 
          "ABaraaapaa", "ABaraaapap", "ABarapaaaa", "ABarapaaap", 
          "ABarapapaa", "ABarapappa", "ABarapappp", "ABprpapppa", "ABprpapppp")

# 2. 定义第二次分组用的 gastrulating2 名单
gastrulating2 <- c("Ea", "Ep", "Z2", "Z3", 
                   "Daa", "Dap", "Dpa", "Dpp", 
                   "Capaa", "Capap", "Cappa", "Cappp", 
                   "Cppaa", "Cppap", "Cpppa", "Cpppp", 
                   "MSaaaa", "MSaaap", "MSapaa", "MSapap", 
                   "MSappa", "MSappp", "MSpaaa", "MSpaap", 
                   "MSppaa", "MSppap", "MSpppa", "MSpppp",
                   "MSaapaa", "MSaapap", "MSaappa", "MSaappp", 
                   "MSpapaa", "MSpapap", "MSpappa", "MSpappp",
                   "ABalaap", "ABalapp", "ABalpppa", "ABalpppp", 
                   "ABarappp", "ABplppap", "ABprpaap", "ABalpaaap", 
                   "ABalpaapp", "ABaraapaa", "ABaraapap", "ABaraappa", "ABaraappp", 
                   "ABarapaap", "ABalpaaaaa", "ABalpaaaap", "ABalpaapaa", 
                   "ABalpaapap", "ABalpappaa", "ABalpapppa", "ABalpapppp", 
                   "ABaraaapaa", "ABaraaapap", "ABarapaaaa", "ABarapaaap", 
                   "ABarapapaa", "ABarapappa", "ABarapappp", "ABprpapppa", "ABprpapppp")

# 3. 定义第三次分组用的 gastrulating3 名单
gastrulating3 <- c("ABarpapa", "ABarpapp", "ABalaapp", "ABalappa", 
                   "ABalaapa", "ABarapppp", "ABarapppa", "ABpraappp", "ABpraappa")

# 4. 第一次分组：为每个 cell 判断是否在 gas1 中
data_split1 <- data %>%
  mutate(final_group = if_else(cell %in% gas1,
                               "1_gastrulating cells",
                               "1_other cells"),
         split = "split1")  # 标记为第一次分组

# 5. 第二次分组：为每个 cell 判断是否在 gastrulating2 中
data_split2 <- data %>%
  mutate(final_group = if_else(cell %in% gastrulating2,
                               "2_gastrulating cells",
                               "2_other cells"),
         split = "split2")  # 标记为第二次分组

# 6. 第三次分组：为每个 cell 判断是否在 gastrulating3 中
data_split3 <- data %>%
  mutate(final_group = if_else(cell %in% gastrulating3,
                               "3_gastrulating cells",
                               "3_other cells"),
         split = "split3")  # 标记为第三次分组

# 7. 合并三份数据（每个 cell 可能会出现多次，对应不同的分组方案）
combined_data <- bind_rows(data_split1, data_split2, data_split3)

# 8. 将 RMSD 的三个指标转换为长格式，同时保证顺序 ce, cb, cn
data_long <- combined_data %>%
  pivot_longer(
    cols = c("ce_movement_mean", "cb_movement_mean", "cn_movement_mean"),
    names_to = "measure",
    values_to = "RMSD"
  ) %>%
  mutate(measure = recode(measure,
                          "ce_movement_mean" = "ce",
                          "cb_movement_mean" = "cb",
                          "cn_movement_mean" = "cn"),
         measure = factor(measure, levels = c("ce", "cb", "cn"))
  )

# 9. 根据 split 列分别提取三个分组的数据
data_long_split1 <- data_long %>% filter(split == "split1")
data_long_split2 <- data_long %>% filter(split == "split2")
data_long_split3 <- data_long %>% filter(split == "split3")

# 10. 定义一个基础绘图函数，方便重复使用相同的图形主题与样式
base_boxplot <- function(data, title_text) {
  # 定义一个位置调整参数
  dodge_val <- position_dodge(width = 0.75)
  
  ggplot(data, aes(x = final_group, y = RMSD, fill = measure)) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.42,
      size = 0.5,
      position = dodge_val    # 添加位置调整
    ) +
    geom_boxplot(
      width = 0.65, 
      outlier.shape = NA,      # 隐藏异常值
      color = "black",         # 箱线边框为黑色
      size = 0.5,
      position = dodge_val     # 添加位置调整，使箱体与误差条对齐
    ) +
    scale_fill_manual(values = c(
      "ce" = "#CD534C", 
      "cb" = "#0073C2", 
      "cn" = "#EFC000"
    )) +
    labs(
      x = "",
      y = "Mean RMSD",
      title = title_text
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black"),  # 显示刻度线
      axis.text = element_text(size = 12),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0.2, 1.0), expand = c(0, 0))
}
# 11. 分别生成三个分组对应的箱线图
p1 <- base_boxplot(data_long_split1, "Group 1")
p2 <- base_boxplot(data_long_split2, "Group 2")
p3 <- base_boxplot(data_long_split3, "Group 3")

# 12. 分别显示图形（每个图形会在独立的绘图窗口或输出区域显示）
print(p1)
print(p2)
print(p3)


# 加载必要的包
library(ggplot2)
library(gridExtra)
library(grid)


# 将三张图拼接成一个对象（这里横向排列，每行一个图，nrow = 1）
combined_plot <- arrangeGrob(p1, p2, p3, nrow = 1)

# 打印查看合并后的图
grid.newpage()        # 清空已有绘图区域
grid.draw(combined_plot)

# 保存合并后的图为 PNG 文件（宽15英寸，高5英寸，300 dpi）
#ggsave("F4_sup6_movement_gas_all.pdf", combined_plot, width = 10, height = 3.5, dpi = 300)







# 设置工作目录，并加载必要的包
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

# 读取数据
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_cellMovement.csv")
colnames(data)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1, 44)
data <- merge(data, stage, by = "cell")
colnames(data)

# 生成 ce_movement_mean, cb_movement_mean, cn_movement_mean 并计算 diff
data <- data %>%
  mutate(ce_movement_mean = round(rowMeans(select(., ce1_norm:ce8_norm), na.rm = TRUE), 2)) %>%
  mutate(cb_movement_mean = round(rowMeans(select(., cb1_norm:cb6_norm), na.rm = TRUE), 2)) %>%
  mutate(cn_movement_mean = round(rowMeans(select(., cn1_norm:cn7_norm), na.rm = TRUE), 2)) %>%
  select(1, 24:26, 23) %>%
  na.omit()
colnames(data)

data <- data %>%
  mutate(diff_cb = cb_movement_mean - ce_movement_mean,
         diff_cn = cn_movement_mean - ce_movement_mean) %>%
  select(1, 6, 7, 5) 

colnames(data)



# 加载所需的包
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggplotify)  # 用于将 grob 转换为 ggplot 对象
library(cowplot)    # 用于组合 ggplot 对象
library(patchwork)

# 定义注释和颜色（与原代码一致）
ann_colors <- list(
  group = c(A = "#F8766D", 
            M = "#A3A500",
            E = "#00BF7D", 
            C = "#00B0F6",
            D = "#E76BF3")
)

# 1. 第一次分组用的 gastrulating cells 名单
gas1 <- c("Ea", "Ep", "ABalapp", "ABalaap", "ABalpppa", "ABalpppp", "ABarappp", "ABprpaap", "ABplppap", "MSaaaa", "MSaaap", "Z2", "Z3", "MSapaa", "MSapap", "MSappa", "MSappp", "MSpaaa", "MSpaap", "MSppaa", "MSppap", "MSpppa", "MSpppp", "Daa", "Dap", "Dpa", "Dpp", "ABalpaapp", "ABaraappa", "ABaraappp", "ABaraapaa", "ABaraapap", "ABarapaap", "Cpppa", "Cpppp", "ABalpaaap", "MSaapaa", "MSaapap", "MSpapaa", "MSpapap", "Cappa", "Cappp", "Cppaa", "Cppap", "MSaappa", "MSaappp", "Capaa", "Capap", "MSpappa", "MSpappp", "ABaraaapaa", "ABaraaapap", "ABalpaapaa", "ABalpaapap", "ABalpappaa", "ABarapappa", "ABarapappp", "ABalpapppa", "ABalpapppp", "ABarapapaa", "ABprpapppa", "ABprpapppp", 
          "ABarapaaaa", "ABarapaaap", "ABalpaaaaa", "ABalpaaaap")

# 2. 第二次分组用的 gastrulating2 名单
gastrulating2 <- c("Ep", "Ea", "ABalaap", "ABalpppp", "ABalpppa", "ABarappp", "ABplppap", "MSpaaa", "MSaaap", "ABprppap", "MSpppp", "MSapap", "Z2", "MSapaa", "MSppap", "MSappp", "MSpaap", "Z3", "MSppaa", "MSappa", "MSaaaa", "MSpppa", "Dpa", "Daa", "Dpp","ABaraappp", "ABaraapap", "Dap", "ABaraapaa", "ABaraaaap", "Cappa", "ABaraappa", "MSaapap", "MSaapaa", "ABarapaap", "Cppap", "Capaa", "MSpapap", "Cappp", "Cppaa", "MSpapaa", "Cpppp", "Capap", "Cpppa", "ABalpaaap", "ABprpapaa", "MSaappp", "MSpappp", "MSpappa", "MSaappa")

# 3. 第三次分组用的 gastrulating3 名单
gastrulating3 <- c("ABarpapa", "ABarpapp", "ABalaapp", "ABalappa", "ABalaapa", "ABarapppp", "ABarapppa", "ABpraappp", "ABpraappa", "ABalpppap", "ABarppppp", "ABarppppa", "ABalpppaa", "ABalppppp","ABalapppa", "ABalappap", "ABalppppa", "ABalapppp")

# 拆分每个列表为前半部分和后半部分
# 注意：这里假设列表长度为偶数，否则可按 floor()/ceiling() 处理
gas1_first    <- gas1[1:(length(gas1)/2)]
gas1_second   <- gas1[(length(gas1)/2 + 1):length(gas1)]

gas2_first    <- gastrulating2[1:(length(gastrulating2)/2)]
gas2_second   <- gastrulating2[(length(gastrulating2)/2 + 1):length(gastrulating2)]

gas3_first    <- gastrulating3[1:(length(gastrulating3)/2)]
gas3_second   <- gastrulating3[(length(gastrulating3)/2 + 1):length(gastrulating3)]

# 4. 第一次分组：为每个 cell 判断是否在 gas1 中
data_split1 <- data %>%
  mutate(final_group = if_else(cell %in% gas1,
                               "1_gastrulating cells",
                               "1_other cells"),
         split = "split1")  # 标记为第一次分组

# 5. 第二次分组：为每个 cell 判断是否在 gastrulating2 中
data_split2 <- data %>%
  mutate(final_group = if_else(cell %in% gastrulating2,
                               "2_gastrulating cells",
                               "2_other cells"),
         split = "split2")  # 标记为第二次分组

# 6. 第三次分组：为每个 cell 判断是否在 gastrulating3 中
data_split3 <- data %>%
  mutate(final_group = if_else(cell %in% gastrulating3,
                               "3_gastrulating cells",
                               "3_other cells"),
         split = "split3")  # 标记为第三次分组

# 7. 合并三份数据（每个 cell 可能会出现多次，对应不同的分组方案）
combined_data <- bind_rows(data_split1, data_split2, data_split3)

# 8. 将 RMSD 的两个指标转换成长格式，同时保证顺序 ce, cb, cn
data_long <- combined_data %>%
  pivot_longer(
    cols = c("diff_cb", "diff_cn"),
    names_to = "measure",
    values_to = "distance"
  ) %>%
  mutate(
    measure = factor(measure, levels = c("diff_cb", "diff_cn")),
    stage = factor(stage, levels = c("Stage_AB4", "Stage_AB8", "Stage_AB16", 
                                     "Stage_AB32", "Stage_AB64", "Stage_AB128", "Stage_AB256"))
  )

# 9. 根据 split 列分别提取三个分组的数据
data_long_split1 <- data_long %>% filter(split == "split1")
data_long_split2 <- data_long %>% filter(split == "split2")
data_long_split3 <- data_long %>% filter(split == "split3")

# 10. 对每个分组数据按照各自列表进行因子设定
# 此处原代码直接设置因子顺序，现在我们在后续根据拆分好的列表分别做子集
data_long_split1 <- data_long_split1 %>% mutate(cell = factor(cell, levels = gas1))
data_long_split2 <- data_long_split2 %>% mutate(cell = factor(cell, levels = gastrulating2))
data_long_split3 <- data_long_split3 %>% mutate(cell = factor(cell, levels = gastrulating3))

# 11. 将每个分组数据拆分为前半部分和后半部分
data_long_split1_first <- data_long_split1 %>% filter(cell %in% gas1_first) %>% 
  mutate(cell = factor(cell, levels = gas1_first))
data_long_split1_second <- data_long_split1 %>% filter(cell %in% gas1_second) %>% 
  mutate(cell = factor(cell, levels = gas1_second))

data_long_split2_first <- data_long_split2 %>% filter(cell %in% gas2_first) %>% 
  mutate(cell = factor(cell, levels = gas2_first))
data_long_split2_second <- data_long_split2 %>% filter(cell %in% gas2_second) %>% 
  mutate(cell = factor(cell, levels = gas2_second))

data_long_split3_first <- data_long_split3 %>% filter(cell %in% gas3_first) %>% 
  mutate(cell = factor(cell, levels = gas3_first))
data_long_split3_second <- data_long_split3 %>% filter(cell %in% gas3_second) %>% 
  mutate(cell = factor(cell, levels = gas3_second))

# 12. 绘图函数（箱线图）
base_boxplot <- function(data, title_text) {
  # 仅保留 final_group 中包含 "gastrulating cells" 的数据
  data <- data %>% filter(grepl("gastrulating cells", final_group))
  
  dodge <- position_dodge(width = 0.65)
  
  ggplot(data, aes(x = measure, y = distance, fill = measure)) +
    geom_violin(width = 0.7, 
                position = dodge, 
                trim = FALSE, 
                color = "black",
                size = 0.3) +
    geom_boxplot(aes(group = interaction(final_group, measure)), 
                 width = 0.25, 
                 position = dodge,
                 outlier.shape = NA, 
                 color = "black", 
                 size = 0.3,
                 fill = "black",  
                 alpha = 0.2) +
    scale_fill_manual(values = c("diff_cb" = "#d9a190", 
                                 "diff_cn" = "#a6c8cf")) +
    labs(x = NULL, y = "Mean movement", title = title_text) +
    theme_minimal() +
    theme(
      panel.background   = element_blank(),
      panel.grid         = element_blank(),
      plot.background    = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line          = element_blank(),
      axis.ticks         = element_line(color = "black"),
      axis.text          = element_text(size = 12),
      legend.position    = "none"
    ) +
    scale_y_continuous(limits = c(-0.092, 0.28), expand = c(0.001, 0.05))+
    scale_x_discrete(expand = c(0.22, 0.22))+
    # 添加显著性比较：比较 "diff_cb" 与 "diff_cn"，采用 t-test 结果以星号显示
    stat_compare_means(comparisons = list(c("diff_cb", "diff_cn")),
                       label = "p.signif", method = "wilcox.test", size = 5)
  
}

# 13. 分别生成每个拆分后的箱线图
p1_first   <- base_boxplot(data_long_split1_first, "Group 1 - First Half")
p1_second  <- base_boxplot(data_long_split1_second, "Group 1 - Second Half")
p2_first   <- base_boxplot(data_long_split2_first, "Group 2 - First Half")
p2_second  <- base_boxplot(data_long_split2_second, "Group 2 - Second Half")
p3_first   <- base_boxplot(data_long_split3_first, "Group 3 - First Half")
p3_second  <- base_boxplot(data_long_split3_second, "Group 3 - Second Half")

# 14. 合并 6 张图：例如每一行代表一个分组（Group 1, Group 2, Group 3），左右分别为 First Half 和 Second Half
combined_plot <- (p1_first | p1_second) / (p2_first | p2_second) / (p3_first | p3_second)

# 15. 显示组合后的图
print(combined_plot)

# 可选：保存最终组合图到 PDF 文件
#ggsave("movement_gas_group.pdf", combined_plot, width = 10, height = 15)








library(ggplot2)
library(dplyr)
library(ggpubr)

base_barplot <- function(data, title_text) {
  # 过滤只保留包含 "gastrulating cells" 的数据
  data <- data %>% filter(grepl("gastrulating cells", final_group))
  
  # 定义统一的位置调整参数
  dodge <- position_dodge(width = 0.65)
  
  ggplot(data, aes(x = cell, y = distance, fill = measure)) +
    # 用 stat_summary 绘制柱状图（显示均值）
    stat_summary(geom = "col", fun = "mean",
                 position = dodge,
                 width = 0.7,
                 color = "black") +
    scale_fill_manual(values = c("diff_cb" = "#d9a190",
                                 "diff_cn" = "#a6c8cf")) +
    labs(x = NULL,
         y = "Mean movement",
         title = title_text) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line       = element_blank(),
      axis.ticks      = element_line(color = "black"),
      axis.text       = element_text(size = 12),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),  # 将 x 轴标签旋转45度
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(-0.08, 0.26))
}

# 调用示例：
p1 <- base_barplot(data_long_split1, "Group 1")
p2 <- base_barplot(data_long_split2, "Group 2")
p3 <- base_barplot(data_long_split3, "Group 3")

print(p1)
print(p2)
print(p3)

library(patchwork)

# 组合图形：这里用 | 操作符将三个图水平排列，
# 如果想垂直排列可以使用 / 操作符，例如：p_4e / p2 / p3
combined_plot <- p1 / p2 / p3

# 显示组合后的图（在 RStudio 中会自动绘制）
print(combined_plot)


# 保存组合图到文件
ggsave(filename = "movement_diff_gasonly_cell.pdf",
       plot = combined_plot,
       width = 16,    # 根据排列方式调整宽度
       height = 8)








setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_cellMovement.csv")
colnames(data)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(data, stage, by = "cell")
colnames(data)

data <- data %>%
  rowwise() %>%
  mutate(
    # 对 ce 系列列（如 ce1_norm ~ ce8_norm）计算均值和标准差，再计算 cv
    ce_mean = mean(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_sd   = sd(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_cv   = ce_sd / ce_mean,         # 若要以百分比表示则乘以100
    
    # 对 cb 系列列（cb1_norm ~ cb6_norm）计算均值和标准差，再计算 cv
    cb_mean = mean(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_sd   = sd(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_cv   = cb_sd / cb_mean,         # 如果需要百分比：cb_sd / cb_mean * 100
    
    # 对 cn 系列列（cn1_norm ~ cn7_norm）计算均值和标准差，再计算 cv
    cn_mean = mean(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_sd   = sd(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_cv   = cn_sd / cn_mean          # 或乘以100得到百分比
  ) %>%
  ungroup()%>%
  select(1, 26, 29, 32, 23) %>%
  na.omit()

colnames(data)


data_long <- data %>%
  pivot_longer(
    cols = c("ce_cv", "cb_cv", "cn_cv"),
    names_to = "measure",
    values_to = "RMSD"
  )
# 2. 修改 factor 标签，使 y 轴显示 ce, cb, cn
data_long$measure <- factor(data_long$measure,
                            levels = c("ce_cv", "cb_cv", "cn_cv"),
                            labels = c("ce", "cb", "cn"))

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

# 绘制热图
ggplot(data_long, aes(x = cell, y = measure, fill = RMSD)) +
  geom_tile(color = "white") +
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  # scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
  #                      midpoint = 0.17, name = "Z-score") +
  scale_fill_gradientn(
    colors = c("#2166AC", '#92C5DE', "#D1E5F0", "#FDDBC7", "#D6604D", "#B2182B"),
    values = c(0, 0.15, 0.35, 0.5, 0.75, 1)
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


# 绘制热图
ggsave(
  filename = "F4_sup4_movement_cv.pdf",  # 文件名，你也可以给出绝对路径
  plot = last_plot(),                  # 要保存的图形对象
  width = 10,                            # 图形宽度（单位默认为英寸）
  height = 2.5,                            # 图形高度
  dpi = 300                              # 分辨率
)



# 设置工作目录，并加载必要的包
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

# 读取数据
data <- read.csv("~/Desktop/cbcn/draft_code/new/con_cellMovement.csv")
colnames(data)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(data, stage, by = "cell")
colnames(data)

data <- data %>%
  rowwise() %>%
  mutate(
    # 对 ce 系列列（如 ce1_norm ~ ce8_norm）计算均值和标准差，再计算 cv
    ce_mean = mean(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_sd   = sd(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_cv   = ce_sd / ce_mean,         # 若要以百分比表示则乘以100
    
    # 对 cb 系列列（cb1_norm ~ cb6_norm）计算均值和标准差，再计算 cv
    cb_mean = mean(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_sd   = sd(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_cv   = cb_sd / cb_mean,         # 如果需要百分比：cb_sd / cb_mean * 100
    
    # 对 cn 系列列（cn1_norm ~ cn7_norm）计算均值和标准差，再计算 cv
    cn_mean = mean(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_sd   = sd(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_cv   = cn_sd / cn_mean          # 或乘以100得到百分比
  ) %>%
  ungroup()%>%
  select(1, 26, 29, 32, 23) %>%
  na.omit()

colnames(data)
tissue <- read.table("~/Desktop/cbcn/linage/AllLineage.tsv", header = TRUE, sep = "\t")
colnames(tissue)
tissue <- tissue %>% select(1, 3)
names(tissue)[1] <- "cell"
colnames(tissue)

data <- merge(data, tissue, by = "cell")
colnames(data)
head(data)

data_long <- data %>%
  pivot_longer(
    cols = c("ce_cv", "cb_cv", "cn_cv"),
    names_to = "measure",
    values_to = "RMSD"
  )
# 2. 修改 factor 标签，使 y 轴显示 ce, cb, cn
data_long$measure <- factor(data_long$measure,
                            levels = c("ce_cv", "cb_cv", "cn_cv"),
                            labels = c("ce", "cb", "cn"))

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
library(ggpubr)  # 用于添加比较显著性标注

# 绘制热图
ggplot(data_long, aes(x = measure, y = RMSD, fill = measure)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  # 添加须线端帽（横线）
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      # 隐藏异常值，避免图形杂乱
               color = "black",         # 箱线边框为黑色
               size = 0.5) +
  facet_grid(. ~ first_letter, scales = "free_x", space = "free") +
  #facet_grid(. ~ CellFate, scales = "free_x", space = "free") +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  #scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(title = "",
       x = "movement_cv",
       y = "Value") +
  scale_fill_manual(values = c("ce" = "#CD534C",   # 自定义 ce_mean 颜色
                               "cb" = "#0073C2",   # 自定义 cb_mean 颜色
                               "cn" = "#EFC000")) +# 自定义 cn_mean 颜色
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),            # 显示刻度线
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  # 添加差异显著性比较，这里比较 "ce"、"cb" 和 "cn" 三组之间的差异
  stat_compare_means(comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
                     label = "p.signif",   # 使用星号显示显著性结果
                     method = "wilcox.test",
                     label.y = c(0.28, 0.28, 0.28)) 

# # 绘制热图
# ggsave(
#   filename = "movement_cv_box_all.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 6,                            # 图形宽度（单位默认为英寸）
#   height = 6,                            # 图形高度
#   dpi = 300                              # 分辨率
# )





# 绘制热图
ggplot(data_long, aes(x = measure, y = RMSD, fill = measure)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  # 添加须线端帽（横线）
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      # 隐藏异常值，避免图形杂乱
               color = "black",         # 箱线边框为黑色
               size = 0.5) +
  facet_grid(. ~ CellFate, scales = "free_x", space = "free") +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  #scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(title = "",
       x = "movement_cv",
       y = "Value") +
  scale_fill_manual(values = c("ce" = "#CD534C",   # 自定义 ce_mean 颜色
                               "cb" = "#0073C2",   # 自定义 cb_mean 颜色
                               "cn" = "#EFC000")) +# 自定义 cn_mean 颜色
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),            # 显示刻度线
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  # 添加差异显著性比较，这里比较 "ce"、"cb" 和 "cn" 三组之间的差异
  stat_compare_means(comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
                     label = "p.signif",   # 使用星号显示显著性结果
                     method = "wilcox.test",
                     label.y = c(0.28, 0.28, 0.28)) 

# # 绘制热图
# ggsave(
#   filename = "movement_cv_box_all.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 6,                            # 图形宽度（单位默认为英寸）
#   height = 6,                            # 图形高度
#   dpi = 300                              # 分辨率
# )

