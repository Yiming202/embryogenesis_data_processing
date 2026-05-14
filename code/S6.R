


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

# -------------------------------
# 简化后的代码：分块计算两两数据间在相同时间点下每个细胞距离向量对的 RMSD
# -------------------------------
# 定义函数：计算给定距离矩阵列表中所有数据集两两组合的 RMSD
calc_pairwise_RMSD <- function(distance_list, n_splits = 5) {
  dataset_names <- names(distance_list)
  all_results <- list()  # 存放所有结果
  result_idx <- 1
  
  # 利用 combn 获取所有两两组合（保证不重复）
  pairs <- combn(dataset_names, 2, simplify = FALSE)
  
  for (pair in pairs) {
    ds1 <- pair[[1]]
    ds2 <- pair[[2]]
    
    # 取出 ds1 与 ds2 各自的时间点
    timepoints1 <- names(distance_list[[ds1]])
    timepoints2 <- names(distance_list[[ds2]])
    common_tp <- intersect(timepoints1, timepoints2)
    if (length(common_tp) == 0) next  # 若无公共时间点，则跳过
    
    # 分块：将公共时间点分为 n_splits 份，避免一次计算过多
    splits <- split(common_tp, cut(seq_along(common_tp), n_splits, labels = FALSE))
    
    for (chunk in splits) {
      for (tp in chunk) {
        m1 <- distance_list[[ds1]][[tp]]
        m2 <- distance_list[[ds2]][[tp]]
        
        # 检查矩阵是否有行名和列名
        if (is.null(rownames(m1)) || is.null(colnames(m1)) ||
            is.null(rownames(m2)) || is.null(colnames(m2))) {
          warning(paste("时间点", tp, "在", ds1, "或", ds2, "缺少行/列名。"))
          next
        }
        
        # 保留两个矩阵共有的细胞
        common_cells <- intersect(rownames(m1), rownames(m2))
        if (length(common_cells) == 0) next
        
        for (cell in common_cells) {
          # 提取对应距离向量并去除自身的值
          v1 <- m1[cell, ]
          v1 <- v1[names(v1) != cell]
          v2 <- m2[cell, ]
          v2 <- v2[names(v2) != cell]
          
          # 对向量取交集，保证比较的对象一致
          common_vec_names <- intersect(names(v1), names(v2))
          if (length(common_vec_names) == 0) next
          v1 <- v1[common_vec_names]
          v2 <- v2[common_vec_names]
          
          # 计算 RMSD：sqrt(mean((v1 - v2)^2))
          rmsd_val <- sqrt(mean((v1 - v2)^2, na.rm = TRUE))
          
          # 保存结果行
          all_results[[result_idx]] <- data.frame(
            Dataset1  = ds1,
            Dataset2  = ds2,
            TimePoint = tp,
            Cell      = cell,
            RMSD      = rmsd_val,
            stringsAsFactors = FALSE
          )
          result_idx <- result_idx + 1
        } # end for each cell
      } # end for each timepoint in chunk
      
      # 输出当前块处理信息
      cat("Completed chunk for", ds1, ds2, "with timepoints:", paste(chunk, collapse = ", "), "\n")
    }
  }
  
  # 合并所有块结果为一个数据框
  if (length(all_results) > 0) {
    final_df <- do.call(rbind, all_results)
  } else {
    final_df <- data.frame()
  }
  
  return(final_df)
}

# ---------------------------------------------------------------------
# 调用示例
# 对于 ce、cb、cn 三类数据分别计算 RMSD
final_pairwise_ce <- calc_pairwise_RMSD(distance_ce, n_splits = 5)
final_pairwise_cb <- calc_pairwise_RMSD(distance_cb, n_splits = 5)
final_pairwise_cn <- calc_pairwise_RMSD(distance_cn, n_splits = 5)

# ---------------------------------------------------------------------
# 对结果进行分组统计（例如按时间点和细胞求均值）
ce_RMSD_mean_cell_tp <- final_pairwise_ce %>%
  group_by(TimePoint, Cell) %>%
  summarise(mean_RMSD = mean(RMSD, na.rm = TRUE), .groups = "drop")

ce_RMSD_mean_cell <- final_pairwise_ce %>%
  group_by(Cell) %>%
  summarise(mean_RMSD = mean(RMSD, na.rm = TRUE)) %>%
  ungroup()

cb_RMSD_mean_cell_tp <- final_pairwise_cb %>%
  group_by(TimePoint, Cell) %>%
  summarise(mean_RMSD = mean(RMSD, na.rm = TRUE), .groups = "drop")

cb_RMSD_mean_cell <- final_pairwise_cb %>%
  group_by(Cell) %>%
  summarise(mean_RMSD = mean(RMSD, na.rm = TRUE)) %>%
  ungroup()

cn_RMSD_mean_cell_tp <- final_pairwise_cn %>%
  group_by(TimePoint, Cell) %>%
  summarise(mean_RMSD = mean(RMSD, na.rm = TRUE), .groups = "drop")

cn_RMSD_mean_cell <- final_pairwise_cn %>%
  group_by(Cell) %>%
  summarise(mean_RMSD = mean(RMSD, na.rm = TRUE)) %>%
  ungroup()

save(ce_RMSD_mean_cell_tp, ce_RMSD_mean_cell,
     cb_RMSD_mean_cell_tp, cb_RMSD_mean_cell,
     cn_RMSD_mean_cell_tp, cn_RMSD_mean_cell,
     file = "RMSD_results.RData")











setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

load("RMSD_results.RData")

ce <- ce_RMSD_mean_cell %>% 
  rename(ce_mean_RMSD = mean_RMSD)
cb <- cb_RMSD_mean_cell %>% 
  rename(cb_mean_RMSD = mean_RMSD)
cn <- cn_RMSD_mean_cell %>% 
  rename(cn_mean_RMSD = mean_RMSD)

combined_df <- ce %>% 
  left_join(cb, by = "Cell") %>%
  left_join(cn, by = "Cell")
names(combined_df)[1] <- "cell"

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(combined_df, stage, by = "cell")

colnames(data)


# 1. 将宽格式数据转换为长格式，并设置因子顺序及显示名称
data_long <- data %>%
  pivot_longer(
    cols = c(ce_mean_RMSD, cb_mean_RMSD, cn_mean_RMSD),
    names_to = "experiment",
    values_to = "RMSD"
  ) %>%
  mutate(experiment = factor(experiment,
                             levels = c("ce_mean_RMSD", "cb_mean_RMSD", "cn_mean_RMSD"),
                             labels = c("ce", "cb", "cn")))

# 2. 绘制 box plot，并使用 scale_fill_manual 设置指定颜色：
p_4b <-ggplot(data_long, aes(x = experiment, y = RMSD, fill = experiment)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  # 添加须线端帽（横线）
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      # 隐藏异常值，避免图形杂乱
               color = "black",         # 箱线边框为黑色
               size = 0.5) +
  scale_fill_manual(values = c("ce" = "#CD534C", 
                               "cb" = "#0073C2", 
                               "cn" = "#EFC000")) +
  labs(
    x = "",
    y = "Mean RMSD",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),            # 显示刻度线
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0.025, 0.127), expand = c(0.001, 0.001)) +
  scale_x_discrete(expand = c(0.2, 0.2))         # 缩小 x 轴两端的空白区域

print(p_4b)
# ggsave(
#   filename = "pv_compare_cell_all.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 7,                            # 图形宽度（单位默认为英寸）
#   height = 6,                            # 图形高度
#   dpi = 300                              # 分辨率
# )




library(reshape2)
# 1. 将数据转换为长格式，保留 cell 与 xorder 列，同时将三个度量值转换为 long 格式
data_long <- melt(data,
                  id.vars = c("cell",'stage'),
                  measure.vars = c("ce_mean_RMSD", "cb_mean_RMSD", "cn_mean_RMSD"))

# 2. 修改 factor 标签，使 y 轴显示 ce, cb, cn
data_long$variable <- factor(data_long$variable,
                             levels = c("cn_mean_RMSD", "cb_mean_RMSD","ce_mean_RMSD" ),
                             labels = c('Cni', "Cbr", "Cel"))

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
ggplot(data_long, aes(x = cell, y = variable, fill = value)) +
  geom_tile() +  # 使用 tile 几何对象绘制热图
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  labs(x = "Cell", y = "Variable", title = "") +
  scale_fill_gradientn(
    colors = c("#2166AC", '#92C5DE', "#D1E5F0", "#FDDBC7", "#D6604D", "#B2182B"),
    values = c(0, 0.15, 0.3, 0.4, 0.75, 1)
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
#   filename = "pv_compare_cell_stage.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 8,                            # 图形宽度（单位默认为英寸）
#   height = 2,                            # 图形高度
#   dpi = 300                              # 分辨率
# )


# 绘制热图
ggplot(data_long, aes(x = cell, y = variable, fill = value)) +
  geom_tile() +  # 使用 tile 几何对象绘制热图
  facet_grid(. ~ first_letter, scales = "free_x", space = "free") +
  labs(x = "Cell", y = "Variable", title = "") +
  scale_fill_gradientn(
    colors = c("#2166AC", '#92C5DE', "#D1E5F0", "#FDDBC7", "#D6604D", "#B2182B"),
    values = c(0, 0.15, 0.3, 0.4, 0.75, 1)
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
#   filename = "pv_compare_cell_lineage.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 8,                            # 图形宽度（单位默认为英寸）
#   height = 2,                            # 图形高度
#   dpi = 300                              # 分辨率
# )














setwd("~/Desktop/cbcn/draft_code/new/")

library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(magrittr)
load("RMSD_results.RData")

colnames(ce_RMSD_mean_cell_tp)


# 定义需要删除的 Cell 值
#bad_cells <- c("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")
bad_cells <- c("ABa",'ABp','EMS','P2')


# --- 处理 ce_RMSD_mean_cell_tp ---
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  # 先转换 TimePoint 为数值型
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

# 得到需要过滤的 TimePoint 值
bad_timepoints_ce <- ce_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

# 过滤数据并将 TimePoint 减 9
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_ce)) %>%
  mutate(TimePoint = TimePoint - 9.54)

# --- 处理 cb_RMSD_mean_cell_tp ---
cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cb <- cb_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cb)) %>%
  mutate(TimePoint = TimePoint + 1.66)

# --- 处理 cn_RMSD_mean_cell_tp ---
cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cn <- cn_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cn)) %>%
  mutate(TimePoint = TimePoint + 2.28)

# 函数：分割细胞名称
trim_cell <- function(Cell){
  last_cell <- str_sub(Cell, -1)
  remaining_cell <- str_sub(Cell, 1, -2)
  return(c(remaining_cell, last_cell))
}

# 函数：计算x偏移量
xoffset <- function(node, root, klist = list(), slist = list()){
  inter_node <- node
  ROUTE <- c()
  
  while (str_detect(inter_node, root) & inter_node != root) {
    cell <- trim_cell(inter_node)
    inter_node <- cell[1]
    tail_cell <- cell[2]
    
    k <- 1
    if(length(klist) > 0 && sum(str_detect(inter_node, names(klist))) > 0){
      k <- prod(unlist(klist[str_detect(inter_node, names(klist))]))
    }
    
    ROUTE <- c(if(tail_cell == 'a') -k else k, ROUTE)
  }
  
  xoff <- sum((.5^(1:length(ROUTE))) * ROUTE)
  
  if(length(slist) > 0 && sum(str_detect(node, names(slist))) > 0){
    xoff <- xoff + sum(unlist(slist[str_detect(node, names(slist))]))
  }
  
  return(xoff)
}

# 常量定义
klist <- list(Za = 2, Zpap = 0.3, Zppa = 0.5, Zpppa = 0.8, Zppp = .8, Zpppp = .3)
slist <- list(Za = -.5, Zp = -.2, Zpa = .22, Zpap = -.08, Zpp = -.04, Zppa = .04, Zppp = -.06, Zpppa = .01, Zpppp = -.03)
id.show.len <- 4
bg.color <- 'darkgray'

# 读取所有谱系数据
AllLineage.tsv <- read.csv("~/Desktop/cbcn/CellID.csv") %>% tibble()

n_death <- sum(AllLineage.tsv$CellFate == "Death", na.rm = TRUE)
print(n_death)

# 函数：处理谱系数据并绘图
process_lineage_and_plot <- function(lineage_data, std_column, time_limit, 
                                     intercept, slope, normfactor, data_type = "cb", expdata) {
  # 处理谱系数据
  lineage <- lineage_data %>%
    as_tibble() %>%
    filter(time <= time_limit) %>%
    mutate(time = (round((time - intercept) / slope, 2))-normfactor) %>%
    group_by(cell) %>%
    summarize(
      Start = min(time),
      End = max(time)
    ) %>%
    rename(CellName = cell) %>%
    left_join(AllLineage.tsv %>% select(CellName, ID, CellFate), by = "CellName") %>%
    rename(Fate = CellFate) 
  
  # 根据数据类型添加根节点
  if (data_type == "ce") {
    lineage %<>%
      add_row(CellName = "P0", Start = -20, End = -5, ID = 'Z', Fate = 'Unspecified') %>%
      add_row(CellName = "AB", Start = -5, End = 0, ID = 'Za', Fate = 'Unspecified') %>%
      add_row(CellName = "P1", Start = -5, End = 0, ID = 'Zp', Fate = 'Unspecified') %>%
      distinct(CellName, .keep_all = TRUE)
  } else {  # cb 或 cn
    lineage %<>%
      add_row(CellName = "P0", Start = -20, End = -5, ID = 'Z', Fate = 'Unspecified') %>%
      add_row(CellName = "AB", Start = -5, End = 0, ID = 'Za', Fate = 'Unspecified') %>%
      add_row(CellName = "P1", Start = -5, End = 0, ID = 'Zp', Fate = 'Unspecified') %>%
      distinct(CellName, .keep_all = TRUE)
  }
  
  # 计算位置和颜色
  p.data <- lineage %>%
    mutate(End = End + 1) %>%
    mutate(k = if_else(str_detect(ID, 'Zppa'), 0.3, 1)) %>%
    rowwise() %>%
    mutate(x = xoffset(ID, 'Z', klist, slist)) %>%
    mutate(parentx = xoffset(str_sub(ID, end = -2), 'Z', klist, slist))
  
  p.data.y <- p.data %>%
    dplyr::select(CellName, ID, x, Start, End, Fate) %>%
    mutate(Mid = (End + Start) / 2) %>%
    melt(id.vars = c('CellName', 'ID', 'x', 'Mid', 'Fate'),
         variable.name = 'SE', value.name = 'Posi')
  
  p.data.x <- p.data %>% 
    dplyr::select(ID, x, Start, parentx) %>%
    melt(id.vars = c('ID', 'Start'), variable.name = 'FT', value.name = 'Posi')
  
  # 生成 Exp 数据，并新建 UK 列
  Exp.csv <- expdata %>%
    tibble() %>%
    rename(
      CellName = Cell,
      TP = TimePoint,
      EXP = mean_RMSD
    ) %>%
    mutate(UK = TP + 0.96)
  
  Exp.csv %>%
    mutate(EXP = if_else(is.na(EXP), 0, EXP)) %>%
    filter(EXP > 0) %>%
    left_join(AllLineage.tsv %>% dplyr::select(CellName, ID), by = join_by(CellName == CellName)) %>%
    mutate(G = paste0(CellName, TP)) %>%
    rowwise() %>%
    mutate(x = xoffset(ID, 'Z', klist, slist)) -> p.data.exp
  
  p.data.exp %<>%
    bind_rows(p.data.exp %>% mutate(TP = TP + 1))
  
  # 创建图表
  p <- p.data.y %>%
    mutate(label = if_else(nchar(ID) < id.show.len, CellName, NA)) %>%
    ggplot() + 
    geom_path(aes(x = x, y = -Posi, group = ID), color = bg.color) +
    geom_text(aes(x = x, y = -Mid, label = label)) +
    # 这里修改：根据 EXP 列值映射颜色，从蓝色到红色形成渐变
    geom_path(data = p.data.exp, aes(x = x, y = -TP, group = G, color = EXP)) +
    geom_path(data = p.data.x, aes(x = Posi, y = -Start, group = ID), color = bg.color) +
    scale_x_continuous(n.breaks = 20) +
    scale_y_continuous(limits = c(-200, 20), labels = function(x) abs(x)) +
    # 添加根据 EXP 渐变的颜色比例尺（你可根据需要修改 low 和 high 参数）
    # values 参数设置三个颜色在梯度中的位置(0:低值区域, 1:高值区域)
    scale_color_gradientn(
      # 定义三个基本颜色，分别对应蓝色、黄色、红色
      colors = c("blue", "yellow","yellow",'#f55042',"#B21F1F"),
      # 在所属区间 [0.03, 0.068] 内，0.03 对应 0，0.045 对应中间值，0.068 对应 1
      values = c(0, 0.15, 0.4,0.81, 1),
      # 自定义的 rescaler：将值小于 0.03 的统一夹紧为 0.03，值大于 0.068的夹紧为 0.068，
      # 然后将区间 [0.03, 0.068] 线性映射到 [0,1]
      rescaler = function(x, to = c(0,1), from = NULL) {
        x_clamped <- ifelse(x < 0.023, 0.023, ifelse(x > 0.1, 0.1, x))
        scales::rescale(x_clamped, to = to, from = c(0.023, 0.1))
      },
      # 超出范围的值（极端情况）按照 squish 规则处理
      oob = scales::squish
    ) +
    # scale_color_gradientn(
    #   colors = c("blue",'blue', "yellow", "red", "red"),
    #   values = scales::rescale(c(0,0.03,0.045, 0.068, 0.09), to = c(0, 1)),
    #   limits = c(0, 0.09)
    # )
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.y = element_line(color = "black")
    )
  
  return(p)
}

# 处理 CE 数据并绘图
p_ce <- process_lineage_and_plot(
  read.csv("~/Desktop/cbcn/CDFile/ce/CD200113plc1p2.csv"),
  "std_len_ce", 200, 7.04, 1.04, 9.54,data_type = "ce",ce_RMSD_mean_cell_tp
)
print(p_ce)


# 处理 CB 数据并绘图
p_cb <- process_lineage_and_plot(
  read.csv("~/Desktop/cbcn/CDFile/she1/CD241202cbhis72p1.csv"),
  "std_len_cb", 160, 8.57, 0.87, -1.66,data_type = "cb", cb_RMSD_mean_cell_tp
)
print(p_cb)


p_cn <- process_lineage_and_plot(
  read.csv("~/Desktop/cbcn/CDFile/cn/CD241207cnhis72p1.csv"), 
  "std_len_cn", 230, 15.57,1.33, -2.28,data_type = "cn", cn_RMSD_mean_cell_tp
)
print(p_cn)

library(patchwork)
# 使用 patchwork 包将三张图纵向拼接
combined_plot <- p_ce / p_cb / p_cn
print(combined_plot)

#ggsave("PV_combined_lineage_plot_2.pdf", combined_plot, width = 12, height = 7.5)





setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

load("RMSD_results.RData")

colnames(ce_RMSD_mean_cell_tp)

# 定义需要删除的 Cell 值
bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")

# --- 处理 ce_RMSD_mean_cell_tp ---
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  # 先转换 TimePoint 为数值型
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

# 得到需要过滤的 TimePoint 值
bad_timepoints_ce <- ce_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

# 过滤数据并将 TimePoint 减 9
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_ce)) %>%
  mutate(TimePoint = TimePoint - 22.04)

# --- 处理 cb_RMSD_mean_cell_tp ---
cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cb <- cb_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cb)) %>%
  mutate(TimePoint = TimePoint - 13.29)

# --- 处理 cn_RMSD_mean_cell_tp ---
cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cn <- cn_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cn)) %>%
  mutate(TimePoint = TimePoint - 9.1)

# 计算 cn_RMSD_mean_cell_tp 中 TimePoint 的最大值
max_cn_time <- max(cn_RMSD_mean_cell_tp$TimePoint, na.rm = TRUE)

# 过滤 ce 数据框中 TimePoint 大于 cn 中最大值的行
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

# 过滤 cb 数据框中 TimePoint 大于 cn 中最大值的行
cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

sorted_tp <- sort(as.numeric(names(table(cb_RMSD_mean_cell_tp$TimePoint))))
print(sorted_tp)

# # 检查每个数据集中唯一 TimePoint 的数量和范围
# all_data %>%
#   group_by(Group) %>%
#   summarise(
#     n_unique = n_distinct(TimePoint),
#     min_time = min(TimePoint),
#     max_time = max(TimePoint)
#   )


# 为每个数据框添加分组信息
ce_RMSD_mean_cell_tp$Group <- "ce"
cb_RMSD_mean_cell_tp$Group <- "cb"
cn_RMSD_mean_cell_tp$Group <- "cn"

# 合并三个数据框
all_data <- bind_rows(ce_RMSD_mean_cell_tp, cb_RMSD_mean_cell_tp, cn_RMSD_mean_cell_tp)

library(dplyr)

result <- all_data %>%
  filter(Group == "ce") %>%                   # 筛选 Group 为 "ce" 的数据
  group_by(TimePoint) %>%                     # 按 TimePoint 分组
  summarise(nCells = n_distinct(Cell)) %>%    # 统计每个 TimePoint 下不同 Cell 的个数
  filter(nCells == 340)                       # 只保留 Cell 数量为 26 的 TimePoint
##28cell, 170 cell, 300cell, 509cell
print(result)
cetimepoints <- result$TimePoint
cbtimepoints <- result$TimePoint
cntimepoints <- result$TimePoint
print(cetimepoints)
print(cbtimepoints)
print(cntimepoints)

# 计算所有数据中 TimePoint 的最小值和最大值，用于统一 x 轴范围
min_time <- min(all_data$TimePoint) - 1
max_time <- max(all_data$TimePoint) + 1

# # 绘制箱线图：
p_4c <- ggplot(all_data, aes(x = TimePoint, y = mean_RMSD, color = Group, fill = Group)) +
  # 绘制箱体：设置 alpha = 0 使填充透明，同时绘制箱型图轮廓和中位数线
  geom_boxplot(aes(group = interaction(TimePoint, Group)),
               outlier.shape = NA,
               width = 1,              # 调整箱体宽度
               position = position_dodge(width = 0),
               alpha = 0,                # 填充设为透明
               coef = 0,                 # 上下 whisker 与箱体边界重合
               size = 0.2) +               # 上下 whisker 与箱体边界重合
  # 第二层：添加 whisker 并设置半透明效果（alpha = 0.4）
  stat_boxplot(geom = "errorbar",
               width = 0.2,
               aes(group = interaction(TimePoint, Group)),
               position = position_dodge(width = 0),
               coef = 1.3,    # 使用默认的 1.5 倍 IQR 计算 whisker 端点
               alpha = 0.4,                 # 上下 whisker 与箱体边界重合
               linewidth = 0.3) +
  # 添加 x=12.5 的竖直虚线
  geom_vline(xintercept = 20.23, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 68.3, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 98.11, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 143.31, linetype = "dashed", color = "black", size = 0.5) +
  labs(x = "TimePoint", y = "RMSD", color = "Group") +
  scale_y_continuous(limits = c(0.02, 0.13), expand = c(0, 0)) +
  scale_x_continuous(limits = c(min_time, max_time),
                     expand = c(0, 0)) +
  # 设置不同组别对应的边框和填充颜色（填充由 alpha=0 达到透明效果）
  scale_color_manual(values = c("ce" = "#CD534C", 
                                "cb" = "#0073C2", 
                                "cn" = "#EFC000")) +
  scale_fill_manual(values = c("ce" = "#CD534C", 
                               "cb" = "#0073C2", 
                               "cn" = "#EFC000")) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),            # 显示刻度线
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) 


print(p_4c)

#ggsave(filename = "pv_box1.pdf", plot = p_4c, device = "pdf", width = 8, height = 7, units = "in")







setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

load("RMSD_results.RData")

colnames(ce_RMSD_mean_cell_tp)

# 定义需要删除的 Cell 值
bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")

# --- 处理 ce_RMSD_mean_cell_tp ---
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  # 先转换 TimePoint 为数值型
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

# 得到需要过滤的 TimePoint 值
bad_timepoints_ce <- ce_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

# 过滤数据并将 TimePoint 减 9
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_ce)) %>%
  mutate(TimePoint = TimePoint - 22.04)

# --- 处理 cb_RMSD_mean_cell_tp ---
cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cb <- cb_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cb)) %>%
  mutate(TimePoint = TimePoint - 13.29)

# --- 处理 cn_RMSD_mean_cell_tp ---
cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cn <- cn_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cn)) %>%
  mutate(TimePoint = TimePoint - 9.1)

# 计算 cn_RMSD_mean_cell_tp 中 TimePoint 的最大值
max_cn_time <- max(cn_RMSD_mean_cell_tp$TimePoint, na.rm = TRUE)

# 过滤 ce 数据框中 TimePoint 大于 cn 中最大值的行
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

# 过滤 cb 数据框中 TimePoint 大于 cn 中最大值的行
cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

sorted_tp <- sort(as.numeric(names(table(cb_RMSD_mean_cell_tp$TimePoint))))
print(sorted_tp)

# # 检查每个数据集中唯一 TimePoint 的数量和范围
# all_data %>%
#   group_by(Group) %>%
#   summarise(
#     n_unique = n_distinct(TimePoint),
#     min_time = min(TimePoint),
#     max_time = max(TimePoint)
#   )


# 为每个数据框添加分组信息
ce_RMSD_mean_cell_tp$Group <- "ce"
cb_RMSD_mean_cell_tp$Group <- "cb"
cn_RMSD_mean_cell_tp$Group <- "cn"

# 合并三个数据框
all_data <- bind_rows(ce_RMSD_mean_cell_tp, cb_RMSD_mean_cell_tp, cn_RMSD_mean_cell_tp)
# 提取 Cell 列的第一个字母，并将字母 "Z" 归为 "P"
all_data <- all_data %>%
  mutate(Facet = if_else(substr(Cell, 1, 1) == "Z", "P", substr(Cell, 1, 1)))%>%
  filter(Facet != "P")%>%
  mutate(Facet = factor(Facet, levels = c("A", "E", "M", "C", "D")))

# 计算所有数据中 TimePoint 的最小值和最大值，用于统一 x 轴范围
min_time <- min(all_data$TimePoint) - 1
max_time <- max(all_data$TimePoint) + 1

# 绘制箱线图，并根据新变量 Facet 进行分面绘制：
p_4c <- ggplot(all_data, aes(x = TimePoint, y = mean_RMSD, color = Group, fill = Group)) +
  # 绘制箱体：设置 alpha = 0 使填充透明，同时绘制箱型图轮廓和中位数线
  geom_boxplot(aes(group = interaction(TimePoint, Group)),
               outlier.shape = NA,
               width = 1,              # 调整箱体宽度
               position = position_dodge(width = 0),
               alpha = 0,              # 填充设为透明
               coef = 0,               # whisker 与箱体边界重合
               size = 0.2) +
  # 第二层：添加 whisker 并设置半透明效果（alpha = 0.4）
  stat_boxplot(geom = "errorbar",
               width = 0.2,
               aes(group = interaction(TimePoint, Group)),
               position = position_dodge(width = 0),
               coef = 1.3,
               alpha = 0.4,
               linewidth = 0.3) +
  #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  # 添加 x=12.5 的竖直虚线
  geom_vline(xintercept = 20.23, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 68.3, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 98.11, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 143.31, linetype = "dashed", color = "black", size = 0.5) +
  labs(x = "TimePoint", y = "RMSD", color = "Group") +
  scale_y_continuous(limits = c(0.02, 0.13), expand = c(0, 0)) +
  scale_x_continuous(limits = c(min_time, max_time), expand = c(0, 0)) +
  # 设置不同组别对应的颜色
  scale_color_manual(values = c("ce" = "#CD534C", 
                                "cb" = "#0073C2", 
                                "cn" = "#EFC000")) +
  scale_fill_manual(values = c("ce" = "#CD534C", 
                               "cb" = "#0073C2", 
                               "cn" = "#EFC000")) +
  theme_classic(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    legend.position = "none",
    strip.background = element_blank()  # 移除分面标题边框
  ) +
  # 根据提取的首字母分面绘制
  facet_wrap(~ Facet, nrow = 1)

# 显示图形
print(p_4c)
#ggsave(filename = "pv_box_lineage.pdf", plot = p_4c, device = "pdf", width = 16, height = 3.8)







# ============= 预处理数据 =============
setwd("~/Desktop/cbcn/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggpubr)

load("RMSD_results.RData")
setwd("~/Desktop/cbcn/draft_code/fig4_all/")

# 检查 ce_RMSD_mean_cell_tp 的列名
colnames(ce_RMSD_mean_cell)

# 定义需要删除的 Cell 值
bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")


# 为每个数据框添加分组信息
ce_RMSD_mean_cell$Group <- "ce"
cb_RMSD_mean_cell$Group <- "cb"
cn_RMSD_mean_cell$Group <- "cn"

# 合并三个数据框
all_data <- bind_rows(ce_RMSD_mean_cell, cb_RMSD_mean_cell, cn_RMSD_mean_cell)
colnames(all_data)

# 将 Group 转换为因子，并按照 ce, cb, cn 的顺序排列
all_data <- all_data %>%
  mutate(Group = factor(Group, levels = c("ce", "cb", "cn")))

colnames(all_data)


# 读取左右分组数据
lr <- read.csv("~/Desktop/cbcn/LR_symmetry_2.csv")
colnames(lr)

lr_long <- lr %>%
  # 为每一行配对生成一个 pair id（例如，1,2,3,...）
  mutate(pair = row_number()) %>%
  # 将 L 与 R 转换成长格式
  pivot_longer(cols = c("L", "R"), names_to = "side", values_to = "Cell") %>%
  # 根据 side 决定 lr 取值：若 side 为 "L" 则为 "l"，若 "R" 则为 "r"
  mutate(lr = if_else(side == "L", "l", "r")) %>%
  # 可以去掉 side 列
  select(-side)

# 2. 将 lr_long 与 all_data 按照 Cell 进行左连接，加入 lr 与 pair 信息
all_data <- all_data %>%
  left_join(lr_long, by = "Cell") %>%
  # 对于未匹配到 lr 信息的行，将 lr 设置为 "n"
  mutate(lr = if_else(is.na(lr), "n", lr))

colnames(all_data)

# ===================== 以下部分为绘图处理 =====================

# 1. 计算pair层面的平均值
#    仅考虑左右配对数据，即 lr 为 "l" 或 "r" 的记录；如果某个 pair 对应多次测量，则对每侧取平均
data_pair <- all_data %>%
  filter(lr %in% c("l", "r")) %>% 
  group_by(Group, pair, lr) 

# 2. 将左右数据转换为宽格式，每个 pair 一行，得到左右侧的平均值
data_pair_wide <- data_pair %>%
  pivot_wider(names_from = lr, values_from = c(mean_RMSD, Cell)) %>%
  # 此时 cols为：Group, pair, avg_RMSD_l, avg_RMSD_r, cell_l, cell_r
  # 仅保留 l 和 r 都存在的 pair
  filter(!is.na(mean_RMSD_l) & !is.na(mean_RMSD_r)) %>%
  mutate(diff = mean_RMSD_l - mean_RMSD_r)

# -------------------- 作图 --------------------

## (A) 作图①：显示每个pair中左右两侧的平均 RMSD 值。
#     用点图展示左右两侧，用线连接同一 pair 内的两侧数据，
#     并按 Group 分面显示
p_pair_box <- ggplot(data_pair, aes(x = lr, y = mean_RMSD)) +
  # stat_boxplot(geom = "errorbar", width = 0.4) +
  # 为 errorbar 添加 color 映射（使用 aes(color = lr)）
  stat_boxplot(aes(color = lr), geom = "errorbar", width = 0.5,size = 0.8) +
  # 绘制点图层，固定所有点为透明的灰色，组内使用 pair 连接
  geom_point(aes(group = factor(pair)), 
             color = "grey", 
             size = 1.2, 
             alpha = 0.4) +
  # 添加连线层，连接同一个 pair 内的数据点
  geom_line(aes(group = factor(pair)), color = "gray", alpha = 0.15) +
  # 绘制箱型图，映射边框颜色到 lr，固定填充颜色为灰色，半透明显示
  geom_boxplot(aes(color = lr),
               fill = "transparent", 
               outlier.shape = NA, 
               width = 0.65, 
               alpha = 0.2,
               size = 0.6) +
  scale_color_manual(values = c("l" = "#e3882f", "r" = "#49548a")) +
  facet_wrap(~ Group) +
  labs(title = "Paired RMSD for each pair (l vs r)",
       x = "Side",
       y = "Average RMSD") +
  theme_classic() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_blank(),  # 移除分面标题边框
        # 删除原本的 xy 轴线
        axis.line = element_blank(),
        # 为每个 facet 添加四边框
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  # 添加左右组（l vs r）统计比较标记
  stat_compare_means(comparisons = list(c("l", "r")),
                     method = "t.test",      # 可选 "t.test" 或 "wilcox.test"
                     label = "p.signif",     # 使用 p.signif 显示星级
                     hide.ns = TRUE,         # 不显示不显著的标记
                     # label.y 设置位置，这里取所有数据的最大值的1.1倍，必要时可调整
                     label.y = max(data_pair$mean_RMSD, na.rm = TRUE) * 0.9)

print(p_pair_box)

#ggsave("pv_lr_pair.pdf", plot = p_pair_box, width = 5, height = 2.8)



# 
# bidirectional_data <- data_pair_wide %>%
#   pivot_longer(
#     cols = c(mean_RMSD_l, mean_RMSD_r), 
#     names_to = "side", 
#     values_to = "RMSD"
#   ) %>%
#   mutate(
#     # 将变量名称转换为 “l” 或 “r”
#     side = if_else(side == "mean_RMSD_l", "l", "r"),
#     # 对左侧（l）的 RMSD 取负，以便在图中显示在下方
#     RMSD = if_else(side == "l", -RMSD, RMSD)
#   )
# 
# # -------------------- 绘制双向柱状图 --------------------
# p_bidirectional <- ggplot(bidirectional_data, aes(x = factor(pair), y = RMSD, fill = side)) +
#   geom_col(width = 0.7) +
#   facet_wrap(~ Group, scales = "free_x") +
#   # 使用绝对值作为 y 轴标签，使得下方数值也以正数显示
#   scale_y_continuous(labels = abs) +
#   labs(
#     title = "Bidirectional RMSD for each pair",
#     x = "Pair ID",
#     y = "Average RMSD"
#   ) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_blank(),   # 隐藏 x 轴刻度文本
#     strip.text = element_text(size = 12)
#   )
# 
# print(p_bidirectional)
# 
# # 保存图形
# ggsave("pv_bidirectional_pair.pdf", plot = p_bidirectional, width = 10, height = 6)




## (B) 作图②：展示每对 pair 的差值 (l - r)
# (B1) 用条形图显示每个 pair 的差值（按 Group 分面）；
p_diff_individual <- ggplot(data_pair_wide, aes(x = factor(pair), y = diff, fill = Group)) +
  geom_col() +
  facet_wrap(~ Group, scales = "free_x") +
  labs(title = "Difference (l - r) for each pair",
       x = "Pair ID",
       y = "Difference (l - r)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),  # 隐藏 x 轴刻度文本
        strip.text = element_text(size = 12))
print(p_diff_individual)
ggsave("pv_lr_pair_value.pdf", plot = p_diff_individual, width = 10, height = 6)

# p_diff_individual_points <- ggplot(data_pair_wide, aes(x = factor(pair), y = diff, color = Group)) +
#   geom_point(size = 1, position = position_jitter(width = 0.2)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
#   facet_wrap(~ Group, scales = "free_x") +
#   labs(title = "Difference (l - r) for each pair",
#        x = "Pair ID",
#        y = "Difference (l - r)") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.text = element_text(size = 12))
# 
# print(p_diff_individual_points)


# (B2) 用箱型图展示各 Group 内所有 pair 差值的分布
p_diff_overall <- ggplot(data_pair_wide, aes(x = Group, y = abs(diff), fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.42, linewidth = 0.5) +  # 添加须线端帽（横线）
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      # 隐藏异常值，避免图形杂乱
               color = "black",         # 箱线边框为黑色
               size = 0.5) +
  scale_fill_manual(values = c("ce" = "#CD534C", 
                               "cb" = "#0073C2", 
                               "cn" = "#EFC000")) +
  labs(title = "Overall difference (l - r) across pairs",
       x = "Group",
       y = "Difference (l - r)") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),            # 显示刻度线
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
    method = "t.test",
    label = "p.signif",    # 以星号标记显著性
    hide.ns = TRUE,
    label.y = max(abs(data_pair_wide$diff), na.rm = TRUE) * 0.57
  )+
  scale_y_continuous(limits = c(0, 0.045), expand = c(0.0005, 0.001))
print(p_diff_overall)

ggsave("pv_lr_pair_diff_box.pdf", plot = p_diff_overall, width = 10, height = 10)

# ---------------- 合并绘图 ----------------
# 将图①和图②进行拼接排列（可以根据实际需求调整排列方式）
combined_diff <- (p_pair | p_diff_individual) / p_diff_overall
print(combined_diff)

# 保存最终图形（宽、高可根据需要调整）
#ggsave("pair_difference_plots.pdf", plot = combined_diff, width = 12, height = 10)







setwd("~/Desktop/cbcn/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

load("RMSD_results.RData")
setwd("~/Desktop/cbcn/draft_code/fig4_all/")

colnames(ce_RMSD_mean_cell_tp)

# 定义需要删除的 Cell 值
bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")

# --- 处理 ce_RMSD_mean_cell_tp ---
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  # 先转换 TimePoint 为数值型
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

# 得到需要过滤的 TimePoint 值
bad_timepoints_ce <- ce_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

# 过滤数据并将 TimePoint 减 9
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_ce)) %>%
  mutate(TimePoint = TimePoint - 22.04)

# --- 处理 cb_RMSD_mean_cell_tp ---
cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cb <- cb_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cb)) %>%
  mutate(TimePoint = TimePoint - 13.29)

# --- 处理 cn_RMSD_mean_cell_tp ---
cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cn <- cn_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cn)) %>%
  mutate(TimePoint = TimePoint - 9.1)

# 计算 cn_RMSD_mean_cell_tp 中 TimePoint 的最大值
max_cn_time <- max(cn_RMSD_mean_cell_tp$TimePoint, na.rm = TRUE)

# 过滤 ce 数据框中 TimePoint 大于 cn 中最大值的行
ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

# 过滤 cb 数据框中 TimePoint 大于 cn 中最大值的行
cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

sorted_tp <- sort(as.numeric(names(table(cb_RMSD_mean_cell_tp$TimePoint))))
print(sorted_tp)

# # 检查每个数据集中唯一 TimePoint 的数量和范围
# all_data %>%
#   group_by(Group) %>%
#   summarise(
#     n_unique = n_distinct(TimePoint),
#     min_time = min(TimePoint),
#     max_time = max(TimePoint)
#   )


# 为每个数据框添加分组信息
ce_RMSD_mean_cell_tp$Group <- "ce"
cb_RMSD_mean_cell_tp$Group <- "cb"
cn_RMSD_mean_cell_tp$Group <- "cn"

# 合并三个数据框
all_data <- bind_rows(ce_RMSD_mean_cell_tp, cb_RMSD_mean_cell_tp, cn_RMSD_mean_cell_tp)
colnames(all_data)
lr    <- read.csv("~/Desktop/cbcn/LR_symmetry_2.csv")
colnames(lr)

all_data <- all_data %>%
  mutate(lr = case_when(
    Cell %in% lr$L ~ "l",
    Cell %in% lr$R ~ "r",
    TRUE ~ "n"
  ))
colnames(all_data)

# 计算所有数据中 TimePoint 的最小值和最大值，用于统一 x 轴范围
min_time <- min(all_data$TimePoint) - 1
max_time <- max(all_data$TimePoint) + 1

# 绘制箱线图，并根据新变量 Facet 进行分面绘制：
p_4c <- ggplot(all_data, aes(x = TimePoint, y = mean_RMSD, color = Group, fill = Group)) +
  # 绘制箱体：设置 alpha = 0 使填充透明，同时绘制箱型图轮廓和中位数线
  geom_boxplot(aes(group = interaction(TimePoint, Group)),
               outlier.shape = NA,
               width = 1,              # 调整箱体宽度
               position = position_dodge(width = 0),
               alpha = 0,              # 填充设为透明
               coef = 0,               # whisker 与箱体边界重合
               size = 0.2) +
  # 第二层：添加 whisker 并设置半透明效果（alpha = 0.4）
  stat_boxplot(geom = "errorbar",
               width = 0.2,
               aes(group = interaction(TimePoint, Group)),
               position = position_dodge(width = 0),
               coef = 1.3,
               alpha = 0.4,
               linewidth = 0.3) +
  labs(x = "TimePoint", y = "RMSD", color = "Group") +
  scale_y_continuous(limits = c(0.02, 0.13), expand = c(0, 0)) +
  scale_x_continuous(limits = c(min_time, max_time), expand = c(0, 0)) +
  # 设置不同组别对应的颜色
  scale_color_manual(values = c("ce" = "#CD534C", 
                                "cb" = "#0073C2", 
                                "cn" = "#EFC000")) +
  scale_fill_manual(values = c("ce" = "#CD534C", 
                               "cb" = "#0073C2", 
                               "cn" = "#EFC000")) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  # 根据提取的首字母分面绘制
  facet_wrap(~ factor(lr,levels = c('l','r','n')))

# 显示图形
print(p_4c)
#ggsave(filename = "pv_box_lineage.pdf", plot = p_4c, device = "pdf", width = 10, height = 6)






# ============= 预处理数据 =============
setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggpubr)

load("RMSD_results.RData")

# 检查 ce_RMSD_mean_cell_tp 的列名
colnames(ce_RMSD_mean_cell)

# 定义需要删除的 Cell 值
bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")


# 为每个数据框添加分组信息
ce_RMSD_mean_cell$Group <- "ce"
cb_RMSD_mean_cell$Group <- "cb"
cn_RMSD_mean_cell$Group <- "cn"

# 合并三个数据框
all_data <- bind_rows(ce_RMSD_mean_cell, cb_RMSD_mean_cell, cn_RMSD_mean_cell)
colnames(all_data)

# 将 Group 转换为因子，并按照 ce, cb, cn 的顺序排列
all_data <- all_data %>%
  mutate(Group = factor(Group, levels = c("ce", "cb", "cn")))

colnames(all_data)


# 读取左右分组数据
lr <- read.csv("~/Desktop/cbcn/LR_symmetry_2.csv")
colnames(lr)

lr_long <- lr %>%
  # 为每一行配对生成一个 pair id（例如，1,2,3,...）
  mutate(pair = row_number()) %>%
  # 将 L 与 R 转换成长格式
  pivot_longer(cols = c("L", "R"), names_to = "side", values_to = "Cell") %>%
  # 根据 side 决定 lr 取值：若 side 为 "L" 则为 "l"，若 "R" 则为 "r"
  mutate(lr = if_else(side == "L", "l", "r")) %>%
  # 可以去掉 side 列
  select(-side)

# 2. 将 lr_long 与 all_data 按照 Cell 进行左连接，加入 lr 与 pair 信息
all_data <- all_data %>%
  left_join(lr_long, by = "Cell") %>%
  # 对于未匹配到 lr 信息的行，将 lr 设置为 "n"
  mutate(lr = if_else(is.na(lr), "n", lr))

colnames(all_data)

# ===================== 以下部分为绘图处理 =====================

# 1. 计算pair层面的平均值
#    仅考虑左右配对数据，即 lr 为 "l" 或 "r" 的记录；如果某个 pair 对应多次测量，则对每侧取平均
data_pair <- all_data %>%
  filter(lr %in% c("l", "r")) %>% 
  group_by(Group, pair, lr) 

# 2. 将左右数据转换为宽格式，每个 pair 一行，得到左右侧的平均值
data_pair_wide <- data_pair %>%
  pivot_wider(names_from = lr, values_from = c(mean_RMSD, Cell)) %>%
  # 此时 cols为：Group, pair, avg_RMSD_l, avg_RMSD_r, cell_l, cell_r
  # 仅保留 l 和 r 都存在的 pair
  filter(!is.na(mean_RMSD_l) & !is.na(mean_RMSD_r)) %>%
  mutate(diff = mean_RMSD_l - mean_RMSD_r)

# -------------------- 作图 --------------------

## (A) 作图①：显示每个pair中左右两侧的平均 RMSD 值。
#     用点图展示左右两侧，用线连接同一 pair 内的两侧数据，
#     并按 Group 分面显示
p_pair_box <- ggplot(data_pair, aes(x = lr, y = mean_RMSD)) +
  # stat_boxplot(geom = "errorbar", width = 0.4) +
  # 为 errorbar 添加 color 映射（使用 aes(color = lr)）
  stat_boxplot(aes(color = lr), geom = "errorbar", width = 0.5,size = 0.8) +
  # 绘制点图层，固定所有点为透明的灰色，组内使用 pair 连接
  geom_point(aes(group = factor(pair)), 
             color = "grey", 
             size = 1.2, 
             alpha = 0.4) +
  # 添加连线层，连接同一个 pair 内的数据点
  geom_line(aes(group = factor(pair)), color = "gray", alpha = 0.15) +
  # 绘制箱型图，映射边框颜色到 lr，固定填充颜色为灰色，半透明显示
  geom_boxplot(aes(color = lr),
               fill = "transparent", 
               outlier.shape = NA, 
               width = 0.65, 
               alpha = 0.2,
               size = 0.6) +
  scale_color_manual(values = c("l" = "#e3882f", "r" = "#49548a")) +
  facet_wrap(~ Group) +
  labs(title = "Paired RMSD for each pair (l vs r)",
       x = "Side",
       y = "Average RMSD") +
  theme_classic() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_blank(),  # 移除分面标题边框
        # 删除原本的 xy 轴线
        axis.line = element_blank(),
        # 为每个 facet 添加四边框
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  # 添加左右组（l vs r）统计比较标记
  stat_compare_means(comparisons = list(c("l", "r")),
                     method = "t.test",      # 可选 "t.test" 或 "wilcox.test"
                     label = "p.signif",     # 使用 p.signif 显示星级
                     hide.ns = TRUE,         # 不显示不显著的标记
                     # label.y 设置位置，这里取所有数据的最大值的1.1倍，必要时可调整
                     label.y = max(data_pair$mean_RMSD, na.rm = TRUE) * 0.9)

print(p_pair_box)

#ggsave("pv_lr_pair.pdf", plot = p_pair_box, width = 5, height = 2.8)





## (B) 作图②：展示每对 pair 的差值 (l - r)
# (B1) 用条形图显示每个 pair 的差值（按 Group 分面）；
p_diff_individual <- ggplot(data_pair_wide, aes(x = factor(pair), y = diff, fill = Group)) +
  geom_col() +
  facet_wrap(~ Group, scales = "free_x") +
  labs(title = "Difference (l - r) for each pair",
       x = "Pair ID",
       y = "Difference (l - r)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),  # 隐藏 x 轴刻度文本
        strip.text = element_text(size = 12))
print(p_diff_individual)
#ggsave("pv_lr_pair_value.pdf", plot = p_diff_individual, width = 10, height = 6)




# (B2) 用箱型图展示各 Group 内所有 pair 差值的分布
p_diff_overall <- ggplot(data_pair_wide, aes(x = Group, y = abs(diff), fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.42, linewidth = 0.5) +  # 添加须线端帽（横线）
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      # 隐藏异常值，避免图形杂乱
               color = "black",         # 箱线边框为黑色
               size = 0.5) +
  scale_fill_manual(values = c("ce" = "#CD534C", 
                               "cb" = "#0073C2", 
                               "cn" = "#EFC000")) +
  labs(title = "Overall difference (l - r) across pairs",
       x = "Group",
       y = "Difference (l - r)") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),            # 显示刻度线
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
    method = "t.test",
    label = "p.signif",    # 以星号标记显著性
    hide.ns = TRUE,
    label.y = max(abs(data_pair_wide$diff), na.rm = TRUE) * 0.57
  )+
  scale_y_continuous(limits = c(0, 0.045), expand = c(0.0005, 0.001))
print(p_diff_overall)

#ggsave("pv_lr_pair_diff_box.pdf", plot = p_diff_overall, width = 10, height = 10)




# ==========================================================
# Figure 5 & Figure 2(new_F2): pv (RMSD mean) vs asynchrony
# ==========================================================

setwd("~/Desktop/cbcn/draft_code/new/")

# ----------------------------------------------------------
# 1. Load libraries
# ----------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)   # for stat_cor()

# ----------------------------------------------------------
# 2. Load asynchrony datasets
# ----------------------------------------------------------
ce <- read.table("~/Desktop/cbcn/output_folder/test/cedivision_asynchrony_phenotype/ce_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t") %>% rename(ce_mean = cell_rmsd_mean)
cb <- read.table("~/Desktop/cbcn/output_folder/test/cbdivision_asynchrony_phenotype/cb_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t") %>% rename(cb_mean = cell_rmsd_mean)
cn <- read.table("~/Desktop/cbcn/output_folder/test/cndivision_asynchrony_phenotype/cn_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t") %>% rename(cn_mean = cell_rmsd_mean)

asy_data <- ce %>%
  left_join(cb, by = "all_cells") %>%
  left_join(cn, by = "all_cells") %>%
  rename(cell = all_cells)

# stage / tissue annotation
stage <- read.csv("~/Desktop/cbcn_draft/file/nor_Length_Time_correct.csv") %>%
  select(cell = 1, stage = 44)
tissue <- read.table("~/Desktop/cbcn_draft/file/linage/AllLineage.tsv", header = TRUE, sep = "\t") %>%
  select(cell = 1, CellFate = 3)

# combine annotation
asy_data <- asy_data %>%
  left_join(stage, by = "cell") %>%
  left_join(tissue, by = "cell") %>%
  select(cell, ce_mean, cb_mean, cn_mean, stage, CellFate)

# ----------------------------------------------------------
# 3. Load pv datasets (RMSD mean per cell)
# ----------------------------------------------------------
load("~/Desktop/cbcn/draft_code/new/RMSD_results.RData")

ce_RMSD_mean_cell <- ce_RMSD_mean_cell %>% rename(ce_mean = mean_RMSD)
cb_RMSD_mean_cell <- cb_RMSD_mean_cell %>% rename(cb_mean = mean_RMSD)
cn_RMSD_mean_cell <- cn_RMSD_mean_cell %>% rename(cn_mean = mean_RMSD)

pv_data <- ce_RMSD_mean_cell %>%
  left_join(cb_RMSD_mean_cell, by = "Cell") %>%
  left_join(cn_RMSD_mean_cell, by = "Cell") %>%
  rename(cell = Cell) %>%
  left_join(stage, by = "cell") %>%
  left_join(tissue, by = "cell")

# ----------------------------------------------------------
# 4. Merge pv and asynchrony
# ----------------------------------------------------------
merged_data <- inner_join(pv_data, asy_data, by = c("cell","stage","CellFate"),
                          suffix = c("_pv","_asy"))

# ----------------------------------------------------------
# 5. Long format
# ----------------------------------------------------------
data_long_pv <- merged_data %>%
  pivot_longer(cols = c(ce_mean_pv, cb_mean_pv, cn_mean_pv),
               names_to = "measure", values_to = "value_pv") %>%
  mutate(measure = recode(measure,
                          "ce_mean_pv"="Cel","cb_mean_pv"="Cbr","cn_mean_pv"="Cni"))

data_long_asy <- merged_data %>%
  pivot_longer(cols = c(ce_mean_asy, cb_mean_asy, cn_mean_asy),
               names_to = "measure", values_to = "value_asy") %>%
  mutate(measure = recode(measure,
                          "ce_mean_asy"="Cel","cb_mean_asy"="Cbr","cn_mean_asy"="Cni"))

final_data <- left_join(data_long_pv, data_long_asy,
                        by = c("cell","stage","CellFate","measure")) %>%
  mutate(measure = factor(measure, levels = c("Cel","Cbr","Cni")))

# ----------------------------------------------------------
# 6. Thresholds (optional, e.g. 85%)
# ----------------------------------------------------------
facet_thresholds <- final_data %>%
  group_by(measure, CellFate) %>%
  summarise(thresh_x = quantile(value_pv, 0.85, na.rm = TRUE),
            thresh_y = quantile(value_asy, 0.85, na.rm = TRUE),
            count_overlap = sum(value_pv >= quantile(value_pv,0.85,na.rm=TRUE) &
                                  value_asy >= quantile(value_asy,0.85,na.rm=TRUE),
                                na.rm=TRUE))


# ==========================================================
# Figure 5: pv (RMSD mean) vs asynchrony (all cells)
# ==========================================================

facet_thresholds_all <- final_data %>%
  group_by(measure) %>%
  summarise(
    thresh_x = quantile(value_pv, 0.85, na.rm = TRUE),
    thresh_y = quantile(value_asy, 0.85, na.rm = TRUE),
    count_overlap = sum(value_pv >= quantile(value_pv,0.85,na.rm=TRUE) &
                          value_asy >= quantile(value_asy,0.85,na.rm=TRUE),
                        na.rm=TRUE)
  )

p_all <- ggplot(final_data, aes(x = value_pv, y = value_asy, color = measure)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_abline(slope = 1, intercept = -2, color = "gray") +
  stat_cor(aes(group = 1), method = "pearson",
           label.x.npc = "left", label.y.npc = "top", size = 4) +
  geom_vline(data = facet_thresholds_all, aes(xintercept = thresh_x),
             color = "black", size = 1, alpha = 0.4) +
  geom_hline(data = facet_thresholds_all, aes(yintercept = thresh_y),
             color = "black", size = 1, alpha = 0.4) +
  geom_text(data = facet_thresholds_all,
            aes(x = thresh_x, y = thresh_y, label = paste0("n=", count_overlap)),
            vjust = -2, hjust = 0, color = "black", size = 5) +
  facet_wrap(~ measure, scales = "free") +
  labs(x="Value in pv data", y="Value in asynchrony",
       title="pv vs asynchrony", color="Group") +
  scale_x_continuous(limits=c(0.03,0.13)) +
  scale_y_continuous(limits=c(0.04,0.34)) +
  scale_color_manual(values=c("Cel"="#CD534C","Cbr"="#0073C2","Cni"="#EFC000")) +
  theme_minimal(base_size=14) +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=1),
        axis.ticks=element_line(color="black"))
print(p_all)
#ggsave("final_scatter_plot.pdf", p_all, width=10, height=3.5)


# ==========================================================
# Figure 2 (new_F2): pv vs asynchrony (specific tissue: Skin)
# ==========================================================

final_data_skin <- final_data %>% filter(CellFate == "Skin")

facet_thresholds_skin <- final_data_skin %>%
  group_by(measure) %>%
  summarise(
    thresh_x = quantile(value_pv, 0.8, na.rm = TRUE),
    thresh_y = quantile(value_asy, 0.8, na.rm = TRUE),
    count_overlap = sum(value_pv >= quantile(value_pv,0.8,na.rm=TRUE) &
                          value_asy >= quantile(value_asy,0.8,na.rm=TRUE),
                        na.rm=TRUE)
  )

p_skin <- ggplot(final_data_skin, aes(x=value_pv, y=value_asy, color=measure)) +
  geom_point(size=2) +
  stat_cor(aes(group=1), method="pearson",
           label.x.npc="left", label.y.npc="top", size=3) +
  geom_smooth(method="lm", se=FALSE,
              color=scales::alpha("black",0.8), size=0.8) +
  facet_wrap(~ measure, nrow=1, scales="free") +
  labs(x="Value in pv data", y="Value in asynchrony",
       title="pv vs asynchrony (Skin)", color="Group") +
  scale_color_manual(values=c("Cel"="#CD534C","Cbr"="#0073C2","Cni"="#EFC000")) +
  theme_minimal(base_size=14) +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=1),
        axis.ticks=element_line(color="black"),
        axis.text.x=element_text(size=8))
print(p_skin)
#ggsave("asychrony_pv_correlation_tissue_skin.pdf", p_skin, width=10, height=3.4)




# ==========================================================
# Figure 2 (new_F2): pv vs asynchrony (specific tissue: Skin)
# ==========================================================

p <- ggplot(final_data, aes(x = value_pv, y = value_asy)) +
  # 按 measure 设置颜色映射
  geom_point(aes(color = measure),size = 0.8) +
  # 添加相关系数
  stat_cor(aes(group = 1),
           method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top",
           size = 2) +
  # 添加回归线：solid 实线，透明黑色（40% 透明）
  geom_smooth(method = "lm", se = FALSE, 
              linetype = "solid", 
              color = scales::alpha("black", 0.8),
              size = 0.8) +
  # 使用 facet_grid 同时按 measure（行）与 Facet（列）分面
  facet_grid(measure ~ CellFate, scales = "free") +
  labs(x = "Value in data_pv", 
       y = "Value in data_asychrony",
       title = "pv VS asychrony",
       color = "Group") +
  scale_x_continuous(limits = c(0.03, 0.13)) +
  scale_y_continuous(limits = c(0.04, 0.34)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    # 设置 x 轴刻度字体大小为 8（可以根据需要调整）
    axis.text.x = element_text(size = 5)
  ) +
  # 根据 measure 赋予固定的颜色
  scale_color_manual(values = c("Cel" = "#CD534C", 
                                "Cbr" = "#0073C2", 
                                "Cni" = "#EFC000"))

# 显示图形
print(p)
#ggsave("asychrony_pv_correlation_tissue.pdf", p, width=10, height=4.5)






