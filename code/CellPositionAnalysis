


setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(reshape2)
library(RColorBrewer)

load("normalized_position.RData")
length_all_results <- read.table("~/Desktop/cbcn/draft_code/new/embryo_lengths.txt", header = TRUE, sep = "\t")

summary(ce_aligned)
head(ce_aligned[[1]], n = 5)
str(ce_aligned)

compute_distance_by_time <- function(df, embryo_name) {
  embryo_length_current <- length_all_results$embryo_length[
    length_all_results$embryo == embryo_name
  ]
  
  if (length(embryo_length_current) == 0) {
    warning(paste("没有找到胚胎", embryo_name, "对应的胚胎长度，使用1作为默认值"))
    embryo_length_current <- 1
  }
  
  time_points <- unique(df$con_time)
  
  distances_list <- list()
  
  for (t in time_points) {
    sub_df <- df %>% filter(con_time == t)
    
    coords <- as.matrix(sub_df %>% select(x, y, z))
    
    dmat <- as.matrix(dist(coords, method = "euclidean"))
    
    dmat <- dmat / embryo_length_current
    
    rownames(dmat) <- sub_df$cell
    colnames(dmat) <- sub_df$cell
    
    distances_list[[as.character(t)]] <- dmat
  }
  
  return(distances_list)
}


distance_ce <- lapply(names(ce_aligned), function(embryo_name) {
  compute_distance_by_time(ce_aligned[[embryo_name]], embryo_name)
})
names(distance_ce) <- names(ce_aligned)

distance_cb <- lapply(names(cb_aligned), function(embryo_name) {
  compute_distance_by_time(cb_aligned[[embryo_name]], embryo_name)
})
names(distance_cb) <- names(cb_aligned)

distance_cn <- lapply(names(cn_aligned), function(embryo_name) {
  compute_distance_by_time(cn_aligned[[embryo_name]], embryo_name)
})
names(distance_cn) <- names(cn_aligned)

summary(distance_ce)
names(distance_ce)
names(distance_ce$ce1)

calc_pairwise_RMSD <- function(distance_list, n_splits = 5) {
  dataset_names <- names(distance_list)
  all_results <- list()  
  result_idx <- 1
  
  
  pairs <- combn(dataset_names, 2, simplify = FALSE)
  
  for (pair in pairs) {
    ds1 <- pair[[1]]
    ds2 <- pair[[2]]
    
    
    timepoints1 <- names(distance_list[[ds1]])
    timepoints2 <- names(distance_list[[ds2]])
    common_tp <- intersect(timepoints1, timepoints2)
    if (length(common_tp) == 0) next 
    
    
    splits <- split(common_tp, cut(seq_along(common_tp), n_splits, labels = FALSE))
    
    for (chunk in splits) {
      for (tp in chunk) {
        m1 <- distance_list[[ds1]][[tp]]
        m2 <- distance_list[[ds2]][[tp]]
        
      
        if (is.null(rownames(m1)) || is.null(colnames(m1)) ||
            is.null(rownames(m2)) || is.null(colnames(m2))) {
          warning(paste("时间点", tp, "在", ds1, "或", ds2, "缺少行/列名。"))
          next
        }
        
        
        common_cells <- intersect(rownames(m1), rownames(m2))
        if (length(common_cells) == 0) next
        
        for (cell in common_cells) {
        
          v1 <- m1[cell, ]
          v1 <- v1[names(v1) != cell]
          v2 <- m2[cell, ]
          v2 <- v2[names(v2) != cell]
          
        
          common_vec_names <- intersect(names(v1), names(v2))
          if (length(common_vec_names) == 0) next
          v1 <- v1[common_vec_names]
          v2 <- v2[common_vec_names]
          
         
          rmsd_val <- sqrt(mean((v1 - v2)^2, na.rm = TRUE))
          
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
      
      cat("Completed chunk for", ds1, ds2, "with timepoints:", paste(chunk, collapse = ", "), "\n")
    }
  }
  
  if (length(all_results) > 0) {
    final_df <- do.call(rbind, all_results)
  } else {
    final_df <- data.frame()
  }
  
  return(final_df)
}


final_pairwise_ce <- calc_pairwise_RMSD(distance_ce, n_splits = 5)
final_pairwise_cb <- calc_pairwise_RMSD(distance_cb, n_splits = 5)
final_pairwise_cn <- calc_pairwise_RMSD(distance_cn, n_splits = 5)


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


data_long <- data %>%
  pivot_longer(
    cols = c(ce_mean_RMSD, cb_mean_RMSD, cn_mean_RMSD),
    names_to = "experiment",
    values_to = "RMSD"
  ) %>%
  mutate(experiment = factor(experiment,
                             levels = c("ce_mean_RMSD", "cb_mean_RMSD", "cn_mean_RMSD"),
                             labels = c("ce", "cb", "cn")))

p_4b <-ggplot(data_long, aes(x = experiment, y = RMSD, fill = experiment)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      
               color = "black",         
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
    axis.ticks = element_line(color = "black"),            
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0.025, 0.127), expand = c(0.001, 0.001)) +
  scale_x_discrete(expand = c(0.2, 0.2))       

print(p_4b)



library(reshape2)
data_long <- melt(data,
                  id.vars = c("cell",'stage'),
                  measure.vars = c("ce_mean_RMSD", "cb_mean_RMSD", "cn_mean_RMSD"))

data_long$variable <- factor(data_long$variable,
                             levels = c("cn_mean_RMSD", "cb_mean_RMSD","ce_mean_RMSD" ),
                             labels = c('Cni', "Cbr", "Cel"))

order_stage  <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                  "Stage_AB64", "Stage_AB128", "Stage_AB256")
order_letter <- c("A", "M", "E", "C", "D", "P")

data_long <- data_long %>%
  mutate(
    stage = factor(stage, levels = order_stage),
    first_letter = factor(substr(cell, 1, 1), levels = order_letter),
    cell_length = nchar(cell)
  ) %>%
  arrange(stage, first_letter, cell_length)

data_long$cell <- factor(data_long$cell, levels = unique(data_long$cell))

ggplot(data_long, aes(x = cell, y = variable, fill = value)) +
  geom_tile() +  
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  labs(x = "Cell", y = "Variable", title = "") +
  scale_fill_gradientn(
    colors = c("#2166AC", '#92C5DE', "#D1E5F0", "#FDDBC7", "#D6604D", "#B2182B"),
    values = c(0, 0.15, 0.3, 0.4, 0.75, 1)
  ) +
  theme_minimal() +  
  theme(
    axis.text.x  = element_blank(),  
    axis.text.y = element_text(size = 10),             
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )


ggplot(data_long, aes(x = cell, y = variable, fill = value)) +
  geom_tile() +  
  facet_grid(. ~ first_letter, scales = "free_x", space = "free") +
  labs(x = "Cell", y = "Variable", title = "") +
  scale_fill_gradientn(
    colors = c("#2166AC", '#92C5DE', "#D1E5F0", "#FDDBC7", "#D6604D", "#B2182B"),
    values = c(0, 0.15, 0.3, 0.4, 0.75, 1)
  ) +
  theme_minimal() + 
  theme(
    axis.text.x  = element_blank(),  
    axis.text.y = element_text(size = 10),              
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )





setwd("~/Desktop/cbcn/draft_code/new/")

library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(magrittr)
load("RMSD_results.RData")

colnames(ce_RMSD_mean_cell_tp)


bad_cells <- c("ABa",'ABp','EMS','P2')


ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_ce <- ce_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_ce)) %>%
  mutate(TimePoint = TimePoint - 9.54)

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cb <- cb_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cb)) %>%
  mutate(TimePoint = TimePoint + 1.66)

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cn <- cn_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cn)) %>%
  mutate(TimePoint = TimePoint + 2.28)

trim_cell <- function(Cell){
  last_cell <- str_sub(Cell, -1)
  remaining_cell <- str_sub(Cell, 1, -2)
  return(c(remaining_cell, last_cell))
}

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

klist <- list(Za = 2, Zpap = 0.3, Zppa = 0.5, Zpppa = 0.8, Zppp = .8, Zpppp = .3)
slist <- list(Za = -.5, Zp = -.2, Zpa = .22, Zpap = -.08, Zpp = -.04, Zppa = .04, Zppp = -.06, Zpppa = .01, Zpppp = -.03)
id.show.len <- 4
bg.color <- 'darkgray'

AllLineage.tsv <- read.csv("~/Desktop/cbcn/CellID.csv") %>% tibble()

n_death <- sum(AllLineage.tsv$CellFate == "Death", na.rm = TRUE)
print(n_death)

process_lineage_and_plot <- function(lineage_data, std_column, time_limit, 
                                     intercept, slope, normfactor, data_type = "cb", expdata) {
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
  
  p <- p.data.y %>%
    mutate(label = if_else(nchar(ID) < id.show.len, CellName, NA)) %>%
    ggplot() + 
    geom_path(aes(x = x, y = -Posi, group = ID), color = bg.color) +
    geom_text(aes(x = x, y = -Mid, label = label)) +
    geom_path(data = p.data.exp, aes(x = x, y = -TP, group = G, color = EXP)) +
    geom_path(data = p.data.x, aes(x = Posi, y = -Start, group = ID), color = bg.color) +
    scale_x_continuous(n.breaks = 20) +
    scale_y_continuous(limits = c(-200, 20), labels = function(x) abs(x)) +
    scale_color_gradientn(
      colors = c("blue", "yellow","yellow",'#f55042',"#B21F1F"),
      values = c(0, 0.15, 0.4,0.81, 1),
      rescaler = function(x, to = c(0,1), from = NULL) {
        x_clamped <- ifelse(x < 0.023, 0.023, ifelse(x > 0.1, 0.1, x))
        scales::rescale(x_clamped, to = to, from = c(0.023, 0.1))
      },
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

p_ce <- process_lineage_and_plot(
  read.csv("~/Desktop/cbcn/CDFile/ce/CD200113plc1p2.csv"),
  "std_len_ce", 200, 7.04, 1.04, 9.54,data_type = "ce",ce_RMSD_mean_cell_tp
)
print(p_ce)


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

bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_ce <- ce_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_ce)) %>%
  mutate(TimePoint = TimePoint - 22.04)

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cb <- cb_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cb)) %>%
  mutate(TimePoint = TimePoint - 13.29)

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cn <- cn_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cn)) %>%
  mutate(TimePoint = TimePoint - 9.1)

max_cn_time <- max(cn_RMSD_mean_cell_tp$TimePoint, na.rm = TRUE)

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

sorted_tp <- sort(as.numeric(names(table(cb_RMSD_mean_cell_tp$TimePoint))))
print(sorted_tp)


ce_RMSD_mean_cell_tp$Group <- "ce"
cb_RMSD_mean_cell_tp$Group <- "cb"
cn_RMSD_mean_cell_tp$Group <- "cn"

all_data <- bind_rows(ce_RMSD_mean_cell_tp, cb_RMSD_mean_cell_tp, cn_RMSD_mean_cell_tp)

library(dplyr)

result <- all_data %>%
  filter(Group == "ce") %>%                  
  group_by(TimePoint) %>%                     
  summarise(nCells = n_distinct(Cell)) %>%    
  filter(nCells == 340)                       
##28cell, 170 cell, 300cell, 509cell
print(result)
cetimepoints <- result$TimePoint
cbtimepoints <- result$TimePoint
cntimepoints <- result$TimePoint
print(cetimepoints)
print(cbtimepoints)
print(cntimepoints)

min_time <- min(all_data$TimePoint) - 1
max_time <- max(all_data$TimePoint) + 1

p_4c <- ggplot(all_data, aes(x = TimePoint, y = mean_RMSD, color = Group, fill = Group)) +
  geom_boxplot(aes(group = interaction(TimePoint, Group)),
               outlier.shape = NA,
               width = 1,              
               position = position_dodge(width = 0),
               alpha = 0,                
               coef = 0,                
               size = 0.2) +              
  stat_boxplot(geom = "errorbar",
               width = 0.2,
               aes(group = interaction(TimePoint, Group)),
               position = position_dodge(width = 0),
               coef = 1.3,    
               alpha = 0.4,                 
               linewidth = 0.3) +
  geom_vline(xintercept = 20.23, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 68.3, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 98.11, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 143.31, linetype = "dashed", color = "black", size = 0.5) +
  labs(x = "TimePoint", y = "RMSD", color = "Group") +
  scale_y_continuous(limits = c(0.02, 0.13), expand = c(0, 0)) +
  scale_x_continuous(limits = c(min_time, max_time),
                     expand = c(0, 0)) +
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

bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_ce <- ce_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_ce)) %>%
  mutate(TimePoint = TimePoint - 22.04)

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cb <- cb_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cb)) %>%
  mutate(TimePoint = TimePoint - 13.29)

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cn <- cn_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cn)) %>%
  mutate(TimePoint = TimePoint - 9.1)

max_cn_time <- max(cn_RMSD_mean_cell_tp$TimePoint, na.rm = TRUE)

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

sorted_tp <- sort(as.numeric(names(table(cb_RMSD_mean_cell_tp$TimePoint))))
print(sorted_tp)


ce_RMSD_mean_cell_tp$Group <- "ce"
cb_RMSD_mean_cell_tp$Group <- "cb"
cn_RMSD_mean_cell_tp$Group <- "cn"

all_data <- bind_rows(ce_RMSD_mean_cell_tp, cb_RMSD_mean_cell_tp, cn_RMSD_mean_cell_tp)
all_data <- all_data %>%
  mutate(Facet = if_else(substr(Cell, 1, 1) == "Z", "P", substr(Cell, 1, 1)))%>%
  filter(Facet != "P")%>%
  mutate(Facet = factor(Facet, levels = c("A", "E", "M", "C", "D")))

min_time <- min(all_data$TimePoint) - 1
max_time <- max(all_data$TimePoint) + 1

p_4c <- ggplot(all_data, aes(x = TimePoint, y = mean_RMSD, color = Group, fill = Group)) +
  geom_boxplot(aes(group = interaction(TimePoint, Group)),
               outlier.shape = NA,
               width = 1,             
               position = position_dodge(width = 0),
               alpha = 0,             
               coef = 0,               
               size = 0.2) +
  stat_boxplot(geom = "errorbar",
               width = 0.2,
               aes(group = interaction(TimePoint, Group)),
               position = position_dodge(width = 0),
               coef = 1.3,
               alpha = 0.4,
               linewidth = 0.3) +
  geom_vline(xintercept = 20.23, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 68.3, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 98.11, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 143.31, linetype = "dashed", color = "black", size = 0.5) +
  labs(x = "TimePoint", y = "RMSD", color = "Group") +
  scale_y_continuous(limits = c(0.02, 0.13), expand = c(0, 0)) +
  scale_x_continuous(limits = c(min_time, max_time), expand = c(0, 0)) +
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
    strip.background = element_blank()  
  ) +
  facet_wrap(~ Facet, nrow = 1)

print(p_4c)
#ggsave(filename = "pv_box_lineage.pdf", plot = p_4c, device = "pdf", width = 16, height = 3.8)



setwd("~/Desktop/cbcn/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggpubr)

load("RMSD_results.RData")
setwd("~/Desktop/cbcn/draft_code/fig4_all/")

colnames(ce_RMSD_mean_cell)

bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")

ce_RMSD_mean_cell$Group <- "ce"
cb_RMSD_mean_cell$Group <- "cb"
cn_RMSD_mean_cell$Group <- "cn"

all_data <- bind_rows(ce_RMSD_mean_cell, cb_RMSD_mean_cell, cn_RMSD_mean_cell)
colnames(all_data)

all_data <- all_data %>%
  mutate(Group = factor(Group, levels = c("ce", "cb", "cn")))

colnames(all_data)


lr <- read.csv("~/Desktop/cbcn/LR_symmetry_2.csv")
colnames(lr)

lr_long <- lr %>%
  mutate(pair = row_number()) %>%
  pivot_longer(cols = c("L", "R"), names_to = "side", values_to = "Cell") %>%
  mutate(lr = if_else(side == "L", "l", "r")) %>%
  select(-side)

all_data <- all_data %>%
  left_join(lr_long, by = "Cell") %>%
  mutate(lr = if_else(is.na(lr), "n", lr))

colnames(all_data)

data_pair <- all_data %>%
  filter(lr %in% c("l", "r")) %>% 
  group_by(Group, pair, lr) 

data_pair_wide <- data_pair %>%
  pivot_wider(names_from = lr, values_from = c(mean_RMSD, Cell)) %>%
  filter(!is.na(mean_RMSD_l) & !is.na(mean_RMSD_r)) %>%
  mutate(diff = mean_RMSD_l - mean_RMSD_r)



p_pair_box <- ggplot(data_pair, aes(x = lr, y = mean_RMSD)) +
  stat_boxplot(aes(color = lr), geom = "errorbar", width = 0.5,size = 0.8) +
  geom_point(aes(group = factor(pair)), 
             color = "grey", 
             size = 1.2, 
             alpha = 0.4) +
  geom_line(aes(group = factor(pair)), color = "gray", alpha = 0.15) +
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
        strip.background = element_blank(),  
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  stat_compare_means(comparisons = list(c("l", "r")),
                     method = "t.test",     
                     label = "p.signif",     
                     hide.ns = TRUE,        
                     label.y = max(data_pair$mean_RMSD, na.rm = TRUE) * 0.9)

print(p_pair_box)



p_diff_individual <- ggplot(data_pair_wide, aes(x = factor(pair), y = diff, fill = Group)) +
  geom_col() +
  facet_wrap(~ Group, scales = "free_x") +
  labs(title = "Difference (l - r) for each pair",
       x = "Pair ID",
       y = "Difference (l - r)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),  
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


p_diff_overall <- ggplot(data_pair_wide, aes(x = Group, y = abs(diff), fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.42, linewidth = 0.5) +  
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,     
               color = "black",         
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
    axis.ticks = element_line(color = "black"),           
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
    method = "t.test",
    label = "p.signif",   
    hide.ns = TRUE,
    label.y = max(abs(data_pair_wide$diff), na.rm = TRUE) * 0.57
  )+
  scale_y_continuous(limits = c(0, 0.045), expand = c(0.0005, 0.001))
print(p_diff_overall)

ggsave("pv_lr_pair_diff_box.pdf", plot = p_diff_overall, width = 10, height = 10)

combined_diff <- (p_pair | p_diff_individual) / p_diff_overall
print(combined_diff)



setwd("~/Desktop/cbcn/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)

load("RMSD_results.RData")
setwd("~/Desktop/cbcn/draft_code/fig4_all/")

colnames(ce_RMSD_mean_cell_tp)

bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_ce <- ce_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_ce)) %>%
  mutate(TimePoint = TimePoint - 22.04)

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cb <- cb_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cb)) %>%
  mutate(TimePoint = TimePoint - 13.29)

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  mutate(TimePoint = as.numeric(as.character(TimePoint)))

bad_timepoints_cn <- cn_RMSD_mean_cell_tp %>%
  filter(Cell %in% bad_cells) %>%
  pull(TimePoint) %>%
  unique()

cn_RMSD_mean_cell_tp <- cn_RMSD_mean_cell_tp %>%
  filter(!(TimePoint %in% bad_timepoints_cn)) %>%
  mutate(TimePoint = TimePoint - 9.1)

max_cn_time <- max(cn_RMSD_mean_cell_tp$TimePoint, na.rm = TRUE)

ce_RMSD_mean_cell_tp <- ce_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

cb_RMSD_mean_cell_tp <- cb_RMSD_mean_cell_tp %>%
  filter(TimePoint <= max_cn_time)

sorted_tp <- sort(as.numeric(names(table(cb_RMSD_mean_cell_tp$TimePoint))))
print(sorted_tp)


ce_RMSD_mean_cell_tp$Group <- "ce"
cb_RMSD_mean_cell_tp$Group <- "cb"
cn_RMSD_mean_cell_tp$Group <- "cn"

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

min_time <- min(all_data$TimePoint) - 1
max_time <- max(all_data$TimePoint) + 1

p_4c <- ggplot(all_data, aes(x = TimePoint, y = mean_RMSD, color = Group, fill = Group)) +
  geom_boxplot(aes(group = interaction(TimePoint, Group)),
               outlier.shape = NA,
               width = 1,            
               position = position_dodge(width = 0),
               alpha = 0,           
               coef = 0,              
               size = 0.2) +
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
  facet_wrap(~ factor(lr,levels = c('l','r','n')))

print(p_4c)





setwd("~/Desktop/cbcn/draft_code/new/")
getwd()

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggpubr)

load("RMSD_results.RData")

colnames(ce_RMSD_mean_cell)

bad_cells <- c('ABa','ABp','EMS','P2',"ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")


ce_RMSD_mean_cell$Group <- "ce"
cb_RMSD_mean_cell$Group <- "cb"
cn_RMSD_mean_cell$Group <- "cn"

all_data <- bind_rows(ce_RMSD_mean_cell, cb_RMSD_mean_cell, cn_RMSD_mean_cell)
colnames(all_data)

all_data <- all_data %>%
  mutate(Group = factor(Group, levels = c("ce", "cb", "cn")))

colnames(all_data)


lr <- read.csv("~/Desktop/cbcn/LR_symmetry_2.csv")
colnames(lr)

lr_long <- lr %>%
  mutate(pair = row_number()) %>%
  pivot_longer(cols = c("L", "R"), names_to = "side", values_to = "Cell") %>%
  mutate(lr = if_else(side == "L", "l", "r")) %>%
  select(-side)

all_data <- all_data %>%
  left_join(lr_long, by = "Cell") %>%
  mutate(lr = if_else(is.na(lr), "n", lr))

colnames(all_data)

data_pair <- all_data %>%
  filter(lr %in% c("l", "r")) %>% 
  group_by(Group, pair, lr) 

data_pair_wide <- data_pair %>%
  pivot_wider(names_from = lr, values_from = c(mean_RMSD, Cell)) %>%
  filter(!is.na(mean_RMSD_l) & !is.na(mean_RMSD_r)) %>%
  mutate(diff = mean_RMSD_l - mean_RMSD_r)


p_pair_box <- ggplot(data_pair, aes(x = lr, y = mean_RMSD)) +
  stat_boxplot(aes(color = lr), geom = "errorbar", width = 0.5,size = 0.8) +
  geom_point(aes(group = factor(pair)), 
             color = "grey", 
             size = 1.2, 
             alpha = 0.4) +
  geom_line(aes(group = factor(pair)), color = "gray", alpha = 0.15) +
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
        strip.background = element_blank(),  
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  stat_compare_means(comparisons = list(c("l", "r")),
                     method = "t.test",     
                     label = "p.signif",    
                     hide.ns = TRUE,        
                     label.y = max(data_pair$mean_RMSD, na.rm = TRUE) * 0.9)

print(p_pair_box)




p_diff_individual <- ggplot(data_pair_wide, aes(x = factor(pair), y = diff, fill = Group)) +
  geom_col() +
  facet_wrap(~ Group, scales = "free_x") +
  labs(title = "Difference (l - r) for each pair",
       x = "Pair ID",
       y = "Difference (l - r)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        strip.text = element_text(size = 12))
print(p_diff_individual)
#ggsave("pv_lr_pair_value.pdf", plot = p_diff_individual, width = 10, height = 6)




p_diff_overall <- ggplot(data_pair_wide, aes(x = Group, y = abs(diff), fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.42, linewidth = 0.5) + 
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      
               color = "black",      
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
    axis.ticks = element_line(color = "black"),         
    axis.text = element_text(size = 16),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
    method = "t.test",
    label = "p.signif", 
    hide.ns = TRUE,
    label.y = max(abs(data_pair_wide$diff), na.rm = TRUE) * 0.57
  )+
  scale_y_continuous(limits = c(0, 0.045), expand = c(0.0005, 0.001))
print(p_diff_overall)




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
  geom_point(aes(color = measure),size = 0.8) +
  stat_cor(aes(group = 1),
           method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top",
           size = 2) +
  geom_smooth(method = "lm", se = FALSE, 
              linetype = "solid", 
              color = scales::alpha("black", 0.8),
              size = 0.8) +
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
    axis.text.x = element_text(size = 5)
  ) +
  scale_color_manual(values = c("Cel" = "#CD534C", 
                                "Cbr" = "#0073C2", 
                                "Cni" = "#EFC000"))

print(p)
#ggsave("asychrony_pv_correlation_tissue.pdf", p, width=10, height=4.5)






