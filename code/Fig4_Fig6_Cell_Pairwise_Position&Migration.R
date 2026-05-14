setwd("~/Desktop/cbcn/draft_code/new/")
getwd()
library(ggplot2);library(dplyr);library(tidyverse);library(ggsci);library(ggpubr);library(reshape2);library(RColorBrewer)

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

data_list <- list(
  ce1 = ce1, ce2 = ce2, ce3 = ce3, ce4 = ce4, ce5 = ce5, ce6 = ce6, ce7 = ce7, ce8 = ce8,
  cb1 = cb1, cb2 = cb2, cb3 = cb3, cb4 = cb4, cb5 = cb5, cb6 = cb6,
  cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cn5 = cn5, cn6 = cn6, cn7 = cn7
)

library(dplyr)

common_cells <- Reduce(intersect, lapply(data_list, function(df) unique(df$cell)))

data_list <- lapply(data_list, function(df) {
  df %>% 
    filter(cell %in% common_cells)
})

data_list <- lapply(data_list, function(df) {
  df %>% group_by(cell) %>% 
    mutate(time_percent = (time - min(time)) / (max(time) - min(time))) %>% 
    ungroup()
})

align_cell <- function(cell_name, data_list) {
  cell_data <- lapply(data_list, function(df) {
    df %>% 
      filter(cell == cell_name) %>% 
      arrange(time_percent)
  })
  
  n_rows <- sapply(cell_data, nrow)
  template_index <- which.min(n_rows)
  template_data <- cell_data[[template_index]]
  template_time <- template_data$time_percent 
  
  aligned_data <- lapply(cell_data, function(df) {
    indices <- sapply(template_time, function(t_val) {
      which.min(abs(df$time_percent - t_val))
    })
    indices <- unlist(indices)
    df[indices, ]
  })
  return(aligned_data)
}


aligned_by_cell <- lapply(common_cells, function(cell) {
  align_cell(cell, data_list)
})
names(aligned_by_cell) <- common_cells


library(purrr)
aligned_data_list <- vector("list", length(data_list))
names(aligned_data_list) <- names(data_list)

for (i in seq_along(data_list)) {
  aligned_data_list[[i]] <- map_dfr(aligned_by_cell, ~ .x[[i]])
}


str(aligned_data_list[[1]])
summary(aligned_data_list)


library(dplyr)
library(purrr)

emblength <- read.table("~/Desktop/cbcn/embryo_lengths.txt", header = TRUE, sep = "\t")

calculate_distance <- function(df, embryo_id) {
  df <- df %>%
    mutate(embryo = embryo_id)
  
  df_distance <- df %>%
    group_by(cell) %>%                           
    arrange(time, .by_group = TRUE) %>%            
    mutate(
      distance = sqrt((lead(x) - x)^2 + (lead(y) - y)^2 + (lead(z) - z)^2),
      distance = if_else(row_number() == 1, NA_real_, distance)
    ) %>%
    summarise(total_distance = sum(distance, na.rm = TRUE), .groups = "drop")
  
  embryo_length_value <- emblength$embryo_length[which(emblength$embryo == embryo_id)]
  
  df_distance <- df_distance %>%
    mutate(normalized_distance = total_distance / embryo_length_value)
  
  return(df_distance)
}

results_list <- mapply(function(embryo_id, df) {
  calculate_distance(df, embryo_id)
}, names(aligned_data_list), aligned_data_list, SIMPLIFY = FALSE)

norm_list <- lapply(names(results_list), function(embryo_id) {
  results_list[[embryo_id]] %>% 
    select(cell, normalized_distance) %>% 
    rename(!!paste0(embryo_id, "_norm") := normalized_distance)
})

combined_norm <- reduce(norm_list, full_join, by = "cell")

print(colnames(combined_norm))

write.csv(combined_norm, file = "con_cellMovement.csv", row.names = FALSE)



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
    ce_mean = mean(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_sd   = sd(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_cv   = ce_sd / ce_mean,         
    
    cb_mean = mean(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_sd   = sd(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_cv   = cb_sd / cb_mean,         
    cn_mean = mean(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_sd   = sd(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_cv   = cn_sd / cn_mean         
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
data_long$measure <- factor(data_long$measure,
                            levels = c("ce_cv", "cb_cv", "cn_cv"),
                            labels = c("ce", "cb", "cn"))

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
library(ggpubr)  

ggplot(data_long, aes(x = measure, y = RMSD, fill = measure)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +  
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,    
               color = "black",     
               size = 0.5) +
  #facet_grid(. ~ first_letter, scales = "free_x", space = "free") +
  #facet_grid(. ~ CellFate, scales = "free_x", space = "free") +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  #scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(title = "",
       x = "movement_cv",
       y = "Value") +
  scale_fill_manual(values = c("ce" = "#CD534C",  
                               "cb" = "#0073C2",   
                               "cn" = "#EFC000")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),            
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_compare_means(comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
                     label = "p.signif",   
                     method = "wilcox.test",
                     label.y = c(0.28, 0.28, 0.28)) 



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
  select(1, 44)
data <- merge(data, stage, by = "cell")
colnames(data)

data <- data %>%
  mutate(ce_movement_mean = round(rowMeans(select(., ce1_norm:ce8_norm), na.rm = TRUE), 2)) %>%
  mutate(cb_movement_mean = round(rowMeans(select(., cb1_norm:cb6_norm), na.rm = TRUE), 2)) %>%
  mutate(cn_movement_mean = round(rowMeans(select(., cn1_norm:cn7_norm), na.rm = TRUE), 2)) %>%
  select(1, 24:26, 23) %>%
  na.omit()
colnames(data)

data_long <- data %>%
  pivot_longer(
    cols = c("ce_movement_mean", "cb_movement_mean", "cn_movement_mean"),
    names_to = "experiment",
    values_to = "movement_mean"
  ) %>%
  mutate(
    experiment = gsub("_movement_mean", "", experiment),
    experiment = factor(experiment, levels = c("ce", "cb", "cn"))
  )

p_4d <-ggplot(data_long, aes(x = experiment, y = movement_mean, fill = experiment)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) + 
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,    
               color = "black",        
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
    axis.ticks = element_line(color = "black"),        
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0.2, 1.0), expand = c(0.01, 0.01)) +
  scale_x_discrete(expand = c(0.2, 0.2))       

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

library(dplyr)
library(tidyr)
library(ggplot2)

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

gastrulating3 <- c("ABarpapa", "ABarpapp", "ABalaapp", "ABalappa", 
                   "ABalaapa", "ABarapppp", "ABarapppa", "ABpraappp", "ABpraappa")

data_split1 <- data %>%
  mutate(final_group = if_else(cell %in% gas1,
                               "1_gastrulating cells",
                               "1_other cells"),
         split = "split1") 

data_split2 <- data %>%
  mutate(final_group = if_else(cell %in% gastrulating2,
                               "2_gastrulating cells",
                               "2_other cells"),
         split = "split2") 
data_split3 <- data %>%
  mutate(final_group = if_else(cell %in% gastrulating3,
                               "3_gastrulating cells",
                               "3_other cells"),
         split = "split3")  

combined_data <- bind_rows(data_split1, data_split2, data_split3)

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

data_long_split1 <- data_long %>% filter(split == "split1")
data_long_split2 <- data_long %>% filter(split == "split2")
data_long_split3 <- data_long %>% filter(split == "split3")

base_boxplot <- function(data, title_text) {
  dodge_val <- position_dodge(width = 0.75)
  
  ggplot(data, aes(x = final_group, y = RMSD, fill = measure)) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.42,
      size = 0.5,
      position = dodge_val  
    ) +
    geom_boxplot(
      width = 0.65, 
      outlier.shape = NA,     
      color = "black",      
      size = 0.5,
      position = dodge_val    
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
      axis.ticks = element_line(color = "black"),  
      axis.text = element_text(size = 12),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0.2, 1.0), expand = c(0, 0))
}
p1 <- base_boxplot(data_long_split1, "Group 1")
p2 <- base_boxplot(data_long_split2, "Group 2")
p3 <- base_boxplot(data_long_split3, "Group 3")

print(p1)
print(p2)
print(p3)

library(ggplot2)
library(gridExtra)
library(grid)


combined_plot <- arrangeGrob(p1, p2, p3, nrow = 1)

grid.newpage()        
grid.draw(combined_plot)


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
  select(1, 44)
data <- merge(data, stage, by = "cell")
colnames(data)

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


library(pheatmap)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggplotify)  
library(cowplot)    
library(patchwork)

ann_colors <- list(
  group = c(A = "#F8766D", 
            M = "#A3A500",
            E = "#00BF7D", 
            C = "#00B0F6",
            D = "#E76BF3")
)

gas1 <- c("Ea", "Ep", "ABalapp", "ABalaap", "ABalpppa", "ABalpppp", "ABarappp", "ABprpaap", "ABplppap", "MSaaaa", "MSaaap", "Z2", "Z3", "MSapaa", "MSapap", "MSappa", "MSappp", "MSpaaa", "MSpaap", "MSppaa", "MSppap", "MSpppa", "MSpppp", "Daa", "Dap", "Dpa", "Dpp", "ABalpaapp", "ABaraappa", "ABaraappp", "ABaraapaa", "ABaraapap", "ABarapaap", "Cpppa", "Cpppp", "ABalpaaap", "MSaapaa", "MSaapap", "MSpapaa", "MSpapap", "Cappa", "Cappp", "Cppaa", "Cppap", "MSaappa", "MSaappp", "Capaa", "Capap", "MSpappa", "MSpappp", "ABaraaapaa", "ABaraaapap", "ABalpaapaa", "ABalpaapap", "ABalpappaa", "ABarapappa", "ABarapappp", "ABalpapppa", "ABalpapppp", "ABarapapaa", "ABprpapppa", "ABprpapppp", 
          "ABarapaaaa", "ABarapaaap", "ABalpaaaaa", "ABalpaaaap")

gastrulating2 <- c("Ep", "Ea", "ABalaap", "ABalpppp", "ABalpppa", "ABarappp", "ABplppap", "MSpaaa", "MSaaap", "ABprppap", "MSpppp", "MSapap", "Z2", "MSapaa", "MSppap", "MSappp", "MSpaap", "Z3", "MSppaa", "MSappa", "MSaaaa", "MSpppa", "Dpa", "Daa", "Dpp","ABaraappp", "ABaraapap", "Dap", "ABaraapaa", "ABaraaaap", "Cappa", "ABaraappa", "MSaapap", "MSaapaa", "ABarapaap", "Cppap", "Capaa", "MSpapap", "Cappp", "Cppaa", "MSpapaa", "Cpppp", "Capap", "Cpppa", "ABalpaaap", "ABprpapaa", "MSaappp", "MSpappp", "MSpappa", "MSaappa")

gastrulating3 <- c("ABarpapa", "ABarpapp", "ABalaapp", "ABalappa", "ABalaapa", "ABarapppp", "ABarapppa", "ABpraappp", "ABpraappa", "ABalpppap", "ABarppppp", "ABarppppa", "ABalpppaa", "ABalppppp","ABalapppa", "ABalappap", "ABalppppa", "ABalapppp")

gas1_first    <- gas1[1:(length(gas1)/2)]
gas1_second   <- gas1[(length(gas1)/2 + 1):length(gas1)]

gas2_first    <- gastrulating2[1:(length(gastrulating2)/2)]
gas2_second   <- gastrulating2[(length(gastrulating2)/2 + 1):length(gastrulating2)]

gas3_first    <- gastrulating3[1:(length(gastrulating3)/2)]
gas3_second   <- gastrulating3[(length(gastrulating3)/2 + 1):length(gastrulating3)]

data_split1 <- data %>%
  mutate(final_group = if_else(cell %in% gas1,
                               "1_gastrulating cells",
                               "1_other cells"),
         split = "split1")  

data_split2 <- data %>%
  mutate(final_group = if_else(cell %in% gastrulating2,
                               "2_gastrulating cells",
                               "2_other cells"),
         split = "split2") 

data_split3 <- data %>%
  mutate(final_group = if_else(cell %in% gastrulating3,
                               "3_gastrulating cells",
                               "3_other cells"),
         split = "split3")  

combined_data <- bind_rows(data_split1, data_split2, data_split3)

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

data_long_split1 <- data_long %>% filter(split == "split1")
data_long_split2 <- data_long %>% filter(split == "split2")
data_long_split3 <- data_long %>% filter(split == "split3")

data_long_split1 <- data_long_split1 %>% mutate(cell = factor(cell, levels = gas1))
data_long_split2 <- data_long_split2 %>% mutate(cell = factor(cell, levels = gastrulating2))
data_long_split3 <- data_long_split3 %>% mutate(cell = factor(cell, levels = gastrulating3))

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

base_boxplot <- function(data, title_text) {
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
    stat_compare_means(comparisons = list(c("diff_cb", "diff_cn")),
                       label = "p.signif", method = "wilcox.test", size = 5)
  
}

p1_first   <- base_boxplot(data_long_split1_first, "Group 1 - First Half")
p1_second  <- base_boxplot(data_long_split1_second, "Group 1 - Second Half")
p2_first   <- base_boxplot(data_long_split2_first, "Group 2 - First Half")
p2_second  <- base_boxplot(data_long_split2_second, "Group 2 - Second Half")
p3_first   <- base_boxplot(data_long_split3_first, "Group 3 - First Half")
p3_second  <- base_boxplot(data_long_split3_second, "Group 3 - Second Half")

combined_plot <- (p1_first | p1_second) / (p2_first | p2_second) / (p3_first | p3_second)

print(combined_plot)
#ggsave("movement_gas_group.pdf", combined_plot, width = 10, height = 15)








library(ggplot2)
library(dplyr)
library(ggpubr)

base_barplot <- function(data, title_text) {
  data <- data %>% filter(grepl("gastrulating cells", final_group))
  
  dodge <- position_dodge(width = 0.65)
  
  ggplot(data, aes(x = cell, y = distance, fill = measure)) +
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

p1 <- base_barplot(data_long_split1, "Group 1")
p2 <- base_barplot(data_long_split2, "Group 2")
p3 <- base_barplot(data_long_split3, "Group 3")

print(p1)
print(p2)
print(p3)

library(patchwork)

combined_plot <- p1 / p2 / p3

print(combined_plot)


ggsave(filename = "movement_diff_gasonly_cell.pdf",
       plot = combined_plot,
       width = 16,  
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
    ce_mean = mean(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_sd   = sd(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_cv   = ce_sd / ce_mean,        
    
    cb_mean = mean(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_sd   = sd(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_cv   = cb_sd / cb_mean,        
    
    cn_mean = mean(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_sd   = sd(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_cv   = cn_sd / cn_mean          
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
data_long$measure <- factor(data_long$measure,
                            levels = c("ce_cv", "cb_cv", "cn_cv"),
                            labels = c("ce", "cb", "cn"))

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

ggplot(data_long, aes(x = cell, y = measure, fill = RMSD)) +
  geom_tile(color = "white") +
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  # scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
  #                      midpoint = 0.17, name = "Z-score") +
  scale_fill_gradientn(
    colors = c("#2166AC", '#92C5DE', "#D1E5F0", "#FDDBC7", "#D6604D", "#B2182B"),
    values = c(0, 0.15, 0.35, 0.5, 0.75, 1)
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
    ce_mean = mean(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_sd   = sd(c_across(matches("^ce\\d+_norm$")), na.rm = TRUE),
    ce_cv   = ce_sd / ce_mean,     
    
    cb_mean = mean(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_sd   = sd(c_across(matches("^cb\\d+_norm$")), na.rm = TRUE),
    cb_cv   = cb_sd / cb_mean,        
    
    cn_mean = mean(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_sd   = sd(c_across(matches("^cn\\d+_norm$")), na.rm = TRUE),
    cn_cv   = cn_sd / cn_mean        
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
data_long$measure <- factor(data_long$measure,
                            levels = c("ce_cv", "cb_cv", "cn_cv"),
                            labels = c("ce", "cb", "cn"))

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
library(ggpubr) 

ggplot(data_long, aes(x = measure, y = RMSD, fill = measure)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,    
               color = "black",     
               size = 0.5) +
  facet_grid(. ~ first_letter, scales = "free_x", space = "free") +
  #facet_grid(. ~ CellFate, scales = "free_x", space = "free") +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  #scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(title = "",
       x = "movement_cv",
       y = "Value") +
  scale_fill_manual(values = c("ce" = "#CD534C",   
                               "cb" = "#0073C2",  
                               "cn" = "#EFC000")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),        
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_compare_means(comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
                     label = "p.signif",   
                     method = "wilcox.test",
                     label.y = c(0.28, 0.28, 0.28)) 

ggplot(data_long, aes(x = measure, y = RMSD, fill = measure)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) + 
  geom_boxplot(width = 0.65, 
               outlier.shape = NA,      
               color = "black",        
               size = 0.5) +
  facet_grid(. ~ CellFate, scales = "free_x", space = "free") +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0, 0)) +
  #scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(title = "",
       x = "movement_cv",
       y = "Value") +
  scale_fill_manual(values = c("ce" = "#CD534C", 
                               "cb" = "#0073C2",   
                               "cn" = "#EFC000")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),         
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  stat_compare_means(comparisons = list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn")),
                     label = "p.signif",  
                     method = "wilcox.test",
                     label.y = c(0.28, 0.28, 0.28)) 

