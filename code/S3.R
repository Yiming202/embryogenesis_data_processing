# ============================================================
# Setup
# ============================================================
setwd("~/Desktop/cbcn/draft_code/new/")

library(readr)
library(dplyr)
library(tidyr)
library(ape)
library(igraph)
library(ggtree)
library(ggplot2)
library(patchwork)

# ============================================================
# Helper: compute per-parent asymmetry for a given group
# ============================================================
process_group <- function(data, lineage_data_top, new_prefix) {
  pattern <- switch(
    new_prefix,
    ce = "corTim_cep",
    cb = "corTim_cbp",
    cn = "corTim_cnp"
  )
  
  data <- data %>%
    rename_with(
      ~ sub(pattern, paste0("time_", new_prefix), .x),
      starts_with(pattern)
    )
  
  data_processed <- data %>%
    mutate(parent = substr(cell, 1, nchar(cell) - 1)) %>%
    filter(
      nchar(cell) > 1,
      !cell %in% c("M", "P2", "P3", "P4", "Z2", "Z3", "EMS", "MS")
    ) %>%
    rename(child = cell) %>%
    bind_rows(lineage_data_top) %>%
    filter(!(parent %in% c("P0", "P1", "P4")))
  
  time_cols <- grep(paste0("^time_", new_prefix), names(data_processed), value = TRUE)
  
  data_processed %>%
    group_by(parent) %>%
    summarise(
      across(
        all_of(time_cols),
        ~ {
          last_letters <- substr(child, nchar(child), nchar(child))
          pos_val <- .x[last_letters %in% c("p", "r", "v", "3", "4")]
          neg_val <- .x[last_letters %in% c("a", "l", "d", "C", "D")]
          if (length(pos_val) == 1 && length(neg_val) == 1) pos_val - neg_val else NA_real_
        },
        .names = "diff_{sub('time_', '', .col)}"
      ),
      .groups = "drop"
    ) %>%
    filter(!if_any(starts_with(paste0("diff_", new_prefix)), is.na))
}

# ============================================================
# CE group processing
# ============================================================
data_ce <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv") %>%
  select(cell, 23:30)

lineage_data_top_ce <- data.frame(
  parent = c("P0", "P0", "P1", "P1", "EMS", "EMS", "P2", "P2", "P3", "P3", "P4", "P4"),
  child  = c("AB", "P1", "EMS", "P2", "MS", "E", "P3", "C", "P4", "D", "Z2", "Z3"),
  time_ce1 = 1, time_ce2 = 1, time_ce3 = 1, time_ce4 = 1,
  time_ce5 = 1, time_ce6 = 1, time_ce7 = 1, time_ce8 = 1
)

all_ce <- process_group(data_ce, lineage_data_top_ce, "ce")

# ============================================================
# CB group processing
# ============================================================
data_cb <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv") %>%
  select(cell, 31:36)

lineage_data_top_cb <- data.frame(
  parent = c("P0", "P0", "P1", "P1", "EMS", "EMS", "P2", "P2", "P3", "P3", "P4", "P4"),
  child  = c("AB", "P1", "EMS", "P2", "MS", "E", "P3", "C", "P4", "D", "Z2", "Z3"),
  time_cb1 = 1, time_cb2 = 1, time_cb3 = 1,
  time_cb4 = 1, time_cb5 = 1, time_cb6 = 1
)

all_cb <- process_group(data_cb, lineage_data_top_cb, "cb")

# ============================================================
# CN group processing
# ============================================================
data_cn <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv") %>%
  select(cell, 37:43)

lineage_data_top_cn <- data.frame(
  parent = c("P0", "P0", "P1", "P1", "EMS", "EMS", "P2", "P2", "P3", "P3", "P4", "P4"),
  child  = c("AB", "P1", "EMS", "P2", "MS", "E", "P3", "C", "P4", "D", "Z2", "Z3"),
  time_cn1 = 1, time_cn2 = 1, time_cn3 = 1,
  time_cn4 = 1, time_cn5 = 1, time_cn6 = 1, time_cn7 = 1
)

all_cn <- process_group(data_cn, lineage_data_top_cn, "cn")

# ============================================================
# Merge CE/CB/CN summaries
# ============================================================
merged_data <- all_ce %>%
  full_join(all_cb, by = "parent") %>%
  full_join(all_cn, by = "parent")

merged_data <- merged_data %>%
  mutate(
    ce_mean = round(apply(select(., diff_ce1:diff_ce8), 1, mean, na.rm = TRUE), 2),
    cb_mean = round(apply(select(., diff_cb1:diff_cb6), 1, mean, na.rm = TRUE), 2),
    cn_mean = round(apply(select(., diff_cn1:diff_cn7), 1, mean, na.rm = TRUE), 2)
  ) %>%
  select(1, 23:25)

data1 <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv") %>%
  select(1, 44)

colnames(merged_data)[1] <- "cell"

data <- merged_data %>%
  inner_join(data1, by = "cell")

# ============================================================
# Prepare plotting table (mean asymmetry)
# ============================================================
data_plot <- data %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "Stage_AB4",
        "Stage_AB8",
        "Stage_AB16",
        "Stage_AB32",
        "Stage_AB64",
        "Stage_AB128",
        "Stage_AB256"
      )
    ),
    cell_letter = substring(cell, 1, 1),
    cell_nchar = nchar(cell),
    cell_letter = factor(cell_letter, levels = c("A", "M", "E", "C", "D", "P"))
  ) %>%
  arrange(stage, cell_letter, cell_nchar)

#write.csv(data_plot, "SuppTable_data_sa_value.csv", row.names = FALSE)

data_plot$cell <- factor(data_plot$cell, levels = unique(data_plot$cell))

y_axis_order <- c("cn_mean", "cb_mean", "ce_mean")

highlight_stages <- c("Stage_AB4", "Stage_AB16", "Stage_AB64", "Stage_AB256")

highlight_rects <- data_plot %>%
  filter(stage %in% highlight_stages) %>%
  group_by(stage) %>%
  summarize(
    xmin = min(as.numeric(cell)),
    xmax = max(as.numeric(cell)),
    ymin = 0.5,
    ymax = length(y_axis_order) + 0.5
  )

data_plot_long <- data_plot %>%
  pivot_longer(
    cols = 2:4,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(variable = factor(variable, levels = y_axis_order))

# ============================================================
# Heatmap visualization
# ============================================================
ggplot(data_plot_long, aes(x = cell, y = variable, fill = value)) +
  geom_tile() +
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  labs(x = "Cell", y = "Variable", title = "") +
  scale_fill_gradientn(
    colors = c("#2166AC", "#92C5DE", "white", "#FDDBC7", "#B2182B"),
    values = c(0, 0.25, 0.6, 0.75, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )

# ggsave(
#   filename = "asymmtry_mean_value.pdf",
#   plot = last_plot(),
#   width = 10,
#   height = 3,
#   dpi = 300
# )



# ==========================================================
# Figure 3: Asymmetry significance (q-value heatmap for CB / CN vs CE)
# ==========================================================


# ----------------------------------------------------------
# Load libraries
# ----------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# ----------------------------------------------------------
# 1. Define processing function (same as before)
# ----------------------------------------------------------
process_group <- function(data, lineage_data_top, new_prefix) {
  pattern <- switch(new_prefix,
                    ce = "corTim_cep",
                    cb = "corTim_cbp",
                    cn = "corTim_cnp")
  
  data <- data %>%
    rename_with(~ sub(pattern, paste0("time_", new_prefix), .),
                starts_with(pattern))
  
  data_processed <- data %>%
    mutate(parent = substr(cell, 1, nchar(cell)-1)) %>%
    filter(nchar(cell) > 1, !cell %in% c("M","P2","P3","P4","Z2","Z3","EMS","MS")) %>%
    rename(child = cell) %>%
    bind_rows(lineage_data_top) %>%
    filter(!(parent %in% c("P0","P1","P4")))
  
  time_cols <- grep(paste0("^time_", new_prefix), names(data_processed), value=TRUE)
  
  data_processed %>%
    group_by(parent) %>%
    summarise(
      across(all_of(time_cols),
             ~ {
               last_letters <- substr(child, nchar(child), nchar(child))
               pos_val <- .x[last_letters %in% c("p","r","v","3","4")]
               neg_val <- .x[last_letters %in% c("a","l","d","C","D")]
               if (length(pos_val)==1 && length(neg_val)==1) pos_val - neg_val else NA_real_
             },
             .names = "diff_{sub('time_', '', .col)}"),
      .groups="drop"
    ) %>%
    filter(!if_any(starts_with(paste0("diff_", new_prefix)), is.na))
}

# ----------------------------------------------------------
# 2. Load raw data
# ----------------------------------------------------------
raw <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")

# CE group
data_ce <- raw %>% select(cell, 23:30)
lineage_data_top_ce <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_ce1=1,time_ce2=1,time_ce3=1,time_ce4=1,
  time_ce5=1,time_ce6=1,time_ce7=1,time_ce8=1
)
all_ce <- process_group(data_ce, lineage_data_top_ce, "ce")

# CB group
data_cb <- raw %>% select(cell, 31:36)
lineage_data_top_cb <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_cb1=1,time_cb2=1,time_cb3=1,time_cb4=1,time_cb5=1,time_cb6=1
)
all_cb <- process_group(data_cb, lineage_data_top_cb, "cb")

# CN group
data_cn <- raw %>% select(cell, 37:43)
lineage_data_top_cn <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_cn1=1,time_cn2=1,time_cn3=1,time_cn4=1,
  time_cn5=1,time_cn6=1,time_cn7=1
)
all_cn <- process_group(data_cn, lineage_data_top_cn, "cn")

# ----------------------------------------------------------
# 3. Merge CE / CB / CN differences
# ----------------------------------------------------------
merged <- all_ce %>%
  full_join(all_cb, by="parent") %>%
  full_join(all_cn, by="parent") %>%
  mutate(
    ce_mean = round(rowMeans(select(., diff_ce1:diff_ce8), na.rm=TRUE), 2),
    cb_mean = round(rowMeans(select(., diff_cb1:diff_cb6), na.rm=TRUE), 2),
    cn_mean = round(rowMeans(select(., diff_cn1:diff_cn7), na.rm=TRUE), 2),
    ce_sd   = round(apply(select(., diff_ce1:diff_ce8), 1, sd, na.rm=TRUE), 2)
  ) %>%
  select(parent, ce_mean, cb_mean, cn_mean, ce_sd) %>%
  rename(cell=parent)

stage <- raw %>% select(cell, stage=44)
data <- inner_join(merged, stage, by="cell")

# ----------------------------------------------------------
# 4. Compute z-scores, p/q-values, significance
# ----------------------------------------------------------
data <- data %>%
  mutate(across(c(cb_mean, cn_mean),
                list(zscore = ~ (. - ce_mean)/ce_sd),
                .names="{.col}_{.fn}"))

# Keep only cell + zscores + stage
data <- data %>% select(cell, cb_mean_zscore, cn_mean_zscore, stage, ce_sd)

# Compute p/q-values and significance (q_bh)
z_cols <- c("cb_mean_zscore","cn_mean_zscore")

for (zname in z_cols) {
  z <- data[[zname]]
  p_val <- 2*pnorm(-abs(z))
  p_dir <- ifelse(is.na(z), NA, ifelse(z<0,-1, ifelse(z>0,1, NA)))
  q_val <- p.adjust(p_val, method="fdr")
  q_bh  <- ifelse(is.na(q_val), NA,
                  ifelse(q_val<0.01 & p_dir==1, 1,
                         ifelse(q_val<0.01 & p_dir==-1, -1, 0)))
  
  data[[paste0("p_",zname)]] <- p_val
  data[[paste0("q_",zname)]] <- q_val
  data[[paste0("q_bh_",zname)]] <- q_bh
  data[[paste0("pdir_",zname)]] <- p_dir
}

# Save tables
write.csv(merged, "SuppTable_sa_qvalue_1.csv", row.names=FALSE)
data_plot <- data %>% select(cell, stage, starts_with("q_bh"))
write.csv(data_plot, "SuppTable_sa_qvalue_2.csv", row.names=FALSE)

# ----------------------------------------------------------
# 5. Prepare data for heatmap
# ----------------------------------------------------------
data_plot <- data_plot %>%
  mutate(stage=factor(stage, levels=c("Stage_AB4","Stage_AB8","Stage_AB16",
                                      "Stage_AB32","Stage_AB64","Stage_AB128","Stage_AB256")),
         cell_letter=factor(substr(cell,1,1), levels=c("A","M","E","C","D","P")),
         cell_nchar=nchar(cell)) %>%
  arrange(stage, cell_letter, cell_nchar)

data_plot$cell <- factor(data_plot$cell, levels=unique(data_plot$cell))

# Reshape to long format
y_axis_order <- c("q_bh_cn_mean_zscore","q_bh_cb_mean_zscore")
data_long <- data_plot %>%
  pivot_longer(cols=starts_with("q_bh"),
               names_to="variable", values_to="value") %>%
  mutate(value=factor(value, levels=c("-1","0","1",NA)),
         variable=factor(variable, levels=y_axis_order))

# ----------------------------------------------------------
# 6. Plot heatmap
# ----------------------------------------------------------
p <- ggplot(data_long, aes(x=cell, y=variable, fill=value)) +
  geom_tile() +
  facet_grid(.~stage, scales="free_x", space="free") +
  scale_fill_manual(values=c("1"="#9c2831","0"="#f5efe7","-1"="navy"),
                    na.value="gray", name="Value") +
  labs(x="Cell", y="") +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=8),
    legend.position="right",
    panel.border=element_rect(color="black", fill=NA, linewidth=1),
    panel.background=element_blank(),
    panel.grid=element_blank(),
    plot.background=element_blank()
  )

print(p)

#ggsave("qua_asymmetry_0.01_cbcn_mean.pdf", plot=p, width=10, height=3, dpi=300)



# ==========================================================
# Comparison of asymmetry mean and variability (CE / CB / CN)
# ==========================================================

setwd("~/Desktop/cbcn/draft_code/new/")

# ----------------------------------------------------------
# Load libraries
# ----------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)   # for stat_compare_means

# ----------------------------------------------------------
# 1. Define function to process CE / CB / CN
# ----------------------------------------------------------
process_group <- function(data, lineage_data_top, new_prefix) {
  pattern <- switch(new_prefix,
                    ce = "corTim_cep",
                    cb = "corTim_cbp",
                    cn = "corTim_cnp")
  
  # Rename time columns
  data <- data %>%
    rename_with(~ sub(pattern, paste0("time_", new_prefix), .),
                starts_with(pattern))
  
  # Preprocess
  data_processed <- data %>%
    mutate(parent = substr(cell, 1, nchar(cell)-1)) %>%
    filter(nchar(cell) > 1,
           !cell %in% c("M","P2","P3","P4","Z2","Z3","EMS","MS")) %>%
    rename(child = cell) %>%
    bind_rows(lineage_data_top) %>%
    filter(!(parent %in% c("P0","P1","P4")))
  
  time_cols <- grep(paste0("^time_", new_prefix), names(data_processed), value=TRUE)
  
  # Compute asymmetry (positive - negative daughter)
  data_processed %>%
    group_by(parent) %>%
    summarise(
      across(all_of(time_cols),
             ~{
               last_letters <- substr(child, nchar(child), nchar(child))
               pos_val <- .x[last_letters %in% c("p","r","v","3","4")]
               neg_val <- .x[last_letters %in% c("a","l","d","C","D")]
               if (length(pos_val)==1 && length(neg_val)==1) pos_val - neg_val else NA_real_
             },
             .names = "diff_{sub('time_', '', .col)}"),
      child = paste(unique(child), collapse=","),
      .groups="drop"
    ) %>%
    filter(!if_any(starts_with(paste0("diff_", new_prefix)), is.na))
}

# ----------------------------------------------------------
# 2. Load raw data
# ----------------------------------------------------------
raw <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")

# CE
data_ce <- raw %>% select(cell, 23:30)
lineage_data_top_ce <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_ce1=1,time_ce2=1,time_ce3=1,time_ce4=1,
  time_ce5=1,time_ce6=1,time_ce7=1,time_ce8=1
)
all_ce <- process_group(data_ce, lineage_data_top_ce, "ce")

# CB
data_cb <- raw %>% select(cell, 31:36)
lineage_data_top_cb <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_cb1=1,time_cb2=1,time_cb3=1,time_cb4=1,time_cb5=1,time_cb6=1
)
all_cb <- process_group(data_cb, lineage_data_top_cb, "cb")

# CN
data_cn <- raw %>% select(cell, 37:43)
lineage_data_top_cn <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_cn1=1,time_cn2=1,time_cn3=1,time_cn4=1,
  time_cn5=1,time_cn6=1,time_cn7=1
)
all_cn <- process_group(data_cn, lineage_data_top_cn, "cn")

# ----------------------------------------------------------
# 3. Merge CE / CB / CN mean values (F2D)
# ----------------------------------------------------------
merged <- all_ce %>%
  full_join(all_cb, by="parent") %>%
  full_join(all_cn, by="parent") %>%
  mutate(
    ce_mean = round(rowMeans(select(., diff_ce1:diff_ce8), na.rm=TRUE),2),
    cb_mean = round(rowMeans(select(., diff_cb1:diff_cb6), na.rm=TRUE),2),
    cn_mean = round(rowMeans(select(., diff_cn1:diff_cn7), na.rm=TRUE),2)
  ) %>%
  select(parent, ce_mean, cb_mean, cn_mean) %>%
  rename(cell=parent)

# Reshape
data_long <- merged %>%
  pivot_longer(cols = c(ce_mean, cb_mean, cn_mean),
               names_to="channel", values_to="value") %>%
  mutate(channel=factor(channel, levels=c("ce_mean","cb_mean","cn_mean")))

comparisons <- list(c("ce_mean","cb_mean"),
                    c("ce_mean","cn_mean"),
                    c("cb_mean","cn_mean"))

# Plot F2D
p_F2D <- ggplot(data_long, aes(x=channel, y=value, fill=channel)) +
  stat_boxplot(geom="errorbar", width=0.42, size=0.5) +
  geom_boxplot(width=0.65, outlier.shape=NA, color="black", size=0.5) +
  labs(x="Group", y="Mean asymmetry value") +
  scale_fill_manual(values=c("ce_mean"="#CD534C","cb_mean"="#0073C2","cn_mean"="#EFC000")) +
  theme_classic(base_size=14) +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=1),
        axis.line=element_blank(),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(size=12),
        legend.position="none") +
  scale_y_continuous(limits=c(-0.2,11), expand=c(0.01,0)) +
  stat_compare_means(comparisons=comparisons, method="wilcox.test",
                     label="p.signif", label.y=rep(5.5, length(comparisons)))
print(p_F2D)


# ----------------------------------------------------------
# 4. Stage-specific variability (F2E, example: AB128)
# ----------------------------------------------------------
merged_sd <- all_ce %>%
  full_join(all_cb, by="parent") %>%
  full_join(all_cn, by="parent") %>%
  mutate(
    ce_sd = round(apply(select(., diff_ce1:diff_ce8),1,sd,na.rm=TRUE),2),
    cb_sd = round(apply(select(., diff_cb1:diff_cb6),1,sd,na.rm=TRUE),2),
    cn_sd = round(apply(select(., diff_cn1:diff_cn7),1,sd,na.rm=TRUE),2)
  ) %>%
  select(parent, ce_sd, cb_sd, cn_sd) %>%
  rename(cell=parent)

stage <- raw %>% select(cell, stage=44)
data_stage <- merged_sd %>%
  inner_join(stage, by="cell") %>%
  filter(stage=="Stage_AB128")

data_long <- pivot_longer(data_stage, cols=c("ce_sd","cb_sd","cn_sd"),
                          names_to="Group", values_to="Value")
data_long$Group <- factor(gsub("_sd","",data_long$Group), levels=c("ce","cb","cn"))

p_F2E <- ggplot(data_long, aes(x=Group, y=Value, fill=Group)) +
  stat_boxplot(geom="errorbar", width=0.42, size=0.5) +
  geom_boxplot(width=0.65, outlier.shape=NA, color="black", size=0.5) +
  labs(x="Group", y="SD of asymmetry (Stage_AB128)") +
  scale_fill_manual(values=c("ce"="#CD534C","cb"="#0073C2","cn"="#EFC000")) +
  theme_classic(base_size=14) +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=1),
        axis.line=element_blank(),
        axis.text=element_text(size=12),
        legend.position="none") +
  scale_y_continuous(limits=c(0,12), expand=c(0.01,0))
print(p_F2E)

# ----------------------------------------------------------
# 4. All-stage variability 
# ----------------------------------------------------------
merged_sd <- all_ce %>%
  full_join(all_cb, by="parent") %>%
  full_join(all_cn, by="parent") %>%
  mutate(
    ce_sd = round(apply(select(., diff_ce1:diff_ce8),1,sd,na.rm=TRUE),2),
    cb_sd = round(apply(select(., diff_cb1:diff_cb6),1,sd,na.rm=TRUE),2),
    cn_sd = round(apply(select(., diff_cn1:diff_cn7),1,sd,na.rm=TRUE),2)
  ) %>%
  select(parent, ce_sd, cb_sd, cn_sd) %>%
  rename(cell=parent)


data_long <- pivot_longer(merged_sd, cols=c("ce_sd","cb_sd","cn_sd"),
                          names_to="Group", values_to="Value")
data_long$Group <- factor(gsub("_sd","",data_long$Group), levels=c("ce","cb","cn"))

p_F2E <- ggplot(data_long, aes(x=Group, y=Value, fill=Group)) +
  stat_boxplot(geom="errorbar", width=0.42, size=0.5) +
  geom_boxplot(width=0.65, outlier.shape=NA, color="black", size=0.5) +
  labs(x="Group", y="SD of asymmetry (Stage_AB128)") +
  scale_fill_manual(values=c("ce"="#CD534C","cb"="#0073C2","cn"="#EFC000")) +
  theme_classic(base_size=14) +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=1),
        axis.line=element_blank(),
        axis.text=element_text(size=12),
        legend.position="none") +
  scale_y_continuous(limits=c(0,6), expand=c(0.1,0.01))
print(p_F2E)

# ----------------------------------------------------------
# 5. Growth vs specification variability (F2F)
# ----------------------------------------------------------
tissue <- read.table("~/Desktop/cbcn/linage/AllLineage.tsv", header=TRUE, sep="\t") %>%
  select(parent=1, CellFate=3)

data_sep <- merged_sd %>%
  inner_join(all_ce %>% select(parent, child), by=c("cell"="parent")) %>%
  separate(child, into=c("child1","child2"), sep=",", remove=FALSE) %>%
  left_join(tissue, by=c("child1"="parent")) %>% rename(fate1=CellFate) %>%
  left_join(tissue, by=c("child2"="parent")) %>% rename(fate2=CellFate) %>%
  mutate(Group=if_else(fate1==fate2,"growth","specification"))

data_long <- pivot_longer(data_sep, cols=c("ce_sd","cb_sd","cn_sd"),
                          names_to="measure", values_to="value") %>%
  mutate(group_measure=paste(Group,measure,sep="_"),
         group_measure=factor(group_measure,
                              levels=c("growth_ce_sd","growth_cb_sd","growth_cn_sd",
                                       "specification_ce_sd","specification_cb_sd","specification_cn_sd")))

# F2F: compare ce/cb/cn inside growth and inside specification
comparisons <- list(
  c("growth_ce_sd", "growth_cb_sd"),
  c("growth_ce_sd", "growth_cn_sd"),
  c("growth_cb_sd", "growth_cn_sd"),
  c("specification_ce_sd", "specification_cb_sd"),
  c("specification_ce_sd", "specification_cn_sd"),
  c("specification_cb_sd", "specification_cn_sd")
)

p_F2F <- ggplot(data_long, aes(x=group_measure, y=value, fill=measure)) +
  stat_boxplot(geom="errorbar", width=0.42, size=0.5) +
  geom_boxplot(width=0.65, outlier.shape=NA, color="black", size=0.5) +
  labs(x="Group and measure", y="SD of asymmetry") +
  scale_fill_manual(values=c("ce_sd"="#CD534C","cb_sd"="#0073C2","cn_sd"="#EFC000")) +
  theme_classic(base_size=14) +
  theme(panel.border=element_rect(color="black", fill=NA, linewidth=1),
        axis.line   = element_blank(),
        axis.text=element_text(size=12),
        legend.position="none") +
  scale_y_continuous(limits=c(0,8), expand=c(0.01,0)) +
  stat_compare_means(comparisons=comparisons, method="wilcox.test",
                     label="p.signif", label.y=rep(7,6))
print(p_F2F)

# ----------------------------------------------------------
# Save plots
# ----------------------------------------------------------
# ggsave("F2D_asymmetry_mean.pdf", p_F2D, width=5, height=5)
# ggsave("F2E_asymmetry_sd_stage128.pdf", p_F2E, width=5, height=5)
# ggsave("F2F_asymmetry_sd_growth_spec.pdf", p_F2F, width=6, height=5)





# ==========================================================
# Asymmetry variability (SD heatmap)
# ==========================================================
# ----------------------------------------------------------
# Load libraries
# ----------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# ----------------------------------------------------------
# 1. Define processing function
# ----------------------------------------------------------
process_group <- function(data, lineage_data_top, new_prefix) {
  pattern <- switch(new_prefix,
                    ce = "corTim_cep",
                    cb = "corTim_cbp",
                    cn = "corTim_cnp")
  
  data <- data %>%
    rename_with(~ sub(pattern, paste0("time_", new_prefix), .),
                starts_with(pattern))
  
  data_processed <- data %>%
    mutate(parent = substr(cell, 1, nchar(cell)-1)) %>%
    filter(nchar(cell) > 1,
           !cell %in% c("M","P2","P3","P4","Z2","Z3","EMS","MS")) %>%
    rename(child = cell) %>%
    bind_rows(lineage_data_top) %>%
    filter(!(parent %in% c("P0","P1","P4")))
  
  time_cols <- grep(paste0("^time_", new_prefix), names(data_processed), value=TRUE)
  
  data_processed %>%
    group_by(parent) %>%
    summarise(
      across(all_of(time_cols),
             ~{
               last_letters <- substr(child, nchar(child), nchar(child))
               pos_val <- .x[last_letters %in% c("p","r","v","3","4")]
               neg_val <- .x[last_letters %in% c("a","l","d","C","D")]
               if (length(pos_val)==1 && length(neg_val)==1) pos_val - neg_val else NA_real_
             },
             .names = "diff_{sub('time_', '', .col)}"),
      .groups="drop"
    ) %>%
    filter(!if_any(starts_with(paste0("diff_", new_prefix)), is.na))
}

# ----------------------------------------------------------
# 2. Load CE / CB / CN data
# ----------------------------------------------------------
raw <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")

# CE
data_ce <- raw %>% select(cell, 23:30)
lineage_data_top_ce <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_ce1=1,time_ce2=1,time_ce3=1,time_ce4=1,
  time_ce5=1,time_ce6=1,time_ce7=1,time_ce8=1
)
all_ce <- process_group(data_ce, lineage_data_top_ce, "ce")

# CB
data_cb <- raw %>% select(cell, 31:36)
lineage_data_top_cb <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_cb1=1,time_cb2=1,time_cb3=1,time_cb4=1,time_cb5=1,time_cb6=1
)
all_cb <- process_group(data_cb, lineage_data_top_cb, "cb")

# CN
data_cn <- raw %>% select(cell, 37:43)
lineage_data_top_cn <- data.frame(
  parent=c("P0","P0","P1","P1","EMS","EMS","P2","P2","P3","P3","P4","P4"),
  child =c("AB","P1","EMS","P2","MS","E","P3","C","P4","D","Z2","Z3"),
  time_cn1=1,time_cn2=1,time_cn3=1,time_cn4=1,
  time_cn5=1,time_cn6=1,time_cn7=1
)
all_cn <- process_group(data_cn, lineage_data_top_cn, "cn")

# ----------------------------------------------------------
# 3. Merge and compute SD
# ----------------------------------------------------------
merged <- all_ce %>%
  full_join(all_cb, by="parent") %>%
  full_join(all_cn, by="parent") %>%
  mutate(
    ce_sd = round(apply(select(., diff_ce1:diff_ce8),1,sd,na.rm=TRUE),2),
    cb_sd = round(apply(select(., diff_cb1:diff_cb6),1,sd,na.rm=TRUE),2),
    cn_sd = round(apply(select(., diff_cn1:diff_cn7),1,sd,na.rm=TRUE),2)
  ) %>%
  select(parent, ce_sd, cb_sd, cn_sd) %>%
  rename(cell=parent)

stage <- raw %>% select(cell, stage=44)
data <- merged %>% inner_join(stage, by="cell")

# ----------------------------------------------------------
# 4. Prepare for plotting
# ----------------------------------------------------------
data_plot <- data %>%
  mutate(stage=factor(stage, levels=c("Stage_AB4","Stage_AB8","Stage_AB16",
                                      "Stage_AB32","Stage_AB64","Stage_AB128","Stage_AB256")),
         cell_letter=substr(cell,1,1),
         cell_nchar=nchar(cell),
         cell_letter=factor(cell_letter, levels=c("A","M","E","C","D","P"))) %>%
  arrange(stage, cell_letter, cell_nchar)

data_plot$cell <- factor(data_plot$cell, levels=unique(data_plot$cell))

y_axis_order <- c("cn_sd","cb_sd","ce_sd")

data_long <- data_plot %>%
  pivot_longer(cols=c("ce_sd","cb_sd","cn_sd"),
               names_to="variable", values_to="value") %>%
  mutate(variable=factor(variable, levels=y_axis_order))

# ----------------------------------------------------------
# 5. Heatmap
# ----------------------------------------------------------
p <- ggplot(data_long, aes(x=cell, y=variable, fill=value)) +
  geom_tile() +
  facet_grid(.~stage, scales="free_x", space="free") +
  scale_fill_gradientn(colors=c("#2166AC","#92C5DE","#D1E5F0",
                                "#FDDBC7","#D6604D","#B2182B"),
                       values=c(0,0.15,0.3,0.4,0.75,1)) +
  labs(x="Cell", y="SD of asymmetry") +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=10),
    legend.position="right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background=element_blank(),
    panel.grid=element_blank(),
    plot.background=element_blank()
  )

print(p)

#ggsave("f3_sup8_asymmetry_sd_inner.pdf", plot=p, width=8, height=2, dpi=300)




