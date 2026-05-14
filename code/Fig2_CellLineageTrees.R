
# =============================================================================
# PLOT_Lineage Tree Visualization_CASE STUDY
# =============================================================================
# =============================================================================
# Cell Lineage Tree Visualization with Asynchrony Coloring
# =============================================================================

library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(magrittr)
library(patchwork)

setwd("~/Desktop/cbcn/")

# =============================================================================
# Helper Functions
# =============================================================================

trim_cell <- function(cell) {
  c(str_sub(cell, 1, -2), str_sub(cell, -1))
}

xoffset <- function(node, root, klist = list(), slist = list()) {
  inter_node <- node
  route <- c()
  
  while (str_detect(inter_node, root) & inter_node != root) {
    cell <- trim_cell(inter_node)
    inter_node <- cell[1]
    tail_cell <- cell[2]
    
    k <- 1
    if (length(klist) > 0 && sum(str_detect(inter_node, names(klist))) > 0) {
      k <- prod(unlist(klist[str_detect(inter_node, names(klist))]))
    }
    
    route <- c(if (tail_cell == "a") -k else k, route)
  }
  
  xoff <- sum((0.5^(1:length(route))) * route)
  
  if (length(slist) > 0 && sum(str_detect(node, names(slist))) > 0) {
    xoff <- xoff + sum(unlist(slist[str_detect(node, names(slist))]))
  }
  
  return(xoff)
}

# =============================================================================
# Configuration
# =============================================================================

KLIST <- list(Za = 2, Zpap = 0.3, Zppa = 0.5, Zpppa = 0.8, Zppp = 0.8, Zpppp = 0.3)
SLIST <- list(Za = -0.5, Zp = -0.2, Zpa = 0.22, Zpap = -0.08, Zpp = -0.04,
              Zppa = 0.04, Zppp = -0.06, Zpppa = 0.01, Zpppp = -0.03)
ID_SHOW_LEN <- 3
BG_COLOR <- "darkgray"

# Cell stages
STAGE_RAW <- c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2")
STAGE_AB4 <- c("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")
STAGE_AB8 <- c("ABala", "ABalp", "ABara", "ABarp", "ABpla", "ABplp", "ABpra", "ABprp",
               "MSa", "MSp", "Ca", "Cp", "P3")
CELLS_TO_KEEP <- c(STAGE_RAW, STAGE_AB4, STAGE_AB8)

# =============================================================================
# Load Data
# =============================================================================

AllLineage <- read.csv("CellID.csv") %>%
  tibble() %>%
  filter(CellName %in% CELLS_TO_KEEP)

cat("Death cells:", sum(AllLineage$CellFate == "Death", na.rm = TRUE), "\n")

df_asynchrony <- read.csv("asynchrony_0.01.csv") %>%
  select(CellName = 1, corTim_cep5, corTim_cbp3, corTim_cnp1) %>%
  mutate(
    corTim_cep5 = case_when(CellName == "P3" ~ 2L, CellName %in% STAGE_AB8 ~ 1L, TRUE ~ 0L),
    corTim_cbp3 = case_when(CellName == "P3" ~ 2L, CellName %in% STAGE_AB8 ~ 1L, TRUE ~ 0L),
    corTim_cnp1 = case_when(CellName == "P3" ~ 2L, CellName %in% STAGE_AB8 ~ 1L, TRUE ~ 0L)
  )

# =============================================================================
# Main Plotting Function
# =============================================================================

process_lineage_and_plot <- function(lineage_data, slope, std_column, time_limit, time_factor, data_type = "cb") {
  
  # Process lineage data
  lineage <- lineage_data %>%
    as_tibble() %>%
    filter(time <= time_limit) %>%
    mutate(time = time * time_factor) %>%
    group_by(cell) %>%
    summarize(Start = min(time) / slope, End = max(time) / slope, .groups = "drop") %>%
    rename(CellName = cell) %>%
    filter(CellName %in% CELLS_TO_KEEP) %>%
    left_join(AllLineage %>% select(CellName, ID, CellFate), by = "CellName") %>%
    rename(Fate = CellFate) %>%
    left_join(df_asynchrony %>% select(CellName, !!sym(std_column)), by = "CellName")
  
  # Add root nodes
  lineage %<>%
    add_row(CellName = "P0", Start = -20, End = -5, ID = "Z", Fate = "Unspecified") %>%
    add_row(CellName = "AB", Start = -5, End = 0, ID = "Za", Fate = "Unspecified") %>%
    add_row(CellName = "P1", Start = -5, End = 0, ID = "Zp", Fate = "Unspecified") %>%
    distinct(CellName, .keep_all = TRUE)
  
  # Calculate positions
  p.data <- lineage %>%
    mutate(
      End = End + 1,
      k = if_else(str_detect(ID, "Zppa"), 0.3, 1)
    ) %>%
    rowwise() %>%
    mutate(
      x = xoffset(ID, "Z", KLIST, SLIST),
      parentx = xoffset(str_sub(ID, end = -2), "Z", KLIST, SLIST)
    ) %>%
    ungroup()
  
  # Assign colors
  p.data$node_colors <- ifelse(p.data[[std_column]] == 1, "blue",
                               ifelse(p.data[[std_column]] == 2, "red", NA))
  
  # Prepare data for vertical lines
  p.data.y <- p.data %>%
    select(CellName, ID, x, Start, End, node_colors, Fate) %>%
    mutate(Mid = Start + 0.03) %>%
    melt(id.vars = c("CellName", "ID", "x", "Mid", "node_colors", "Fate"),
         variable.name = "SE", value.name = "Posi")
  
  # Prepare data for horizontal lines
  p.data.x <- p.data %>%
    select(ID, x, Start, parentx) %>%
    melt(id.vars = c("ID", "Start"), variable.name = "FT", value.name = "Posi")
  
  # Create plot
  p <- p.data.y %>%
    mutate(label = if_else(nchar(ID) < ID_SHOW_LEN, CellName, NA)) %>%
    ggplot() +
    geom_path(aes(x = x, y = -Posi, group = ID), color = BG_COLOR) +
    geom_text(aes(x = x, y = -Mid, label = label)) +
    geom_path(data = p.data.x, aes(x = Posi, y = -Start, group = ID), color = BG_COLOR) +
    geom_segment(
      data = subset(p.data.y, SE == "End" & !is.na(node_colors)),
      aes(x = x - 0.03, xend = x + 0.03, y = -Posi, yend = -Posi, color = node_colors),
      linewidth = 0.8
    ) +
    scale_y_continuous(limits = c(-65, 20), labels = function(x) abs(x)) +
    guides(color = "none") +
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

# =============================================================================
# Generate Plots
# =============================================================================

p_ce <- process_lineage_and_plot(
  read.csv("CDFile/ce/CD200322plc1p2.csv"),
  0.95, "corTim_cep5", 195, 1.44, "ce"
)
print(p_ce)

p_cb <- process_lineage_and_plot(
  read.csv("CDFile/she1/CD240731cbhis72p3.csv"),
  0.9, "corTim_cbp3", 170, 1.57, "cb"
)
print(p_cb)

p_cn <- process_lineage_and_plot(
  read.csv("CDFile/cn/CD241202cnhis72p1.csv"),
  1.35, "corTim_cnp1", 235, 1.58, "cn"
)
print(p_cn)

# Combine plots
p_3f <- p_ce + p_cb + p_cn
print(p_3f)

# ggsave("lineage_plot_P3_rep.pdf", p_3f, width = 8, height = 5)






# =============================================================================
# ASYNCHRONY ANALYSIS_SEP 
# =============================================================================
# =============================================================================
# Asynchrony Analysis: Division Timing Matrix & Phenotype Detection
# =============================================================================
setwd("~/Desktop/cbcn/")

library(dplyr)
library(readr)
library(tidyr)
library(scales)
library(stats)
library(stringr)
library(pracma)
library(ggplot2)

# =============================================================================
# Configuration
# =============================================================================

CTR_timing_matrix_path <- "~/Desktop/cbcn/output_folder/0114/CTR_division_timing_matrix"
RNAi_timing_matrix_path <- "~/Desktop/cbcn/output_folder/0114/RNAi_division_timing_matrix"
Division_asynchrony_phenotype_folder <- "~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype"

if (!dir.exists(CTR_timing_matrix_path)) dir.create(CTR_timing_matrix_path)
if (!dir.exists(RNAi_timing_matrix_path)) dir.create(RNAi_timing_matrix_path)
if (!dir.exists(Division_asynchrony_phenotype_folder)) dir.create(Division_asynchrony_phenotype_folder)

# Cell Stages
Stage_AB4 <- c("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C")
Stage_AB8 <- c("ABala", "ABalp", "ABara", "ABarp", "ABpla", "ABplp", "ABpra", "ABprp", "MSa", "MSp", "Ca", "Cp", "P3")
Stage_AB16 <- c("ABalaa", "ABalap", "ABalpa", "ABalpp", "ABaraa", "ABarap", "ABarpa", "ABarpp", "ABplaa", "ABplap", "ABplpa", "ABplpp", "ABpraa", "ABprap", "ABprpa", "ABprpp", "MSaa", "MSap", "MSpa", "MSpp", "Ea", "Ep", "D")
Stage_AB32 <- c("ABalaaa", "ABalaap", "ABalapa", "ABalapp", "ABalpaa", "ABalpap", "ABalppa", "ABalppp", "ABaraaa", "ABaraap", "ABarapa", "ABarapp", "ABarpaa", "ABarpap", "ABarppa", "ABarppp", "ABplaaa", "ABplaap", "ABplapa", "ABplapp", "ABplpaa", "ABplpap", "ABplppa", "ABplppp", "ABpraaa", "ABpraap", "ABprapa", "ABprapp", "ABprpaa", "ABprpap", "ABprppa", "ABprppp", "MSaaa", "MSaap", "MSapa", "MSapp", "MSpaa", "MSpap", "MSppa", "MSppp", "Eal", "Ear", "Epl", "Epr", "Caa", "Cap", "Cpa", "Cpp", "Da", "Dp", "P4")
Stage_AB64 <- c("ABalaaaa", "ABalaaap", "ABalaapa", "ABalaapp", "ABalapaa", "ABalapap", "ABalappa", "ABalappp", "ABalpaaa", "ABalpaap", "ABalpapa", "ABalpapp", "ABalppaa", "ABalppap", "ABalpppa", "ABalpppp", "ABaraaaa", "ABaraaap", "ABaraapa", "ABaraapp", "ABarapaa", "ABarapap", "ABarappa", "ABarappp", "ABarpaaa", "ABarpaap", "ABarpapa", "ABarpapp", "ABarppaa", "ABarppap", "ABarpppa", "ABarpppp", "ABplaaaa", "ABplaaap", "ABplaapa", "ABplaapp", "ABplapaa", "ABplapap", "ABplappa", "ABplappp", "ABplpaaa", "ABplpaap", "ABplpapa", "ABplpapp", "ABplppaa", "ABplppap", "ABplpppa", "ABplpppp", "ABpraaaa", "ABpraaap", "ABpraapa", "ABpraapp", "ABprapaa", "ABprapap", "ABprappa", "ABprappp", "ABprpaaa", "ABprpaap", "ABprpapa", "ABprpapp", "ABprppaa", "ABprppap", "ABprpppa", "ABprpppp", "MSaaaa", "MSaaap", "MSaapa", "MSaapp", "MSapaa", "MSapap", "MSappa", "MSappp", "MSpaaa", "MSpaap", "MSpapa", "MSpapp", "MSppaa", "MSppap", "MSpppa", "MSpppp", "Caaa", "Caap", "Capa", "Capp", "Cpaa", "Cpap", "Cppa", "Cppp")
Stage_AB128 <- c("ABalaaaal", "ABalaaaar", "ABalaaapa", "ABalaaapp", "ABalaapaa", "ABalaapap", "ABalaappa", "ABalaappp", "ABalapaaa", "ABalapaap", "ABalapapa", "ABalapapp", "ABalappaa", "ABalappap", "ABalapppa", "ABalapppp", "ABalpaaaa", "ABalpaaap", "ABalpaapa", "ABalpaapp", "ABalpapaa", "ABalpapap", "ABalpappa", "ABalpappp", "ABalppaaa", "ABalppaap", "ABalppapa", "ABalppapp", "ABalpppaa", "ABalpppap", "ABalppppa", "ABalppppp", "ABaraaaaa", "ABaraaaap", "ABaraaapa", "ABaraaapp", "ABaraapaa", "ABaraapap", "ABaraappa", "ABaraappp", "ABarapaaa", "ABarapaap", "ABarapapa", "ABarapapp", "ABarappaa", "ABarappap", "ABarapppa", "ABarapppp", "ABarpaaaa", "ABarpaaap", "ABarpaapa", "ABarpaapp", "ABarpapaa", "ABarpapap", "ABarpappa", "ABarpappp", "ABarppaaa", "ABarppaap", "ABarppapa", "ABarppapp", "ABarpppaa", "ABarpppap", "ABarppppa", "ABarppppp", "ABplaaaaa", "ABplaaaap", "ABplaaapa", "ABplaaapp", "ABplaapaa", "ABplaapap", "ABplaappa", "ABplaappp", "ABplapaaa", "ABplapaap", "ABplapapa", "ABplapapp", "ABplappaa", "ABplappap", "ABplapppa", "ABplapppp", "ABplpaaaa", "ABplpaaap", "ABplpaapa", "ABplpaapp", "ABplpapaa", "ABplpapap", "ABplpappa", "ABplpappp", "ABplppaaa", "ABplppaap", "ABplppapa", "ABplppapp", "ABplpppaa", "ABplpppap", "ABplppppa", "ABplppppp", "ABpraaaaa", "ABpraaaap", "ABpraaapa", "ABpraaapp", "ABpraapaa", "ABpraapap", "ABpraappa", "ABpraappp", "ABprapaaa", "ABprapaap", "ABprapapa", "ABprapapp", "ABprappaa", "ABprappap", "ABprapppa", "ABprapppp", "ABprpaaaa", "ABprpaaap", "ABprpaapa", "ABprpaapp", "ABprpapaa", "ABprpapap", "ABprpappa", "ABprpappp", "ABprppaaa", "ABprppaap", "ABprppapa", "ABprppapp", "ABprpppaa", "ABprpppap", "ABprppppa", "ABprppppp", "MSaaaaa", "MSaaaap", "MSaaapa", "MSaapaa", "MSaapap", "MSapaaa", "MSapapa", "MSapapp", "MSpaaaa", "MSpaaap", "MSpaapa", "MSpapaa", "MSpapap", "MSppaaa", "MSppapa", "MSppapp", "Eala", "Ealp", "Eara", "Earp", "Epla", "Eplp", "Epra", "Eprp", "Caaaa", "Caaap", "Caapp", "Capaa", "Capap", "Cappa", "Cappp", "Cpaaa", "Cpaap", "Cpapa", "Cpapp", "Cppaa", "Cppap", "Cpppa", "Cpppp", "Daa", "Dap", "Dpa", "Dpp")
Stage_AB256 <- c("ABalaaappr", "ABalaapppa", "ABalaapppp", "ABalapaapp", "ABalapappa", "ABalappapp", "ABalapppap", "ABalappppa", "ABalappppp", "ABalpaapaa", "ABalpaapap", "ABalpaappa", "ABalpaappp", "ABalpapaaa", "ABalpapaap", "ABalpapapa", "ABalpapapp", "ABalpappaa", "ABalpappap", "ABalpapppp", "ABalppapaa", "ABalppapap", "ABalppappa", "ABalppappp", "ABalpppapa", "ABalpppapp", "ABalppppaa", "ABalppppap", "ABalpppppa", "ABalpppppp", "ABaraaapaa", "ABaraaapap", "ABaraaappa", "ABaraaappp", "ABaraapaaa", "ABaraapaap", "ABaraapapa", "ABaraapapp", "ABaraappaa", "ABaraappap", "ABaraapppa", "ABaraapppp", "ABarapaapp", "ABarapapap", "ABarapappp", "ABarappaap", "ABarappapa", "ABarappapp", "ABarapppaa", "ABarapppap", "ABarappppa", "ABarappppp", "ABarpapaap", "ABplaaaaap", "ABplaapaap", "ABplaapapa", "ABplaapapp", "ABplapaaaa", "ABplapaaap", "ABplapappp", "ABplapppap", "ABplpaaaaa", "ABplpaaaap", "ABplpaaapa", "ABplpaaapp", "ABplpaapaa", "ABplpaapap", "ABplpaappa", "ABplpaappp", "ABplpapaaa", "ABplpapaap", "ABplpapapa", "ABplpapapp", "ABplpappaa", "ABplpapppa", "ABplppaaaa", "ABplppaapa", "ABplppaapp", "ABplppapaa", "ABplppapap", "ABplppappa", "ABplppappp", "ABplpppaaa", "ABplpppaap", "ABplpppapa", "ABplppppaa", "ABplppppap", "ABplpppppa", "ABplpppppp", "ABpraaappp", "ABpraapaap", "ABpraapapp", "ABprapaaap", "ABprpaaaaa", "ABprpaaaap", "ABprpaaapa", "ABprpaaapp", "ABprpaapaa", "ABprpaapap", "ABprpaappa", "ABprpaappp", "ABprpapaaa", "ABprpapaap", "ABprpapapa", "ABprpapapp", "ABprpappaa", "ABprpappap", "ABprpapppa", "ABprpapppp", "ABprppaaaa", "ABprppaapp", "ABprppapaa", "ABprppapap", "ABprppappa", "ABprppappp", "ABprpppaaa", "ABprpppaap", "ABprpppapa", "ABprppppaa", "ABprppppap", "ABprpppppa", "ABprpppppp", "Caapa", "Capaaa", "Capaap", "Capapa", "Capapp", "Cappaa", "Cappap", "Capppa", "Capppp", "Cppaaa", "Cppaap", "Cppapa", "Cppapp", "Cpppaa", "Cpppap", "Cppppa", "Cppppp", "Daaa", "Daap", "Dapa", "Dapp", "Dpaa", "Dpap", "Dppa", "Dppp", "MSaaapp", "MSaappa", "MSaappp", "MSapaap", "MSappaa", "MSappap", "MSapppa", "MSapppp", "MSpappa", "MSpappp", "MSppaap", "MSpppaa", "MSpppap", "MSppppa", "MSppppp")

stages <- list(Stage_AB4 = Stage_AB4, Stage_AB8 = Stage_AB8, Stage_AB16 = Stage_AB16,
               Stage_AB32 = Stage_AB32, Stage_AB64 = Stage_AB64, Stage_AB128 = Stage_AB128, Stage_AB256 = Stage_AB256)
all_cells <- unique(unlist(stages))

# =============================================================================
# Load Data
# =============================================================================

data <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
CTR_celldivision <- data %>% select(1, 23:30)
RNAi_celldivision <- data %>% select(1, 31:43)

# =============================================================================
# Step 1: Create CTR Timing Matrices
# =============================================================================

for (item in colnames(CTR_celldivision)[2:ncol(CTR_celldivision)]) {
  emb <- CTR_celldivision[, c("cell", item)]
  division_time_dist <- matrix(NA, nrow = length(all_cells), ncol = length(all_cells))
  rownames(division_time_dist) <- colnames(division_time_dist) <- all_cells
  
  for (stage_name in list(Stage_AB4, Stage_AB8, Stage_AB16, Stage_AB32, Stage_AB64, Stage_AB128, Stage_AB256)) {
    cell_list <- stage_name
    cell_value <- emb[emb$cell %in% cell_list, ]
    
    if (nrow(cell_value) == length(cell_list)) {
      stage_time_dist <- matrix(0, nrow = length(cell_list), ncol = length(cell_list))
      rownames(stage_time_dist) <- colnames(stage_time_dist) <- cell_list
      
      for (z in seq_along(cell_list)) {
        cell <- cell_list[z]
        value <- cell_value[cell_value$cell == cell, 2]
        values <- data.frame(cell = cell_value$cell, values = cell_value[, 2] - value)
        values_df <- data.frame(key = values[, 1], value = values[, 2])
        stage_time_dist <- cbind(stage_time_dist, values_df[order(rownames(stage_time_dist)), -1])
      }
      stage_time_dist <- stage_time_dist[, !colnames(stage_time_dist) %in% cell_list]
      max_value <- max(stage_time_dist, na.rm = TRUE)
      stage_time_dist <- stage_time_dist / max_value
      colnames(stage_time_dist) <- rownames(stage_time_dist)
    } else {
      stage_time_dist <- matrix(NA, nrow = length(cell_list), ncol = length(cell_list))
      rownames(stage_time_dist) <- colnames(stage_time_dist) <- cell_list
    }
    
    for (row in rownames(stage_time_dist)) {
      for (col in colnames(stage_time_dist)) {
        if (!is.na(stage_time_dist[row, col])) {
          division_time_dist[row, col] <- stage_time_dist[row, col]
        }
      }
    }
  }
  output_file <- file.path(CTR_timing_matrix_path, paste(item, "_time_Matrix_normlized.txt", sep = ""))
  write.table(division_time_dist, file = output_file, sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
}

# =============================================================================
# Step 2: Create RNAi Timing Matrices
# =============================================================================

for (item in colnames(RNAi_celldivision)[2:ncol(RNAi_celldivision)]) {
  emb <- RNAi_celldivision[, c("cell", item)]
  division_time_dist <- matrix(NA, nrow = length(all_cells), ncol = length(all_cells))
  rownames(division_time_dist) <- colnames(division_time_dist) <- all_cells
  
  for (stage_name in list(Stage_AB4, Stage_AB8, Stage_AB16, Stage_AB32, Stage_AB64, Stage_AB128, Stage_AB256)) {
    cell_list <- stage_name
    cell_value <- emb[emb$cell %in% cell_list, ]
    
    if (nrow(cell_value) == length(cell_list)) {
      stage_time_dist <- matrix(0, nrow = length(cell_list), ncol = length(cell_list))
      rownames(stage_time_dist) <- colnames(stage_time_dist) <- cell_list
      
      for (z in seq_along(cell_list)) {
        cell <- cell_list[z]
        value <- cell_value[cell_value$cell == cell, 2]
        values <- data.frame(cell = cell_value$cell, values = cell_value[, 2] - value)
        values_df <- data.frame(key = values[, 1], value = values[, 2])
        stage_time_dist <- cbind(stage_time_dist, values_df[order(rownames(stage_time_dist)), -1])
      }
      stage_time_dist <- stage_time_dist[, !colnames(stage_time_dist) %in% cell_list]
      max_value <- max(stage_time_dist, na.rm = TRUE)
      stage_time_dist <- stage_time_dist / max_value
      colnames(stage_time_dist) <- rownames(stage_time_dist)
    } else {
      stage_time_dist <- matrix(NA, nrow = length(cell_list), ncol = length(cell_list))
      rownames(stage_time_dist) <- colnames(stage_time_dist) <- cell_list
    }
    
    for (row in rownames(stage_time_dist)) {
      for (col in colnames(stage_time_dist)) {
        if (!is.na(stage_time_dist[row, col])) {
          division_time_dist[row, col] <- stage_time_dist[row, col]
        }
      }
    }
  }
  output_file <- file.path(RNAi_timing_matrix_path, paste(item, "_time_Matrix_normlized.txt", sep = ""))
  write.table(division_time_dist, file = output_file, sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
}

# =============================================================================
# Step 3: Asynchrony Phenotype Analysis
# =============================================================================

Embryo_dist_mean_WT <- matrix(NA, nrow = length(all_cells), ncol = length(all_cells))
rownames(Embryo_dist_mean_WT) <- colnames(Embryo_dist_mean_WT) <- all_cells

RNAi_cellposition_zscore <- data.frame(all_cells = unlist(all_cells))
rownames(RNAi_cellposition_zscore) <- RNAi_cellposition_zscore[, 1]

WT_single_cell_RMSD_all <- data.frame(all_cells = unlist(all_cells), cell_rmsd_mean = NA, cell_rmsd_sd = NA)
rownames(WT_single_cell_RMSD_all) <- WT_single_cell_RMSD_all[, 1]

RNAi_cellposition_qvalue_binary <- data.frame(all_cells = unlist(all_cells))
rownames(RNAi_cellposition_qvalue_binary) <- RNAi_cellposition_qvalue_binary[, 1]

RNAi_cellposition_qvalue <- data.frame(all_cells = unlist(all_cells))
rownames(RNAi_cellposition_qvalue) <- RNAi_cellposition_qvalue[, 1]

RNAi_single_cell_RMSD_all <- data.frame(all_cells = unlist(all_cells))
rownames(RNAi_single_cell_RMSD_all) <- RNAi_single_cell_RMSD_all[, 1]

for (stage in stages) {
  cell_list <- stage
  setwd(CTR_timing_matrix_path)
  
  File_list <- list.files(CTR_timing_matrix_path)
  Total_matrix <- NULL
  
  for (file in File_list) {
    Embryo_dist <- read.table(file, header = TRUE, sep = "\t", row.names = 1)
    Embryo_dist <- Embryo_dist[cell_list, cell_list]
    Total_matrix <- if (is.null(Total_matrix)) Embryo_dist else cbind(Total_matrix, Embryo_dist)
  }
  
  cell_rmsd_mean <- list()
  cell_rmsd_sd <- list()
  Embryo_dist_mean <- data.frame()
  
  for (i in seq_along(cell_list)) {
    cell <- cell_list[i]
    matching_columns <- grep(cell, colnames(Total_matrix))
    cell_vector_matrix <- Total_matrix[, matching_columns, drop = FALSE]
    cell_vector_mean <- rowMeans(cell_vector_matrix, na.rm = TRUE)
    
    cell_vector_matrix <- cell_vector_matrix[-which(rownames(cell_vector_matrix) == cell), , drop = FALSE]
    cell_rmsd <- dist(t(cell_vector_matrix)) / sqrt(ncol(t(cell_vector_matrix)))
    cell_rmsd <- cell_rmsd[!is.nan(cell_rmsd)]
    cell_rmsd <- cell_rmsd[!is.na(cell_rmsd) & !is.nan(cell_rmsd)]
    
    cell_rmsd_mean[cell] <- mean(cell_rmsd)
    cell_rmsd_sd[cell] <- sd(cell_rmsd)
    
    Embryo_dist_mean <- if (is.null(dim(Embryo_dist_mean)) || all(dim(Embryo_dist_mean) == 0)) {
      data.frame(cell_vector_mean)
    } else {
      cbind(Embryo_dist_mean, cell_vector_mean)
    }
  }
  colnames(Embryo_dist_mean) <- rownames(Embryo_dist_mean)
  
  for (row in rownames(Embryo_dist_mean)) {
    for (col in colnames(Embryo_dist_mean)) {
      if (!is.na(Embryo_dist_mean[row, col])) {
        Embryo_dist_mean_WT[row, col] <- Embryo_dist_mean[row, col]
      }
    }
  }
  
  key_dict_mean <- data.frame(key = names(cell_rmsd_mean), value = unlist(cell_rmsd_mean))
  for (row in rownames(key_dict_mean)) {
    if (!is.na(key_dict_mean[row, 2])) {
      WT_single_cell_RMSD_all[row, 2] <- key_dict_mean[row, 2]
    }
  }
  
  key_dict_sd <- data.frame(key = names(cell_rmsd_sd), value = unlist(cell_rmsd_sd))
  for (row in rownames(key_dict_sd)) {
    if (!is.na(key_dict_sd[row, 2])) {
      WT_single_cell_RMSD_all[row, 3] <- key_dict_sd[row, 2]
    }
  }
  
  # Read RNAi data
  setwd(RNAi_timing_matrix_path)
  embryo_list <- list()
  single_cell_rmsd_matrix <- NULL
  a <- 0
  
  for (file in list.files(RNAi_timing_matrix_path)) {
    Embryo <- strsplit(file, "_time_")[[1]][1]
    embryo_list <- append(embryo_list, Embryo)
    Embryo_dist <- read.table(file.path(RNAi_timing_matrix_path, file), header = TRUE, sep = "\t", row.names = 1)
    Embryo_dist <- Embryo_dist[cell_list, cell_list]
    
    single_cell_RMSD <- list()
    for (i in seq_along(cell_list)) {
      cell <- cell_list[i]
      cell_vactor_temp <- Embryo_dist[cell, , drop = FALSE]
      cell_vactor_mean <- Embryo_dist_mean[cell, , drop = FALSE]
      
      cell_vactor <- data.frame(cell_vactor_temp = unlist(cell_vactor_temp), cell_vactor_mean = unlist(cell_vactor_mean))
      cell_vactor <- cell_vactor[-which(rownames(cell_vactor) == cell), , drop = FALSE]
      cell_vactor <- cell_vactor[complete.cases(cell_vactor), ]
      
      if (nrow(cell_vactor) > 0) {
        cell_vactor_t <- t(cell_vactor)
        cell_rmsd <- as.vector(dist(cell_vactor_t) / sqrt(ncol(cell_vactor_t)))
        single_cell_RMSD[cell] <- cell_rmsd
      } else {
        single_cell_RMSD[cell] <- NA
      }
    }
    
    single_cell_rmsd_matrix <- if (is.list(single_cell_rmsd_matrix)) {
      as.data.frame(single_cell_RMSD)
    } else {
      cbind(single_cell_rmsd_matrix, unlist(single_cell_RMSD))
    }
    a <- a + 1
    print(a)
  }
  
  rownames(single_cell_rmsd_matrix) <- cell_list
  colnames(single_cell_rmsd_matrix) <- embryo_list
  
  for (row in rownames(single_cell_rmsd_matrix)) {
    for (col in colnames(single_cell_rmsd_matrix)) {
      if (!is.na(single_cell_rmsd_matrix[row, col])) {
        RNAi_single_cell_RMSD_all[row, col] <- single_cell_rmsd_matrix[row, col]
      }
    }
  }
  
  # Calculate z-score, p-value, q-value
  control_mean_values <- unlist(cell_rmsd_mean)
  control_std_values <- unlist(cell_rmsd_sd)
  
  diff_matrix <- sweep(single_cell_rmsd_matrix, 1, control_mean_values, FUN = "-")
  z_score <- sweep(diff_matrix, 1, control_std_values, FUN = "/")
  
  for (row in rownames(z_score)) {
    for (col in colnames(z_score)) {
      if (!is.na(z_score[row, col])) {
        RNAi_cellposition_zscore[row, col] <- z_score[row, col]
      }
    }
  }
  
  p_value <- pnorm(abs(z_score), lower.tail = FALSE)
  
  q_value <- p_value
  for (col in colnames(p_value)) {
    q_value[, col] <- p.adjust(p_value[, col], method = "fdr")
  }
  
  for (row in rownames(q_value)) {
    for (col in colnames(q_value)) {
      if (!is.na(q_value[row, col])) {
        RNAi_cellposition_qvalue[row, col] <- q_value[row, col]
      }
    }
  }
  
  q_value_binary <- matrix(0, nrow = nrow(q_value), ncol = ncol(q_value))
  rownames(q_value_binary) <- rownames(q_value)
  colnames(q_value_binary) <- colnames(q_value)
  
  for (i in 1:nrow(q_value)) {
    for (j in 1:ncol(q_value)) {
      if (!is.na(q_value[i, j]) && q_value[i, j] < 0.01) {
        q_value_binary[i, j] <- ifelse(z_score[i, j] >= 0, 1, -1)
      } else {
        q_value_binary[i, j] <- 0
      }
    }
  }
  
  for (row in rownames(q_value_binary)) {
    for (col in colnames(q_value_binary)) {
      if (!is.na(q_value_binary[row, col])) {
        RNAi_cellposition_qvalue_binary[row, col] <- q_value_binary[row, col]
      }
    }
  }
  
  cat("Stage", length(stage), "cells processing complete.\n")
}

# =============================================================================
# Step 4: Save Results
# =============================================================================

write.table(Embryo_dist_mean_WT, file = file.path(Division_asynchrony_phenotype_folder, "CTR_division_timing_vector_mean.txt"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(WT_single_cell_RMSD_all, file = file.path(Division_asynchrony_phenotype_folder, "CTR_division_asynchrony_RMSD_mean.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(RNAi_single_cell_RMSD_all, file = file.path(Division_asynchrony_phenotype_folder, "RNAi_division_asynchrony_RMSD.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(RNAi_cellposition_qvalue, file = file.path(Division_asynchrony_phenotype_folder, "RNAi_division_asynchrony_Q-value.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(RNAi_cellposition_qvalue_binary, file = file.path(Division_asynchrony_phenotype_folder, "RNAi_division_asynchrony_Q-value_binary_0.01.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(RNAi_cellposition_zscore, file = file.path(Division_asynchrony_phenotype_folder, "RNAi_division_asynchrony_Z-score.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


# =============================================================================
# Step 5: Visualization - Q-value Heatmap
# =============================================================================

setwd("~/Desktop/cbcn/draft_code/")

cedata <- read.table("~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/CTR_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t")
cbcndata <- read.table("~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/RNAi_division_asynchrony_RMSD.txt", header = TRUE, sep = "\t")
stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv") %>% select(1, 44)

data <- merge(cbcndata, cedata, by = "all_cells") %>%
  mutate(cb_mean = rowMeans(select(., corTim_cbp1:corTim_cbp6), na.rm = TRUE),
         cn_mean = rowMeans(select(., corTim_cnp1:corTim_cnp7), na.rm = TRUE)) %>%
  select(1, 17, 18, 15, 16)
names(data)[1] <- "cell"
data <- merge(data, stage, by = "cell")

# CB statistics
data$cb_z_score <- (data$cb_mean - data$cell_rmsd_mean) / data$cell_rmsd_sd
data$cb_p_value <- pnorm(abs(data$cb_z_score), lower.tail = FALSE)
data$cb_q_value <- p.adjust(data$cb_p_value, method = "fdr")
data$cb_q_value_binary <- ifelse(data$cb_q_value < 0.01 & data$cb_z_score >= 0, 1,
                                 ifelse(data$cb_q_value < 0.01 & data$cb_z_score < 0, -1, 0))

# CN statistics
data$cn_z_score <- (data$cn_mean - data$cell_rmsd_mean) / data$cell_rmsd_sd
data$cn_p_value <- pnorm(abs(data$cn_z_score), lower.tail = FALSE)
data$cn_q_value <- p.adjust(data$cn_p_value, method = "fdr")
data$cn_q_value_binary <- ifelse(data$cn_q_value < 0.01 & data$cn_z_score >= 0, 1,
                                 ifelse(data$cn_q_value < 0.01 & data$cn_z_score < 0, -1, 0))

# Prepare plot data
data <- data %>% select(1, 10, 14, 6)
names(data)[2] <- "cb_mean"
names(data)[3] <- "cn_mean"

data_plot <- data %>%
  mutate(
    stage = factor(stage, levels = c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32", "Stage_AB64", "Stage_AB128", "Stage_AB256")),
    cell_letter = substring(cell, 1, 1),
    cell_nchar = nchar(cell),
    cell_letter = factor(cell_letter, levels = c("A", "M", "E", "C", "D", "P"))
  ) %>%
  arrange(stage, cell_letter, cell_nchar)

data_plot$cell <- factor(data_plot$cell, levels = unique(data_plot$cell))

y_axis_order <- c("cn_mean", "cb_mean")
data_plot_long <- data_plot %>%
  pivot_longer(cols = 2:3, names_to = "variable", values_to = "value") %>%
  mutate(
    value = factor(value, levels = c("-1", "0", "1", NA)),
    variable = factor(variable, levels = y_axis_order)
  )

# Plot
ggplot(data_plot_long, aes(x = cell, y = variable, fill = value)) +
  geom_tile() +
  facet_grid(. ~ stage, scales = "free_x", space = "free") +
  scale_fill_manual(
    values = c("1" = "#9c2831", "0" = "#f5efe7", "-1" = "#024163"),
    na.value = "gray",
    name = "Value"
  ) +
  labs(x = "Cell", y = "Variable", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )

# ggsave("f3_sup7_qua_mean_0.01.pdf", last_plot(), width = 10, height = 2.5, dpi = 300)






# ============================================================
# Plot_significant cell counts by lineage & tissue + Venn diagram
# ============================================================
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(ggVennDiagram)
library(ggplotify)
library(grid)
library(VennDiagram)

setwd("~/Desktop/cbcn/draft_code/")

# ============================================================
# Load input tables
# ============================================================
cedata <- read.table(
  "~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/CTR_division_asynchrony_RMSD_mean.txt",
  header = TRUE, sep = "\t"
)

cbcndata <- read.table(
  "~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/RNAi_division_asynchrony_RMSD.txt",
  header = TRUE, sep = "\t"
)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv") %>%
  select(cell = 1, stage = 44)

tissue <- read.table(
  "~/Desktop/cbcn/linage/AllLineage.tsv",
  header = TRUE, sep = "\t"
) %>%
  select(cell = 1, CellFate = 3)

# ============================================================
# Compute CB / CN mean, z-score, p-value, q-value, and binary calls
# ============================================================
data <- merge(cbcndata, cedata, by = "all_cells") %>%
  as_tibble() %>%
  mutate(
    cb_mean = rowMeans(select(., corTim_cbp1:corTim_cbp6), na.rm = TRUE),
    cn_mean = rowMeans(select(., corTim_cnp1:corTim_cnp7), na.rm = TRUE)
  ) %>%
  select(
    cell = all_cells,
    cb_mean,
    cn_mean,
    cell_rmsd_mean,
    cell_rmsd_sd
  ) %>%
  merge(stage, by = "cell")

data <- data %>%
  mutate(
    cb_z_score = (cb_mean - cell_rmsd_mean) / cell_rmsd_sd,
    cb_p_value = pnorm(abs(cb_z_score), lower.tail = FALSE),
    cb_q_value = p.adjust(cb_p_value, method = "fdr"),
    cb_q_value_binary = case_when(
      cb_q_value < 0.01 & cb_z_score >= 0 ~ 1,
      cb_q_value < 0.01 & cb_z_score < 0 ~ -1,
      TRUE ~ 0
    ),
    cn_z_score = (cn_mean - cell_rmsd_mean) / cell_rmsd_sd,
    cn_p_value = pnorm(abs(cn_z_score), lower.tail = FALSE),
    cn_q_value = p.adjust(cn_p_value, method = "fdr"),
    cn_q_value_binary = case_when(
      cn_q_value < 0.01 & cn_z_score >= 0 ~ 1,
      cn_q_value < 0.01 & cn_z_score < 0 ~ -1,
      TRUE ~ 0
    )
  )

colnames(data)

# ============================================================
# Generic bar-plot function
# ============================================================
create_barplot <- function(count_data, x_var, y_limit, colors) {
  ggplot(count_data, aes(x = .data[[x_var]], y = Count, fill = channel)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = Count), vjust = -0.3, size = 3) +
    facet_wrap(~ channel) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0.8)) +
    labs(title = "", x = x_var, y = "Significant Number") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_blank(),
      legend.position = "none"
    )
}

# ============================================================
# Figure 1: counts by lineage (cell initial)
# ============================================================
data_lineage <- data %>%
  mutate(
    cell_letter = substring(cell, 1, 1),
    cell_letter = factor(cell_letter, levels = c("A", "M", "E", "C", "D", "P"))
  ) %>%
  select(cell, cell_letter, cb_q_value_binary, cn_q_value_binary)

data_long_lineage <- data_lineage %>%
  pivot_longer(
    cols = c(cb_q_value_binary, cn_q_value_binary),
    names_to = "channel",
    values_to = "value"
  ) %>%
  mutate(channel = recode(channel,
                          cb_q_value_binary = "cb_mean",
                          cn_q_value_binary = "cn_mean")) %>%
  filter(value == 1)

count_data_lineage <- data_long_lineage %>%
  count(channel, cell_letter, name = "Count")

p1 <- create_barplot(
  count_data_lineage,
  x_var = "cell_letter",
  y_limit = 195,
  colors = c("cb_mean" = "#0073C2", "cn_mean" = "#EFC000")
)
print(p1)
# ggsave("asychrony_signifi_cbcn_lineage_number.pdf", plot = p1, width = 7, height = 4)

# ============================================================
# Figure 2: counts by tissue (CellFate)
# ============================================================
data_tissue <- data %>%
  left_join(tissue, by = "cell") %>%
  filter(!is.na(CellFate)) %>%
  select(CellFate, cb_q_value_binary, cn_q_value_binary)

data_long_tissue <- data_tissue %>%
  pivot_longer(
    cols = c(cb_q_value_binary, cn_q_value_binary),
    names_to = "channel",
    values_to = "value"
  ) %>%
  mutate(channel = recode(channel,
                          cb_q_value_binary = "cb_mean",
                          cn_q_value_binary = "cn_mean")) %>%
  filter(value == 1)

count_data_tissue <- data_long_tissue %>%
  count(channel, CellFate, name = "Count")

p2 <- create_barplot(
  count_data_tissue,
  x_var = "CellFate",
  y_limit = 180,
  colors = c("cb_mean" = "#0073C2", "cn_mean" = "#EFC000")
)
print(p2)
# ggsave("asychrony_signifi_cbcn_tissue.pdf", plot = p2, width = 8, height = 5)

# ============================================================
# Venn diagram of significant CB vs CN cells
# ============================================================
cb_cells <- data$cell[data$cb_q_value_binary == 1]
cn_cells <- data$cell[data$cn_q_value_binary == 1]

n_cb <- length(cb_cells)
n_cn <- length(cn_cells)
n_overlap <- length(intersect(cb_cells, cn_cells))

cat("Number of cells with cb_mean = 1:", n_cb, "\n")
cat("Number of cells with cn_mean = 1:", n_cn, "\n")
cat("Number of cells with both = 1:", n_overlap, "\n")

grid::grid.newpage()
venn_plot <- draw.pairwise.venn(
  area1 = n_cb,
  area2 = n_cn,
  cross.area = n_overlap,
  category = c("cb_mean = 1", "cn_mean = 1"),
  fill = c("#0073C240", "#EFC00040"),
  col = c("#0073C2", "#EFC000"),
  cex = 1.5,
  cat.cex = 1.5
)
# ggsave("venn_plot.pdf", plot = venn_plot, width = 8, height = 8)











# =============================================================================
# ASYNCHRONY ANALYSIS_MEAN (inner group calculation first)
# =============================================================================
setwd("~/Desktop/cbcn/")
library(dplyr); library(readr); library(tidyr); library(scales)
library(stats); library(stringr); library(pracma)

# --------------------------- Shared configuration ----------------------------
data <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stages <- list(
  Stage_AB4   = c('ABal','ABar','ABpl','ABpr','MS','E','C'),
  Stage_AB8   = c('ABala','ABalp','ABara','ABarp','ABpla','ABplp','ABpra','ABprp','MSa','MSp','Ca','Cp','P3'),
  Stage_AB16  = c('ABalaa','ABalap','ABalpa','ABalpp','ABaraa','ABarap','ABarpa','ABarpp','ABplaa','ABplap',
                  'ABplpa','ABplpp','ABpraa','ABprap','ABprpa','ABprpp','MSaa','MSap','MSpa','MSpp','Ea','Ep','D'),
  Stage_AB32  = c('ABalaaa','ABalaap','ABalapa','ABalapp','ABalpaa','ABalpap','ABalppa','ABalppp','ABaraaa','ABaraap',
                  'ABarapa','ABarapp','ABarpaa','ABarpap','ABarppa','ABarppp','ABplaaa','ABplaap','ABplapa','ABplapp',
                  'ABplpaa','ABplpap','ABplppa','ABplppp','ABpraaa','ABpraap','ABprapa','ABprapp','ABprpaa','ABprpap',
                  'ABprppa','ABprppp','MSaaa','MSaap','MSapa','MSapp','MSpaa','MSpap','MSppa','MSppp','Eal','Ear',
                  'Epl','Epr','Caa','Cap','Cpa','Cpp','Da','Dp','P4'),
  Stage_AB64  = c('ABalaaaa','ABalaaap','ABalaapa','ABalaapp','ABalapaa','ABalapap','ABalappa','ABalappp','ABalpaaa',
                  'ABalpaap','ABalpapa','ABalpapp','ABalppaa','ABalppap','ABalpppa','ABalpppp','ABaraaaa','ABaraaap',
                  'ABaraapa','ABaraapp','ABarapaa','ABarapap','ABarappa','ABarappp','ABarpaaa','ABarpaap','ABarpapa',
                  'ABarpapp','ABarppaa','ABarppap','ABarpppa','ABarpppp','ABplaaaa','ABplaaap','ABplaapa','ABplaapp',
                  'ABplapaa','ABplapap','ABplappa','ABplappp','ABplpaaa','ABplpaap','ABplpapa','ABplpapp','ABplppaa',
                  'ABplppap','ABplpppa','ABplpppp','ABpraaaa','ABpraaap','ABpraapa','ABpraapp','ABprapaa','ABprapap',
                  'ABprappa','ABprappp','ABprpaaa','ABprpaap','ABprpapa','ABprpapp','ABprppaa','ABprppap','ABprpppa',
                  'ABprpppp','MSaaaa','MSaaap','MSaapa','MSaapp','MSapaa','MSapap','MSappa','MSappp','MSpaaa','MSpaap',
                  'MSpapa','MSpapp','MSppaa','MSppap','MSpppa','MSpppp','Caaa','Caap','Capa','Capp','Cpaa','Cpap','Cppa','Cppp'),
  Stage_AB128 = c('ABalaaaal','ABalaaaar','ABalaaapa','ABalaaapp','ABalaapaa','ABalaapap','ABalaappa','ABalaappp',
                  'ABalapaaa','ABalapaap','ABalapapa','ABalapapp','ABalappaa','ABalappap','ABalapppa','ABalapppp',
                  'ABalpaaaa','ABalpaaap','ABalpaapa','ABalpaapp','ABalpapaa','ABalpapap','ABalpappa','ABalpappp',
                  'ABalppaaa','ABalppaap','ABalppapa','ABalppapp','ABalpppaa','ABalpppap','ABalppppa','ABalppppp',
                  'ABaraaaaa','ABaraaaap','ABaraaapa','ABaraaapp','ABaraapaa','ABaraapap','ABaraappa','ABaraappp',
                  'ABarapaaa','ABarapaap','ABarapapa','ABarapapp','ABarappaa','ABarappap','ABarapppa','ABarapppp',
                  'ABarpaaaa','ABarpaaap','ABarpaapa','ABarpaapp','ABarpapaa','ABarpapap','ABarpappa','ABarpappp',
                  'ABarppaaa','ABarppaap','ABarppapa','ABarppapp','ABarpppaa','ABarpppap','ABarppppa','ABarppppp',
                  'ABplaaaaa','ABplaaaap','ABplaaapa','ABplaaapp','ABplaapaa','ABplaapap','ABplaappa','ABplaappp',
                  'ABplapaaa','ABplapaap','ABplapapa','ABplapapp','ABplappaa','ABplappap','ABplapppa','ABplapppp',
                  'ABplpaaaa','ABplpaaap','ABplpaapa','ABplpaapp','ABplpapaa','ABplpapap','ABplpappa','ABplpappp',
                  'ABplppaaa','ABplppaap','ABplppapa','ABplppapp','ABplpppaa','ABplpppap','ABplppppa','ABplppppp',
                  'ABpraaaaa','ABpraaaap','ABpraaapa','ABpraaapp','ABpraapaa','ABpraapap','ABpraappa','ABpraappp',
                  'ABprapaaa','ABprapaap','ABprapapa','ABprapapp','ABprappaa','ABprappap','ABprapppa','ABprapppp',
                  'ABprpaaaa','ABprpaaap','ABprpaapa','ABprpaapp','ABprpapaa','ABprpapap','ABprpappa','ABprpappp',
                  'ABprppaaa','ABprppaap','ABprppapa','ABprppapp','ABprpppaa','ABprpppap','ABprppppa','ABprppppp',
                  'MSaaaaa','MSaaaap','MSaaapa','MSaapaa','MSaapap','MSapaaa','MSapapa','MSapapp','MSpaaaa','MSpaaap',
                  'MSpaapa','MSpapaa','MSpapap','MSppaaa','MSppapa','MSppapp','Eala','Ealp','Eara','Earp','Epla','Eplp',
                  'Epra','Eprp','Caaaa','Caaap','Caapp','Capaa','Capap','Cappa','Cappp','Cpaaa','Cpaap','Cpapa','Cpapp',
                  'Cppaa','Cppap','Cpppa','Cpppp','Daa','Dap','Dpa','Dpp'),
  Stage_AB256 = c("ABalaaappr","ABalaapppa","ABalaapppp","ABalapaapp","ABalapappa","ABalappapp","ABalapppap",
                  "ABalappppa","ABalappppp","ABalpaapaa","ABalpaapap","ABalpaappa","ABalpaappp","ABalpapaaa","ABalpapaap",
                  "ABalpapapa","ABalpapapp","ABalpappaa","ABalpappap","ABalpapppp","ABalppapaa","ABalppapap","ABalppappa",
                  "ABalppappp","ABalpppapa","ABalpppapp","ABalppppaa","ABalppppap","ABalpppppa","ABalpppppp","ABaraaapaa",
                  "ABaraaapap","ABaraaappa","ABaraaappp","ABaraapaaa","ABaraapaap","ABaraapapa","ABaraapapp","ABaraappaa",
                  "ABaraappap","ABaraapppa","ABaraapppp","ABarapaapp","ABarapapap","ABarapappp","ABarappaap","ABarappapa",
                  "ABarappapp","ABarapppaa","ABarapppap","ABarappppa","ABarappppp","ABarpapaap","ABplaaaaap","ABplaapaap",
                  "ABplaapapa","ABplaapapp","ABplapaaaa","ABplapaaap","ABplapappp","ABplapppap","ABplpaaaaa","ABplpaaaap",
                  "ABplpaaapa","ABplpaaapp","ABplpaapaa","ABplpaapap","ABplpaappa","ABplpaappp","ABplpapaaa","ABplpapaap",
                  "ABplpapapa","ABplpapapp","ABplpappaa","ABplpapppa","ABplppaaaa","ABplppaapa","ABplppaapp","ABplppapaa",
                  "ABplppapap","ABplppappa","ABplppappp","ABplpppaaa","ABplpppaap","ABplpppapa","ABplppppaa","ABplppppap",
                  "ABplpppppa","ABplpppppp","ABpraaappp","ABpraapaap","ABpraapapp","ABprapaaap","ABprpaaaaa","ABprpaaaap",
                  "ABprpaaapa","ABprpaaapp","ABprpaapaa","ABprpaapap","ABprpaappa","ABprpaappp","ABprpapaaa","ABprpapaap",
                  "ABprpapapa","ABprpapapp","ABprpappaa","ABprpappap","ABprpapppa","ABprpapppp","ABprppaaaa","ABprppaapp",
                  "ABprppapaa","ABprppapap","ABprppappa","ABprppappp","ABprpppaaa","ABprpppaap","ABprpppapa","ABprppppaa",
                  "ABprppppap","ABprpppppa","ABprpppppp","Caapa","Capaaa","Capaap","Capapa","Capapp","Cappaa","Cappap",
                  "Capppa","Capppp","Cppaaa","Cppaap","Cppapa","Cppapp","Cpppaa","Cpppap","Cppppa","Cppppp","Daaa","Daap",
                  "Dapa","Dapp","Dpaa","Dpap","Dppa","Dppp","MSaaapp","MSaappa","MSaappp","MSapaap","MSappaa","MSappap",
                  "MSapppa","MSapppp","MSpappa","MSpappp","MSppaap","MSpppaa","MSpppap","MSppppa","MSppppp")
)
all_cells <- unique(unlist(stages))

# make sure output path exist
ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)

# input
write_wt_matrices <- function(celldivision_df, stages, all_cells, out_dir) {
  ensure_dir(out_dir)
  for (item in colnames(celldivision_df)[-1]) {
    emb <- celldivision_df[, c("cell", item)]
    division_time_dist <- matrix(NA_real_, length(all_cells), length(all_cells),
                                 dimnames = list(all_cells, all_cells))
    for (stage_cells in stages) {
      cell_value <- emb[emb$cell %in% stage_cells, ]
      if (nrow(cell_value) == length(stage_cells)) {
        stage_time_dist <- matrix(0, length(stage_cells), length(stage_cells),
                                  dimnames = list(stage_cells, stage_cells))
        for (cell in stage_cells) {
          value <- cell_value[cell_value$cell == cell, 2]
          values <- data.frame(cell = cell_value$cell, values = cell_value[, 2] - value)
          values_df <- data.frame(key = values[, 1], value = values[, 2])
          stage_time_dist <- cbind(stage_time_dist,
                                   values_df[order(rownames(stage_time_dist)), -1, drop = FALSE])
        }
        stage_time_dist <- stage_time_dist[, !colnames(stage_time_dist) %in% stage_cells, drop = FALSE]
        max_value <- max(stage_time_dist, na.rm = TRUE)
        stage_time_dist <- stage_time_dist / max_value
        colnames(stage_time_dist) <- rownames(stage_time_dist)
      } else {
        stage_time_dist <- matrix(NA_real_, length(stage_cells), length(stage_cells),
                                  dimnames = list(stage_cells, stage_cells))
      }
      idx <- !is.na(stage_time_dist)
      division_time_dist[rownames(stage_time_dist), colnames(stage_time_dist)][idx] <- stage_time_dist[idx]
    }
    write.table(
      division_time_dist,
      file = file.path(out_dir, paste0(item, "_time_Matrix_normlized.txt")),
      sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE
    )
  }
}

# asychrony analysis
compute_wt_asynchrony <- function(timing_dir, stages, all_cells, out_dir, prefix) {
  ensure_dir(out_dir)
  Embryo_dist_mean_WT <- matrix(NA_real_, length(all_cells), length(all_cells),
                                dimnames = list(all_cells, all_cells))
  WT_single_cell_RMSD_all <- data.frame(
    all_cells = all_cells,
    cell_rmsd_mean = NA_real_,
    cell_rmsd_sd = NA_real_,
    cell_rmsd_n = NA_real_,
    row.names = all_cells
  )
  file_list <- list.files(timing_dir)
  for (stage_cells in stages) {
    Total_matrix <- NULL
    for (file in file_list) {
      Embryo_dist <- read.table(file.path(timing_dir, file), header = TRUE, sep = "\t", row.names = 1)
      Embryo_dist <- Embryo_dist[stage_cells, stage_cells]
      Total_matrix <- if (is.null(Total_matrix)) Embryo_dist else cbind(Total_matrix, Embryo_dist)
    }
    cell_rmsd_mean <- cell_rmsd_sd <- cell_rmsd_n <- list()
    Embryo_dist_mean <- data.frame()
    for (cell in stage_cells) {
      matching_columns <- grep(cell, colnames(Total_matrix))
      cell_vector_matrix <- Total_matrix[, matching_columns, drop = FALSE]
      cell_vector_mean <- rowMeans(cell_vector_matrix, na.rm = TRUE)
      cell_vector_matrix <- cell_vector_matrix[-which(rownames(cell_vector_matrix) == cell), , drop = FALSE]
      cell_rmsd <- dist(t(cell_vector_matrix)) / sqrt(ncol(t(cell_vector_matrix)))
      cell_rmsd <- cell_rmsd[!is.na(cell_rmsd)]
      cell_rmsd_mean[[cell]] <- mean(cell_rmsd)
      cell_rmsd_sd[[cell]] <- sd(cell_rmsd)
      cell_rmsd_n[[cell]] <- length(cell_rmsd)
      Embryo_dist_mean <- if (ncol(Embryo_dist_mean) == 0) data.frame(cell_vector_mean) else
        cbind(Embryo_dist_mean, cell_vector_mean)
    }
    Embryo_dist_mean <- as.matrix(Embryo_dist_mean)
    colnames(Embryo_dist_mean) <- rownames(Embryo_dist_mean)
    idx <- !is.na(Embryo_dist_mean)
    Embryo_dist_mean_WT[rownames(Embryo_dist_mean), colnames(Embryo_dist_mean)][idx] <- Embryo_dist_mean[idx]
    WT_single_cell_RMSD_all[names(cell_rmsd_mean), "cell_rmsd_mean"] <- unlist(cell_rmsd_mean)
    WT_single_cell_RMSD_all[names(cell_rmsd_sd),   "cell_rmsd_sd"]   <- unlist(cell_rmsd_sd)
    WT_single_cell_RMSD_all[names(cell_rmsd_n),    "cell_rmsd_n"]    <- unlist(cell_rmsd_n)
    cat("Stage", paste(head(stage_cells, 3), collapse = "/"), "... processed.\n")
  }
  write.table(
    Embryo_dist_mean_WT,
    file = file.path(out_dir, paste0(prefix, "_division_timing_vector_mean.txt")),
    sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE
  )
  write.table(
    WT_single_cell_RMSD_all,
    file = file.path(out_dir, paste0(prefix, "_division_asynchrony_RMSD_mean.txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# --------------------------- Pipeline for ce/cb/cn ---------------------------
jobs <- list(
  ce = list(cols = 23:30, timing_dir = "~/Desktop/cbcn/output_folder/test/ce_division_timing_matrix",
            out_dir = "~/Desktop/cbcn/output_folder/test/cedivision_asynchrony_phenotype"),
  cb = list(cols = 31:36, timing_dir = "~/Desktop/cbcn/output_folder/test/cb_division_timing_matrix",
            out_dir = "~/Desktop/cbcn/output_folder/test/cbdivision_asynchrony_phenotype"),
  cn = list(cols = 37:43, timing_dir = "~/Desktop/cbcn/output_folder/test/cn_division_timing_matrix",
            out_dir = "~/Desktop/cbcn/output_folder/test/cndivision_asynchrony_phenotype")
)

for (prefix in names(jobs)) {
  cat("==== Processing", prefix, "====\n")
  job <- jobs[[prefix]]
  celldivision <- data[, c(1, job$cols)]
  colnames(celldivision)[1] <- "cell"
  write_wt_matrices(celldivision, stages, all_cells, job$timing_dir)
  compute_wt_asynchrony(job$timing_dir, stages, all_cells, job$out_dir, prefix)
}













# 设置工作目录
setwd("~/Desktop/cbcn/draft_code/")

library(dplyr);library(readr);library(tidyr);library(scales);library(stats);library(stringr);library(pracma)
library(ggplot2)
library(ggpattern)

# 读取数据
ce <- read.table("~/Desktop/cbcn/output_folder/test/cedivision_asynchrony_phenotype/ce_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t")
cb <- read.table("~/Desktop/cbcn/output_folder/test/cbdivision_asynchrony_phenotype/cb_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t")
cn <- read.table("~/Desktop/cbcn/output_folder/test/cndivision_asynchrony_phenotype/cn_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t")

# 查看各数据框的列名（可选）
colnames(cn)

ce <- ce %>% 
  rename(ce_mean_RMSD = cell_rmsd_mean)
cb <- cb %>% 
  rename(cb_mean_RMSD = cell_rmsd_mean)
cn <- cn %>% 
  rename(cn_mean_RMSD = cell_rmsd_mean)
colnames(ce)

all_data <- ce %>% 
  left_join(cb, by = "all_cells") %>%
  left_join(cn, by = "all_cells")
names(all_data)[1] <- "cell"

colnames(all_data)


stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
stage <- stage %>%
  select(1,44)
data <- merge(all_data, stage, by = "cell")

colnames(data)
data <- data %>%
  select(1,2,5,8,11)
colnames(data)

tissue <- read.table("~/Desktop/cbcn/linage/AllLineage.tsv", header=TRUE, sep="\t") %>%
  select(parent=1, CellFate=3)
colnames(tissue)
names(tissue)[1] <- "cell"
data <- merge(data, tissue, by = "cell")
colnames(data)


# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

# -------------------------------------------
# 假设数据已存储在 data 数据框中，
# 列名包括："cell", "ce_mean_RMSD", "cb_mean_RMSD", "cn_mean_RMSD", "stage"
# -------------------------------------------

# 1. 将 stage 列转换为因子，并按指定顺序排序
order_stage <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                 "Stage_AB64", "Stage_AB128", "Stage_AB256")
data$stage <- factor(data$stage, levels = order_stage)

# 2. 将数据从宽格式转换为长格式
#    将 ce_mean_RMSD、cb_mean_RMSD、cn_mean_RMSD 三列整合到一起
data_long <- pivot_longer(
  data,
  cols = c("ce_mean_RMSD", "cb_mean_RMSD", "cn_mean_RMSD"),
  names_to = "Metric",
  values_to = "RMSD_value"
)

# 3. 根据要求生成辅助变量并对数据排序：
#    - 强制 stage 按指定顺序排序
#    - 提取 cell 的首字母，并转换为因子（指定顺序）
#    - 计算 cell 名称的长度
#    - 按照 stage、首字母及 cell 排序，再将排序后的 cell 转为因子，确保 x 轴显示顺序正确
order_letter <- c("A", "M", "E", "C", "D", "P")
data_long <- data_long %>%
  mutate(
    stage = factor(stage, levels = order_stage),
    first_letter = factor(substr(cell, 1, 1), levels = order_letter),
    cell_length = nchar(cell)
  ) %>%
  arrange(stage, first_letter, cell,CellFate)
data_long$cell <- factor(data_long$cell, levels = unique(data_long$cell))

# 4. 计算每个 cell 的新 x 坐标，新坐标在每个 stage 内连续，并在不同 stage 之间增加间隔
gap_size <- 2   # 定义每个 stage 之间增加的间隔大小
# 提取各 cell 信息（cell及其所属 stage），并计算 stage 内序号
cell_order <- data_long %>% 
  distinct(cell, stage) %>%
  group_by(stage) %>%
  arrange(cell) %>%
  mutate(cell_index = row_number()) %>%
  ungroup()

# 4.1 计算每个 stage 的 cell 数量及前序累计量
stage_counts <- cell_order %>% 
  group_by(stage) %>% 
  summarize(count = n()) %>%
  arrange(factor(stage, levels = order_stage)) %>%
  mutate(cum = lag(cumsum(count), default = 0))
cell_order <- left_join(cell_order, stage_counts, by = "stage")

# 4.2 计算新的 x 坐标：
#  新坐标 = 累计的 cell 数量(来自之前 stage) + 当前 stage 内序号 + 间隔总量( gap_size*(stage序号-1) )
cell_order <- cell_order %>%
  mutate(new_x = cum + cell_index + gap_size * (as.integer(stage) - 1))

# 将新 x 坐标合并回 data_long
data_long <- left_join(data_long, cell_order %>% select(cell, new_x), by = "cell")

# 5. 根据新的 x 坐标，计算指定 stage 的背景区域（灰色背景，仅针对部分阶段）
# rect_data <- cell_order %>%
#   filter(stage %in% c("Stage_AB8", "Stage_AB32", "Stage_AB128")) %>%
#   group_by(stage) %>%
#   summarize(xmin = min(new_x) - 0.5,
#             xmax = max(new_x) + 0.5) %>%
#   ungroup()

# 6. 对每个 cell 计算三个指标的最小值、中间值和最大值及对应的指标名称，
#    同时记录该 cell 在 x 轴上的新坐标（new_x）
fill_df <- data_long %>%
  group_by(cell) %>%
  arrange(RMSD_value) %>%  # 升序，第一行为最小值，第二行为中间值，第三行为最大值
  summarize(
    min_value  = first(RMSD_value),
    Metric_min = first(Metric),
    mid_value  = nth(RMSD_value, 2),
    Metric_mid = nth(Metric, 2),
    max_value  = last(RMSD_value),
    Metric_max = last(Metric),
    x = first(new_x)
  ) %>%
  ungroup()

# 7. 定义指标对应的颜色
metric_colors <- c("ce_mean_RMSD" = "#CD534C", 
                   "cb_mean_RMSD" = "#0073C2", 
                   "cn_mean_RMSD" = "#EFC000")
# 将每个分段对应的颜色加入 fill_df
fill_df <- fill_df %>%
  mutate(
    fill_min = metric_colors[Metric_min],
    fill_mid = metric_colors[Metric_mid],
    fill_max = metric_colors[Metric_max]
  )

# 8. 根据 cell_order 生成每个 stage 的边框数据
stage_border <- cell_order %>%
  group_by(stage) %>%
  summarize(
    xmin = min(new_x) - 0.5,
    xmax = max(new_x) + 0.5
  ) %>%
  ungroup()

# 9. 绘图
p_main <-ggplot(data_long, aes(x = new_x, y = RMSD_value, 
                               group = interaction(Metric, stage), color = Metric)) +
  # # 添加部分 stage 的灰色背景
  # geom_rect(data = rect_data,
  #           aes(xmin = xmin, xmax = xmax),
  #           ymin = -Inf, ymax = Inf,
  #           fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # 添加每个 stage 的黑色边框
  geom_rect(data = stage_border,
            aes(xmin = xmin, xmax = xmax),
            ymin = -Inf, ymax = Inf,
            fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
  # 第一步：每个 cell下方区域从 y = 0 到最小值，填充对应最小值指标的颜色
  geom_rect(data = fill_df,
            aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 0, ymax = min_value, fill = fill_min),
            inherit.aes = FALSE, alpha = 0.4) +
  # 第二步：每个 cell 中间区域从最小值到中间值，填充对应中间值指标的颜色
  geom_rect(data = fill_df,
            aes(xmin = x - 0.5, xmax = x + 0.5, ymin = min_value, ymax = mid_value, fill = fill_mid),
            inherit.aes = FALSE, alpha = 0.4) +
  # 第三步：每个 cell 上方区域从中间值到最大值，填充对应最大值指标的颜色
  geom_rect(data = fill_df,
            aes(xmin = x - 0.5, xmax = x + 0.5, ymin = mid_value, ymax = max_value, fill = fill_max),
            inherit.aes = FALSE, alpha = 0.4) +
  # 绘制折线（和顺序：先填充再绘制折线）
  geom_line(size = 0.4) +
  # 设置折线颜色
  scale_color_manual(values = metric_colors) +
  # 利用 scale_fill_identity() 根据数据中的颜色值直接设置填充颜色
  scale_fill_identity() +
  labs(
    title = "",
    x = "Cell",
    y = "Mean RMSD"
  ) +
  # 设置 x 轴刻度使用 cell 对应的新 x 坐标，并显示 cell 标签
  scale_x_continuous(
    expand = c(0, 0),
    breaks = cell_order$new_x,
    labels = cell_order$cell
  ) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_blank(),  # 不显示 x 轴文本标签
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
print(p_main)
#ggsave("asychrony_compare_line.pdf",p_main, width = 15, height = 4)


# ---------------------
# 新增：生成 treemap 数据及图形（使 treemap 每个 stage 的宽度和下面主图一致）
library(treemapify)
library(patchwork)
library(purrr)  # 用于 map 和 split 操作

# 计算 treemap 数据：统计每个 cell 中最大值对应指标（Metric_max）在各 stage 下的数量
tree_data <- fill_df %>% 
  left_join(cell_order %>% select(cell, stage), by = "cell") %>%
  group_by(stage, Metric_max) %>%
  summarize(count = n(), .groups = "drop")

# 计算每个stage宽度
stage_border <- cell_order %>%
  group_by(stage) %>%
  summarize(
    xmin = min(new_x) - 0.5,
    xmax = max(new_x) + 0.5
  ) %>%
  ungroup() %>%
  mutate(width = xmax - xmin)

stage_widths <- stage_border$width
names(stage_widths) <- stage_border$stage

# 拆分treemap数据，分别绘制
library(purrr)
treemap_list <- map(order_stage, function(stg) {
  df_sub <- tree_data %>% filter(stage == stg)
  ggplot(df_sub, aes(area = count, fill = Metric_max, label = paste(count))) +
    #geom_treemap() +
    geom_treemap(alpha = 0.6) +  # 设置透明度为40%
    geom_treemap_text(color = "white", place = "centre", grow = FALSE, size = 10, min.size = 0) +
    scale_fill_manual(values = metric_colors) +
    labs(title = stg) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),  # 去掉标题
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(1, 1, 1, 2.5)  # 去除图形边距
    )
})

# 横向拼接treemap，按stage宽度比例
library(patchwork)
p_tree_combined <- wrap_plots(treemap_list, nrow = 1, widths = stage_widths)

# 竖向拼接treemap和主图
p_3g <- p_tree_combined / p_main + plot_layout(heights = c(0.5, 3))

print(p_3g)
#ggsave("f3_3.pdf", p_3g, width = 15, height = 5)




# =========================
# Boxplot 需要的预处理
# =========================
data_long2 <- data_long %>%
  mutate(
    Metric = factor(Metric, levels = c("ce_mean_RMSD", "cb_mean_RMSD", "cn_mean_RMSD")),
    first_letter = factor(substr(cell, 1, 1), levels = order_letter)
  ) %>%
  filter(!is.na(RMSD_value), !is.na(Metric))

# （可选）更好看的指标名字
metric_labels <- c(
  ce_mean_RMSD = "cel",
  cb_mean_RMSD = "cbr",
  cn_mean_RMSD = "cni"
)

# =========================
# 1) 所有 cell 汇总：一张 box plot（按 Metric 分组）
# =========================
p_box_all <- ggplot(data_long2, aes(x = Metric, y = RMSD_value, fill = Metric)) +
  stat_boxplot(geom="errorbar", width=0.42, size=0.5) +
  geom_boxplot(width=0.65, outlier.shape=NA, color="black", size=0.5) +
  scale_fill_manual(values = metric_colors, labels = metric_labels) +
  scale_x_discrete(labels = metric_labels) +
  labs(x = NULL, y = "Mean RMSD") +
  theme_classic(base_size=14) +
  theme(
    panel.grid = element_blank(),
    axis.line=element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )+
  scale_y_continuous(limits=c(0,0.3), expand=c(0,0))

print(p_box_all)
#ggsave("f2_1.pdf", plot=p_box_all, width=5, height=6, dpi=300)
 
# library(ggpubr)
# 
# comparisons <- list(
#   c("ce_mean_RMSD","cb_mean_RMSD"),
#   c("ce_mean_RMSD","cn_mean_RMSD"),
#   c("cb_mean_RMSD","cn_mean_RMSD")
# )
# 
# p_box_all_sig <- ggplot(data_long2, aes(x = Metric, y = RMSD_value, fill = Metric)) +
#   stat_boxplot(geom="errorbar", width=0.42, linewidth=0.5) +
#   geom_boxplot(width=0.65, outlier.shape=NA, color="black", linewidth=0.5) +
#   scale_fill_manual(values = metric_colors, labels = metric_labels) +
#   scale_x_discrete(labels = metric_labels) +
#   labs(x = NULL, y = "Mean RMSD") +
#   theme_classic(base_size=14) +
#   theme(
#     panel.grid = element_blank(),
#     axis.line = element_blank(),
#     legend.position = "none",
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
#   ) +
#   scale_y_continuous(expand = c(0, 0)) +
#   coord_cartesian(ylim = c(0, 0.3)) +   # 关键：用它替代 limits
#   stat_compare_means(
#     comparisons = comparisons,
#     method = "wilcox.test",
#     p.adjust.method = "BH",
#     label = "p.signif",
#     label.y = c(0.22, 0.245, 0.27),     # 往下调：自行微调这三个数
#     tip.length = 0.01
#   )
# 
# print(p_box_all_sig)

# =========================
# 2) 按 first_letter 分组：一张 box plot
#    展示方式：x=first_letter，每个 letter 内用 Metric 分组（并排箱线图）
# =========================
# 确保 first_letter 和 Metric 都是 factor（很关键）
data_long2 <- data_long2 %>%
  mutate(
    first_letter = factor(first_letter, levels = c("A","M","E","C","D","P")),
    Metric = factor(Metric, levels = c("ce_mean_RMSD","cb_mean_RMSD","cn_mean_RMSD"))
  )

p_box_letter <- ggplot(
  data_long2,
  aes(x = first_letter, y = RMSD_value, fill = Metric)
) +
  stat_boxplot(
    aes(group = interaction(first_letter, Metric)),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.42,
    linewidth = 0.5
  ) +
  geom_boxplot(
    aes(group = interaction(first_letter, Metric)),
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = metric_colors, labels = metric_labels) +
  labs(x = "First letter", y = "Mean RMSD", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0))

print(p_box_letter)
ggsave("f2_2.pdf", plot=p_box_letter, width=15, height=6, dpi=300)


# =========================
# 3) 按 CellFate 分组：一张 box plot
#    注意：CellFate 类别可能很多，建议翻转坐标让标签可读
# =========================
data_long2 <- data_long2 %>%
  mutate(
    Metric   = factor(Metric, levels = c("ce_mean_RMSD","cb_mean_RMSD","cn_mean_RMSD")),
    CellFate = factor(CellFate)  # 如果你有想要的顺序，也可以在这里指定 levels
  )

p_box_fate <- ggplot(
  data_long2,
  aes(x = CellFate, y = RMSD_value, fill = Metric)
) +
  stat_boxplot(
    aes(group = interaction(CellFate, Metric)),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.42,
    linewidth = 0.5
  ) +
  geom_boxplot(
    aes(group = interaction(CellFate, Metric)),
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = metric_colors, labels = metric_labels) +
  labs(x = "Cell fate", y = "Mean RMSD", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1)  # fate 多的话更清楚
  ) +
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0))

print(p_box_fate)
ggsave("f2_3.pdf", plot=p_box_fate, width=15, height=6, dpi=300)








# 设置工作目录
setwd("~/Desktop/cbcn/draft_code/fig5_all/")

library(dplyr);library(readr);library(tidyr);library(scales);library(stats);library(stringr);library(pracma)
library(ggplot2)
library(ggpubr)

# 1. 读取数据
data_cb <- read.csv("~/Desktop/cbcn/draft_code/fig5_all/code/analysis/cb_mean_asychrony.csv")
data_cn <- read.csv("~/Desktop/cbcn/draft_code/fig5_all/code/analysis/cn_mean_asychrony.csv")


# cedata <- read.table("~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/CTR_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t")
# cbcndata <- read.table("~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/RNAi_division_asynchrony_RMSD.txt", header = TRUE, sep = "\t")
# stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
# colnames(stage)
# stage <- stage %>%
#   select(1,44)
# colnames(stage)
# colnames(cedata)
# colnames(cbcndata)
# data <- merge(cbcndata, cedata, by = "all_cells")
# colnames(data)
# data <- data %>%
#   mutate(cb_mean = rowMeans(select(., corTim_cbp1:corTim_cbp6), na.rm = TRUE),
#          cn_mean = rowMeans(select(., corTim_cnp1:corTim_cnp7), na.rm = TRUE))%>%
#   select(1,17,18,15,16)
# names(data)[1] <- "cell"
# data <- merge(data, stage, by = "cell")
# colnames(data)
# 
# 
# # -----------------------------
# # 1. 计算 CB 部分的 z-score, p-value, q-value, q_value_binary
# # z-score = (cb_mean - cell_rmsd_mean) / cell_rmsd_sd
# data$cb_z_score <- (data$cb_mean - data$cell_rmsd_mean) / data$cell_rmsd_sd
# 
# # 计算双侧 p-value，注意使用绝对值 z-score
# data$cb_p_value <- pnorm(abs(data$cb_z_score), lower.tail = FALSE)
# 
# # FDR 校正得到 q-value
# data$cb_q_value <- p.adjust(data$cb_p_value, method = "fdr")
# 
# # 将 q_value 二值化：若 q < 0.01 且 z_score >= 0 则赋值为 1；若 q < 0.01 且 z_score < 0 则赋值为 -1；其他情况赋值为 0
# data$cb_q_value_binary <- ifelse(data$cb_q_value < 0.01 & data$cb_z_score >= 0, 1,
#                                  ifelse(data$cb_q_value < 0.01 & data$cb_z_score < 0, -1, 0))
# 
# # -----------------------------
# # 2. 计算 CN 部分的 z-score, p-value, q-value, q_value_binary
# # z-score = (cn_mean - cell_rmsd_mean) / cell_rmsd_sd
# data$cn_z_score <- (data$cn_mean - data$cell_rmsd_mean) / data$cell_rmsd_sd
# 
# # 计算双侧 p-value
# data$cn_p_value <- pnorm(abs(data$cn_z_score), lower.tail = FALSE)
# 
# # FDR 校正得到 q-value
# data$cn_q_value <- p.adjust(data$cn_p_value, method = "fdr")
# 
# # 二值化 q_value：若 q < 0.01 且 z_score >= 0 则赋值为 1; 若 q < 0.01 且 z_score < 0 则赋值为 -1; 否则赋值为 0
# data$cn_q_value_binary <- ifelse(data$cn_q_value < 0.01 & data$cn_z_score >= 0, 1,
#                                  ifelse(data$cn_q_value < 0.01 & data$cn_z_score < 0, -1, 0))
# 
# # -----------------------------
# # 3. 将所有结果组合到一个数据框中
# # 提取需要的列：cell, control 的均值与 sd, CB 以及 CN 的测量均值及计算指标，stage 信息
# final_results <- data[, c("cell", 
#                           "cell_rmsd_mean", "cell_rmsd_sd", 
#                           "cb_mean", "cb_z_score", "cb_p_value", "cb_q_value", "cb_q_value_binary",
#                           "cn_mean", "cn_z_score", "cn_p_value", "cn_q_value", "cn_q_value_binary",
#                           "stage")]
# 
# # 查看最终结果
# print(head(final_results))
# # 假设你要保存的数据框名为 DataFrameName，
# # 将数据保存到当前工作目录下的 "DataFrameName.csv" 文件中（不保存行名）
# #write.csv(final_results, file = "cbcn_comparece_mean_inter.csv", row.names = FALSE)




# ----------------------------------------------------------
# 1. Load input data
# ----------------------------------------------------------
data_cb <- read.table("~/Desktop/cbcn/output_folder/test/cbdivision_asynchrony_phenotype/cb_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t")
data_cn <- read.table("~/Desktop/cbcn/output_folder/test/cndivision_asynchrony_phenotype/cn_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t")
cbcn    <- read.csv("~/Desktop/cbcn/output_folder/0114/cbcn_comparece_mean_inter.csv")


# ----------------------------------------------------------
# 2. Merge CE reference with CB / CN
# ----------------------------------------------------------
# Extract CE values for merging
cecn <- cbcn %>% select(cell, cn_mean, stage)
cecb <- cbcn %>% select(cell, cb_mean, stage)

merged_cn <- merge(data_cn, cecn, by.x = "all_cells", by.y = "cell")
merged_cb <- merge(data_cb, cecb, by.x = "all_cells", by.y = "cell")

# ----------------------------------------------------------
# 3. Add tissue and lineage annotation
# ----------------------------------------------------------
tissue <- read.table("~/Desktop/cbcn/linage/AllLineage.tsv", header = TRUE, sep = "\t") %>%
  select(all_cells = 1, CellFate = 3)

merged_cb <- merged_cb %>%
  merge(tissue, by = "all_cells") %>%
  mutate(lineage = if_else(substr(all_cells, 1, 1) == "Z", "P", substr(all_cells, 1, 1)))

merged_cn <- merged_cn %>%
  merge(tissue, by = "all_cells") %>%
  mutate(lineage = if_else(substr(all_cells, 1, 1) == "Z", "P", substr(all_cells, 1, 1)))

# ----------------------------------------------------------
# 4. Define scatter plot function
# ----------------------------------------------------------
plot_scatter <- function(df, intra_col, point_color) {
  ggplot(df, aes_string(x = intra_col, y = "cell_rmsd_mean")) +
    geom_point(color = point_color, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1, alpha = 0.6) +
    stat_cor(method = "pearson", label.x = 0.1, label.y = 0.30) +
    labs(x = "Intra (CB or CN)", y = "Inter (CE)") +
    # facet_wrap(~ factor(stage, levels = c("Stage_AB4", "Stage_AB8", "Stage_AB16",
    #                                       "Stage_AB32", "Stage_AB64", "Stage_AB128", "Stage_AB256")),
    #            nrow = 1)+
    #facet_wrap(~ factor(CellFate),nrow = 1)+
    scale_x_continuous(limits = c(0, 0.8)) +
    scale_y_continuous(limits = c(0, 0.33)) +
    theme_minimal() +
    theme(
      axis.text.y  = element_text(size = 10),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_blank(),
      axis.ticks = element_line(color = "black")
    )
}

# ----------------------------------------------------------
# 5. Plot CB vs CE and CN vs CE
# ----------------------------------------------------------
p_cb <- plot_scatter(merged_cb, "cb_mean", "#0073C2")
p_cn <- plot_scatter(merged_cn, "cn_mean", "#EFC000")

# ----------------------------------------------------------
# 6. Combine plots and save
# ----------------------------------------------------------
combined_plot <- p_cb + p_cn + plot_layout(ncol = 1)

print(combined_plot)

# ggsave("inter_intra_all.pdf", plot = combined_plot, width = 4, height = 8)
# ggsave("inter_intra_stage.pdf", plot = combined_plot, width = 10, height = 4.6)
# ggsave("inter_intra_tissue.pdf", plot = combined_plot, width = 10, height = 4.6)
# ggsave("inter_intra_lineage.pdf", plot = combined_plot, width = 10, height = 4.6)





setwd("~/Desktop/cbcn/draft_code/")
library(ggplot2);library(dplyr);library(tidyverse);library(ggsci);library(ggpubr);library(reshape2);library(RColorBrewer)

################################################################################
############################### mean value ##################################### 
cedata <- read.table("~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/CTR_division_asynchrony_RMSD_mean.txt", header = TRUE, sep = "\t")
cbcndata <- read.table("~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/RNAi_division_asynchrony_RMSD.txt", header = TRUE, sep = "\t")
stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv")
colnames(stage)
stage <- stage %>%
  select(1,44)
colnames(stage)
colnames(cedata)
colnames(cbcndata)
data <- merge(cbcndata, cedata, by = "all_cells")
colnames(data)
data <- data %>%
  mutate(cb_mean = rowMeans(select(., corTim_cbp1:corTim_cbp6), na.rm = TRUE),
         cn_mean = rowMeans(select(., corTim_cnp1:corTim_cnp7), na.rm = TRUE))%>%
  select(1,17,18,15,16)
names(data)[1] <- "cell"
data <- merge(data, stage, by = "cell")
colnames(data)
names(data)[4] <- "ce_mean"
data <- data %>%
  select(1,4,2,3,6)
colnames(data)
library(dplyr)

data <- data %>%
  mutate(
    cb_mean = cb_mean - ce_mean,
    cn_mean = cn_mean - ce_mean
  )
colnames(data)
data <- data %>%
  select(1,3,4,5)
colnames(data)
library(dplyr)

data <- data %>%
  mutate(stage = ifelse(stage %in% c("Stage_AB4", "Stage_AB8"), "Stage_AB16", stage))

# 查看修改后的数据
head(data)

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)

# 指定 stage 的顺序
stage_levels <- c( "Stage_AB16", 
                   "Stage_AB32", "Stage_AB64", "Stage_AB128", "Stage_AB256")

# 将 stage 列转换为因子并指定级别顺序
data <- data %>%
  mutate(stage = factor(stage, levels = stage_levels))

# 将 ce_mean, cb_mean, cn_mean 三个变量转换成长格式，同时指定因子顺序（确保顺序：ce, cb, cn）
data_long <- data %>%
  pivot_longer(cols = c(cb_mean, cn_mean), 
               names_to = "measurement", 
               values_to = "value") %>%
  mutate(measurement = factor(measurement, levels = c("cb_mean", "cn_mean")))


# 计算散点在每个 (measurement, stage) 分组内按 value 排序后的新 x 坐标
# 箱线图的箱宽设置为 0.65，因此散点在该范围内均匀分布
data_long2 <- data_long %>%
  group_by(measurement, stage) %>%
  arrange(value, .by_group = TRUE) %>%
  mutate(
    x_offset = seq(-0.65/2, 0.65/2, length.out = n()),
    x_new = as.numeric(stage) + x_offset
  ) %>%
  ungroup()

p <- ggplot() +
  # 绘制 errorbar（须线），这里显式设定 group = stage
  stat_boxplot(
    data = subset(data_long, measurement == "cb_mean"),
    aes(x = as.numeric(stage), y = value, group = stage),
    geom = "errorbar", width = 0.42, size = 0.5,
    color = "black"
  ) +
  stat_boxplot(
    data = subset(data_long, measurement == "cn_mean"),
    aes(x = as.numeric(stage), y = value, group = stage),
    geom = "errorbar", width = 0.42, size = 0.5,
    color = "black"
  ) +
  # 绘制箱线图（箱体填充为白色、无透明度），添加 group = stage 保证各组数据合并
  geom_boxplot(
    data = subset(data_long, measurement == "cb_mean"),
    aes(x = as.numeric(stage), y = value, group = stage),
    width = 0.65, outlier.shape = NA, fill = "white", alpha = 1,
    color = "black", size = 0.5
  ) +
  geom_boxplot(
    data = subset(data_long, measurement == "cn_mean"),
    aes(x = as.numeric(stage), y = value, group = stage),
    width = 0.65, outlier.shape = NA, fill = "white", alpha = 1,
    color = "black", size = 0.5
  ) +
  # 绘制散点，使用预先计算好的 x_new 保证散点在对应箱宽范围中且按 value 顺序排列
  geom_point(
    data = data_long2, 
    aes(x = x_new, y = value, color = stage),
    size = 1.5, alpha = 0.4
  ) +
  # 按 measurement 分面，每个面板显示各自的 y 轴范围
  facet_wrap(~ measurement, ncol = 2, scales = "free_y") +
  scale_x_continuous(
    breaks = unique(as.numeric(data_long2$stage)),
    labels = levels(data_long2$stage),
    expand = c(0.1, 0.1)
  ) +
  # 限定 y 轴显示范围
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.1, 0.2)) +
  theme_bw() +
  labs(
    title = "Boxplot of cb_mean and cn_mean by Stage",
    x = "Stage", 
    y = "Value"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    # 去除分面小标题的边框背景
    strip.background = element_blank()
  ) +
  # 使用 Set1 调色板确保同一 stage 的散点颜色一致
  scale_color_brewer(palette = "Set1")

# 显示图形
print(p)


# # 绘制热图
# ggsave(
#   filename = "copyscience.pdf",  # 文件名，你也可以给出绝对路径
#   plot = last_plot(),                  # 要保存的图形对象
#   width = 10,                            # 图形宽度（单位默认为英寸）
#   height = 6,                            # 图形高度
#   dpi = 300                              # 分辨率
# )






# ============================================================
# Setup
# ============================================================
setwd("~/Desktop/cbcn/draft_code/new/")

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(ggVennDiagram)
library(ggplotify)
library(grid)
library(ggbreak)

# ============================================================
# Load input tables
# ============================================================
cedata <- read.table(
  "~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/CTR_division_asynchrony_RMSD_mean.txt",
  header = TRUE,
  sep = "\t"
)

cbcndata <- read.table(
  "~/Desktop/cbcn/output_folder/0114/division_asynchrony_phenotype/RNAi_division_asynchrony_RMSD.txt",
  header = TRUE,
  sep = "\t"
)

stage <- read.csv("~/Desktop/cbcn/nor_Length_Time_correct.csv") %>%
  select(cell = 1, stage = 44)

tissue <- read.table(
  "~/Desktop/cbcn/linage/AllLineage.tsv",
  header = TRUE,
  sep = "\t"
) %>%
  select(cell = 1, CellFate = 3)

# ============================================================
# Compute CB / CN mean, z-score, p-value, q-value, and binary calls
# ============================================================
data <- merge(cbcndata, cedata, by = "all_cells") %>%
  as_tibble() %>%
  mutate(
    cb_mean = rowMeans(select(., corTim_cbp1:corTim_cbp6), na.rm = TRUE),
    cn_mean = rowMeans(select(., corTim_cnp1:corTim_cnp7), na.rm = TRUE)
  ) %>%
  select(
    cell = all_cells,
    cb_mean,
    cn_mean,
    cell_rmsd_mean,
    cell_rmsd_sd
  ) %>%
  merge(stage, by = "cell")

data <- data %>%
  mutate(
    cb_z_score = (cb_mean - cell_rmsd_mean) / cell_rmsd_sd,
    cb_p_value = pnorm(abs(cb_z_score), lower.tail = FALSE),
    cb_q_value = p.adjust(cb_p_value, method = "fdr"),
    cb_q_value_binary = case_when(
      cb_q_value < 0.01 & cb_z_score >= 0 ~ 1,
      cb_q_value < 0.01 & cb_z_score < 0 ~ -1,
      TRUE ~ 0
    ),
    cn_z_score = (cn_mean - cell_rmsd_mean) / cell_rmsd_sd,
    cn_p_value = pnorm(abs(cn_z_score), lower.tail = FALSE),
    cn_q_value = p.adjust(cn_p_value, method = "fdr"),
    cn_q_value_binary = case_when(
      cn_q_value < 0.01 & cn_z_score >= 0 ~ 1,
      cn_q_value < 0.01 & cn_z_score < 0 ~ -1,
      TRUE ~ 0
    )
  )

colnames(data)
data <- data %>% select(1, 10, 14)

names(data)[2] <- "cb_mean"
names(data)[3] <- "cn_mean"
colnames(data)

data <- data %>%
  mutate(
    cell_letter = substring(cell, 1, 1),
    cell_letter = factor(cell_letter, levels = c("A", "M", "E", "C", "D", "P"))
  )

# ============================================================
# Percentage of significant cells per cell-letter (CB vs CN)
# ============================================================
data_long <- data %>%
  pivot_longer(
    cols = c(cb_mean, cn_mean),
    names_to = "channel",
    values_to = "value"
  ) %>%
  filter(value == 1)

count_data <- data_long %>%
  group_by(channel, cell_letter) %>%
  summarise(Count = n(), .groups = "drop")

denom <- data %>%
  group_by(cell_letter) %>%
  summarise(total_cells = n(), .groups = "drop")

count_data <- count_data %>%
  left_join(denom, by = "cell_letter") %>%
  mutate(percentage = Count / total_cells * 100)

p2 <- ggplot(count_data, aes(x = cell_letter, y = percentage, fill = channel)) +
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 3
  ) +
  facet_wrap(~channel) +
  scale_fill_manual(values = c("cb_mean" = "#0073C2", "cn_mean" = "#EFC000")) +
  theme_minimal() +
  labs(
    title = "",
    x = "Cell Letter",
    y = "Percentage (%)"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.8)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_blank(),
    legend.position = "none"
  )

print(p2)
# ggsave("asychrony_signifi_cbcn_lineage_percent.pdf", plot = p2, width = 7, height = 4)

# ============================================================
# Counts by cell-letter for CB/CN-specific/common significance
# ============================================================
data_new <- data %>%
  filter(cb_mean == 1 | cn_mean == 1) %>%
  mutate(Condition = case_when(
    cb_mean == 1 & cn_mean == 0 ~ "Cbr specific",
    cb_mean == 0 & cn_mean == 1 ~ "Cni specific",
    cb_mean == 1 & cn_mean == 1 ~ "Common",
    TRUE ~ NA_character_
  ))

count_data <- data_new %>%
  group_by(Condition, cell_letter) %>%
  summarise(Count = n(), .groups = "drop")

p1 <- ggplot(count_data, aes(x = cell_letter, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Count), vjust = -0.3, size = 3) +
  facet_wrap(~Condition) +
  scale_fill_manual(values = c(
    "Cbr specific" = "#0073C2",
    "Cni specific" = "#EFC000",
    "Common" = "#c6dfc8"
  )) +
  scale_y_break(c(35, 130), expand = expansion(add = c(0, 5))) +
  theme_minimal() +
  labs(
    title = "",
    x = "CellFate",
    y = "Significant Number"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_blank(),
    legend.position = "none"
  )

print(p1)
# ggsave("asychrony_signifi_overlay_lineage_number.pdf", plot = p1, width = 7, height = 4)

# ============================================================
# Percentages by cell-letter for CB/CN-specific/common significance
# ============================================================
denom <- data %>%
  group_by(cell_letter) %>%
  summarise(Total = n(), .groups = "drop")

count_data <- count_data %>%
  left_join(denom, by = "cell_letter") %>%
  mutate(Percentage = Count / Total * 100)

p3 <- ggplot(count_data, aes(x = cell_letter, y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage)),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 3
  ) +
  facet_wrap(~Condition) +
  scale_fill_manual(values = c(
    "Cbr specific" = "#0073C2",
    "Cni specific" = "#EFC000",
    "Common" = "#c6dfc8"
  )) +
  theme_minimal() +
  labs(
    title = "",
    x = "CellFate",
    y = "Significant Percentage (%)"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.8)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_blank(),
    legend.position = "none"
  )

print(p3)
# ggsave("asychrony_signifi_overlay_lineage_percent.pdf", plot = p3, width = 7, height = 4)

# ============================================================
# Merge with tissue annotation and repeat CB/CN-specific/common analyses
# ============================================================
colnames(data)
colnames(tissue)
data <- merge(data, tissue, by = "cell")
colnames(data)

data_new <- data %>%
  filter(cb_mean == 1 | cn_mean == 1) %>%
  mutate(Condition = case_when(
    cb_mean == 1 & cn_mean == 0 ~ "Cbr specific",
    cb_mean == 0 & cn_mean == 1 ~ "Cni specific",
    cb_mean == 1 & cn_mean == 1 ~ "Common",
    TRUE ~ NA_character_
  ))

count_data <- data_new %>%
  group_by(Condition, CellFate) %>%
  summarise(Count = n(), .groups = "drop")

p1 <- ggplot(count_data, aes(x = CellFate, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Count), vjust = -0.3, size = 3) +
  facet_wrap(~Condition) +
  scale_fill_manual(values = c(
    "Cbr specific" = "#0073C2",
    "Cni specific" = "#EFC000",
    "Common" = "#c6dfc8"
  )) +
  scale_y_break(c(28, 130), expand = expansion(add = c(0, 5))) +
  theme_minimal() +
  labs(
    title = "",
    x = "CellFate",
    y = "Significant Number"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_blank(),
    legend.position = "none"
  )

print(p1)
# ggsave("asychrony_signifi_overlay_tissue.pdf", plot = p1, device = "pdf", width = 8, height = 5)






