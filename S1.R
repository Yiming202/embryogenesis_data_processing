# =============================================================================
# Nematode Embryonic Raw Cell Cycle Analysis
# -----------------------------------------------------------------------------
# Purpose:
#   1. Read tracking CSV files from multiple embryos and species
#   2. Compute, for each cell:
#        - Cell cycle length (number of frames × time resolution)
#        - Division time (last observed frame × time resolution)
#   3. Merge all embryos and annotate cells with developmental stages
#
# Inputs:
#   - CSV files 
#   - Predefined mapping of cell names to developmental stages (`stages`)
#
# Outputs:
#   - All_LengthTim_Absolute.csv
#       Raw length and division time for all cells and all embryos
#   - Length_Time_Absolute.csv
#       Same as above, but:
#         * filtered to remove a few founder cells
#         * each cell annotated with a developmental stage label
# =============================================================================

library(ggplot2)
library(dplyr)
library(tibble)
library(purrr)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Configuration: Define parameters for all datasets
# -----------------------------------------------------------------------------
# Each row in `datasets_config` describes one embryo:
#   - group      : species/group label (e.g.: ce, af, cb, cn)
#   - id         : embryo index within the group
#   - file_path  : path to the tracking CSV
#   - time_limit : last editing/tracking time point
#   - multiplier : time resolution for each embryo (mins per frame)
# =============================================================================

datasets_config <- tribble(
  ~group, ~id, ~file_path,                                    ~time_limit, ~multiplier,
  # CE group
  "ce", 1, "CDFile/ce/CD191108plc1p1.csv",                   205, 1.43,
  "ce", 2, "CDFile/ce/CD200109plc1p1.csv",                   205, 1.43,
  "ce", 3, "CDFile/ce/CD200113plc1p3.csv",                   195, 1.44,
  "ce", 4, "CDFile/ce/CD200113plc1p2.csv",                   205, 1.44,
  "ce", 5, "CDFile/ce/CD200322plc1p2.csv",                   195, 1.44,
  "ce", 6, "CDFile/ce/CD200323plc1p1.csv",                   185, 1.44,
  "ce", 7, "CDFile/ce/CD200326plc1p3.csv",                   220, 1.44,
  "ce", 8, "CDFile/ce/CD200326plc1p4.csv",                   195, 1.44,
  # AF group
  "af", 1, "CDFile/AF16/CD210519ZZY0874p1.csv",              170, 1.29,
  "af", 2, "CDFile/AF16/CD210519ZZY0874p5.csv",              175, 1.29,
  # CB group
  "cb", 1, "CDFile/she1/CD240731cbhis72p1.csv",              165, 1.57,
  "cb", 2, "CDFile/she1/CD240731cbhis72p2.csv",              175, 1.57,
  "cb", 3, "CDFile/she1/CD240731cbhis72p3.csv",              170, 1.57,
  "cb", 4, "CDFile/she1/CD241202cbhis72p1.csv",              160, 1.58,
  "cb", 5, "CDFile/she1/CD241202cbhis72p2.csv",              165, 1.58,
  "cb", 6, "CDFile/she1/CD241202cbhis72p4.csv",              180, 1.58,
  # CN group
  "cn", 1, "CDFile/cn/CD241202cnhis72p1.csv",                235, 1.58,
  "cn", 2, "CDFile/cn/CD240712cnhis72p1.csv",                235, 1.60,
  "cn", 3, "CDFile/cn/CD240712cnhis72p2.csv",                235, 1.60,
  "cn", 4, "CDFile/cn/CD240712cnhis72p3.csv",                235, 1.60,
  "cn", 5, "CDFile/cn/CD241207cnhis72p1.csv",                230, 1.65,
  "cn", 6, "CDFile/cn/CD241207cnhis72p3.csv",                230, 1.65,
  "cn", 7, "CDFile/cn/CD241202cnhis72p2.csv",                200, 1.58
)

# =============================================================================
# 2. Define developmental stages
# =============================================================================

stages <- list(
  Stage_AB4 = c("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C"),
  
  Stage_AB8 = c("ABala", "ABalp", "ABara", "ABarp", "ABpla", "ABplp", "ABpra", "ABprp",
                "MSa", "MSp", "Ca", "Cp", "P3"),
  
  Stage_AB16 = c("ABalaa", "ABalap", "ABalpa", "ABalpp", "ABaraa", "ABarap", "ABarpa", "ABarpp",
                 "ABplaa", "ABplap", "ABplpa", "ABplpp", "ABpraa", "ABprap", "ABprpa", "ABprpp",
                 "MSaa", "MSap", "MSpa", "MSpp", "Ea", "Ep", "D"),
  
  Stage_AB32 = c("ABalaaa", "ABalaap", "ABalapa", "ABalapp", "ABalpaa", "ABalpap", "ABalppa", "ABalppp",
                 "ABaraaa", "ABaraap", "ABarapa", "ABarapp", "ABarpaa", "ABarpap", "ABarppa", "ABarppp",
                 "ABplaaa", "ABplaap", "ABplapa", "ABplapp", "ABplpaa", "ABplpap", "ABplppa", "ABplppp",
                 "ABpraaa", "ABpraap", "ABprapa", "ABprapp", "ABprpaa", "ABprpap", "ABprppa", "ABprppp",
                 "MSaaa", "MSaap", "MSapa", "MSapp", "MSpaa", "MSpap", "MSppa", "MSppp",
                 "Eal", "Ear", "Epl", "Epr", "Caa", "Cap", "Cpa", "Cpp", "Da", "Dp", "P4"),
  
  Stage_AB64 = c("ABalaaaa", "ABalaaap", "ABalaapa", "ABalaapp", "ABalapaa", "ABalapap", "ABalappa", "ABalappp",
                 "ABalpaaa", "ABalpaap", "ABalpapa", "ABalpapp", "ABalppaa", "ABalppap", "ABalpppa", "ABalpppp",
                 "ABaraaaa", "ABaraaap", "ABaraapa", "ABaraapp", "ABarapaa", "ABarapap", "ABarappa", "ABarappp",
                 "ABarpaaa", "ABarpaap", "ABarpapa", "ABarpapp", "ABarppaa", "ABarppap", "ABarpppa", "ABarpppp",
                 "ABplaaaa", "ABplaaap", "ABplaapa", "ABplaapp", "ABplapaa", "ABplapap", "ABplappa", "ABplappp",
                 "ABplpaaa", "ABplpaap", "ABplpapa", "ABplpapp", "ABplppaa", "ABplppap", "ABplpppa", "ABplpppp",
                 "ABpraaaa", "ABpraaap", "ABpraapa", "ABpraapp", "ABprapaa", "ABprapap", "ABprappa", "ABprappp",
                 "ABprpaaa", "ABprpaap", "ABprpapa", "ABprpapp", "ABprppaa", "ABprppap", "ABprpppa", "ABprpppp",
                 "MSaaaa", "MSaaap", "MSaapa", "MSaapp", "MSapaa", "MSapap", "MSappa", "MSappp",
                 "MSpaaa", "MSpaap", "MSpapa", "MSpapp", "MSppaa", "MSppap", "MSpppa", "MSpppp",
                 "Caaa", "Caap", "Capa", "Capp", "Cpaa", "Cpap", "Cppa", "Cppp"),
  
  Stage_AB128 = c("ABalaaaal", "ABalaaaar", "ABalaaapa", "ABalaaapp", "ABalaapaa", "ABalaapap", "ABalaappa", "ABalaappp",
                  "ABalapaaa", "ABalapaap", "ABalapapa", "ABalapapp", "ABalappaa", "ABalappap", "ABalapppa", "ABalapppp",
                  "ABalpaaaa", "ABalpaaap", "ABalpaapa", "ABalpaapp", "ABalpapaa", "ABalpapap", "ABalpappa", "ABalpappp",
                  "ABalppaaa", "ABalppaap", "ABalppapa", "ABalppapp", "ABalpppaa", "ABalpppap", "ABalppppa", "ABalppppp",
                  "ABaraaaaa", "ABaraaaap", "ABaraaapa", "ABaraaapp", "ABaraapaa", "ABaraapap", "ABaraappa", "ABaraappp",
                  "ABarapaaa", "ABarapaap", "ABarapapa", "ABarapapp", "ABarappaa", "ABarappap", "ABarapppa", "ABarapppp",
                  "ABarpaaaa", "ABarpaaap", "ABarpaapa", "ABarpaapp", "ABarpapaa", "ABarpapap", "ABarpappa", "ABarpappp",
                  "ABarppaaa", "ABarppaap", "ABarppapa", "ABarppapp", "ABarpppaa", "ABarpppap", "ABarppppa", "ABarppppp",
                  "ABplaaaaa", "ABplaaaap", "ABplaaapa", "ABplaaapp", "ABplaapaa", "ABplaapap", "ABplaappa", "ABplaappp",
                  "ABplapaaa", "ABplapaap", "ABplapapa", "ABplapapp", "ABplappaa", "ABplappap", "ABplapppa", "ABplapppp",
                  "ABplpaaaa", "ABplpaaap", "ABplpaapa", "ABplpaapp", "ABplpapaa", "ABplpapap", "ABplpappa", "ABplpappp",
                  "ABplppaaa", "ABplppaap", "ABplppapa", "ABplppapp", "ABplpppaa", "ABplpppap", "ABplppppa", "ABplppppp",
                  "ABpraaaaa", "ABpraaaap", "ABpraaapa", "ABpraaapp", "ABpraapaa", "ABpraapap", "ABpraappa", "ABpraappp",
                  "ABprapaaa", "ABprapaap", "ABprapapa", "ABprapapp", "ABprappaa", "ABprappap", "ABprapppa", "ABprapppp",
                  "ABprpaaaa", "ABprpaaap", "ABprpaapa", "ABprpaapp", "ABprpapaa", "ABprpapap", "ABprpappa", "ABprpappp",
                  "ABprppaaa", "ABprppaap", "ABprppapa", "ABprppapp", "ABprpppaa", "ABprpppap", "ABprppppa", "ABprppppp",
                  "MSaaaaa", "MSaaaap", "MSaaapa", "MSaapaa", "MSaapap", "MSapaaa", "MSapapa", "MSapapp",
                  "MSpaaaa", "MSpaaap", "MSpaapa", "MSpapaa", "MSpapap", "MSppaaa", "MSppapa", "MSppapp",
                  "Eala", "Ealp", "Eara", "Earp", "Epla", "Eplp", "Epra", "Eprp",
                  "Caaaa", "Caaap", "Caapp", "Capaa", "Capap", "Cappa", "Cappp",
                  "Cpaaa", "Cpaap", "Cpapa", "Cpapp", "Cppaa", "Cppap", "Cpppa", "Cpppp",
                  "Daa", "Dap", "Dpa", "Dpp"),
  
  Stage_AB256 = c("ABalaaappr", "ABalaapppa", "ABalaapppp", "ABalapaapp", "ABalapappa", "ABalappapp", "ABalapppap",
                  "ABalappppa", "ABalappppp", "ABalpaapaa", "ABalpaapap", "ABalpaappa", "ABalpaappp", "ABalpapaaa",
                  "ABalpapaap", "ABalpapapa", "ABalpapapp", "ABalpappaa", "ABalpappap", "ABalpapppp", "ABalppapaa",
                  "ABalppapap", "ABalppappa", "ABalppappp", "ABalpppapa", "ABalpppapp", "ABalppppaa", "ABalppppap",
                  "ABalpppppa", "ABalpppppp", "ABaraaapaa", "ABaraaapap", "ABaraaappa", "ABaraaappp", "ABaraapaaa",
                  "ABaraapaap", "ABaraapapa", "ABaraapapp", "ABaraappaa", "ABaraappap", "ABaraapppa", "ABaraapppp",
                  "ABarapaapp", "ABarapapap", "ABarapappp", "ABarappaap", "ABarappapa", "ABarappapp", "ABarapppaa",
                  "ABarapppap", "ABarappppa", "ABarappppp", "ABarpapaap", "ABplaaaaap", "ABplaapaap", "ABplaapapa",
                  "ABplaapapp", "ABplapaaaa", "ABplapaaap", "ABplapappp", "ABplapppap", "ABplpaaaaa", "ABplpaaaap",
                  "ABplpaaapa", "ABplpaaapp", "ABplpaapaa", "ABplpaapap", "ABplpaappa", "ABplpaappp", "ABplpapaaa",
                  "ABplpapaap", "ABplpapapa", "ABplpapapp", "ABplpappaa", "ABplpapppa", "ABplppaaaa", "ABplppaapa",
                  "ABplppaapp", "ABplppapaa", "ABplppapap", "ABplppappa", "ABplppappp", "ABplpppaaa", "ABplpppaap",
                  "ABplpppapa", "ABplppppaa", "ABplppppap", "ABplpppppa", "ABplpppppp", "ABpraaappp", "ABpraapaap",
                  "ABpraapapp", "ABprapaaap", "ABprpaaaaa", "ABprpaaaap", "ABprpaaapa", "ABprpaaapp", "ABprpaapaa",
                  "ABprpaapap", "ABprpaappa", "ABprpaappp", "ABprpapaaa", "ABprpapaap", "ABprpapapa", "ABprpapapp",
                  "ABprpappaa", "ABprpappap", "ABprpapppa", "ABprpapppp", "ABprppaaaa", "ABprppaapp", "ABprppapaa",
                  "ABprppapap", "ABprppappa", "ABprppappp", "ABprpppaaa", "ABprpppaap", "ABprpppapa", "ABprppppaa",
                  "ABprppppap", "ABprpppppa", "ABprpppppp", "Caapa", "Capaaa", "Capaap", "Capapa", "Capapp",
                  "Cappaa", "Cappap", "Capppa", "Capppp", "Cppaaa", "Cppaap", "Cppapa", "Cppapp", "Cpppaa",
                  "Cpppap", "Cppppa", "Cppppp", "Daaa", "Daap", "Dapa", "Dapp", "Dpaa", "Dpap", "Dppa", "Dppp",
                  "MSaaapp", "MSaappa", "MSaappp", "MSapaap", "MSappaa", "MSappap", "MSapppa", "MSapppp",
                  "MSpappa", "MSpappp", "MSppaap", "MSpppaa", "MSpppap", "MSppppa", "MSppppp")
)

# =============================================================================
# 3. Processing functions
# =============================================================================

# Calculate raw/absolute cell cycle length (cell count × time resolution)
calc_length <- function(data, time_limit, multiplier) {
  data %>%
    filter(time <= time_limit) %>%
    group_by(cell) %>%
    tally() %>%
    mutate(value = n * multiplier) %>%
    select(cell, value)
}

# Calculate raw/absolute cell division time (max time × time resolution)
calc_time <- function(data, time_limit, multiplier) {
  data %>%
    filter(time <= time_limit) %>%
    group_by(cell) %>%
    summarise(value = max(time) * multiplier, .groups = "drop")
}

# Process a single dataset, calculating both length and time
process_dataset <- function(config_row) {
  data <- read.csv(config_row$file_path, header = TRUE)
  col_prefix <- paste0(config_row$group, "p", config_row$id)
  
  len_data <- calc_length(data, config_row$time_limit, config_row$multiplier) %>%
    rename(!!paste0("Len_", col_prefix) := value)
  
  tim_data <- calc_time(data, config_row$time_limit, config_row$multiplier) %>%
    rename(!!paste0("Tim_", col_prefix) := value)
  
  list(length = len_data, time = tim_data)
}

# Merge multiple data frames by cell
merge_all_by_cell <- function(df_list) {
  reduce(df_list, ~full_join(.x, .y, by = "cell"))
}

# Assign developmental stage based on cell name
assign_stage <- function(cell_name, stages_list) {
  for (stage_name in names(stages_list)) {
    if (cell_name %in% stages_list[[stage_name]]) {
      return(stage_name)
    }
  }
  return(NA_character_)
}

# =============================================================================
# 4. Main processing pipeline
# =============================================================================

# Process all datasets
cat("Processing data...\n")
results <- datasets_config %>%
  rowwise() %>%
  group_split() %>%
  map(process_dataset)

# Merge raw/absolute length and time data separately
mergeLen <- results %>% map("length") %>% merge_all_by_cell() %>% na.omit()
mergeTim <- results %>% map("time") %>% merge_all_by_cell() %>% na.omit()

# Merge all data
merge_all <- merge(mergeLen, mergeTim, by = "cell", all = TRUE) %>%
  filter(!cell %in% c("ABa", "ABp", "EMS", "P2"))

# Save intermediate results
write.csv(merge_all, file = "All_LengthTim_Absolute.csv", row.names = FALSE)
cat("Saved: All_LengthTim_Absolute.csv\n")

# Assign developmental stages
merge_all <- merge_all %>%
  mutate(stage = sapply(cell, assign_stage, stages_list = stages)) %>%
  filter(!is.na(stage))

# Save final results
write.csv(merge_all, file = "Length_Time_Absolute.csv", row.names = FALSE)
cat("Saved: Length_Time_Absolute.csv\n")
cat("Processing complete!\n")






# =============================================================================
# Nematode Normalization Cell Cycle Analysis
# -----------------------------------------------------------------------------
# Purpose:
#  Use elegans as a reference/template to normalize division times/lengths of other species.
#
# Steps:
#   1) Compute template mean length/time across ce1–ce8 for each cell
#   2) For each embryo of ce, cb, cn:
#        - Fit linear model: sample_time ~ template_time
#        - Use slope and intercept to normalize:
#            normalized_time = (time - intercept) / slope
#            normalized_length = length / slope
#   3) Merge normalized values across all embryos and attach stage labels
#
# Inputs:
#   - Length_Time_Absolute.csv
#
# Outputs:
#   - model_parameters.csv (slope & intercept per embryo)
#   - nor_Length_Time_correct.csv (normalized length & time per cell)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tibble)
library(purrr)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Configuration: Define sample groups
# =============================================================================

# Define all samples with their group and index
samples_config <- tribble(
  ~group, ~id,
  # CE group (template reference)
  "ce", 1, "ce", 2, "ce", 3, "ce", 4, "ce", 5, "ce", 6, "ce", 7, "ce", 8,
  # CB group
  "cb", 1, "cb", 2, "cb", 3, "cb", 4, "cb", 5, "cb", 6,
  # CN group
  "cn", 1, "cn", 2, "cn", 3, "cn", 4, "cn", 5, "cn", 6, "cn", 7
)

# =============================================================================
# 2. Load data and calculate template means
# =============================================================================

data <- read.csv("Length_Time_Absolute.csv", header = TRUE)

# Calculate template means from CE group (ce1-ce8)
template <- data %>%
  mutate(
    temLen_mean = round(rowMeans(select(., Len_cep1:Len_cep8), na.rm = TRUE), 2),
    temTim_mean = round(rowMeans(select(., Tim_cep1:Tim_cep8), na.rm = TRUE), 2)
  ) %>%
  select(cell, temLen_mean, temTim_mean)

# =============================================================================
# 3. Processing functions
# =============================================================================

# Fit linear model and create regression plot
fit_linear_model <- function(data, template, tim_col) {
  # Prepare data for regression
  reg_data <- data %>%
    select(cell, !!sym(tim_col)) %>%
    merge(template, by = "cell", all = TRUE) %>%
    na.omit() %>%
    select(-temLen_mean)
  
  # Fit linear model: sample_time ~ template_time
  model <- lm(reg_data[[tim_col]] ~ temTim_mean, data = reg_data)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  
  # Create regression plot
  plot <- ggplot(reg_data, aes(x = temTim_mean, y = .data[[tim_col]])) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "Linear Regression", x = "Template", y = tim_col) +
    xlim(0, 300) + ylim(0, 300) +
    annotate("text", x = 0, y = 300, hjust = 0, vjust = 1, size = 6,
             label = paste0("Intercept = ", round(intercept, 2), "\n",
                            "Slope = ", round(slope, 2))) +
    theme_bw()
  
  list(
    data = reg_data,
    intercept = intercept,
    slope = slope,
    plot = plot
  )
}

# Process a single sample: fit model and calculate normalized values
process_sample <- function(group, id, data, template) {
  tim_col <- paste0("Tim_", group, "p", id)
  len_col <- paste0("Len_", group, "p", id)
  
  # Fit linear regression model
  model_result <- fit_linear_model(data, template, tim_col)
  
  # Calculate normalized division time: (value - intercept) / slope
  nor_tim <- model_result$data %>%
    mutate(!!paste0("corTim_", group, "p", id) :=
             (.data[[tim_col]] - model_result$intercept) / model_result$slope) %>%
    select(cell, last_col())
  
  # Calculate normalized cell cycle length: value / slope
  nor_len <- data %>%
    select(cell, !!sym(len_col)) %>%
    mutate(!!paste0("corLen_", group, "p", id) := .data[[len_col]] / model_result$slope) %>%
    select(cell, last_col())
  
  list(
    group = group,
    id = id,
    model = model_result,
    nor_tim = nor_tim,
    nor_len = nor_len
  )
}

# =============================================================================
# 4. Main processing pipeline
# =============================================================================

cat("Processing samples...\n")

# Process all samples
results <- samples_config %>%
  rowwise() %>%
  group_split() %>%
  map(~process_sample(.x$group, .x$id, data, template))

# Display all regression plots (optional - uncomment to view)
# walk(results, ~print(.x$model$plot))

# =============================================================================
# 5. Extract and save model parameters (slope & intercept)
# =============================================================================

model_params <- results %>%
  map_dfr(~tibble(
    embryo = paste0(.x$group, "p", .x$id),
    group = .x$group,
    intercept = round(.x$model$intercept, 4),
    slope = round(.x$model$slope, 4)
  ))

# Display model parameters
cat("\n=== Model Parameters (Slope & Intercept) ===\n")
print(model_params, n = Inf)

# Save model parameters to CSV
write.csv(model_params, file = "model_parameters.csv", row.names = FALSE)
cat("\nSaved: model_parameters.csv\n")

# =============================================================================
# 6. Merge normalized values
# =============================================================================

# Merge all normalized division times
Nor_Tim_All <- results %>%
  map("nor_tim") %>%
  reduce(~full_join(.x, .y, by = "cell"))

# Merge all normalized cell cycle lengths
Nor_Len_All <- results %>%
  map("nor_len") %>%
  reduce(~full_join(.x, .y, by = "cell"))

# Combine length and time data
nor_Length_time_correct <- Nor_Len_All %>%
  full_join(Nor_Tim_All, by = "cell")

# Add stage column from original data
nor_Length_time_correct <- nor_Length_time_correct %>%
  left_join(data %>% select(cell, stage), by = "cell")

# =============================================================================
# 7. Save results
# =============================================================================

write.csv(nor_Length_time_correct, file = "nor_Length_Time_correct.csv", row.names = FALSE)
cat("Saved: nor_Length_Time_correct.csv\n")
cat("Processing complete!\n")




# =============================================================================
#  PLOT_PCA Module
# -----------------------------------------------------------------------------
# Purpose:
#   Perform PCA on absolute cell cycle measures (length/time) across embryos
#   and visualize embryo clustering in PC1–PC2 space with k-means.
#
# Input:
#   - Length_Time_Absolute.csv (first 24 columns are Len_*/Tim_* per embryo)
#
# Output:
#   - PCA scatter plot with k-means clusters and group-based colors
# =============================================================================

library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(ggalt)
library(stringr)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Load and prepare data
# =============================================================================

data <- read.csv("Length_Time_Absolute.csv", header = TRUE) %>%
  select(1:24)

# Set cell names as row names, then remove cell column
rownames(data) <- data$cell
data$cell <- NULL

# Transpose data: samples (columns) become rows, features become columns
sample_data <- t(data)

# =============================================================================
# 2. Perform PCA analysis
# =============================================================================

pca_result <- prcomp(sample_data, center = TRUE, scale. = TRUE)
summary(pca_result)

# Convert PCA scores to data frame with sample names
pca_scores <- pca_result$x %>%
  as.data.frame() %>%
  mutate(sample = rownames(.))

# Extract variance explained by PC1 and PC2
pca_summary <- summary(pca_result)
pc1_var <- round(pca_summary$importance[2, 1] * 100, 1)
pc2_var <- round(pca_summary$importance[2, 2] * 100, 1)

# =============================================================================
# 3. K-means clustering
# =============================================================================

set.seed(123)
k_result <- kmeans(pca_scores[, c("PC1", "PC2")], centers = 3)
pca_scores$close_group <- as.factor(k_result$cluster)

# =============================================================================
# 4. Define color scheme
# =============================================================================

# Color mapping function for sample groups
get_sample_color <- function(sample) {
  case_when(
    str_detect(sample, "^Len_cep") ~ "#CD534C",   # CE group: red
    str_detect(sample, "^Len_cbp") ~ "#0073C2",   # CB group: blue
    str_detect(sample, "^Len_afp") ~ "#1b7c3d",   # AF group: green
    str_detect(sample, "^Len_cnp") ~ "#EFC000",   # CN group: yellow
    TRUE ~ "black"                                 # Default: black
  )
}

# =============================================================================
# 5. Create PCA plot
# =============================================================================

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, label = sample)) +
  # Encircle k-means clusters
  geom_encircle(
    aes(group = close_group),
    color = "black",
    expand = 0.09,
    size = 2,
    s_shape = 1,
    alpha = 0.5
  ) +
  # Plot points with group-specific colors
  geom_point(aes(color = get_sample_color(sample)), size = 3) +
  # Axis labels with variance explained
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  ggtitle("PCA of Samples") +
  # Use color values directly
  scale_color_identity() +
  # Expand axis ranges for encircle visibility
  scale_x_continuous(expand = expansion(mult = c(0.13, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  # Theme settings
  
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

# Display plot
print(p_pca)

# Save plot
# ggsave("PCA_plot.pdf", plot = p_pca, device = "pdf", width = 5, height = 5)







# =============================================================================
#  PLOT_CUMMULATIVE FREQUENCY
# =============================================================================
# =============================================================================
#  Cell Cycle Length - Cumulative Frequency Analysis
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Load data and calculate group means
# =============================================================================

data <- read.csv("Length_Time_Absolute.csv") %>%
  mutate(
    ceLen_mean = round(rowMeans(select(., Len_cep1:Len_cep8), na.rm = TRUE), 2),
    cbLen_mean = round(rowMeans(select(., Len_cbp1:Len_cbp6), na.rm = TRUE), 2),
    cnLen_mean = round(rowMeans(select(., Len_cnp1:Len_cnp7), na.rm = TRUE), 2)
  ) %>%
  select(cell, ceLen_mean, cbLen_mean, cnLen_mean)

# =============================================================================
# 2. Transform to long format and calculate cumulative frequency
# =============================================================================

# Type name mapping
type_labels <- c(ceLen_mean = "ce", cbLen_mean = "cb", cnLen_mean = "cn")

data_long <- data %>%
  pivot_longer(
    cols = c(ceLen_mean, cbLen_mean, cnLen_mean),
    names_to = "type",
    values_to = "length"
  ) %>%
  mutate(type = recode(type, !!!type_labels)) %>%
  group_by(type) %>%
  arrange(length, .by_group = TRUE) %>%
  mutate(rel_cum_freq = row_number() / n() * 100) %>%
  ungroup()

# =============================================================================
# 3. Plot cumulative relative frequency
# =============================================================================

# Color scheme for groups
group_colors <- c(ce = "#CD534C", cb = "#0073C2", cn = "#EFC000")

p_cum <- ggplot(data_long, aes(x = length, y = rel_cum_freq, color = type)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = seq(0, 140, by = 30)) +
  labs(
    x = "Length",
    y = "Cumulative Relative Frequency (%)",
    title = "Raw Cell Length",
    color = "Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

print(p_cum)

# Save plot
# ggsave("Cumulative_RelFreq_length.pdf", plot = p_cum, width = 9, height = 8)




# =============================================================================
#  PLOT_ABSOLUTE LENGTH/DIVISDION TIME PAIRWISE COMPARISON
# =============================================================================
# =============================================================================
#  Pairwise Group Comparisons (Length & Division Time)
# =============================================================================

library(ggplot2)
library(dplyr)
library(gridExtra)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Common theme and plotting function
# =============================================================================

scatter_theme <- theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

create_scatter_plot <- function(data, x_var, y_var, title, axis_limits,
                                label_data = NULL, label_hjust = -0.1, label_vjust = -1) {
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = "#8da4c3", size = 1.5, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    scale_x_continuous(limits = c(0, axis_limits)) +
    scale_y_continuous(limits = c(0, axis_limits)) +
    labs(title = title, x = x_var, y = y_var) +
    scatter_theme
  
  if (!is.null(label_data) && nrow(label_data) > 0) {
    p <- p + geom_text(data = label_data, aes(label = cell),
                       hjust = label_hjust, vjust = label_vjust,
                       color = "black", size = 4)
  }
  return(p)
}

# =============================================================================
# 2. Cell Cycle LENGTH Comparison
# =============================================================================

data_len <- read.csv("Length_Time_Absolute.csv") %>%
  select(cell, starts_with("Len_cep"), starts_with("Len_cbp"), starts_with("Len_cnp"), stage) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  rowwise() %>%
  mutate(
    ce = mean(c_across(starts_with("Len_cep")), na.rm = TRUE),
    cb = mean(c_across(starts_with("Len_cbp")), na.rm = TRUE),
    cn = mean(c_across(starts_with("Len_cnp")), na.rm = TRUE)
  ) %>%
  ungroup()

# Length plots
furthest_len <- data_len %>%
  filter(cb > ce) %>%
  mutate(distance = (cb - ce) / sqrt(2)) %>%
  slice_max(distance, n = 1)

p_len_ce_cb <- create_scatter_plot(data_len, "ce", "cb", "Length: ce vs cb", 135,
                                   label_data = furthest_len, label_vjust = -1)

p_len_ce_cn <- create_scatter_plot(data_len, "ce", "cn", "ce vs cn", 135)

p_len_cb_cn <- create_scatter_plot(data_len, "cb", "cn", "cb vs cn", 135,
                                   label_data = data_len %>% filter(cn < cb),
                                   label_vjust = 2)

print(p_len_ce_cb)
print(p_len_ce_cn)
print(p_len_cb_cn)
grid.arrange(p_len_ce_cb, p_len_ce_cn, p_len_cb_cn, ncol = 3)

# ggsave(filename = "draft_code/fig2_all/f2_sup2.pdf",
#        plot = grid.arrange(p_len_ce_cb, p_len_ce_cn, p_len_cb_cn, ncol = 3),
#        width = 10, height = 3.2)

# =============================================================================
# 3. Division TIME Comparison
# =============================================================================

data_tim <- read.csv("Length_Time_Absolute.csv") %>%
  select(cell, starts_with("Tim_cep"), starts_with("Tim_cbp"), starts_with("Tim_cnp"), stage) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  rowwise() %>%
  mutate(
    ce = mean(c_across(starts_with("Tim_cep")), na.rm = TRUE),
    cb = mean(c_across(starts_with("Tim_cbp")), na.rm = TRUE),
    cn = mean(c_across(starts_with("Tim_cnp")), na.rm = TRUE)
  ) %>%
  ungroup()

# Time plots
furthest_tim <- data_tim %>%
  filter(cb > ce) %>%
  mutate(distance = (cb - ce) / sqrt(2)) %>%
  slice_max(distance, n = 1)

p_tim_ce_cb <- create_scatter_plot(data_tim, "ce", "cb", "Division Time: ce vs cb", 280,
                                   label_data = furthest_tim, label_vjust = -1)

p_tim_ce_cn <- create_scatter_plot(data_tim, "ce", "cn", "ce vs cn", 360)

p_tim_cb_cn <- create_scatter_plot(data_tim, "cb", "cn", "cb vs cn", 360,
                                   label_data = data_tim %>% filter(cn < cb),
                                   label_vjust = 2)

print(p_tim_ce_cb)
print(p_tim_ce_cn)
print(p_tim_cb_cn)
grid.arrange(p_tim_ce_cb, p_tim_ce_cn, p_tim_cb_cn, ncol = 3)

# Save combined plot (uncomment to use)
# ggsave(filename = "Dotplots_time_ce_cb_cn.pdf",
#        plot = grid.arrange(p_tim_ce_cb, p_tim_ce_cn, p_tim_cb_cn, ncol = 3),
#        height = 6, width = 22)



# =============================================================================
#  PLOT_DEVELOPMENTAL RATE
# =============================================================================
# =============================================================================
#  Developmental Rate Analysis (Length vs Division Time)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpmisc)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Load and prepare data
# =============================================================================
library(dplyr)

# data <- read.csv("Length_Time_Absolute.csv", header = TRUE) %>%
#   mutate(
#     xorder = substr(cell, 1, 1),
#     xorder = ifelse(xorder %in% c("P", "Z"), "P", xorder),
#     xorder = recode(xorder, "A" = "AB", "M" = "MS"),
#     xorder = factor(xorder, levels = c("AB", "MS", "E", "C", "D", "P"))
#   ) %>%
#   na.omit()

data <- read.csv("Length_Time_Absolute.csv", header = TRUE) %>%
  mutate(
    xorder = substr(cell, 1, 1),
    xorder = factor(ifelse(xorder %in% c("P", "Z"), "P", xorder),
                    levels = c("A", "M", "E", "C", "D", "P"))
  ) %>%
  na.omit()



# data <- data %>%
#   mutate(
#     celLen_mean = round(rowMeans(select(., corLen_cep1:corLen_cep8), na.rm = TRUE), 2),
#     cbrLen_mean = round(rowMeans(select(., corLen_cbp1:corLen_cbp6), na.rm = TRUE), 2),
#     cniLen_mean = round(rowMeans(select(., corLen_cnp1:corLen_cnp7), na.rm = TRUE), 2)
#   ) %>%
#   rowwise() %>%
#   mutate(
#     celLen_cv = sd(c_across(corLen_cep1:corLen_cep8)) / celLen_mean,
#     cbrLen_cv = sd(c_across(corLen_cbp1:corLen_cbp6)) / cbrLen_mean,
#     cniLen_cv = sd(c_across(corLen_cnp1:corLen_cnp7)) / cniLen_mean
#   ) %>%
#   ungroup() 
# 
# write.csv(data, file = "Absolute_Len_Time_260302.csv", row.names = FALSE)

# data <- read.csv("nor_Length_Time_correct.csv", header = TRUE) %>%
#   mutate(
#     xorder = substr(cell, 1, 1),
#     xorder = ifelse(xorder %in% c("P", "Z"), "P", xorder),
#     xorder = recode(xorder, "A" = "AB", "M" = "MS"),
#     xorder = factor(xorder, levels = c("AB", "MS", "E", "C", "D", "P"))
#   ) %>%
#   na.omit()
# 
# write.csv(data, file = "Nor_Len_Time_260302.csv", row.names = FALSE)


# =============================================================================
# 2. Transform to long format for plotting
# =============================================================================

plotdata <- data %>%
  pivot_longer(
    cols = matches("^(Len|Tim)_"),
    names_to = c("Measure", "Group", "Rep"),
    names_pattern = "^(Len|Tim)_([a-z]+)(\\d+)$",
    values_to = "Value"
  ) %>%
  pivot_wider(names_from = Measure, values_from = Value) %>%
  filter(Group %in% c("cep", "cbp", "cnp")) %>%
  mutate(
    Group = recode(Group, cep = "ce", cbp = "cb", cnp = "cn"),
    Group = factor(Group, levels = c("ce", "cb", "cn"))
  )

# =============================================================================
# 3. Define color scheme
# =============================================================================

group_colors <- c(ce = "#CD534C", cb = "#0073C2", cn = "#EFC000")

# =============================================================================
# 4. Create faceted scatter plot with regression lines
# =============================================================================

p_grid <- ggplot(plotdata, aes(x = Tim, y = Len, color = Group, fill = Group)) +
  # Scatter points
  geom_point(size = 1, alpha = 0.6, shape = 16) +
  # Black regression line (background)
  geom_smooth(
    method = "lm",
    linewidth = 1.5,
    se = FALSE,
    color = adjustcolor("black", alpha.f = 0.5),
    show.legend = FALSE
  ) +
  # R² annotation
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(..rr.label..)),
    parse = TRUE,
    x.npc = 0.05,
    y.npc = 0.95,
    size = 6
  ) +
  # Color scales
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  # Facet by group
  facet_wrap(~Group, nrow = 1) +
  # Fixed aspect ratio
  coord_fixed(ratio = 2.6) +
  # Axis settings
  scale_x_continuous(
    name = "time of division (min)",
    limits = c(0, 400),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_continuous(
    name = "cell length (min)",
    limits = c(0, 160),
    breaks = seq(0, 160, by = 40),
    expand = expansion(mult = c(0, 0.05))
  ) +
  # Theme
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p_grid)

# Save plot
# ggsave("DF1_DevelopmentalRate_Lineage_abs_all.pdf", plot = p_grid, height = 6, width = 8)





# =============================================================================
#  PLOT_CV BOX PLOT
# =============================================================================
# =============================================================================
#  Cell Cycle Length - Coefficient of Variation (CV) Analysis
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Load data and calculate CV (std/mean) for each group
# =============================================================================

data <- read.csv("Length_Time_Absolute.csv") %>%
  mutate(
    ceLen_mean = round(rowMeans(select(., Len_cep1:Len_cep8), na.rm = TRUE), 2),
    cbLen_mean = round(rowMeans(select(., Len_cbp1:Len_cbp6), na.rm = TRUE), 2),
    cnLen_mean = round(rowMeans(select(., Len_cnp1:Len_cnp7), na.rm = TRUE), 2)
  ) %>%
  rowwise() %>%
  mutate(
    ce = sd(c_across(Len_cep1:Len_cep8)) / ceLen_mean,
    cb = sd(c_across(Len_cbp1:Len_cbp6)) / cbLen_mean,
    cn = sd(c_across(Len_cnp1:Len_cnp7)) / cnLen_mean
  ) %>%
  ungroup() %>%
  select(cell, ce, cb, cn, stage)

# =============================================================================
# 2. Transform to long format
# =============================================================================

data_long <- data %>%
  pivot_longer(
    cols = c(ce, cb, cn),
    names_to = "measurement",
    values_to = "value"
  ) %>%
  mutate(measurement = factor(measurement, levels = c("ce", "cb", "cn")))

# =============================================================================
# 3. Statistical comparison (Wilcoxon test)
# =============================================================================

my_comparisons <- list(c("ce", "cb"), c("ce", "cn"), c("cb", "cn"))

comparison_results <- compare_means(
  formula = value ~ measurement,
  data = data_long,
  method = "wilcox.test",
  comparisons = my_comparisons
)

# Replace p-values < 1e-04 with 1e-04 to avoid symnum errors
comparison_results$p <- pmax(comparison_results$p, 1e-04)

# Calculate significance symbols
comparison_results$p.signif <- symnum(
  comparison_results$p,
  corr = FALSE,
  cutpoints = c(0.0001, 0.001, 0.01, 0.05, 1),
  symbols = c("***", "**", "*", "ns")
)

print(comparison_results)

# =============================================================================
# 4. Filter outliers for plotting
# =============================================================================

data_long_filtered <- data_long %>%
  group_by(measurement) %>%
  mutate(
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1
  ) %>%
  filter(value >= (Q1 - 1.5 * IQR), value <= (Q3 + 1.5 * IQR)) %>%
  ungroup() %>%
  select(-Q1, -Q3, -IQR)

# Set y positions for significance annotations
y_max <- max(data_long_filtered$value, na.rm = TRUE)
comparison_results$y.position <- c(y_max * 1.05, y_max * 1.10, y_max * 1.15)

# =============================================================================
# 5. Create boxplot
# =============================================================================

group_colors <- c(ce = "#CD534C", cb = "#0073C2", cn = "#EFC000")

p_box <- ggplot(data_long_filtered, aes(x = measurement, y = value, fill = measurement)) +
  stat_boxplot(geom = "errorbar", width = 0.4, size = 0.5) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    color = "black",
    size = 0.5
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = "Measurement",
    y = "Value",
    title = "Raw Cell Length (CV)"
  ) +
  coord_cartesian(ylim = c(0, 0.17)) +
  scale_y_continuous(expand = c(0.005, 0)) +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

print(p_box)

# Save plot (uncomment to use)
# ggsave("cv_boxplot_all.pdf", plot = p_box, height = 8, width = 5)



# =============================================================================
#  PLOT_LINE PLOT AND COUNT
# =============================================================================
# =============================================================================
# C. elegans Cell Cycle Length - Stage-wise Comparison with Treemap
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(treemapify)
library(patchwork)
library(purrr)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Load and prepare data
# =============================================================================

data <- read.csv("All_LengthTim_Absolute.csv") %>%
  mutate(
    ceLen_mean = round(rowMeans(select(., Len_cep1:Len_cep8), na.rm = TRUE), 2),
    cbLen_mean = round(rowMeans(select(., Len_cbp1:Len_cbp6), na.rm = TRUE), 2),
    cnLen_mean = round(rowMeans(select(., Len_cnp1:Len_cnp7), na.rm = TRUE), 2)
  ) %>%
  select(cell, ceLen_mean, cbLen_mean, cnLen_mean)

# Add stage information
stage <- read.csv("nor_Length_Time_correct.csv") %>% select(cell, stage)
data <- merge(data, stage, by = "cell")

# =============================================================================
# 2. Define ordering and color scheme
# =============================================================================

order_stage <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                 "Stage_AB64", "Stage_AB128", "Stage_AB256")
order_letter <- c("A", "M", "E", "C", "D", "P")
metric_colors <- c(ceLen_mean = "#CD534C", cbLen_mean = "#0073C2", cnLen_mean = "#EFC000")

# =============================================================================
# 3. Transform to long format and calculate x coordinates
# =============================================================================

data$stage <- factor(data$stage, levels = order_stage)

data_long <- data %>%
  pivot_longer(
    cols = c(ceLen_mean, cbLen_mean, cnLen_mean),
    names_to = "Metric",
    values_to = "Len_mean"
  ) %>%
  mutate(
    stage = factor(stage, levels = order_stage),
    first_letter = factor(substr(cell, 1, 1), levels = order_letter)
  ) %>%
  arrange(stage, first_letter, cell)

data_long$cell <- factor(data_long$cell, levels = unique(data_long$cell))

# Calculate new x coordinates with gaps between stages
gap_size <- 2

cell_order <- data_long %>%
  distinct(cell, stage) %>%
  group_by(stage) %>%
  mutate(cell_index = row_number()) %>%
  ungroup()

stage_counts <- cell_order %>%
  group_by(stage) %>%
  summarize(count = n(), .groups = "drop") %>%
  arrange(factor(stage, levels = order_stage)) %>%
  mutate(cum = lag(cumsum(count), default = 0))

cell_order <- cell_order %>%
  left_join(stage_counts, by = "stage") %>%
  mutate(new_x = cum + cell_index + gap_size * (as.integer(stage) - 1))

data_long <- left_join(data_long, cell_order %>% select(cell, new_x), by = "cell")

# =============================================================================
# 4. Calculate min/mid/max values for fill regions
# =============================================================================

fill_df <- data_long %>%
  group_by(cell) %>%
  arrange(Len_mean) %>%
  summarize(
    min_value = first(Len_mean), Metric_min = first(Metric),
    mid_value = nth(Len_mean, 2), Metric_mid = nth(Metric, 2),
    max_value = last(Len_mean), Metric_max = last(Metric),
    x = first(new_x),
    .groups = "drop"
  ) %>%
  mutate(
    fill_min = metric_colors[Metric_min],
    fill_mid = metric_colors[Metric_mid],
    fill_max = metric_colors[Metric_max]
  )

# =============================================================================
# 5. Calculate stage borders
# =============================================================================

stage_border <- cell_order %>%
  group_by(stage) %>%
  summarize(xmin = min(new_x) - 0.5, xmax = max(new_x) + 0.5, .groups = "drop") %>%
  mutate(width = xmax - xmin)

# =============================================================================
# 6. Create main plot
# =============================================================================

p_main <- ggplot(data_long, aes(x = new_x, y = Len_mean,
                                group = interaction(Metric, stage), color = Metric)) +
  # Stage borders
  geom_rect(data = stage_border,
            aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
            fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
  # Fill regions (bottom to top: min, mid, max)
  geom_rect(data = fill_df, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 0, ymax = min_value, fill = fill_min),
            inherit.aes = FALSE, alpha = 0.4) +
  geom_rect(data = fill_df, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = min_value, ymax = mid_value, fill = fill_mid),
            inherit.aes = FALSE, alpha = 0.4) +
  geom_rect(data = fill_df, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = mid_value, ymax = max_value, fill = fill_max),
            inherit.aes = FALSE, alpha = 0.4) +
  # Line plot
  geom_line(size = 0.4) +
  scale_color_manual(values = metric_colors) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0, 0), breaks = cell_order$new_x, labels = cell_order$cell) +
  labs(x = "Cell", y = "Len_mean") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(p_main)

# =============================================================================
# 7. Create treemap for each stage
# =============================================================================

tree_data <- fill_df %>%
  left_join(cell_order %>% select(cell, stage), by = "cell") %>%
  group_by(stage, Metric_max) %>%
  summarize(count = n(), .groups = "drop")

stage_widths <- setNames(stage_border$width, stage_border$stage)

treemap_list <- map(order_stage, function(stg) {
  df_sub <- tree_data %>% filter(stage == stg)
  ggplot(df_sub, aes(area = count, fill = Metric_max, label = paste(count))) +
    geom_treemap(alpha = 0.6) +
    geom_treemap_text(color = "white", place = "centre", grow = FALSE, size = 10, min.size = 0) +
    scale_fill_manual(values = metric_colors) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(1, 1, 1, 2.5)
    )
})

# =============================================================================
# 8. Combine treemap and main plot
# =============================================================================

p_tree_combined <- wrap_plots(treemap_list, nrow = 1, widths = stage_widths)
p_sup <- p_tree_combined / p_main + plot_layout(heights = c(0.5, 3))

print(p_sup)

# Save plot (uncomment to use)
# ggsave("draft_code/fig2_all/f2_sup1.pdf", p_sup, width = 10, height = 3)





# =============================================================================
#  PLOT_DEVELPOPMENTAL RATE FACETED_LINEAR REGRESSION
# =============================================================================
# =============================================================================
# Nematode Developmental Rate Analysis - Faceted by Lineage and Stage
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpmisc)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Load and prepare data
# =============================================================================

data <- read.csv("Length_Time_Absolute.csv", header = TRUE) %>%
  mutate(
    Len_ce = rowMeans(select(., Len_cep1:Len_cep8), na.rm = TRUE),
    Len_cb = rowMeans(select(., Len_cbp1:Len_cbp6), na.rm = TRUE),
    Len_cn = rowMeans(select(., Len_cnp1:Len_cnp7), na.rm = TRUE),
    Tim_ce = rowMeans(select(., Tim_cep1:Tim_cep8), na.rm = TRUE),
    Tim_cb = rowMeans(select(., Tim_cbp1:Tim_cbp6), na.rm = TRUE),
    Tim_cn = rowMeans(select(., Tim_cnp1:Tim_cnp7), na.rm = TRUE)
  ) %>%
  select(cell, Len_ce, Len_cb, Len_cn, Tim_ce, Tim_cb, Tim_cn, stage) %>%
  mutate(
    xorder = factor(ifelse(substr(cell, 1, 1) %in% c("P", "Z"), "P", substr(cell, 1, 1)),
                    levels = c("A", "M", "E", "C", "D", "P")),
    stage = factor(stage, levels = c("Stage_AB4", "Stage_AB8", "Stage_AB16",
                                     "Stage_AB32", "Stage_AB64", "Stage_AB128", "Stage_AB256"))
  ) %>%
  na.omit()

# =============================================================================
# 2. Transform to long format
# =============================================================================

plotdata <- data %>%
  pivot_longer(
    cols = matches("^(Len|Tim)_"),
    names_to = c("Measure", "Group"),
    names_pattern = "^(Len|Tim)_([a-z]+)$",
    values_to = "Value"
  ) %>%
  pivot_wider(names_from = Measure, values_from = Value) %>%
  filter(Group %in% c("ce", "cb", "cn")) %>%
  mutate(Group = factor(Group, levels = c("ce", "cb", "cn")))

# =============================================================================
# 3. Define color scheme and common plot elements
# =============================================================================

group_colors <- c(ce = "#CD534C", cb = "#0073C2", cn = "#EFC000")

base_theme <- theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

# Function to create faceted developmental rate plot
create_faceted_plot <- function(data, facet_formula) {
  ggplot(data, aes(x = Tim, y = Len, color = Group, fill = Group)) +
    geom_point(size = ifelse(deparse(facet_formula) == "~xorder", 2, 1), alpha = 0.6, shape = 16) +
    geom_smooth(method = "lm", linewidth = 1, se = FALSE,
                color = adjustcolor("black", alpha.f = 0.5), show.legend = FALSE) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..rr.label..)),
                 parse = TRUE, x.npc = 0.05, y.npc = 0.95, size = 6) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    facet_grid(facet_formula) +
    coord_fixed(ratio = 2.6) +
    scale_x_continuous(name = "time of division (min)", limits = c(0, 400),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(name = "cell length (min)", limits = c(0, 160),
                       breaks = seq(0, 160, by = 40), expand = expansion(mult = c(0, 0.05))) +
    base_theme
}

# =============================================================================
# 4. Create plots
# =============================================================================

# Plot 1: Faceted by Group (rows) x Lineage (columns)
p_lineage <- ggplot(plotdata, aes(x = Tim, y = Len, color = Group, fill = Group)) +
  geom_point(size = 2, alpha = 0.6, shape = 16) +
  geom_smooth(method = "lm", linewidth = 1, se = FALSE,
              color = adjustcolor("black", alpha.f = 0.5), show.legend = FALSE) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(..rr.label..)),
               parse = TRUE, x.npc = 0.05, y.npc = 0.95, size = 6) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_grid(Group ~ xorder) +
  coord_fixed(ratio = 2.6) +
  scale_x_continuous(name = "time of division (min)", limits = c(0, 400),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(name = "cell length (min)", limits = c(0, 160),
                     breaks = seq(0, 160, by = 40), expand = expansion(mult = c(0, 0.05))) +
  base_theme

print(p_lineage)
# ggsave("f2_sup3_lineage.pdf", plot = p_lineage, height = 4, width = 10)

# Plot 2: Faceted by Group (rows) x Stage (columns)
p_stage <- ggplot(plotdata, aes(x = Tim, y = Len, color = Group, fill = Group)) +
  geom_point(size = 1, alpha = 0.6, shape = 16) +
  geom_smooth(method = "lm", linewidth = 1, se = FALSE,
              color = adjustcolor("black", alpha.f = 0.5), show.legend = FALSE) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(..rr.label..)),
               parse = TRUE, x.npc = 0.05, y.npc = 0.95, size = 6) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_grid(Group ~ stage) +
  coord_fixed(ratio = 2.6) +
  scale_x_continuous(name = "time of division (min)", limits = c(0, 400),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(name = "cell length (min)", limits = c(0, 160),
                     breaks = seq(0, 160, by = 40), expand = expansion(mult = c(0, 0.05))) +
  base_theme

print(p_stage)
# ggsave("Raw_length_Tim_stage.pdf", plot = p_stage, height = 4, width = 10)




# =============================================================================
#  PLOT_TREE
# =============================================================================
# =============================================================================
# Nematode Cell Lineage Tree Visualization with CV Coloring
# =============================================================================

library(dplyr)    
library(ggplot2)
library(reshape2)
library(stringr)
library(magrittr)
library(patchwork)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Helper functions for tree layout
# =============================================================================

# Split cell name into parent and last character
trim_cell <- function(cell) {
  c(str_sub(cell, 1, -2), str_sub(cell, -1))
}

# Calculate x-offset for tree positioning
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
# 2. Configuration
# =============================================================================

# Tree layout adjustment parameters
KLIST <- list(Za = 2, Zpap = 0.3, Zppa = 0.5, Zpppa = 0.8, Zppp = 0.8, Zpppp = 0.3)
SLIST <- list(Za = -0.5, Zp = -0.2, Zpa = 0.22, Zpap = -0.08, Zpp = -0.04,
              Zppa = 0.04, Zppp = -0.06, Zpppa = 0.01, Zpppp = -0.03)

# Display settings
ID_SHOW_LEN <- 5
BG_COLOR <- "darkgray"

# CV color breaks and colors
CV_BREAKS <- c(-Inf, 0.07, 0.10, Inf)
CV_COLORS <- c(red = "#fff5f5", green = "#d16d5b", blue = "red")

# =============================================================================
# 3. Load reference data
# =============================================================================

# Cell ID mapping
AllLineage <- read.csv("CellID.csv") %>% tibble()

# Calculate CV for each species
df_cv <- read.csv("Length_Time_Absolute.csv") %>%
  mutate(
    ceLen_mean = round(rowMeans(select(., Len_cep1:Len_cep8), na.rm = TRUE), 2),
    cbLen_mean = round(rowMeans(select(., Len_cbp1:Len_cbp6), na.rm = TRUE), 2),
    cnLen_mean = round(rowMeans(select(., Len_cnp1:Len_cnp7), na.rm = TRUE), 2)
  ) %>%
  rowwise() %>%
  mutate(
    cv_len_ce = sd(c_across(Len_cep1:Len_cep8)) / ceLen_mean,
    cv_len_cb = sd(c_across(Len_cbp1:Len_cbp6)) / cbLen_mean,
    cv_len_cn = sd(c_across(Len_cnp1:Len_cnp7)) / cnLen_mean
  ) %>%
  ungroup() %>%
  select(CellName = cell, cv_len_ce, cv_len_cb, cv_len_cn)

# =============================================================================
# 4. Main plotting function
# =============================================================================

process_lineage_and_plot <- function(lineage_file, cv_column, time_limit, time_factor, data_type) {
  
  # Load and process lineage data
  lineage <- read.csv(lineage_file) %>%
    as_tibble() %>%
    filter(time <= time_limit) %>%
    mutate(time = time * time_factor) %>%
    group_by(cell) %>%
    summarize(Start = min(time), End = max(time), .groups = "drop") %>%
    rename(CellName = cell) %>%
    left_join(AllLineage %>% select(CellName, ID, CellFate), by = "CellName") %>%
    rename(Fate = CellFate) %>%
    left_join(df_cv %>% select(CellName, !!sym(cv_column)), by = "CellName")
  
  # Add root nodes based on species
  root_start <- if (data_type == "ce") -15 else -20
  lineage %<>%
    add_row(CellName = "P0", Start = root_start, End = if (data_type == "ce") 0 else -5,
            ID = "Z", Fate = "Unspecified") %>%
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
  
  # Assign colors based on CV
  p.data$node_colors <- cut(p.data[[cv_column]],
                            breaks = CV_BREAKS,
                            labels = names(CV_COLORS),
                            right = FALSE)
  
  # Prepare data for vertical lines
  p.data.y <- p.data %>%
    select(CellName, ID, x, Start, End, node_colors, Fate) %>%
    mutate(Mid = (End + Start) / 2) %>%
    melt(id.vars = c("CellName", "ID", "x", "Mid", "node_colors", "Fate"),
         variable.name = "SE", value.name = "Posi")
  
  # Prepare data for horizontal lines
  p.data.x <- p.data %>%
    select(ID, x, Start, parentx) %>%
    melt(id.vars = c("ID", "Start"), variable.name = "FT", value.name = "Posi")
  
  # Create plot
  p <- p.data.y %>%
    mutate(label = if_else(nchar(ID) < ID_SHOW_LEN, CellName, NA_character_)) %>%
    ggplot() +
    # Vertical lines (cell lifespan)
    geom_path(aes(x = x, y = -Posi, group = ID), color = BG_COLOR, size = 0.3) +
    # Cell name labels
    geom_text(aes(x = x, y = -Mid, label = label)) +
    # Horizontal lines (connections)
    geom_path(data = p.data.x, aes(x = Posi, y = -Start, group = ID),
              color = BG_COLOR, size = 0.3) +
    # CV-colored points at cell birth
    geom_point(data = subset(p.data.y, SE == "Start" & !is.na(node_colors)),
               aes(x = x, y = -Posi, color = node_colors), size = 1.5) +
    scale_color_manual(values = CV_COLORS) +
    # Death markers (triangles)
    # Death markers (triangles)
    geom_point(data = subset(p.data.y, SE == "End" & Fate == "Death"),
               aes(x = x, y = -Posi),
               shape = 24, fill = "black", color = "black", size = 1) +
    
    # geom_point(data = subset(p.data.y, SE == "End" & Fate == "Death" &
    # #geom_point(data = subset(p.data.y,  Fate == "Death" &
    #                            !(Posi %in% c(281.8, 260.05, 372.3))),
    #            aes(x = x, y = -Posi),
    #            shape = 24, fill = "black", color = "black", size = 1) +
    # Axis settings
    scale_x_continuous(n.breaks = 20) +
    scale_y_continuous(limits = c(-400, 20), labels = function(x) abs(x)) +
    guides(color = "none") +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.y = element_line(color = "black")
    )
  
  # # 统计死亡细胞数量
  # death_count <- p.data.y %>%
  #   filter(SE == "End" & Fate == "Death") %>%
  #   nrow()
  # 
  # cat(sprintf("%s: %d death cells (triangles)\n", data_type, death_count))
  
  # 统计死亡细胞
  death_cells <- p.data.y %>%
    filter(SE == "End" & Fate == "Death")
  
  death_count <- nrow(death_cells)
  cat(sprintf("%s: %d death cells\n", data_type, death_count))
  
  # 返回plot和统计信息
  list(plot = p, death_count = death_count, death_cells = death_cells)
  
  
  #return(p)
  
}

# =============================================================================
# 5. Generate plots for each species
# =============================================================================

result_ce <- process_lineage_and_plot("CDFile/ce/CD200113plc1p3.csv", "cv_len_ce", 195, 1.44, "ce")
result_cb <- process_lineage_and_plot("CDFile/she1/CD240731cbhis72p1.csv", "cv_len_cb", 165, 1.57, "cb")
result_cn <- process_lineage_and_plot("CDFile/cn/CD241202cnhis72p1.csv", "cv_len_cn", 235, 1.58, "cn")

# 查看统计
cat("\nSummary:\n")
cat("CE:", result_ce$death_count, "triangles\n")
cat("CB:", result_cb$death_count, "triangles\n")
cat("CN:", result_cn$death_count, "triangles\n")

# 合并图
combined_plot <- result_ce$plot / result_cb$plot / result_cn$plot
print(result_ce$plot)
print(combined_plot)

# 查看死亡细胞名称
print(result_ce$death_cells %>% select(CellName, Posi) %>% distinct())
# 查看死亡细胞名称
print(result_cb$death_cells %>% select(CellName, Posi) %>% distinct())
# 查看死亡细胞名称
print(result_cn$death_cells %>% select(CellName, Posi) %>% distinct())
# 获取各物种死亡细胞名称
ce_death <- result_ce$death_cells %>% 
  filter(SE == "End") %>% 
  pull(CellName) %>% 
  unique()

cb_death <- result_cb$death_cells %>% 
  filter(SE == "End") %>% 
  pull(CellName) %>% 
  unique()

cn_death <- result_cn$death_cells %>% 
  filter(SE == "End") %>% 
  pull(CellName) %>% 
  unique()

# dead cells that cb don't have and ce have
cb_missing <- setdiff(ce_death, cb_death)
cat("\n=== CB missing (in CE but not in CB):", length(cb_missing), "===\n")
print(cb_missing)

# dead cells that cn don't have and ce have
cn_missing <- setdiff(ce_death, cn_death)
cat("\n=== CN missing (in CE but not in CN):", length(cn_missing), "===\n")
print(cn_missing)

# common dead cells among three nematodes
common_death <- Reduce(intersect, list(ce_death, cb_death, cn_death))
cat("\n=== Common death cells (all 3 species):", length(common_death), "===\n")
print(common_death)


# C. elegans
p_ce <- process_lineage_and_plot(
  "CDFile/ce/CD200113plc1p3.csv", "cv_len_ce", 195, 1.44, "ce"
)
print(p_ce)

# C. briggsae
p_cb <- process_lineage_and_plot(
  "CDFile/she1/CD240731cbhis72p1.csv", "cv_len_cb", 165, 1.57, "cb"
)
print(p_cb)

# C. nigoni
p_cn <- process_lineage_and_plot(
  "CDFile/cn/CD241202cnhis72p1.csv", "cv_len_cn", 235, 1.58, "cn"
)
print(p_cn)

# =============================================================================
# 6. Combine plots
# =============================================================================

combined_plot <- p_ce / p_cb / p_cn
print(combined_plot)

# Save plot (uncomment to use)
ggsave("all_tree.pdf", combined_plot, width = 10, height = 6.5)



# =============================================================================
#  PLOT_NORMALIZED FREQUENCY
# =============================================================================
# =============================================================================
# Nematode Cell Cycle Length - Normalized Cumulative Frequency Analysis
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Desktop/cbcn/")

# =============================================================================
# 1. Load data and calculate group means
# =============================================================================

data <- read.csv("nor_Length_Time_correct.csv") %>%
  mutate(
    ceLen_mean = round(rowMeans(select(., corLen_cep1:corLen_cep8), na.rm = TRUE), 2),
    cbLen_mean = round(rowMeans(select(., corLen_cbp1:corLen_cbp6), na.rm = TRUE), 2),
    cnLen_mean = round(rowMeans(select(., corLen_cnp1:corLen_cnp7), na.rm = TRUE), 2)
  ) %>%
  select(cell, ceLen_mean, cbLen_mean, cnLen_mean)

# =============================================================================
# 2. Transform to long format and calculate cumulative frequency
# =============================================================================

type_labels <- c(ceLen_mean = "ce", cbLen_mean = "cb", cnLen_mean = "cn")

data_long <- data %>%
  pivot_longer(
    cols = c(ceLen_mean, cbLen_mean, cnLen_mean),
    names_to = "type",
    values_to = "length"
  ) %>%
  mutate(type = recode(type, !!!type_labels)) %>%
  group_by(type) %>%
  arrange(length, .by_group = TRUE) %>%
  mutate(rel_cum_freq = row_number() / n() * 100) %>%
  ungroup()

# =============================================================================
# 3. Plot cumulative relative frequency
# =============================================================================

group_colors <- c(ce = "#CD534C", cb = "#0073C2", cn = "#EFC000")

p_3a <- ggplot(data_long, aes(x = length, y = rel_cum_freq, color = type)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = seq(0, 140, by = 30)) +
  labs(
    x = "Length",
    y = "Cumulative Relative Frequency (%)",
    title = "",
    color = "Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

print(p_3a)

# Save plot (uncomment to use)
# ggsave("fig3a.pdf", plot = p_3a, height = 5, width = 5)




# =============================================================================
#  PLOT_SIGNIFICANCE BAR PLOT
# =============================================================================
# =============================================================================
# Z-Score Analysis and Significance Bar Plot
#
# Description: Calculate z-scores for cb and cn relative to ce (template),
#              perform FDR correction, and visualize significant differences
# =============================================================================
library(ggplot2)
library(dplyr)

setwd("~/Desktop/cbcn/")

# Color scheme
GROUP_COLORS <- c(cb = "#0073C2", cn = "#EFC000")

# -----------------------------------------------------------------------------
# Data Loading and Mean Calculation
# -----------------------------------------------------------------------------

data <- read.csv("nor_Length_Time_correct.csv") %>%
  select(cell = 1, corLen_cep1:corLen_cep8, corLen_cbp1:corLen_cbp6,
         corLen_cnp1:corLen_cnp7, stage) %>%
  mutate(across(where(is.numeric), ~ round(., 1))) %>%
  filter(!is.na(stage)) %>%
  mutate(
    ce_mean = round(rowMeans(select(., corLen_cep1:corLen_cep8), na.rm = TRUE), 2),
    ce_sd   = round(apply(select(., corLen_cep1:corLen_cep8), 1, sd, na.rm = TRUE), 2),
    cb_mean = round(rowMeans(select(., corLen_cbp1:corLen_cbp6), na.rm = TRUE), 2),
    cn_mean = round(rowMeans(select(., corLen_cnp1:corLen_cnp7), na.rm = TRUE), 2)
  ) %>%
  select(cell, ce_mean, ce_sd, cb_mean, cn_mean, stage)

# -----------------------------------------------------------------------------
# Z-Score Calculation
# -----------------------------------------------------------------------------

# Calculate z-scores: (species_mean - ce_mean) / ce_sd
data <- data %>%
  mutate(
    zscore_cb = (cb_mean - ce_mean) / ce_sd,
    zscore_cn = (cn_mean - ce_mean) / ce_sd
  )

# -----------------------------------------------------------------------------
# Statistical Testing with FDR Correction
# -----------------------------------------------------------------------------

# Function to calculate p-value, q-value (FDR), and significance direction
calc_significance <- function(zscore) {
  # Two-tailed p-value
  p_val <- 2 * pnorm(-abs(zscore))
  
  # Direction: -1 if zscore < 0, 1 if zscore > 0
  
  p_dir <- case_when(
    is.na(zscore) ~ NA_real_,
    zscore < 0 ~ -1,
    zscore > 0 ~ 1,
    TRUE ~ NA_real_
  )
  
  # FDR-adjusted q-value
  q_val <- p.adjust(p_val, method = "fdr")
  
  # Significance classification: 1/-1 if q < 0.05, else 0
  q_bh <- case_when(
    is.na(q_val) ~ NA_real_,
    q_val < 0.05 & p_dir == 1 ~ 1,
    q_val < 0.05 & p_dir == -1 ~ -1,
    TRUE ~ 0
  )
  
  list(p_val = p_val, q_val = q_val, q_bh = q_bh, p_dir = p_dir)
}

# Apply to cb and cn
cb_stats <- calc_significance(data$zscore_cb)
cn_stats <- calc_significance(data$zscore_cn)

data <- data %>%
  mutate(
    q_bh_cb = cb_stats$q_bh,
    q_bh_cn = cn_stats$q_bh
  )

# -----------------------------------------------------------------------------
# Prepare Data for Bar Plot
# -----------------------------------------------------------------------------

# Count significant differences (q_bh = 1 or -1) for each species
plot_data <- data.frame(
  zscore = rep(c("1", "-1"), each = 2),
  type   = rep(c("cb", "cn"), times = 2),
  count  = c(
    sum(data$q_bh_cb == 1, na.rm = TRUE),   # cb longer than ce
    sum(data$q_bh_cn == 1, na.rm = TRUE),   # cn longer than ce
    sum(data$q_bh_cb == -1, na.rm = TRUE),  # cb shorter than ce
    sum(data$q_bh_cn == -1, na.rm = TRUE)   # cn shorter than ce
  )
)

# -----------------------------------------------------------------------------
# Generate Bar Plot
# -----------------------------------------------------------------------------

p_3b <- ggplot(plot_data, aes(x = zscore, y = count, fill = type)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.6),
           width = 0.6,
           color = "black") +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.5),
            vjust = -0.5,
            size = 5) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(title = "", x = "", y = "", fill = "") +
  scale_y_continuous(limits = c(0, 130), expand = c(0.003, 0.01)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

print(p_3b)
# ggsave("zscore_significance_barplot.pdf", plot = p_3b, width = 5, height = 5)




# =============================================================================
#  PLOT_CV BOX PLOT
# =============================================================================
# =============================================================================
# CV Boxplot with Pairwise Significance Test
#
# Description: Calculate coefficient of variation (CV) for normalized cell
#              cycle length, perform pairwise Wilcoxon tests, and visualize
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Desktop/cbcn/")

# Color scheme
GROUP_COLORS <- c(ce = "#CD534C", cb = "#0073C2", cn = "#EFC000")

# -----------------------------------------------------------------------------
# Data Loading and CV Calculation
# -----------------------------------------------------------------------------

data <- read.csv("nor_Length_Time_correct.csv") %>%
  mutate(
    ceLen_mean = round(rowMeans(select(., corLen_cep1:corLen_cep8), na.rm = TRUE), 2),
    cbLen_mean = round(rowMeans(select(., corLen_cbp1:corLen_cbp6), na.rm = TRUE), 2),
    cnLen_mean = round(rowMeans(select(., corLen_cnp1:corLen_cnp7), na.rm = TRUE), 2)
  ) %>%
  rowwise() %>%
  mutate(
    cv_len_ce = sd(c_across(corLen_cep1:corLen_cep8)) / ceLen_mean,
    cv_len_cb = sd(c_across(corLen_cbp1:corLen_cbp6)) / cbLen_mean,
    cv_len_cn = sd(c_across(corLen_cnp1:corLen_cnp7)) / cnLen_mean
  ) %>%
  ungroup() %>%
  select(cell, cv_len_ce, cv_len_cb, cv_len_cn, stage)

# Convert to long format
data_long <- data %>%
  pivot_longer(
    cols = starts_with("cv_len_"),
    names_to = "group",
    names_prefix = "cv_len_",
    values_to = "cv_length"
  ) %>%
  mutate(group = factor(group, levels = c("ce", "cb", "cn")))

# -----------------------------------------------------------------------------
# Pairwise Significance Test
# -----------------------------------------------------------------------------

# Perform pairwise Wilcoxon test with FDR correction
pairwise_results <- pairwise.wilcox.test(
  data_long$cv_length,
  data_long$group,
  p.adjust.method = "BH"
)

# Format results with significance symbols
format_pairwise_results <- function(pairwise_result) {
  p_matrix <- pairwise_result$p.value
  
  data.frame(
    comparison = c("ce vs cb", "ce vs cn", "cb vs cn"),
    p_value = c(p_matrix["cb", "ce"], p_matrix["cn", "ce"], p_matrix["cn", "cb"])
  ) %>%
    mutate(
      significance = case_when(
        p_value < 0.0001 ~ "****",
        p_value < 0.001  ~ "***",
        p_value < 0.01   ~ "**",
        p_value < 0.05   ~ "*",
        TRUE             ~ "ns"
      ),
      p_formatted = ifelse(p_value < 0.0001,
                           format(p_value, scientific = TRUE, digits = 2),
                           round(p_value, 4))
    )
}

results_table <- format_pairwise_results(pairwise_results)
print(results_table)

# -----------------------------------------------------------------------------
# CV Boxplot
# -----------------------------------------------------------------------------

p_3c <- ggplot(data_long, aes(x = group, y = cv_length, fill = group)) +
  stat_boxplot(geom = "errorbar", width = 0.42, size = 0.5) +
  geom_boxplot(width = 0.65,
               outlier.shape = NA,
               color = "black",
               size = 0.5) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(x = "Measurement", y = "Value", title = "") +
  scale_y_continuous(limits = c(0, 0.13), expand = c(0.01, 0)) +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

print(p_3c)
# ggsave("fig3c.pdf", plot = p_3c, height = 5, width = 5)





# =============================================================================
#  PLOT_Linear Regression Comparison Plot
# =============================================================================
# =============================================================================
# Linear Regression Comparison Plot
#
# Description: Compare division time of each embryo against ce template mean
#              with linear regression lines for each species group
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

setwd("~/Desktop/cbcn/")

# -----------------------------------------------------------------------------
# Data Loading and Template Calculation
# -----------------------------------------------------------------------------

data <- read.csv("Length_Time_Absolute.csv", header = TRUE)

# Calculate template mean from C. elegans
template <- data %>%
  mutate(
    temLen_mean = round(rowMeans(select(., Len_cep1:Len_cep8), na.rm = TRUE), 2),
    temTim_mean = round(rowMeans(select(., Tim_cep1:Tim_cep8), na.rm = TRUE), 2)
  ) %>%
  select(cell, temLen_mean, temTim_mean)

# -----------------------------------------------------------------------------
# Helper Function for Data Preparation and Plotting
# -----------------------------------------------------------------------------

create_regression_plot <- function(data, template, col_prefix, n_samples, title) {
  # Generate column names
  cols <- paste0("Tim_", col_prefix, 1:n_samples)
  
  # Reshape to long format and merge with template
  plot_data <- data %>%
    select(cell, all_of(cols)) %>%
    pivot_longer(
      cols = all_of(cols),
      names_to = "variable",
      values_to = "time_value"
    ) %>%
    left_join(template, by = "cell") %>%
    na.omit()
  
  # Create scatter plot with regression lines
  ggplot(plot_data, aes(x = temTim_mean, y = time_value, color = variable)) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = title, x = "Template Time Mean", y = "Time Value") +
    xlim(0, 300) +
    ylim(0, 300) +
    theme_bw() +
    theme(legend.title = element_blank())
}

# -----------------------------------------------------------------------------
# Generate Plots for Each Species Group
# -----------------------------------------------------------------------------

plot_ce <- create_regression_plot(data, template, "cep", 8, "CE Group")
plot_cb <- create_regression_plot(data, template, "cbp", 6, "CB Group")
plot_cn <- create_regression_plot(data, template, "cnp", 7, "CN Group")

# Combine plots horizontally
grid.arrange(plot_ce, plot_cb, plot_cn, nrow = 1)
# ggsave("regression_comparison.pdf", grid.arrange(plot_ce, plot_cb, plot_cn, nrow = 1), width = 15, height = 5)






# =============================================================================
#  PLOT_ZSCORE_HEATMAP
# =============================================================================
# =============================================================================
# Z-Score Significance Heatmap by Developmental Stage
#
# Description: Calculate z-scores for cb and cn relative to ce, perform FDR
#              correction, and visualize significance patterns across stages
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Desktop/cbcn/")

# Configuration
ORDER_STAGE <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                 "Stage_AB64", "Stage_AB128", "Stage_AB256")
ORDER_LETTER <- c("A", "M", "E", "C", "D", "P")
HEATMAP_COLORS <- c("1" = "#9c2831", "0" = "#f5efe7", "-1" = "navy")

# -----------------------------------------------------------------------------
# Data Loading and Mean/SD Calculation
# -----------------------------------------------------------------------------

data <- read.csv("nor_Length_Time_correct.csv") %>%
  select(cell = 1, corLen_cep1:corLen_cep8, corLen_cbp1:corLen_cbp6,
         corLen_cnp1:corLen_cnp7, stage) %>%
  mutate(across(where(is.numeric), ~ round(., 1))) %>%
  filter(!is.na(stage)) %>%
  mutate(
    ce_mean = round(rowMeans(select(., corLen_cep1:corLen_cep8), na.rm = TRUE), 2),
    ce_sd   = round(apply(select(., corLen_cep1:corLen_cep8), 1, sd, na.rm = TRUE), 2),
    cb_mean = round(rowMeans(select(., corLen_cbp1:corLen_cbp6), na.rm = TRUE), 2),
    cn_mean = round(rowMeans(select(., corLen_cnp1:corLen_cnp7), na.rm = TRUE), 2)
  ) %>%
  select(cell, ce_mean, ce_sd, cb_mean, cn_mean, stage)

# -----------------------------------------------------------------------------
# Z-Score and Significance Calculation
# -----------------------------------------------------------------------------

# Calculate z-scores
data <- data %>%
  mutate(
    zscore_cb = (cb_mean - ce_mean) / ce_sd,
    zscore_cn = (cn_mean - ce_mean) / ce_sd
  )

# Function to calculate significance with FDR correction
calc_qbh <- function(zscore, threshold = 0.01) {
  p_val <- 2 * pnorm(-abs(zscore))
  p_dir <- case_when(is.na(zscore) ~ NA_real_, zscore < 0 ~ -1, zscore > 0 ~ 1, TRUE ~ NA_real_)
  q_val <- p.adjust(p_val, method = "fdr")
  case_when(
    is.na(q_val) ~ NA_real_,
    q_val < threshold & p_dir == 1 ~ 1,
    q_val < threshold & p_dir == -1 ~ -1,
    TRUE ~ 0
  )
}

data <- data %>%
  mutate(
    q_bh_cb = calc_qbh(zscore_cb),
    q_bh_cn = calc_qbh(zscore_cn)
  )

# -----------------------------------------------------------------------------
# Data Preparation for Heatmap
# -----------------------------------------------------------------------------

data_plot <- data %>%
  select(cell, q_bh_cb, q_bh_cn, stage) %>%
  mutate(
    stage = factor(stage, levels = ORDER_STAGE),
    cell_letter = factor(substr(cell, 1, 1), levels = ORDER_LETTER),
    cell_nchar = nchar(cell)
  ) %>%
  arrange(stage, cell_letter, cell_nchar) %>%
  mutate(cell = factor(cell, levels = unique(cell)))

# Convert to long format
data_plot_long <- data_plot %>%
  pivot_longer(
    cols = c(q_bh_cb, q_bh_cn),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    value = factor(value, levels = c("-1", "0", "1")),
    variable = recode(variable, q_bh_cb = "q_cb", q_bh_cn = "q_cn"),
    variable = factor(variable, levels = c("q_cn", "q_cb"))
  )

# -----------------------------------------------------------------------------
# Calculate Summary Statistics for Facet Labels
# -----------------------------------------------------------------------------

summary_counts <- data_plot_long %>%
  group_by(stage, variable, value) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = value, values_from = n, values_fill = 0)

facet_summary <- summary_counts %>%
  group_by(stage) %>%
  summarize(
    cb_red   = sum(if_else(variable == "q_cb", `1`, 0L)),
    cb_blue  = sum(if_else(variable == "q_cb", `-1`, 0L)),
    cb_white = sum(if_else(variable == "q_cb", `0`, 0L)),
    cn_red   = sum(if_else(variable == "q_cn", `1`, 0L)),
    cn_blue  = sum(if_else(variable == "q_cn", `-1`, 0L)),
    cn_white = sum(if_else(variable == "q_cn", `0`, 0L)),
    facet_label = paste0("cb: ", cb_red, "/", cb_blue, "/", cb_white, "\n",
                         "cn: ", cn_red, "/", cn_blue, "/", cn_white)
  ) %>%
  ungroup()

stage_labels <- setNames(facet_summary$facet_label, facet_summary$stage)

# -----------------------------------------------------------------------------
# Generate Heatmap
# -----------------------------------------------------------------------------

p_heatmap <- ggplot(data_plot_long, aes(x = cell, y = variable, fill = value)) +
  geom_tile() +
  facet_grid(. ~ stage, scales = "free_x", space = "free",
             labeller = labeller(stage = stage_labels)) +
  scale_fill_manual(values = HEATMAP_COLORS, na.value = "gray", name = "Value") +
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

print(p_heatmap)
# ggsave("f3_sup5_qvalue_0.01.pdf", plot = p_heatmap, width = 10, height = 3, dpi = 300)





# =============================================================================
#  PLOT_LINEMAP AND COUNT
# =============================================================================
# =============================================================================
# Stage-wise Cell Cycle Length Comparison with Treemap
#
# Description: Compare normalized cell cycle length across species (ce, cb, cn)
#              with line plots and treemaps showing which species is longest
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(treemapify)
library(patchwork)

setwd("~/Desktop/cbcn/")

# Configuration
ORDER_STAGE <- c("Stage_AB4", "Stage_AB8", "Stage_AB16", "Stage_AB32",
                 "Stage_AB64", "Stage_AB128", "Stage_AB256")
ORDER_LETTER <- c("A", "M", "E", "C", "D", "P")
METRIC_COLORS <- c(ceLen_mean = "#CD534C", cbLen_mean = "#0073C2", cnLen_mean = "#EFC000")
GAP_SIZE <- 2  # Gap between stages on x-axis

# -----------------------------------------------------------------------------
# Data Loading and Mean Calculation
# -----------------------------------------------------------------------------

data <- read.csv("nor_Length_Time_correct.csv") %>%
  mutate(
    ceLen_mean = round(rowMeans(select(., corLen_cep1:corLen_cep8), na.rm = TRUE), 2),
    cbLen_mean = round(rowMeans(select(., corLen_cbp1:corLen_cbp6), na.rm = TRUE), 2),
    cnLen_mean = round(rowMeans(select(., corLen_cnp1:corLen_cnp7), na.rm = TRUE), 2)
  ) %>%
  select(cell, ceLen_mean, cbLen_mean, cnLen_mean, stage) %>%
  mutate(stage = factor(stage, levels = ORDER_STAGE))

# -----------------------------------------------------------------------------
# Convert to Long Format and Calculate X Coordinates
# -----------------------------------------------------------------------------

# Convert to long format with cell ordering
data_long <- data %>%
  pivot_longer(
    cols = c(ceLen_mean, cbLen_mean, cnLen_mean),
    names_to = "Metric",
    values_to = "Len_mean"
  ) %>%
  mutate(
    first_letter = factor(substr(cell, 1, 1), levels = ORDER_LETTER),
    cell_length = nchar(cell)
  ) %>%
  arrange(stage, first_letter, cell) %>%
  mutate(cell = factor(cell, levels = unique(cell)))

# Calculate new x coordinates with gaps between stages
cell_order <- data_long %>%
  distinct(cell, stage) %>%
  group_by(stage) %>%
  mutate(cell_index = row_number()) %>%
  ungroup()

stage_counts <- cell_order %>%
  group_by(stage) %>%
  summarize(count = n(), .groups = "drop") %>%
  arrange(stage) %>%
  mutate(cum = lag(cumsum(count), default = 0))

cell_order <- cell_order %>%
  left_join(stage_counts, by = "stage") %>%
  mutate(new_x = cum + cell_index + GAP_SIZE * (as.integer(stage) - 1))

data_long <- left_join(data_long, cell_order %>% select(cell, new_x), by = "cell")

# -----------------------------------------------------------------------------
# Calculate Fill Regions (Min/Mid/Max for Each Cell)
# -----------------------------------------------------------------------------

fill_df <- data_long %>%
  group_by(cell) %>%
  arrange(Len_mean) %>%
  summarize(
    min_value  = first(Len_mean),
    Metric_min = first(Metric),
    mid_value  = nth(Len_mean, 2),
    Metric_mid = nth(Metric, 2),
    max_value  = last(Len_mean),
    Metric_max = last(Metric),
    x = first(new_x),
    .groups = "drop"
  ) %>%
  mutate(
    fill_min = METRIC_COLORS[Metric_min],
    fill_mid = METRIC_COLORS[Metric_mid],
    fill_max = METRIC_COLORS[Metric_max]
  )

# Calculate stage borders
stage_border <- cell_order %>%
  group_by(stage) %>%
  summarize(
    xmin = min(new_x) - 0.5,
    xmax = max(new_x) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(width = xmax - xmin)

# -----------------------------------------------------------------------------
# Main Line Plot
# -----------------------------------------------------------------------------

p_main <- ggplot(data_long, aes(x = new_x, y = Len_mean,
                                group = interaction(Metric, stage), color = Metric)) +
  # Stage borders
  geom_rect(data = stage_border,
            aes(xmin = xmin, xmax = xmax),
            ymin = -Inf, ymax = Inf,
            fill = NA, color = "black", size = 0.3, inherit.aes = FALSE) +
  # Fill regions: bottom (min), middle, top (max)
  geom_rect(data = fill_df,
            aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 0, ymax = min_value, fill = fill_min),
            inherit.aes = FALSE, alpha = 0.4) +
  geom_rect(data = fill_df,
            aes(xmin = x - 0.5, xmax = x + 0.5, ymin = min_value, ymax = mid_value, fill = fill_mid),
            inherit.aes = FALSE, alpha = 0.4) +
  geom_rect(data = fill_df,
            aes(xmin = x - 0.5, xmax = x + 0.5, ymin = mid_value, ymax = max_value, fill = fill_max),
            inherit.aes = FALSE, alpha = 0.4) +
  # Lines
  geom_line(size = 0.4) +
  scale_color_manual(values = METRIC_COLORS) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0, 0), breaks = cell_order$new_x, labels = cell_order$cell) +
  labs(title = "", x = "Cell", y = "Len_mean") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# -----------------------------------------------------------------------------
# Treemap (Count of Longest Species per Stage)
# -----------------------------------------------------------------------------

# Calculate treemap data
tree_data <- fill_df %>%
  left_join(cell_order %>% select(cell, stage), by = "cell") %>%
  group_by(stage, Metric_max) %>%
  summarize(count = n(), .groups = "drop")

# Create individual treemaps for each stage
treemap_list <- map(ORDER_STAGE, function(stg) {
  df_sub <- tree_data %>% filter(stage == stg)
  ggplot(df_sub, aes(area = count, fill = Metric_max, label = paste(count))) +
    geom_treemap(alpha = 0.6) +
    geom_treemap_text(color = "white", place = "centre", grow = FALSE, size = 10, min.size = 0) +
    scale_fill_manual(values = METRIC_COLORS) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(1, 1, 1, 2.5)
    )
})

# Combine treemaps with widths proportional to stage cell counts
stage_widths <- setNames(stage_border$width, stage_border$stage)
p_tree_combined <- wrap_plots(treemap_list, nrow = 1, widths = stage_widths)

# -----------------------------------------------------------------------------
# Final Combined Plot
# -----------------------------------------------------------------------------

p_sup <- p_tree_combined / p_main + plot_layout(heights = c(0.5, 3))

print(p_sup)
# ggsave("f3_sup6.pdf", p_sup, width = 10, height = 3)





# =============================================================================
#  PLOT_GASTRULATION ORDER
# =============================================================================
# =============================================================================
# Consensus Cell Order Ranking (Median Rank across All Replicates)
# Cells grouped only when CE, CB, CN median ranks are ALL identical
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Desktop/cbcn/")

# =============================================================================
# Configuration
# =============================================================================

GAS1 <- c("Ea", "Ep", "Z2", "Z3", "Daa", "Dap", "Dpa", "Dpp", "Capaa", "Capap",
          "Cappa", "Cappp", "Cppaa", "Cppap", "Cpppa", "Cpppp",
          "MSaaaa", "MSaaap", "MSapaa", "MSapap", "MSappa", "MSappp",
          "MSpaaa", "MSpaap", "MSppaa", "MSppap", "MSpppa", "MSpppp",
          "MSaapaa", "MSaapap", "MSaappa", "MSaappp",
          "MSpapaa", "MSpapap", "MSpappa", "MSpappp",
          "ABalaap", "ABalapp", "ABalpppa", "ABalpppp", "ABarappp", "ABplppap", "ABprpaap",
          "ABalpaaap", "ABalpaapp", "ABaraapaa", "ABaraapap", "ABaraappa", "ABaraappp",
          "ABarapaap", "ABalpaaaaa", "ABalpaaaap", "ABalpaapaa", "ABalpaapap",
          "ABalpappaa", "ABalpapppa", "ABalpapppp", "ABaraaapaa", "ABaraaapap",
          "ABarapaaaa", "ABarapaaap", "ABarapapaa", "ABarapappa", "ABarapappp",
          "ABprpapppa", "ABprpapppp")

EGRESSION <- c("ABarpapa", "ABarpapp", "ABalaapp", "ABalappa", "ABalaapa", "ABarapppp", "ABarapppa", "ABpraappp", "ABpraappa", "ABalpppap", "ABarppppp", "ABarppppa", "ABalpppaa", "ABalppppp","ABalapppa", "ABalappap", "ABalppppa", "ABalapppp")

CE_FILES <- list(
  ce1 = "CDFile/ce/CD191108plc1p1.csv", ce2 = "CDFile/ce/CD200109plc1p1.csv",
  ce3 = "CDFile/ce/CD200113plc1p3.csv", ce4 = "CDFile/ce/CD200113plc1p2.csv",
  ce5 = "CDFile/ce/CD200322plc1p2.csv", ce6 = "CDFile/ce/CD200323plc1p1.csv",
  ce7 = "CDFile/ce/CD200326plc1p3.csv", ce8 = "CDFile/ce/CD200326plc1p4.csv")

CB_FILES <- list(
  cb1 = "CDFile/she1/CD240731cbhis72p1.csv", cb2 = "CDFile/she1/CD240731cbhis72p2.csv",
  cb3 = "CDFile/she1/CD240731cbhis72p3.csv", cb4 = "CDFile/she1/CD241202cbhis72p1.csv",
  cb5 = "CDFile/she1/CD241202cbhis72p2.csv", cb6 = "CDFile/she1/CD241202cbhis72p4.csv")

CN_FILES <- list(
  cn1 = "CDFile/cn/CD241202cnhis72p1.csv", cn2 = "CDFile/cn/CD240712cnhis72p1.csv",
  cn3 = "CDFile/cn/CD240712cnhis72p2.csv", cn4 = "CDFile/cn/CD240712cnhis72p3.csv",
  cn5 = "CDFile/cn/CD241207cnhis72p1.csv", cn6 = "CDFile/cn/CD241207cnhis72p3.csv",
  cn7 = "CDFile/cn/CD241202cnhis72p2.csv")

# =============================================================================
# Functions
# =============================================================================

calc_median_rank <- function(file_list, cell_list) {
  rank_list <- lapply(names(file_list), function(name) {
    read.csv(file_list[[name]], header = TRUE) %>%
      select(1:3) %>%
      group_by(cell) %>%
      summarise(firstTime = min(time, na.rm = TRUE), .groups = "drop") %>%
      filter(cell %in% cell_list) %>%
      mutate(rank = dense_rank(firstTime)) %>%
      select(cell, rank) %>%
      rename(!!name := rank)
  })
  
  result <- data.frame(cell = cell_list, stringsAsFactors = FALSE)
  for (df in rank_list) result <- left_join(result, df, by = "cell")
  
  result %>%
    rowwise() %>%
    mutate(median_rank = median(c_across(-cell), na.rm = TRUE)) %>%
    ungroup() %>%
    select(cell, median_rank)
}


run_analysis <- function(cell_list, name) {
  cat(sprintf("\n========== %s (%d cells) ==========\n", name, length(cell_list)))
  
  # Calculate median ranks
  ce <- calc_median_rank(CE_FILES, cell_list) %>% rename(ce_median = median_rank)
  cb <- calc_median_rank(CB_FILES, cell_list) %>% rename(cb_median = median_rank)
  cn <- calc_median_rank(CN_FILES, cell_list) %>% rename(cn_median = median_rank)
  
  # Combined table
  combined <- ce %>%
    left_join(cb, by = "cell") %>%
    left_join(cn, by = "cell") %>%
    arrange(ce_median)
  
  cat("\nCombined Table:\n")
  print(combined, n = nrow(combined))
  
  # Grouped table
  grouped <- combined %>%
    group_by(ce_median, cb_median, cn_median) %>%
    summarise(cell_group = paste(cell, collapse = "/"), n_cells = n(), .groups = "drop") %>%
    arrange(ce_median, cb_median, cn_median) %>%
    mutate(order = row_number())
  
  cat("\nGrouped Table:\n")
  print(grouped, n = nrow(grouped))
  
  # 创建颜色调色板
  n_groups <- nrow(grouped)
  color_palette <- colorRampPalette(c("#F2F7FB", "#9CB6DD","darkblue"))(n_groups)
  
  # Heatmap
  df_plot <- grouped %>%
    pivot_longer(c(ce_median, cb_median, cn_median), names_to = "species", values_to = "median_rank") %>%
    group_by(species) %>%
    arrange(species, median_rank) %>%
    mutate(
      color_rank = row_number(),
      fill_color = color_palette[color_rank]
    ) %>%
    ungroup() %>%
    mutate(
      species = recode(species, ce_median = "cel", cb_median = "cbr", cn_median = "cni"),
      species = factor(species, levels = c("cel", "cbr", "cni")),
      cell_group = factor(cell_group, levels = grouped$cell_group))
  
  p <- ggplot(df_plot, aes(x = order, y = species, fill = color_rank)) +
    geom_tile(color = "black", size = 0.1) +
    scale_fill_gradientn(
      colors = c("#F2F7FB", "#9CB6DD", "darkblue"),
      name = "Order Rank"
    ) +
    scale_x_continuous(breaks = grouped$order, labels = grouped$cell_group, expand = c(0, 0)) +
    labs(x = "Cell Group (Ordered by CE)", y = "Species", title = paste0(name, " - Cell Order Comparison")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.position = "right")  # 添加这一行
  # p <- ggplot(df_plot, aes(x = order, y = species, fill = fill_color)) +
  #   geom_tile(color = "black", size = 0.1) +
  #   scale_fill_identity() +
  #   scale_x_continuous(breaks = grouped$order, labels = grouped$cell_group, expand = c(0, 0)) +
  #   labs(x = "Cell Group (Ordered by CE)", y = "Species", title = paste0(name, " - Cell Order Comparison")) +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
  #         panel.grid = element_blank(),
  #         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  
  print(p)
  list(combined = combined, grouped = grouped, plot = p)
}

# =============================================================================
# Run Analysis
# =============================================================================

results_gas1 <- run_analysis(GAS1, "Gas1")
results_egression <- run_analysis(EGRESSION, "Egression")

# Save results (uncomment to use)
# write.csv(results_gas1$combined, "gas1_median_rank_table.csv", row.names = FALSE)
# write.csv(results_gas1$grouped, "gas1_grouped_table.csv", row.names = FALSE)
#ggsave("gas1_heatmap.pdf", results_gas1$plot, width = 14, height = 4)

# write.csv(results_egression$combined, "egression_median_rank_table.csv", row.names = FALSE)
# write.csv(results_egression$grouped, "egression_grouped_table.csv", row.names = FALSE)
#ggsave("egression_heatmap.pdf", results_egression$plot, width = 8, height = 4)





