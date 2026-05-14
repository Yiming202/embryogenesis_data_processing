library(dplyr)
library(tidyr)
library(DESeq2)
library(tibble)
library(tidyplots)
library(org.Ce.eg.db)
library(patchwork)
library(ggplot2)

options(scipen = 999)
# Define a helper function to read the file and rename the 7th column
read_counts <- function(file_path) {
  df <- read.table(file_path,
                   header = TRUE,
                   sep = "\t",
                   comment.char = "#",
                   stringsAsFactors = FALSE)
  colnames(df)[7] <- "Counts"
  return(df)
}

# Specify the file paths
file1 <- "/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/counts/AF_a_counts.txt"
file2 <- "/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/counts/AF_b_counts.txt"
file3 <- "/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/counts/AF_c_counts.txt"
file4 <- "/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/counts/Ju_a_counts.txt"
file5 <- "/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/counts/Ju_b_counts.txt"
file6 <- "/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/counts/Ju_c_counts.txt"
file7 <- "/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/counts/Ce_a_counts.txt"
file8 <- "/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/counts/Ce_b_counts.txt"

# Use the helper function to read in each file
df_AF16_rep_1 <- read_counts(file1)
df_AF16_rep_2 <- read_counts(file2)
df_AF16_rep_3 <- read_counts(file3)
df_JU1421_rep_1 <- read_counts(file4)
df_JU1421_rep_2 <- read_counts(file5)
df_JU1421_rep_3 <- read_counts(file6)
df_Ce_rep_1 <- read_counts(file7) %>% mutate(Geneid = sub("^Gene:", "", Geneid))
df_Ce_rep_2 <- read_counts(file8) %>% mutate(Geneid = sub("^Gene:", "", Geneid))

columns(org.Ce.eg.db)
ce_map <- AnnotationDbi::select(
  org.Ce.eg.db,
  keys    = keys(org.Ce.eg.db, keytype = "SYMBOL"),   # all gene symbols
  columns = c("SYMBOL", "WORMBASE"),                  # <- correct name
  keytype = "SYMBOL"
)

df_mutual_cn_cb <- read.table("/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/mutual_best_hit/AF16_JU1421_mutual_best_hits.txt",
                              header = FALSE,
                              sep = "\t",
                              stringsAsFactors = FALSE,
                              col.names = c("Cbr_gene_name", "Cni_gene_name"))

df_mutual_cb_ce <- read.table("/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/mutual_best_hit/AF16_c_elegans_mutual_best_hits.txt",
                              header = FALSE,
                              sep = "\t",
                              stringsAsFactors = FALSE,
                              col.names = c("Cbr_gene_name", "Cel_gene_name"))

df_mutual_cn_ce <- read.table("/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/mutual_best_hit/JU1421_c_elegans_mutual_best_hits.txt",
                              header = FALSE,
                              sep = "\t",
                              stringsAsFactors = FALSE,
                              col.names = c("Cni_gene_name", "Cel_gene_name"))

length(unique(df_mutual_cn_cb[[1]]))

transcript2gene <- read.table("/Users/jeffrey/Desktop/Yiming_hybrid_lineaging_MS/RNAseq_analysis/mutual_best_hit/transcript2gene.tsv",
                              header = FALSE,
                              sep = "\t",
                              stringsAsFactors = FALSE,
                              col.names = c("transcript", "gene"))

df_mutual_cb_ce <- df_mutual_cb_ce %>%
  # Remove the "Transcript:" prefix
  mutate(Cel_gene_name = sub("^Transcript:", "", Cel_gene_name)) %>%
  # Join with the transcript2gene table to get the corresponding gene ID
  left_join(transcript2gene, by = c("Cel_gene_name" = "transcript")) %>%
  # Replace the Cel_gene_name with the gene column obtained from the join
  mutate(Cel_gene_name = gene) %>%
  # Remove the extra column from the join
  select(-gene)

# For df_mutual_cn_ce: Update the Cel_gene_name column in a similar way
df_mutual_cn_ce <- df_mutual_cn_ce %>%
  mutate(Cel_gene_name = sub("^Transcript:", "", Cel_gene_name)) %>%
  left_join(transcript2gene, by = c("Cel_gene_name" = "transcript")) %>%
  mutate(Cel_gene_name = gene) %>%
  select(-gene)

df_cn_cb_count_table <- df_mutual_cn_cb %>%
  # Join replicates for Cbr (AF16)
  left_join(
    df_AF16_rep_1 %>%
      dplyr::select(Geneid, Counts) %>%
      dplyr::rename(Cbr_counts_rep1 = Counts),
    by = c("Cbr_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_AF16_rep_2 %>%
      dplyr::select(Geneid, Counts) %>%
      dplyr::rename(Cbr_counts_rep2 = Counts),
    by = c("Cbr_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_AF16_rep_3 %>%
      dplyr::select(Geneid, Counts) %>%
      dplyr::rename(Cbr_counts_rep3 = Counts),
    by = c("Cbr_gene_name" = "Geneid")
  ) %>%
  # Join replicates for Cni (JU1421)
  left_join(
    df_JU1421_rep_1 %>%
      dplyr::select(Geneid, Counts) %>%
      dplyr::rename(Cni_counts_rep1 = Counts),
    by = c("Cni_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_JU1421_rep_2 %>%
      dplyr::select(Geneid, Counts) %>%
      dplyr::rename(Cni_counts_rep2 = Counts),
    by = c("Cni_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_JU1421_rep_3 %>%
      dplyr::select(Geneid, Counts) %>%
      dplyr::rename(Cni_counts_rep3 = Counts),
    by = c("Cni_gene_name" = "Geneid")
  )
# Include the gene names when selecting columns
count_matrix <- df_cn_cb_count_table %>%
  dplyr::select(Cbr_gene_name, Cbr_counts_rep1, Cbr_counts_rep2, Cbr_counts_rep3,
                Cni_counts_rep1, Cni_counts_rep2, Cni_counts_rep3)

# Set the gene names as row names and remove the gene name column from the data matrix
mat <- as.data.frame(count_matrix)
rownames(mat) <- mat$Cbr_gene_name  # set row names
mat <- mat[, -1]                   # drop the now redundant gene name column

# Create the sample info object as before
sample_info <- data.frame(
  condition = rep(c("CBR", "CNI"), each = 3)
)
rownames(sample_info) <- colnames(mat)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = as.matrix(mat),
                              colData = sample_info,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

# Retrieve normalized counts (rownames are gene names)
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)

# (Optional) Turn rownames into an explicit column if you want:
norm_counts_df <- rownames_to_column(norm_counts_df, var = "Cbr_gene_name")
norm_counts_df <- left_join(norm_counts_df,
                            df_cn_cb_count_table %>% 
                            dplyr::select(Cbr_gene_name, Cni_gene_name), by = "Cbr_gene_name")%>% 
                            dplyr::select(Cbr_gene_name, Cni_gene_name, everything())%>%
  mutate(
    Mean_CBR = rowMeans(dplyr::select(., starts_with("Cbr_counts")), na.rm = TRUE),
    Mean_CNI = rowMeans(dplyr::select(., starts_with("Cni_counts")), na.rm = TRUE)
  )

res <- results(dds)
res_df <- as.data.frame(res) %>% rownames_to_column(var = "Cbr_gene_name")

# Here, we join the DESeq2 results with our normalized counts table
# Then, we add a new column to flag significantly deviated genes based on an adjusted p-value threshold (e.g., < 0.05)
norm_counts_df <- norm_counts_df %>%
  left_join(res_df %>% dplyr::select(Cbr_gene_name, padj), by = "Cbr_gene_name") %>%
  mutate(Significance = if_else(!is.na(padj) & padj < 0.05, "Significant", "Not Significant")) %>%
  dplyr::select(-padj) 

norm_counts_df <- norm_counts_df %>%
  mutate(
    log2_Mean_CBR = log2(Mean_CBR + 1),
    log2_Mean_CNI = log2(Mean_CNI + 1)
  )

norm_counts_df <- norm_counts_df %>%
  left_join(
    df_mutual_cb_ce %>%
      dplyr::select(Cbr_gene_name, Cel_gene_name) %>%
      distinct(),
    by = "Cbr_gene_name"
  ) %>%
  relocate(Cel_gene_name, .after = Cni_gene_name) %>%
  mutate(elegans_gene_name = if_else(
    !is.na(Cel_gene_name),
    # Look up the SYMBOL that corresponds to Cel_gene_name by matching to ce_map$WORMBASE.
    as.character(ce_map$SYMBOL[match(Cel_gene_name, ce_map$WORMBASE)]),
    NA_character_
  ))%>%
  relocate(elegans_gene_name, .after = Cel_gene_name)

# Build cb_ce_count_table using the mutual-best hit info (CBR vs. CE)
df_cb_ce_count_table <- df_mutual_cb_ce %>%
  # Join replicates for CBR (AF16)
  left_join(
    df_AF16_rep_1 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cbr_counts_rep1 = Counts),
    by = c("Cbr_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_AF16_rep_2 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cbr_counts_rep2 = Counts),
    by = c("Cbr_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_AF16_rep_3 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cbr_counts_rep3 = Counts),
    by = c("Cbr_gene_name" = "Geneid")
  ) %>%
  # Join replicates for Cel (C. elegans)
  left_join(
    df_Ce_rep_1 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cel_counts_rep1 = Counts),
    by = c("Cel_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_Ce_rep_2 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cel_counts_rep2 = Counts),
    by = c("Cel_gene_name" = "Geneid")
  )

# Create count matrix for CB vs. CE
count_matrix_cb_ce <- df_cb_ce_count_table %>%
  select(Cbr_gene_name, Cbr_counts_rep1, Cbr_counts_rep2, Cbr_counts_rep3,
         Cel_counts_rep1, Cel_counts_rep2)

mat_cb_ce <- as.data.frame(count_matrix_cb_ce)
rownames(mat_cb_ce) <- mat_cb_ce$Cbr_gene_name
mat_cb_ce <- mat_cb_ce[, -1]

# Create sample information for CB vs. CE (3 AF16 replicates vs. 2 Ce replicates)
sample_info_cb_ce <- data.frame(
  condition = c(rep("CBR", 3), rep("CE", 2))
)
rownames(sample_info_cb_ce) <- colnames(mat_cb_ce)

# Build DESeqDataSet and run DE analysis for CB vs. CE
dds_cb_ce <- DESeqDataSetFromMatrix(countData = as.matrix(mat_cb_ce),
                                    colData = sample_info_cb_ce,
                                    design = ~ condition)
dds_cb_ce <- dds_cb_ce[rowSums(counts(dds_cb_ce)) > 1, ]
dds_cb_ce <- DESeq(dds_cb_ce)

# Retrieve normalized counts and DE results
norm_counts_cb_ce <- counts(dds_cb_ce, normalized = TRUE)
norm_counts_cb_ce_df <- as.data.frame(norm_counts_cb_ce) %>%
  rownames_to_column(var = "Cbr_gene_name")
res_cb_ce <- results(dds_cb_ce)
res_cb_ce_df <- as.data.frame(res_cb_ce) %>% rownames_to_column(var = "Cbr_gene_name")

# Merge and add significance flag and compute means/log2 transformation
norm_counts_cb_ce_df <- norm_counts_cb_ce_df %>%
  left_join(res_cb_ce_df %>% dplyr::select(Cbr_gene_name, padj), by = "Cbr_gene_name") %>%
  mutate(Significance = if_else(!is.na(padj) & padj < 0.05, "Significant", "Not Significant")) %>%
  dplyr::select(-padj) %>%
  mutate(
    Mean_CBR = rowMeans(select(., starts_with("Cbr_counts")), na.rm = TRUE),
    Mean_Cel = rowMeans(select(., starts_with("Cel_counts")), na.rm = TRUE),
    log2_Mean_CBR = log2(Mean_CBR + 1),
    log2_Mean_Cel = log2(Mean_Cel + 1)
  )

# Build cn_ce_count_table using the mutual-best hit info (CNI vs. CE)
df_cn_ce_count_table <- df_mutual_cn_ce %>%
  # Join replicates for CNI (JU1421)
  left_join(
    df_JU1421_rep_1 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cni_counts_rep1 = Counts),
    by = c("Cni_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_JU1421_rep_2 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cni_counts_rep2 = Counts),
    by = c("Cni_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_JU1421_rep_3 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cni_counts_rep3 = Counts),
    by = c("Cni_gene_name" = "Geneid")
  ) %>%
  # Join replicates for Cel (C. elegans)
  left_join(
    df_Ce_rep_1 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cel_counts_rep1 = Counts),
    by = c("Cel_gene_name" = "Geneid")
  ) %>%
  left_join(
    df_Ce_rep_2 %>% dplyr::select(Geneid, Counts) %>% dplyr::rename(Cel_counts_rep2 = Counts),
    by = c("Cel_gene_name" = "Geneid")
  )

# Create count matrix for CN vs. CE
count_matrix_cn_ce <- df_cn_ce_count_table %>%
  dplyr::select(Cni_gene_name, Cni_counts_rep1, Cni_counts_rep2, Cni_counts_rep3,
                Cel_counts_rep1, Cel_counts_rep2)

mat_cn_ce <- as.data.frame(count_matrix_cn_ce)
rownames(mat_cn_ce) <- mat_cn_ce$Cni_gene_name
mat_cn_ce <- mat_cn_ce[, -1]

# Create sample information for CN vs. CE (3 JU1421 replicates vs. 2 Ce replicates)
sample_info_cn_ce <- data.frame(
  condition = c(rep("CNI", 3), rep("CE", 2))
)
rownames(sample_info_cn_ce) <- colnames(mat_cn_ce)

# Build DESeqDataSet and run DE analysis for CN vs. CE
dds_cn_ce <- DESeqDataSetFromMatrix(countData = as.matrix(mat_cn_ce),
                                    colData = sample_info_cn_ce,
                                    design = ~ condition)
dds_cn_ce <- dds_cn_ce[rowSums(counts(dds_cn_ce)) > 1, ]
dds_cn_ce <- DESeq(dds_cn_ce)

# Retrieve normalized counts and DE results
norm_counts_cn_ce <- counts(dds_cn_ce, normalized = TRUE)
norm_counts_cn_ce_df <- as.data.frame(norm_counts_cn_ce) %>%
  rownames_to_column(var = "Cni_gene_name")
res_cn_ce <- results(dds_cn_ce)
res_cn_ce_df <- as.data.frame(res_cn_ce) %>% rownames_to_column(var = "Cni_gene_name")

norm_counts_cn_ce_df <- norm_counts_cn_ce_df %>%
  left_join(res_cn_ce_df %>% dplyr::select(Cni_gene_name, padj), by = "Cni_gene_name") %>%
  mutate(Significance = if_else(!is.na(padj) & padj < 0.05, "Significant", "Not Significant")) %>%
  dplyr::select(-padj) %>%
  mutate(
    Mean_CNI = rowMeans(select(., starts_with("Cni_counts")), na.rm = TRUE),
    Mean_Cel = rowMeans(select(., starts_with("Cel_counts")), na.rm = TRUE),
    log2_Mean_CNI = log2(Mean_CNI + 1),
    log2_Mean_Cel = log2(Mean_Cel + 1)
  ) %>%
  # Add the Cel gene name back by joining with the original count table
  left_join(
    df_cn_ce_count_table %>%
      dplyr::select(Cni_gene_name, Cel_gene_name) %>%
      distinct(),
    by = "Cni_gene_name"
  ) %>%
  # Move Cel_gene_name to the second column (after Cni_gene_name)
  relocate(Cel_gene_name, .after = Cni_gene_name) %>%
  mutate(elegans_gene_name = if_else(
    !is.na(Cel_gene_name),
    # Look up the SYMBOL that corresponds to Cel_gene_name by matching to ce_map$WORMBASE.
    as.character(ce_map$SYMBOL[match(Cel_gene_name, ce_map$WORMBASE)]),
    NA_character_
  ))%>%
  relocate(elegans_gene_name, .after = Cel_gene_name)

genes <- c("hsp-90", "mys-1", "trr-1", "dpy-22",
           "hmg-1.2", "din-1", "egl-27")
Cn_higher_gene <- c("din-1", "dpy-22")


genetic_modifier_genes <- c(
  "aph-1", "car-1", "cdc-37", "cdc-42", "emb-30", "fat-2", "lag-1", "lsy-22",
  "mel-26", "mel-28", "mel-32", "mex-3", "mom-2", "mom-5", "nmy-2", "par-1",
  "par-2", "par-3", "par-4", "par-5", "par-6", "pkc-3", "pos-1", "rfc-3",
  "rpn-9", "rpn-10", "skn-1", "skr-2", "sur-6"
)
Cn_higher_genetic_modifier_genes <- c("mex-3", "par-1","par-4", "rpn-10")
Cn_higher_cb_genetic_modifier_genes <- c("mex-3", "par-3","fat-2", "pos-1", "aph-1")

df_test <-norm_counts_df %>%
  filter(elegans_gene_name %in% genes)
df_test <-norm_counts_df %>%
  filter(elegans_gene_name %in% genetic_modifier_genes)

"#FBB4AE" "#B3CDE3" "#CCEBC5" "#DECBE4" "#FED9A6" "#FFFFCC"
"#E5D8BD" "#FDDAEC" "#F2F2F2"

# Build the scatter plot with log-transformed mean read counts and gene labeling.
scatter_plot <- norm_counts_df |>
  tidyplot(x = log2_Mean_CBR, y = log2_Mean_CNI, dodge_width = 0) |>
  # Overlay unsignificant genes in grey
  add_data_points(
    data = filter_rows(Significance == "Not Significant"),
    rasterize = FALSE,
    size = 0.5,
    color = "#F2F2F2"
  ) |>
  # Layer for significant genes (default color)
  add_data_points(
    data = filter_rows(Significance == "Significant"),
    rasterize = FALSE,
    size = 0.5,
    color = "#B3CDE3"
  ) |>
  # Finally, add the layer for gene "CBR06396" so it's drawn on top
  add_data_points(
    data = filter_rows(elegans_gene_name %in% Cn_higher_cb_genetic_modifier_genes),
    rasterize = FALSE,
    size = 1,
    color = "red"  # or any color to highlight the gene
  ) |>
  # (Optional) Add a label for the gene if desired
  add_data_labels_repel(
    data = filter_rows(elegans_gene_name %in% Cn_higher_cb_genetic_modifier_genes),
    label = elegans_gene_name,
    color = "black"
  ) |>
  remove_padding() |>
  adjust_title("C. nigoni vs. C. briggsae")|>
  adjust_x_axis_title("C. briggsae log2 (mean reads count)") |>
  adjust_y_axis_title("C. nigoni log2 (mean reads count)") |>
  # Set x-axis breaks and labels at 0, 5, 10, 15
  add(ggplot2::scale_x_continuous(breaks = c(0, 5, 10, 15))) |>
  # Set y-axis breaks and labels at 0, 5, 10, 15
  add(ggplot2::scale_y_continuous(breaks = c(0, 5, 10, 15)))|>
  adjust_legend_title("My new legend title")

print(scatter_plot)

ggsave(filename = "scatter_Cn_Cb.pdf", 
       plot = scatter_plot, 
       device = "pdf", 
       width = 4, 
       height = 4)

df_test <-norm_counts_cn_ce_df %>%
  filter(elegans_gene_name %in% genes)
df_test <-norm_counts_cn_ce_df %>%
  filter(elegans_gene_name %in% genetic_modifier_genes)

scatter_cn_ce <- norm_counts_cn_ce_df |>
  tidyplot(
    x = log2_Mean_Cel,                # x-axis: C. nigoni
    y = log2_Mean_CNI,                 # y-axis: C. elegans
    dodge_width = 0                   # no horizontal jitter
  ) |>
  add_data_points(
    data       = filter_rows(Significance == "Not Significant"),
    color      = "#F2F2F2",
    size       = 0.5,
    rasterize  = FALSE
  ) |>
  add_data_points(
    data       = filter_rows(Significance == "Significant"),
    color      = "#B3CDE3",
    size       = 0.5,
    rasterize  = FALSE
  ) |>
  add_data_points(
    data       = filter_rows(elegans_gene_name %in% Cn_higher_gene),
    color      = "orange",
    size       = 1,
    rasterize  = FALSE
  ) |>
  add_data_points(
    data       = filter_rows(elegans_gene_name %in% Cn_higher_genetic_modifier_genes),
    color      = "red",
    size       = 1,
    rasterize  = FALSE
  ) |>
  add_data_labels_repel(
    data   = filter_rows(elegans_gene_name %in% Cn_higher_gene|elegans_gene_name %in% Cn_higher_genetic_modifier_genes),
    label  = elegans_gene_name,
    color  = "black",
    size   = 3
  ) |>
  remove_padding() |>
  adjust_title("C. nigoni vs. C. elegans")|>
  adjust_x_axis_title("C. elegans log2 (mean reads count)") |>
  adjust_y_axis_title("C. nigoni log2 (mean reads count)") 

print(scatter_cn_ce)

ggsave(filename = "scatter_Cn_Ce.pdf", 
       plot = scatter_cn_ce, 
       device = "pdf", 
       width = 4, 
       height = 4)

combined_plot <- wrap_elements(full = scatter_cn_ce) + wrap_elements(full = scatter_plot) 
  
print(combined_plot)

ggsave(filename = "scatter_combined.pdf", 
       plot = combined_plot, 
       device = "pdf", 
       width = 6, 
       height = 3)


