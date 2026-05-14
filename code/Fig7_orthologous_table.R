library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(org.Ce.eg.db)
library(patchwork)

# Othorgroup analysis among the three species
# 1. Read the data while keeping the "Orthogroup" column along with numeric columns
df_data <- read.table(
  file = "/Users/jeffrey/Desktop/genome_file/orthogroup/Final/Intra_all_Orthogroups.GeneCount.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  dplyr::select(Orthogroup, JU1421, AF16, N2) %>%
  
  # Remove rows where all numeric columns are 0
  filter(!if_all(c(JU1421, AF16, N2), ~ . == 0)) %>%
  
  # Remove rows where all numeric columns are 1
  filter(!if_all(c(JU1421, AF16, N2), ~ . == 1)) %>%
  
  # Keep rows where:
  # (sum >= 3) OR (exactly one value is 2 and the others are 0)
  filter(
    rowSums(across(c(JU1421, AF16, N2))) >= 3 |
      (rowSums(across(c(JU1421, AF16, N2), ~ . == 2)) == 1 &
         rowSums(across(c(JU1421, AF16, N2), ~ . == 0)) == 2)
  ) %>%
  mutate(RowID = row_number())

df_data %>%
  summarise(total_N2 = sum(N2, na.rm = TRUE))

df_data %>%
  summarise(total_AF16 = sum(AF16, na.rm = TRUE))

df_data %>%
  summarise(total_JU1421 = sum(JU1421, na.rm = TRUE))

# View the filtered data with the original Orthogroup column and new RowID column
print(df_data)

df_data %>%
  summarise(total = sum(N2))

# Count the rows with JU1421 == 0, AF16 == 0 and N2 > 0
n_rows <- df_data %>%
  filter(JU1421 == 0, AF16 == 0, N2 > 0) %>%
  nrow()

# Print the number of rows meeting the criteria
print(n_rows)

# Convert data from wide to long format for all three species, exclude rows where the gene count is 0,
# and bin the gene counts.
df_long <- df_data %>%
  pivot_longer(
    cols = c(JU1421, AF16, N2),
    names_to = "Species",
    values_to = "GeneCount"
  ) %>%
  # Keep only GeneCount > 1
  filter(GeneCount > 1) %>%
  mutate(
    Bin = cut(GeneCount,
              breaks = c(2, 5, 10, 100, Inf),
              labels = c("2-5", "5-10", "10-100", ">100"),
              include.lowest = TRUE,
              right = FALSE),
    # Reorder species
    Species = factor(Species, levels = c("N2", "AF16", "JU1421"))
  )


# Create a 100% stacked bar plot with reversed stacking order, custom theme,
# and numeric labels (number of orthogroups) for each segment.
stacked_bar <- ggplot(df_long, aes(x = Species, fill = Bin)) +
  geom_bar(width = 0.7, 
           position = position_fill(reverse = TRUE), 
           color = NA,  # Removed the white borders so segments connect directly.
           size = 0.7) +
  # Label each segment with the number of orthogroups (rows) that fall into that bin
  geom_text(stat = "count", 
            aes(label = paste0("n=", ..count..)), 
            position = position_fill(vjust = 0.5, reverse = TRUE), 
            color = "black", size = 3) +
  scale_y_continuous(
    limits = c(0, 1.2),  # Extend the y-axis to 120% (1.2)
    breaks = c(0, 0.25, 0.50, 0.75, 1),  # Tick marks at 0%, 25%, 50%, 75%, and 100%
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Percentage of Orthogroups per Gene Count Range for Each Species",
    x = "Species",
    y = "Percentage",
    fill = "Gene Count Range"
  ) +
  scale_fill_manual(
    values = c(
      "2-5"    = "#E69F00",  # custom color for 1-5
      "5-10"   = "#56B4E9",  # custom color for 5-10
      "10-100" = "#009E73",  # custom color for 10-100
      ">100"   = "#F0E442"   # custom color for >100
    ),
    guide = guide_legend(reverse = TRUE)  # Reverse the order in the legend
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Graph border
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.ticks.x = element_line(color = "black"), # Add x-axis ticks
    axis.ticks.y = element_line(color = "black")  # Add y-axis ticks
  )

# Display the plot
print(stacked_bar)

ggsave(filename = "Orthogroup_bar_percentage.pdf", 
       plot = stacked_bar, 
       device = "pdf", 
       width = 6, 
       height = 4)

# 2. Determine for each row which strain has the maximum value using the new RowID for reference.

df_data_with_max <- df_data %>%
  rowwise() %>%
  mutate(
    # Gather the values for the three strains in the current row
    max_value = max(c_across(c(JU1421, AF16, N2))),
    num_max = sum(c_across(c(JU1421, AF16, N2)) == max_value),
    MaxStrain = ifelse(
      num_max == 1,
      c("JU1421", "AF16", "N2")[which.max(c_across(c(JU1421, AF16, N2)))],
      NA_character_
    )
  ) %>%
  ungroup()

# 3. Count the number of rows where the maximum numeric value belongs to each strain
max_counts <- df_data_with_max %>%
  dplyr::count(MaxStrain)

# Print the counts
print(max_counts)

max_counts_no_na <- max_counts %>% 
  filter(!is.na(MaxStrain))

# Ensure the factor levels are correctly set
max_counts_no_na$MaxStrain <- factor(max_counts_no_na$MaxStrain, levels = c("JU1421", "AF16", "N2"))

# Create a pie chart using ggplot2:
pie_chart <- ggplot(max_counts_no_na, aes(x = "", y = n, fill = MaxStrain)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 1.5) +
  coord_polar(theta = "y", start = 0) +
  labs(title = "Distribution of Maximum Strain Counts", fill = "Strain") +
  theme_void() +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "JU1421" = "#EFC000",  # Custom orange
    "AF16"   = "#0073C2",  # Custom light blue
    "N2"     = "#CD534C"   # Custom green
  ))

# Display the pie chart
print(pie_chart)

ggsave(filename = "Orthogroup_piechart.pdf", 
       plot = pie_chart, 
       device = "pdf", 
       width = 4, 
       height = 4)

df_filtered <- df_data %>%
  filter(rowSums(across(c(JU1421, AF16, N2))) >= 20)

# For each row in df_filtered, determine the dominant strain
df_data_with_dominant <- df_filtered %>%
  rowwise() %>%
  mutate(
    # Gather the counts for the three strains in the current row
    dominant_value = max(c_across(c(JU1421, AF16, N2))),
    count_dominant = sum(c_across(c(JU1421, AF16, N2)) == dominant_value),
    dominant_strain = ifelse(
      count_dominant == 1,
      c("JU1421", "AF16", "N2")[which.max(c_across(c(JU1421, AF16, N2)))],
      NA_character_
    )
  ) %>%
  ungroup()

# Count the number of rows for each dominant strain (ignoring ties)
strain_distribution <- df_data_with_dominant %>%
  dplyr::count(dominant_strain)

# Print the counts
print(strain_distribution)

# Remove rows with NA in dominant_strain
strain_distribution_no_na <- strain_distribution %>% 
  filter(!is.na(dominant_strain))

# Ensure the factor levels are correctly set
strain_distribution_no_na$dominant_strain <- factor(strain_distribution_no_na$dominant_strain, levels = c("JU1421", "AF16", "N2"))

# Create a pie chart using ggplot2:
pie_chart <- ggplot(strain_distribution_no_na, aes(x = "", y = n, fill = dominant_strain)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 1.5) +
  coord_polar(theta = "y", start = 0) +
  labs(title = "Distribution of Dominant Strain Counts", fill = "Strain") +
  theme_void() +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  # Custom color palette: change the colors below to your desired choices
  scale_fill_manual(values = c(
    "JU1421" = "#EFC000",  # Custom orange
    "AF16"   = "#0073C2",  # Custom light blue
    "N2"     = "#CD534C"   # Custom green
  ))

# Display the pie chart
print(pie_chart)

ggsave(filename = "Orthogroup_piechart_2.pdf", 
       plot = pie_chart, 
       device = "pdf", 
       width = 4, 
       height = 4)

################################################
#checking right interval for data visualization
# Filter rows where the total gene count is between 10 and 20 (inclusive)
df_interval <- df_data %>%
  filter(
    rowSums(across(c(JU1421, AF16, N2))) >= 20
  )

# For each row in the filtered data,
# determine the dominant strain if the maximum count is unique.
df_interval_with_max <- df_interval %>%
  rowwise() %>%
  mutate(
    max_value = max(c_across(c(JU1421, AF16, N2))),
    num_max = sum(c_across(c(JU1421, AF16, N2)) == max_value),
    MaxStrain = ifelse(
      num_max == 1,
      c("JU1421", "AF16", "N2")[which.max(c_across(c(JU1421, AF16, N2)))],
      NA_character_
    )
  ) %>%
  ungroup()

# Count the number of rows (orthogroups) for each dominant strain in this interval
max_counts_interval <- df_interval_with_max %>%
  dplyr::filter(!is.na(MaxStrain)) %>%  # Exclude ties
  dplyr::count(MaxStrain)

# Print the counts for the customized interval 10~20
print(max_counts_interval)
###############################################################################



# ==============================================================================
# 1. Start with your original data: df_data with columns "JU1421", "AF16", "N2"
# ==============================================================================

# (Assuming df_data is already read in and available.)

# ==============================================================================
# 2. Add a total count and assign each row an interval based on total counts.
#    New intervals:
#      - <5
#      - 5-25
#      - >25
# ==============================================================================

df_data_with_interval <- df_data %>%
  mutate(
    total = rowSums(across(c(JU1421, AF16, N2))),
    interval = case_when(
      total > 25              ~ ">25",
      total >= 5 & total <= 25  ~ "5-25",
      total < 5               ~ "<5"
    )
  )

# ==============================================================================
# 3. For each row, determine the dominant strain.
#    Only assign a dominant name if the maximum count is unique.
# ==============================================================================

df_data_with_dominant <- df_data_with_interval %>%
  rowwise() %>%
  mutate(
    dominant_value = max(c_across(c(JU1421, AF16, N2))),
    count_dominant = sum(c_across(c(JU1421, AF16, N2)) == dominant_value),
    dominant_strain = ifelse(
      count_dominant == 1,
      # Select the strain corresponding to the unique maximum.
      c("JU1421", "AF16", "N2")[which.max(c_across(c(JU1421, AF16, N2)))],
      NA_character_
    )
  ) %>%
  ungroup()

# ==============================================================================
# 4. Count rows (ignoring ties where dominant_strain is NA) for each interval
#    and for each dominant_strain.
# ==============================================================================

strain_distribution_intervals <- df_data_with_dominant %>%
  dplyr::filter(!is.na(dominant_strain)) %>%
  dplyr::count(interval, dominant_strain)

# ==============================================================================
# 5. Compute the proportion for each dominant strain within each interval.
#    This ensures that each facet (interval) represents a whole (100%).
# ==============================================================================

strain_distribution_intervals <- strain_distribution_intervals %>%
  group_by(interval) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ==============================================================================
# 6. Set factor levels so that facets and legend appear in the desired order.
# ==============================================================================

strain_distribution_intervals$dominant_strain <- factor(
  strain_distribution_intervals$dominant_strain, 
  levels = c("JU1421", "AF16", "N2")
)

strain_distribution_intervals$interval <- factor(
  strain_distribution_intervals$interval, 
  levels = c(">25", "5-25", "<5")
)

# ==============================================================================
# 7. Create the pie chart.
#    Each facet corresponds to one interval (with each facet's bars summing to 1)
#    and the labels now show the actual count (n) rather than percentages.
# ==============================================================================

pie_chart <- ggplot(strain_distribution_intervals, aes(x = "", y = prop, fill = dominant_strain)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 1.5) +
  coord_polar(theta = "y", start = 0) +
  facet_wrap(~ interval) +
  labs(
    title = "Distribution of Dominant Strain Counts\nby Total Count Intervals",
    fill = "Dominant Strain"
  ) +
  theme_void() +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(
    "JU1421" = "#EFC000",  # Custom orange
    "AF16"   = "#0073C2",  # Custom light blue
    "N2"     = "#CD534C"   # Custom red (adjust as desired)
  ))

# ==============================================================================
# 8. Display and (optionally) save the pie chart.
# ==============================================================================

print(pie_chart)

ggsave(filename = "Dominant_Strain_Distribution_Intervals.pdf", 
       plot = pie_chart, 
       device = "pdf", 
       width = 8, 
       height = 4)




df_intra_table<- read.table(
  file = "/Users/jeffrey/Desktop/genome_file/orthogroup/Final/Intra_all_Orthogroups.tsv", 
  header = TRUE, 
  sep = "\t", 
  stringsAsFactors = FALSE
) 

df_N2_genes <- df_intra_table %>%
  dplyr::select(Orthogroup, N2) %>%
  separate_rows(N2, sep = ",\\s*") %>%
  dplyr::rename(gene_name = N2) %>%
  filter(!is.na(gene_name) & gene_name != "") 

JU1421_orthogroups <- df_data_with_max %>%
  filter(MaxStrain == "JU1421") %>%
  pull(Orthogroup)

AF16_orthogroups <- df_data_with_max %>%
  filter(MaxStrain == "AF16") %>%
  pull(Orthogroup)

N2_orthogroups <- df_data_with_max %>%
  filter(MaxStrain == "N2") %>%
  pull(Orthogroup)

df_pfam <- read.csv(
  file             = "/Users/jeffrey/Desktop/genome_file/orthogroup/Final/Pfam.txt",
  header           = FALSE,        # file itself has no header row
  sep              = "\t",
  stringsAsFactors = FALSE,
  col.names        = c("gene_name", "pfam", "Desc",
                       "start", "end", "Interpro", "strain")
)

table(df_pfam$strain)
head(df_pfam)

df_gene <- read.table(
  file             = "/Users/jeffrey/Desktop/genome_file/orthogroup/Final/Orthogroups.gene.txt",
  header           = FALSE,        # file itself has no header row
  sep              = "\t",
  stringsAsFactors = FALSE,
  col.names        = c("gene_name", "orthogroup")
)

df_pfam <- df_pfam %>%
  left_join(df_gene, by = "gene_name")

JU1421_gene_names <- df_pfam %>%
  filter(strain == "JU1421", orthogroup %in% JU1421_orthogroups) %>%
  pull(gene_name)

AF16_gene_names <- df_pfam %>%
  filter(strain == "AF16", orthogroup %in% AF16_orthogroups) %>%
  pull(gene_name)

N2_gene_names <- df_N2_genes%>%
  filter(Orthogroup %in% N2_orthogroups) %>%
  pull(gene_name)

# Extract GO IDs
go_ids <- keys(GO.db)

# Create a data frame to store results
go_data <- tibble(
  GOID = go_ids,
  TERM = Term(GOTERM[go_ids]),
  ONTOLOGY = Ontology(GOTERM[go_ids]),
  DEFINITION = Definition(GOTERM[go_ids]),
  stringsAsFactors = FALSE
)
go_data

read.table('/Users/jeffrey/Desktop/genome_file/orthogroup/Final/GO/JU1421.GO.txt', col.names = c('Gene','GO')) %>%
  tibble() %>%
  dplyr::select(GO,Gene) -> JU1421.T2G

JU1421.T2G

read.table('/Users/jeffrey/Desktop/genome_file/orthogroup/Final/GO/AF16.GO.txt', col.names = c('Gene','GO')) %>%
  tibble() %>%
  dplyr::select(GO,Gene) -> AF16.T2G

AF16.T2G


JU1421.T2G %>%
  left_join(go_data,by=join_by(GO == GOID)) %>%
  filter(ONTOLOGY == 'BP') %>%
  dplyr::select(GO,Gene) -> JU1421.T2G.BP
JU1421.T2G.BP

JU1421.T2G %>%
  left_join(go_data,by=join_by(GO == GOID)) %>%
  filter(ONTOLOGY == 'MF') %>%
  dplyr::select(GO,Gene) -> JU1421.T2G.MF
JU1421.T2G.MF

JU1421.T2G %>%
  left_join(go_data,by=join_by(GO == GOID)) %>%
  filter(ONTOLOGY == 'CC') %>%
  dplyr::select(GO,Gene) -> JU1421.T2G.CC
JU1421.T2G.CC

AF16.T2G %>%
  left_join(go_data,by=join_by(GO == GOID)) %>%
  filter(ONTOLOGY == 'BP') %>%
  dplyr::select(GO,Gene) -> AF16.T2G.BP
AF16.T2G.BP

AF16.T2G %>%
  left_join(go_data,by=join_by(GO == GOID)) %>%
  filter(ONTOLOGY == 'MF') %>%
  dplyr::select(GO,Gene) -> AF16.T2G.MF
AF16.T2G.MF

AF16.T2G %>%
  left_join(go_data,by=join_by(GO == GOID)) %>%
  filter(ONTOLOGY == 'CC') %>%
  dplyr::select(GO,Gene) -> AF16.T2G.CC
AF16.T2G.CC

go_data %>% 
  dplyr::select(GOID,TERM) -> T2N

JU1421_GO <- enricher(JU1421_gene_names, TERM2GENE = JU1421.T2G.BP , TERM2NAME = T2N,
                      pvalueCutoff  = 0.05,         # set the p-value cutoff
                      pAdjustMethod = "BH",         # set the p-value adjustment method (e.g., "BH", "bonferroni")
                      qvalueCutoff  = 0.01)
remove_ids <- c("GO:0015074","GO:0001707","GO:0003007","GO:0000723","GO:0006281","GO:0001708")
JU1421_GO@result <- JU1421_GO@result[ !JU1421_GO@result$ID %in% remove_ids, ]
test <- as.data.frame(JU1421_GO@result)

AF16_GO <- enricher(AF16_gene_names, TERM2GENE = AF16.T2G.BP , TERM2NAME = T2N,
                      pvalueCutoff  = 0.05,         # set the p-value cutoff
                      pAdjustMethod = "BH",         # set the p-value adjustment method (e.g., "BH", "bonferroni")
                      qvalueCutoff  = 0.01)
remove_ids <- c("GO:0015074","GO:0006310","GO:0007219")
AF16_GO@result <- AF16_GO@result[ !AF16_GO@result$ID %in% remove_ids, ]


test <- as.data.frame(AF16_GO@result)


N2_GO <- enrichGO(
  gene          = N2_gene_names,
  OrgDb         = org.Ce.eg.db,
  keyType       = "WORMBASE",  # change to "ENTREZID" if your IDs are Entrez
  ont           = "BP",        # specify 'BP' for Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE         # converts gene IDs to gene symbols if possible
)
remove_ids <- c("GO:0007600","GO:0050906","GO:0050911","GO:0050907","GO:0007608", "GO:0051606","GO:0043207", 
"GO:0006952","GO:0006955","GO:0002376")
N2_GO@result <- N2_GO@result[ !N2_GO@result$ID %in% remove_ids, ]
test <- as.data.frame(N2_GO@result)

summary(N2_GO)
N2_GO_details <- as.data.frame(N2_GO)
JU1421_GO_details <- as.data.frame(JU1421_GO)
AF16_GO_details <- as.data.frame(AF16_GO)

common_ids <- intersect(JU1421_GO_details$ID, N2_GO_details $ID)

# View the first few rows of the detailed results

head(go_details)
# Optionally, check its structure
str(go_details)

# Plot for N2_GO
p_N2 <- clusterProfiler::dotplot(N2_GO, 
                                 x = "Count",                 
                                 showCategory = 10,           
                                 font.size = 12,              
                                 title = "Enriched GO Terms (N2)",
                                 color = "p.adjust") +
  scale_color_gradient(low = "#3466A6", high = "#D77760") +
  # Customize the size legend for GeneRatio (bubbles) to only show three reference values.
  scale_size_continuous(
    # Adjust these break values according to your GeneRatio range.
    breaks = c(0.06, 0.08, 0.1),
    labels = c("0.06", "0.08", "0.1")
  ) +
  theme_minimal() +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor  = element_blank()
  )


# Display the plot
print(p_N2)
# Plot for AF16_GO
p_AF16 <- clusterProfiler::dotplot(AF16_GO, 
                                   x = "Count",                
                                   showCategory = 10,          
                                   font.size = 12,             
                                   title = "Enriched GO Terms (AF16)",
                                   color = "p.adjust") +
  scale_color_gradient(low = "#3466A6", high = "#D77760") +
  theme_minimal() +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor  = element_blank()
  )
# Display the plot
print(p_AF16)


# Plot for JU1421_GO
p_JU1421 <- clusterProfiler::dotplot(JU1421_GO, 
                                     x = "Count",                 
                                     showCategory = 10,  # note: different number than others
                                     font.size = 12,              
                                     title = "Enriched GO Terms (JU1421)",
                                     color = "p.adjust") +
  scale_color_gradient(low = "#3466A6", high = "#D77760") +
  theme_minimal() +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor   = element_blank()
  )
# Display the plot
print(p_JU1421)

combined_plot <- p_N2 + p_AF16 + p_JU1421 + plot_layout(ncol = 3)
print(combined_plot)

ggsave(filename = "Orthogroup_GO.pdf", 
       plot = combined_plot, 
       device = "pdf", 
       width = 17, 
       height = 3.5)



