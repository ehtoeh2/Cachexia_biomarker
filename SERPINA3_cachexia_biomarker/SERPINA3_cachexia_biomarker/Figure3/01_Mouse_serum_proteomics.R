# ============================================================================
# Figure 3: Mouse Serum Proteomics Validation
# ============================================================================
# This script performs:
#   - Figure 3G: Volcano plot
#   - Figure 3H: Mouse serum proteomics heatmap
#   - Figure 3I: Concordant biomarker lollipop plots
#
# Author: Dae-Hwan Kim, Ji beom Ko
# Paper: "An integrative multi-omics framework identifies SERPINA3 as a 
#         circulating biomarker for cancer cachexia"
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(patchwork)

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment
#
# Required input files:
#   - Mouse serum proteomics data (Mass spectrometry results)
#   - Per-sample proteomics abundance data
#   - Previously generated: res_df (DESeq2 results), survivor lists

DATA_DIR <- "path/to/your/data"
OUTPUT_DIR <- file.path(DATA_DIR, "results")

# Load data files
# Mouse_protein: Mouse serum proteomics data with proDA analysis results
# ms_prot_sample: Per-sample protein abundance data
Mouse_protein <- read.xlsx(file.path(DATA_DIR, "Mouse_serum_proteomics.xlsx"))
ms_prot_sample <- read.xlsx(file.path(DATA_DIR, "proteomics_per_sample.xlsx"))

# Load previously generated data
res_df <- readRDS(file.path(OUTPUT_DIR, "DESeq2_results.rds"))
survivor_list <- readRDS(file.path(OUTPUT_DIR, "CellType_survivors.rds"))

# Extract survivor gene lists
FAP_gc2_survivors <- survivor_list$FAP_GC2
FAP_gc3_survivors <- survivor_list$FAP_GC3
FAP_gc4_survivors <- survivor_list$FAP_GC4
Myonuclei_GC2_survivors <- survivor_list$Myo_GC2
Myonuclei_GC3_survivors <- survivor_list$Myo_GC3
Myonuclei_GC4_survivors <- survivor_list$Myo_GC4

# ----------------------------------------------------------------------------
# 2. Figure 3G: Proteomics Volcano plot
# ----------------------------------------------------------------------------

# See Figure 4B volcano plot

# ----------------------------------------------------------------------------
# 3. Figure 3H: Proteomics Heatmap by Cell Type Origin
# ----------------------------------------------------------------------------
cat("\nCreating proteomics heatmap by cell type origin...\n")

# Combine all survivor genes by cell type
all_fap_genes <- unique(c(FAP_gc2_survivors, FAP_gc3_survivors, FAP_gc4_survivors))
all_myo_genes <- unique(c(Myonuclei_GC2_survivors, Myonuclei_GC3_survivors, Myonuclei_GC4_survivors))

# Filter proteomics data
target_prot2 <- ms_prot_sample %>%
  filter(Gene %in% c(all_fap_genes, all_myo_genes)) %>%
  distinct(Gene, .keep_all = TRUE)

# Assign cell type origin (FAP takes priority if shared)
gene_origins <- character(nrow(target_prot2))
names(gene_origins) <- target_prot2$Gene

for (g in target_prot2$Gene) {
  if (g %in% all_fap_genes) {
    gene_origins[g] <- "FAP Origin"
  } else {
    gene_origins[g] <- "Myonuclei Origin"
  }
}

# Prepare matrix
rownames(target_prot2) <- target_prot2$Gene
mat_prot2 <- target_prot2 %>% 
  dplyr::select(-Gene) %>% 
  as.matrix()

# Z-score transformation
scaled_prot2 <- t(scale(t(mat_prot2)))

# Sample order
valid_samples2 <- intersect(sample_order, colnames(scaled_prot2))
scaled_prot2 <- scaled_prot2[, valid_samples2]

# Update row split vector
final_genes <- rownames(scaled_prot2)
row_split_vec2 <- gene_origins[final_genes]

# Column annotation
sample_cond2 <- ifelse(grepl("CX", colnames(scaled_prot2)), "Cachexia", "Control")
col_anno2 <- HeatmapAnnotation(
  Condition = sample_cond2)

# Draw heatmap
Heatmap(
  scaled_prot2,
  name = "Z-score"),
  
  # Row settings (split by origin)
  row_split = factor(row_split_vec2, levels = c("FAP Origin", "Myonuclei Origin")),
  cluster_row_slices = FALSE,
  row_title_rot = 0,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  
  # Column settings
  top_annotation = col_anno2,
  column_order = valid_samples2,
  column_split = factor(sample_cond2, levels = c("Control", "Cachexia")),
  cluster_columns = FALSE)


# ----------------------------------------------------------------------------
# 4. Figure 3I: Multi-omics Integration (RNA-seq + Proteomics)
# ----------------------------------------------------------------------------
# Prepare origin data frame
origin_df <- data.frame(
  Gene = names(gene_origins), 
  Origin = gene_origins, 
  stringsAsFactors = FALSE
)

# Clean RNA-seq data
res_clean <- res_df %>%
  rownames_to_column("Gene") %>%
  dplyr::select(Gene, RNA_logFC = log2FoldChange, RNA_padj = padj, RNA_pval = pvalue)

# Clean proteomics data
prot_clean <- Mouse_protein %>%
  dplyr::select(Gene, Prot_logFC = proDA_logFC, Prot_padj = proDA_padj, Prot_pval = proDA_p)

# Merge all three datasets
merged_df <- list(res_clean, prot_clean, origin_df) %>%
  purrr::reduce(inner_join, by = "Gene")

# Filter for concordant direction and remove NAs
concordant_df <- merged_df %>%
  filter(!is.na(RNA_logFC) & !is.na(Prot_logFC)) %>%
  filter(!is.na(RNA_pval) & !is.na(Prot_pval)) %>%
  filter(sign(RNA_logFC) == sign(Prot_logFC))  # Same direction


# ----------------------------------------------------------------------------
# 5. Stouffer's Method for Meta-analysis
# ----------------------------------------------------------------------------
# Function to calculate combined Z-score
calculate_stouffer <- function(pvals, logfcs) {
  # Convert p-values to Z-scores (with direction from logFC)
  z_scores <- qnorm(1 - pvals / 2) * sign(logfcs)
  
  # Combine Z-scores (Stouffer's formula)
  combined_z <- sum(z_scores) / sqrt(length(z_scores))
  
  return(combined_z)
}

# Apply to each gene
concordant_df$Stouffer_Z <- mapply(
  function(p1, fc1, p2, fc2) {
    calculate_stouffer(c(p1, p2), c(fc1, fc2))
  },
  concordant_df$RNA_pval, 
  concordant_df$RNA_logFC, 
  concordant_df$Prot_pval, 
  concordant_df$Prot_logFC
)

# Sort by absolute Z-score
concordant_df <- concordant_df %>%
  arrange(desc(abs(Stouffer_Z)))


# ----------------------------------------------------------------------------
# 6. Lollipop Plots for Top Biomarkers
# ----------------------------------------------------------------------------

# Add direction column
lollipop_data <- concordant_df %>%
  mutate(Direction = ifelse(Stouffer_Z > 0, "Up-regulated", "Down-regulated")) %>%
  arrange(desc(abs(Stouffer_Z)))

# Function to create lollipop plot
create_custom_lollipop <- function(data) {
  ggplot(data, aes(x = Stouffer_Z, y = Gene)) +
    geom_segment(aes(x = 0, xend = Stouffer_Z, y = Gene, yend = Gene)) +
    geom_point()}

# Prepare data by cell type origin
fap_ready <- lollipop_data %>%
  filter(Origin == "FAP Origin") %>%
  head(10) %>%
  mutate(Gene = factor(Gene, levels = Gene[order(Stouffer_Z)]))

myo_ready <- lollipop_data %>%
  filter(Origin == "Myonuclei Origin") %>%
  head(10) %>%
  mutate(Gene = factor(Gene, levels = Gene[order(Stouffer_Z)]))

# Create individual plots
plot_fap <- create_custom_lollipop(fap_ready, "FAP Origin Biomarkers")
plot_myo <- create_custom_lollipop(myo_ready, "Myonuclei Origin Biomarkers")


# ----------------------------------------------------------------------------
# 7. Save Session Info
# ----------------------------------------------------------------------------

cat("\n============================================\n")
cat("Figure 3 analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
