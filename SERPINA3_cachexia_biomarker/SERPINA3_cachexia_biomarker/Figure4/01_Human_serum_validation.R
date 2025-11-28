# ============================================================================
# Human Serum Proteomics Validation (PXD035832)
# ============================================================================
# This script performs:
#   - Volcano plot analysis of human serum proteomics
#   - Cross-species validation scatter plot (Mouse vs Human LogFC)
#   - PCA analysis with batch effect correction (ComBat)
#
# Data source: ProteomeXchange PXD035832
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(sva)

set.seed(1234)

# 1. Data Paths Configuration
# NOTE: Modify these paths according to your local environment
DATA_DIR <- "path/to/your/data"
OUTPUT_DIR <- file.path(DATA_DIR, "results")

# Required input files:
#   - proDa_result.xlsx: Differential expression analysis results
#   - human_prot_per_sample.xlsx: Sample-level protein abundance matrix
#   - fap_ready, myo_ready: Mouse proteomics data (from previous analysis)

# Load data
human_protein <- read.xlsx(file.path(DATA_DIR, "proDa_result.xlsx"))
human_protein_per_sample <- read.xlsx(file.path(DATA_DIR, "human_prot_per_sample.xlsx"))

# Load mouse data (from Figure 2/3 analysis)
# fap_ready and myo_ready should contain: Gene, Prot_logFC columns
# Mouse_protein for reference

# ============================================================================
# Figure 4A: PCA Analysis with Batch Effect Correction
# ============================================================================

# Data Preparation

raw_data <- human_protein_per_sample
rownames(raw_data) <- raw_data$Protein
raw_data <- raw_data %>% dplyr::select(-Protein)

# Remove missing values (required for ComBat)
clean_data <- na.omit(raw_data)


# Create Metadata
samples <- colnames(clean_data)
condition <- substr(samples, 1, 2)  # CT = Control, CX = Cachexia

# Define batch information (modify according to your experimental design)
batch <- c(rep(1, 3), rep(1, 3), rep(2, 3), rep(2, 3))

meta_data <- data.frame(
  Sample = samples,
  condition = factor(condition, levels = c("CT", "CX")),
  batch = factor(batch)
)


# Batch Effect Correction using ComBat
# Model matrix preserving condition effect
mod_combat <- model.matrix(~ condition, data = meta_data)

# Run ComBat
combat_mat <- ComBat(
  dat = as.matrix(clean_data),
  batch = meta_data$batch,
  mod = mod_combat,
  par.prior = TRUE,
  prior.plots = FALSE
)

# PCA Calculation
# Note: scale. = FALSE for log2-transformed data (variance already stabilized)
pca_res <- prcomp(t(combat_mat), scale. = FALSE)
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# Prepare PCA data for plotting
pcaData <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  condition = meta_data$condition,
  batch = meta_data$batch,
  name = meta_data$Sample
)


# ----------------------------------------------------------------------------
# PCA Plot
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point()



# ============================================================================
# Figure 4B: Volcano Plot - Human Serum Proteomics
# ============================================================================

# Target genes identified in Figures 2 and 3
my_target_genes <- c("SERPINA3", "LBP", "APOD", "C4B", "CFH", "LUM", 
                     "CPQ", "CTSL", "IFI30", "THBS4")

# Significance thresholds
PVAL_CUTOFF <- 0.05
FC_CUTOFF <- 0.58  # log2(1.5)

# Prepare volcano plot data
human_volcano <- human_protein %>%
  mutate(
    Significance = case_when(
      adj_pval < PVAL_CUTOFF & diff > FC_CUTOFF  ~ "Up",
      adj_pval < PVAL_CUTOFF & diff < -FC_CUTOFF ~ "Down",
      TRUE ~ "NS"
    ),
    # Label: significant genes with |diff| > 1
    Label = ifelse(adj_pval < PVAL_CUTOFF & abs(diff) > 1, Gene, NA)
  )

# Create volcano plot
volcano_plot <- ggplot(human_volcano, aes(x = diff, y = -log10(adj_pval))) +
  geom_point()



# ============================================================================
# Figure 4C: Cross-Species Validation Scatter Plot
# ============================================================================

# Combine mouse target data
mouse_targets <- bind_rows(fap_ready, myo_ready)

# Map mouse gene names to human orthologs
mouse_targets$Human_Gene <- toupper(mouse_targets$Gene)

# Handle specific mouse-human ortholog mappings
mouse_targets$Human_Gene[mouse_targets$Gene == "Serpina3n"] <- "SERPINA3"
mouse_targets$Human_Gene[mouse_targets$Gene == "Serpina3m"] <- "SERPINA3"

# Merge mouse and human data
cross_species_df <- mouse_targets %>%
  inner_join(human_protein, by = c("Human_Gene" = "Gene")) %>%
  dplyr::select(
    Mouse_Gene = Gene,
    Human_Gene,
    Origin,
    Mouse_LogFC = Prot_logFC,
    Human_LogFC = diff,
    Human_AdjP = adj_pval
  )

# Determine concordance between species
plot_data <- cross_species_df %>%
  mutate(
    Concordance = case_when(
      Mouse_LogFC > 0 & Human_LogFC > 0 ~ "Consistent Up",
      Mouse_LogFC < 0 & Human_LogFC < 0 ~ "Consistent Down",
      TRUE ~ "Inconsistent"
    )
  )

# Create scatter plot
scatter_plot <- ggplot(plot_data, aes(x = Mouse_LogFC, y = Human_LogFC)) 
  
# Points
  geom_point(aes(color = Concordance, shape = Origin), size = 4, alpha = 0.8) +
  # Labels for concordant genes
  geom_text_repel(
    data = subset(plot_data, Concordance != "Inconsistent"),
    aes(label = Human_Gene)) 



# ============================================================================
# Summary Statistics
# ============================================================================
cat("\n============================================\n")
cat("Human Serum Validation Analysis Complete\n")
cat("============================================\n")
