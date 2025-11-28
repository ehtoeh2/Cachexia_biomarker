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
# ----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(sva)

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# 2. Target Genes for Visualization
# ----------------------------------------------------------------------------
my_target_genes <- c("SERPINA3", "LBP", "APOD", "C4B", "CFH", "LUM", 
                     "CPQ", "CTSL", "IFI30", "THBS4")

# ============================================================================
# SECTION A: Volcano Plot - Human Serum Proteomics
# ============================================================================

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
volcano_plot <- ggplot(human_volcano, 
                       aes(x = diff, y = -log10(adj_pval), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), 
             linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(PVAL_CUTOFF), 
             linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Up" = "#F5A946", 
                                "Down" = "#1A3664", 
                                "NS" = "gray80")) +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20) +
  labs(
    title = "Human Serum Proteomics (Cachexia vs Control)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_bw()

print(volcano_plot)
ggsave(file.path(OUTPUT_DIR, "Human_volcano_plot.pdf"), 
       volcano_plot, width = 8, height = 6)

# ============================================================================
# SECTION B: Cross-Species Validation Scatter Plot
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

# Calculate Pearson correlation
cor_val <- cor(plot_data$Mouse_LogFC, plot_data$Human_LogFC, method = "pearson")

cat("Cross-species correlation (Pearson r):", round(cor_val, 3), "\n")

# Create scatter plot
scatter_plot <- ggplot(plot_data, aes(x = Mouse_LogFC, y = Human_LogFC)) +
  # Quadrant background
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, 
           fill = "#E41A1C", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           fill = "#377EB8", alpha = 0.05) +
  # Points
  geom_point(aes(color = Concordance, shape = Origin), size = 4, alpha = 0.8) +
  # Labels for concordant genes
  geom_text_repel(
    data = subset(plot_data, Concordance != "Inconsistent"),
    aes(label = Human_Gene),
    size = 4.5, fontface = "bold", box.padding = 0.5
  ) +
  # Color scale
  scale_color_manual(values = c("Consistent Up" = "#E41A1C",
                                "Consistent Down" = "#377EB8",
                                "Inconsistent" = "gray70")) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # Labels
  labs(
    title = "Cross-Species Validation (LogFC Correlation)",
    subtitle = paste0("Mouse Serum Protein vs Human Serum Protein (R = ", 
                      round(cor_val, 3), ")"),
    x = "Mouse Serum Log2 Fold Change",
    y = "Human Serum Log2 Fold Change"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(scatter_plot)
ggsave(file.path(OUTPUT_DIR, "Cross_species_scatter.pdf"), 
       scatter_plot, width = 8, height = 7)

# ============================================================================
# SECTION C: PCA Analysis with Batch Effect Correction
# ============================================================================

# ----------------------------------------------------------------------------
# C1. Data Preparation
# ----------------------------------------------------------------------------
raw_data <- human_protein_per_sample
rownames(raw_data) <- raw_data$Protein
raw_data <- raw_data %>% dplyr::select(-Protein)

# Remove missing values (required for ComBat)
clean_data <- na.omit(raw_data)

cat("Proteins after NA removal:", nrow(clean_data), "\n")

# ----------------------------------------------------------------------------
# C2. Create Metadata
# ----------------------------------------------------------------------------
samples <- colnames(clean_data)
condition <- substr(samples, 1, 2)  # CT = Control, CX = Cachexia

# Define batch information (modify according to your experimental design)
batch <- c(rep(1, 3), rep(1, 3), rep(2, 3), rep(2, 3))

meta_data <- data.frame(
  Sample = samples,
  condition = factor(condition, levels = c("CT", "CX")),
  batch = factor(batch)
)

print(meta_data)

# ----------------------------------------------------------------------------
# C3. Batch Effect Correction using ComBat
# ----------------------------------------------------------------------------
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

cat("ComBat correction completed.\n")

# ----------------------------------------------------------------------------
# C4. PCA Calculation
# ----------------------------------------------------------------------------
# Note: scale. = FALSE for log2-transformed data (variance already stabilized)
pca_res <- prcomp(t(combat_mat), scale. = FALSE)

# Calculate variance explained
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

cat("Variance explained - PC1:", percentVar[1], "%, PC2:", percentVar[2], "%\n")

# Prepare PCA data for plotting
pcaData <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  condition = meta_data$condition,
  batch = meta_data$batch,
  name = meta_data$Sample
)

# Rename conditions for visualization
pcaData$condition <- ifelse(pcaData$condition == "CT", "Control", "Cachexia")
pcaData$condition <- factor(pcaData$condition, levels = c("Control", "Cachexia"))

# ----------------------------------------------------------------------------
# C5. PCA Plot
# ----------------------------------------------------------------------------
pca_colors <- c("Control" = "#F8766D", "Cachexia" = "#00BFC4")
pca_fill <- c("Control" = "lightgray", "Cachexia" = "skyblue")

pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, fill = condition)) +
  # Confidence ellipse
  stat_ellipse(
    geom = "polygon",
    level = 0.95,
    alpha = 0.2,
    show.legend = FALSE
  ) +
  # Points
  geom_point(size = 3, alpha = 0.9) +
  # Colors
  scale_color_manual(values = pca_colors) +
  scale_fill_manual(values = pca_fill) +
  # Axis labels with variance explained
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  # Labels and theme
  labs(
    title = "PCA Plot (Post-ComBat)",
    subtitle = "Human Serum Proteomics Data",
    color = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-7, 7), ylim = c(-8, 8))

print(pca_plot)
ggsave(file.path(OUTPUT_DIR, "Human_PCA_plot.pdf"), 
       pca_plot, width = 7, height = 6)

# ============================================================================
# Summary Statistics
# ============================================================================
cat("\n============================================\n")
cat("Human Serum Validation Analysis Complete\n")
cat("============================================\n")

# Volcano plot summary
cat("\nVolcano Plot Summary:\n")
cat("  Total proteins:", nrow(human_protein), "\n")
cat("  Upregulated (adj.P < 0.05, log2FC > 0.58):", 
    sum(human_volcano$Significance == "Up"), "\n")
cat("  Downregulated (adj.P < 0.05, log2FC < -0.58):", 
    sum(human_volcano$Significance == "Down"), "\n")

# Cross-species summary
cat("\nCross-Species Validation:\n")
cat("  Total matched proteins:", nrow(plot_data), "\n")
cat("  Concordant (same direction):", 
    sum(plot_data$Concordance != "Inconsistent"), "\n")
cat("  Pearson correlation:", round(cor_val, 3), "\n")

# PCA summary
cat("\nPCA Analysis:\n")
cat("  Proteins used:", nrow(combat_mat), "\n")
cat("  PC1 variance:", percentVar[1], "%\n")
cat("  PC2 variance:", percentVar[2], "%\n")

cat("\n============================================\n")

# Session info
sessionInfo()
