# ============================================================================
# Figure 3: Mouse Serum Proteomics Validation
# ============================================================================
# This script performs:
#   - Figure 3F-G: Serum proteomics heatmap (by GC group and cell type origin)
#   - Figure 3H: Multi-omics integration (Stouffer's method)
#   - Figure 3I: Concordant biomarker lollipop plots
#
# Author: Dae-Hwan Kim, Jibeom Ko
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
# 2. Figure 3F: Proteomics Heatmap by GC Group
# ----------------------------------------------------------------------------
cat("Creating proteomics heatmap by GC group...\n")

# Create gene grouping list
gene_groups_list <- list(
  "GC2_Inflammatory" = unique(c(FAP_gc2_survivors, Myonuclei_GC2_survivors)),
  "GC3_Catabolic"    = unique(c(FAP_gc3_survivors, Myonuclei_GC3_survivors)),
  "GC4_Structural"   = unique(c(FAP_gc4_survivors, Myonuclei_GC4_survivors))
)

# Convert to data frame for annotation
gene_meta <- stack(gene_groups_list)
colnames(gene_meta) <- c("Gene", "GC_Group")

# Filter proteomics data for target genes
target_prot <- ms_prot_sample %>%
  filter(Gene %in% gene_meta$Gene) %>%
  distinct(Gene, .keep_all = TRUE)

cat("Target survivors:", nrow(gene_meta), "\n")
cat("Serum detected:", nrow(target_prot), "\n")

# Prepare matrix
rownames(target_prot) <- target_prot$Gene
mat_prot <- target_prot %>% 
  dplyr::select(-Gene) %>% 
  as.matrix()

# Z-score transformation (row-wise)
scaled_prot <- t(scale(t(mat_prot)))

# Define sample order (Control -> Cachexia)
sample_order <- c(
  "F_CT1", "F_CT2", "M_CT1", "M_CT2", "M_CT3",  # Control
  "F_CX4", "F_CX5", "M_CX6", "M_CX8", "M_CX9"   # Cachexia
)
valid_samples <- intersect(sample_order, colnames(scaled_prot))
scaled_prot <- scaled_prot[, valid_samples]

# Row annotation (GC group)
row_split_vec <- gene_meta$GC_Group[match(rownames(scaled_prot), gene_meta$Gene)]

# Column annotation (Condition)
sample_condition <- ifelse(grepl("CX", colnames(scaled_prot)), "Cachexia", "Control")
col_anno <- HeatmapAnnotation(
  Condition = sample_condition,
  col = list(Condition = c("Control" = "#BEBEBE", "Cachexia" = "#1A3664"))
)

# Draw heatmap
pdf(file.path(OUTPUT_DIR, "Figure3F_Proteomics_GC_heatmap.pdf"), width = 8, height = 10)
Heatmap(
  scaled_prot,
  name = "Z-score",
  col = colorRamp2(c(-1.5, 0, 1.5), c("#113F77", "white", "#751B3E")),
  
  # Row settings (split by GC group)
  row_split = factor(row_split_vec, 
                     levels = c("GC2_Inflammatory", "GC3_Catabolic", "GC4_Structural")),
  cluster_row_slices = FALSE,
  row_title_rot = 0,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  
  # Column settings
  top_annotation = col_anno,
  column_order = valid_samples,
  column_split = factor(sample_condition, levels = c("Control", "Cachexia")),
  cluster_columns = FALSE,
  column_names_rot = 45,
  
  # Design
  rect_gp = gpar(col = "white", lwd = 1),
  border = TRUE,
  column_title = "Serum Proteomics Validation (by GC Group)"
)
dev.off()

# ----------------------------------------------------------------------------
# 3. Figure 3G: Proteomics Heatmap by Cell Type Origin
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

cat("Gene origins:\n")
print(table(gene_origins))

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
  Condition = sample_cond2,
  col = list(Condition = c("Control" = "#BEBEBE", "Cachexia" = "#1A3664"))
)

# Cell size settings
cell_h <- unit(3, "mm")
cell_w <- unit(5, "mm")

# Draw heatmap
pdf(file.path(OUTPUT_DIR, "Figure3G_Proteomics_CellType_heatmap.pdf"), width = 8, height = 10)
Heatmap(
  scaled_prot2,
  name = "Z-score",
  col = colorRamp2(c(-1.5, 0, 1.5), c("#113F77", "white", "#751B3E")),
  
  # Fixed cell size
  width = ncol(scaled_prot2) * cell_w,
  height = nrow(scaled_prot2) * cell_h,
  
  # Row settings (split by origin)
  row_split = factor(row_split_vec2, levels = c("FAP Origin", "Myonuclei Origin")),
  cluster_row_slices = FALSE,
  row_title_rot = 0,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  
  # Column settings
  top_annotation = col_anno2,
  column_order = valid_samples2,
  column_split = factor(sample_cond2, levels = c("Control", "Cachexia")),
  cluster_columns = FALSE,
  column_names_rot = 45,
  
  # Design
  rect_gp = gpar(col = "white", lwd = 1),
  border = TRUE,
  column_title = "Serum Proteomics Validation (by Cell Type Origin)"
)
dev.off()

# ----------------------------------------------------------------------------
# 4. Figure 3H: Multi-omics Integration (RNA-seq + Proteomics)
# ----------------------------------------------------------------------------
cat("\nPerforming multi-omics integration...\n")

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

cat("Concordant genes:", nrow(concordant_df), "\n")

# ----------------------------------------------------------------------------
# 5. Stouffer's Method for Meta-analysis
# ----------------------------------------------------------------------------
# Function to calculate combined Z-score
calculate_stouffer <- function(pvals, logfcs) {
  # Safety check for very small p-values
  pvals[pvals < 1e-300] <- 1e-300
  
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

cat("\nTop 10 concordant biomarkers:\n")
print(head(concordant_df, 10))

# Save concordance results
write.xlsx(concordant_df, file.path(OUTPUT_DIR, "Concordant_biomarkers.xlsx"))

# ----------------------------------------------------------------------------
# 6. Figure 3I: Lollipop Plots for Top Biomarkers
# ----------------------------------------------------------------------------
cat("\nCreating lollipop plots...\n")

# Add direction column
lollipop_data <- concordant_df %>%
  mutate(Direction = ifelse(Stouffer_Z > 0, "Up-regulated", "Down-regulated")) %>%
  arrange(desc(abs(Stouffer_Z)))

# Function to create lollipop plot
create_custom_lollipop <- function(data, title_text, limit_range = c(-10, 10)) {
  
  ggplot(data, aes(x = Stouffer_Z, y = Gene, color = Direction)) +
    
    # Lollipop stick
    geom_segment(aes(x = 0, xend = Stouffer_Z, y = Gene, yend = Gene), 
                 size = 1.2, alpha = 0.8) +
    
    # Lollipop head
    geom_point(size = 3) +
    
    # Reference line
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    
    # Colors
    scale_color_manual(values = c("Up-regulated" = "#751B3E", 
                                  "Down-regulated" = "#113F77")) +
    
    # X-axis range
    coord_cartesian(xlim = limit_range) +
    
    # Labels
    labs(
      title = title_text,
      x = "Combined Z-score (Directional)",
      y = NULL
    ) +
    
    # Theme
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(face = "bold.italic", color = "black"),
      axis.text.x = element_text(face = "bold", color = "black"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

# Prepare data by cell type origin
fap_ready <- lollipop_data %>%
  filter(Origin == "FAP Origin") %>%
  head(10) %>%
  mutate(Gene = factor(Gene, levels = Gene[order(Stouffer_Z)]))

myo_ready <- lollipop_data %>%
  filter(Origin == "Myonuclei Origin") %>%
  head(10) %>%
  mutate(Gene = factor(Gene, levels = Gene[order(Stouffer_Z)]))

cat("FAP top genes:", nrow(fap_ready), "\n")
cat("Myonuclei top genes:", nrow(myo_ready), "\n")

# X-axis range
my_xlim <- c(-10, 10)

# Create individual plots
plot_fap <- create_custom_lollipop(fap_ready, "FAP Origin Biomarkers", my_xlim)
plot_myo <- create_custom_lollipop(myo_ready, "Myonuclei Origin Biomarkers", my_xlim)

# Combine plots
final_fig <- (plot_fap | plot_myo) +
  plot_annotation(
    title = 'Top Concordant Biomarkers in Serum & Tissue',
    subtitle = '(Integration of RNA-seq & Proteomics)',
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

print(final_fig)

# Save
ggsave(file.path(OUTPUT_DIR, "Figure3I_Lollipop_concordant.pdf"), 
       final_fig, width = 14, height = 8)

# Individual plots
ggsave(file.path(OUTPUT_DIR, "Figure3I_Lollipop_FAP.pdf"), 
       plot_fap, width = 7, height = 6)
ggsave(file.path(OUTPUT_DIR, "Figure3I_Lollipop_Myo.pdf"), 
       plot_myo, width = 7, height = 6)

# ----------------------------------------------------------------------------
# 7. Save Session Info
# ----------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "Figure3_session_info.txt"))

cat("\n============================================\n")
cat("Figure 3 analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
cat("\nKey outputs:\n")
cat("  - Figure3F_Proteomics_GC_heatmap.pdf\n")
cat("  - Figure3G_Proteomics_CellType_heatmap.pdf\n")
cat("  - Figure3I_Lollipop_concordant.pdf\n")
cat("  - Concordant_biomarkers.xlsx\n")
cat("============================================\n")
