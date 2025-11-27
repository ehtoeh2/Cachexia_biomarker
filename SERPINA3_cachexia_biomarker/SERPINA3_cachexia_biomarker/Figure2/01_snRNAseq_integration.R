# ============================================================================
# Figure 2: Single-Nucleus RNA-seq Integration and Cell Type Analysis
# ============================================================================
# This script performs:
#   - Figure 2A: UMAP visualization
#   - Figure 2B: Marker gene dot plot for cell type annotation
#   - Figure 2C: Cell type proportion analysis
#
# Author: Dae-Hwan Kim, Jibeom Ko
# Paper: "An integrative multi-omics framework identifies SERPINA3 as a 
#         circulating biomarker for cancer cachexia"
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(scales)

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment
#
# Public datasets used:
#   - GSE272085 (Zhang et al., Cell Reports 2024)
#   - Zenodo 11090497 (Agca et al., JCSM 2024)
#
# Required: 10X Genomics format directories with:
#   - matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz

DATA_DIR <- "path/to/your/data"
OUTPUT_DIR <- file.path(DATA_DIR, "results")

# Define paths to 10X data directories
CONTROL1_PATH <- file.path(DATA_DIR, "Con_filtered_feature_bc_matrix")
CONTROL2_PATH <- file.path(DATA_DIR, "ctrl_feature_bc_matrix")
CACHEXIA1_PATH <- file.path(DATA_DIR, "cac_filtered_feature_bc_matrix")
CACHEXIA2_PATH <- file.path(DATA_DIR, "ccx_feature_bc_matrix")

# Create output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Load 10X Data
# ----------------------------------------------------------------------------
cat("Loading 10X data...\n")

control.data1 <- Read10X(data.dir = CONTROL1_PATH)
control.data2 <- Read10X(data.dir = CONTROL2_PATH)
cac.data1 <- Read10X(data.dir = CACHEXIA1_PATH)
cac.data2 <- Read10X(data.dir = CACHEXIA2_PATH)

# ----------------------------------------------------------------------------
# 3. Create Seurat Objects
# ----------------------------------------------------------------------------
Control1 <- CreateSeuratObject(control.data1, project = "Control")
Control2 <- CreateSeuratObject(control.data2, project = "Control")
Cachexia1 <- CreateSeuratObject(cac.data1, project = "Cachexia")
Cachexia2 <- CreateSeuratObject(cac.data2, project = "Cachexia")

cat("Initial cell counts:\n")
cat("  Control1:", ncol(Control1), "\n")
cat("  Control2:", ncol(Control2), "\n")
cat("  Cachexia1:", ncol(Cachexia1), "\n")
cat("  Cachexia2:", ncol(Cachexia2), "\n")

# ----------------------------------------------------------------------------
# 4. Calculate QC Metrics
# ----------------------------------------------------------------------------
# Mitochondrial content
Control1[["percent.mt"]] <- PercentageFeatureSet(Control1, pattern = "^mt-")
Control2[["percent.mt"]] <- PercentageFeatureSet(Control2, pattern = "^mt-")
Cachexia1[["percent.mt"]] <- PercentageFeatureSet(Cachexia1, pattern = "^mt-")
Cachexia2[["percent.mt"]] <- PercentageFeatureSet(Cachexia2, pattern = "^mt-")

# Ribosomal content
Control1[["percent.ribo"]] <- PercentageFeatureSet(Control1, pattern = "^Rps|^Rpl")
Control2[["percent.ribo"]] <- PercentageFeatureSet(Control2, pattern = "^Rps|^Rpl")
Cachexia1[["percent.ribo"]] <- PercentageFeatureSet(Cachexia1, pattern = "^Rps|^Rpl")
Cachexia2[["percent.ribo"]] <- PercentageFeatureSet(Cachexia2, pattern = "^Rps|^Rpl")

# ----------------------------------------------------------------------------
# 5. Quality Control Filtering
# ----------------------------------------------------------------------------
# Filter criteria:
#   - nFeature_RNA > 200 and < 5000
#   - percent.mt < 5%

Control1 <- subset(Control1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
Control2 <- subset(Control2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
Cachexia1 <- subset(Cachexia1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
Cachexia2 <- subset(Cachexia2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

cat("\nAfter QC filtering:\n")
cat("  Control1:", ncol(Control1), "\n")
cat("  Control2:", ncol(Control2), "\n")
cat("  Cachexia1:", ncol(Cachexia1), "\n")
cat("  Cachexia2:", ncol(Cachexia2), "\n")

# ----------------------------------------------------------------------------
# 6. Doublet Detection and Removal (scDblFinder)
# ----------------------------------------------------------------------------
# Helper function for doublet removal
run_scDbl <- function(sobj, dbr = 0.07, seed = 1234) {
  set.seed(seed)
  sce <- as.SingleCellExperiment(sobj)
  sce <- scDblFinder(sce, dbr = dbr)
  sobj$scDblFinder.class <- colData(sce)$scDblFinder.class
  sobj$scDblFinder.score <- colData(sce)$scDblFinder.score
  subset(sobj, subset = scDblFinder.class == "singlet")
}

cat("\nRemoving doublets...\n")
Control1 <- run_scDbl(Control1, dbr = 0.07)
Control2 <- run_scDbl(Control2, dbr = 0.07)
Cachexia1 <- run_scDbl(Cachexia1, dbr = 0.07)
Cachexia2 <- run_scDbl(Cachexia2, dbr = 0.07)

cat("After doublet removal:\n")
cat("  Control1:", ncol(Control1), "\n")
cat("  Control2:", ncol(Control2), "\n")
cat("  Cachexia1:", ncol(Cachexia1), "\n")
cat("  Cachexia2:", ncol(Cachexia2), "\n")

# ----------------------------------------------------------------------------
# 7. Merge Datasets
# ----------------------------------------------------------------------------
# Merge by condition
Control_merged <- merge(Control1, y = Control2, add.cell.ids = c("Control1", "Control2"))
Cachexia_merged <- merge(Cachexia1, y = Cachexia2, add.cell.ids = c("Cachexia1", "Cachexia2"))

# Final merge
combined <- merge(Control_merged, y = Cachexia_merged, add.cell.ids = c("Control", "Cachexia"))

cat("\nTotal cells after merging:", ncol(combined), "\n")

# ----------------------------------------------------------------------------
# 8. SCTransform Normalization
# ----------------------------------------------------------------------------
options(future.globals.maxSize = 1e15)  # Increase memory limit
DefaultAssay(combined) <- "RNA"

cat("\nRunning SCTransform...\n")
combined <- SCTransform(combined, vars.to.regress = "percent.mt", verbose = FALSE)

# ----------------------------------------------------------------------------
# 9. Dimensionality Reduction (PCA)
# ----------------------------------------------------------------------------
cat("Running PCA...\n")
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)

# Visualize elbow plot
pdf(file.path(OUTPUT_DIR, "ElbowPlot.pdf"), width = 6, height = 4)
ElbowPlot(combined, ndims = 50)
dev.off()

# ----------------------------------------------------------------------------
# 10. Batch Correction with Harmony
# ----------------------------------------------------------------------------
cat("Running Harmony batch correction...\n")
combined <- RunHarmony(combined, group.by.vars = "orig.ident", dims = 1:10)

# ----------------------------------------------------------------------------
# 11. Clustering
# ----------------------------------------------------------------------------
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:10)
combined <- FindClusters(combined, resolution = 0.4)

# ----------------------------------------------------------------------------
# 12. UMAP Visualization
# ----------------------------------------------------------------------------
cat("Running UMAP...\n")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:10)

# Set condition factor levels
combined$orig.ident <- factor(combined$orig.ident, levels = c("Control", "Cachexia"))

# UMAP plot split by condition
umap_plot <- DimPlot(
  combined, 
  reduction = "umap", 
  label = TRUE, 
  split.by = "orig.ident", 
  label.size = 4.0, 
  pt.size = 0.3
) + plot_annotation(title = "Integrated snRNA-seq Data")

print(umap_plot)
ggsave(file.path(OUTPUT_DIR, "Figure2A_UMAP_split.pdf"), umap_plot, width = 12, height = 6)

# ----------------------------------------------------------------------------
# 13. Cell Type Annotation
# ----------------------------------------------------------------------------
# Prepare for marker identification
combined_marker <- PrepSCTFindMarkers(combined)

# Find all markers
all_markers <- FindAllMarkers(
  combined_marker, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

write.xlsx(all_markers, file.path(OUTPUT_DIR, "cluster_all_markers.xlsx"))

# ----------------------------------------------------------------------------
# 14. Figure 2B: Marker Gene Dot Plot
# ----------------------------------------------------------------------------
# Define canonical marker genes for cell type annotation
marker_genes <- c(
  # Type IIb myonuclei
  "Myh4", "Myh1", "Actn3",
  # Type IIa/IIx myonuclei
  "Myh2",
  # FAPs (Fibro-adipogenic progenitors)
  "Pdgfra", "Cd34",
  # MuSCs (Muscle stem cells)
  "Pax7", "Calcr",
  # Endothelial cells
  "Pecam1", "Cdh5", "Vwf",
  # Pericytes
  "Pdgfrb", "Rgs5", "Abcc9",
  # NMJ (Neuromuscular junction)
  "Chrne", "Musk", "Lrp4",
  # Macrophages
  "Adgre1", "Csf1r", "Mrc1",
  # MTJ (Myotendinous junction)
  "Col22a1", "Ankrd1", "Slc24a2"
)

DefaultAssay(combined) <- "SCT"

dotplot <- DotPlot(combined, features = marker_genes, cols = c("#FFFFCC", "blue3")) + 
  RotatedAxis() +
  labs(title = "Marker Gene Expression by Cluster")

print(dotplot)
ggsave(file.path(OUTPUT_DIR, "Figure2B_DotPlot_markers.pdf"), dotplot, width = 12, height = 8)

# ----------------------------------------------------------------------------
# 15. Assign Cell Type Identities
# ----------------------------------------------------------------------------
# Create copy for annotation
combined2 <- combined

# Define cluster to cell type mapping
# NOTE: Adjust this mapping based on your clustering results
new.cluster.ids <- c(
  "IIb", "IIb", "IIb", "FAPs", "IIa/IIx", "IIa/IIx", "FAPs", "IIb", 
  "IIa/IIx", "IIb", "IIa/IIx", "IIb", "Endothelial", "MTJ", "Pericyte", 
  "Macrophage", "MuSCs", "NMJ"
)

names(new.cluster.ids) <- levels(combined2)
combined2 <- RenameIdents(combined2, new.cluster.ids)

# Set cell type order for visualization
new_order <- c("IIb", "IIa/IIx", "FAPs", "MuSCs", "Endothelial", "Pericyte", "NMJ", "Macrophage", "MTJ")
combined2 <- SetIdent(combined2, value = factor(Idents(combined2), levels = new_order))

# Update metadata
combined2$orig.ident <- factor(combined2$orig.ident, levels = c("Control", "Cachexia"))
combined2$CellType <- Idents(combined2)

# Dot plot with annotated cell types
dotplot_annotated <- DotPlot(combined2, features = marker_genes, cols = c("#FFFFCC", "blue3")) + 
  RotatedAxis() +
  labs(title = "Marker Gene Expression by Cell Type")

print(dotplot_annotated)
ggsave(file.path(OUTPUT_DIR, "Figure2B_DotPlot_annotated.pdf"), dotplot_annotated, width = 12, height = 6)

# ----------------------------------------------------------------------------
# 16. Figure 2A: Final UMAP with Cell Type Labels
# ----------------------------------------------------------------------------
# Color scheme for cell types
cluster_colors <- c(
  "IIb"         = "#1f78b4",
  "IIa/IIx"     = "#a6cee3",
  "FAPs"        = "#4daf4a",
  "MuSCs"       = "#bebada",
  "Endothelial" = "#984ea3",
  "Pericyte"    = "#8dd3c7",
  "NMJ"         = "#bc80bd",
  "Macrophage"  = "#66c2a5",
  "MTJ"         = "#8c96c6"
)

# UMAP with cell type colors
umap_celltype <- DimPlot(
  combined2, 
  reduction = "umap",
  label = TRUE, 
  label.size = 3.5, 
  pt.size = 0.4,
  cols = cluster_colors
) + 
  theme(legend.position = "right") +
  labs(title = "Cell Type Annotation")

print(umap_celltype)
ggsave(file.path(OUTPUT_DIR, "Figure2A_UMAP_celltype.pdf"), umap_celltype, width = 8, height = 6)

# UMAP split by condition
umap_split <- DimPlot(
  combined2, 
  reduction = "umap", 
  split.by = "orig.ident",
  label = FALSE, 
  label.size = 3.5, 
  pt.size = 0.4,
  cols = cluster_colors
) + 
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

print(umap_split)
ggsave(file.path(OUTPUT_DIR, "Figure2A_UMAP_split_celltype.pdf"), umap_split, width = 12, height = 5)

# ----------------------------------------------------------------------------
# 17. Figure 2C: Cell Type Proportion Analysis
# ----------------------------------------------------------------------------
# Calculate cell type proportions
cell_stats <- combined2@meta.data %>%
  dplyr::select(orig.ident) %>%
  mutate(CellType = Idents(combined2)) %>%
  group_by(orig.ident, CellType) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(
    Total = sum(Count),
    Prop = Count / Total
  ) %>%
  ungroup()

# Set cell type order (reverse for stacking)
cell_stats$CellType <- factor(cell_stats$CellType, levels = rev(new_order))

# Stacked bar plot
proportion_plot <- ggplot(cell_stats, aes(x = orig.ident, y = Prop, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, color = "black", size = 0.2) +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  labs(
    title = "Cell Type Composition",
    x = NULL, 
    y = "Proportion (%)", 
    fill = "Cell Type"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

print(proportion_plot)
ggsave(file.path(OUTPUT_DIR, "Figure2C_CellType_proportion.pdf"), proportion_plot, width = 6, height = 6)

# ----------------------------------------------------------------------------
# 18. Feature Plots for Key Genes
# ----------------------------------------------------------------------------
# Visualize candidate biomarker genes
key_genes <- c("Apod", "Serpina3n", "Lox")

feature_plot <- FeaturePlot(
  combined2, 
  features = key_genes, 
  split.by = "orig.ident", 
  order = TRUE
) + 
  theme(legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "FeaturePlot_key_genes.pdf"), feature_plot, width = 12, height = 10)

# ----------------------------------------------------------------------------
# 19. Save Processed Data
# ----------------------------------------------------------------------------
cat("\nSaving processed Seurat object...\n")
saveRDS(combined2, file.path(OUTPUT_DIR, "snRNAseq_integrated.rds"))

# Save cell type proportions
write.xlsx(cell_stats, file.path(OUTPUT_DIR, "CellType_proportions.xlsx"))

# Session info
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "snRNAseq_session_info.txt"))

cat("\n============================================\n")
cat("Figure 2A-C analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
cat("\nCell counts by type:\n")
print(table(Idents(combined2)))
cat("============================================\n")
