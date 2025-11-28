# ============================================================================
# Figure 2D-E: Cell Type-Specific DEG Analysis and NicheNet Interaction
# ============================================================================
# This script performs:
#   - Figure 2D: Cell type-specific DEG heatmap for secreted candidates
#   - Figure 2E: NicheNet ligand-receptor interaction analysis
#
# Author: Dae-Hwan Kim, Ji beom Ko
# Paper: "An integrative multi-omics framework identifies SERPINA3 as a 
#         circulating biomarker for cancer cachexia"
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(nichenetr)
library(openxlsx)

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Load Data
# ----------------------------------------------------------------------------
# NOTE: Modify paths according to your environment

DATA_DIR <- "path/to/your/data"
OUTPUT_DIR <- file.path(DATA_DIR, "results")

# Load integrated snRNA-seq object from Figure 2A-C analysis
combined2 <- readRDS(file.path(OUTPUT_DIR, "snRNAseq_integrated.rds"))

# Load secretome lists from Figure 1F analysis
# These should contain genes from GC2, GC3, GC4 clusters
GC2_genes <- readRDS(file.path(OUTPUT_DIR, "secretome_list.rds"))$GC2_secreted
GC3_genes <- readRDS(file.path(OUTPUT_DIR, "secretome_list.rds"))$GC3_secreted
GC4_genes <- readRDS(file.path(OUTPUT_DIR, "secretome_list.rds"))$GC4_secreted

# ----------------------------------------------------------------------------
# 2. Subset to Target Cell Types
# ----------------------------------------------------------------------------
# Focus on major populations: Type II myonuclei (IIb, IIa/IIx) and FAPs

# Store cell type in metadata
combined2$CellType <- Idents(combined2)
Idents(combined2) <- "CellType"

# Subset to target populations
target_cells <- c("IIb", "IIa/IIx", "FAPs")
subset_obj <- subset(combined2, idents = target_cells)


# ----------------------------------------------------------------------------
# 3. Cell Type-Specific Differential Expression Analysis
# ----------------------------------------------------------------------------
# Function to run DEG analysis for a specific cell type
run_celltype_DEG <- function(seurat_obj, celltype_name) {
  # Subset to cell type
  cell_obj <- subset(seurat_obj, idents = celltype_name)
  
  # Prepare for FindMarkers
  DefaultAssay(cell_obj) <- "RNA"
  cell_obj <- JoinLayers(cell_obj)
  cell_obj <- NormalizeData(cell_obj)
  
  # Set condition as identity for comparison
  Idents(cell_obj) <- "orig.ident"
  
  # Find DEGs: Cachexia vs Control
  markers <- FindMarkers(
    cell_obj, 
    ident.1 = "Cachexia", 
    ident.2 = "Control",
    assay = "RNA",
    logfc.threshold = 1, 
    min.pct = 0.01
  )
  
  return(markers)
}


# Run for each cell type
IIb_markers <- run_celltype_DEG(subset_obj, "IIb")
IIa_IIx_markers <- run_celltype_DEG(subset_obj, "IIa/IIx")
FAP_markers <- run_celltype_DEG(subset_obj, "FAPs")

# ----------------------------------------------------------------------------
# 4. Filter Significant DEGs
# ----------------------------------------------------------------------------
# Apply stringent criteria: adjusted P < 0.05, |log2FC| > 1

IIb_sig <- IIb_markers %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)

IIa_IIx_sig <- IIa_IIx_markers %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)

FAP_sig <- FAP_markers %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)


# ----------------------------------------------------------------------------
# 5. Identify Cell Type-Specific Secreted Candidates (Survivors)
# ----------------------------------------------------------------------------
# Combine myonuclei signatures
myonuclei_sig_genes <- unique(c(rownames(IIb_sig), rownames(IIa_IIx_sig)))
fap_sig_genes <- rownames(FAP_sig)

# Intersect with GC clusters
# GC2 survivors
FAP_gc2_survivors <- intersect(GC2_genes, fap_sig_genes)
Myonuclei_GC2_survivors <- intersect(GC2_genes, myonuclei_sig_genes)

# GC3 survivors
FAP_gc3_survivors <- intersect(GC3_genes, fap_sig_genes)
Myonuclei_GC3_survivors <- intersect(GC3_genes, myonuclei_sig_genes)

# GC4 survivors
FAP_gc4_survivors <- intersect(GC4_genes, fap_sig_genes)
Myonuclei_GC4_survivors <- intersect(GC4_genes, myonuclei_sig_genes)


# Save survivor lists
survivor_list <- list(
  FAP_GC2 = FAP_gc2_survivors,
  FAP_GC3 = FAP_gc3_survivors,
  FAP_GC4 = FAP_gc4_survivors,
  Myo_GC2 = Myonuclei_GC2_survivors,
  Myo_GC3 = Myonuclei_GC3_survivors,
  Myo_GC4 = Myonuclei_GC4_survivors
)
saveRDS(survivor_list, file.path(OUTPUT_DIR, "CellType_survivors.rds"))

# ----------------------------------------------------------------------------
# 6. Figure 2D: Cell Type-Specific Heatmap
# ----------------------------------------------------------------------------

# Create gene grouping information
gene_group_vec <- character(length(all_survivors))
names(gene_group_vec) <- all_survivors

gene_group_vec[c(FAP_gc2_survivors, Myonuclei_GC2_survivors)] <- "GC2"
gene_group_vec[c(FAP_gc3_survivors, Myonuclei_GC3_survivors)] <- "GC3"
gene_group_vec[c(FAP_gc4_survivors, Myonuclei_GC4_survivors)] <- "GC4"

# Create group ID for heatmap
subset_obj$group_id <- paste(Idents(subset_obj), subset_obj$orig.ident, sep = "-")

# Calculate average expression
DefaultAssay(subset_obj) <- "RNA"
avg_exp <- AverageExpression(
  subset_obj, 
  features = all_survivors, 
  group.by = "group_id", 
  assay = "RNA", 
  layer = "data"
)$RNA

# Z-score normalization
scaled_mat <- t(scale(t(avg_exp)))
scaled_mat[is.na(scaled_mat)] <- 0

# Define column order
col_order <- c(
  "IIb-Control", "IIb-Cachexia",
  "IIa/IIx-Control", "IIa/IIx-Cachexia",
  "FAPs-Control", "FAPs-Cachexia"
)
valid_cols <- intersect(col_order, colnames(scaled_mat))
scaled_mat_final <- scaled_mat[, valid_cols]

# Column split vector
split_vec <- c("IIb", "IIb", "IIa/IIx", "IIa/IIx", "FAPs", "FAPs")

# Cell size settings
cell_width <- unit(9, "mm")
cell_height <- unit(3, "mm")

#minimal heatmap 
Heatmap(
  scaled_mat_final,
  row_split = gene_group_vec[rownames(scaled_mat_final)],
  column_split = split_vec
)

# ----------------------------------------------------------------------------
# 7. Figure 2E: NicheNet Ligand-Receptor Interaction Analysis
# ----------------------------------------------------------------------------
# Load mouse ligand-receptor network
# Download from: https://zenodo.org/record/7074291
lr_network_path <- file.path(DATA_DIR, "lr_network_mouse_21122021.rds")

if (file.exists(lr_network_path)) {
  lr_network <- readRDS(lr_network_path)
} else {
  cat("\nWARNING: Ligand-receptor network file not found.\n")
  cat("Please download from Zenodo and update the path.\n")
  cat("Skipping NicheNet analysis...\n")
}

# --- Part A: FAP to Myonuclei signaling ---

# FAP-derived ligands
fap_ligands_all <- unique(c(FAP_gc2_survivors, FAP_gc3_survivors, FAP_gc4_survivors))

# Get expressed receptors in Myonuclei
Idents(subset_obj) <- "CellType"
myo_obj <- subset(subset_obj, idents = c("IIb", "IIa/IIx"))
DefaultAssay(myo_obj) <- "RNA"
myo_obj <- JoinLayers(myo_obj)

myo_data <- GetAssayData(myo_obj, assay = "RNA", layer = "data")
expressed_receptors_myo <- rownames(myo_data)[rowSums(myo_data > 0) > (ncol(myo_obj) * 0.05)]

# Filter valid interactions
valid_interactions_fap_myo <- lr_network %>%
  filter(from %in% fap_ligands_all) %>%
  filter(to %in% expressed_receptors_myo) %>%
  distinct(from, to)

# Add known interaction (Serpina3n -> Lrp1)
manual_pairs <- data.frame(from = "Serpina3n", to = "Lrp1")
final_interactions_fap_myo <- bind_rows(valid_interactions_fap_myo[, c("from", "to")], manual_pairs) %>%
  distinct()

cat("\nFAP -> Myonuclei interactions:", nrow(final_interactions_fap_myo), "\n")

# --- Part B: Myonuclei to FAP signaling ---

# Myonuclei-derived ligands
myo_ligands_all <- unique(c(Myonuclei_GC2_survivors, Myonuclei_GC3_survivors, Myonuclei_GC4_survivors))

# Get expressed receptors in FAPs
fap_obj <- subset(subset_obj, idents = "FAPs")
DefaultAssay(fap_obj) <- "RNA"
fap_obj <- JoinLayers(fap_obj)

fap_data <- GetAssayData(fap_obj, assay = "RNA", layer = "data")
expressed_receptors_fap <- rownames(fap_data)[rowSums(fap_data > 0) > (ncol(fap_obj) * 0.05)]

# Filter valid interactions
valid_interactions_myo_fap <- lr_network %>%
  filter(from %in% myo_ligands_all) %>%
  filter(to %in% expressed_receptors_fap) %>%
  distinct(from, to)

cat("Myonuclei -> FAP interactions:", nrow(valid_interactions_myo_fap), "\n")

# ----------------------------------------------------------------------------
# 8. Create Chord Diagrams (Figure 2E)
# ----------------------------------------------------------------------------
# Ligand information with GC groups
ligand_info_fap <- bind_rows(
  data.frame(gene = FAP_gc2_survivors, group = "GC2 (Inflammatory)"),
  data.frame(gene = FAP_gc3_survivors, group = "GC3 (Catabolic)"),
  data.frame(gene = FAP_gc4_survivors, group = "GC4 (Structural)")
) %>% distinct(gene, .keep_all = TRUE)

# Add Serpina3n if not present
if (!"Serpina3n" %in% ligand_info_fap$gene) {
  ligand_info_fap <- bind_rows(
    ligand_info_fap,
    data.frame(gene = "Serpina3n", group = "GC2 (Inflammatory)")
  ) %>% distinct(gene, .keep_all = TRUE)
}

# Color scheme
gc_colors <- c(
  "GC2 (Inflammatory)" = "#EEA2AD",
  "GC3 (Catabolic)" = "#FF5B5B",
  "GC4 (Structural)" = "#87CEEB"
)
receptor_color <- "#4DAF4A"

# --- FAP to Myonuclei Chord Diagram ---
plot_df_fap_myo <- final_interactions_fap_myo

# Prepare colors
all_genes_plot <- unique(c(plot_df_fap_myo$from, plot_df_fap_myo$to))
grid_col <- structure(rep(receptor_color, length(all_genes_plot)), names = all_genes_plot)

for (i in 1:nrow(ligand_info_fap)) {
  gene_name <- ligand_info_fap$gene[i]
  if (gene_name %in% names(grid_col)) {
    grid_col[gene_name] <- gc_colors[ligand_info_fap$group[i]]
  }
}

# Sector order
ordered_ligands <- ligand_info_fap$gene[ligand_info_fap$gene %in% plot_df_fap_myo$from]
ordered_receptors <- unique(plot_df_fap_myo$to)
sector_order <- c(ordered_ligands, ordered_receptors)

# Draw chord diagram
pdf(file.path(OUTPUT_DIR, "Figure2E_ChordDiagram_FAP_to_Myo.pdf"), width = 10, height = 10)

circos.clear()
gaps <- c(rep(1, length(ordered_ligands) - 1), 10, rep(1, length(ordered_receptors) - 1), 10)
circos.par(gap.after = gaps, start.degree = 90)

chordDiagram(
  plot_df_fap_myo,
  order = sector_order,
  grid.col = grid_col,
  transparency = 0.2,
  annotationTrack = "grid",
  preAllocateTracks = 1,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow"
)

# Add gene names
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

# Add labels
text(-0.9, 0.5, "FAP Ligands\n(By GC Group)", cex = 1.2, font = 2)
text(0.9, 0.5, "Myonuclei\nReceptors", cex = 1.2, font = 2)

# Legend
legend("bottomleft", legend = names(gc_colors), fill = gc_colors, 
       title = "Ligand Groups", bty = "n", cex = 0.8)

circos.clear()
dev.off()

# --- Myonuclei to FAP Chord Diagram ---
if (nrow(valid_interactions_myo_fap) > 0) {
  
  # Ligand info for myonuclei
  myo_ligand_info <- bind_rows(
    data.frame(gene = Myonuclei_GC2_survivors, group = "GC2 (Inflammatory)"),
    data.frame(gene = Myonuclei_GC3_survivors, group = "GC3 (Catabolic)"),
    data.frame(gene = Myonuclei_GC4_survivors, group = "GC4 (Structural)")
  ) %>% distinct(gene, .keep_all = TRUE)
  
  plot_df_myo_fap <- valid_interactions_myo_fap[, c("from", "to")]
  
  # Colors
  receptor_color_fap <- "#FFAE42"
  all_genes_rev <- unique(c(plot_df_myo_fap$from, plot_df_myo_fap$to))
  grid_col_rev <- structure(rep(receptor_color_fap, length(all_genes_rev)), names = all_genes_rev)
  
  for (i in 1:nrow(myo_ligand_info)) {
    gene_name <- myo_ligand_info$gene[i]
    if (gene_name %in% names(grid_col_rev)) {
      grid_col_rev[gene_name] <- gc_colors[myo_ligand_info$group[i]]
    }
  }
  
  ordered_ligands_rev <- myo_ligand_info$gene[myo_ligand_info$gene %in% plot_df_myo_fap$from]
  ordered_receptors_rev <- unique(plot_df_myo_fap$to)
  sector_order_rev <- c(ordered_ligands_rev, ordered_receptors_rev)
  
  pdf(file.path(OUTPUT_DIR, "Figure2E_ChordDiagram_Myo_to_FAP.pdf"), width = 10, height = 10)
  
  circos.clear()
  gaps_rev <- c(rep(1, length(ordered_ligands_rev) - 1), 10, rep(1, length(ordered_receptors_rev) - 1), 10)
  circos.par(gap.after = gaps_rev, start.degree = 90)
  
  chordDiagram(
    plot_df_myo_fap,
    order = sector_order_rev,
    grid.col = grid_col_rev,
    transparency = 0.2,
    annotationTrack = "grid",
    preAllocateTracks = 1,
    directional = 1,
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow"
  )
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
  }, bg.border = NA)
  
  text(-0.9, 0.5, "Myonuclei Ligands\n(By GC Group)", cex = 1.2, font = 2)
  text(0.9, 0.5, "FAP\nReceptors", cex = 1.2, font = 2)
  
  legend("bottomleft", legend = names(gc_colors), fill = gc_colors, 
         title = "Ligand Groups", bty = "n", cex = 0.8)
  
  circos.clear()
  dev.off()
}

# ----------------------------------------------------------------------------
# 9. Save Results
# ----------------------------------------------------------------------------
# Save interaction data
write.xlsx(
  list(
    FAP_to_Myo = final_interactions_fap_myo,
    Myo_to_FAP = valid_interactions_myo_fap
  ),
  file.path(OUTPUT_DIR, "NicheNet_interactions.xlsx")
)

cat("\n============================================\n")
cat("Figure 2D-E analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
