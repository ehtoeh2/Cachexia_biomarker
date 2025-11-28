# ============================================================================
# Figure 1: Integrated Bulk RNA-seq Analysis for Cancer Cachexia
# ============================================================================
# This script performs:
#   - Figure 1B: PCA plot
#   - Figure 1C: Volcano plot  
#   - Figure 1D: K-means clustering heatmap (GC1-GC5)
#   - Figure 1E: GO Over-Representation Analysis (ORA) dot plot
#
# Author: Dae-Hwan Kim, Jibeom Ko
# Paper: "An integrative multi-omics framework identifies SERPINA3 as a 
#         circulating biomarker for cancer cachexia"
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
library(dplyr)
library(ggrepel)
library(tximport)
library(DESeq2)
library(pheatmap)
library(openxlsx)
library(sva)
library(enrichplot)
library(ggplot2)
library(tibble)
library(cluster)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment
#
# Required data files:
#   - TXNAME: Transcript ID file (tab-delimited, no header)
#   - SYMBOL: Gene symbol file (tab-delimited, no header)
#   - Sample .tsv files: Kallisto quantification outputs
#
# Public datasets used (GEO accessions):
#   - GSE65936, GSE123310, GSE138464, GSE142455

DATA_DIR <- "path/to/your/data"  # Change this to your data directory
TXNAME_PATH <- file.path(DATA_DIR, "TXNAME")
SYMBOL_PATH <- file.path(DATA_DIR, "SYMBOL")
KALLISTO_DIR <- file.path(DATA_DIR, "kallisto_outputs")
OUTPUT_DIR <- file.path(DATA_DIR, "results")

# Create output directory if not exists
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Load Transcript-to-Gene Mapping
# ----------------------------------------------------------------------------
TXNAME <- read.delim(TXNAME_PATH, header = FALSE)
SYMBOL <- read.delim(SYMBOL_PATH, header = FALSE)
tx2gene <- data.frame(TXNAME, SYMBOL)
colnames(tx2gene) <- c("TXNAME", "SYMBOL")

# ----------------------------------------------------------------------------
# 3. Filter for Protein-Coding Genes Only
# ----------------------------------------------------------------------------
# Connect to Ensembl database (requires internet connection)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Retrieve gene biotype information
gene_info <- getBM(
  attributes = c("external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = unique(tx2gene$SYMBOL),
  mart = mart
)

# Filter for protein-coding genes
coding_genes <- gene_info %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(external_gene_name)

tx2gene_coding <- tx2gene %>%
  filter(SYMBOL %in% coding_genes)

cat("Transcripts before filtering:", nrow(tx2gene), "\n")
cat("Transcripts after filtering (protein-coding only):", nrow(tx2gene_coding), "\n")

# ----------------------------------------------------------------------------
# 4. Define Sample Information
# ----------------------------------------------------------------------------
sample_names <- c(
  # --- Control samples (n=20) ---
  'Data_A_Cont1.tsv','Data_A_Cont2.tsv','Data_A_Cont3.tsv','Data_A_Cont4.tsv',
  "Data_B_Cont1.tsv","Data_B_Cont2.tsv","Data_B_Cont3.tsv","Data_B_Cont4.tsv",
  "Data_C_Cont1.tsv","Data_C_Cont2.tsv","Data_C_Cont3.tsv",
  "Data_D_Cont1.tsv", "Data_D_Cont2.tsv", "Data_D_Cont3.tsv",
  "Data_E_Cont1.tsv","Data_E_Cont2.tsv","Data_E_Cont3.tsv",
  "Data_E_Cont4.tsv","Data_E_Cont5.tsv","Data_E_Cont6.tsv",
  
  # --- Cachexia samples (n=19) ---
  "Data_A_Cachexia1.tsv","Data_A_Cachexia2.tsv","Data_A_Cachexia3.tsv","Data_A_Cachexia4.tsv",
  "Data_B_Cachexia1.tsv","Data_B_Cachexia2.tsv","Data_B_Cachexia3.tsv","Data_B_Cachexia4.tsv",
  "Data_C_Cachexia1.tsv","Data_C_Cachexia2.tsv","Data_C_Cachexia3.tsv",
  "Data_D_Cachexia1.tsv","Data_D_Cachexia2.tsv","Data_D_Cachexia3.tsv","Data_D_Cachexia4.tsv",
  "Data_E_Cachexia1.tsv","Data_E_Cachexia2.tsv","Data_E_Cachexia3.tsv", "Data_E_Cachexia4.tsv"
)

# Define file paths
files <- file.path(KALLISTO_DIR, sample_names)
names(files) <- sample_names

# Define batch information for batch effect correction
batch <- c(
  # Control batches
  rep("Batch_A", 4), rep("Batch_B", 4), rep("Batch_C", 3),
  rep("Batch_D", 3), rep("Batch_E", 6),
  # Cachexia batches
  rep("Batch_A", 4), rep("Batch_B", 4), rep("Batch_C", 3),
  rep("Batch_D", 4), rep("Batch_E", 4)
)

# Create sample table
sampleTable <- data.frame(
  condition = factor(rep(c("Control", "Cachexia"), times = c(20, 19))),
  batch = factor(batch)
)

# ----------------------------------------------------------------------------
# 5. Import Kallisto Quantification Data
# ----------------------------------------------------------------------------
txi.kallisto <- tximport(
  files, 
  type = 'kallisto', 
  tx2gene = tx2gene_coding, 
  ignoreAfterBar = TRUE, 
  ignoreTxVersion = TRUE
)

rownames(sampleTable) <- colnames(txi.kallisto$counts)

# ----------------------------------------------------------------------------
# 6. DESeq2 Analysis with Batch Correction
# ----------------------------------------------------------------------------
# Create DESeq2 dataset with batch + condition design
dds_batch <- DESeqDataSetFromTximport(
  txi.kallisto, 
  sampleTable,
  design = ~ batch + condition
)

# Filter low-count genes (keep genes with >= 10 counts in >= 19 samples)
smallestGroupSize <- 19
keep <- rowSums(counts(dds_batch) >= 10) >= smallestGroupSize
dds_batch <- dds_batch[keep, ]

# Set Control as reference level
dds_batch$condition <- relevel(dds_batch$condition, ref = "Control")

# Run DESeq2
deseq2.dds <- DESeq(dds_batch)
deseq2.res <- results(deseq2.dds)
deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]

# Convert to data frame
res_df <- data.frame(deseq2.res)

# Save DESeq2 results
write.xlsx(res_df, file.path(OUTPUT_DIR, "DESeq2_results.xlsx"), rowNames = TRUE)

# ----------------------------------------------------------------------------
# 7. Figure 1B: PCA Plot (with ComBat batch correction for visualization)
# ----------------------------------------------------------------------------
# VST normalization
vsd <- vst(dds_batch, blind = FALSE)
vst_mat <- assay(vsd)

# ComBat batch correction for visualization
mod_combat <- model.matrix(~ condition, data = sampleTable)
batch_var <- sampleTable$batch

combat_mat <- ComBat(
  dat = vst_mat, 
  batch = batch_var, 
  mod = mod_combat, 
  par.prior = TRUE, 
  prior.plots = FALSE
)

# PCA calculation
pca_res <- prcomp(t(combat_mat))
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# Create PCA data frame
pcaData <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  condition = sampleTable$condition,
  batch = sampleTable$batch,
  name = colnames(dds_batch)
)

# Color scheme
pca_colors <- c("Control" = "A", "Cachexia" = "B") 
pca_fill_colors <- c("Control" = "C", "Cachexia" = "D")

# Plot PCA (Figure 1B)
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, fill = condition)) 
print(pca_plot)
ggsave(file.path(OUTPUT_DIR, "Figure1B_PCA.pdf"), pca_plot, width = 8, height = 6)

# ----------------------------------------------------------------------------
# 8. Figure 1C: Volcano Plot
# ----------------------------------------------------------------------------
# Prepare data for volcano plot
plot_data <- res_df %>%
  as.data.frame() %>%
  rownames_to_column("gene_symbol") %>%
  filter(!is.na(padj)) %>%
  mutate(gene = case_when(
    log2FoldChange >= 1 & padj < 0.05 ~ "Upregulated",
    log2FoldChange <= -1 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Non significant"
  ))


# Color scheme
volcano_colors <- c(
  "Downregulated" = 'A',
  "Upregulated" = 'B',
  "Non significant" = 'C'
)

# Plot Volcano (Figure 1C)

volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)")


# ----------------------------------------------------------------------------
# 9. Figure 1D: K-means Clustering Heatmap
# ----------------------------------------------------------------------------
# Filter significant DEGs
sig_genes <- deseq2.res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
  pull(gene)


# Subset and scale expression matrix
mat_deg <- combat_mat[intersect(rownames(combat_mat), sig_genes), , drop = FALSE]
keep_sd <- apply(mat_deg, 1, sd) > 0
mat_deg <- mat_deg[keep_sd, , drop = FALSE]

# Z-score normalization
z <- t(scale(t(mat_deg)))

# Determine optimal k using elbow 
wss <- sapply(2:10, function(k) {
  kmeans(z, centers = k, nstart = 100, iter.max = 100)$tot.withinss
})

elbow_df <- data.frame(
  k = 2:10,
  wss = wss)

# K-means clustering (k = 5)
set.seed(1234)
k <- 5
km <- kmeans(z, centers = k, nstart = 200, iter.max = 100)

# Create cluster mapping table
cluster_tbl <- tibble(
  gene = names(km$cluster),
  cluster_num = km$cluster,
  cluster_label = paste0("GC", km$cluster)
) %>%
  arrange(cluster_num) %>%
  mutate(cluster_label = factor(cluster_label, levels = paste0("GC", 1:k)))

# Order expression matrix by cluster
z_ord <- z[cluster_tbl$gene, , drop = FALSE]

# Prepare column annotation
ann_col <- sampleTable %>% 
  dplyr::select(condition, batch) %>%
  mutate(
    condition = factor(condition, levels = c("Control", "Cachexia")),
    batch = factor(batch)
  ) %>%
  arrange(condition, batch)

z_ord <- z_ord[, rownames(ann_col), drop = FALSE]

# Prepare row annotation
ann_row <- data.frame(Cluster = cluster_tbl$cluster_label)
rownames(ann_row) <- cluster_tbl$gene


# Heatmap (Figure 1D) - minimal 
pheatmap(
  z_ord,
  scale = "none"
)


# ----------------------------------------------------------------------------
# 10. Figure 1E: GO Over-Representation Analysis (ORA)
# ----------------------------------------------------------------------------
# Run compareCluster for all gene clusters
ora_all <- compareCluster(
  gene ~ cluster_label,
  data = cluster_tbl,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Save full ORA results
df_ora <- as.data.frame(ora_all)

# Select representative GO terms for visualization
picked <- tribble(
  ~Cluster, ~Description,
  "GC1", "axon guidance",
  "GC1", "neuron projection guidance",
  "GC1", "amide metabolic process",
  "GC2", "acute inflammatory response",
  "GC2", "leukocyte migration",
  "GC2", "regulation of inflammatory response",
  "GC3", "regulation of cellular catabolic process",
  "GC3", "regulation of protein ubiquitination",
  "GC3", "regulation of autophagy",
  "GC4", "muscle system process",
  "GC4", "muscle contraction",
  "GC4", "muscle cell differentiation",
  "GC5", "extracellular matrix organization",
  "GC5", "external encapsulating structure organization",
  "GC5", "extracellular structure organization"
)

# Filter ORA results for selected terms
cmpdf_sub <- df_ora %>%
  inner_join(picked, by = c("Cluster", "Description"))

target_cluster_order <- c("GC1", "GC2", "GC3", "GC4", "GC5")
cmpdf_sub$Cluster <- factor(cmpdf_sub$Cluster, levels = target_cluster_order)

real_desc_order <- picked$Description[picked$Description %in% cmpdf_sub$Description]
cmpdf_sub$Description <- factor(cmpdf_sub$Description, levels = rev(real_desc_order))

# Create subset object for plotting
cmp_sub <- ora_all
cmp_sub@compareClusterResult <- cmpdf_sub

# Minimal ORA Dotplot (Figure 1E)
dotplot(cmp_sub)

# ----------------------------------------------------------------------------
# Session info
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))

cat("\n============================================\n")
cat("Figure 1 analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
