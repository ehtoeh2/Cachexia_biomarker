# ============================================================================
# Figure 1F: Secretome Prediction and Candidate Biomarker Identification
# ============================================================================
# This script performs:
#   - Generation of secreted protein candidate heatmaps (Figure 1F)
#   - Signal peptide prediction using SignalP 5.0 
#   - Transmembrane domain filtering using DeepTMHMM 
#
# Author: Dae-Hwan Kim, Ji beom Ko
# Paper: "An integrative multi-omics framework identifies SERPINA3 as a 
#         circulating biomarker for cancer cachexia"
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Load Required Packages
# ----------------------------------------------------------------------------
library(dplyr)
library(biomaRt)
library(Biostrings)
library(openxlsx)
library(pheatmap)
library(tibble)

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment
#
# Required input:
#   - k_means_cluster: K-means clustering results from Figure 1 analysis
#   - z_ord: Z-score normalized expression matrix from Figure 1 analysis
#
# External tools required:
#   - SignalP 5.0 (https://services.healthtech.dtu.dk/service.php?SignalP-5.0)
#   - DeepTMHMM (https://dtu.biolib.com/DeepTMHMM)

DATA_DIR <- "path/to/your/data"
OUTPUT_DIR <- file.path(DATA_DIR, "results")

# ----------------------------------------------------------------------------
# 2. Filter Target Gene Clusters (GC2, GC3, GC4)
# ----------------------------------------------------------------------------
# Focus on biologically relevant clusters:
#   - GC2: Inflammatory (upregulated)
#   - GC3: Catabolic (upregulated)
#   - GC4: Muscle structural (downregulated)

target_genes <- k_means_cluster %>% 
  filter(cluster_label %in% c("GC2", "GC3", "GC4"))

gene_list <- unique(target_genes$gene)

# ----------------------------------------------------------------------------
# 3. Retrieve Protein Sequences from BioMart
# ----------------------------------------------------------------------------
# Connect to Ensembl database (for mouse; use "hsapiens_gene_ensembl" for human)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Query protein sequences based on gene symbols
seq_data <- getBM(
  attributes = c("external_gene_name", "peptide"),
  filters = "external_gene_name",
  values = gene_list,
  mart = mart
)

# Quality control: Remove empty or too short sequences
seq_data_clean <- seq_data %>%
  filter(peptide != "") %>%
  filter(nchar(peptide) > 30)  # SignalP requires minimum length

# Keep only the longest isoform per gene (canonical sequence)
seq_data_final <- seq_data_clean %>%
  group_by(external_gene_name) %>%
  dplyr::slice(which.max(nchar(peptide))) %>%
  ungroup()

# ----------------------------------------------------------------------------
# 4. Generate FASTA File for SignalP Analysis
# ----------------------------------------------------------------------------
# Create AAStringSet object for FASTA export
aa_set <- AAStringSet(seq_data_final$peptide)
names(aa_set) <- seq_data_final$external_gene_name

# Save FASTA file for SignalP web server upload
fasta_path <- file.path(OUTPUT_DIR, "GC234_candidates.fasta")
writeXStringSet(aa_set, filepath = fasta_path)

cat("\n============================================\n")
cat("FASTA file saved to:", fasta_path, "\n")
cat("============================================\n")
cat("\nNext steps:\n")
cat("1. Upload this FASTA file to SignalP 5.0 web server\n")
cat("   (https://services.healthtech.dtu.dk/service.php?SignalP-5.0)\n")
cat("2. Download the prediction results\n")
cat("3. Continue with Section 5 below\n")
cat("============================================\n")

# ----------------------------------------------------------------------------
# 5. Process SignalP Results
# ----------------------------------------------------------------------------
# NOTE: After running SignalP 5.0, download results and update the path below

# Load SignalP results (update path to your downloaded file)
sp_result_path <- file.path(OUTPUT_DIR, "SignalP_prediction_results.txt")

# For demonstration, assuming sp_result is loaded:
secreted_candidates <- sp_result %>%
  filter(prediction_class == "SP") %>%
  dplyr::select(gene, prediction_class)

# ----------------------------------------------------------------------------
# 6. Generate FASTA for TMHMM Analysis (Transmembrane Domain Filtering)
# ----------------------------------------------------------------------------
# Filter sequences for TMHMM analysis (only SP-positive proteins)
tmhmm_input_df <- seq_data_final %>%
  inner_join(secreted_candidates, by = c("external_gene_name" = "gene"))

tmhmm_aa_set <- AAStringSet(tmhmm_input_df$peptide)
names(tmhmm_aa_set) <- tmhmm_input_df$external_gene_name

# Save FASTA for TMHMM
tmhmm_fasta_path <- file.path(OUTPUT_DIR, "TMHMM_input_candidates.fasta")
writeXStringSet(tmhmm_aa_set, filepath = tmhmm_fasta_path)

cat("\n============================================\n")
cat("TMHMM input FASTA saved to:", tmhmm_fasta_path, "\n")
cat("============================================\n")
cat("\nNext steps:\n")
cat("1. Upload this FASTA file to DeepTMHMM\n")
cat("   (https://dtu.biolib.com/DeepTMHMM)\n")
cat("2. Download the prediction results\n")
cat("3. Continue with Section 7 below\n")
cat("============================================\n")

# ----------------------------------------------------------------------------
# 7. Process TMHMM Results
# ----------------------------------------------------------------------------
# NOTE: After running DeepTMHMM, download results and update the path below

# Load TMHMM results (update path to your downloaded file)
tmhmm_result_path <- file.path(OUTPUT_DIR, "TMHMM_result.txt")

# Filter for proteins WITHOUT transmembrane domains (PredHel = 0)
real_secreted_proteins <- tmhmm_df %>%
  filter(pred_hel_count == 0)

cat("Final secreted proteins (no TM domains):", nrow(real_secreted_proteins), "\n")

# Clean gene names (remove any trailing characters) if needed

# ----------------------------------------------------------------------------
# 8. Create Cluster Information for Final Candidates
# ----------------------------------------------------------------------------
# Match secreted candidates with cluster assignments
cluster_info <- k_means_cluster %>%
  filter(gene %in% target_secreted_genes_clean) %>%
  filter(cluster_label %in% c("GC2", "GC3", "GC4")) %>%
  mutate(cluster_label = factor(cluster_label, levels = c("GC2", "GC3", "GC4"))) %>%
  arrange(cluster_label, gene)


# ----------------------------------------------------------------------------
# 9. Figure 1F: Secreted Candidates Heatmap
# ----------------------------------------------------------------------------
# Subset expression matrix for secreted candidates
z_final <- z_ord[cluster_info$gene, , drop = FALSE]

# Row annotation
ann_row_final <- data.frame(Cluster = cluster_info$cluster_label)
rownames(ann_row_final) <- cluster_info$gene

# Column annotation (assuming ann_col2 is available from Figure 1 analysis)
# ann_col2 should have condition and batch information

# Function to plot individual cluster heatmaps (minimal version)
plot_cluster_heatmap <- function(target_cluster, z_matrix, cluster_df, output_path) {
  current_genes <- cluster_df %>%
    filter(cluster_label == target_cluster) %>%
    pull(gene)
  
  # Subset matrix
  z_sub <- z_matrix[current_genes, , drop = FALSE]
  
  # Plot minimal heatmap
  pdf(output_path)
  pheatmap(
    z_sub,
    scale = "none"
  )
  dev.off()
}

# Generate heatmaps for each cluster (GC2, GC3, GC4)
for (cluster in c("GC2", "GC3", "GC4")) {
  output_file <- file.path(OUTPUT_DIR, paste0("Figure1F_", cluster, "_heatmap.pdf"))
  plot_cluster_heatmap(
    target_cluster = cluster,
    z_matrix = z_final,
    cluster_df = cluster_info,
    output_path = output_file
  )
}



# ----------------------------------------------------------------------------
# 10. Save Final Secretome Lists
# ----------------------------------------------------------------------------
# Save secreted candidates by cluster
secretome_list <- list(
  GC2_secreted = cluster_info %>% filter(cluster_label == "GC2") %>% pull(gene),
  GC3_secreted = cluster_info %>% filter(cluster_label == "GC3") %>% pull(gene),
  GC4_secreted = cluster_info %>% filter(cluster_label == "GC4") %>% pull(gene)
)


cat("\n============================================\n")
cat("Figure 1F analysis completed!\n")
cat("============================================\n")
cat("\nSecretome candidates summary:\n")
cat("  GC2 (Inflammatory, Up):", length(secretome_list$GC2_secreted), "genes\n")
cat("  GC3 (Catabolic, Up):", length(secretome_list$GC3_secreted), "genes\n")
cat("  GC4 (Structural, Down):", length(secretome_list$GC4_secreted), "genes\n")
cat("  Total:", nrow(cluster_info), "secreted protein candidates\n")
cat("============================================\n")
