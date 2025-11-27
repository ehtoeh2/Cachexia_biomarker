# ============================================================================
# Figure 4: Human Serum Proteomics Validation (Cross-Species Translation)
# ============================================================================
# This script performs:
#   - Figure 4A: Human gastric cancer cachexia serum proteomics analysis
#   - Figure 4B: Composite biomarker score calculation
#   - Figure 4C: ROC curve analysis for biomarker panel validation
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
library(ggplot2)
library(pROC)
library(readxl)
library(stringr)
library(openxlsx)

set.seed(1234)

# ----------------------------------------------------------------------------
# 1. Data Paths Configuration
# ----------------------------------------------------------------------------
# NOTE: Modify these paths according to your local environment
#
# Required input:
#   - Human gastric cancer (GC) serum proteomics data
#   - Mouse prot_long data (for direction determination)

DATA_DIR <- "path/to/your/data"
OUTPUT_DIR <- file.path(DATA_DIR, "results")

# Load mouse proteomics data for direction (rho) calculation
# prot_long should have: Gene, sample, intensity, muscle_weight_mg columns
prot_long <- readRDS(file.path(OUTPUT_DIR, "mouse_prot_long.rds"))

# ----------------------------------------------------------------------------
# 2. Define Final Biomarker Panel (Mouse -> Human Orthologs)
# ----------------------------------------------------------------------------
# Mouse biomarker panel from multi-omics integration
panel_mouse <- c(
  "Ccl8", "Cfh", "Clec3b", "Cpq", "Ctsd", "Ctsl",
  "Fetub", "Fgl1", "Hp", "Itih4", "Lbp", "Lum",
  "Postn", "Prg4", "Saa1", "Selp", "Serpina3m",
  "Serpina3n", "Tfrc"
)

# Mouse to Human ortholog mapping
mouse2human <- c(
  Ccl8       = "CCL8",
  Cfh        = "CFH",
  Clec3b     = "CLEC3B",
  Cpq        = "CPQ",
  Ctsd       = "CTSD",
  Ctsl       = "CTSL",
  Fetub      = "FETUB",
  Fgl1       = "FGL1",
  Hp         = "HP",
  Itih4      = "ITIH4",
  Lbp        = "LBP",
  Lum        = "LUM",
  Postn      = "POSTN",
  Prg4       = "PRG4",
  Saa1       = "SAA1",
  Selp       = "SELP",
  Serpina3m  = "SERPINA3",  # Both map to SERPINA3
  Serpina3n  = "SERPINA3",
  Tfrc       = "TFRC"
)

cat("Mouse biomarkers:", length(panel_mouse), "\n")
cat("Unique human orthologs:", length(unique(mouse2human)), "\n")

# ----------------------------------------------------------------------------
# 3. Calculate Direction (Sign) from Mouse Data
# ----------------------------------------------------------------------------
# Convert to numeric
prot_long <- prot_long %>%
  mutate(
    intensity = as.numeric(intensity),
    muscle_weight_mg = as.numeric(muscle_weight_mg)
  )

# Filter for panel genes
prot_panel <- prot_long %>%
  filter(Gene %in% panel_mouse)

# Calculate Spearman correlation with muscle weight for each gene
candidates <- prot_panel %>%
  group_by(Gene) %>%
  summarise(
    n = sum(!is.na(intensity) & !is.na(muscle_weight_mg)),
    rho_mw = suppressWarnings(
      cor(intensity, muscle_weight_mg,
          method = "spearman", use = "complete.obs")
    ),
    .groups = "drop"
  ) %>%
  # Assign direction: negative rho means gene increases with cachexia
  mutate(sign = if_else(rho_mw < 0, 1, -1))

cat("\nDirection assignment based on muscle weight correlation:\n")
print(candidates)

# Create human panel with directions
panel_sign_mouse <- candidates %>%
  filter(Gene %in% names(mouse2human)) %>%
  mutate(human_gene = mouse2human[Gene]) %>%
  dplyr::select(mouse_gene = Gene, human_gene, sign) %>%
  distinct(human_gene, .keep_all = TRUE)  # Remove duplicates (SERPINA3)

cat("\nHuman panel with directions:\n")
print(panel_sign_mouse)

# ----------------------------------------------------------------------------
# 4. Load Human Gastric Cancer Proteomics Data
# ----------------------------------------------------------------------------
# Load TMT-labeled proteomics data
# Format: Excel file with ratio columns (e.g., 127/126, 128/126, etc.)

human_data_path <- file.path(DATA_DIR, "Human_GC_proteomics.xlsx")

raw <- read_excel(human_data_path, sheet = 1, skip = 2)

# Clean column names
hdr <- as.character(raw[1, ])
dat <- raw[-1, ]
colnames(dat) <- hdr

# Define sample columns (TMT ratio channels)
sample_cols <- c("126/126", "127/126", "128/126",
                 "129/126", "130/126", "131/126")

# Extract gene names from Description column
dat_num <- dat %>%
  mutate(
    across(all_of(sample_cols), as.numeric),
    Gene = str_extract(Description, "GN=[^ ]+"),
    Gene = str_remove(Gene, "^GN=")
  )

cat("\nHuman proteomics data:\n")
cat("Total proteins:", nrow(dat_num), "\n")

# Filter for panel genes present in human data
panel_use <- panel_sign_mouse %>%
  filter(human_gene %in% dat_num$Gene)

cat("Panel genes detected in human serum:", nrow(panel_use), "\n")
print(panel_use)

# ----------------------------------------------------------------------------
# 5. Prepare Data in Long Format
# ----------------------------------------------------------------------------
# Filter and reshape data
gaca_long <- dat_num %>%
  filter(Gene %in% panel_use$human_gene) %>%
  dplyr::select(Gene, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "channel",
    values_to = "ratio"
  ) %>%
  mutate(
    # Map channels to sample names
    sample = dplyr::recode(channel,
      "126/126" = "early4",
      "127/126" = "early5",
      "128/126" = "early6",
      "129/126" = "adv8",
      "130/126" = "adv11",
      "131/126" = "adv12"
    ),
    # Assign groups: Early stage vs Advanced stage
    group = if_else(str_detect(sample, "^early"), "Early", "Advanced"),
    group = factor(group, levels = c("Early", "Advanced"))
  )

cat("\nSamples per group:\n")
print(table(gaca_long$group) / nrow(panel_use))

# Convert to wide format (sample × gene matrix)
gaca_wide <- gaca_long %>%
  dplyr::select(sample, group, Gene, ratio) %>%
  distinct() %>%
  pivot_wider(
    names_from = Gene,
    values_from = ratio
  )

print(gaca_wide)

# ----------------------------------------------------------------------------
# 6. Calculate Composite Biomarker Score
# ----------------------------------------------------------------------------
# Extract expression matrix
expr_mat_h <- gaca_wide %>%
  dplyr::select(all_of(panel_use$human_gene)) %>%
  as.matrix()

# Z-score normalization (column-wise)
z_mat_h <- apply(expr_mat_h, 2, function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  (x - mu) / sdv
})
z_mat_h <- as.matrix(z_mat_h)

# Match sign vector order to matrix columns
panel_use_matched <- panel_use %>%
  arrange(match(human_gene, colnames(z_mat_h)))

sign_vec_h <- panel_use_matched$sign
K_total_h <- length(sign_vec_h)

# Apply sign to Z-scores
signed_z_h <- sweep(z_mat_h, 2, sign_vec_h, `*`)

# Calculate composite score (handling NA values)
K_i_h <- apply(!is.na(signed_z_h), 1, sum)
sum_i_h <- apply(signed_z_h, 1, function(x) sum(x, na.rm = TRUE))

# Final biomarker score formula: sum(signed_Z) / sqrt(K)
score_vec_h <- sum_i_h / sqrt(K_i_h)
score_vec_h[K_i_h == 0] <- NA

# Add scores to data
gaca_score <- gaca_wide %>%
  mutate(biomarker_score = score_vec_h)

cat("\nBiomarker scores by group:\n")
print(gaca_score %>% dplyr::select(sample, group, biomarker_score))

# ----------------------------------------------------------------------------
# 7. Figure 4B: Violin/Box Plot of Biomarker Scores
# ----------------------------------------------------------------------------
# Statistical test
tt <- t.test(biomarker_score ~ group, data = gaca_score)
p_lab <- paste0("P = ", signif(tt$p.value, 3))

cat("\nT-test result:\n")
print(tt)

# Create violin plot
p_violin_h <- ggplot(gaca_score,
                     aes(x = group, y = biomarker_score, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.1, outlier.size = 0.8) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("Early" = "#BEBEBE", "Advanced" = "#1A3664")) +
  labs(
    x = "",
    y = "Composite Biomarker Score",
    title = "Human Gastric Cancer Serum Validation"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# Add p-value annotation
y_max <- max(gaca_score$biomarker_score, na.rm = TRUE)
y_pos <- y_max + 0.15 * diff(range(gaca_score$biomarker_score, na.rm = TRUE))

p_violin_final <- p_violin_h +
  annotate("text",
    x = 1.5,
    y = y_pos,
    label = p_lab,
    size = 5,
    fontface = "bold"
  )

print(p_violin_final)
ggsave(file.path(OUTPUT_DIR, "Figure4B_Violin_human.pdf"), 
       p_violin_final, width = 5, height = 6)

# ----------------------------------------------------------------------------
# 8. Figure 4C: ROC Curve Analysis
# ----------------------------------------------------------------------------
# Filter valid samples
gaca_roc_df <- gaca_score %>%
  filter(!is.na(biomarker_score))

# Calculate ROC curve
roc_obj_h <- roc(
  response = gaca_roc_df$group,
  predictor = gaca_roc_df$biomarker_score,
  levels = c("Early", "Advanced"),
  direction = "<"  # Higher score = Advanced stage
)

# Get AUC and confidence interval
auc_val <- auc(roc_obj_h)
auc_ci <- ci.auc(roc_obj_h)

cat("\nROC Analysis:\n")
cat("AUC:", round(auc_val, 3), "\n")
cat("95% CI:", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), "\n")

# Create ROC plot with ggplot2
roc_df <- data.frame(
  sensitivity = roc_obj_h$sensitivities,
  specificity = roc_obj_h$specificities
)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(color = "#1A3664", size = 1.5) +
  geom_abline(linetype = "dashed", color = "gray50") +
  annotate("text", 
    x = 0.6, y = 0.3,
    label = paste0("AUC = ", round(auc_val, 3), 
                   "\n95% CI: ", round(auc_ci[1], 3), "-", round(auc_ci[3], 3)),
    size = 5, fontface = "bold"
  ) +
  labs(
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    title = "ROC Curve: Composite Biomarker Panel"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    aspect.ratio = 1
  ) +
  coord_equal()

print(roc_plot)
ggsave(file.path(OUTPUT_DIR, "Figure4C_ROC_human.pdf"), 
       roc_plot, width = 6, height = 6)

# Alternative: Base R ROC plot
pdf(file.path(OUTPUT_DIR, "Figure4C_ROC_base.pdf"), width = 6, height = 6)
plot(roc_obj_h, 
     print.auc = TRUE, 
     print.auc.x = 0.4, 
     print.auc.y = 0.2,
     col = "#1A3664",
     lwd = 2,
     main = "ROC Curve: Composite Biomarker Panel")
dev.off()

# ----------------------------------------------------------------------------
# 9. Individual Biomarker Analysis
# ----------------------------------------------------------------------------
cat("\nIndividual biomarker analysis:\n")

# Function to analyze individual biomarkers
analyze_biomarker <- function(gene_name, data_long, data_wide) {
  
  if (!gene_name %in% colnames(data_wide)) {
    return(NULL)
  }
  
  # Extract data
  gene_data <- data_wide %>%
    dplyr::select(sample, group, !!sym(gene_name)) %>%
    filter(!is.na(!!sym(gene_name)))
  
  colnames(gene_data)[3] <- "value"
  
  # T-test
  tt <- t.test(value ~ group, data = gene_data)
  
  # ROC analysis
  roc_result <- roc(
    response = gene_data$group,
    predictor = gene_data$value,
    levels = c("Early", "Advanced"),
    direction = "<"
  )
  
  return(list(
    gene = gene_name,
    p_value = tt$p.value,
    auc = as.numeric(auc(roc_result)),
    mean_early = mean(gene_data$value[gene_data$group == "Early"], na.rm = TRUE),
    mean_advanced = mean(gene_data$value[gene_data$group == "Advanced"], na.rm = TRUE)
  ))
}

# Analyze all panel genes
individual_results <- lapply(panel_use$human_gene, function(g) {
  analyze_biomarker(g, gaca_long, gaca_wide)
})

individual_df <- do.call(rbind, lapply(individual_results, as.data.frame))

cat("\nIndividual biomarker performance:\n")
print(individual_df %>% arrange(p_value))

# Save individual results
write.xlsx(individual_df, file.path(OUTPUT_DIR, "Individual_biomarker_analysis.xlsx"))

# ----------------------------------------------------------------------------
# 10. Combined Figure
# ----------------------------------------------------------------------------
# Combine violin and ROC plots
combined_fig <- p_violin_final + roc_plot +
  plot_annotation(
    title = "Human Gastric Cancer Serum Validation",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "Figure4_combined.pdf"), 
       combined_fig, width = 11, height = 5)

# ----------------------------------------------------------------------------
# 11. Save Results
# ----------------------------------------------------------------------------
# Save processed data
write.xlsx(
  list(
    Panel_info = panel_use,
    Scores = gaca_score,
    Individual = individual_df
  ),
  file.path(OUTPUT_DIR, "Human_validation_results.xlsx")
)

# Session info
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "Figure4_session_info.txt"))

cat("\n============================================\n")
cat("Figure 4 analysis completed!\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("============================================\n")
cat("\nKey findings:\n")
cat("  - Biomarker panel size:", nrow(panel_use), "genes\n")
cat("  - Group comparison P-value:", signif(tt$p.value, 3), "\n")
cat("  - Composite biomarker AUC:", round(auc_val, 3), "\n")
cat("============================================\n")
