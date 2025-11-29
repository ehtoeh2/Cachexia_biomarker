# SERPINA3 Cachexia Biomarker Analysis

## An integrative multi-omics framework identifies SERPINA3 as a circulating biomarker for cancer cachexia

**Authors:** Dae-Hwan Kim, Jibeom Ko, Yongju Kim, Kwan Hyi Lee, Jihoon Kim

**Affiliation:** Korea Institute of Science and Technology (KIST)

---

## Overview

This repository contains analysis code for identifying and validating circulating biomarkers for cancer cachexia using an integrative multi-omics approach. The workflow combines:

1. **Bulk RNA-seq analysis** - Meta-analysis of 4 independent cachexia studies
2. **Secretome prediction** - SignalP/TMHMM-based filtering for secreted proteins
3. **Single-nucleus RNA-seq** - Cell type-specific expression analysis
4. **NicheNet analysis** - Ligand-receptor interaction networks
5. **Mouse serum proteomics** - Mass spectrometry validation
6. **Human serum proteomics** - Cross-species translation and clinical validation

---

## Repository Structure

```
SERPINA3_cachexia_biomarker/
│
├── README.md                          # This file
│
├── Figure1/                           # Bulk RNA-seq & Secretome Analysis
│   ├── 01_BulkRNAseq_DESeq2_analysis.R   # Fig 1B-E: PCA, Volcano, Heatmap, ORA
│   └── 02_Secretome_prediction.R          # Fig 1F: SignalP + TMHMM filtering
│
├── Figure2/                           # Single-nucleus RNA-seq Analysis
│   ├── 01_snRNAseq_integration.R         # Fig 2A-C: UMAP, DotPlot, Proportions
│   └── 02_DEG_NicheNet_analysis.R        # Fig 2D-E: DEG heatmap, NicheNet
│
├── Figure3/                           # Mouse Serum Proteomics
│   └── 01_Mouse_serum_proteomics.R       # Fig 3H-I: Heatmaps, Stouffer, Lollipop
│
├── Figure4/                           # Human Serum Validation
│   └── 01_Human_serum_validation.R       # Fig 4A-C: PCA plot, volcano plot
│
└── data/                              # Data directory (not included)
    └── README_data.md                 # Data availability information
```

---

## Figure-to-Code Mapping

| Figure | Panel | Description | Script |
|--------|-------|-------------|--------|
| **1** | A | Schematic workflow | N/A (Illustration) |
| | B | PCA plot | `Figure1/01_BulkRNAseq_DESeq2_analysis.R` |
| | C | Volcano plot | `Figure1/01_BulkRNAseq_DESeq2_analysis.R` |
| | D | K-means clustering heatmap | `Figure1/01_BulkRNAseq_DESeq2_analysis.R` |
| | E | GO ORA dot plot | `Figure1/01_BulkRNAseq_DESeq2_analysis.R` |
| | F | Secreted candidates heatmap | `Figure1/02_Secretome_prediction.R` |
| **2** | A | UMAP visualization | `Figure2/01_snRNAseq_integration.R` |
| | B | Marker gene dot plot | `Figure2/01_snRNAseq_integration.R` |
| | C | Cell type proportions | `Figure2/01_snRNAseq_integration.R` |
| | D | Cell type-specific DEG heatmap | `Figure2/02_DEG_NicheNet_analysis.R` |
| | E | NicheNet chord diagram | `Figure2/02_DEG_NicheNet_analysis.R` |
| **3** | G | Proteomics volcanoplot | `Figure3/01_Mouse_serum_proteomics.R` |
| | H | Multi-omics integration | `Figure3/01_Mouse_serum_proteomics.R` |
| | I | Concordant biomarker lollipop | `Figure3/01_Mouse_serum_proteomics.R` |
| **4** | A | Human proteomics PCA analysis | `Figure4/01_Human_serum_validation.R` |
| | B | Human proteomics volcano plot | `Figure4/01_Human_serum_validation.R` |
| | C | Scatter plot for cross-species concordant | `Figure4/01_Human_serum_validation.R` |

---

## Data Availability

### Public Datasets

| Dataset | GEO/Zenodo | Description | Reference |
|---------|------------|-------------|-----------|
| GSE65936 | GEO | Bulk RNA-seq, Lewis lung carcinoma | Judge et al., JCI 2014 |
| GSE123310 | GEO | Bulk RNA-seq, Pancreatic cancer | Tseng et al., JNCI 2020 |
| GSE138464 | GEO | Bulk RNA-seq, KPC model | Rupert et al., JEM 2021 |
| GSE142455 | GEO | Bulk RNA-seq, C26 model | Fontes-Oliveira et al., EMBO 2022 |
| GSE272085 | GEO | snRNA-seq, Cachexia muscle | Zhang et al., Cell Reports 2024 |
| Zenodo 11090497 | Zenodo | snRNA-seq, Cachexia muscle | Agca et al., JCSM 2024 |

### Proteomics Data

Mouse serum proteomics data will be provided PRIDE after publication.
Human serum proteomics data | PXD035832 (Public Datasets)


### External Resources

- **SignalP 5.0**: https://services.healthtech.dtu.dk/service.php?SignalP-5.0
- **DeepTMHMM**: https://dtu.biolib.com/DeepTMHMM
- **NicheNet**: https://github.com/saeyslab/nichenetr
- **Mouse L-R network**: https://zenodo.org/record/7074291

---

## Software Requirements

### R Version
- R >= 4.2.0

### Required R Packages

```r
# Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "tximport",
  "biomaRt",
  "clusterProfiler",
  "org.Mm.eg.db",
  "ComplexHeatmap",
  "SingleCellExperiment",
  "scDblFinder",
  "enrichplot"
))

# CRAN packages
install.packages(c(
  "Seurat",
  "harmony",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "circlize",
  "patchwork",
  "openxlsx",
  "sva",
  "pROC",
  "readxl",
  "stringr",
  "tibble",
  "scales"
))

# NicheNet (GitHub)
devtools::install_github("saeyslab/nichenetr")
```

---

## Usage

### 1. Clone Repository
```bash
git clone https://github.com/[username]/SERPINA3_cachexia_biomarker.git
cd SERPINA3_cachexia_biomarker
```

### 2. Prepare Data
- Download public datasets from GEO/Zenodo
- Place raw data in appropriate directories
- Update data paths in each script

### 3. Run Analysis (Sequential Order)

```r
# Step 1: Bulk RNA-seq analysis
source("Figure1/01_BulkRNAseq_DESeq2_analysis.R")

# Step 2: Secretome prediction
source("Figure1/02_Secretome_prediction.R")
# Note: Requires manual upload to SignalP/TMHMM web servers

# Step 3: snRNA-seq integration
source("Figure2/01_snRNAseq_integration.R")

# Step 4: DEG and NicheNet analysis
source("Figure2/02_DEG_NicheNet_analysis.R")

# Step 5: Mouse proteomics validation
source("Figure3/01_Mouse_serum_proteomics.R")

# Step 6: Human validation
source("Figure4/01_Human_serum_validation.R")
```

---

## Key Analysis Parameters

### Bulk RNA-seq
- **Batch correction**: ComBat (sva package)
- **DESeq2 design**: `~ batch + condition`
- **Gene filtering**: ≥10 counts in ≥19 samples
- **DEG criteria**: |log2FC| > 1, adjusted P < 0.05
- **K-means clustering**: k = 5, nstart = 25

### snRNA-seq
- **Quality control**: 200 < nFeature < 5000, percent.mt < 5%
- **Doublet detection**: scDblFinder (dbr = 0.07)
- **Normalization**: SCTransform (regress percent.mt)
- **Batch correction**: Harmony (dims = 1:10)
- **Clustering**: resolution = 0.4

### Multi-omics Integration
- **Concordance**: Same direction in RNA-seq and proteomics
- **Meta-analysis**: Stouffer's Z-score method
- **Combined Z**: sum(signed_Z) / sqrt(K)

---

## Output Files

Each analysis script generates output in the `results/` directory:

- `*.pdf` - Publication-ready figures
- `*.xlsx` - Data tables and gene lists
- `*.rds` - R objects for downstream analysis
- `*_session_info.txt` - Package versions

---

## Citation

If you use this code, please cite:

```
Kim DH, Ko J, Kim Y, Lee KH, Kim J. An integrative multi-omics framework 
identifies SERPINA3 as a circulating biomarker for cancer cachexia. 
[Journal Name]. 2025.
```

---

## Contact

For questions or issues, please contact:
- Dae-Hwan Kim: [email]
- Jihoon Kim: [email]

Korea Institute of Science and Technology (KIST)

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.
