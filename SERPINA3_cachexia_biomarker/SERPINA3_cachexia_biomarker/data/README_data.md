# Data Directory

This directory should contain the input data files for the analysis.

## Required Data Structure

```
data/
├── bulk_rnaseq/
│   ├── kallisto_outputs/          # Kallisto quantification results (.tsv)
│   ├── TXNAME                     # Transcript ID mapping file
│   └── SYMBOL                     # Gene symbol mapping file
│
├── snrnaseq/
│   ├── Con_filtered_feature_bc_matrix/   # Control 10X data
│   ├── ctrl_feature_bc_matrix/           # Control 10X data (dataset 2)
│   ├── cac_filtered_feature_bc_matrix/   # Cachexia 10X data
│   └── ccx_feature_bc_matrix/            # Cachexia 10X data (dataset 2)
│
├── proteomics/
│   ├── Mouse_serum_proteomics.xlsx       # Mouse mass spec results
│   ├── proteomics_per_sample.xlsx        # Per-sample abundance
│   └── Human_GC_proteomics.xlsx          # Human gastric cancer serum data
│
└── external/
    └── lr_network_mouse_21122021.rds     # NicheNet L-R network
```

## Downloading Public Data

### Bulk RNA-seq (GEO)

```bash
# Download raw FASTQ files from GEO
# GSE65936, GSE123310, GSE138464, GSE142455

# Process with Kallisto
kallisto quant -i transcriptome.idx -o output sample_R1.fastq sample_R2.fastq
```

### snRNA-seq (GEO/Zenodo)

```bash
# GSE272085 - Download from GEO
# Zenodo 11090497 - Download from Zenodo

# Files should be in 10X Genomics format:
# - matrix.mtx.gz
# - features.tsv.gz  
# - barcodes.tsv.gz
```

### NicheNet Resources

```r
# Download from Zenodo
# https://zenodo.org/record/7074291

# Or use nichenetr package
library(nichenetr)
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
```

## Data Not Included

The following data are not publicly available but can be requested:
- Mouse serum proteomics (Mass spectrometry)
- Human gastric cancer serum proteomics (TMT-labeled)

Contact the corresponding authors for data access.

## File Format Details

### Kallisto Output (.tsv)
```
target_id    length    eff_length    est_counts    tpm
ENSMUST...    1500      1200          150.0         12.5
```

### 10X Genomics Format
Standard CellRanger output format with:
- Sparse matrix (matrix.mtx.gz)
- Gene features (features.tsv.gz)
- Cell barcodes (barcodes.tsv.gz)

### Proteomics Excel Format
- Gene column: Gene symbols
- Sample columns: Abundance values (log2 or ratio)
- Statistics: proDA_logFC, proDA_padj, proDA_p
