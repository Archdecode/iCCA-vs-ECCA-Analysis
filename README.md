#iCCA-vs-ECCA-Analysis
Stage-specific transcriptional programs in intrahepatic (iCCA) vs extrahepatic (eCCA) cholangiocarcinoma.

## Overview

This project identifies shared and divergent stage-specific transcriptional programs between iCCA and eCCA using TCGA (CHOL) and GEO cohorts. We work at the gene-expression level (RNA-seq counts) and map results to pathways (MSigDB Hallmark) to define transcriptional programs, comparing early vs late stage within each subtype before contrasting between subtypes.

## Pipeline
The analysis is split into four sequential modules.
| Step | Script | Description |
|------|--------|-------------|
| 1 | `Master_Sample.ipynb` | Builds the unified sample metadata table (`sample_meta_master`) across TCGA and GEO — sample IDs, subtype, stage, tumour/normal status, dataset |
| 2 | `QC AND NORMALIZATION.R` | Imports raw TCGA STAR counts and GEO expression matrices, aligns to metadata, filters and normalizes via DESeq2, builds `dds`/`vsd` objects |
| 3 | `Stage_specificDE` | Fits DESeq2 interaction model (`~ subtype + stage_group + subtype:stage_group`), extracts stage contrasts within iCCA and eCCA, classifies genes as shared or subtype-specific |
| 4 | `GSEA_and_Survival.R` | Runs fgsea on ranked gene lists per subtype, classifies pathways as shared/divergent, generates figures, links pathway activity (ssGSEA) to overall survival via Cox regression |

## Data
Raw data, intermediate `.rds`/`.csv` objects, and figures are stored in a shared Google Drive (not committed to this repo due to file size). This repo contains only the analysis code.

Required inputs (not included here):
- TCGA STAR count matrix + clinical/CDR survival data
- GEO expression matrix + sample metadata
- `sample_meta_master.csv` (output of step 1)

## Requirements

R ≥ 4.3, with the following packages:

```r
# CRAN
install.packages(c("tidyverse", "data.table", "janitor", "msigdbr",
                    "ggrepel", "pheatmap", "RColorBrewer",
                    "survminer", "readxl"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", "fgsea", "GSVA"))
```

**Key finding:** iCCA and eCCA show almost entirely divergent stage-associated transcriptional programs, with E2F target activity linked to survival in iCCA and interferon-alpha response activity linked to survival in eCCA.

## Reproducibility

Each script sets `set.seed(1234)` and ends by writing `sessionInfo()` to its results folder for version tracking.
