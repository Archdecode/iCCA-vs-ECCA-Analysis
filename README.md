# iCCA-vs-ECCA-Analysis

Stage-specific transcriptional programs in intrahepatic (iCCA) vs extrahepatic (eCCA) cholangiocarcinoma using TCGA-CHOL and GEO cohorts.

## Overview

This project identifies shared and divergent stage-specific transcriptional programs between iCCA and eCCA. Analysis is performed at the gene-expression level (TCGA STAR counts) and mapped to pathways (MSigDB Hallmark), comparing early-stage (I/II) vs late-stage (III/IV) within each subtype before contrasting between subtypes. An independent GEO cohort (GSE107943, n=30 iCCA) was used for external validation. The pipeline addresses all six gaps required for Q2 publication, including immune deconvolution, prognostic signature validation, hub gene analysis, mutation overlay, and power analysis.

## Scripts

| # | Script | What it does |
|---|--------|-------------|
| 1 | `Master_Sample.ipynb` | Builds the unified sample metadata table (`sample_meta_master.csv`) across TCGA and GEO with sample IDs, subtype, stage, tumour/normal label, dataset, and STAR count availability flags |
| 2 | `QC AND NORMALIZATION.R` | Loads 43 TCGA STAR count samples, strips ENSEMBL version suffixes, builds a DESeq2 `dds` object, filters low-count genes, runs `vst()` normalization, processes GEO microarray expression, saves `tcga_dds.rds`, `tcga_vsd.rds`, `tcga_meta.rds`, `geo_expr_filtered.rds` and CSV equivalents |
| 3 | `Stage_specificDE` | Pulls TCGA-CDR clinical data, annotates iCCA/eCCA from histology, collapses TNM staging into Early (I/II) and Late (III/IV), fits a DESeq2 interaction model (`~ subtype + stage_group + subtype:stage_group`), extracts stage contrasts within each subtype and the interaction term, classifies all genes as `Shared_concordant`, `Shared_opposite`, `iCCA_specific`, `eCCA_specific`, or `Not_significant` |
| 4 | `GSEA_and_Survival.R` | Maps ENSEMBL IDs to gene symbols, runs fgsea on ranked gene lists (log2FC, Late vs Early) per subtype against MSigDB Hallmark gene sets, classifies pathways using the same shared/divergent logic as the gene-level analysis, generates NES scatter plot, pathway heatmap and enrichment plots, runs ssGSEA to score per-sample pathway activity, links scores to overall survival via Cox regression and Kaplan-Meier |
| 5 | `EXTERNAL_VALIDATION.R` | Downloads GSE107943 from GEO, aligns RPKM expression matrix to gene symbols via ENSEMBL mapping, bridges internal sample IDs to GSM accessions, checks concordance of the top 20 TCGA iCCA stage genes in the independent cohort (Late vs Early), generates concordance barplot |
| 6a | `Prognostic_signature.R` | Selects top 5 iCCA-specific DE genes by adjusted p-value, builds a composite z-score risk signature, applies median split in TCGA-CHOL iCCA samples, tests overall survival separation with Kaplan-Meier and log-rank test |
| 6b | `Signature_validation.R` | Validates the iCCA prognostic signature in the independent GSE107943 cohort, applies median split on signature score, generates KM curve and log-rank summary |
| 7 | `immune_Deconvolution.R` | Runs quanTIseq immune deconvolution on the TCGA-CHOL VST matrix and the GSE107943 expression matrix, tests CD8+ T cell, NK cell, macrophage M1/M2, and dendritic cell fractions for stage-associated differences within each subtype using Wilcoxon tests, generates boxplots and heatmaps for both cohorts, correlates immune scores with the prognostic signature score |
| 8 | `GENE_SYMBOLS_HUBS.R` | Builds the master ENSEMBL-to-symbol mapping table shared across all scripts, extracts the top 50 iCCA-specific and eCCA-specific genes by log2FC for STRING network submission, identifies top 5 hub genes per subtype from STRING node degree output, regenerates all volcano, scatter, and heatmap figures with HGNC gene symbols |
| 9 | `Mutation_Overlay.R` | Downloads TCGA-CHOL masked somatic mutation MAF files via TCGAbiolinks, builds a binary mutation matrix for 10 CCA driver genes (TP53, KRAS, IDH1, IDH2, FGFR2, ARID1A, BAP1, SMAD4, PBRM1, CDKN2A), generates an oncoprint annotated by subtype and stage, runs Fisher's exact tests for mutation-stage association within iCCA |
| 10 | `POWER_SENSITIVITY.R` | Computes minimum detectable effect size (Cohen's d = 0.99, large effect) given n=34 at alpha=0.05, power=0.80; runs 100-iteration leave-two-out bootstrap sensitivity analysis showing 99.1% of top 10 iCCA-specific genes remain significant after random sample removal |

## Key Findings

- iCCA and eCCA show almost entirely divergent stage-specific transcriptional programs — no shared concordant pathways were identified across 50 Hallmark gene sets
- E2F targets are the top iCCA-specific stage pathway and associate with overall survival
- Interferon alpha response is the top eCCA-specific stage pathway and associates with overall survival
- One shared-opposite pathway was identified, representing the strongest point of biological divergence between subtypes
- A 5-gene iCCA prognostic signature stratifies overall survival in both TCGA-CHOL and the independent GSE107943 cohort
- 9 of 20 top iCCA stage genes showed concordant direction of expression change in GSE107943, with ABCB1 and ARSJ showing the strongest cross-cohort signal
- Immune deconvolution validated pathway-level immune findings at the cell-composition level in both cohorts
- Power analysis confirmed the study detects large effect sizes (Cohen's d ≥ 0.99) and 99.1% of top findings are stable under leave-two-out resampling

## Data

Raw data, processed objects, DE results, GSEA outputs, and figures are stored in a shared Google Drive. This repository contains analysis code only.

Required inputs to reproduce the analysis:
- TCGA-CHOL STAR count matrix and clinical CDR supplemental table
- GEO expression matrix and sample metadata
- `sample_meta_master.csv` (output of `Master_Sample.ipynb`)
- `gene_symbol_map_master.csv` (output of `GENE_SYMBOLS_HUBS.R`, required by scripts 6a, 6b, 7)
- STRING node degree TSV files (downloaded manually from string-db.org, required by `GENE_SYMBOLS_HUBS.R`)

## Run order

Scripts must be run in the numbered order above. Scripts 6a, 6b, 7, 8, 9, and 10 can be run in parallel once scripts 1-5 are complete.

## Requirements

R ≥ 4.3

```r
# CRAN
install.packages(c("tidyverse", "data.table", "janitor", "msigdbr",
                   "ggrepel", "pheatmap", "RColorBrewer", "survminer",
                   "readxl", "pwr", "maftools"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi",
                       "fgsea", "GSVA", "GEOquery", "TCGAbiolinks",
                       "ComplexHeatmap", "immunedeconv"))
```

## Reproducibility

All scripts set `set.seed(1234)`. Session info is written to the results folder at the end of each run.
