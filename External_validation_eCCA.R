# SECTION 0 — Install and load packages

cran_pkgs <- c("tidyverse", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer")
new_cran  <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_cran) > 0) install.packages(new_cran)

bioc_pkgs <- c("GEOquery", "limma", "org.Hs.eg.db", "AnnotationDbi", "Biobase")
new_bioc  <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(new_bioc)
}

library(GEOquery)
library(limma)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Biobase)
library(ggplot2)
library(ggrepel)
library(pheatmap)

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
rename    <- dplyr::rename

set.seed(1234)

base_dir  <- "C:/Users/archi/Desktop/RVCE/EL/sem_6/OMICS_PBL"
out_dir   <- file.path(base_dir, "results/eCCA_GSE132305")
de_file   <- file.path(base_dir, "res_eCCA_stage.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Downloading GSE132305...\n")
gse       <- getGEO("GSE132305", GSEMatrix = TRUE)
geo_expr  <- exprs(gse[[1]])        # probes x samples matrix
geo_pheno <- pData(gse[[1]])        # GEO sample annotation

cat("Expression matrix dimensions:", nrow(geo_expr), "probes x",
    ncol(geo_expr), "samples\n")
cat("Available annotation columns:\n")
print(colnames(geo_pheno))
cat("\nTissue labels:\n")
print(table(geo_pheno$`tissue:ch1`))

cat("\nBuilding metadata table...\n")
geo_meta <- geo_pheno %>%
  rownames_to_column("sample_id") %>%
  mutate(
    group   = case_when(
      grepl("extrahepatic cholangiocarcinoma", `tissue:ch1`,
            ignore.case = TRUE) ~ "Tumour",
      grepl("non-tumor|non-tumour|bile duct|normal",
            `tissue:ch1`, ignore.case = TRUE) ~ "Normal",
      TRUE ~ NA_character_
    ),
    subtype  = "eCCA",
    dataset  = "GSE132305",
    platform = platform_id,
    # no stage available in this dataset — record explicitly
    stage_group = NA_character_,
    OS.time     = NA_real_,
    OS          = NA_real_
  ) %>%
  dplyr::select(sample_id, group, subtype, dataset, platform,
                stage_group, OS.time, OS,
                tissue_label = `tissue:ch1`)

cat("Tumour samples:", sum(geo_meta$group == "Tumour", na.rm = TRUE), "\n")
cat("Normal samples:", sum(geo_meta$group == "Normal", na.rm = TRUE), "\n")

write.csv(geo_meta,
          file.path(out_dir, "ecca_gse132305_meta.csv"),
          row.names = FALSE)
cat("Metadata saved.\n")

cat("\nMapping probes to gene symbols...\n")

# Get the platform annotation table
platform_data <- gse[[1]]@featureData@data
cat("Platform annotation columns:", paste(colnames(platform_data), collapse = ", "), "\n")
cat("First 3 rows of platform data:\n")
print(head(platform_data, 3))

symbol_col <- "Gene Symbol"   

probe_map <- platform_data %>%
  rownames_to_column("probe_id") %>%
  dplyr::select(probe_id, gene_symbol = all_of(symbol_col)) %>%
  filter(!is.na(gene_symbol), gene_symbol != "", gene_symbol != "---")

cat("Probes with valid gene symbols:", nrow(probe_map), "\n")
cat("Unique genes:", length(unique(probe_map$gene_symbol)), "\n")


cat("\nChecking expression value range...\n")
cat("Min:", round(min(geo_expr, na.rm = TRUE), 2),
    "Max:", round(max(geo_expr, na.rm = TRUE), 2),
    "Mean:", round(mean(geo_expr, na.rm = TRUE), 2), "\n")

# If max > 100, data is not log-transformed — log2 transform it
if (max(geo_expr, na.rm = TRUE) > 100) {
  cat("Values appear to be raw intensity — applying log2 transformation\n")
  geo_expr <- log2(geo_expr + 1)
} else {
  cat("Values appear to be log2-transformed already — no transformation needed\n")
}

# Collapse probes to gene symbols by median
cat("\nCollapsing probes to gene symbols...\n")
geo_expr_sym <- geo_expr %>%
  as.data.frame() %>%
  rownames_to_column("probe_id") %>%
  left_join(probe_map, by = "probe_id") %>%
  filter(!is.na(gene_symbol)) %>%
  dplyr::select(-probe_id) %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  column_to_rownames("gene_symbol")

cat("Genes after collapse:", nrow(geo_expr_sym), "\n")
cat("Samples:", ncol(geo_expr_sym), "\n")

write.csv(geo_expr_sym,
          file.path(out_dir, "ecca_gse132305_expr_norm.csv"),
          row.names = TRUE)
cat("Normalised expression matrix saved.\n")

cat("\nRunning limma Tumour vs Normal DE...\n")

# Align expression matrix columns to metadata
common_samples <- intersect(colnames(geo_expr_sym), geo_meta$sample_id)
cat("Samples in both expression and metadata:", length(common_samples), "\n")

expr_aligned <- geo_expr_sym[, common_samples]
meta_aligned <- geo_meta %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

stopifnot(identical(colnames(expr_aligned), meta_aligned$sample_id))

# Build design matrix
group_factor <- factor(meta_aligned$group, levels = c("Normal", "Tumour"))
design       <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Normal", "Tumour")

# Fit limma model
fit    <- lmFit(as.matrix(expr_aligned), design)
contr  <- makeContrasts(Tumour - Normal, levels = design)
fit2   <- contrasts.fit(fit, contr)
fit2   <- eBayes(fit2)

# Extract results
de_results <- topTable(fit2, number = Inf, adjust.method = "BH") %>%
  rownames_to_column("gene_symbol") %>%
  as_tibble() %>%
  rename(log2FC = logFC, padj = adj.P.Val, pval = P.Value) %>%
  mutate(
    direction = case_when(
      padj < 0.05 & log2FC > 1  ~ "UP",
      padj < 0.05 & log2FC < -1 ~ "DOWN",
      TRUE                       ~ "NS"
    )
  ) %>%
  arrange(padj)

cat("\n--- Tumour vs Normal DE summary ---\n")
cat("Upregulated in tumour (padj<0.05, log2FC>1):",
    sum(de_results$direction == "UP"), "\n")
cat("Downregulated in tumour (padj<0.05, log2FC<-1):",
    sum(de_results$direction == "DOWN"), "\n")

write.csv(de_results,
          file.path(out_dir, "ecca_gse132305_DE_results.csv"),
          row.names = FALSE)
cat("DE results saved.\n")

cat("\nRunning eCCA gene concordance analysis...\n")

# Load your team's TCGA eCCA DE results
res_eCCA_df <- read.csv(de_file)

# Map ENSEMBL IDs to gene symbols
symbol_map_tcga <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res_eCCA_df$gene_id),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)
symbol_map_tcga <- symbol_map_tcga[!is.na(symbol_map_tcga$SYMBOL), ]
symbol_map_tcga <- symbol_map_tcga[!duplicated(symbol_map_tcga$ENSEMBL), ]

# Top 20 eCCA stage genes by adjusted p-value
top20_eCCA <- res_eCCA_df %>%
  left_join(symbol_map_tcga, by = c("gene_id" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL), !is.na(padj)) %>%
  arrange(padj) %>%
  head(20) %>%
  dplyr::select(gene_id, gene_symbol = SYMBOL, log2FC, padj)

cat("Top 20 eCCA stage genes:\n")
print(top20_eCCA$gene_symbol)

# Check overlap with GSE132305
overlap_genes <- intersect(top20_eCCA$gene_symbol, de_results$gene_symbol)
cat("\nGenes found in GSE132305:", length(overlap_genes), "/ 20\n")

# Check directional concordance
# TCGA direction: log2FC from Late vs Early (positive = up in Late stage)
# GSE132305 direction: log2FC from Tumour vs Normal (positive = up in tumour)
# These aren't identical comparisons but same-direction is still informative
tcga_dir <- top20_eCCA %>%
  filter(gene_symbol %in% overlap_genes) %>%
  mutate(tcga_direction = ifelse(log2FC > 0, "up", "down")) %>%
  dplyr::select(gene_symbol, log2FC_tcga = log2FC, tcga_direction)

geo_dir <- de_results %>%
  filter(gene_symbol %in% overlap_genes) %>%
  dplyr::select(gene_symbol, log2FC_geo = log2FC, pval,
                geo_direction = direction) %>%
  mutate(geo_direction = ifelse(log2FC_geo > 0, "up", "down"))

concordance <- tcga_dir %>%
  left_join(geo_dir, by = "gene_symbol") %>%
  mutate(concordant = tcga_direction == geo_direction)

cat("\n--- Concordance summary ---\n")
cat("Concordant:", sum(concordance$concordant), "/",
    nrow(concordance), "\n")
print(concordance %>%
        dplyr::select(gene_symbol, tcga_direction,
                      geo_direction, log2FC_geo, pval, concordant))

write.csv(concordance,
          file.path(out_dir, "ecca_concordance_results.csv"),
          row.names = FALSE)

p_conc <- concordance %>%
  mutate(
    gene        = reorder(gene_symbol, log2FC_geo),
    concordance = ifelse(concordant, "Concordant", "Discordant"),
    sig_label   = ifelse(pval < 0.1, "*", "")
  ) %>%
  ggplot(aes(x = gene, y = log2FC_geo, fill = concordance)) +
  geom_col() +
  geom_text(aes(label = sig_label), hjust = -0.3, size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Concordant" = "#1B7837",
                               "Discordant" = "#D6604D")) +
  coord_flip() +
  labs(
    title    = "GSE132305 validation: eCCA Tumour vs Normal",
    subtitle = "Direction compared to TCGA eCCA stage genes | * = p < 0.1",
    x        = "Gene",
    y        = "GSE132305 log2FC (Tumour vs Normal)",
    fill     = "Concordance with TCGA"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "ecca_concordance_barplot.png"),
       p_conc, width = 8, height = 7, dpi = 300)
cat("Figure saved.\n")

processed_object <- list(
  expr     = geo_expr_sym,       # gene symbols x samples, normalised
  meta     = geo_meta,           # sample metadata in team format
  de       = de_results,         # tumour vs normal DE results
  concordance = concordance      # concordance with TCGA eCCA stage genes
)

saveRDS(processed_object,
        file.path(out_dir, "ecca_gse132305_processed.rds"))

cat("\n--- GSE132305 processing complete ---\n")
cat("All outputs saved to:", out_dir, "\n")
cat("Files:\n")
print(list.files(out_dir))

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_eCCA_processing.txt"))