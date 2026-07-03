# =============================================================================
# Gene Symbol Conversion + Hub Gene Analysis
# iCCA vs eCCA Stage Transcriptomics Pipeline
# =============================================================================
# Description:
# - Converts all Ensembl IDs across DE result files to HGNC gene symbols
# - Builds a master symbol mapping table for use by all team members
# - Extracts top 50 subtype-specific genes for STRING network analysis
# - Identifies top 5 hub genes per subtype from STRING node degree output
# - Regenerates volcano, scatter, and heatmap figures with gene symbols
#
# Inputs:
#   - res_iCCA_stage.csv         : DESeq2 results for iCCA (Late vs Early)
#   - res_eCCA_stage.csv         : DESeq2 results for eCCA (Late vs Early)
#   - res_merged_classified.csv  : merged DE results with subtype categories
#   - string_node_degrees.tsv    : node degree file exported from STRING DB
#   - results/processed/tcga_vsd.rds   : variance-stabilized DESeq2 object
#   - results/processed/tcga_meta.rds  : aligned TCGA sample metadata
#
# Outputs (saved to Downloads/):
#   - gene_symbol_map_master.csv   : Ensembl -> gene symbol mapping table
#   - top50_iCCA_genes.txt         : top 50 iCCA-specific genes for STRING
#   - top50_eCCA_genes.txt         : top 50 eCCA-specific genes for STRING
#   - hub_genes_iCCA.csv           : top 5 hub genes for iCCA (by degree)
#   - hub_genes_eCCA.csv           : top 5 hub genes for eCCA (by degree)
#   - figures/volcano_iCCA.pdf
#   - figures/volcano_eCCA.pdf
#   - figures/scatter_iCCA_vs_eCCA.pdf
#   - figures/heatmap_top_DE_genes.pdf
#
# Note: STRING node degree files must be downloaded manually from
#       string-db.org after uploading each gene list, then re-run
#       the hub gene section for iCCA and eCCA separately.
# =============================================================================

# --- Libraries ----------------------------------------------------------------

library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(DESeq2)

# --- File paths ---------------------------------------------------------------
# Update these paths if your files are stored elsewhere

INPUT_DIR   <- "C:/Users/navni/Downloads"
OUTPUT_DIR  <- "C:/Users/navni/Downloads"
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load DE result files --------------------------------------------------

res_iCCA   <- read.csv(file.path(INPUT_DIR, "res_iCCA_stage.csv"))
res_eCCA   <- read.csv(file.path(INPUT_DIR, "res_eCCA_stage.csv"))
res_merged <- read.csv(file.path(INPUT_DIR, "res_merged_classified.csv"))

# --- 2. Build master gene symbol mapping table --------------------------------

all_ensg <- unique(c(res_iCCA$gene_id, res_eCCA$gene_id, res_merged$gene_id))

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = all_ensg,
  columns  = c("ENSEMBL", "SYMBOL"),
  keytype  = "ENSEMBL"
) %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  rename(gene_id = ENSEMBL, gene_symbol = SYMBOL)

# Diagnostics
cat("Total genes mapped:", nrow(symbol_map), "\n")
cat("Genes with no symbol (NA):", sum(is.na(symbol_map$gene_symbol)), "\n")

# Save — share this file with all team members immediately
write.csv(symbol_map,
          file.path(OUTPUT_DIR, "gene_symbol_map_master.csv"),
          row.names = FALSE)

# --- 3. Extract top 50 subtype-specific genes for STRING ----------------------

# iCCA-specific genes ranked by |log2FC|
top50_iCCA <- res_merged %>%
  left_join(symbol_map, by = "gene_id") %>%
  filter(category == "iCCA_specific", !is.na(gene_symbol)) %>%
  arrange(desc(abs(log2FC_iCCA))) %>%
  head(50) %>%
  pull(gene_symbol)

# eCCA-specific genes ranked by |log2FC|
top50_eCCA <- res_merged %>%
  left_join(symbol_map, by = "gene_id") %>%
  filter(category == "eCCA_specific", !is.na(gene_symbol)) %>%
  arrange(desc(abs(log2FC_eCCA))) %>%
  head(50) %>%
  pull(gene_symbol)

cat("iCCA top 50 genes extracted:", length(top50_iCCA), "\n")
cat("eCCA top 50 genes extracted:", length(top50_eCCA), "\n")

writeLines(top50_iCCA, file.path(OUTPUT_DIR, "top50_iCCA_genes.txt"))
writeLines(top50_eCCA, file.path(OUTPUT_DIR, "top50_eCCA_genes.txt"))

# NOTE: Upload each .txt file to string-db.org -> Multiple Proteins ->
# Homo sapiens -> Continue -> Exports -> download "protein node degrees" TSV
# Run the hub gene section below separately for iCCA and eCCA

# --- 4. Hub gene identification from STRING node degrees ----------------------
# Run this section twice:
#   First:  after uploading iCCA list to STRING and downloading node degrees
#   Second: after uploading eCCA list to STRING and downloading node degrees

# iCCA hub genes
node_degrees_iCCA <- read.table(
  file.path(INPUT_DIR, "string_node_degrees.tsv"),
  header = TRUE, sep = "\t"
)
colnames(node_degrees_iCCA) <- c("gene", "accession", "degree")
node_degrees_iCCA <- node_degrees_iCCA %>% arrange(desc(degree))

cat("Top 10 iCCA hub genes:\n")
print(head(node_degrees_iCCA, 10))

hub_genes_iCCA <- node_degrees_iCCA %>% head(5)
write.csv(hub_genes_iCCA,
          file.path(OUTPUT_DIR, "hub_genes_iCCA.csv"),
          row.names = FALSE)

# eCCA hub genes
# Re-download string_node_degrees.tsv from STRING after uploading eCCA list
node_degrees_eCCA <- read.table(
  file.path(INPUT_DIR, "string_node_degrees.tsv"),
  header = TRUE, sep = "\t"
)
colnames(node_degrees_eCCA) <- c("gene", "accession", "degree")
node_degrees_eCCA <- node_degrees_eCCA %>% arrange(desc(degree))

cat("Top 10 eCCA hub genes:\n")
print(head(node_degrees_eCCA, 10))

hub_genes_eCCA <- node_degrees_eCCA %>% head(5)
write.csv(hub_genes_eCCA,
          file.path(OUTPUT_DIR, "hub_genes_eCCA.csv"),
          row.names = FALSE)

# --- 5. Add gene symbols to result files for figures -------------------------

res_iCCA_sym <- res_iCCA %>%
  left_join(symbol_map, by = "gene_id") %>%
  mutate(label = ifelse(sig == TRUE & !is.na(gene_symbol), gene_symbol, NA))

res_eCCA_sym <- res_eCCA %>%
  left_join(symbol_map, by = "gene_id") %>%
  mutate(label = ifelse(sig == TRUE & !is.na(gene_symbol), gene_symbol, NA))

res_merged_sym <- res_merged %>%
  left_join(symbol_map, by = "gene_id") %>%
  filter(!is.na(gene_symbol)) %>%
  mutate(label = ifelse(
    category %in% c("iCCA_specific", "eCCA_specific") &
      (abs(log2FC_iCCA) > quantile(abs(log2FC_iCCA), 0.95, na.rm = TRUE) |
       abs(log2FC_eCCA) > quantile(abs(log2FC_eCCA), 0.95, na.rm = TRUE)),
    gene_symbol, NA
  ))

# --- 6. Volcano plot — iCCA ---------------------------------------------------

ggplot(res_iCCA_sym, aes(x = log2FC, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  na.rm = TRUE) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  labs(title = "iCCA: Late vs Early Stage",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  theme_classic()

ggsave(file.path(FIGURES_DIR, "volcano_iCCA.pdf"), width = 8, height = 6)

# --- 7. Volcano plot — eCCA ---------------------------------------------------

ggplot(res_eCCA_sym, aes(x = log2FC, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  na.rm = TRUE) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  labs(title = "eCCA: Late vs Early Stage",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  theme_classic()

ggsave(file.path(FIGURES_DIR, "volcano_eCCA.pdf"), width = 8, height = 6)

# --- 8. Scatter plot — iCCA vs eCCA log2FC ------------------------------------

ggplot(res_merged_sym, aes(x = log2FC_iCCA, y = log2FC_eCCA, color = category)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "iCCA vs eCCA log2FC comparison",
       x = "log2FC iCCA (Late vs Early)",
       y = "log2FC eCCA (Late vs Early)") +
  theme_classic()

ggsave(file.path(FIGURES_DIR, "scatter_iCCA_vs_eCCA.pdf"), width = 8, height = 7)

# --- 9. Heatmap — top DE genes with gene symbols ------------------------------

vsd      <- readRDS("C:/Users/navni/Downloads/tcga_vsd.rds")
tcga_meta <- readRDS("C:/Users/navni/Downloads/tcga_meta.rds")

# Top 30 significant genes per subtype
top_genes <- bind_rows(
  res_iCCA_sym %>% filter(sig == TRUE) %>% arrange(padj) %>% head(30),
  res_eCCA_sym %>% filter(sig == TRUE) %>% arrange(padj) %>% head(30)
) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_id, .keep_all = TRUE)

# Extract and label VST matrix
vsd_mat <- assay(vsd)
rownames(vsd_mat) <- sub("\\..*", "", rownames(vsd_mat))
heat_mat <- vsd_mat[rownames(vsd_mat) %in% top_genes$gene_id, ]
rownames(heat_mat) <- top_genes$gene_symbol[match(rownames(heat_mat), top_genes$gene_id)]
heat_mat <- heat_mat[!is.na(rownames(heat_mat)), ]

# Z-score per gene
heat_scaled <- t(scale(t(heat_mat)))

# Column annotation
ann_col <- data.frame(
  label = tcga_meta$label,
  row.names = tcga_meta$sample_id
)

pheatmap(heat_scaled,
         annotation_col  = ann_col,
         show_colnames   = FALSE,
         fontsize_row    = 8,
         main            = "Top DE genes — iCCA & eCCA",
         filename        = file.path(FIGURES_DIR, "heatmap_top_DE_genes.pdf"),
         width           = 10,
         height          = 12)

cat("Gap 4 complete. All outputs saved to:", OUTPUT_DIR, "\n")
