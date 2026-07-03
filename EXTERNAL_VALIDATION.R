# External Validation — GSE107943 (iCCA)
# Member 4 | validate top 20 iCCA stage genes in independent cohort

library(GEOquery)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

base_dir  <- "C:/Users/archi/Desktop/RVCE/EL/sem_6/OMICS_PBL"
val_dir   <- file.path(base_dir, "results/validation")
rpkm_file <- file.path(val_dir, "GSE107943/GSE107943_RPKM.txt.gz")
de_file   <- file.path(base_dir, "res_iCCA_stage.csv")

dir.create(val_dir, recursive = TRUE, showWarnings = FALSE)


# Metadata

gse       <- getGEO("GSE107943", GSEMatrix = TRUE)
geo_pheno <- pData(gse[[1]])

geo_meta <- geo_pheno %>%
  rownames_to_column("sample_id") %>%
  mutate(
    group = case_when(
      grepl("Tumor",    `tissue:ch1`, ignore.case = TRUE) ~ "Tumour",
      grepl("adjacent", `tissue:ch1`, ignore.case = TRUE) ~ "Normal",
      TRUE ~ NA_character_
    ),
    stage_raw   = `stageajcc:ch1`,
    stage_group = case_when(
      stage_raw %in% c("I", "II")           ~ "Early",
      stage_raw %in% c("III", "IVA", "IVB") ~ "Late",
      TRUE ~ NA_character_
    ),
    OS.time = as.numeric(`survival(mo):ch1`),
    OS      = as.numeric(`death:ch1`)
  )


# Download expression file if not already present

if (!file.exists(rpkm_file)) {
  getGEOSuppFiles("GSE107943", baseDir = val_dir)
}


# Load and clean expression matrix

geo_expr <- read.table(
  rpkm_file,
  header      = TRUE,
  sep         = "\t",
  row.names   = NULL,
  check.names = FALSE
)

geo_expr_clean <- geo_expr %>%
  filter(!is.na(Ensenble), Ensenble != "") %>%
  column_to_rownames("Ensenble") %>%
  dplyr::select(-No, -Chr, -Genesymbol, -Start, -Stop)


# Map ENSEMBL IDs to gene symbols

row_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = rownames(geo_expr_clean),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)
row_map <- row_map[!is.na(row_map$SYMBOL), ]
row_map <- row_map[!duplicated(row_map$ENSEMBL), ]

keep         <- rownames(geo_expr_clean) %in% row_map$ENSEMBL
geo_expr_sym <- geo_expr_clean[keep, ]
new_symbols  <- row_map$SYMBOL[match(rownames(geo_expr_sym), row_map$ENSEMBL)]

geo_expr_sym <- geo_expr_sym[!duplicated(new_symbols) & !is.na(new_symbols), ]
new_symbols  <- new_symbols[!duplicated(new_symbols) & !is.na(new_symbols)]
rownames(geo_expr_sym) <- new_symbols


# Bridge internal sample IDs to GSM accessions

expr_cols      <- colnames(geo_expr_sym)
expr_cols      <- expr_cols[expr_cols != "CodingLength"]
expr_ids_clean <- gsub("^s_", "",    expr_cols)
expr_ids_clean <- gsub("_RPKM$", "", expr_ids_clean)

bridge <- data.frame(
  expr_col = expr_cols,
  clean_id = expr_ids_clean,
  stringsAsFactors = FALSE
) %>%
  left_join(
    geo_pheno %>%
      rownames_to_column("sample_id") %>%
      dplyr::select(sample_id, description),
    by = c("clean_id" = "description")
  )

geo_expr_final <- geo_expr_sym %>%
  dplyr::select(-CodingLength)

colnames(geo_expr_final) <- bridge$sample_id[
  match(colnames(geo_expr_final), bridge$expr_col)
]

geo_meta_aligned <- geo_meta %>%
  filter(sample_id %in% colnames(geo_expr_final))


# Top 20 iCCA stage genes from TCGA

res_iCCA_df <- read.csv(de_file)

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res_iCCA_df$gene_id),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)
symbol_map <- symbol_map[!is.na(symbol_map$SYMBOL), ]
symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), ]

top20_iCCA <- res_iCCA_df %>%
  left_join(symbol_map, by = c("gene_id" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL), !is.na(padj)) %>%
  arrange(padj) %>%
  head(20) %>%
  dplyr::select(gene_id, gene_symbol = SYMBOL, log2FC, padj)

overlap_genes <- intersect(top20_iCCA$gene_symbol, rownames(geo_expr_final))
cat("Genes found in GSE107943:", length(overlap_genes), "/ 20\n")


# Late vs Early concordance

early_samples <- geo_meta_aligned$sample_id[
  !is.na(geo_meta_aligned$stage_group) & geo_meta_aligned$stage_group == "Early"
]
late_samples <- geo_meta_aligned$sample_id[
  !is.na(geo_meta_aligned$stage_group) & geo_meta_aligned$stage_group == "Late"
]
early_samples <- early_samples[!is.na(early_samples)]
late_samples  <- late_samples[!is.na(late_samples)]

stage_concordance <- map_dfr(overlap_genes, function(g) {
  early_vals <- as.numeric(geo_expr_final[g, early_samples])
  late_vals  <- as.numeric(geo_expr_final[g, late_samples])
  if (length(early_vals) < 2 | length(late_vals) < 2) return(NULL)
  test <- t.test(late_vals, early_vals)
  tibble(
    gene                = g,
    geo_mean_early      = mean(early_vals, na.rm = TRUE),
    geo_mean_late       = mean(late_vals,  na.rm = TRUE),
    geo_stage_pval      = test$p.value,
    geo_stage_direction = ifelse(mean(late_vals) > mean(early_vals), "up", "down")
  )
})

tcga_direction <- top20_iCCA %>%
  mutate(tcga_direction = ifelse(log2FC > 0, "up", "down")) %>%
  dplyr::select(gene_symbol, log2FC, tcga_direction)

final_stage_concordance <- stage_concordance %>%
  left_join(tcga_direction, by = c("gene" = "gene_symbol")) %>%
  mutate(concordant = geo_stage_direction == tcga_direction)

cat("Concordant:", sum(final_stage_concordance$concordant), "/",
    nrow(final_stage_concordance), "\n")


# Save outputs

write.csv(final_stage_concordance,
          file.path(val_dir, "concordance_late_vs_early.csv"), row.names = FALSE)
write.csv(top20_iCCA,
          file.path(val_dir, "top20_iCCA_genes.csv"), row.names = FALSE)
write.csv(as.data.frame(geo_expr_final),
          file.path(val_dir, "GSE107943_expr_symbols.csv"), row.names = TRUE)
write.csv(geo_meta_aligned,
          file.path(val_dir, "GSE107943_meta.csv"), row.names = FALSE)


# Concordance figure

library(ggplot2)

p_conc <- final_stage_concordance %>%
  mutate(
    gene        = reorder(gene, geo_mean_late - geo_mean_early),
    concordance = ifelse(concordant, "Concordant", "Discordant"),
    sig_label   = ifelse(geo_stage_pval < 0.1, "*", "")
  ) %>%
  ggplot(aes(x = gene, y = geo_mean_late - geo_mean_early, fill = concordance)) +
  geom_col() +
  geom_text(aes(label = sig_label),
            hjust = -0.3, size = 5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Concordant" = "#1B7837", "Discordant" = "#D6604D")) +
  coord_flip() +
  labs(
    title    = "External validation: GSE107943 (iCCA, n=30)",
    subtitle = "Late vs Early stage | * = p < 0.1 | Green = concordant with TCGA direction",
    x        = "Gene",
    y        = "GEO: mean(Late) - mean(Early)",
    fill     = "Concordance with TCGA"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(val_dir, "concordance_stage_barplot.png"),
       p_conc, width = 8, height = 7, dpi = 300)
cat("Figure saved.\n")

p_conc <- final_stage_concordance %>%
  mutate(
    gene        = reorder(gene, geo_mean_late - geo_mean_early),
    concordance = ifelse(concordant, "Concordant", "Discordant"),
    diff_capped = pmax(pmin(geo_mean_late - geo_mean_early, 5), -5),
    sig_label   = ifelse(geo_stage_pval < 0.1, "*", "")
  ) %>%
  ggplot(aes(x = gene, y = diff_capped, fill = concordance)) +
  geom_col() +
  geom_text(aes(label = sig_label), hjust = -0.3, size = 5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Concordant" = "#1B7837", "Discordant" = "#D6604D")) +
  coord_flip() +
  labs(
    title    = "External validation: GSE107943 (iCCA, n=30)",
    subtitle = "Late vs Early stage | * = p < 0.1 | Bars capped at \u00b15 for visibility",
    x        = "Gene",
    y        = "GEO: mean(Late) - mean(Early) [capped at \u00b15]",
    fill     = "Concordance with TCGA"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(val_dir, "concordance_stage_barplot2.png"),
       p_conc, width = 8, height = 7, dpi = 300)
