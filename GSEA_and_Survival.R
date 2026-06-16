# Fix the namespace conflict permanently for this session
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise

cran_pkgs <- c(
  "tidyverse", "data.table", "janitor",
  "msigdbr", "ggrepel", "pheatmap",
  "RColorBrewer", "survminer", "readxl"
)
new_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_cran) > 0) {
  message("Installing CRAN packages: ", paste(new_cran, collapse = ", "))
  install.packages(new_cran)
}
# Bioconductor packages
bioc_pkgs <- c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", "fgsea", "GSVA")
new_bioc  <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(new_bioc, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(new_bioc)
}

cat("All packages installed.\n")

library(DESeq2)
library(tidyverse)
library(data.table)
library(janitor)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)
library(fgsea)
library(GSVA)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(survival)
library(survminer)
library(readxl)

set.seed(1234)
options(stringsAsFactors = FALSE)

results_dir <- "C:/Users/archi/Desktop/RVCE/EL/sem_6/OMICS PBL"

# Input paths — all files are directly in OMICS PBL
de_dir   <- results_dir
proc_dir <- results_dir
cdr_path <- file.path(results_dir, "TCGA-CDR-SupplementalTableS1.xlsx")

# Output paths — R will create these fresh for your outputs
gsea_dir   <- file.path(results_dir, "results/gsea")
fig_dir    <- file.path(results_dir, "results/figures")
enrich_dir <- file.path(results_dir, "results/figures/gsea_enrichment_plots")

# Create output folders
dir.create(gsea_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

cat("Directories ready.\n")

cat("Loading Member 3 DE results...\n")
res_merged  <- read.csv(file.path(de_dir, "res_merged_classified.csv"))
res_iCCA_df <- read.csv(file.path(de_dir, "res_iCCA_stage.csv"))
res_eCCA_df <- read.csv(file.path(de_dir, "res_eCCA_stage.csv"))

cat("Loading Member 2 VSD object...\n")
vsd       <- readRDS(file.path(proc_dir, "tcga_vsd.rds"))
tcga_meta <- readRDS(file.path(proc_dir, "tcga_meta.rds"))

cat("Rows in res_merged:", nrow(res_merged), "\n")
cat("Gene categories from Member 3:\n")
print(table(res_merged$category))

cat("\nMapping ENSEMBL IDs to gene symbols...\n")

all_ensg <- unique(res_merged$gene_id)

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = all_ensg,
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
) %>%
  as_tibble() %>%
  rename(gene_id = ENSEMBL, gene_symbol = SYMBOL) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_id, .keep_all = TRUE)

cat("Mapped", nrow(symbol_map), "of", length(all_ensg), "genes to symbols.\n")

res_merged_sym <- res_merged %>%
  left_join(symbol_map, by = "gene_id") %>%
  filter(!is.na(gene_symbol))

res_iCCA_sym <- res_iCCA_df %>%
  left_join(symbol_map, by = "gene_id") %>%
  filter(!is.na(gene_symbol))

res_eCCA_sym <- res_eCCA_df %>%
  left_join(symbol_map, by = "gene_id") %>%
  filter(!is.na(gene_symbol))

cat("\nFetching MSigDB Hallmark gene sets...\n")

msig_h <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  group_by(gs_name) %>%
  summarise(gene_symbols = list(gene_symbol)) %>%
  deframe()

cat("Hallmark gene sets loaded:", length(msig_h), "pathways.\n")


cat("\nBuilding rank vectors...\n")

rank_iCCA <- res_merged_sym %>%
  filter(!is.na(log2FC_iCCA)) %>%
  arrange(desc(log2FC_iCCA)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  { setNames(.$log2FC_iCCA, .$gene_symbol) }

rank_eCCA <- res_merged_sym %>%
  filter(!is.na(log2FC_eCCA)) %>%
  arrange(desc(log2FC_eCCA)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  { setNames(.$log2FC_eCCA, .$gene_symbol) }

cat("iCCA rank vector length:", length(rank_iCCA), "\n")
cat("eCCA rank vector length:", length(rank_eCCA), "\n")

cat("\nRunning fgsea for iCCA...\n")
fg_iCCA <- fgsea(
  pathways    = msig_h,
  stats       = rank_iCCA,
  nPermSimple = 10000,
  minSize     = 15,
  maxSize     = 500
)

cat("Running fgsea for eCCA...\n")
fg_eCCA <- fgsea(
  pathways    = msig_h,
  stats       = rank_eCCA,
  nPermSimple = 10000,
  minSize     = 15,
  maxSize     = 500
)

fg_iCCA_df <- as_tibble(fg_iCCA) %>% mutate(subtype = "iCCA")
fg_eCCA_df <- as_tibble(fg_eCCA) %>% mutate(subtype = "eCCA")

cat("\niCCA significant pathways (padj<0.05):",
    sum(fg_iCCA_df$padj < 0.05, na.rm = TRUE), "\n")
cat("eCCA significant pathways (padj<0.05):",
    sum(fg_eCCA_df$padj < 0.05, na.rm = TRUE), "\n")

write.csv(fg_iCCA_df %>% select(-leadingEdge),
          file.path(gsea_dir, "fg_iCCA.csv"), row.names = FALSE)
write.csv(fg_eCCA_df %>% select(-leadingEdge),
          file.path(gsea_dir, "fg_eCCA.csv"), row.names = FALSE)
cat("Individual GSEA tables saved.\n")

cat("\nClassifying pathways...\n")

fg_wide <- fg_iCCA_df %>%
  select(pathway, nes_iCCA = NES, padj_iCCA = padj) %>%
  full_join(
    fg_eCCA_df %>% select(pathway, nes_eCCA = NES, padj_eCCA = padj),
    by = "pathway"
  ) %>%
  mutate(category = case_when(
    padj_iCCA < 0.05 & padj_eCCA < 0.05 & sign(nes_iCCA) == sign(nes_eCCA) ~ "Shared_concordant",
    padj_iCCA < 0.05 & padj_eCCA < 0.05 & sign(nes_iCCA) != sign(nes_eCCA) ~ "Shared_opposite",
    padj_iCCA < 0.05 & (is.na(padj_eCCA) | padj_eCCA >= 0.05)              ~ "iCCA_specific",
    padj_eCCA < 0.05 & (is.na(padj_iCCA) | padj_iCCA >= 0.05)              ~ "eCCA_specific",
    TRUE                                                                      ~ "Not_significant"
  ))

cat("\n--- Pathway classification ---\n")
print(table(fg_wide$category))

write.csv(fg_wide, file.path(gsea_dir, "fg_merged_classified.csv"),
          row.names = FALSE)
cat("Classified pathway table saved.\n")

category_colors <- c(
  "Shared_concordant" = "#1B7837",
  "Shared_opposite"   = "#762A83",
  "iCCA_specific"     = "#2166AC",
  "eCCA_specific"     = "#D6604D",
  "Not_significant"   = "grey80"
)


# ---- 7A: NES scatter plot ---------------------------------------------------

cat("\nGenerating NES scatter plot...\n")

label_pathways <- fg_wide %>%
  filter(category != "Not_significant") %>%
  group_by(category) %>%
  arrange(desc(abs(nes_iCCA) + abs(nes_eCCA))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(pathway)

p_nes_scatter <- fg_wide %>%
  mutate(label = ifelse(pathway %in% label_pathways,
                        gsub("HALLMARK_", "", pathway), NA)) %>%
  ggplot(aes(nes_iCCA, nes_eCCA, color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = label), size = 3,
                  max.overlaps = 25, na.rm = TRUE,
                  segment.color = "grey60") +
  scale_color_manual(values = category_colors) +
  labs(
    title = "Pathway NES: iCCA vs eCCA (Late vs Early stage)",
    x     = "NES — iCCA (Late vs Early)",
    y     = "NES — eCCA (Late vs Early)",
    color = "Category"
  ) +
  theme_classic(base_size = 13)

ggsave(file.path(fig_dir, "pathway_nes_scatter.png"),
       p_nes_scatter, width = 8, height = 7, dpi = 300)
cat("NES scatter saved.\n")

---- 7B: Pathway heatmap ---------------------------------------------------
cat("Generating pathway heatmap...\n")

key_pathways <- fg_wide %>%
  filter(category %in% c("Shared_concordant", "Shared_opposite",
                         "iCCA_specific", "eCCA_specific")) %>%
  mutate(abs_sum = abs(nes_iCCA) + abs(nes_eCCA)) %>%
  arrange(desc(abs_sum)) %>%
  head(40) %>%
  pull(pathway)

if (length(key_pathways) >= 2) {
  mat_nes <- fg_wide %>%
    filter(pathway %in% key_pathways) %>%
    select(pathway, nes_iCCA, nes_eCCA) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  
  rownames(mat_nes) <- gsub("HALLMARK_", "", rownames(mat_nes))
  
  pathway_ann <- fg_wide %>%
    filter(pathway %in% key_pathways) %>%
    mutate(pathway_label = gsub("HALLMARK_", "", pathway)) %>%
    select(pathway_label, Category = category) %>%
    column_to_rownames("pathway_label")
  
  ann_colors <- list(
    Category = category_colors[names(category_colors) != "Not_significant"]
  )
  
  pheatmap(mat_nes,
           annotation_row    = pathway_ann,
           annotation_colors = ann_colors,
           cluster_cols      = FALSE,
           cluster_rows      = TRUE,
           color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(50),
           breaks            = seq(-3, 3, length.out = 51),
           fontsize_row      = 8,
           fontsize_col      = 11,
           main              = "Stage-specific transcriptional programs (NES)",
           labels_col        = c("iCCA", "eCCA"),
           filename          = file.path(fig_dir, "pathway_heatmap.png"),
           width = 10, height = 12)
  cat("Pathway heatmap saved.\n")
} else {
  cat("Not enough significant pathways for heatmap.\n")
}

# ---- 7C: Individual enrichment plots ---------------------------------------

cat("Generating individual enrichment plots...\n")

top_shared    <- fg_wide %>% filter(category == "Shared_concordant") %>%
  arrange(desc(abs(nes_iCCA))) %>% head(3) %>% pull(pathway)
top_iCCA_spec <- fg_wide %>% filter(category == "iCCA_specific") %>%
  arrange(desc(abs(nes_iCCA))) %>% head(2) %>% pull(pathway)
top_eCCA_spec <- fg_wide %>% filter(category == "eCCA_specific") %>%
  arrange(desc(abs(nes_eCCA))) %>% head(2) %>% pull(pathway)

plot_pathways <- unique(c(top_shared, top_iCCA_spec, top_eCCA_spec))

for (pw in plot_pathways) {
  if (!pw %in% names(msig_h)) next
  pw_label <- gsub("HALLMARK_", "", pw)
  
  png(file.path(enrich_dir, paste0("enrich_iCCA_", pw_label, ".png")),
      width = 800, height = 500, res = 120)
  print(plotEnrichment(msig_h[[pw]], rank_iCCA) +
          labs(title = paste("iCCA —", pw_label)) + theme_classic())
  dev.off()
  
  png(file.path(enrich_dir, paste0("enrich_eCCA_", pw_label, ".png")),
      width = 800, height = 500, res = 120)
  print(plotEnrichment(msig_h[[pw]], rank_eCCA) +
          labs(title = paste("eCCA —", pw_label)) + theme_classic())
  dev.off()
}
cat("Enrichment plots saved.\n")

# All significant pathways
sig_pathways_table <- fg_wide %>%
  filter(category != "Not_significant") %>%
  dplyr::select(pathway, nes_iCCA, nes_eCCA, padj_iCCA, padj_eCCA, category) %>%
  arrange(category)

print(sig_pathways_table)
write.csv(sig_pathways_table,
          file.path(gsea_dir, "significant_pathways.csv"),
          row.names = FALSE)
cat("Saved to:", file.path(gsea_dir, "significant_pathways.csv"), "\n")

# Shared opposite pathway
shared_opposite_table <- fg_wide %>%
  filter(category == "Shared_opposite") %>%
  dplyr::select(pathway, nes_iCCA, nes_eCCA, padj_iCCA, padj_eCCA)

print(shared_opposite_table)
write.csv(shared_opposite_table,
          file.path(gsea_dir, "shared_opposite_pathway.csv"),
          row.names = FALSE)
cat("Saved to:", file.path(gsea_dir, "shared_opposite_pathway.csv"), "\n")

cat("\nRunning ssGSEA on VSD expression matrix...\n")

expr_vsd  <- assay(vsd)
ensg_rows <- rownames(expr_vsd)

row_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = ensg_rows,
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
) %>%
  as_tibble() %>%
  rename(gene_id = ENSEMBL, gene_symbol = SYMBOL) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_id, .keep_all = TRUE)

keep_idx       <- ensg_rows %in% row_map$gene_id
expr_vsd_sym   <- expr_vsd[keep_idx, ]
rownames(expr_vsd_sym) <- row_map$gene_symbol[
  match(rownames(expr_vsd_sym), row_map$gene_id)
]
expr_vsd_sym <- expr_vsd_sym[!duplicated(rownames(expr_vsd_sym)), ]

cat("Expression matrix for ssGSEA:", nrow(expr_vsd_sym), "genes x",
    ncol(expr_vsd_sym), "samples\n")

# GSVA >= 1.50 API; falls back to legacy call if older version installed
gsva_scores <- tryCatch({
  gsva_param <- ssgseaParam(exprData = expr_vsd_sym, geneSets = msig_h, minSize = 5)
  gsva(gsva_param)
}, error = function(e) {
  message("ssgseaParam not found, using legacy gsva() call.")
  gsva(expr_vsd_sym, msig_h, method = "ssgsea", min.sz = 5)
})

gsva_scores_df <- t(gsva_scores) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

cat("ssGSEA complete:", nrow(gsva_scores_df), "samples x",
    ncol(gsva_scores_df) - 1, "pathways\n")

write.csv(gsva_scores_df,
          file.path(gsea_dir, "gsva_scores.csv"), row.names = FALSE)
cat("ssGSEA scores saved.\n")

cat("\nLoading CDR clinical data...\n")

cdr      <- as.data.frame(read_excel(cdr_path))
cdr_chol <- cdr %>%
  filter(type == "CHOL") %>%
  select(bcr_patient_barcode, ajcc_pathologic_tumor_stage,
         histological_type, OS, OS.time) %>%
  mutate(
    subtype = case_when(
      grepl("intrahepatic",               histological_type, ignore.case = TRUE) ~ "iCCA",
      grepl("hilar|perihilar|distal",     histological_type, ignore.case = TRUE) ~ "eCCA",
      TRUE ~ NA_character_
    ),
    OS      = as.numeric(OS),
    OS.time = as.numeric(OS.time)
  )

# Match 15-char VSD barcodes to 12-char patient IDs
gsva_scores_df <- gsva_scores_df %>%
  mutate(patient_id = substr(sample_id, 1, 12))

surv_df <- gsva_scores_df %>%
  left_join(cdr_chol, by = c("patient_id" = "bcr_patient_barcode")) %>%
  filter(!is.na(OS), !is.na(OS.time), !is.na(subtype))

cat("Samples with survival data:", nrow(surv_df), "\n")


# ---- 9A: Cox regression for every significant pathway ----------------------

run_cox_pathways <- function(df, subtype_label, pathways_to_test) {
  df_sub <- df %>% filter(subtype == subtype_label)
  map_dfr(pathways_to_test, function(pw) {
    if (!pw %in% colnames(df_sub)) return(NULL)
    cox_data <- df_sub %>%
      select(OS.time, OS, score = all_of(pw)) %>%
      filter(!is.na(score), !is.na(OS.time), !is.na(OS))
    if (nrow(cox_data) < 5) return(NULL)
    fit <- tryCatch(
      coxph(Surv(OS.time, OS) ~ score, data = cox_data),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    tibble(
      pathway  = pw,
      subtype  = subtype_label,
      HR       = s$coefficients[1, "exp(coef)"],
      CI_lower = s$conf.int[1, "lower .95"],
      CI_upper = s$conf.int[1, "upper .95"],
      p_value  = s$coefficients[1, "Pr(>|z|)"]
    )
  })
}

sig_pathways <- fg_wide %>%
  filter(category != "Not_significant") %>%
  pull(pathway) %>% unique()

cat("Running Cox models for", length(sig_pathways), "pathways per subtype...\n")

cox_iCCA <- run_cox_pathways(surv_df, "iCCA", sig_pathways)
cox_eCCA <- run_cox_pathways(surv_df, "eCCA", sig_pathways)

cox_all <- bind_rows(cox_iCCA, cox_eCCA) %>%
  group_by(subtype) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(padj)

cat("\nTop 10 pathway-survival associations:\n")
print(head(cox_all %>% select(pathway, subtype, HR, p_value, padj), 10))

write.csv(cox_all,
          file.path(gsea_dir, "cox_pathway_survival.csv"), row.names = FALSE)


# ---- 9B: KM plots for top pathway per subtype ------------------------------

km_plot <- function(df, subtype_label, pathway_name, out_path) {
  df_sub <- df %>%
    dplyr::filter(subtype == subtype_label, !is.na(.data[[pathway_name]])) %>%
    dplyr::mutate(activity_group = ifelse(
      .data[[pathway_name]] >= median(.data[[pathway_name]]),
      "High", "Low"
    ))
  if (nrow(df_sub) < 4) {
    cat("Not enough samples for KM plot:", subtype_label, pathway_name, "\n")
    return(invisible(NULL))
  }
  km_fit <- survfit(Surv(OS.time, OS) ~ activity_group, data = df_sub)
  p <- ggsurvplot(km_fit, data = df_sub,
                  pval        = TRUE,
                  conf.int    = TRUE,
                  risk.table  = TRUE,
                  palette     = c("#D6604D", "#2166AC"),
                  title       = paste(subtype_label, "—",
                                      gsub("HALLMARK_", "", pathway_name)),
                  xlab        = "Time (days)",
                  ylab        = "Overall survival",
                  legend.labs = c("High activity", "Low activity"))
  
  # correct way to save ggsurvplot — use png() directly instead of ggsave
  png(out_path, width = 1000, height = 800, res = 120)
  print(p)
  dev.off()
  
  cat("KM plot saved:", out_path, "\n")
}

top_iCCA_surv <- if (nrow(cox_iCCA) > 0)
  cox_iCCA %>% arrange(p_value) %>% slice(1) %>% pull(pathway) else character(0)
top_eCCA_surv <- if (nrow(cox_eCCA) > 0)
  cox_eCCA %>% arrange(p_value) %>% slice(1) %>% pull(pathway) else character(0)

if (length(top_iCCA_surv) > 0)
  km_plot(surv_df, "iCCA", top_iCCA_surv,
          file.path(fig_dir, paste0("survival_km_iCCA_",
                                    gsub("HALLMARK_", "", top_iCCA_surv), ".png")))

if (length(top_eCCA_surv) > 0)
  km_plot(surv_df, "eCCA", top_eCCA_surv,
          file.path(fig_dir, paste0("survival_km_eCCA_",
                                    gsub("HALLMARK_", "", top_eCCA_surv), ".png")))
writeLines(capture.output(sessionInfo()),
           file.path(gsea_dir, "sessionInfo_member4.txt"))

cat("\n--- Member 4 complete ---\n")
cat("GSEA outputs:\n");     print(list.files(gsea_dir))
cat("Figures:\n");          print(list.files(fig_dir))
cat("Enrichment plots:\n"); print(list.files(enrich_dir))
