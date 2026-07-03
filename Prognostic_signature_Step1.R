# Gap 2 — iCCA-only prognostic signature + KM curve
# Inputs: tcga_vsd.rds, res_iCCA_stage.csv, gene_symbol_map_master.csv, tcga_master_labels.csv
# Outputs: saved to working directory via setwd()

library(tidyverse)
library(survival)
library(survminer)
library(tibble)

select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate
arrange <- dplyr::arrange

# -----------------------------
# INPUT FILE PATHS — EDIT THESE
# -----------------------------
tcga_vsd_file    <- "C:/Users/kaavy/Downloads/tcga_vsd.rds"
res_icca_file    <- "C:/Users/kaavy/Downloads/res_iCCA_stage.csv"
symbol_map_file  <- "C:/Users/kaavy/Downloads/gene_symbol_map_master.csv"
tcga_labels_file <- "C:/Users/kaavy/Downloads/tcga_master_labels.csv"

# -----------------------------
# WORKING DIRECTORY — outputs save here
# -----------------------------
setwd("C:/Users/kaavy/OneDrive/Documents/OMICS_PBL/Prognostic_signature_results")

# -----------------------------
# STEP 1 — Load inputs
# -----------------------------
library(SummarizedExperiment)

expr_vsd_sym <- readRDS(tcga_vsd_file)
res_iCCA_df  <- read.csv(res_icca_file)
symbol_map   <- read.csv(symbol_map_file)
tcga_master  <- read.csv(tcga_labels_file)

# -----------------------------
# STEP 1b — Extract matrix + convert ENSEMBL rownames to gene symbols
# -----------------------------
expr_matrix <- assay(expr_vsd_sym)

symbol_map_clean <- symbol_map %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_id, .keep_all = TRUE)

gene_lookup <- symbol_map_clean$gene_symbol
names(gene_lookup) <- symbol_map_clean$gene_id

matched_symbols <- gene_lookup[rownames(expr_matrix)]
keep_rows <- !is.na(matched_symbols)

expr_matrix_sym <- expr_matrix[keep_rows, , drop = FALSE]
rownames(expr_matrix_sym) <- matched_symbols[keep_rows]

# collapse duplicate gene symbols by averaging (if multiple ENSEMBL IDs map to the same symbol)
expr_matrix_sym <- rowsum(expr_matrix_sym, group = rownames(expr_matrix_sym)) /
  as.vector(table(rownames(expr_matrix_sym))[rownames(rowsum(expr_matrix_sym, group = rownames(expr_matrix_sym)))])

expr_vsd_sym <- expr_matrix_sym

cat("Converted matrix dimensions:", dim(expr_vsd_sym), "\n")
cat("Sample of converted rownames:", paste(head(rownames(expr_vsd_sym)), collapse = ", "), "\n")

intersect(top5_iCCA$gene_symbol, rownames(expr_vsd_sym))
# -----------------------------
# STEP 2 — Clean symbol map
# -----------------------------
symbol_map <- symbol_map %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_id, .keep_all = TRUE)

# -----------------------------
# STEP 3 — Select top iCCA genes
# -----------------------------
res_iCCA_sym <- res_iCCA_df %>%
  left_join(symbol_map, by = "gene_id") %>%
  filter(!is.na(gene_symbol))

if ("sig" %in% colnames(res_iCCA_sym)) {
  top5_iCCA <- res_iCCA_sym %>%
    filter(sig == TRUE, !is.na(padj)) %>%
    arrange(padj) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    head(5)
} else {
  top5_iCCA <- res_iCCA_sym %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    head(5)
}

write.csv(top5_iCCA, "prognostic_signature_genes_iCCA.csv", row.names = FALSE)

# -----------------------------
# STEP 4 — Build iCCA signature
# -----------------------------
build_signature_score <- function(expr_mat, gene_list) {
  keep_genes <- intersect(gene_list, rownames(expr_mat))
  
  if (length(keep_genes) == 0) {
    stop("None of the selected genes were found in expr_vsd_sym rownames.")
  }
  
  sub_expr <- expr_mat[keep_genes, , drop = FALSE]
  
  if (is.data.frame(sub_expr)) {
    sub_expr <- as.matrix(sub_expr)
  }
  
  sub_expr <- apply(sub_expr, 2, as.numeric)
  rownames(sub_expr) <- keep_genes
  
  z_scores <- t(scale(t(sub_expr)))
  score <- colMeans(z_scores, na.rm = TRUE)
  
  list(score = score, genes_used = keep_genes)
}

sig_obj <- build_signature_score(expr_vsd_sym, top5_iCCA$gene_symbol)

signature_df <- tibble(
  sample_id = names(sig_obj$score),
  iCCA_sig_score = as.numeric(sig_obj$score)
)

write.csv(signature_df, "signature_scores_iCCA_all_samples.csv", row.names = FALSE)

# -----------------------------
# STEP 5 — Merge with survival metadata
# -----------------------------
surv_sig_df <- signature_df %>%
  left_join(
    tcga_master %>%
      select(SAMPLE_ID, PATIENT_ID, subtype, OS_days, OS_event),
    by = c("sample_id" = "SAMPLE_ID")
  ) %>%
  filter(
    subtype == "iCCA",
    !is.na(OS_days),
    !is.na(OS_event),
    !is.na(iCCA_sig_score)
  ) %>%
  mutate(
    OS.time = OS_days,
    OS = OS_event
  )

write.csv(surv_sig_df, "survival_input_iCCA.csv", row.names = FALSE)

# -----------------------------
# STEP 6 — Median split
# -----------------------------
median_cut <- median(surv_sig_df$iCCA_sig_score, na.rm = TRUE)

surv_iCCA <- surv_sig_df %>%
  mutate(
    risk_group = ifelse(iCCA_sig_score >= median_cut, "High", "Low"),
    risk_group = factor(risk_group, levels = c("Low", "High"))
  )

write.csv(surv_iCCA, "survival_input_iCCA_mediansplit.csv", row.names = FALSE)

# -----------------------------
# STEP 7 — Kaplan–Meier
# -----------------------------
km_fit_iCCA <- survfit(Surv(OS.time, OS) ~ risk_group, data = surv_iCCA)

km_plot <- ggsurvplot(
  km_fit_iCCA,
  data = surv_iCCA,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#1B7837", "#D6604D"),
  title = "iCCA 5-gene prognostic signature",
  xlab = "Overall survival (days)",
  ylab = "Survival probability",
  legend.title = "Risk group",
  legend.labs = c("Low", "High")
)

png("KM_iCCA_signature_with_risktable.png", width = 8, height = 8, units = "in", res = 300)
print(km_plot)
dev.off()
# -----------------------------
# STEP 8 — Log-rank summary
# -----------------------------
logrank_test <- survdiff(Surv(OS.time, OS) ~ risk_group, data = surv_iCCA)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)

summary_df <- tibble(
  n_total_iCCA = nrow(surv_iCCA),
  n_low = sum(surv_iCCA$risk_group == "Low"),
  n_high = sum(surv_iCCA$risk_group == "High"),
  median_cutoff = median_cut,
  n_signature_genes_selected = nrow(top5_iCCA),
  n_signature_genes_used = length(sig_obj$genes_used),
  genes_used = paste(sig_obj$genes_used, collapse = "; "),
  logrank_chisq = unname(logrank_test$chisq),
  logrank_pvalue = logrank_p
)

write.csv(summary_df, "KM_iCCA_summary.csv", row.names = FALSE)

cat("Done. Outputs saved to:", getwd(), "\n")
cat("Top selected genes:", paste(top5_iCCA$gene_symbol, collapse = ", "), "\n")
cat("Genes present in expr_vsd_sym:", paste(sig_obj$genes_used, collapse = ", "), "\n")
cat("iCCA samples used for survival:", nrow(surv_iCCA), "\n")
cat("Median cutoff:", median_cut, "\n")
cat("Log-rank p-value:", logrank_p, "\n")