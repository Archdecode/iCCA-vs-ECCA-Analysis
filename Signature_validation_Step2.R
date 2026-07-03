# Gap 2, Step 4 — External validation of iCCA prognostic signature
# Cohort: GSE107943 (independent iCCA cohort)

library(tidyverse)
library(survival)
library(survminer)

# =====================================================
# EDIT THESE PATHS TO MATCH YOUR SETUP
# =====================================================

# Input: expression matrix from GSE107943 (gene symbols x samples)
expr_file <- "C:/Users/kaavy/Downloads/GSE107943_expr_symbols.csv"

# Input: sample metadata with stage_group, OS.time, OS columns
meta_file <- "C:/Users/kaavy/Downloads/GSE107943_meta.csv"

# Input: your TCGA iCCA signature genes from Gap 2 Step 1
sig_genes_file <- "C:/Users/kaavy/OneDrive/Documents/OMICS_PBL/Prognostic_signature_results/prognostic_signature_genes_iCCA.csv"

# Output: folder where all results/plots from THIS script will be saved
out_dir <- "C:/Users/kaavy/OneDrive/Documents/OMICS_PBL/Signature_validation"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# Load external cohort data
# =====================================================
geo_expr_final <- read.csv(expr_file, row.names = 1, check.names = FALSE)
geo_meta       <- read.csv(meta_file)
top5_iCCA      <- read.csv(sig_genes_file)

signature_genes <- top5_iCCA$gene_symbol

genes_in_geo <- intersect(signature_genes, rownames(geo_expr_final))
cat("Signature genes found in GSE107943:", length(genes_in_geo), "/", length(signature_genes), "\n")

if (length(genes_in_geo) < 2) {
  stop("Too few signature genes present in external cohort — cannot build composite score.")
}

# =====================================================
# Build composite z-score signature
# =====================================================
build_signature_score <- function(expr_mat, gene_list) {
  keep_genes <- intersect(gene_list, rownames(expr_mat))
  sub_expr <- expr_mat[keep_genes, , drop = FALSE]
  sub_expr <- apply(sub_expr, 2, as.numeric)
  rownames(sub_expr) <- keep_genes
  z_scores <- t(scale(t(sub_expr)))
  colMeans(z_scores, na.rm = TRUE)
}

geo_score <- build_signature_score(geo_expr_final, genes_in_geo)

geo_signature_df <- tibble(
  sample_id = names(geo_score),
  sig_score = as.numeric(geo_score)
)

# =====================================================
# Merge with survival + stage metadata
# =====================================================
geo_surv_df <- geo_signature_df %>%
  left_join(geo_meta %>% select(sample_id, stage_group, OS.time, OS), by = "sample_id") %>%
  filter(!is.na(OS.time), !is.na(OS), !is.na(sig_score))

cat("Samples with valid survival + score:", nrow(geo_surv_df), "\n")

# =====================================================
# Median split within external cohort
# =====================================================
geo_cutoff <- median(geo_surv_df$sig_score, na.rm = TRUE)

geo_surv_df <- geo_surv_df %>%
  mutate(
    risk_group = ifelse(sig_score >= geo_cutoff, "High", "Low"),
    risk_group = factor(risk_group, levels = c("Low", "High"))
  )

write.csv(geo_surv_df,
          file.path(out_dir, "external_validation_survival_input.csv"),
          row.names = FALSE)

# =====================================================
# KM curve — external validation
# =====================================================
km_fit_geo <- survfit(Surv(OS.time, OS) ~ risk_group, data = geo_surv_df)

km_plot_geo <- ggsurvplot(
  km_fit_geo,
  data = geo_surv_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#1B7837", "#D6604D"),
  title = "External validation: GSE107943 (iCCA signature)",
  subtitle = paste0(length(genes_in_geo), "-gene signature, median-split, n = ", nrow(geo_surv_df)),
  xlab = "Overall survival (months)",
  ylab = "Survival probability",
  legend.title = "Risk group",
  legend.labs = c("Low", "High")
)

ggsave(
  filename = file.path(out_dir, "KM_external_GSE107943.png"),
  plot = km_plot_geo$plot,
  width = 7, height = 6, dpi = 300
)

ggsave(
  filename = file.path(out_dir, "KM_external_GSE107943_with_risktable.png"),
  plot = arrange_ggsurvplots(list(km_plot_geo), print = FALSE, ncol = 1, nrow = 1),
  width = 8, height = 8, dpi = 300
)

# =====================================================
# Log-rank summary
# =====================================================
logrank_geo <- survdiff(Surv(OS.time, OS) ~ risk_group, data = geo_surv_df)
logrank_p_geo <- 1 - pchisq(logrank_geo$chisq, df = length(logrank_geo$n) - 1)

geo_summary <- tibble(
  n_total = nrow(geo_surv_df),
  n_low = sum(geo_surv_df$risk_group == "Low"),
  n_high = sum(geo_surv_df$risk_group == "High"),
  median_cutoff = geo_cutoff,
  n_genes_used = length(genes_in_geo),
  genes_used = paste(genes_in_geo, collapse = ", "),
  logrank_chisq = unname(logrank_geo$chisq),
  logrank_pvalue = logrank_p_geo
)

write.csv(geo_summary,
          file.path(out_dir, "KM_external_GSE107943_summary.csv"),
          row.names = FALSE)

cat("\nExternal validation complete.\n")
cat("Genes used:", paste(genes_in_geo, collapse = ", "), "\n")
cat("n =", nrow(geo_surv_df), "| Log-rank p =", logrank_p_geo, "\n")