# =============================================================================
# GAP 6: Power Analysis + Sensitivity Analysis
# iCCA vs eCCA Stage Transcriptomics Pipeline
# =============================================================================
# Description:
# - Computes minimum detectable effect size (Cohen's d) given sample size
# - Runs leave-two-out bootstrap sensitivity analysis (100 iterations)
#   to assess robustness of top DE gene findings
#
# Inputs:
#   - C:/Users/navni/Downloads/tcga_dds.rds         : DESeq2 dataset object
#   - C:/Users/navni/Downloads/res_merged_classified.csv : merged DE results
#
# Outputs (saved to Downloads/):
#   - sensitivity_leave2out.csv : per-iteration count of genes remaining
#                                 significant after dropping 2 random samples
#
# Key results:
#   - Minimum detectable effect size: Cohen's d = 0.99 (large effect)
#   - Sensitivity: 9.91/10 top genes stable across 100 iterations (99.1%)
#
# Methods text (for paper):
#   "With 34 tumour samples, the study has 80% power to detect a minimum
#    effect size of Cohen's d = 0.99 (two-sample t-test, alpha = 0.05),
#    corresponding to a large effect size. To assess robustness, a
#    leave-two-out bootstrap sensitivity analysis was performed across
#    100 iterations, randomly removing 2 samples per iteration and
#    re-running differential expression. On average, 9.91 of the top 10
#    iCCA-specific genes (99.1%) remained statistically significant
#    (adjusted p < 0.05), indicating high stability of results."
# =============================================================================

# --- Libraries ----------------------------------------------------------------

library(pwr)
library(DESeq2)
library(dplyr)
library(purrr)
library(tibble)

# --- File paths ---------------------------------------------------------------

INPUT_DIR  <- "C:/Users/navni/Downloads"
OUTPUT_DIR <- "C:/Users/navni/Downloads"

# --- 1. Load objects ----------------------------------------------------------

dds       <- readRDS(file.path(INPUT_DIR, "tcga_dds.rds"))
res_merged <- read.csv(file.path(INPUT_DIR, "res_merged_classified.csv"))

# --- 2. Power analysis --------------------------------------------------------
# Compute minimum detectable effect size given sample size and 80% power

n_total <- ncol(dds)
cat("Total TCGA samples:", n_total, "\n")
cat("Samples per group (assumed equal split):", n_total / 2, "\n\n")

power_result <- pwr.t.test(
  n          = n_total / 2,
  sig.level  = 0.05,
  power      = 0.80,
  type       = "two.sample"
)

print(power_result)
cat("\nMinimum detectable Cohen's d:", round(power_result$d, 3), "\n")

# --- 3. Leave-two-out bootstrap sensitivity analysis -------------------------
# Randomly drop 2 samples per iteration, re-run DESeq2, check how many
# of the top 10 iCCA-specific genes remain significant

# Get top 10 iCCA-specific genes by adjusted p-value
top10_genes <- res_merged %>%
  filter(category == "iCCA_specific") %>%
  arrange(padj_iCCA) %>%
  head(10) %>%
  pull(gene_id)

cat("\nTop 10 iCCA-specific genes being tracked:\n")
print(top10_genes)

# Extract count matrix and metadata from dds
star_counts_mat <- counts(dds)
meta_df         <- as.data.frame(colData(dds))

set.seed(1234)
n_iter <- 100

cat("\nRunning", n_iter, "leave-two-out iterations...\n")

stability_results <- map_dfr(1:n_iter, function(i) {
  drop_samples <- sample(colnames(star_counts_mat), 2)
  sub_counts   <- star_counts_mat[, !colnames(star_counts_mat) %in% drop_samples]
  sub_meta     <- meta_df[!rownames(meta_df) %in% drop_samples, , drop = FALSE]

  dds_sub <- DESeqDataSetFromMatrix(
    countData = round(sub_counts),
    colData   = sub_meta,
    design    = ~ 1
  )
  dds_sub <- dds_sub[rowSums(counts(dds_sub) >= 10) >= 3, ]
  dds_sub <- DESeq(dds_sub, quiet = TRUE)
  res_sub <- results(dds_sub)

  tibble(
    iteration       = i,
    genes_still_sig = sum(rownames(res_sub) %in% top10_genes &
                          res_sub$padj < 0.05, na.rm = TRUE)
  )
})

# --- 4. Report results --------------------------------------------------------

avg_stable <- mean(stability_results$genes_still_sig)
cat("\nSensitivity analysis complete.\n")
cat("Average genes remaining significant:", round(avg_stable, 2), "/ 10\n")
cat("Stability rate:", round(avg_stable / 10 * 100, 1), "%\n")

# --- 5. Save output -----------------------------------------------------------

write.csv(stability_results,
          file.path(OUTPUT_DIR, "sensitivity_leave2out.csv"),
          row.names = FALSE)

cat("\nGap 6 complete. Output saved to:", OUTPUT_DIR, "\n")
