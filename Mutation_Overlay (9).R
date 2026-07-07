# ══════════════════════════════════════════════════════════
# GAP 5 — Somatic Mutation Overlay (Driver Gene Oncoprint)
# TCGA-CHOL (n=34)
# ══════════════════════════════════════════════════════════

library(DESeq2)
library(tidyverse)
library(TCGAbiolinks)
library(maftools)
library(ComplexHeatmap)

dl  <- "C:/Users/Subhashis/Downloads"
out <- "C:/Users/Subhashis/Desktop/OMICS_PBL"
setwd(dl)

dds_sub <- readRDS(file.path(out, "results/de/tcga_dds_interaction.rds"))
complete_samples <- as.data.frame(colData(dds_sub))
complete_samples$sample_id  <- rownames(complete_samples)
complete_samples$patient_id <- substr(complete_samples$sample_id, 1, 12)

# ── Download and prepare MAF ────────────────────────────────
maf_dir <- file.path(dl, "GDC_MAF")
dir.create(maf_dir, showWarnings = FALSE)

maf_query <- GDCquery(
  project = "TCGA-CHOL",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query, directory = maf_dir)
maf_data <- GDCprepare(maf_query, directory = maf_dir)

maf_data$patient_id <- substr(maf_data$Tumor_Sample_Barcode, 1, 12)
maf_sub <- maf_data %>% filter(patient_id %in% complete_samples$patient_id)

# ── Driver gene panel ────────────────────────────────────────
driver_genes <- c("TP53", "KRAS", "IDH1", "IDH2", "FGFR2",
                  "ARID1A", "BAP1", "SMAD4", "PBRM1", "CDKN2A")
maf_driver <- maf_sub %>% filter(Hugo_Symbol %in% driver_genes)

# ── Binary mutation matrix (gene x patient) ─────────────────
mut_matrix <- maf_driver %>%
  distinct(Hugo_Symbol, patient_id, Variant_Classification) %>%
  mutate(present = Variant_Classification) %>%
  select(Hugo_Symbol, patient_id, present) %>%
  pivot_wider(names_from = patient_id, values_from = present,
              values_fn = function(x) paste(unique(x), collapse = ";"), values_fill = "") %>%
  column_to_rownames("Hugo_Symbol") %>% as.matrix()

missing_genes <- setdiff(driver_genes, rownames(mut_matrix))
if (length(missing_genes) > 0) {
  empty_rows <- matrix("", nrow = length(missing_genes), ncol = ncol(mut_matrix),
                       dimnames = list(missing_genes, colnames(mut_matrix)))
  mut_matrix <- rbind(mut_matrix, empty_rows)
}
mut_matrix <- mut_matrix[driver_genes, , drop = FALSE]

missing_patients <- setdiff(complete_samples$patient_id, colnames(mut_matrix))
if (length(missing_patients) > 0) {
  na_cols <- matrix("", nrow = nrow(mut_matrix), ncol = length(missing_patients),
                    dimnames = list(rownames(mut_matrix), missing_patients))
  mut_matrix <- cbind(mut_matrix, na_cols)
}
mut_matrix <- mut_matrix[, complete_samples$patient_id, drop = FALSE]

# ── Oncoprint ────────────────────────────────────────────────
ann_onco <- complete_samples %>% select(patient_id, subtype, stage_group) %>%
  distinct(patient_id, .keep_all = TRUE)
ann_onco <- ann_onco[match(colnames(mut_matrix), ann_onco$patient_id), ]

top_ann <- HeatmapAnnotation(
  Subtype = ann_onco$subtype, Stage = ann_onco$stage_group,
  col = list(Subtype = c("iCCA" = "#2166AC", "eCCA" = "#D6604D"),
             Stage   = c("Early" = "#92C5DE", "Late" = "#B2182B"))
)
alter_fun <- list(
  background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "grey90", col = NA)),
  Mutated    = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#08519C", col = NA))
)
mut_matrix_simple <- ifelse(mut_matrix == "", "", "Mutated")

png(file.path(out, "results/figures/oncoprint_CCA.png"), width = 12, height = 6, units = "in", res = 300)
oncoPrint(mut_matrix_simple, alter_fun = alter_fun, col = c(Mutated = "#08519C"),
          top_annotation = top_ann, column_title = "Driver gene mutations — TCGA-CHOL cohort",
          row_names_side = "left", pct_side = "right", show_column_names = FALSE)
dev.off()

# ── Fisher's exact test — mutation vs stage (iCCA only) ─────
fisher_results <- data.frame()
for (gene in driver_genes) {
  gene_mut <- ifelse(mut_matrix[gene, ] != "", 1, 0)
  df_test <- data.frame(patient_id = colnames(mut_matrix), mutated = gene_mut) %>%
    left_join(ann_onco, by = "patient_id") %>% filter(subtype == "iCCA")
  if (length(unique(df_test$stage_group)) == 2 && length(unique(df_test$mutated)) > 1) {
    tab <- table(df_test$mutated, df_test$stage_group)
    ft  <- tryCatch(fisher.test(tab), error = function(e) NULL)
    fisher_results <- rbind(fisher_results, data.frame(
      gene = gene, p_value = if (!is.null(ft)) ft$p.value else NA,
      n_mutated_Early = sum(df_test$mutated == 1 & df_test$stage_group == "Early"),
      n_mutated_Late  = sum(df_test$mutated == 1 & df_test$stage_group == "Late"),
      n_total_Early   = sum(df_test$stage_group == "Early"),
      n_total_Late    = sum(df_test$stage_group == "Late")
    ))
  }
}
fisher_results <- fisher_results %>% arrange(p_value)
write.csv(fisher_results, file.path(out, "results/de/mutation_stage_fisher.csv"), row.names = FALSE)
write.csv(as.data.frame(mut_matrix), file.path(out, "results/de/mutation_matrix.csv"))

# ── Mutation frequency summary + bar chart ──────────────────
mut_freq_summary <- data.frame(
  gene = driver_genes,
  n_mutated = sapply(driver_genes, function(g) sum(mut_matrix[g, ] != "")),
  pct_mutated = round(100 * sapply(driver_genes, function(g) sum(mut_matrix[g, ] != "")) / ncol(mut_matrix), 1)
) %>% arrange(desc(n_mutated))
write.csv(mut_freq_summary, file.path(out, "results/de/mutation_frequency_summary.csv"), row.names = FALSE)

p_freq <- ggplot(mut_freq_summary, aes(x = reorder(gene, -n_mutated), y = pct_mutated)) +
  geom_col(fill = "#2166AC") +
  geom_text(aes(label = paste0(n_mutated, "/", ncol(mut_matrix))), vjust = -0.5, size = 3.5) +
  labs(title = "Driver gene mutation frequency — TCGA-CHOL cohort", x = NULL, y = "% mutated") +
  theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out, "results/figures/mutation_frequency_bar.png"), p_freq, width = 7, height = 5, dpi = 300)

cat("GAP 5 COMPLETE.\n")
