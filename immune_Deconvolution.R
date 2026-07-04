# ══════════════════════════════════════════════════════════
# GAP 3 — Immune Microenvironment Deconvolution
# TCGA-CHOL (n=34) discovery + GSE107943 (n=30) validation
# ══════════════════════════════════════════════════════════

library(DESeq2)
library(tidyverse)
library(readxl)
library(pheatmap)
library(immunedeconv)

set.seed(1234)
dl  <- "C:/Users/Subhashis/Downloads"
out <- "C:/Users/Subhashis/Desktop/OMICS_PBL"
setwd(dl)

# ── PART A: TCGA discovery cohort (n=34) ───────────────────
dds_sub  <- readRDS(file.path(out, "results/de/tcga_dds_interaction.rds"))
vsd_sub  <- vst(dds_sub, blind = FALSE)
expr_vsd <- assay(vsd_sub)

complete_samples <- as.data.frame(colData(dds_sub))
complete_samples$sample_id  <- rownames(complete_samples)
complete_samples$patient_id <- substr(complete_samples$sample_id, 1, 12)

gene_map <- read.csv(file.path(dl, "gene_symbol_map_master.csv"))
colnames(gene_map) <- c("gene_id", "gene_symbol")

expr_symbol <- as.data.frame(expr_vsd) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_map, by = "gene_id") %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  column_to_rownames("gene_symbol") %>%
  select(-gene_id) %>%
  as.matrix()

immune_scores <- deconvolute(expr_symbol, method = "quantiseq")
write.csv(immune_scores, file.path(out, "results/de/immune_deconv_scores.csv"), row.names = FALSE)

immune_long <- immune_scores %>%
  pivot_longer(-cell_type, names_to = "sample_id", values_to = "fraction") %>%
  left_join(complete_samples %>% select(sample_id, subtype, stage_group), by = "sample_id")

all_tests <- immune_long %>%
  group_by(cell_type, subtype) %>%
  filter(n_distinct(stage_group) == 2) %>%
  summarise(
    p_value    = tryCatch(wilcox.test(fraction ~ stage_group)$p.value, error = function(e) NA_real_),
    n_Early    = sum(stage_group == "Early"),
    n_Late     = sum(stage_group == "Late"),
    mean_Early = mean(fraction[stage_group == "Early"], na.rm = TRUE),
    mean_Late  = mean(fraction[stage_group == "Late"],  na.rm = TRUE),
    direction  = ifelse(mean_Late > mean_Early, "UP_in_Late", "DOWN_in_Late"),
    .groups = "drop"
  ) %>%
  arrange(p_value)
all_tests$padj_iCCA_only <- p.adjust(ifelse(all_tests$subtype == "iCCA", all_tests$p_value, NA), method = "BH")
write.csv(all_tests, file.path(out, "results/de/immune_wilcoxon_results.csv"), row.names = FALSE)

p_immune <- immune_long %>%
  filter(cell_type %in% c("T cell CD8+", "NK cell", "Macrophage M1", "Macrophage M2", "Myeloid dendritic cell")) %>%
  ggplot(aes(x = stage_group, y = fraction, fill = stage_group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  facet_grid(cell_type ~ subtype, scales = "free_y") +
  scale_fill_manual(values = c("Early" = "#2166AC", "Late" = "#D6604D")) +
  labs(title = "TCGA (n=34): Immune cell infiltration by stage and subtype",
       x = "Stage", y = "Estimated cell fraction", fill = "Stage") +
  theme_classic(base_size = 11) + theme(strip.text = element_text(size = 8), legend.position = "bottom")
ggsave(file.path(out, "results/figures/immune_boxplot.png"), p_immune, width = 8, height = 10, dpi = 300)

immune_wide <- immune_scores %>% as.data.frame() %>% remove_rownames() %>%
  column_to_rownames("cell_type") %>% as.matrix()
ann_col_immune <- complete_samples %>% select(sample_id, subtype, stage_group) %>%
  remove_rownames() %>% column_to_rownames("sample_id")
common_samples <- intersect(colnames(immune_wide), rownames(ann_col_immune))
immune_wide    <- immune_wide[, common_samples]
ann_col_immune <- ann_col_immune[common_samples, , drop = FALSE]

pheatmap(immune_wide, annotation_col = ann_col_immune, show_colnames = FALSE,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("white", "#2166AC", "#08306B"))(50),
         main = "Immune cell fractions — TCGA-CHOL (n=34)",
         filename = file.path(out, "results/figures/immune_heatmap.png"), width = 10, height = 6)

# ── PART B: GEO GSE107943 validation cohort (n=30) ─────────
geo_expr <- read.csv(file.path(dl, "GSE107943_expr_symbols.csv"), row.names = 1) %>% as.matrix()
geo_meta <- read.csv(file.path(dl, "external_validation_survival_input.csv"))
common   <- intersect(colnames(geo_expr), geo_meta$sample_id)
geo_expr_sub <- geo_expr[, common]

geo_immune <- deconvolute(geo_expr_sub, method = "quantiseq")
write.csv(geo_immune, file.path(out, "results/de/geo_immune_deconv_scores.csv"), row.names = FALSE)

geo_immune_long <- geo_immune %>%
  pivot_longer(-cell_type, names_to = "sample_id", values_to = "fraction") %>%
  left_join(geo_meta %>% select(sample_id, stage_group, sig_score, risk_group), by = "sample_id")

geo_stage_tests <- geo_immune_long %>%
  group_by(cell_type) %>%
  filter(n_distinct(stage_group) == 2) %>%
  summarise(
    p_value = tryCatch(wilcox.test(fraction ~ stage_group)$p.value, error = function(e) NA_real_),
    n_Early = sum(stage_group == "Early"), n_Late = sum(stage_group == "Late"),
    mean_Early = mean(fraction[stage_group == "Early"], na.rm = TRUE),
    mean_Late  = mean(fraction[stage_group == "Late"],  na.rm = TRUE),
    direction  = ifelse(mean_Late > mean_Early, "UP_in_Late", "DOWN_in_Late"),
    .groups = "drop"
  ) %>% mutate(padj = p.adjust(p_value, method = "BH")) %>% arrange(p_value)
write.csv(geo_stage_tests, file.path(out, "results/de/geo_immune_stage_wilcoxon.csv"), row.names = FALSE)

geo_corr_tests <- geo_immune_long %>%
  group_by(cell_type) %>%
  summarise(
    spearman_r = cor(fraction, sig_score, method = "spearman"),
    p_value    = cor.test(fraction, sig_score, method = "spearman")$p.value,
    .groups = "drop"
  ) %>% mutate(padj = p.adjust(p_value, method = "BH")) %>% arrange(p_value)
write.csv(geo_corr_tests, file.path(out, "results/de/geo_immune_riskscore_correlation.csv"), row.names = FALSE)

p_geo_box <- geo_immune_long %>%
  filter(cell_type %in% c("T cell CD8+", "Myeloid dendritic cell", "NK cell", "Macrophage M1", "Macrophage M2")) %>%
  ggplot(aes(x = stage_group, y = fraction, fill = stage_group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Early" = "#2166AC", "Late" = "#D6604D")) +
  labs(title = "GEO validation (n=30): Immune cell infiltration by stage",
       x = "Stage", y = "Estimated cell fraction", fill = "Stage") +
  theme_classic(base_size = 11) + theme(strip.text = element_text(size = 8), legend.position = "bottom")
ggsave(file.path(out, "results/figures/geo_immune_boxplot.png"), p_geo_box, width = 12, height = 4, dpi = 300)

immune_wide_geo <- geo_immune %>% as.data.frame() %>% remove_rownames() %>%
  column_to_rownames("cell_type") %>% as.matrix()
ann_col_geo <- geo_meta %>% select(sample_id, stage_group, risk_group) %>%
  remove_rownames() %>% column_to_rownames("sample_id")
common_geo <- intersect(colnames(immune_wide_geo), rownames(ann_col_geo))
immune_wide_geo <- immune_wide_geo[, common_geo]
ann_col_geo     <- ann_col_geo[common_geo, , drop = FALSE]

pheatmap(immune_wide_geo, annotation_col = ann_col_geo,
         annotation_colors = list(stage_group = c("Early" = "#92C5DE", "Late" = "#B2182B"),
                                  risk_group  = c("Low" = "#1B7837", "High" = "#762A83")),
         show_colnames = FALSE, clustering_method = "ward.D2",
         color = colorRampPalette(c("white", "#2166AC", "#08306B"))(50),
         main = "Immune cell fractions — GEO GSE107943 (n=30)",
         filename = file.path(out, "results/figures/geo_immune_heatmap.png"), width = 10, height = 6)

cat("GAP 3 COMPLETE — TCGA discovery + GEO validation done.\n")
