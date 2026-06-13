library(DESeq2)
library(data.table)
library(dplyr)

dir.create("results/processed", recursive = TRUE, showWarnings = FALSE)

sample_meta_master <- fread("C:/Users/navni/Downloads/sample_meta_master (2).csv")

star_counts <- fread("C:/Users/navni/Downloads/star_counts_5samples_gene_by_sample (1).csv")

star_counts <- star_counts %>% 
  filter(!gene_id %in% c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous"))

star_counts$gene_id <- sub("\\..*", "", star_counts$gene_id)

if (any(duplicated(star_counts$gene_id))) {
  star_counts <- star_counts %>% 
    group_by(gene_id) %>% 
    summarise(across(everything(), sum)) %>% 
    ungroup()
}

gene_ids <- star_counts$gene_id
tcga_counts <- as.matrix(star_counts %>% select(-gene_id))
rownames(tcga_counts) <- gene_ids

# strip vial letter: TCGA-XX-XXXX-01A -> TCGA-XX-XXXX-01
colnames(tcga_counts) <- substr(colnames(tcga_counts), 1, 15)

dim(tcga_counts)
head(tcga_counts[, 1:5])
colnames(tcga_counts)
rownames(tcga_counts)[1:10]

# check before re-running, or just verify directly:
length(gene_ids)
length(unique(gene_ids))

tcga_meta <- sample_meta_master %>% 
  filter(dataset == "TCGA", has_star_counts == 1)

common_samples <- intersect(colnames(tcga_counts), tcga_meta$sample_id)
tcga_counts <- tcga_counts[, common_samples]
tcga_meta <- tcga_meta[match(common_samples, tcga_meta$sample_id), ]

stopifnot(all(colnames(tcga_counts) == tcga_meta$sample_id))
stopifnot(!any(is.na(tcga_meta$sample_id)))
stopifnot(length(common_samples) > 0)

tcga_meta$label <- factor(tcga_meta$label)

dim(tcga_counts)
tcga_meta

dds <- DESeqDataSetFromMatrix(
  countData = tcga_counts,
  colData = as.data.frame(tcga_meta),
  design = ~ 1
)

keep <- rowSums(counts(dds) >= 10) >= (0.1 * ncol(dds))
dds <- dds[keep, ]

dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

dim(dds)
sizeFactors(dds)
dim(vsd)

geo_expr <- fread("C:/Users/navni/Downloads/geo_expression_labelled.csv")

colnames(geo_expr)[1:5]

head(geo_expr[, 1:5])
dim(geo_expr)

geo_expr <- fread("C:/Users/navni/Downloads/geo_expression_labelled.csv", header = TRUE)

# Extract the embedded label row (row 1 of data, since header was already consumed)
geo_labels_row <- geo_expr[1, ]
geo_expr <- geo_expr[-1, ]  # remove label row from the data

# Set probe IDs as rownames
probe_ids <- geo_expr[[1]]
geo_matrix <- as.matrix(geo_expr[, -1])
rownames(geo_matrix) <- probe_ids

# Convert to numeric (in case it's still character type)
geo_matrix <- apply(geo_matrix, 2, as.numeric)
rownames(geo_matrix) <- probe_ids

# Build a labels data frame from the extracted row
geo_sample_labels <- data.frame(
  sample_id = colnames(geo_matrix),
  label = as.character(geo_labels_row[1, -1])
)

dim(geo_matrix)
head(geo_sample_labels)

geo_meta <- sample_meta_master %>% filter(dataset == "GEO")

common_geo_samples <- intersect(colnames(geo_matrix), geo_meta$sample_id)
length(common_geo_samples)

geo_matrix <- geo_matrix[, common_geo_samples]
geo_meta <- geo_meta[match(common_geo_samples, geo_meta$sample_id), ]

stopifnot(all(colnames(geo_matrix) == geo_meta$sample_id))
stopifnot(!any(is.na(geo_meta$sample_id)))
stopifnot(length(common_geo_samples) > 0)

dim(geo_matrix)

keep_probes <- rowMeans(geo_matrix) > 4
geo_matrix_filtered <- geo_matrix[keep_probes, ]

dim(geo_matrix)
dim(geo_matrix_filtered)

boxplot(geo_matrix_filtered[, 1:10], main = "Check distributions across samples")

saveRDS(dds, "results/processed/tcga_dds.rds")
saveRDS(vsd, "results/processed/tcga_vsd.rds")
saveRDS(tcga_meta, "results/processed/tcga_meta.rds")
saveRDS(geo_matrix_filtered, "results/processed/geo_expr_filtered.rds")
saveRDS(geo_meta, "results/processed/geo_meta.rds")

fwrite(as.data.frame(counts(dds, normalized = TRUE)), 
       "results/processed/tcga_normalized_counts.csv", row.names = TRUE)

fwrite(as.data.frame(assay(vsd)), 
       "results/processed/tcga_vst_matrix.csv", row.names = TRUE)

fwrite(as.data.frame(geo_matrix_filtered), 
       "results/processed/geo_expr_filtered.csv", row.names = TRUE)

list.files("results/processed/")

getwd()