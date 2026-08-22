library(readr)
library(dplyr)

iCCA_deg <- read_csv("C:/Users/kaavy/Downloads/res_iCCA_stage (1).csv")
eCCA_deg <- read_csv("C:/Users/kaavy/Downloads/res_eCCA_stage (1).csv")
symbol_map <- read_csv("C:/Users/kaavy/Downloads/gene_symbol_map_master (1).csv")

#inspecting column names
names(iCCA_deg)
names(eCCA_deg)
names(symbol_map)

#standardize ENSEMBL IDs, remove suffixes
strip_ensembl_version <- function(x) {
  sub("\\..*$", "", as.character(x))
}

iCCA_deg <- iCCA_deg %>%
  mutate(gene_id_clean = strip_ensembl_version(gene_id))

eCCA_deg <- eCCA_deg %>%
  mutate(gene_id_clean = strip_ensembl_version(gene_id))

symbol_map <- symbol_map %>%
  mutate(gene_id_clean = strip_ensembl_version(gene_id))

#remove duplicates & missing symbols
symbol_map_clean <- symbol_map %>%
  dplyr::select(gene_id_clean, gene_symbol) %>%
  dplyr::filter(
    !is.na(gene_id_clean),
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  dplyr::distinct(gene_id_clean, .keep_all = TRUE)

iCCA_deg_symbol <- iCCA_deg %>%
  left_join(symbol_map_clean, by = "gene_id_clean")
eCCA_deg_symbol <- eCCA_deg %>%
  left_join(symbol_map_clean, by = "gene_id_clean")

#map gene symbols to ENSEMBL IDs
data.frame(
  dataset = c("iCCA", "eCCA"),
  total_genes = c(
    nrow(iCCA_deg_symbol),
    nrow(eCCA_deg_symbol)
  ),
  mapped_genes = c(
    sum(!is.na(iCCA_deg_symbol$gene_symbol)),
    sum(!is.na(eCCA_deg_symbol$gene_symbol))
  )
)

mean(!is.na(iCCA_deg_symbol$gene_symbol)) * 100
mean(!is.na(eCCA_deg_symbol$gene_symbol)) * 100

names(iCCA_deg_symbol)
names(eCCA_deg_symbol)

#finding significant genes
iCCA_sig <- iCCA_deg_symbol %>%
  dplyr::filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    !is.na(padj),
    !is.na(log2FC),
    padj < 0.05,
    abs(log2FC) > 1
  )

eCCA_sig <- eCCA_deg_symbol %>%
  dplyr::filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    !is.na(padj),
    !is.na(log2FC),
    padj < 0.05,
    abs(log2FC) > 1
  )

nrow(iCCA_sig)
nrow(eCCA_sig)
table(iCCA_sig$direction)
table(eCCA_sig$direction)


unique(iCCA_deg_symbol$direction)
unique(eCCA_deg_symbol$direction)

table(iCCA_deg_symbol$direction, useNA = "always")
table(eCCA_deg_symbol$direction, useNA = "always")

#top genes list
get_top_genes <- function(deg, direction_value, n = 50) {
  deg %>%
    dplyr::filter(
      .data$direction == direction_value,
      !is.na(gene_symbol),
      gene_symbol != "",
      !is.na(padj),
      padj < 0.05,
      abs(log2FC) > 1
    ) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FC))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(gene_symbol)
}

iCCA_up <- get_top_genes(iCCA_deg_symbol, "UP")
iCCA_down <- get_top_genes(iCCA_deg_symbol, "DOWN")

eCCA_up <- get_top_genes(eCCA_deg_symbol, "UP")
eCCA_down <- get_top_genes(eCCA_deg_symbol, "DOWN")
length(iCCA_up)
length(iCCA_down)
length(eCCA_up)
length(eCCA_down)
writeLines(iCCA_up, "top_iCCA_UP.txt")
writeLines(iCCA_down, "top_iCCA_DOWN.txt")

writeLines(eCCA_up, "top_eCCA_UP.txt")
writeLines(eCCA_down, "top_eCCA_DOWN.txt")

iCCA_up <- readLines("top_iCCA_UP.txt")
iCCA_down <- readLines("top_iCCA_DOWN.txt")
eCCA_up <- readLines("top_eCCA_UP.txt")
eCCA_down <- readLines("top_eCCA_DOWN.txt")

length(iCCA_up)
length(iCCA_down)
length(eCCA_up)
length(eCCA_down)

head(iCCA_up)

clean_gene_list <- function(x) {
  x <- trimws(x)
  x <- x[x != ""]
  unique(x)
}

#input files for ChEA3
iCCA_up <- clean_gene_list(iCCA_up)
iCCA_down <- clean_gene_list(iCCA_down)
eCCA_up <- clean_gene_list(eCCA_up)
eCCA_down <- clean_gene_list(eCCA_down)

install.packages("rChEA3")
library(rChEA3)

tf_iCCA_up <- queryChEA3(
  iCCA_up,
  query_name = "iCCA_UP"
)

tf_iCCA_down <- queryChEA3(
  iCCA_down,
  query_name = "iCCA_DOWN"
)

tf_eCCA_up <- queryChEA3(
  eCCA_up,
  query_name = "eCCA_UP"
)

tf_eCCA_down <- queryChEA3(
  eCCA_down,
  query_name = "eCCA_DOWN"
)
names(tf_iCCA_up)
names(tf_iCCA_down)
str(tf_iCCA_up, max.level = 1)
saveRDS(tf_iCCA_up, "ChEA3_iCCA_UP_raw.rds")
saveRDS(tf_iCCA_down, "ChEA3_iCCA_DOWN_raw.rds")
saveRDS(tf_eCCA_up, "ChEA3_eCCA_UP_raw.rds")
saveRDS(tf_eCCA_down, "ChEA3_eCCA_DOWN_raw.rds")

names(tf_iCCA_up)
#ranking queried results
tf_iCCA_up_table <- as.data.frame(
  tf_iCCA_up[["Integrated--meanRank"]]
)

tf_iCCA_down_table <- as.data.frame(
  tf_iCCA_down[["Integrated--meanRank"]]
)

tf_eCCA_up_table <- as.data.frame(
  tf_eCCA_up[["Integrated--meanRank"]]
)

tf_eCCA_down_table <- as.data.frame(
  tf_eCCA_down[["Integrated--meanRank"]]
)

write.csv(
  tf_iCCA_up_table,
  "ChEA3_iCCA_UP_integrated.csv",
  row.names = FALSE
)

write.csv(
  tf_iCCA_down_table,
  "ChEA3_iCCA_DOWN_integrated.csv",
  row.names = FALSE
)

write.csv(
  tf_eCCA_up_table,
  "ChEA3_eCCA_UP_integrated.csv",
  row.names = FALSE
)

write.csv(
  tf_eCCA_down_table,
  "ChEA3_eCCA_DOWN_integrated.csv",
  row.names = FALSE
)

head(tf_iCCA_up_table, 20)
head(tf_iCCA_down_table, 20)
head(tf_eCCA_up_table, 20)
head(tf_eCCA_down_table, 20)

names(tf_iCCA_up_table)

#saving top 20 TFs
top_TFs <- dplyr::bind_rows(
  tf_iCCA_up_table %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(subtype = "iCCA", direction = "UP"),
  
  tf_iCCA_down_table %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(subtype = "iCCA", direction = "DOWN"),
  
  tf_eCCA_up_table %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(subtype = "eCCA", direction = "UP"),
  
  tf_eCCA_down_table %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(subtype = "eCCA", direction = "DOWN")
)

#final results
write.csv(top_TFs, "top_TFs_all_four_lists.csv", row.names = FALSE)

