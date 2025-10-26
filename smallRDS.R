# Script to create a much smaller RDS for the Shiny app.
# It extracts only the pieces needed at runtime:
#  - gene_lookup: id_to_symbol and symbol_to_id
#  - volcano: per-contrast, per-cellline reduced results tables
#  - counts: normalized counts in long format + sample metadata (condition,treatment,cellline)
#  - stats: within-cellline comparisons used by the plot annotations (group1,group2,p.adj)
#  - gene_ids: vector of gene IDs present
#
# Usage:
# 1) Place this script in the project root (or edit paths below).
# 2) Ensure you have the original large RDS files:
#      Processed Data/fDds_final.RDS
#      Processed Data/gene_lookup.RDS
# 3) Run locally where you have enough memory:
#      Rscript make_small_rds.R
#
# The produced file is: Processed Data/fDds_small.RDS

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# helper used in the app too
revString <- function(text){
  paste(rev(unlist(strsplit(text, NULL))), collapse = "")
}

# Paths (edit if needed)
full_fDds_path <- "Processed Data/fDds_final.RDS"
full_gene_lookup_path <- "Processed Data/gene_lookup.RDS"
out_path <- "Processed Data/fDds_small.RDS"

# Optional filters to reduce size further:
# If you want to keep only genes with baseMean >= MIN_BASEMEAN, set MIN_BASEMEAN > 0
MIN_BASEMEAN <- 0    # set >0 to filter, 0 keeps all genes

# Load full objects
cat("Loading full DESeq2 object...\n")
fDds <- readRDS(full_fDds_path)
cat("Loading gene lookup...\n")
gene_lookup <- readRDS(full_gene_lookup_path)

# Keep vector of gene ids (rownames)
gene_ids_all <- rownames(fDds)

# Optionally filter by baseMean
if (MIN_BASEMEAN > 0) {
  cat("Filtering genes by baseMean >= ", MIN_BASEMEAN, " ...\n", sep = "")
  res_base <- results(fDds)                       # default results to get baseMean
  keep_genes <- rownames(res_base)[res_base$baseMean >= MIN_BASEMEAN]
} else {
  keep_genes <- gene_ids_all
}
cat("Number of genes retained:", length(keep_genes), "\n")

# 1) counts: normalized counts (counts(..., normalized=TRUE))
cat("Extracting normalized counts and sample metadata...\n")
norm_counts <- as.data.frame(counts(fDds, normalized = TRUE)[keep_genes, , drop = FALSE])
counts_long <- norm_counts %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "count")

coldat <- as.data.frame(colData(fDds))
if (!all(c("condition","treatment","cellline") %in% colnames(coldat))) {
  stop("Expected columns condition,treatment,cellline in colData(fDds). Adjust script if different.")
}
coldat$sample <- rownames(coldat)

counts_long <- counts_long %>%
  left_join(coldat %>% select(sample, condition, treatment, cellline), by = "sample") %>%
  mutate(count = as.numeric(count))

# 2) Precompute reduced volcano tables for each contrast/cellline
cat("Computing reduced volcano tables...\n")
contrast_map <- list(
  "Control vs Hypoxia" = list(m = c("VHm","VNm"), t = c("VHt","VNt")),
  "Control vs Fulvestrant" = list(m = c("FHm","VHm"), t = c("FHt","VHt")),
  "Hypoxia vs Hypoxia + Fulvestrant" = list(m = c("FHm","VHm"), t = c("FHt","VHt"))
)

volcano <- list()
for (cname in names(contrast_map)) {
  cm <- contrast_map[[cname]]
  # MCF7 (m)
  res_m <- results(fDds, contrast = c("CoI", cm$m[1], cm$m[2]))
  dfm <- as.data.frame(res_m[keep_genes, , drop = FALSE])
  dfm$gene_id <- rownames(dfm)
  dfm$symbol <- gene_lookup$id_to_symbol[dfm$gene_id]
  dfm$negLog10Padj <- -log10(dfm$padj)
  dfm$meets_fc <- abs(dfm$log2FoldChange) >= 2
  dfm$meets_padj <- dfm$padj < 0.05
  dfm$cellline <- "MCF7"
  dfm_small <- dfm %>% select(gene_id, symbol, log2FoldChange, padj, negLog10Padj, meets_fc, meets_padj, cellline)
  
  # T47D (t)
  res_t <- results(fDds, contrast = c("CoI", cm$t[1], cm$t[2]))
  dft <- as.data.frame(res_t[keep_genes, , drop = FALSE])
  dft$gene_id <- rownames(dft)
  dft$symbol <- gene_lookup$id_to_symbol[dft$gene_id]
  dft$negLog10Padj <- -log10(dft$padj)
  dft$meets_fc <- abs(dft$log2FoldChange) >= 2
  dft$meets_padj <- dft$padj < 0.05
  dft$cellline <- "T47D"
  dft_small <- dft %>% select(gene_id, symbol, log2FoldChange, padj, negLog10Padj, meets_fc, meets_padj, cellline)
  
  volcano[[cname]] <- list(m = dfm_small, t = dft_small)
}

# 3) Stats: within-cellline comparisons used in the app (pairs)
# Added NF vs NV so we can annotate Normoxia Vehicle vs Normoxia Fulvestrant
cat("Computing annotation stats for within-cellline comparisons...\n")
pairs <- list(c("HF","HV"), c("HF","NV"), c("HV","NV"), c("NF","NV"))
stats <- list()
for (cell_tag in c("m", "t")) {
  out_rows <- list()
  for (pair in pairs) {
    group1 <- pair[1]; group2 <- pair[2]
    lvl1 <- paste0(revString(group1), cell_tag)
    lvl2 <- paste0(revString(group2), cell_tag)
    res_pair <- results(fDds, contrast = c("CoI", lvl1, lvl2))
    res_df <- as.data.frame(res_pair[keep_genes, , drop = FALSE])
    res_df$gene_id <- rownames(res_df)
    out_rows[[paste0(group1,"_vs_",group2)]] <- res_df %>%
      select(gene_id, padj) %>%
      rename(p.adj = padj) %>%
      mutate(group1 = group1, group2 = group2) %>%
      select(gene_id, group1, group2, p.adj)
  }
  stats[[cell_tag]] <- bind_rows(out_rows)
  rownames(stats[[cell_tag]]) <- NULL
}

# Build and save the reduced object
small_obj <- list(
  gene_lookup = gene_lookup,
  volcano = volcano,
  counts = counts_long,
  stats = stats,
  gene_ids = keep_genes
)

cat("Saving smaller RDS to:", out_path, "\n")
saveRDS(small_obj, file = out_path)
cat("Done.\n")