library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

fDds <- readRDS("Processed Data/fDds_final.RDS")

gene_ids <- rownames(fDds)

id_to_symbol_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENTREZID",
  columns = c("SYMBOL")
)

# Entrez -> symbol
id_to_symbol_tbl <- id_to_symbol_raw %>%
  group_by(ENTREZID) %>%
  summarise(SYMBOL = SYMBOL[1], .groups = "drop")

id_to_symbol <- id_to_symbol_tbl$SYMBOL
names(id_to_symbol) <- id_to_symbol_tbl$ENTREZID

# symbol -> Entrez
symbol_to_id_tbl <- id_to_symbol_raw %>%
  filter(!is.na(SYMBOL), !is.na(ENTREZID)) %>%
  group_by(SYMBOL) %>%
  summarise(ENTREZID = ENTREZID[1], .groups = "drop")

symbol_to_id <- symbol_to_id_tbl$ENTREZID
names(symbol_to_id) <- symbol_to_id_tbl$SYMBOL

gene_lookup <- list(
  id_to_symbol = id_to_symbol,   # names = Entrez IDs, values = Symbols
  symbol_to_id = symbol_to_id    # names = Symbols,   values = Entrez IDs
)

saveRDS(gene_lookup, file = "Processed Data/gene_lookup.RDS")
