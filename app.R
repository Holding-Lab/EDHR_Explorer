library(shiny)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(plotly)

############################################################
## Helper functions (unchanged behaviour)
############################################################

mapIDsToSymbols <- function(ids) {
  ids <- as.character(ids)
  syms <- gene_lookup$id_to_symbol[ids]
  # fallback: if we don't know the symbol, just show the ID
  syms[is.na(syms)] <- ids[is.na(syms)]
  unname(syms)
}

revString <- function(text){
  paste(rev(unlist(strsplit(text,NULL))),collapse="")
}

# Pretty labels used for plotting / faceting
CoI_labels <- c(
  'HV'='Hypoxia & Vehicle',
  'HF'='Hypoxia & Fulvestrant',
  'NV'='Normoxia & Vehicle',
  'NF'='Normoxia & Fulvestrant'
)

cellline_labels <- c('m'='MCF7','t'='T47D')

########################
# Volcano helpers
########################

# Interactive volcano with POI on top (in red)
makeInteractiveVolcano <- function(df, poi_gene_query){
  
  # mark poi / sig / bg
  poi_mask <- (df$gene_id == poi_gene_query) | (df$symbol == poi_gene_query)
  sig_mask <- (!poi_mask) & df$meets_fc & df$meets_padj
  bg_mask  <- (!poi_mask) & !sig_mask
  
  # hover text
  safe_padj <- ifelse(is.na(df$padj), "NA", signif(df$padj,3))
  safe_neglog10 <- ifelse(is.na(df$negLog10Padj), "NA", signif(df$negLog10Padj,3))
  hover_txt <- function(idx_vec) {
    paste0(
      "Gene: ",
      ifelse(is.na(df$symbol[idx_vec]), df$gene_id[idx_vec], df$symbol[idx_vec]), "<br>",
      "log2FC: ", signif(df$log2FoldChange[idx_vec],3), "<br>",
      "adj.P: ", safe_padj[idx_vec], "<br>",
      "-log10(adj.P): ", safe_neglog10[idx_vec]
    )
  }
  
  p <- plot_ly()
  
  # background (light grey)
  if (any(bg_mask, na.rm=TRUE)) {
    idx <- which(bg_mask)
    p <- add_trace(
      p,
      x = df$log2FoldChange[idx],
      y = df$negLog10Padj[idx],
      type = "scatter",
      mode = "markers",
      text = hover_txt(idx),
      hoverinfo = "text",
      marker = list(
        size = 6,
        color = "lightgrey"
      ),
      showlegend = FALSE,
      name = "background"
    )
  }
  
  # significant (black)
  if (any(sig_mask, na.rm=TRUE)) {
    idx <- which(sig_mask)
    p <- add_trace(
      p,
      x = df$log2FoldChange[idx],
      y = df$negLog10Padj[idx],
      type = "scatter",
      mode = "markers",
      text = hover_txt(idx),
      hoverinfo = "text",
      marker = list(
        size = 6,
        color = "black"
      ),
      showlegend = FALSE,
      name = "significant"
    )
  }
  
  # point of interest (red, bigger, last so always on top)
  if (any(poi_mask, na.rm=TRUE)) {
    idx <- which(poi_mask)
    p <- add_trace(
      p,
      x = df$log2FoldChange[idx],
      y = df$negLog10Padj[idx],
      type = "scatter",
      mode = "markers",
      text = hover_txt(idx),
      hoverinfo = "text",
      marker = list(
        size = 8,
        color = "red",
        line = list(width = 1, color = "red")
      ),
      showlegend = FALSE,
      name = "POI"
    )
  }
  
  # cutoff lines: |LFC| >= 2 and padj < 0.05
  p <- add_lines(
    p,
    x = c(-2,-2),
    y = c(0, max(df$negLog10Padj, na.rm=TRUE)),
    inherit = FALSE,
    line = list(dash = "dash", width = 1, color = "black"),
    showlegend = FALSE,
    hoverinfo = "none"
  )
  
  p <- add_lines(
    p,
    x = c(2,2),
    y = c(0, max(df$negLog10Padj, na.rm=TRUE)),
    inherit = FALSE,
    line = list(dash = "dash", width = 1, color = "black"),
    showlegend = FALSE,
    hoverinfo = "none"
  )
  
  p <- add_lines(
    p,
    x = c(min(df$log2FoldChange, na.rm=TRUE),
          max(df$log2FoldChange, na.rm=TRUE)),
    y = rep(-log10(0.05), 2),
    inherit = FALSE,
    line = list(dash = "dash", width = 1, color = "black"),
    showlegend = FALSE,
    hoverinfo = "none"
  )
  
  p <- layout(
    p,
    title = unique(df$cellline),
    xaxis = list(title = "log2 Fold Change"),
    yaxis = list(title = "-log10(adj p)"),
    showlegend = FALSE
  )
  
  p
}

########################
# Gene expression / stats helpers (adapted to use the smaller RDS)
########################

# Pull per-sample counts + metadata for one gene from the small counts table
getGeneCountsDf <- function(geneID) {
  d <- counts_df %>% dplyr::filter(gene_id == as.character(geneID))
  d$combined <- paste0(d$condition, d$treatment)
  d
}

# Compute DESeq2 padj for a set of within-cell-line comparisons and place them vertically
# now uses precomputed stats_store (long format)
getDeSeqStat_oneCell <- function(geneID, cell_tag, ypos_base) {
  pairs <- list(
    c("HF","HV"),
    c("NF","NV"),
    c("HV","NV")
  )
  
  stat_long <- stats_store[[cell_tag]] %>% dplyr::filter(gene_id == as.character(geneID))
  
  out <- Map(function(pair, i){
    group1 <- pair[1]
    group2 <- pair[2]
    
    p.adj.val <- stat_long %>%
      dplyr::filter(group1 == !!group1, group2 == !!group2) %>%
      dplyr::pull(p.adj)
    if (length(p.adj.val) == 0) p.adj.val <- NA
    
    data.frame(
      group1 = group1,
      group2 = group2,
      p.adj = signif(p.adj.val, 2),
      y.position = ypos_base * (1 + 0.1*(i-1)),
      cellline = cell_tag,
      stringsAsFactors = FALSE
    )
  }, pairs, seq_along(pairs))
  
  do.call(rbind, out)
}

# Plot one cell line (either m or t) with optional fixed y cap and with that
# cell line's stats table (unchanged)
plotGene_oneCell <- function(d_sub, geneID, cell_tag, ymax_override = NULL, stats_df = NULL) {
  
  # 1. figure out how tall the data go
  data_max <- max(d_sub$count, na.rm = TRUE)
  
  # 2. figure out how tall the annotations want to go
  if (!is.null(stats_df) && nrow(stats_df) > 0) {
    anno_max <- max(stats_df$y.position, na.rm = TRUE)
  } else {
    anno_max <- data_max
  }
  
  # 3. pick a target max (before global override)
  # a little padding, e.g. +10%
  local_ceiling <- max(data_max, anno_max) * 1.1
  
  # 4. if caller provided an override (shared scale mode), use that instead
  final_ceiling <- if (!is.null(ymax_override)) ymax_override else local_ceiling
  
  p <- ggplot(d_sub, aes(x=combined, y=count)) +
    geom_boxplot(aes(fill=condition), linewidth = 0.5) +
    ggtitle(paste0(mapIDsToSymbols(geneID), " (", cellline_labels[[cell_tag]], ")")) +
    theme_pubr() +
    scale_x_discrete(labels=CoI_labels) +
    labs(y= "Normalised Counts", x = "Condition & Treatment") +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, final_ceiling))
  
  if (!is.null(stats_df) && nrow(stats_df) > 0) {
    p <- p + stat_pvalue_manual(
      stats_df,
      size = 2.8,
      unit = "pt",
      label = "p = {p.adj}"
    )
  }
  
  p +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

# Make the combined 2-panel plot with synced or free y-scales
plotGeneCombined <- function(geneID, lock_scale = TRUE) {
  
  d_all <- getGeneCountsDf(geneID)
  
  d_m <- d_all %>% dplyr::filter(cellline == "m")
  d_t <- d_all %>% dplyr::filter(cellline == "t")
  
  # per-cell-line maxima
  ymax_m <- max(d_m$count, na.rm=TRUE) * 1.2
  ymax_t <- max(d_t$count, na.rm=TRUE) * 1.2
  
  # choose per-plot y-limits
  if (isTRUE(lock_scale)) {
    common_max <- max(ymax_m, ymax_t)
    stats_m <- getDeSeqStat_oneCell(geneID, "m", ypos_base = common_max)
    stats_t <- getDeSeqStat_oneCell(geneID, "t", ypos_base = common_max)
    ylim_m <- common_max*1.2
    ylim_t <- common_max*1.2
  } else {
    stats_m <- getDeSeqStat_oneCell(geneID, "m", ypos_base = ymax_m)
    stats_t <- getDeSeqStat_oneCell(geneID, "t", ypos_base = ymax_t)
    ylim_m <- NULL
    ylim_t <- NULL
  }
  
  
  p_m <- plotGene_oneCell(
    d_sub = d_m,
    geneID = geneID,
    cell_tag = "m",
    ymax_override = ylim_m,
    stats_df = stats_m
  )
  
  p_t <- plotGene_oneCell(
    d_sub = d_t,
    geneID = geneID,
    cell_tag = "t",
    ymax_override = ylim_t,
    stats_df = stats_t
  )
  
  ggpubr::ggarrange(p_m, p_t, nrow = 1)
}

############################################################
## Data load (smaller precomputed object)
############################################################

# NEW: read the reduced RDS prepared by the helper script (see make_small_rds.R)
small <- readRDS("Processed Data/fDds_small.RDS")

# components expected in the small object:
# small$gene_lookup: list with id_to_symbol and symbol_to_id
# small$volcano: named list of contrasts, each a list with $m and $t data.frames (columns same as before)
# small$counts: long data.frame with columns gene_id,sample,count,condition,treatment,cellline
# small$stats: list with $m and $t long data.frames columns: gene_id, group1, group2, p.adj
# small$gene_ids: character vector of gene IDs present
gene_lookup <- small$gene_lookup
volcano_store <- small$volcano
counts_df <- small$counts
stats_store <- small$stats
gene_ids <- small$gene_ids

# biological -> CoI mapping for each cell line (unchanged)
contrast_map <- list(
  "Control vs Hypoxia" = list(
    m = c("VHm","VNm"),  # MCF7: Vehicle Hypoxia vs Vehicle Normoxia
    t = c("VHt","VNt")   # T47D: Vehicle Hypoxia vs Vehicle Normoxia
  ),
  "Control vs Fulvestrant" = list(
    m = c("FHm","VHm"),  # Fulvestrant Hypoxia vs Vehicle Hypoxia
    t = c("FHt","VHt")   # Fulvestrant Hypoxia vs Vehicle Hypoxia
  ),
  "Hypoxia vs Hypoxia + Fulvestrant" = list(
    m = c("FHm","VHm"),  # hypoxia +/- fulvestrant
    t = c("FHt","VHt")
  )
)

############################################################
## UI (unchanged)
############################################################

ui <- fluidPage(
  titlePanel("EDHR explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "contrast_choice",
        "Select contrast:",
        choices = names(contrast_map),
        selected = "Control vs Hypoxia"
      ),
      textInput(
        "gene_search",
        "Gene (symbol or Entrez ID):",
        value = "SCNN1B"
      ),
      fluidRow(
        column(
          width = 7,
          checkboxInput(
            "lock_scale",
            "Match y-axis",
            value = TRUE
          )
        ),
        column(
          width = 5,
          actionButton("go_gene", "Plot gene")
        )
      ),
      helpText("Click 'Plot gene' to refresh highlight & expression. Tick 'Match y-axis' to force both panels to same Y range.")
    ),
    mainPanel(
      fluidRow(
        column(
          width = 6,
          plotlyOutput("volcano_m")
        ),
        column(
          width = 6,
          plotlyOutput("volcano_t")
        )
      ),
      hr(),
      h3(textOutput("gene_title")),
      plotOutput("gene_boxplot", height = "500px")
    )
  )
)

############################################################
## Server (adapted to use precomputed, smaller data)
############################################################

server <- function(input, output, session){
  
  # freeze the requested gene at the moment "Plot gene" is clicked
  selected_gene <- eventReactive(input$go_gene, {
    input$gene_search
  }, ignoreInit = FALSE)  # will run once on load
  
  # recompute volcano data whenever contrast changes (use precomputed store)
  volcano_data <- reactive({
    # take the precomputed tables for this named contrast
    vs <- volcano_store[[ input$contrast_choice ]]
    # vs is expected to be a list with $m and $t data.frames
    list(m = vs$m, t = vs$t)
  })
  
  # volcano plots (highlight based on last confirmed gene)
  output$volcano_m <- renderPlotly({
    makeInteractiveVolcano(
      volcano_data()$m,
      poi_gene_query = selected_gene()
    )
  })
  
  output$volcano_t <- renderPlotly({
    makeInteractiveVolcano(
      volcano_data()$t,
      poi_gene_query = selected_gene()
    )
  })
  
  # title above gene stats plot
  output$gene_title <- renderText({
    paste0("Per-gene expression for: ", selected_gene())
  })
  
  # boxplot + stats: now uses combined two-panel plot with proper scaling logic
  output$gene_boxplot <- renderPlot({
    gene_query <- selected_gene()
    
    # Accept Entrez ID rowname or gene symbol (no DB call)
    target_row <- NULL
    gene_query_clean <- gene_query
    
    # Case 1: user typed an Entrez ID that matches gene_ids
    if (gene_query_clean %in% gene_ids) {
      target_row <- gene_query_clean
    } else {
      # Case 2: user typed a gene symbol that we know
      if (gene_query_clean %in% names(gene_lookup$symbol_to_id)) {
        candidate <- gene_lookup$symbol_to_id[[gene_query_clean]]
        if (candidate %in% gene_ids) {
          target_row <- candidate
        }
      }
    }
    
    if (is.null(target_row)) {
      validate(
        need(FALSE,
             paste0("Couldn't find gene '",gene_query,
                    "' in this dataset. Try Entrez ID or another symbol?"))
      )
    }
    
    plotGeneCombined(
      geneID = target_row,
      lock_scale = input$lock_scale
    )
  })
}

############################################################
## Launch (unchanged)
############################################################

shinyApp(ui, server)