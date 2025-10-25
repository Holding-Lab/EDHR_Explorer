##############################
## app.R
##############################

library(shiny)
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)
library(ggpubr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(plotly)

############################################################
## Helper functions
############################################################

mapIDsToSymbols <- function(ids) {
  eg2symmap <- as.list(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])
  ids <- as.character(ids)
  mapped_symbols <- eg2symmap[ids]
  names(mapped_symbols) <- ids
  symbols <- sapply(mapped_symbols, function(x) if (is.null(x)) NA else x[1])
  setNames(unlist(symbols), names(mapped_symbols))
}

revString <- function(text){
  paste(rev(unlist(strsplit(text,NULL))),collapse="")
}

plotGene <- function(fDds,geneID) {
  facets_labels <- c('m'='MCF7','t'='T47D')
  CoI_labels <- c('HV'='Hypoxia & Vehicle',
                  'HF'='Hypoxia & Fulvestrant',
                  'NV'='Normoxia & Vehicle',
                  'NF'='Normoxia & Fulvestrant')
  
  d <- plotCounts(
    fDds,
    gene = as.character(geneID),
    intgroup = c("CoI","cellline","condition","treatment"),
    returnData = TRUE,
    transform = FALSE
  )
  
  d$combined <- paste0(d$condition,d$treatment)
  
  ggplot(d, aes(x=combined, y=count)) +
    geom_boxplot(aes(fill=condition), linewidth = 0.5) +
    facet_grid(cellline ~ ., labeller = as_labeller(facets_labels), scales="free") +
    ggtitle(mapIDsToSymbols(geneID)) +
    theme_pubr() +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(labels=CoI_labels) +
    labs(y= "Normalised Counts", x = "Condition & Treatment") +
    theme(legend.position = "none")
}

getDeSeqStatCompare<-function(fDds, cellline, group1, group2, geneID, ypos) {
  stat.res <- results(
    fDds,
    contrast = c(
      "CoI",
      paste0(revString(group1),cellline),
      paste0(revString(group2),cellline)
    )
  )
  data.frame(
    cellline = cellline,
    group1 = group1,
    group2 = group2,
    p.adj = signif(stat.res[as.character(geneID),'padj'],2),
    y.position = ypos
  )
}

getDeSeqStat<-function(fDds,geneID,ypos) {
  rbind(
    getDeSeqStatCompare(fDds, 'm',"HF","HV",geneID,ypos),
    getDeSeqStatCompare(fDds, 'm',"HF","NV",geneID,ypos*1.1),
    getDeSeqStatCompare(fDds, 'm',"HV","NV",geneID,ypos*1.2),
    getDeSeqStatCompare(fDds, 't',"HF","HV",geneID,ypos),
    getDeSeqStatCompare(fDds, 't',"HF","NV",geneID,ypos*1.1),
    getDeSeqStatCompare(fDds, 't',"HV","NV",geneID,ypos*1.2)
  )
}

plotGeneStats<-function(fDds,geneID,ypos=1000) {
  plotGene(fDds, geneID) +
    stat_pvalue_manual(
      getDeSeqStat(fDds,geneID,ypos),
      size=2.8,
      unit="pt",
      label = "p = {p.adj}"
    )
}

# Build volcano data frame for a given contrast level in a given cell line
getVolcanoTable <- function(fDds, coi_A, coi_B, which_cellline){
  res <- results(fDds, contrast=c("CoI", coi_A, coi_B))
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df$symbol <- mapIDsToSymbols(df$gene_id)
  df$negLog10Padj <- -log10(df$padj)
  df$cellline <- which_cellline
  
  df$meets_fc   <- abs(df$log2FoldChange) >= 2
  df$meets_padj <- df$padj < 0.05
  
  df
}

# Interactive volcano with POI on top
makeInteractiveVolcano <- function(df, poi_gene_query){
  
  # mark poi / sig / bg
  poi_mask <- (df$gene_id == poi_gene_query) | (df$symbol == poi_gene_query)
  sig_mask <- (!poi_mask) & df$meets_fc & df$meets_padj
  bg_mask  <- (!poi_mask) & !sig_mask
  
  # safe hover text
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
  
  # point of interest (red, bigger, last)
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
  
  # cutoff lines
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

############################################################
## Data load + DESeq2 setup
############################################################

mcf7 <- readRDS("Processed Data/countMatrix-MCF7.RDS")
t47d <- readRDS("Processed Data/countMatrix-t47D.RDS")

mcf7Counts <- mcf7$counts
t47dCounts <- t47d$counts

mcf7Samplenames <- matrix(nrow=2, unlist(strsplit(colnames(mcf7Counts),"\\.")))[1,]

t47dSamplenames <- sub(
  "_2","",
  matrix(nrow=2, unlist(strsplit(colnames(t47dCounts),"\\.")))[1,]
)
t47dSamplenames <- paste0(
  substr(t47dSamplenames,3,3),
  substr(t47dSamplenames,1,2),
  substr(t47dSamplenames,4,4)
)

colnames(mcf7Counts) <- paste0(mcf7Samplenames,"_m")
colnames(t47dCounts) <- paste0(t47dSamplenames,"_t")

generateMetadata <- function(samplenames,cellline) {
  data.frame(
    treatment = substr(samplenames,1,1),
    condition = substr(samplenames,nchar(samplenames),nchar(samplenames)),
    rep = gsub("[^0-9.-]", "", samplenames),
    cellline = cellline,
    row.names = paste0(samplenames,"_",cellline)
  )
}

coldata <- rbind(
  generateMetadata(mcf7Samplenames,'m'),
  generateMetadata(t47dSamplenames,'t')
)

counts <- cbind(mcf7Counts,t47dCounts)
coldata <- coldata[colnames(counts),,drop=FALSE]

coldata$CoI <- paste0(
  coldata$treatment,
  coldata$condition,
  coldata$cellline
)

removeSamples <- c(grep("4",colnames(counts), value=TRUE),"V21H_t")
keep <- !colnames(counts) %in% removeSamples
filteredCounts <- counts[, keep, drop=FALSE]
filteredColdata <- coldata[keep, , drop=FALSE]

fDds <- DESeqDataSetFromMatrix(
  countData = filteredCounts,
  colData = filteredColdata,
  design = ~ CoI
)
fDds <- DESeq(fDds)

# biological -> CoI mapping for each cell line
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
## UI
############################################################

ui <- fluidPage(
  titlePanel("MCF7 vs T47D dual-volcano explorer"),
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
      actionButton("go_gene", "Update gene"),
      helpText("Click 'Update gene' to refresh highlight & expression.")
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
## Server
############################################################

server <- function(input, output, session){
  
  # freeze the requested gene at the moment "Update gene" is clicked
  selected_gene <- eventReactive(input$go_gene, {
    input$gene_search
  }, ignoreInit = FALSE)  # so it works on first load
  
  # recompute volcano data when contrast changes
  volcano_data <- reactive({
    cm <- contrast_map[[ input$contrast_choice ]]
    
    vm <- getVolcanoTable(
      fDds,
      coi_A = cm$m[1],
      coi_B = cm$m[2],
      which_cellline = "MCF7"
    )
    
    vt <- getVolcanoTable(
      fDds,
      coi_A = cm$t[1],
      coi_B = cm$t[2],
      which_cellline = "T47D"
    )
    
    list(m = vm, t = vt)
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
  
  # title updates with last confirmed gene
  output$gene_title <- renderText({
    paste0("Per-gene expression for: ", selected_gene())
  })
  
  # boxplot + stats also use last confirmed gene
  output$gene_boxplot <- renderPlot({
    gene_query <- selected_gene()
    
    # match Entrez ID or symbol
    target_row <- NULL
    if (gene_query %in% rownames(fDds)) {
      target_row <- gene_query
    } else {
      sym2entrez <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = gene_query,
        keytype = "SYMBOL",
        columns = c("ENTREZID","SYMBOL")
      )
      if (!is.null(sym2entrez$ENTREZID)) {
        match_id <- sym2entrez$ENTREZID[1]
        if (match_id %in% rownames(fDds)) {
          target_row <- match_id
        }
      }
    }
    
    if (is.null(target_row)) {
      validate(
        need(FALSE,
             paste0("Couldn't find gene '",gene_query,
                    "' in this DESeqDataSet. Try Entrez ID?"))
      )
    }
    
    dtmp <- plotCounts(
      fDds,
      gene = target_row,
      intgroup = c("CoI","cellline","condition","treatment"),
      returnData = TRUE,
      transform = FALSE
    )
    yguess <- max(dtmp$count, na.rm=TRUE) * 1.2
    
    plotGeneStats(fDds, target_row, ypos = yguess)
  })
}

############################################################
## Launch
############################################################

shinyApp(ui, server)
