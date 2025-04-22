setwd("../proteomics_GSEA_GO_Reactome/")

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(enrichplot)
library(readxl)
library(writexl)
library(forcats)
library(ggstance)
library(msigdbr)
library(here)
library(future)
library(tidyverse)

# Define gene set configurations
gene_sets_config <- list(
  C5_GO_BP = list(category = "C5", subcategory = "GO:BP"),
  C2_REACTOME = list(category = "C2", subcategory = "CP:REACTOME")
)

# Get all CSV files in the folder
csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file
for (file in csv_files) {
  sample_name <- tools::file_path_sans_ext(basename(file))
  cat("Processing:", sample_name, "\n")
  
  # Load data
  deg <- read.csv(file)
  
  # Check column names exist
  if (!all(c("AVG.Log2.Ratio", "Qvalue", "Genes") %in% colnames(deg))) {
    warning(paste("Skipping", sample_name, "- required columns missing"))
    next
  }
  
  # Filter genes
  deg_filtered <- deg %>%
    filter(abs(AVG.Log2.Ratio) > 0.25 & Qvalue < 0.05)
  
  if (nrow(deg_filtered) == 0) {
    cat("  No genes passed the filter for:", sample_name, "\n")
    next
  }
  
  # Create ranked gene list
  geneList <- deg_filtered$AVG.Log2.Ratio
  names(geneList) <- deg_filtered$Genes
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Run GSEA for each gene set
  for (set_name in names(gene_sets_config)) {
    config <- gene_sets_config[[set_name]]
    cat("  Running GSEA for:", set_name, "\n")
    
    gene_sets <- msigdbr(
      species = "Homo sapiens", 
      category = config$category, 
      subcategory = config$subcategory
    )
    
    geneset <- gene_sets %>%
      select(gs_name, gene_symbol) %>%
      distinct()
    
    gsea_result <- GSEA(
      geneList,
      TERM2GENE = geneset,
      minGSSize = 5,
      pvalueCutoff = 0.05,
      verbose = FALSE,
      eps = 0
    )
    
    result_df <- gsea_result@result
    result_df$ID <- NULL
    
    write.csv(result_df, paste0(sample_name, "_GSEA_", set_name, ".csv"), row.names = FALSE)
    
    # Plot top pathways
    top_pathways <- result_df %>%
      mutate(ordering = abs(NES)) %>%
      arrange(desc(ordering)) %>%
      group_by(sign(NES)) %>%
      slice_max(order_by = ordering, n = 20, with_ties = FALSE) %>%
      ungroup()
    
    pdf(paste0("GSEA_", sample_name, "_", set_name, ".pdf"), width = 10, height = 6)
    print(
      ggplot(top_pathways, aes(x = NES, y = fct_reorder(Description, NES), fill = p.adjust)) + 
        geom_bar(stat = 'identity') +
        scale_fill_continuous(low = '#440154FF', high = '#1F968BFF', guide = guide_colorbar(reverse = TRUE)) + 
        theme_minimal() + 
        labs(y = NULL, x = "NES", title = paste("GSEA -", sample_name, "-", set_name))
    )
    dev.off()
  }
}
