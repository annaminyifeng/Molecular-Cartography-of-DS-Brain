setwd("..")

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
  C5_GO_CC = list(category = "C5", subcategory = "GO:CC"),
  C5_GO_MF = list(category = "C5", subcategory = "GO:MF"),
  C2_REACTOME = list(category = "C2", subcategory = "CP:REACTOME"),
  C2_CP = list(category = "C2", subcategory = "CP")
)

# Get all CSV files in the folder
csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file
for (file in csv_files) {
  # Extract a short name from the file (without extension)
  sample_name <- tools::file_path_sans_ext(basename(file))
  cat("Processing:", sample_name, "\n")
  
  # Load DGE results
  deg <- read.csv(file)
  geneList <- deg$avg_log2FC
  names(geneList) <- deg$gene
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Loop through gene set configurations
  for (set_name in names(gene_sets_config)) {
    config <- gene_sets_config[[set_name]]
    cat("  Running GSEA for:", set_name, "\n")
    
    # Get gene sets
    gene_sets <- msigdbr(
      species = "Homo sapiens", 
      category = config$category, 
      subcategory = config$subcategory
    )
    
    geneset <- gene_sets %>%
      dplyr::select(gs_name, gene_symbol) %>%
      distinct()
    
    # Run GSEA
    gsea_result <- GSEA(
      geneList,
      TERM2GENE = geneset,
      minGSSize = 5,
      pvalueCutoff = 0.05,
      verbose = FALSE,
      eps = 0
    )
    
    # Save results
    result_df <- gsea_result@result
    result_df$ID <- NULL
    write.csv(result_df, paste0(sample_name, "_GSEA_", set_name, ".csv"), row.names = FALSE)
    
    # Plot top pathways
    top_pathways <- result_df %>%
      mutate(ordering = abs(NES)) %>%
      arrange(desc(ordering)) %>%
      group_by(sign(NES)) %>%
      slice_max(order_by = ordering, n = 20)
    
    pdf(paste0("GSEA_", sample_name, "_", set_name, ".pdf"), width = 10, height = 6)
    print(
      ggplot(top_pathways, aes(x = NES, y = fct_reorder(Description, NES), fill = p.adjust)) + 
        geom_bar(stat = 'identity', orientation = 'y') +
        scale_fill_continuous(low = '#440154FF', high = '#1F968BFF', guide = guide_colorbar(reverse = TRUE)) + 
        theme_minimal() + 
        labs(y = NULL, title = paste("GSEA -", sample_name, "-", set_name))
    )
    dev.off()
  }
}
