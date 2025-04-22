library(Seurat)
library(ggplot2)
library(dplyr)

setwd("..")

# Load your Seurat object
seurat_obj <- readRDS("../prenatal_final.rds")

# Load your differential expression results (replace with actual dataframe)
de_results <- read.csv("../proteomics.csv")

# Define thresholds
log2FC_threshold <- 0.25
qvalue_threshold <- 0.05

# Identify upregulated and downregulated genes
upregulated_genes <- de_results %>%
  filter(AVG.Log2.Ratio >= log2FC_threshold & Qvalue <= qvalue_threshold) %>%
  pull(Genes)

downregulated_genes <- de_results %>%
  filter(AVG.Log2.Ratio <= -log2FC_threshold & Qvalue <= qvalue_threshold) %>%
  pull(Genes)

# Check which genes are present in the Seurat object
valid_upregulated_genes <- upregulated_genes[upregulated_genes %in% rownames(seurat_obj@assays$RNA@data)]
valid_downregulated_genes <- downregulated_genes[downregulated_genes %in% rownames(seurat_obj@assays$RNA@data)]

library(Seurat)
library(ggplot2)

# Ensure unique factor levels
seurat_obj$cellsubtype <- factor(seurat_obj$cellsubtype, levels = unique(seurat_obj$cellsubtype))

# Remove genes found in multiple assays
valid_upregulated_genes <- valid_upregulated_genes[valid_upregulated_genes %in% rownames(seurat_obj)]

# Plot and save DotPlot for upregulated genes
if (length(valid_upregulated_genes) > 0) {
  pdf("DotPlot_Upregulated_Genes.pdf", width = 17, height = 6)
  dotplot_up <- DotPlot(seurat_obj, features = valid_upregulated_genes, cols="RdBu", group.by = "cellsubtype") +
    theme_minimal() +
    ggtitle("DotPlot of Upregulated Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(dotplot_up)
  dev.off()
}

# Remove genes found in multiple assays for downregulated genes
valid_downregulated_genes <- valid_downregulated_genes[valid_downregulated_genes %in% rownames(seurat_obj)]

# Plot and save DotPlot for downregulated genes
if (length(valid_downregulated_genes) > 0) {
  pdf("DotPlot_Downregulated_Genes.pdf", width = 15, height = 6)
  dotplot_down <- DotPlot(seurat_obj, features = valid_downregulated_genes, cols="RdBu", group.by = "cellsubtype") +
    theme_minimal() +
    ggtitle("DotPlot of Downregulated Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(dotplot_down)
  dev.off()
}

