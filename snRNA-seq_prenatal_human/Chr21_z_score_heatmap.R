library("dendsort")
library("seriation")
library("gplots")
library("RColorBrewer")
library(devtools)
library("ComplexHeatmap")
library(colorRamp2)

NucSeq.atac <- readRDS("../prenatal_final.rds")
DefaultAssay(NucSeq.atac) <- 'SCT'
#reset idents
Idents(NucSeq.atac) <- factor(as.character(NucSeq.atac$cellsubtype), levels=unique(as.character(NucSeq.atac$cellsubtype))[order(unique(as.character(NucSeq.atac$cellsubtype)))])


genes <- c("ABCG1", "ADAMTS5", "ADARB1", "AGPAT3", "APP", "B3GALT5", "BACH1", "BTG3", "CBR1", "CCT8", 
           "CHAF1B", "COL18A1", "COL5A2", "COL6A1", "CRYZL1", "CSTB", "CXADR", "DSCAM", "DYRK1A", 
           "ERG", "ETS2", "FTCD", "GABPA", "GART", "GRIK1", "HSF2BP", "IFNAR1", "IFNAR2", "IFNGR2", 
           "IGSF5", "IL10RB", "KCNE1", "KCNJ6", "LSS", "MCM3", "MX1", "NCAM2", "NDUFV3", "NRIP1", 
           "OLIG2", "PCBP3", "PCNT", "PCP4", "PDE9A", "PDXK", "PFKL", "PKNOX1", "RBM11", "RCAN1", "RUNX1", "S100B", 
           "SH3BGR", "SLC19A1", "SLC37A1", "SLC5A3", "SOD1", "SON", "SYNJ1", "TCP10L", 
           "TIAM1", "TRPC7", "TTC3", "U2AF1", "UBE2G2", "USP16", "USP25", "WDR4")
gene_anno_list <- genes

# Define a function to calculate average expression for each gene across each cell subtype
calculate_average_expression <- function(seurat_object, genes) {
  # Initialize an empty list to store average expression matrices for each cell subtype
  average_expression_list <- list()
  
  # Loop through each unique cell subtype
  for (subtype in unique(seurat_object$cellsubtype)) {
    # Subset data for the current subtype
    cur_seurat_subtype <- subset(seurat_object, cellsubtype == subtype)
    
    # Calculate pseudo-expression for the filtered genes
    pseudoexpression_matrix <- AverageExpression(cur_seurat_subtype, assay = 'SCT', slot = 'counts', features = genes)
    
    # Store the pseudo-expression matrix in the list with the subtype as the key
    average_expression_list[[subtype]] <- pseudoexpression_matrix
  }
  
  return(average_expression_list)
}

# Calculate average expression matrices for each cell subtype
average_expression_list <- calculate_average_expression(NucSeq.atac, genes)

# Combine average expression matrices horizontally
combined_matrix <- Reduce(cbind, lapply(average_expression_list, function(mat) as.matrix(as.data.frame(mat))))

# Change column names to match cell subtype names
colnames(combined_matrix) <- names(average_expression_list)

# Calculate z-scores for combined matrix
matrix_z <- apply(combined_matrix, 1, zScore) %>% t()
matrix_z <- matrix_z[, order(colnames(matrix_z))]

# Set up row annotations
#gene_anno_list <- genes %>% top_n(1, wt = avg_log2FC) %>% .$gene %>% unique
#gene_anno_list <- c(gene_anno_list, more_gene_list) %>% unique
#gene_anno_list <- gene_anno_list[gene_anno_list %in% rownames(matrix_z)]
#ha <- rowAnnotation(foo = anno_mark(at = unlist(lapply(gene_anno_list, function(gene) {
which(rownames(matrix_z) == gene)
#})), labels = gene_anno_list))

# Hierarchical clustering
row_dend <- dendsort(hclust(dist(matrix_z)))
col_dend <- dendsort(hclust(dist(t(matrix_z))))

pdf('heatmap_z_score_cell_subtypes.pdf', width = 8, height = 10)
ComplexHeatmap::Heatmap(
  matrix_z, show_column_names = FALSE, show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 8),  # Adjust font size of column names
  cluster_columns = FALSE,
  bottom_annotation = HeatmapAnnotation(
    text = anno_text(colnames(matrix_z), rot = 90, location = unit(1, "npc"), just = "right"),
    annotation_height = max_text_width(colnames(matrix_z))
  ),
  right_annotation = ha,
  col = colorRamp2(c(min(matrix_z, na.rm = TRUE), 0, max(matrix_z, na.rm = TRUE)), c("#440154FF", "white", "#1F968BFF"))
)
dev.off()

