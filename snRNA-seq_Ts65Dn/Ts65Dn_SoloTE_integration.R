library(Seurat)
library(future)
options(future.globals.maxSize = 1024^3)  
future::plan("multisession", workers = 4)

setwd("../")

seurat_list <- readRDS("combined_list.rds")

# Merge the objects instead of using FindIntegrationAnchors
seurat_merged <- merge(seurat_list[[1]], y = seurat_list[-1])

# Normalize, find variable features, and scale
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged)
seurat_merged <- ScaleData(seurat_merged)

# Run PCA before Harmony
seurat_merged <- RunPCA(seurat_merged)

# Run Harmony
seurat_merged <- RunHarmony(seurat_merged, group.by.vars = "orig.ident")

# Run UMAP and clustering on Harmony embeddings
seurat_merged <- RunUMAP(seurat_merged, reduction = "harmony", dims = 1:20)
seurat_merged <- FindNeighbors(seurat_merged, reduction = "harmony", dims = 1:20)
seurat_merged <- FindClusters(seurat_merged)

DefaultAssay(seurat_merged) <- "RNA"

saveRDS(seurat_merged, "SoloTE_Ts65Dn.rds")
