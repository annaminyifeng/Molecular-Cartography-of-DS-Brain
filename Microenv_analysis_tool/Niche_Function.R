library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
library(tidyverse)

BuildNicheAssay <- function(object, fov, group.by, assay = "niche", neighbors.k = 20, 
                            niches.k = 4) 
{
  
  coords <- GetTissueCoordinates(object[[fov]], which = "centroids")
  cells <- coords$cell
  rownames(coords) <- cells
  coords <- as.matrix(coords[, c("x", "y")])
  neighbors <- FindNeighbors(coords, k.param = neighbors.k)
  
  # neighbors$nn <- neighbors$nn[Cells(object), Cells(object)]
  
  ct.mtx <- matrix(data = 0, nrow = length(cells), ncol = length(unlist(unique(object[[group.by]]))))
  rownames(ct.mtx) <- cells
  colnames(ct.mtx) <- unique(unlist(object[[group.by]]))
  cts <- object[[group.by]]
  for (i in 1:length(cells)) {
    ct <- as.character(cts[cells[[i]], ])
    ct.mtx[cells[[i]], ct] <- 1
  }
  sum.mtx <- as.matrix(neighbors$nn %*% ct.mtx)
  niche.assay <- CreateAssayObject(counts = t(sum.mtx))
  object[[assay]] <- niche.assay
  DefaultAssay(object) <- assay
  object <- ScaleData(object)
  results <- kmeans(x = t(object[[assay]]@scale.data), centers = niches.k, 
                    nstart = 30)
  object$niches <- results[["cluster"]]
  return(object)
}