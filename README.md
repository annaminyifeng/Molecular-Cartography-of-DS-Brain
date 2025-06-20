# Molecular Cartography of the Human and Mouse Down Syndrome Brain

This GitHub repository contains the custom code for microenvironment analysis, aimed at quantifying changes in local cell type compositions between Down Syndrome (DS) and euploid samples. This repository also includes all relevant codes for the snRNA-seq, Slide-seq, and MERFISH data analyses from the associated publication:

**"Molecular Cartography of the Human and Mouse Down Syndrome Brain"**.

# Abstract:

Down syndrome (DS, or Trisomy 21) is one of the most common genetic causes of intellectual disability. DS results in both abnormal neurodevelopment and accelerated neurodegeneration, but the molecular mechanisms underlying abnormal corticogenesis are incompletely understood. To gain molecular insight into the prenatal neurobiology of DS, we performed single-nucleus sequencing, spatial transcriptomics, and proteomics on mid-gestational prenatal human cortex. We captured altered expression dynamics of lineage commitment genes and pronounced de-repression of transposable elements in DS neural progenitor cells, which suggest changes to the fate and functionality of neuronal and glial cells. Given the importance of linking human and model system pathobiology, we also performed highly multiplexed RNA in situ spatial transcriptomics on a well-established trisomic mouse model (Ts65Dn) to study the cellular landscape of the trisomic brain during early development and maturation. We profiled the spatial transcriptome of > 240,000 cells in the mouse brain and identified trisomy-associated gene expression patterns in the molecular control of neurogenesis and gliogenesis. Together, our study provides an extensive resource for understanding of the complex multicellular processes underlying DS neurodevelopment.

## Reproducing the Microenvironment Analysis

To reproduce the microenvironment analysis, we have uploaded a small MERFISH dataset of the P0 CC ("P0_cc_final.rds") into the `Microenv_analysis_tool` folder. The expected output will match our results from **Figure 6d** in the paper.

Please refer to the **Microenvironment Analysis Methods** section for detailed instructions on how to use the tool.

## Package Requirements

This code is compataible with **4.2.1** and was tested in a Code Ocean container. To run the analysis, you will need the following R packages: 

| Package           | Version   |
|------------------|-----------|
| tidyverse         | 2.0.0     |
| ggplot2           | 3.5.2     |
| dplyr             | 1.1.4     |
| tidyr             | 1.3.1     |
| reshape2          | 1.4.4     |
| future            | 1.49.0    |
| purrr             | 1.0.4     |
| nlme              | 3.1-168   |
| Seurat            | 5.1.0     |
| ggrepel           | 0.9.6     |
| nloptr            | latest    |
| Rfast             | latest    |
| sctransform       | latest    |
| SeuratObject      | latest    |
| Matrix            | 1.6.5     |
| stringi           | 1.4.6     |
| devtools          | latest    |
| SingleCellExperiment | latest via Bioconductor |
