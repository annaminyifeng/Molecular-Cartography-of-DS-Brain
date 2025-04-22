```{r}
library(fgsea)
library(readxl)
library(data.table)
library(ggplot2)
library(Seurat)
library(tidyverse)
```

```{r}
input_file <-"../6mo_Ts65Dn_TE_celltypes.rds"
te_obj <- readRDS(input_file)
te_obj <- ScaleData(te_obj, assay = "RNA")
DefaultAssay(te_obj) <- "RNA"


data <- read_excel("../TE_families_modified.xlsx")
grouped_data <- aggregate(all_features ~ `Classification (Final)`, data, function(x) list(x))
te_list <- with(grouped_data, setNames(lapply(all_features, unlist), `Classification (Final)`))
te_list <- te_list[names(te_list) != "NA"]
print(te_list)

all_features <- Features(te_obj)
solo_TE_features <- all_features[grepl("^SoloTE-", all_features)]
my_te_set <- te_list


```

```{r}
group.by <- "celltype"
celltypes <- unique(te_obj@meta.data[[group.by]])
Idents(te_obj) <- group.by
gsea_results <- list()
de_results <- list()

# Loop through cell types for MAST
for(cell_type in celltypes){
  sub_te_obj <- subset(te_obj, idents = cell_type)
  
  Idents(sub_te_obj) <- "Condition"
  DEG <- FindMarkers(
    sub_te_obj,
    ident.1 = "DS",
    ident.2 = "CTL",
    test.use = "MAST",
    features = solo_TE_features
  )
  
  DEG$gene <- rownames(DEG)
  DEG$celltype <- cell_type
  DEG <- DEG %>% arrange(desc(avg_log2FC))
  de_results[[cell_type]] <- DEG
}

# Loop through cell types for GSEA & Volcano Plot
for(cell_type in celltypes){
  FC <- de_results[[cell_type]]$avg_log2FC
  de_results[[cell_type]]$gene_name <- sub("^SoloTE-", "", de_results[[cell_type]]$gene)
  names(FC) <- de_results[[cell_type]]$gene_name
  
  # Add small noise to resolve ties
  set.seed(123) # For reproducibility
  FC <- FC + rnorm(length(FC), mean = 0, sd = 1e-5)
  
  # Run fgsea
  subset_result <- fgseaMultilevel(
    pathways = my_te_set,
    stats = FC,
    minSize = 1,
    eps = 0,
    maxSize = 500,
    nPermSimple = 10000
  )
  
  subset_result$celltype <- cell_type
  gsea_results[[cell_type]] <- subset_result
  
  # Save DEG results as CSV
  safe_comparison_name <- gsub("[^a-zA-Z0-9]", "_", cell_type)
  write.csv(de_results[[cell_type]], file.path(output_dir, paste0(safe_comparison_name, ".csv")), row.names = FALSE)
  
  # Prepare data for volcano plot
  df <- de_results[[cell_type]]
  df$gene_label <- rownames(df)
  df$group <- with(df, ifelse(avg_log2FC > 0.25 & p_val_adj < 0.05, "up_sig", 
                              ifelse(avg_log2FC > 0.25 & p_val_adj >= 0.05, "up_non_sig",
                                     ifelse(avg_log2FC <= -0.25 & p_val_adj < 0.05, "down_sig", "down_non_sig"))))
  
  # Generate Volcano Plot
  volcano <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) +
    geom_point(alpha = 0.8) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle(cell_type) +
    ylab("-log10(p_val_adj)") +
    xlab("log2 Fold Change") +
    scale_color_manual(values = c("up_sig" = "#1F968BFF", "down_sig" = "#440154FF", "up_non_sig" = "gray", "down_non_sig" = "gray")) +
    geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.5, linetype = "dashed", color = "gray")
  
  # Save Volcano Plot
  ggsave(file.path(output_dir, paste0(safe_comparison_name, "_volcano.pdf")), volcano, width = 6, height = 4)}

```



```{r}

library(dplyr)
library(tidyr)

# Assuming 'df_list' is your list of data frames and 'cell_types' is a vector of cell type names
combined_df <- bind_rows(gsea_results, .id = "celltype") 


# Creating the heat map with significance markers
plot1 <- ggplot(combined_df, aes(x = celltype, y = pathway, fill = NES)) + 
  geom_tile() +  # This creates the tiles for the heat map
  scale_fill_gradient2(low = "#440154FF", high = "#1F968BFF", mid = "#f7f7f7", midpoint = 0, name = "NES") +
  geom_text(aes(label = ifelse(padj < 0.05, "*", "")), color = "black", size = 5) +  # Add asterisks for significant values
  theme_minimal() +  # Minimal theme for the plot
  labs(fill = "NES", x = "Cell Type", y = "Pathway")  + # Labeling the axes and legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(face = "bold"),
        legend.title.align = 0) 
# Save the first plot
ggsave("../broad_celltype.pdf", plot1, width = 10, height = 4)


```

```{r}
#Plot retrotransposons from interneuron celltype
library(dplyr)
library(ggplot2)

# Read data from CSV
file_path <- "../Inhibitory_Neurons.csv"
df <- read.csv(file_path)

# List of specific genes to highlight
genes_to_highlight <- c("LTR86B1", "RLTR13D6", "RLTR13D",
                        "MLT1F-int", "RLTR4-MM-int", "LTR85a", "IAPEY5-LTR", "MamGypLTR3")
# Prepare data
df$group <- with(df, ifelse(avg_log2FC > 0.25 & p_val_adj < 0.05, "up_sig", 
                            ifelse(avg_log2FC > 0.25 & p_val_adj >= 0.05, "up_non_sig",
                                   ifelse(avg_log2FC <= -0.25 & p_val_adj < 0.05, "down_sig", "down_non_sig"))))

# Highlight specific genes
df$highlight <- ifelse(df$gene_name %in% genes_to_highlight, "highlight", "no_highlight")

# Filter highlighted data
highlight_df <- subset(df, highlight == "highlight")

# Plotting
volcano <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = c("up_sig" = "#1F968BFF", "down_sig" = "#440154FF", 
                                "up_non_sig" = "gray", "down_non_sig" = "gray")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Inhibitory Neurons Volcano Plot") +
  ylab("-log10(p_val_adj)") +
  xlab("log2 Fold Change") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.5, linetype = "dashed", color = "gray") +
  
  # Add gene labels and connecting lines
  geom_text(data = highlight_df, aes(label = gene_name),
            color = "black", vjust = -0.5, hjust = 0.5, size = 3) +
  
  geom_segment(data = highlight_df,
               aes(x = avg_log2FC, xend = avg_log2FC,
                   y = -log10(p_val_adj) + 0.5, yend = -log10(p_val_adj)),
               color = "black", linetype = "dotted")

# Save Volcano Plot
ggsave("../inhibitory_neurons_volcano_retrotransposons.pdf", 
       volcano, width = 5, height = 5)

```
