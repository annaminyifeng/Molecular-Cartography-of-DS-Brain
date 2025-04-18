source("../Niche_Function.R")

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(future)
library(purrr)
library(nlme)

# > packageVersion("tidyverse")
# [1] ‘2.0.0’
# > packageVersion("ggplot2")
# [1] ‘3.5.1’
# > packageVersion("dplyr")
# [1] ‘1.1.4’
# > packageVersion("tidyr")
# [1] ‘1.3.1’
# > packageVersion("reshape2")
# [1] ‘1.4.4’
# > packageVersion("future")
# [1] ‘1.34.0’
# > packageVersion("purrr")
# [1] ‘1.0.2’
# > packageVersion("nlme")
# [1] ‘3.1.164’

plan("multisession", workers = 8)

setwd("../Microenv")
saveDir <- "../Microenv/"
input <- "../P0_cc_final.rds"
cell.group = "celltype"  
comparason.group = "Condition" #There must be two conditions in this column
replicate.column = "Replicate" #Column that identifies replicate
baseline_condition <- "WT" #Indicate your baseline
neighbors.k = 100 #Can change to what number you desire to use for microenv analysis
num_niches = 6 #Check Seurat documentation for Niches
skip_threshold_percentage <- 90  #You can adjust this percentage as needed. Currently, it means that if a query cell type is absent around a particular center cell type 10% of the time, that comparison will be skipped.

obj <- readRDS(input)
all_conditions <- unique(obj[[comparason.group]])
disease_condition <- all_conditions[all_conditions != baseline_condition]
neighbor_df <- list()
if (!dir.exists(saveDir)) {dir.create(saveDir, recursive = TRUE)}
fovs <- names(obj@images)
for(fov in fovs){
  obj <- BuildNicheAssay(object = obj, fov = fov, group.by = cell.group, niches.k = num_niches, neighbors.k = 100)
  DefaultAssay(obj) <- 'niche'
  table <- as.data.frame(obj@assays$niche@counts)
  table <- t(table)
  neighbor_df[[fov]] <- table
  print(paste0("Finished Processing: ", fov))
}

combined_df <- do.call(rbind, neighbor_df)
combined_df <- as.data.frame(combined_df)
celltype.vector <- colnames(combined_df)
combined_df$cell_id <- rownames(combined_df)
metadata_df <- obj@meta.data
metadata_df$cell_id <- rownames(metadata_df)
merged_df <- merge(combined_df, metadata_df, by = "cell_id")
write.csv(merged_df, file = paste0(saveDir, "combined_df.csv"))

data <- merged_df
data$Sample <- data[[replicate.column]]
data$celltype <- data[[cell.group]]
data$celltype <- as.factor(gsub("_", "a", data$celltype))
colnames(data) <- gsub("//.", "-", colnames(data))

data$Condition <- factor(data$Condition, levels = c(baseline_condition, disease_condition))
unique_cell_types <- unique(data$celltype)

results_df <- tibble(
  CentroidCellType = character(),
  QueryCellType = character(),
  TestType = character(),
  PValue_model = numeric(),
  PValue_wilcox = numeric(),
  MeanDifference = numeric(),
  Condition1Median = numeric(),
  Condition2Median = numeric()
)

for (centroid_cell_type in unique_cell_types) {
  for (query_cell_type in unique_cell_types) {
    if (!query_cell_type %in% names(data)) { next }
    
    query_data <- data %>%
      dplyr::filter(celltype == centroid_cell_type) %>%
      dplyr::select(cell_id, Condition, Sample, all_of(query_cell_type)) %>%
      drop_na()
    # Skip if not enough data
    if (nrow(query_data) < 2 || length(unique(query_data$Condition)) < 2) { next }
    # Calculate the percentage of zero counts for the query cell type
    zero_count_percentage <- sum(query_data[[query_cell_type]] == 0) / nrow(query_data) * 100
    # Filtration step: skip if more than the threshold percentage of counts are zero
    if (zero_count_percentage > skip_threshold_percentage) { next }
    # Additional filtration: ensure there are at least two non-zero observations for each condition
    condition_counts <- query_data %>% group_by(Condition) %>% 
      summarise(non_zero_count = sum(.data[[query_cell_type]] > 0))
    if (any(condition_counts$non_zero_count < 2)) { next }
    # Additional check: ensure we have enough data to fit a mixed-effects model
    if (length(unique(query_data$Sample)) < 2) { next }  # We need at least 2 groups for random effects
    
    
    query_data_control <- query_data[query_data$Condition == baseline_condition,]
    query_data_disease <- query_data[query_data$Condition == disease_condition,]
    wilcox_result <- wilcox.test(query_data_control[[query_cell_type]], query_data_disease[[query_cell_type]])
    p_value_wilcox <- wilcox_result$p.value
    # Prepare and fit the full linear mixed-effects model using ML instead of REML
    full_formula <- as.formula(paste(query_cell_type, "~ Condition", sep = ""))
    full_model <- tryCatch({
      lme(full_formula, random = ~1 | Sample, data = query_data, method = "ML")
    }, error = function(e) return(NULL))
    
    # Fit the reduced model (without the Condition effect) using ML
    reduced_formula <- as.formula(paste(query_cell_type, "~ 1", sep = ""))
    reduced_model <- tryCatch({
      lme(reduced_formula, random = ~1 | Sample, data = query_data, method = "ML")
    }, error = function(e) return(NULL))
    
    # Only continue if both models were successfully fitted
    if (is.null(full_model) || is.null(reduced_model)) { next }
    
    # Perform the likelihood ratio test
    lrt_result <- anova(full_model, reduced_model)
    p_value <- lrt_result$`p-value`[2] # Get the p-value for the condition effect
    
    condition1_median <- median(query_data %>% filter(Condition == levels(query_data$Condition)[1]) %>% pull(query_cell_type))
    condition2_median <- median(query_data %>% filter(Condition == levels(query_data$Condition)[2]) %>% pull(query_cell_type))
    median_diff <- condition1_median - condition2_median
    
    # Extract fixed effects from the full model
    fixed_effects <- fixef(full_model)
    
    # Assuming the second fixed effect is the one of interest
    condition_effect <- fixed_effects[2]
    
    # Add results to the dataframe
    results_df <- rbind(results_df, tibble(
      CentroidCellType = centroid_cell_type,
      QueryCellType = query_cell_type,
      TestType = "Mixed-Effects Model",
      PValue_model = p_value,
      PValue_condition = p_value_wilcox,
      MeanDifference = condition_effect,  # or median_diff if more appropriate
      Condition1Median = condition1_median,
      Condition2Median = condition2_median
    ))
  }
}

results_df$PValueAdjusted_BH <- p.adjust(results_df$PValue_condition, method = "BH")

results_df <- results_df %>%
  mutate(Comparison_Significance = case_when(PValueAdjusted_BH < 0.05 ~ "*", TRUE ~ ""))
results_df <- results_df %>%
  mutate(Model_Significance = case_when(PValue_model < 0.05 ~ "*", TRUE ~ ""))

results_df <- results_df %>%
  mutate(Significance = case_when(Comparison_Significance == "*" & Model_Significance == "*" ~ "*", TRUE ~ ""))

write_csv(results_df, paste0(saveDir, "results_df.csv"))


pdf("Microenv_CC.pdf")
ggplot(results_df, aes(x = QueryCellType, y = CentroidCellType, fill = MeanDifference)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Significance), vjust = 0.8, color = "black", size = 6) +
  scale_fill_gradient2(low = "#440154FF", high = "#1F968BFF", mid = "white", midpoint = 0, 
                       limits = c(-12.65, 16), name = "Log2 Fold Change") +
  theme_minimal() +
  labs(title = paste0("MicroEnv Change in ", disease_condition , ' vs ', baseline_condition),
       x = "Surrounding Cell Type", y = "Center Cell Type", fill = "Log2FC") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(face = "bold"),
        legend.title.align = 0)
dev.off()
