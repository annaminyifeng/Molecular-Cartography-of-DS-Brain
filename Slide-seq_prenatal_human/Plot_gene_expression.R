# Define colors for plotting including gray for 'singlet'
cell_color <- c(
  "NSC_1"= "#e1a9f6",
  "NSC_2"= "#ead7f5",
  "NSC_3"= "#d1bbf9",
  "NSC_4"= "#cbd9f5",
  "oRG"= "#ffcae5",
  "vRG"= "#c6c7ff",
  "Cycling_progenitors"= "#f9abe7",
  "IP"= "#ef87a9",
  "Excitatory_neuroblast" = "#d24459",
  "IP_IN"= "#bbeeff",
  "TAP_1" ="#7eedff",
  "TAP_2" ="#88bdf2",
  "Inhibitory_neuroblast"= "#6699ed",
  "Mig_L2-3_ExN"="#E05A88",
  "Mig_L2-4_ExN"= "#D12D66",
  "L2-4_ExN"="#be49dd",
  "Layer_5-6_ExN"="#5732a8",
  "Piriform_1"= "#FF4FA7",
  "Piriform_2"= "#FF1D8E",
  "Mig_Layer_2-3_ExN"= "#E05A88",
  "Mig_Layer_2-4_ExN"="#D12D66",
  "Layer_1-6a_ExN"= "#BE1B88",
  "Layer_2-6a_ExN"= "#9b00ae",
  "Layer_1_ExN"=  "#F57373" ,
  "Layer_2-3_ExN" = "#E55EA8",
  "Layer_4-5_ExN"= "#D449DD", 
  "Layer_4-6_ExN"= "#863BBC",
  "Layer_2-6a_ExN"= "#963EC3",
  "Layer_5-6a_ExN"= "#7738B6",
  "Layer_2-4_ExN"= "#CB37A4",
  "Layer_2-5_ExN"=  "#a049dd",
  "Layer_4-6a_ExN"=  "#8A0DAD",
  "Layer_6b_ExN"= "#580A70",
  "Subplate" = "#5723E7",
  "Interneuron"="#3366ff" ,
  "Interneuron_1" = "#3366ff" ,
  "Interneuron_2" = "#31AED3",
  "Interneuron_3" =  "#0096FF",
  "Interneuron_4" ="#00dbff",
  "Interneuron_S1"="#0b4ade",
  "Interneuron_S2"= "#032dff",
  "Interneuron_M1"= "#2c7fb8",
  "Interneuron_M2"="#1d547a",
  "Myeloid"="#d9feeb",
  "Microglia_1"="#75c8be",
  "Immune"="#75c8be",
  "Microglia_2"="#008893",
  "Astrocyte_progenitors"= "#c7e9b4",
  "Astrocyte_1"= "darkolivegreen1",
  "Astrocyte_2"="#79d33d",
  "Astrocyte_3"="#009000", 
  "Astrocyte"="#009000",
  "Astrocyte_4"="#02783D",
  "Fibroblast"= "#ffbf7f",
  "Fibroblast_1"= "#ffbf7f",
  "Vascular"="#FCE205",
  "Smooth_muscle_1"="#FCE205",
  "Smooth_muscle_2"= "#ffce00",
  "Ependymal"="#ffbbb1",
  "OPC"= "#ff6f4b",
  "OPC_1"= "#ff6f4b",
  "OPC_2"= "#ff4000",
  "OL"="#fa1b27",
  # Add gray for 'singlet'
  "singlet" = "gray"
)


###
combined<-readRDS("../Slide-seq.rds")
setwd("../")

library(Seurat)
library(ggplot2)

pdf('CTL_gene_expression.pdf', width=8, height=2)

# Create plots for each feature
plots <- lapply(c("FOS", "EGR1", "IER2", "DNMT3A", "HMGB2"), function(feature) {
  p <- SpatialFeaturePlot(object = brain, 
                          features = feature, 
                          images = "slide_seq.3", 
                          min.cutoff = "q10", 
                          max.cutoff = "q90", 
                          alpha = c(0.1, 1), 
                          stroke = NA) + NoLegend()
  print(p)  # <-- You need to print it
  return(p)
})

# Close PDF device
dev.off()
