library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(grDevices)

# Load the data
df <- read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/KEGG-KOs/KO-pathways-abundances.csv")
df = df %>%
  dplyr::filter(Polymer != "PS")

df$Order_for_heatmap <- as.numeric(as.character(df$Order_for_heatmap))

# Now arrange the data frame by 'Order_for_heatmap'
df_ordered <- df %>%
  arrange(Order_for_heatmap) %>%
  mutate(Enzyme = factor(Enzyme, levels = unique(Enzyme)))

df_ordered$Enzyme = paste(df_ordered$KO.number, ":",df_ordered$Enzyme, sep="")
enzyme_names = df_ordered$Enzyme
pathway=df_ordered$Pathway

# Convert to matrix
data_matrix <- as.matrix(cbind(df_ordered$High_MPs,df_ordered$Low_MPs))
data_matrix = sqrt(data_matrix)
colnames(data_matrix) <- c("High MPs", "Low MPs")
rownames(data_matrix) <- enzyme_names

# Create a color vector for the conditions
pathway_condtions = df_ordered %>% 
  select(Pathway)
  
rownames(pathway_condtions) = df_ordered$Enzyme
annotation_colors <- list(Pathway = c("Styrene degradation" = "darkgreen", "Polycyclic aromatic hydrocarbon degradation" = "yellow"))
# Generate heatmap

pheatmap(data_matrix, cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("white", "darkred"))(100), 
         display_numbers = TRUE, # Optionally display numbers
         main = "Plastic polymer Pathway Enzymes", # Optional title
         fontsize_row = 10,
         number_color="black",
         border_color = "black",
         annotation_row = pathway_condtions,
         angle_col = 90)







