
"Deseq2: Eggnog PCoA"

#Plot PCA: all locations
counts = read.csv('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Publication/tables/Fields shotgun/merged_COGs_file.csv',row.names = 1)
meta = read.csv('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Publication/tables/Fields shotgun/Meta_all_samples.csv',row.names = 1, header=TRUE)
colnames(meta)[1] = "condition"

# Check if all rownames of meta are in the column names of counts
if(all(rownames(meta) %in% colnames(counts))) {
  # Reorder columns of counts to match the order of rownames in meta
  counts <- counts[, rownames(meta)]
} else {
  stop("Not all row names in meta are column names in counts.")
}

library(DESeq2)
library(ggplot2)
library(dplyr)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ condition)

vst_data <- vst(dds, blind=TRUE)
pca_data <- plotPCA(vst_data, intgroup=c("condition"), returnData=TRUE)

#PLOT AVERAGE
# Plotting with ggplot2
library(ggpubr)

ggplot(pca_data, aes(x=PC1, y=PC2, color=condition, label=condition)) +
  geom_point(size=4) +
  theme_pubr() +
  scale_fill_manual(values =  c("Bacterioplankton" = "#b3669e", "Daphnia" = "#98984d")) +
  labs(x=paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
       y=paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance")) +
  theme(legend.position="bottom")+
  ggtitle("PCA for COGs composition: Microbiome vs Bacterioplankton")

# Load the required library
library(stats)

# Fit a multivariate linear model (MANOVA)
manova_model <- manova(cbind(PC1, PC2) ~ condition, data = pca_data)

# Perform the MANOVA test
manova_test <- summary(manova_model) #2.2e-16 ***

