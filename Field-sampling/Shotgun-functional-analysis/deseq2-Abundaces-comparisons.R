"
Deseq2
"

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)  # For theme_pubr()

#Load data and packages:
counts = read.csv('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/14_All_together_datasets/COGs/COGs_all_samples_raw_counts.tsv',row.names = 1,sep="\t")
meta = read.csv('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/metadata.csv',row.names = 1, header=TRUE)
#colnames(meta)[2] = "condition"
meta <- meta[colnames(counts),]

##remove empty COGs
counts_tidy <- counts %>%
  dplyr::filter(rowSums(.) > 0) #proportion

#take proportions
counts_tidy = data.frame(t(counts_tidy))
counts_total = rowSums(counts_tidy)
proportions_counts = sweep(counts_tidy, 1, counts_total, "/")

bray <- proportions_counts %>%
  vegdist(method = "bray")

PCoA= capscale(bray~1)
pcoa_scores <- scores(PCoA, display = "sites")
pcoa_df <- as.data.frame(pcoa_scores)
pcoa_df$Type = meta$Type
pcoa_df$Category = meta$Category

# Plot using ggplot
pca_plot_customized <- ggplot(pcoa_df, aes(x = MDS1, y = MDS2, color = Type, shape = Category)) +  # Adjust axis names as needed
  geom_point(size = 5) +  # Adjust size as needed
  scale_color_manual(values = c("Daphnia" = "#98984d", "Bacterioplankton" = "#b3669e")) +
  scale_shape_manual(values= c("High"= 10, "Low"=16))+
  theme_pubr() +  # Apply a publication-ready theme
  labs(x = paste("PCoA 1:",round(PCoA$CA$eig[1]/sum(PCoA$CA$eig),2)*100, "%"), y = paste("PCoA 2:",round(PCoA$CA$eig[2]/sum(PCoA$CA$eig),2)*100,"%"), color = "Sample Type") +  # Adjust labels as needed
  theme(legend.title = element_text(face = "bold"))  # Bold legend titles

# Display the plot
print(pca_plot_customized)

#PERMANOVAS
##Between bacterioplankton and Daphnia

# Perform betadisper to assess homogeneity of variances
dispersion_test <- betadisper(bray, pcoa_df$Type)

# To see the results
print(dispersion_test)

# To visualize the dispersion
plot(dispersion_test)
anova(dispersion_test) #0.01427 *

##Between High vs Low - Daphnia
meta_daphnia = meta %>%
  filter(Type=="Daphnia")

proportions_counts_daphnia <- proportions_counts %>%
  filter(rownames(proportions_counts) %in% rownames(meta_daphnia))

bray_daphnia <- proportions_counts_daphnia %>%
  vegdist(method = "bray")

# Perform betadisper to assess homogeneity of variances
dispersion_test <- betadisper(bray_daphnia, meta_daphnia$Category)

# To see the results
print(dispersion_test)

# To visualize the dispersion
plot(dispersion_test)
anova(dispersion_test) #0.0825 .

#adonis
# Performing PERMANOVA using the adonis function from the vegan package
permanova_result <- adonis2(bray_daphnia~Category,data=meta_daphnia, method = "bray")

# Viewing the results
print(permanova_result) #0.005 **

##Between High vs Low - Bacterioplankton
meta_bac = meta %>%
  filter(Type=="Bacterioplankton")

proportions_counts_bac <- proportions_counts %>%
  filter(rownames(proportions_counts) %in% rownames(meta_bac))

bray_bac <- proportions_counts_bac %>%
  vegdist(method = "bray")

# Perform betadisper to assess homogeneity of variances
dispersion_test <- betadisper(bray_bac, meta_bac$Category)

# To see the results
print(dispersion_test)

# To visualize the dispersion
plot(dispersion_test)
anova(dispersion_test) #0.01433 *

###End of permanovas

#Divide into Bacterioplankton / Daphnia

###DAPHNIA
ps=phyloseq(otu_table(counts,taxa_are_rows = TRUE), sample_data(meta))
ps_dap=subset_samples(ps,Type == "Daphnia")
counts_dap=data.frame(otu_table(ps_dap))
meta_dap=data.frame(sample_data(ps_dap))

dds <- DESeqDataSetFromMatrix(countData = counts_dap,
                              colData = meta_dap,
                              design = ~ Category)
dds <- DESeq(dds)
res <- results(dds)
res

#Remove low count COGs
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
res = results(dds)
res05 = results(dds,alpha=0.05) #p-value <0.05
summary(res05)

"
out of 16960 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 6093, 36%
LFC < 0 (down)     : 2943, 17%
outliers [1]       : 0, 0%
low counts [2]     : 699, 4.1%
(mean count < 1)

"

###BACTEROLANKTON
ps=phyloseq(otu_table(counts,taxa_are_rows = TRUE), sample_data(meta))
ps_bac=subset_samples(ps,Type == "Bacterioplankton")
counts_bac=data.frame(otu_table(ps_bac))
meta_bac=data.frame(sample_data(ps_bac))

dds <- DESeqDataSetFromMatrix(countData = counts_bac,
                              colData = meta_bac,
                              design = ~ Category)
dds <- DESeq(dds)
res <- results(dds)
res

#Remove low count COGs
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
res = results(dds)
res05 = results(dds,alpha=0.05) #p-value <0.05
summary(res05)

"
out of 54642 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 5730, 10%
LFC < 0 (down)     : 7548, 14%
outliers [1]       : 0, 0%
low counts [2]     : 7554, 14%
(mean count < 2)

"