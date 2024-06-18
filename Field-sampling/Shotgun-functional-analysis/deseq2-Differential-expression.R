"Deseq2: differential expression"

counts = read.csv('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/COGs/COGs_all_samples_raw_counts.tsv',row.names = 1,sep="\t")
meta = read.csv('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/Eggnog/all_samples/meta.csv',row.names = 1, header=TRUE)
#meta = data.frame(meta_all$V2)
#rownames(meta)=rownames(meta_all)

counts <- counts[, rownames(meta)]
all(rownames(meta) == colnames(counts))
colnames(meta)[2] = "condition"

#Remove low-counts COGs
threshold <- 10
n_samples <- 3

# Apply the threshold
keep <- rowSums(counts >= threshold) >= n_samples

# Filter the countData matrix
filtered_countData <- counts[keep, ]

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = filtered_countData,
                              colData = meta,
                              design = ~ condition)
#run DESeq:
dds <- DESeq(dds)

resultsNames(dds)
# Shrinking log fold changes using apeglm:
res_shrink <- lfcShrink(dds, coef="condition_Natural.lake_vs_Artifical.pond", type="apeglm")

# Subset significant genes
sig_genes <- subset(res_shrink, padj < 0.01)
sig_genes = data.frame(sig_genes)

sig_gene_sorted <- rownames(sig_genes[order(sig_genes$log2FoldChange,decreasing = TRUE), ]>0)

# Filter COGs that begin with "COG"
cogs_start_with_cog <- grep("^COG", sig_gene_sorted, value = TRUE)

# Save the filtered COGs to a CSV file
#write.csv(cogs_start_with_cog, "/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/Eggnog/Diff-exps-COGs_filtered.csv", row.names = FALSE)
