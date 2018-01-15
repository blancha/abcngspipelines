#!/usr/bin/env Rscript

#########
# Plots #
#########

library(amap)
library(DESeq2)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(FactoMineR)

load("RData/data.RData")

# Create plot output directory
outputDirectoryPlots <- file.path(outputDirectory, "plots")
dir.create(outputDirectoryPlots)

# Hierarchical clustering.
pdf(file.path(outputDirectoryPlots, "hierachical_clustering_samples.pdf"))
h<-hcluster(t(counts(dds, normalized=TRUE)), method="spearman")
plot(h)
dev.off()

# Get 35 genes with highest FPKM values
most.expressed.genes <- fpkm.counts.annotated %>%
    mutate(sum=rowSums(.[samples])) %>% # Calculate row sums for sample counts.
    arrange(desc(sum)) %>% # Arrange in decreasing order
    slice(1:35) %>% # Select most expressed genes.
    select_(.dots=c("external_gene_name", samples))

# Heatmap of the 35 highest expressed genes
pdf(file.path(outputDirectoryPlots, "heatmap_35_most_expressed_genes_log2_fpkm.pdf"), height=8)
heatmap.2( as.matrix(log2(most.expressed.genes[,2:length(colnames(most.expressed.genes))])), labRow=most.expressed.genes$external_gene_name,
	      scale="none", trace="none", dendrogram="column",
              col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margin=c(13, 16))
dev.off()

# FactoMineR PCA plot
rld <- rlog(dds)
assay.rld.t <- t(assay(rld))
pca <- PCA(assay.rld.t, graph=FALSE)

samples <- rownames(colData(dds))
conditions <- colData(dds)$conditions
summary.conditions <- summary(conditions)

# Compute colors for PCA plot
colors.brewer <- brewer.pal(n=length(levels(conditions)), name="Set1")
colors <- unlist(lapply(seq_along(summary.conditions), function(index) {
		return(rep(colors.brewer[index], summary.conditions[index]))
	      }))

pdf(file.path(outputDirectoryPlots, "PCA_with_colors.pdf"))
plot.PCA(pca, habillage="ind", col.hab=colors)
dev.off()

##########
# MAplot #
##########
res <- results(dds)
resMLE <- results(dds, addMLE=TRUE)
df.MLE <- data.frame(mean = resMLE$baseMean, lfc = resMLE$lfcMLE, isDE = ifelse(is.na(resMLE$padj), FALSE, resMLE$padj < 0.1)) 

pdf(file.path(outputDirectory, "MAplot_unshrunken_log2_fold_change.pdf"))
#par(mfcol=c(1,2))
plotMA(df.MLE, main="unshrunken log2 fold changes", ylim=c(-2,2))
dev.off()

pdf(file.path(outputDirectory, "MAplot_deseq2.pdf"))
#par(mfcol=c(1,2))
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()
