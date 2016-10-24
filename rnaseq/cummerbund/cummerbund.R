#!/usr/bin/env Rscript
print("Loading packages ...")
suppressPackageStartupMessages(library(cummeRbund))
print("Loaded packages.")

inputDirectory <- "../../results/cuffdiff"
outputDirectoryPlots <- "../../results/cummerbund/plots"

# Create ouput directory for plots
dir.create(outputDirectoryPlots, recursive=TRUE, showWarnings=FALSE) 

# Read the Cufflinks data
cuff <- readCufflinks(inputDirectory, rebuild=TRUE)

# Information about the run
cuff
runInfo(cuff)
replicates(cuff)
gene.fpkm<-fpkm(genes(cuff))
head(gene.fpkm)
gene.repFpkm<-repFpkm(genes(cuff))
head(gene.repFpkm)

# Dispersion plot
disp <- dispersionPlot(genes(cuff))
ggsave(filename=file.path(outputDirectoryPlots, "dispersionPlot.pdf"), plot=disp)
ggsave(filename=file.path(outputDirectoryPlots, "dispersionPlot.png"), plot=disp)

# Square coefficient of variation
genes.scv <- fpkmSCVPlot(genes(cuff))
ggsave(filename=file.path(outputDirectoryPlots, "fpkmSCVPlot_genes.pdf"), plot=genes.scv)

isoforms.scv <- fpkmSCVPlot(isoforms(cuff))
ggsave(filename=file.path(outputDirectoryPlots, "fpkmSCVPlot_isoforms.pdf"), plot=isoforms.scv)

# csDensity
dens <- csDensity(genes(cuff))
ggsave(filename=file.path(outputDirectoryPlots, "density.pdf"), plot=dens)

densRep <- csDensity(genes(cuff),replicates=TRUE)
ggsave(filename=file.path(outputDirectoryPlots, "density_with_replicates.pdf"), plot=densRep)

# Box plots
b <- csBoxplot(genes(cuff))
ggsave(filename=file.path(outputDirectoryPlots, "boxPlot.pdf"), plot=b)

bRep <- csBoxplot(genes(cuff), replicates=TRUE)
ggsave(filename=file.path(outputDirectoryPlots, "boxPlot_with_replicates.pdf"), plot=bRep)

# Scatter matrix
s <- csScatterMatrix(genes(cuff))
ggsave(filename=file.path(outputDirectoryPlots, "scatterMatrix.pdf"), plot=s)
ggsave(filename=file.path(outputDirectoryPlots, "scatterMatrix.png"), plot=s)

# Individual scatter plot
#s<-csScatter(genes(cuff),"Ctrl","sleep_dep",smooth=T)
#ggsave("../../cummerbund/scatter.pdf")

# Dendrogram
pdf(file.path(outputDirectoryPlots, "dendrogram.pdf"))
d <- csDendro(genes(cuff), replicates=TRUE)
dev.off()

# Volcano matrix
v <- csVolcanoMatrix(genes(cuff))
ggsave(filename=file.path(outputDirectoryPlots, "volcanoMatrix.pdf"), plot=v)
ggsave(filename=file.path(outputDirectoryPlots, "volcanoMatrix.png"), plot=v)

# Pairwise comparisons. The names must be specified.
#v<-csVolcano(genes(cuff), "Ctrl", "sleep_dep")
#ggsave("../../cummerbund/volcanoPlot.pdf")

# Individual gene analysis
#myGeneId <- "ENSMUSG00000000142"
#myGene<-getGene(cuff,myGeneId)
#gl<-expressionPlot(myGene)
#ggsave("../../cummerbund/individual_gene.pdf")

#gl.rep<-expressionPlot(myGene, replicates=TRUE)
#ggsave("../../cummerbund/individual_gene_with_replicates.pdf")

# Individual gene analysis with isoforms
#gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=T)
#ggsave("../../cummerbund/individual_gene_with_isoforms.pdf")
