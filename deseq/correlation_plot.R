library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

load("RData/data.RData")

outputDirectoryPlots <- file.path(outputDirectory, "plots")

fpkm <- as.data.frame(fpkm(dds))
#conditions <- colData(dds)$condition
fpkm$gene <- rownames(fpkm)

fpkm.tidy <- gather(fpkm, sample, fpkm, -gene)
fpkm.tidy$condition <- gsub("_[1-3]", "", fpkm.tidy$sample)

fpkmscvplot <- function(fpkm.tidy, outputDirectoryPlots) {
    pop.var <- function(x) var(x) * (length(x)-1) / length(x)
    pop.sd <- function(x) sqrt(pop.var(x))

    # Calculate mean and sd for each gene. Remove genes with a mean count of 0.
    fpkm.means.sd <- fpkm.tidy %>% group_by(gene, condition) %>% 
	summarize(mean=mean(fpkm), sd=sd(fpkm), cv=(pop.sd(fpkm)/mean), cv2=cv^2) %>% 
	filter(mean > 1) %>%
	mutate(log10mean=log10(mean))

    FPKMLowerBound = 1

    pdf(file = file.path(outputDirectoryPlots, "fpkmSCVPlot.pdf"))
    p <- ggplot(fpkm.means.sd, aes(x=log10mean, y=cv2)) + # geom_point(aes(color=condition), size=1) +
	xlab(bquote(~log[10]~ " FPKM")) +
	scale_y_continuous(name=bquote(CV^2)) +
	stat_smooth(aes(color=condition, fill=condition), na.rm=TRUE, method="auto", fullrange=TRUE) +
	theme_bw() + 
	ggplot2::xlim(c(log10(FPKMLowerBound), max(fpkm.means.sd$log10mean)))
    p
    dev.off()
}
