#!/usr/bin/env Rscript

options(echo=TRUE)

library(gtools)
library(xlsx)

#setwd("/sb/project/afb-431/BIF/RIC/RIC-BIF-P1/results/methylkit/poly17_vs_ctl17/cpg/methylation_differences/bases")

file <- Sys.glob("methylation_differences*.tsv")[1]

# Sort chromosomes in numeric mixedorder
methylation.differences <- read.table(file, sep="\t", header=TRUE, comment.char="", quote="")

# Keep only rows with p-values under 1%
methylation.differences.1.percent <- methylation.differences[(methylation.differences$pvalue<0.01),]

# Sort chromosomes in numerical mixedorder
methylation.differences.1.percent <- methylation.differences.1.percent[mixedorder(methylation.differences.1.percent$chr),]

# Generate BED file in numerical mixedorder
methylation.differences.1.percent.homer.bed <- data.frame(methylation.differences.1.percent$chr,
														methylation.differences.1.percent$start,
														methylation.differences.1.percent$end,
														paste(methylation.differences.1.percent$chr, ":", methylation.differences.1.percent$start, "_", methylation.differences.1.percent$end, sep=""),
														"\t",
														methylation.differences.1.percent$strand
														)
# Write BED file in Homer format							
write.table(methylation.differences.1.percent.homer.bed, file="methylation_differences_p_value_1_percent.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Annotate BED file with Homer annotatePeaks.pl
system("annotatePeaks.pl methylation_differences_p_value_1_percent.bed mm9 > homer_results.tsv")

# Read Homer results
homer.results <- read.table("homer_results.tsv", sep="\t", header=TRUE, comment.char="", quote="")

# Sort Homer results in numerical mixedorder
homer.results <- homer.results[mixedorder(homer.results$PeakID..cmd.annotatePeaks.pl.methylation_differences_p_value_1_percent.bed.mm9.),]

#Merge Homer results and methylkit results
methylation.differences.1.percent.annotated <- cbind(methylation.differences.1.percent, homer.results)

# Remove columns of no interest
methylation.differences.1.percent.annotated <-  methylation.differences.1.percent.annotated[,c(1:8, 15:ncol(methylation.differences.1.percent.annotated))]

# Write data to Excel file.
write.xlsx2(methylation.differences.1.percent.annotated, file="methylation_differences_p_value_1_percent_annotated.xlsx", row.names=FALSE)