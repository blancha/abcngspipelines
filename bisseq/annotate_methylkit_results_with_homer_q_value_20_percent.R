#!/usr/bin/env Rscript

options(echo=TRUE)

library(gtools)
library(xlsx)

file <- Sys.glob("methylation_differences*.tsv")[1]

# Sort chromosomes in numeric mixedorder
methylation.differences <- read.table(file, sep="\t", header=TRUE, comment.char="", quote="")

# Keep only rows with p-values under 20%
methylation.differences.q.20.percent <- methylation.differences[(methylation.differences$qvalue<0.2),]

# Sort chromosomes in numerical mixedorder
methylation.differences.q.20.percent <- methylation.differences.q.20.percent[mixedorder(methylation.differences.q.20.percent$chr),]

# Generate BED file in numerical mixedorder
if (nrow(methylation.differences.q.20.percent) != 0) {
    methylation.differences.q.20.percent.homer.bed <- data.frame(methylation.differences.q.20.percent$chr,
                                                            methylation.differences.q.20.percent$start,
                                                            methylation.differences.q.20.percent$end,
                                                            paste(methylation.differences.q.20.percent$chr, ":", methylation.differences.q.20.percent$start, "_", methylation.differences.q.20.percent$end, sep=""),
                                                            "\t",
                                                            methylation.differences.q.20.percent$strand
                                                            )} else { # Empty dataframe
    methylation.differences.q.20.percent.homer.bed <- NA
}

# Write BED file in Homer format							
write.table(methylation.differences.q.20.percent.homer.bed, file="methylation_differences_q_value_20_percent.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Annotate BED file with Homer annotatePeaks.pl
system("annotatePeaks.pl methylation_differences_q_value_20_percent.bed mm9 > homer_results_q_value_20_percent.tsv")

# Read Homer results
homer.results <- read.table("homer_results_q_value_20_percent.tsv", sep="\t", header=TRUE, comment.char="", quote="")

# Sort Homer results in numerical mixedorder
homer.results <- homer.results[mixedorder(homer.results$PeakID..cmd.annotatePeaks.pl.methylation_differences_q_value_20_percent.bed.mm9.),]

#Merge Homer results and methylkit results
methylation.differences.q.20.percent.annotated <- cbind(methylation.differences.q.20.percent, homer.results)

# Remove columns of no interest
methylation.differences.q.20.percent.annotated <-  methylation.differences.q.20.percent.annotated[,c(1:8, 15:ncol(methylation.differences.q.20.percent.annotated))]

# Write data to TSV file
write.table(methylation.differences.q.20.percent.annotated, file="methylation_differences_q_value_20_percent_annotated.tsv", row.names=FALSE, quote=FALSE, sep="\t")

# Determine if the folder is CHG, CHH or CPG
# For CPG, generate Excel files.
# For CHH, generated TSV files.
#if (!grepl(pattern="cpg", x=getwd())) {
#  print("cpg directory")
  # Write data to Excel file.
#  write.xlsx2(methylation.differences.q.20.percent.annotated, file="methylation_differences_q_value_20_percent_annotated.xlsx", row.names=FALSE)
#} else {
#  print("Not cpg directory")
#  write.table(methylation.differences.q.20.percent.annotated, file="methylation_differences_q_value_20_percent_annotated.tsv", row.names=FALSE, quote=FALSE, sep="\t")
#}
