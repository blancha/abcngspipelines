#!/usr/bin/env Rscript

############################################################
# Calculation of differential methylation with methylKit.  #
#                                                          #
# Author: Alexis Blanchet-Cohen                            #
# Date: July 25th 2014                                     #
# Version 1.0                                              #
############################################################

# Genomes folder with annotation files.
genomes.folder <- "/sb/project/afb-431/genomes"

suppressPackageStartupMessages(library(argparse))
options(echo=TRUE)
options(verbose=TRUE)

###########################################################
# Command line argument parsing.                          #
# To run in terminal directly, comment out this section,  # 
# and just specify the values of the variables directly.  #
###########################################################

# Create parser object
parser <- ArgumentParser()

# Specify desired options
parser$add_argument("-g", "--genome", help="Genome. E.g. mm9, mm10, hg19", default="mm9")
parser$add_argument("-c", "--context", help="Context. CpG, CHG, or CHH", default="CpG")
parser$add_argument("-t", "--treatment_vector", help="Treatment vector E.g. 1,1,0,0", required=TRUE)
parser$add_argument("-o", "--output_directory", help="Output directory", default="../results/methylkit")
parser$add_argument("-n", "--number_cores", help="Number of cores to use by the function calculateDiffMeth.", default=3)
parser$add_argument("-s", "--sample_ids", nargs="*", help="sample IDs E.g. mut1,mut2,wt1,wt2", required=TRUE)
parser$add_argument("-f", "--files", nargs="*", help="Bismark SAM files. E.g. ../data/mut1.sam,../data/mut2.sam", required=TRUE)
parser$add_argument("-m", "--mincov", help="Minimum coverage", default=10)
parser$add_argument("-q", "--minqual", help="Minimum quality", default=20)

# Get command line options. If help option encountered, print help and exit.
# Otherwise, if options not found on command line then set defaults.
args <- parser$parse_args()

# Process command line arguments
genome <- args$genome
context <- tolower(args$context)
treatment.vector <- as.numeric(unlist(strsplit(args$treatment_vector, split=",")))
output.directory <- file.path(args$output_directory, tolower(context))
number.cores <- args$number_cores
sample.ids <- as.list(unlist(strsplit(args$sample_ids, split=",")))
bismark.files <- as.list(unlist(strsplit(args$files, split=",")))
min.cov <- args$mincov
min.qual <- args$minqual

# To run in terminal
# setwd("~/Downloads/MAS-BIF-P5/scripts/methylkit")
# genome <- "mm9"
# context <- "CpG"
# treatment.vector <- c(1,0)
# output.directory <- "../../results/methylkit"
# number.cores <- 3
# sample.ids <- list("CBP_SD", "CBP_Ctl")
# bismark.files <- list("../../results/bismark/CBP_SD/CBP_SD.sam", "../../results/bismark/CBP_Ctl/CBP_Ctl.sam")

library(methylKit)

# Create output directory
dir.create(output.directory, recursive=TRUE)

######################
# Read Bismark files #
######################
                
print("Reading Bismark SAM files ...")
if (context == "chh") context.methylkit.format="CHH" else if (context == "chg") context.methylkit.format="CHG" else if (context == "cpg") context.methylkit.format="CpG"
# Read files
myMethylRawList.before.normalization = read.bismark(location=bismark.files,
                sample.id=sample.ids,assembly=genome,
                save.folder=file.path(output.directory,"methylation_files"),save.context=c(context.methylkit.format),read.context=context.methylkit.format,
                nolap=FALSE,mincov=min.cov,minqual=min.qual,treatment=treatment.vector)              
print("Finished reading Bismark SAM files.")

print("Normalizing coverage ...")             
# Normalize coverage
myMethylRawList.after.normalization <- normalizeCoverage(myMethylRawList.before.normalization)
print("Coverage normalized.")

# To free up memory, remove myMethylRawList.before.normalization from workspace after normalisation
rm(myMethylRawList.before.normalization)

#######################
# Generate statistics #
#######################

# Function to calculate statistics for each sample in list of samples
getStatistics <- function(myMethylRaw) {
	sample <- myMethylRaw@sample.id
	
	##########################
	# Methylation statistics #
	##########################
	# Create output directory
	methylation.statistics.directory <- file.path(output.directory, "methylation_statistics")
	graphsDirectory <- file.path(methylation.statistics.directory, "graphs")
	statisticsDirectory <- file.path(methylation.statistics.directory, "statistics")
	dir.create(graphsDirectory, recursive=TRUE, showWarnings=FALSE)
	dir.create(statisticsDirectory, showWarnings=FALSE)	
	# direct output to a file 
	sink(file.path(statisticsDirectory,paste(sample, "_methylation_statistics.txt", sep="")), append=FALSE, split=FALSE)
	print(sample)
	getMethylationStats(myMethylRaw,plot=F,both.strands=F)
	# return output to the terminal 
	sink()
	pdf(file.path(graphsDirectory, paste(sample, "_methylation_graph.pdf", sep="")))
	getMethylationStats(myMethylRaw,plot=T,both.strands=F)
	dev.off()

	#######################
	# Coverage statistics #
	#######################
	# Create output directory
	coverage.statistics.directory <- file.path(output.directory, "coverage_statistics")
	graphsDirectory <- file.path(coverage.statistics.directory, "graphs")
	statisticsDirectory <- file.path(coverage.statistics.directory, "statistics")
	dir.create(graphsDirectory, recursive=TRUE, showWarnings=FALSE)
	dir.create(statisticsDirectory, showWarnings=FALSE)	
	# direct output to a file 
	sink(file.path(statisticsDirectory,paste(sample, "_coverage_statistics.txt", sep="")), append=FALSE, split=FALSE)
	print(sample)
	getCoverageStats(myMethylRaw,plot=F,both.strands=F)
	# return output to the terminal 
	sink()
	pdf(file.path(graphsDirectory, paste(sample, "_coverage_graph.pdf", sep="")))
	getCoverageStats(myMethylRaw,plot=T,both.strands=F)
	dev.off()	
}

lapply(myMethylRawList.after.normalization, getStatistics) 

# The "unite" function unites a methylRawList object so that only bases with coverage from all samples are retained. The resulting object is a class of methylBase
myMethylBase.after.normalisation <- unite(myMethylRawList.after.normalization, destrand=FALSE)

# To free up memory, remove myMethylRawList.after.normalization from workspace after uniting all samples together.
rm(myMethylRawList.after.normalization)

# Save myMethylBase.after.normalisation, to be able to restart from this point if necessary.
save(myMethylBase.after.normalisation, file="myMethylBase_after_normalisation.RData")

###############
# Correlation #
###############
correlation.statistics.directory <- file.path(output.directory, "correlation_statistics")
graphsDirectory <- file.path(correlation.statistics.directory, "graphs")
statisticsDirectory <- file.path(correlation.statistics.directory, "statistics")
dir.create(graphsDirectory, recursive=TRUE, showWarnings=FALSE)
dir.create(statisticsDirectory, showWarnings=FALSE)	
# direct output to a file 
sink(file.path(statisticsDirectory, "correlation_statistics.txt"), append=FALSE, split=FALSE)
getCorrelation(myMethylBase.after.normalisation,plot=FALSE)
# return output to the terminal 
sink()
png(file.path(graphsDirectory, "correlation_graph.png"), width=960, height=960)
getCorrelation(myMethylBase.after.normalisation, plot=TRUE)
dev.off()	

######################
# Clustering samples #
######################
clustering.statistics.directory <- file.path(output.directory, "clustering_statistics")
graphsDirectory <- file.path(clustering.statistics.directory, "graphs")
statisticsDirectory <- file.path(clustering.statistics.directory, "statistics")
dir.create(graphsDirectory, recursive=TRUE, showWarnings=FALSE)
dir.create(statisticsDirectory, showWarnings=FALSE)	
# direct output to a file 
sink(file.path(statisticsDirectory, "clustering_statistics.txt"), append=FALSE, split=FALSE)
clustersamples(myMethylBase.after.normalisation, dist="correlation", method="ward", plot=FALSE)
# return output to the terminal 
sink()
png(file.path(graphsDirectory, "clustering_graph.png"), width=960, height=960)
clustersamples(myMethylBase.after.normalisation, dist="correlation", method="ward", plot=TRUE)
dev.off()	

#######
# PCA #
#######
PCADirectory <- file.path(output.directory, "PCA")
graphsDirectory <- file.path(PCADirectory, "graphs")
dir.create(graphsDirectory, recursive=TRUE, showWarnings=FALSE)
pdf(file.path(graphsDirectory, "PCA_screeplot_graph.pdf"))
PCAsamples(myMethylBase.after.normalisation, screeplot=TRUE)
dev.off()	
pdf(file.path(graphsDirectory, "PCA_graph.pdf"))
PCAsamples(myMethylBase.after.normalisation)
dev.off()	

####################################
# Methylation differences of bases #
####################################

differences.bases.directory <- file.path(output.directory, "methylation_differences", "bases")
dir.create(differences.bases.directory, recursive=TRUE)

print("Calculating methylation differences ...")
myMethylDiff.bases <- calculateDiffMeth(myMethylBase.after.normalisation, num.cores=number.cores)
print("Finished calculating methylation differences.")

print("Printing methylation differences to file ...")
write.table(myMethylDiff.bases, file=file.path(differences.bases.directory,  "methylation_differences.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
print("Methylation differences printed to file.")

# Save myMethylDiff.bases here, to be able to restart from this point if necessary.
save(myMethylDiff.bases, file="myMethylDiff_bases.RData")

# To free up memory, remove myMethylDiff.bases
rm(myMethylDiff.bases)

###################
# 100 bases tiles #
###################

differences.tiles.100.directory <- file.path(output.directory, "methylation_differences", "tiles", "100")
dir.create(differences.tiles.100.directory, recursive=TRUE)

# Summarize methylated/unmethylated base counts over tiling windows accross genome.
myMethylBase.tiled.100 <- tileMethylCounts(object=myMethylBase.after.normalisation,win.size=100, step.size=100,cov.bases=0)

print("Calculating methylation differences for 100 bases tiles ...")
myMethylDiff.tiled.100 <- calculateDiffMeth(myMethylBase.tiled.100, num.cores=number.cores)
print("Finished calculating methylation differences for 100 bases tiles ...")

print("Writing 100 bases tiled methylation differences to file ...")
write.table(myMethylDiff.tiled.100, file=file.path(differences.tiles.100.directory,  "methylation_differences_tiles_100_bases.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
print("Finished writing 100 bases tiled methylation differences to file ...")

# Save myMethylDiff.tiled.100 here, to be able to restart from this point if necessary.
save(myMethylDiff.tiled.100, file="myMethylDiff_tiled_100.RData")

# To free up memory, remove myMethylBase.tiled.100 and myMethylDiff.tiled.100
rm(myMethylBase.tiled.100, myMethylDiff.tiled.100)

####################
# 1000 bases tiles #
####################

differences.tiles.1000.directory <- file.path(output.directory, "methylation_differences", "tiles", "1000")
dir.create(differences.tiles.1000.directory, recursive=TRUE)

# Summarize methylated/unmethylated base counts over tiling windows accross genome.
myMethylBase.tiled.1000 <- tileMethylCounts(object=myMethylBase.after.normalisation,win.size=1000, step.size=1000,cov.bases=0)

print("Calculating methylation differences for 1000 bases tiles ...")
myMethylDiff.tiled.1000 <- calculateDiffMeth(myMethylBase.tiled.1000, num.cores=number.cores)
print("Finished calculating methylation differences for 1000 bases tiles ...")

print("Writing 1000 bases tiled methylation differences to file ...")
write.table(myMethylDiff.tiled.1000, file=file.path(differences.tiles.1000.directory,  "methylation_differences_tiles_1000_bases.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
print("Finished writing 1000 bases tiled methylation differences to file ...")

# Save myMethylDiff.tiled.1000 here, to be able to restart from this point if necessary.
save(myMethylDiff.tiled.1000, file="myMethylDiff_tiled_1000.RData")