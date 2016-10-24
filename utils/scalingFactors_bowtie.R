#!/usr/bin/env Rscript

###################################
# Process command-line arguments. #
###################################
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

# Create parser object
parser <- ArgumentParser()

parser$add_argument("-m", "--min_library", help="Specify size of minimum library. If left unspecified, TopHat align.summary files will be read to identify library with minimum size.", default="Unspecified")
parser$add_argument("-o", "--output_file", help="Specify output file. DEFAULT=scaling_factors.txt", default="scaling_factors.txt")

args <- parser$parse_args()

min_library <- args$min_library 
output_file <- args$output_file

inputDirectory <- "bowtie"

# Get samples
samplesFile = fread("samples.txt")
samples = samplesFile$sample
sample = samples[1]
results <- lapply(samples, function(sample) {
    logFile <- file.path(inputDirectory, paste0("bowtie_", sample, ".sh.log"))
    if (!file.exists(logFile)) {
	print("Error!")
	print("The specified file does not exist.")
	print(logFile)
	quit(status=1)
    }
    aligned.1.time <- system(paste("grep 'concordantly exactly 1 time'", logFile), intern=TRUE)
    aligned.1.time <- as.integer(strsplit(trimws(aligned.1.time), split=" " )[[1]][1])
    aligned.many.times <- system(paste("grep 'concordantly >1 times'", logFile), intern=TRUE)
    aligned.many.times <- as.integer(strsplit(trimws(aligned.many.times), split=" " )[[1]][1])
    mapped <- aligned.1.time + aligned.many.times
    return(c(name=sample, number.of.reads=mapped))
})

scaling.factors <- as.data.frame(do.call(rbind, results), stringsAsFactors=FALSE)
scaling.factors$number.of.reads <- as.numeric(scaling.factors$number.of.reads)

if (min_library == "Unspecified") minimum = min(scaling.factors$number.of.reads) else 
    minimum = as.numeric(min_library)

scaling.factors$scaling_factor <- round(minimum/scaling.factors$number.of.reads, digits=2)

write.table(scaling.factors, output_file, sep="\t", quote=FALSE, row.names=FALSE) 
cat(paste0("Wrote file ", output_file, ".\n"))
