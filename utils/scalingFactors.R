#!/usr/bin/env Rscript

###################################
# Process command-line arguments. #
###################################
suppressPackageStartupMessages(library(argparse))

# Create parser object
parser <- ArgumentParser()

parser$add_argument("-m", "--min_library", help="Specify size of minimum library. If left unspecified, TopHat align.summary files will be read to identify library with minimum size.", default="Unspecified")
parser$add_argument("-o", "--output_file", help="Specify output file. DEFAULT=scaling_factors.txt", default="scaling_factors.txt")

args <- parser$parse_args()

min_library <- args$min_library 
output_file <- args$output_file

inputDirectory <- "../results/tophat"

# List all samples in TopHat directory
samples  = list.files(inputDirectory)

results <- lapply(samples, function(sample) {
    align.summary.file <- file.path(inputDirectory, sample, "align_summary.txt")
    if (!file.exists(align.summary.file)) {
	print("Error!")
	print("The specified file does not exist.")
	print(align.summary.file)
	quit(status=1)
    }
    mapped.line <- readLines(file.path(inputDirectory, sample, "align_summary.txt"))[11]
    mapped <- strsplit(mapped.line, "\\s+")[[1]][3]
    return(c(name=sample, number.of.reads=mapped))
})

scaling.factors <- as.data.frame(do.call(rbind, results), stringsAsFactors=FALSE)
scaling.factors$number.of.reads <- as.numeric(scaling.factors$number.of.reads)

if (min_library == "Unspecified") minimum = min(scaling.factors$number.of.reads) else 
    minimum = as.numeric(min_library)

scaling.factors$scaling_factor <- round(minimum/scaling.factors$number.of.reads, digits=2)

write.table(scaling.factors, output_file, sep="\t", quote=FALSE, row.names=FALSE) 
cat(paste0("Wrote file ", output_file, ".\n"))
