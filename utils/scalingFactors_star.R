#!/usr/bin/env Rscript

###################################
# Process command-line arguments. #
###################################
suppressPackageStartupMessages(library(argparse))

# Create parser object
parser <- ArgumentParser()

parser$add_argument("-m", "--min_library", help="Specify size of minimum library. If left unspecified, Star align.summary files will be read to identify library with minimum size.", default="Unspecified")
parser$add_argument("-o", "--output_file", help="Specify output file. DEFAULT=scaling_factors.txt", default="scaling_factors.txt")

args <- parser$parse_args()

min_library <- args$min_library 
output_file <- args$output_file

inputDirectory <- "../results/star"

# List all samples in Star directory
samples  = list.files(inputDirectory)

results <- lapply(samples, function(sample) {
    align.summary.file <- file.path(inputDirectory, sample, paste0(sample, "Log.final.out"))
    if (!file.exists(align.summary.file)) {
	print("Error!")
	print("The specified file does not exist.")
	print(align.summary.file)
	quit(status=1)
    }
    lines <- readLines(align.summary.file)
    number.of.input.reads <- lines[grep("Number of input reads", lines)]
    number.of.input.reads <- as.numeric(strsplit(number.of.input.reads, "\t")[[1]][2])

    mismatches <- lines[grep("too many mismatches", lines)]
    mismatches <- strsplit(mismatches, "\t")[[1]][2]
    mismatches <- as.numeric(gsub("%", "", mismatches))
    short <- lines[grep("too short", lines)]
    short <- strsplit(short, "\t")[[1]][2]
    short <- as.numeric(gsub("%", "", short))
    other <- lines[grep("unmapped: other", lines)]
    other  <- strsplit(other, "\t")[[1]][2]
    other <- as.numeric(gsub("%", "", other))

    unmapped.percentage <- mismatches + short + other
    unmapped <-  (number.of.input.reads * unmapped.percentage) / 100
    mapped <- round(number.of.input.reads - unmapped)

    return(c(name=sample, number.of.reads=mapped))
})

scaling.factors <- as.data.frame(do.call(rbind, results), stringsAsFactors=FALSE)
scaling.factors$number.of.reads <- as.numeric(scaling.factors$number.of.reads)

if (min_library == "Unspecified") minimum = min(scaling.factors$number.of.reads) else 
    minimum = as.numeric(min_library)

scaling.factors$scaling_factor <- round(minimum/scaling.factors$number.of.reads, digits=2)

write.table(scaling.factors, output_file, sep="\t", quote=FALSE, row.names=FALSE) 
cat(paste0("Wrote file ", output_file, ".\n"))
