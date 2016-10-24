#!/usr/bin/env Rscript

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 05/09/2014

samples <- read.table("../samples.txt", header=TRUE)

f <- function(x, output) {
 sample <- x[3]
 filename <- paste(x[3], ".txt", sep="")
 condition <- x[4]
 write(paste(sample, filename, condition, sep="\t"), file= output, append = TRUE)
}


write("sample\tfilename\tcondition", "deseqDesign.txt")

apply(samples, 1, f, output="deseqDesign.txt")
