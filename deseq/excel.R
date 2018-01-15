#!/usr/bin/env Rscript

library(data.table)
library(openxlsx)

inputDirectory  <- "../../results/deseq/tsv_files"

file.names <- list.files(inputDirectory, recursive=TRUE, full.names=TRUE)

excel <- function(file.name) {
    data <- fread(file.name)
    directory.name <- gsub("tsv_files", "excel_files", dirname(file.name))
    dir.create(directory.name, recursive=TRUE, showWarnings=FALSE)
    file.name <- gsub(".tsv", ".xlsx", basename(file.name), fixed=TRUE)
    write.xlsx(data, file=file.path(directory.name, file.name))
}

log <- lapply(file.names, excel)
