#!/usr/bin/env Rscript

#####################################################################
# Formats the CSV file downloaded from Génome Québec's Nanuq server #
# into a format that can be included in the LaTex report.           #
#####################################################################

inputDirectory <- "../../data/FASTQ_files/untrimmed"

nanuq.file <- dir(inputDirectory, pattern="*.csv", full.names=TRUE)

# Check if Nanuq file exists
if (length(nanuq.file) == 0) {
    print(paste("There is no csv file in", inputDirectory))
    quit(status = 1)    
}


# Read file
print(paste("Reading", nanuq.file))
nanuq <- read.table(nanuq.file, sep=",", header=TRUE, check.names=FALSE, colClasses="character")

# Convert all column names to lower case.
colnames(nanuq) <- tolower(colnames(nanuq))

# Select desired columns, and order rows by alphabetical order
nanuqFormatted <- nanuq[order(nanuq$name),
                  c("name","run","run type","library type","number of reads","number of cycles","% duplicate")]

# Give names that will appear in report.
colnames(nanuqFormatted) <- c("Name","Run","Run type","Library type","Reads","Cycles","Duplicates")
nanuqFormatted$Duplicates <- format(as.numeric(nanuqFormatted$Duplicates), digits=2)
nanuqFormatted$Duplicates <- paste(nanuqFormatted$Duplicates, "%")

# Write file that will be used to create table in report.
print("Generating nanuq.csv file, to be included in the report as a table.")
write.table(nanuqFormatted, file="nanuq.csv", row.names=FALSE, quote=TRUE, sep=",")

