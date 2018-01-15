#!/usr/bin/env Rscript
print("Loading packages ...")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
print("Packages loaded.")

# Comparisons
configuration = readLines("../configuration_project.txt")
comparisons = configuration[grep("comparisons", configuration)]
comparisons = str_trim(strsplit(strsplit(comparisons, "comparisons=")[[1]][2], ",")[[1]])
comparisons = lapply(comparisons, function(comparison) {return(strsplit(comparison, " vs ")[[1]])})

inputDirectory <- "../../results/featurecounts_norrna_s2r"
outputDirectory <- "../../results/deseq"
dir.create(outputDirectory)

print("Reading samples file ...")
# Obtain the list of samples, and their condition
samplesFile <- fread("../samples.txt")
# Only keep samples being compared.
samplesFile <- samplesFile[samplesFile$Condition %in% unlist(comparisons),]
conditions <- factor(samplesFile$Condition)
samples <- samplesFile$Sample
print("Samples file read.")

print("Reading featureCounts files ...")
# Read the featureCount files.
countData.list <- lapply(samples, function(sample) {
    file <- Sys.glob(file.path(inputDirectory, paste0(sample, "*txt")))
    print(file)
    dd <- fread(file, skip=1)
    dd <- dd %>% select(Geneid, Length, contains("accepted_hits.bam"))
    setnames(dd, 3, sample)
    return(dd)
})
print("FeatureCounts file read.")

# Convert the list of dataframes to one dataframe
countData <- as.data.frame(do.call(cbind, countData.list))

# Remove the duplicate columns
countData <- countData[, unique(colnames(countData))]

colData <- data.frame(conditions, row.names=colnames(countData)[3:ncol(countData)])

# Create a DESeqDataSet object.
dds <- DESeqDataSetFromMatrix(countData=subset(countData, select=-Length), colData=colData, design= ~ conditions, tidy=TRUE)
dds <- DESeq(dds)

# Give the gene lengths to DESeq
mcols(dds)$basepairs <- countData$Length

# Read the annotation
annotation <- fread("annotation.txt", data.table=FALSE)

# Create output directories
dir.create(file.path(outputDirectory, "tsv_files", "counts_normalized_relative_to_library_size"), recursive=TRUE)
dir.create(file.path(outputDirectory, "tsv_files", "fpkm_counts"))

# Calculate and save normalized counts
normalized.counts <- counts(dds, normalized=TRUE)
normalized.counts.annotated <- merge(annotation, normalized.counts, by.x="ensembl_gene_id", by.y="row.names", all=TRUE)
if (length(levels(conditions)) >= 2) {
    write.table(normalized.counts.annotated, file.path(outputDirectory, "tsv_files", "counts_normalized_relative_to_library_size", "counts_normalized_relative_to_library_size.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
}

# Calculate and save FPKM counts
fpkm.counts <- fpkm(dds)
fpkm.counts.annotated <- merge(annotation, fpkm.counts, by.x="ensembl_gene_id", by.y="row.names", all=TRUE)
if(length(levels(conditions)) >= 2) {
    write.table(fpkm.counts.annotated, file.path(outputDirectory, "tsv_files", "fpkm_counts", "fpkm_counts.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
}

#####################################
# Differential expression analysis. # 
#####################################

differential.expression <- function(comparison) {
    condition.numerator <- comparison[1]
    condition.denominator <- comparison[2]
    if (substr(condition.numerator, 1, 2) == "TS") {
	condition.numerator.match <- paste0(strsplit(condition.numerator, "_KO_")[[1]][1], "_KO_[12]_", strsplit(condition.numerator, "_")[[1]][3])
    }
    else {
	condition.numerator.match <- paste0(strsplit(condition.numerator, "_")[[1]][1], "_[12]_", strsplit(condition.numerator, "_")[[1]][2])
    }
    if (substr(condition.denominator, 1, 2) == "TS") {
        condition.denominator.match <- paste0(strsplit(condition.denominator, "_KO_")[[1]][1], "_KO_[12]_", strsplit(condition.denominator, "_")[[1]][3])
    }
    else {
        condition.denominator.match <- paste0(strsplit(condition.denominator, "_")[[1]][1], "_[12]_", strsplit(condition.denominator, "_")[[1]][2])
    }
    # Check order if levels, putting the denonimator first, if necessary.
    levels <- levels(colData(dds)$condition)
    if(which(levels==condition.denominator) > which(levels==condition.numerator)) {
	colData(dds)$condition <- relevel(colData(dds)$condition, condition.numerator)
	dds <- DESeq(dds)
    }
    deseq.results <- results(dds, contrast=list(paste0("conditions", condition.numerator), paste0("conditions", condition.denominator)))
    # Add annotation and gene counts (only for samples being compared)
    normalized.counts.annotated.filtered <- normalized.counts.annotated %>% select(ensembl_gene_id, external_gene_name, gene_biotype, description, matches(condition.numerator.match), matches(condition.denominator.match))
    normalized.counts.with.deseq.results.annotated <- merge(normalized.counts.annotated.filtered, as.data.frame(deseq.results), by.x="ensembl_gene_id", by.y="row.names")
    fpkm.counts.annotated.filtered <- fpkm.counts.annotated %>% select(ensembl_gene_id, external_gene_name, gene_biotype, description, matches(condition.numerator.match), matches(condition.denominator.match))
    fpkm.counts.with.deseq.results.annotated <- merge(fpkm.counts.annotated.filtered, as.data.frame(deseq.results), by.x="ensembl_gene_id", by.y="row.names")
    # Order by adjusted p-value and p-value
    normalized.counts.with.deseq.results.annotated <- normalized.counts.with.deseq.results.annotated %>% arrange(padj, pvalue)
    fpkm.counts.with.deseq.results.annotated <- fpkm.counts.with.deseq.results.annotated %>% arrange(padj, pvalue)
    # Write to file
    write.table(normalized.counts.with.deseq.results.annotated, 
		file.path(outputDirectory, "tsv_files", "counts_normalized_relative_to_library_size", paste0(condition.numerator, "_vs_", condition.denominator, "_counts_normalized_relative_to_library_size.tsv")),
		quote=FALSE, row.names=FALSE, sep="\t")
    write.table(fpkm.counts.with.deseq.results.annotated, 
		file.path(outputDirectory, "tsv_files", "fpkm_counts", paste0(condition.numerator, "_vs_", condition.denominator, "_fpkm_counts.tsv")),
                quote=FALSE, row.names=FALSE, sep="\t")
}

log <- lapply(comparisons, differential.expression)

dir.create("RData")
save.image("RData/data.RData")
