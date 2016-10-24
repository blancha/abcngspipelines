#!/usr/bin/env Rscript

library("biomaRt")

ensembl_biomart = useMart('ensembl',dataset='mmusculus_gene_ensembl')
gene.annotation=getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description'), mart=ensembl_biomart)

gene.annotation <- gene.annotation[order(gene.annotation$gene_biotype),]

write.table(gene.annotation, "annotation.txt", sep="\t", row.names=FALSE, quote=FALSE)
