//!/usr/bin/env nextflow

params.design_file="/mnt/mydatalocal/Groupe10_Hackaton/etiquetage_echantillon.csv"
fichier_reads = Channel.of(params.design_file)

// Analyse des données

process dataAnalysis {

publishDir "results/count_output/"
  
input:
path reads from FileReads
path metadata from fichier_mutant

output:
file "histogramm.png" into result_deseq

script :
"""
#!/usr/bin/env Rscript
library("DESeq2")
counts <- read.table("${reads}", row.names = 1, header=T)
counts <- counts[,-1:-5]
counts <- counts[rowSums(counts)>0,]
counts <- as.matrix(counts)
coldata = read.csv("${metadata}", row.names = 1, header=TRUE)
cond <- factor(coldata$etat_echantillon)
dds <- DESeqDataSetFromMatrix(counts, DataFrame(cond), ~ cond)
dds <- DESeq(dds)
res <- results(dds)
png("histogramm.png")
hist(res$pvalue, col="blue", main="")
dev.off()
"""
}
