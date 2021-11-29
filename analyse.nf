//!/usr/bin/env nextflow

params.design_file="/mnt/mydatalocal/Groupe10_Hackaton/etiquetage_echantillon.csv"
params.reads_file="/mnt/mydatalocal/Groupe10_Hackaton/count_file"
fichier_meta = Channel.of(params.design_file)
fichier_reads = Channel.of(params.reads_file)

// Analyse des données

process dataAnalysis {

	publishDir "results/plots/"
  
	input:
	path reads from fichier_reads
	path metadata from fichier_meta

	output:
	file "histogramm.png" into histo
	file "histogramm.png" into MA
	file "plot_counts.png" into plots_counts
	file "dispersion.png" into disp
	file "deseq_res_1.csv" into object_res_1
	file "deseq_res_2.csv" into object_res_2

	script :
	"""
	#!/usr/bin/env Rscript
	library("DESeq2")
	counts <- read.table("${reads}", sep="\t", row.names = 1, header=T)
	coldata = read.csv("${metadata}", row.names = 1, header=TRUE)
	counts <- counts[,-1:-5]
	counts <- counts[,row.names(coldata)]
	counts <- counts[rowSums(counts)>0,]
	counts <- as.matrix(counts)
	cond <- factor(coldata[,"etat_echantillon"])
	dds <- DESeqDataSetFromMatrix(counts, DataFrame(cond), ~ cond)
	dds <- DESeq(dds)
	res <- results(dds)
	png("histogramm.png")
	hist(res[,"pvalue"], col="blue", main="")
	dev.off()
	png("MA_plot.png")
	plotMA(res,ylim=c(-5,5))
	dev.off()
	png("plot_counts.png")
	plotCounts(dds, gene=which.min(res[,"padj"]), intgroup="cond")
	dev.off()
	png("dispersion.png")
	plotDispEsts(dds)
	dev.off()
	res <- as.data.frame(res)
	write.csv(res, file="deseq_res_1.csv")
	res2 <- results(dds, tidy=TRUE)
	res2 <- as.data.frame(res2)
	write.csv(res2, file="deseq_res_2.csv")
	"""
	}

// plot du Volcano_plot

process plotVolcano {

	publishDir "results/plots/"
  
	input:
	path res from object_res_2

	output:
	file "volcano.png"

	script :
	"""
	#!/usr/bin/env Rscript
	library(ggplot2)
	res = read.csv("${res}", header=TRUE, sep=",")
	res2.condition1 <- subset(res2, padj<.1)
	res2.condition2 <- subset(res2, (padj<.005 & -log10(pvalue) > 6))
	plot.volcano <- ggplot(data = res2.condition2, aes(x = log2FoldChange, y = -log10(pvalue), label = row)) +
	  geom_point(data = res2, aes(x = log2FoldChange, y = -log10(pvalue)), color = "black") +
	  geom_point(data = res2.condition1, aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue") +
	  geom_point(data = res2.condition2, aes(x = log2FoldChange, y = -log10(pvalue)), color = "red") +
	  geom_text(label=res2.condition2$row, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T, size = 2) +
	  xlim(-8, 8)
	ggsave("volcano.png", plot.volcano)
	"""
	}
