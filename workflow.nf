//!/usr/bin/env nextflow

// définition des variables
samples = Channel.from("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589")

// Liste des chromosomes à télécharger
liste_chromosomes = Channel.from("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16" , "17", "18", "19", "20", "21", "22", "MT", "X", "Y")

// Téléchargement des données

process DownloadFastqFiles {

	publishDir "results/samples_sra/"

        input:
        val sample from samples

        output:
        tuple val (sample), file ("${sample}.sra") into sra_files

        script:
        """
        wget -O ${sample}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${sample}/${sample}.1
        """
        }
	
process sraToFastqZip {

	publishDir "results/samples_fastq/"

        input:
        tuple val (sample), file (samp_sra) from sra_files

        output:
        tuple val (sample), file ("${sample}_1.fastq.gz"), file ("${sample}_2.fastq.gz") into fastq_files

        script:
        """
	fastq-dump --gzip --split-files ${samp_sra}
        """
        }


// Téléchargement du génome humain


process DownloadChromosomes {

	publishDir "results/chr/"

	input:
	val chromo from liste_chromosomes

	output:
	file "${chromo}.fa.gz" into genome_humain_zip
	"""
	wget -O ${chromo}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chromo}.fa.gz
	"""
	}

// Concaténation des séquences des chromosomes en un seul fichier
process MergeChr {

	publishDir "results/genome/"

	input:
	file gen from genome_humain_zip.collect()

	output:
	file "ref.fa" into genome_merge

	script :
	"""
	gunzip -c ${gen} > ref.fa
	"""
}


// Indexation du génome humain
process indexGen {

	publishDir "results/gen_index/"

	input :
	file genome from genome_merge.collect()

	output :
	path "ref" into star_index

	script :
	"""
	STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ${genome}
	"""
}
// Télécharger les annotations du génome

process getGenomic_features{

	publishDir "results/gen_features/"

	output:
	file "*.gtf" into annotation

	script:
	"""
	wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
	gzip -d  Homo_sapiens.GRCh38.101.chr.gtf.gz
	"""
}


// Mapping des fichiers FASTQ compressés avec le génome index (alignement des séquences avec le génome humain)

process mapping_Fastq {

	publishDir "results/bam_files/"

	input:
	tuple val(sample), file(fastq1), file(fastq2) from fastq_files
	path chemin from star_index

	output:
	file "${sample}.bam" into mapped_fastq_1, mapped_fastq_2

	script:
	"""
	STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir ${chemin} \
        --readFilesIn <(gunzip -c ${fastq1}) <(gunzip -c ${fastq2}) \
        --runThreadN ${task.cpus} \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM ${task.memory.toBytes()} \
	 >${sample}.bam
	"""
}

// Indexation des fichiers bam

process indexBam {

	publishDir "results/bam_files/"

	input:
	file bam from mapped_fastq_1

	output:
	file "*.bai" into sam_fastq_files

	script:
	"""
	samtools index *.bam
	"""
}


// Comptage des reads
process getCountReads_feature {

	publishDir "results/count_output/"

	input:
	file gtf from annotation
	file bam from mapped_fastq_2.collect()

	output:
	file "output.counts" into FileReads
	file "output.counts.summary" into logsFileReads

	script:
	"""
	featureCounts -p -t gene -g gene_id -s 0 -a ${gtf} -o output.counts ${bam}
	"""
}

params.design_file="/mnt/mydatalocal/Groupe10_Hackaton/etiquetage_echantillon.csv"
fichier_meta = Channel.of(params.design_file)

// Analyse des données

process dataAnalysis {

	publishDir "results/plots/"

	input:
	path reads from FileReads
	path metadata from fichier_meta

	output:
	file "histogramm.png" into histo
	file "MA_plot.png" into MA
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
	res2 = read.csv("${res}", header=TRUE, sep=",")
	res2.condition1 <- subset(res2, padj<.1)
	res2.condition2 <- subset(res2, (padj<.005 & -log10(pvalue) > 6))
	plot.volcano <- ggplot(data = res2.condition2, aes(x = log2FoldChange, y = -log10(pvalue), label = row)) +
	  geom_point(data = res2, aes(x = log2FoldChange, y = -log10(pvalue)), color = "black") +
	  geom_point(data = res2.condition1, aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue") +
	  geom_point(data = res2.condition2, aes(x = log2FoldChange, y = -log10(pvalue)), color = "red") +
	  geom_text(label=res2.condition2[,"row"], nudge_x = 0.25, nudge_y = 0.25, check_overlap = T, size = 2) +
	  xlim(-8, 8)
	ggsave("volcano.png", plot.volcano)
	"""
	}
