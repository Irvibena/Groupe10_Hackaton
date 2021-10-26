//!/usr/bin/env nextflow

// définition des variables

samples = Channel.from("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589")

// Téléchargement des données

process DownloadFastqFiles {
    input:
    val sample from samples

	output:
	file "${sample}_1.fastq.gz" into fastq_files

	script:
	"""
	wget -O ${sample}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample.substring(0,6)}/${sample}/${sample}_1.fastq.gz
	"""
	}
