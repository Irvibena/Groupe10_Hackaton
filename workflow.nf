//!/usr/bin/env nextflow

// définition des variables

SAMPLES=["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]

// Téléchargement des données

process DownloadFastaFiles {
	output:
	file "srrFiles/{sample}.1" into fasta_files

	script:
	"""
	wget -O {output} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR628/SRR628582/SRR628582_1.fastq.gz
	"""
	}
