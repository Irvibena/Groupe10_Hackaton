//!/usr/bin/env nextflow

// définition des variables

samples = Channel.from("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589")

liste_chromosomes = Channel.from("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16" , "17", "18", "19", "20", "21", "22", "MT", "X", "Y")


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


// Téléchargement du génome humain


process DownloadChromosomes {

	input:
	val chromo from liste_chromosomes

	output:
	file "*.fa.gz" into genome_humain 


	script:
	"""

	wget -O chromo${chromo}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chromo}.fa.gz
	"""
	}
	
process MergeChromosomes {

	input:
	file all_chromosome from genome_humain.collect()
	
	output:
	file "ref.fa" into genome_merge
	
	script:
	"""
	gunzip -c ${all_chromosome} > ref.fa
	"""
	}


