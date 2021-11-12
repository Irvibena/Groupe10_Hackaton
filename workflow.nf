//!/usr/bin/env nextflow

// définition des variables

samples = Channel.from("SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589")

// Téléchargement des données

process DownloadFastqFiles {
        input:
        val sample from samples

        output:
        tuple var (sample), file ("${sample}_1.fastq.gz"), file ("${sample}_2.fastq.gz") into fastq_files

        script:
        """
        wget -O ${sample}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample.substring(0,6)}/${sample}/${sample}_1.fastq.gz
        wget -O ${sample}_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample.substring(0,6)}/${sample}/${sample}_2.fastq.gz
        """
        }


// Téléchargement du génome humain


process DownloadChromosomes {

	input:
	val chromo from liste_chromosomess

	output:
	file "*.fa.gz" into genome_humain_zip


	script:
	"""
	wget -O chromo${chromo}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chromo}.fa.gz
	"""
	}

// Dezipe les chromosomes


process DezipChromosomes {

	input:
	file "gen_zip" from genome_humain_zip

	output:
	file "*.fa" into genome_humain


	script:
	"""
	gunzip chromo${chromo}.fa.gz
	"""
	}

// Concaténation des séquences des chromosomes en un seul fichier
process MergeChr {
    input:
    file "gen" from genome_humain_zip.collect()

    output:
    file "ref.fa" into genome_merge

    script :
    """
    gunzip -c *.fa.gz > ref.fa
    """
}


// Indexation du génome humain
process indexGen {
    input :
    file genome from genome_merge

    output :
    file "ref" into star_index

    script :
    """
    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir ref \
    --genomeFastaFiles ${genome}
    """
}
// Télécharger les annotations du génome

process getGenomic_features{
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

	publishDir "bam_files/"

	input:
	tuple var (sample), file fastq1, file fastq2 from fastq_files
	path chemin from star_index.first()

	output:
	file "${sample}.bam" into mapped_fastq

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
    publishDir "bam_files/"

    input:
    file bam from mapped_fastq

    output:
    file "*.bai" into sam_fastq_files

    script:
    """
    samtools index *.bam
    """
}


// Comptage des reads
process getCountReads_feature {

    publishDir "count_output/"
    input:
    file gtf from annotation
    file bam from mapped_fastq.collect()

    output:
    file "output.counts" into FileReads
    file "output.counts.summary" into logsFileReads

    script:
    """
     featureCounts -p -t gene -g gene_id -s 0 -a ${gtf} -o output.counts ${bam}
    """
}
