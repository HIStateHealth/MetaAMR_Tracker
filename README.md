# MetaAMR_Tracker
MetaAMR-Tracker is a comprehensive metagenomics pipeline designed for tracking Antimicrobial Resistance (AMR) genes and fungi genes in environmental wastewater samples.
From quality control to contig assembly and resistance gene identification, this pipeline integrates state-of-the-art tools like FastQC, MEGAHIT, MetaPhlAn, and ResFinder, 
to deliver a holistic view of the microbial resistome.

![image](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/f834b920-49aa-4d3e-9dec-4e2fcb19e25c)



![image](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/df06b4d7-76af-40c0-814e-65d983f45890)










![pipeline](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/ca9602d0-0d8f-44e6-ad3d-5764e41adeb2)




![pipeline_drawing (1)](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/acb5b2af-0660-414d-a10e-365b4a611294)



# Metagenomic Analysis Pipeline with Nextflow

## Introduction
This pipeline serves as a comprehensive solution for metagenomic data analysis. It is designed to run in a Nextflow environment, leveraging Docker containers for seamless software management. The pipeline incorporates a variety of bioinformatics tools for tasks ranging from quality control to specialized analyses like antibiotic resistance and viral identification.

## Prerequisites
Before you begin, ensure that you have installed the following:
- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/products/docker-desktop)
- Python 3.0 or higher
- git
- 
- All bioinformatics tools listed in the Nextflow script

## Setup
1. Clone the GitHub repository.
    ```bash
    git clone https://github.com/HIStateHealth/MetaAMR_Tracker.git
    ```

2. Navigate into the project directory.
    ```bash
    cd [I will put my Project Directory here]
    ```

## Example Command
Here's an example command to help you get started. Replace the paths with those corresponding to your own setup.

```bash
nextflow run tracker_test_pipeline_2.nf \
--first_fastq /path/to/first_fastq_file \
--second_fastq /path/to/second_fastq_file \
--output_dir /path/to/output_directory \
--adapters /path/to/adapters_file \
--kraken2_db /path/to/kraken2_db \
--metaphlan_db /path/to/metaphlan_db \
--virsorter2_db /path/to/virsorter2_db

## Example Workflow

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.first_fastq = ''
params.second_fastq = ''
params.output_dir = ''
params.adapters = ''
params.kraken2_db = ''
params.metaphlan_db = ''
params.virsorter2_db = ''


process fastQC {

    input:
    path first_fastq
    path second_fastq

    script:
    """
    mkdir -p ${params.output_dir}/fastqc_output
    fastqc -o ${params.output_dir}/fastqc_output ${first_fastq} ${second_fastq}
    """
}

process Trimmomatic {

    input:
    path first_fastq
    path second_fastq
    path adapters

    output:
    path("trimmomatic_output/forward_paired.fq.gz"), emit: forward_paired
    path("trimmomatic_output/reverse_paired.fq.gz"), emit: reverse_paired

    script:
    """
    mkdir -p trimmomatic_output
    trimmomatic PE ${first_fastq} ${second_fastq} \\
    trimmomatic_output/forward_paired.fq.gz trimmomatic_output/forward_unpaired.fq.gz \\
    trimmomatic_output/reverse_paired.fq.gz trimmomatic_output/reverse_unpaired.fq.gz \\
    ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    cp -r trimmomatic_output ${params.output_dir}/
    """
}

process Kraken2 {

    input:
    path forward_paired
    path reverse_paired
    path kraken2_db

    script:
    """
    mkdir -p ${params.output_dir}/kraken2_output
    kraken2 --db ${kraken2_db} --paired ${forward_paired} ${reverse_paired} --report ${params.output_dir}/kraken2_output/kraken2_report.txt
    """
}

process MetaPhlAn4{

	input:
	path forward_paired
	path reverse_paired
	path metaphlan_db
	
	script:
	"""
	mkdir -p ${params.output_dir}/metaphlan4_output
        metaphlan ${forward_paired},${reverse_paired} --input_type fastq --bowtie2db ${metaphlan_db} --bowtie2out ${params.output_dir}/metaphlan4_output/metaphlan_bowtie2.txt > ${params.output_dir}/metaphlan4_output/profiled_metagenome.txt
	"""
}



process MegaHit {
	
	input:
	path forward_paired
	path reverse_paired



	output:
	path("megahit_output_${task.process.hashCode()}/final.contigs.fa"), emit: contigs



	script:
	"""
	mkdir -p ${params.output_dir}/megahit_output
	megahit -1 ${forward_paired} -2 ${reverse_paired} -o megahit_output_${task.process.hashCode()}
	cp megahit_output_${task.process.hashCode()}/final.contigs.fa ${params.output_dir}/megahit_output/
	"""
	}
	
process MetaBAT2 {

    	input:
    	path contigs

    	output:
    	path("./metabat2_output"), emit: bins

    	script:
    	"""
    	mkdir -p metabat2_output
    	metabat2 -i ${contigs} -o ./metabat2_output/bin
    	cp -r ./metabat2_output ${params.output_dir}/
    	"""
}

process QUAST {

    input:
    path contigs

    output:
    path("./quast_output"), emit: quast_report

    script:
    """
    mkdir -p quast_output
    quast ${contigs} -o ./quast_output
    cp -r ./quast_output ${params.output_dir}/
    """
}

process AMRFinderPlus {
	container 'staphb/ncbi-amrfinderplus:latest' 
	
	publishDir "${params.output_dir}/amr_finder_output", mode: 'copy'
	
	input:
	path contigs
	
	output:
	path("amr_finder_output/*"), emit: amr_results
	
	script:
	"""
	mkdir -p amr_finder_output
	amrfinder -n ${contigs} -o amr_finder_output/amr_report.txt
	"""
}


process plasmidFinder{
	container 'staphb/plasmidfinder:latest'
	
	publishDir "${params.output_dir}/plasmid_finder_output", mode: 'copy'
	
	input:
	path contigs
	
	output:
	path("plasmid_finder_output/*"), emit: plasmidFinder_results
	
	script:
	"""
	mkdir -p plasmid_finder_output
	plasmidfinder.py -i ${contigs} -o plasmid_finder_output
	"""
}

process VirSorter2 {
	container 'staphb/virsorter2:latest'
	publishDir "${params.output_dir}/virsorter2_output", mode: 'copy'

	input:
	path contigs
	path virsorter2_db

	output:
	path("virsorter2_output/*"), emit: virsorter2_results


	script:
	"""
	mkdir -p virsorter2_output
	virsorter run -i ${contigs} -w virsorter2_output --provirus-off --max-orf-per-seq 20 all -j 12 --db-dir ${virsorter2_db}
	"""
}
workflow {
    fastQC(params.first_fastq, params.second_fastq)
    trimmed_reads = Trimmomatic(params.first_fastq, params.second_fastq, params.adapters)
    Kraken2(trimmed_reads.forward_paired, trimmed_reads.reverse_paired, params.kraken2_db)
    MetaPhlAn4(trimmed_reads.forward_paired, trimmed_reads.reverse_paired, params.metaphlan_db)
    assembly = MegaHit(trimmed_reads.forward_paired, trimmed_reads.reverse_paired)
    MetaBAT2(assembly.contigs)
    QUAST(assembly.contigs)
    AMRFinderPlus(assembly.contigs)
    plasmidFinder(assembly.contigs)
    VirSorter2(assembly.contigs, params.virsorter2_db)
}












