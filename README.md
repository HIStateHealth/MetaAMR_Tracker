# MetaAMR_Tracker

MetaAMR-Tracker is a comprehensive metagenomics pipeline designed for tracking Antimicrobial Resistance (AMR) genes and fungi genes in environmental wastewater samples. This multifaceted pipeline supports processing of Illumina short reads, Nanopore long reads, and hybrid reads, catering to a diverse range of sequencing technologies.

Users must specify the type of reads (Illumina, Nanopore, or Hybrid) they wish to process in the UI frontend interface. This input is crucial for the pipeline to select the appropriate path and tools for analysis.

From quality control to contig assembly and resistance gene identification, MetaAMR-Tracker integrates state-of-the-art tools like FastQC, MEGAHIT, MetaPhlAn, and ResFinder. The pipeline delivers a holistic view of the microbial resistome, providing valuable insights into the prevalence and distribution of AMR genes in environmental samples.


![image](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/f834b920-49aa-4d3e-9dec-4e2fcb19e25c)



![image](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/df06b4d7-76af-40c0-814e-65d983f45890)







![pipeline](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/ca9602d0-0d8f-44e6-ad3d-5764e41adeb2)



![pipeline_drawing_Hybrid drawio (2)](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/c56b1427-08b0-4aba-86fa-69d50eb9554c)





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
- All bioinformatics tools l![Uploading pipeline_drawing_Hybrid.drawio (2).png…]()
isted in the Nextflow script

## Setup
1. Clone the GitHub repository.
    ```bash
    git clone https://github.com/HIStateHealth/MetaAMR_Tracker.git
    ```

2. Navigate into the project directory.
    ```bash
    cd [I will put my Project Directory here]
    ```
## Pipeline Steps ILLUMINA SIDE****

### FastQC
- **Function**: Performs quality control checks on raw sequence data.
- **Container**: `staphb/fastqc:latest`

### Trimmomatic
- **Function**: Trims and removes adapters from sequencing data.
- **Container**: `staphb/trimmomatic:latest`

### MegaHit
- **Function**: Assembles high-quality reads into contigs.
- **Container**: `nanozoo/megahit:latest`

### QUAST
- **Function**: Assesses the quality of genome assemblies.
- **Container**: `staphb/quast:latest`

### Kraken2
- **Function**: Taxonomically classifies sequencing reads.
- **Script**: `kraken2 --db ${kraken2_db} --paired ${forward_paired} ${reverse_paired} --report ${params.output_dir}/kraken2_output/kraken2_report.txt`

### MetaPhlAn
- **Function**: Profiles microbial community composition.
- **Container**: `stang/metaphlan4:v1`

### MetaWRAP (Binning)
- **Function**: Bins contigs into genome bins.
- **Container**: `pnatarajan84/metawrap_1.3.2_binning_fastq_gz_modified:latest`

### DAS Tool
- **Function**: Refines and selects genome bins.
- **Container**: `bladerunner2945/dastools1.1.6:latest`

### CheckM
- **Function**: Evaluates the quality of genome bins.
- **Container**: `nanozoo/checkm:1.1.3--c79a047`

### GTDB-Tk
- **Function**: Taxonomically classifies bins using the GTDB.
- **Container**: `nanozoo/gtdbtk:2.3.2--641ec99`

### PrependBinNames
- **Function**: Prepares bins for downstream analysis.
- **Container**: `bladerunner2945/prepend_bin_names:latest`

### VSEARCH
- **Function**: Classifies fungal taxa.
- **Container**: `manutamminen/vsearch:v1.28`

### AMRfinder
- **Function**: Identifies antimicrobial resistance genes.
- **Container**: `staphb/ncbi-amrfinderplus:latest`

### PlasmidFinder
- **Function**: Identifies plasmids.
- **Container**: `staphb/plasmidfinder:latest`

### VirSorter2
- **Function**: Identifies and classifies viral sequences.
- **Container**: `staphb/virsorter2:latest`



******### <u>This Markdown documentation is formatted to clearly present each step of your pipeline, including the function and container used. The "Running the Pipeline" section provides users with instructions on how to execute the pipeline.</u>


# Nextflow Pipeline Guide

## Overview
This guide provides instructions on setting up and running a Nextflow pipeline for metagenomic analysis. The pipeline uses Docker containers for each process, ensuring consistency and reproducibility.

## Pipeline Setup

### 1. Nextflow Configuration
Create a `nextflow.config` file in your pipeline's root directory with the following content to enable Docker and set the DSL version:

```nextflow
nextflow.enable.dsl=2

docker {
    enabled = true
}


## Running the Pipeline

Ensure Docker or Singularity is installed and configured.

To run the pipeline, use the following command structure, replacing `<placeholders>` with appropriate paths or values:

```bash
nextflow run [Your Pipeline Script Name].nf --first_fastq <path_to_first_fastq> --second_fastq <path_to_second_fastq> --output_dir <path_to_output_directory> --adapters <path_to_adapters_file> --kraken2_db <path_to_kraken2_db> --utax_reference_db <path_to_utax_db> --virsorter2_db <path_to_virsorter2_db> --gtdbtk_db <path_to_gtdbtk_db>


To run the pipeline, use the following command structure, replacing <placeholders> with appropriate paths or values:

// Define your processes here
// Example:
process fastQC {
    container 'staphb/fastqc:latest'
    // ... rest of the fastQC process
}
// Add other processes similarly

// Define the workflow
workflow {
    // Your workflow steps go here
}


## Important Note on Nextflow DSL Version

This pipeline is written using Nextflow DSL 2 syntax. To run this pipeline, your environment needs to support Nextflow DSL 2.

Please ensure that your Nextflow configuration includes the line `nextflow.enable.dsl=2`. This line is necessary to enable DSL 2 features in the script. If you are using an older version of Nextflow, you might need to update it to a version that supports DSL 2.

To update or confirm the DSL version in your Nextflow environment, you can create or modify the `nextflow.config` file in your project directory with the following content:

```nextflow
nextflow.enable.dsl=2

```Also you need to make sure to create a 


## Example Command
Here's an example command to help you get started. @@@####****REPLACE THE PATHS WITH THOSE
CORRESPONDING TO YOUR OWN SETUP, OR IT WILL NOT WORK....*****####

```nextflow
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.first_fastq = ''
params.second_fastq = ''
params.output_dir = ''
params.adapters = ''
params.kraken2_db = ''
params.utax_reference_db = ''
params.virsorter2_db = ''
params.gtdbtk_db = ''

process fastQC {
    container 'staphb/fastqc:latest'
    publishDir "${params.output_dir}/fastqc_output", mode: 'copy'

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
    container 'staphb/trimmomatic:latest'

    input:
    path first_fastq
    path second_fastq
    path adapters

    output:
    path("trimmomatic_output/paired_R1.fastq.gz"), emit: forward_paired
    path("trimmomatic_output/paired_R2.fastq.gz"), emit: reverse_paired

    script:
    """
    mkdir -p trimmomatic_output
    trimmomatic PE ${first_fastq} ${second_fastq} \\
    trimmomatic_output/paired_R1.fastq.gz trimmomatic_output/forward_unpaired.fq.gz \\
    trimmomatic_output/paired_R2.fastq.gz trimmomatic_output/reverse_unpaired.fq.gz \\
    ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    cp -r trimmomatic_output ${params.output_dir}/
    """
}

process MegaHit {
    container 'nanozoo/megahit:latest'

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

process QUAST {
    container 'staphb/quast:latest'

    input:
    path contigs

    output:
    path("./quast_output"), emit: quast_report

    script:
    """
    mkdir -p quast_output
    metaquast.py ${contigs} -o ./quast_output
    cp -r ./quast_output ${params.output_dir}/
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

process Metaphlan_reads {
    container 'stang/metaphlan4:v1'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path forward_paired
    path reverse_paired

    output:
    path "metaphlan_profile_reads.txt", emit: metaphlan_profile_reads
    path "metaphlan_bowtie2.txt", emit: bowtie2_output

    script:
    """
    metaphlan ${forward_paired},${reverse_paired} --input_type fastq --bowtie2out metaphlan_bowtie2.txt --nproc 10 > metaphlan_profile_reads.txt
    """
}

process metawrapBinning {
    container 'pnatarajan84/metawrap_1.3.2_binning_fastq_gz_modified:latest'
    publishDir "${params.output_dir}/metawrap_output", mode: 'copy'

    input:
    path contigs
    path forward_paired
    path reverse_paired

    output:
    path "metawrap_output/*", emit: metawrap_output_files
    path "metawrap_output/concoct_bins/*_bins_contigs.tsv", emit: concoct_bins_tsv
    path "metawrap_output/maxbin2_bins/*_bins_contigs.tsv", emit: maxbin2_bins_tsv
    path "metawrap_output/metabat2_bins/*_bins_contigs.tsv", emit: metabat2_bins_tsv

    script:
    """
    metawrap binning -o metawrap_output -t 40 -a ${contigs} --metabat2 --maxbin2 --concoct ${forward_paired} ${reverse_paired}
    for dir in metawrap_output/concoct_bins metawrap_output/maxbin2_bins metawrap_output/metabat2_bins; do
        if [ -d "\$dir" ]; then
           bin_tool=\$(basename \$dir)
           tsv_file="\$dir/\${bin_tool}_bins_contigs.tsv"
           for bin_file in "\$dir"/*.fa; do
               bin_name=\$(basename \$bin_file .fa);
               grep "^>" \$bin_file | sed "s/>//" | awk -v bin="\$bin_name" '{print \$1 "\t" bin}' >> \$tsv_file;
           done;
       fi;
    done
    """
}

process DAStool {
    container 'bladerunner2945/dastools1.1.6:latest'
    publishDir "${params.output_dir}/dastool_output", mode: 'copy'

    input:
    path concoct_bins_tsv
    path maxbin2_bins_tsv
    path metabat2_bins_tsv
    path contigs 

    output:
    path "dastool_output_DASTool_bins", emit: dastool_output_bins

    script:
    """
    # Ensure the output directory exists
    mkdir -p dastool_output

    # Run DAS_Tool
    DAS_Tool -i ${concoct_bins_tsv},${maxbin2_bins_tsv},${metabat2_bins_tsv} -l concoct,maxbin2,metabat2 -c ${contigs} -o dastool_output --write_bin_evals --threads 20 --write_bins --write_unbinned

    # Debugging: Check the output files
    echo "Checking output files in dastool_output directory:"


workflow {
    fastQC(params.first_fastq, params.second_fastq)
    trimmed_reads = Trimmomatic(params.first_fastq, params.second_fastq, params.adapters)
    assembly = MegaHit(trimmed_reads.forward_paired, trimmed_reads.reverse_paired)
    quast_report = QUAST(assembly.contigs)
    kraken2_report = Kraken2(trimmed_reads.forward_paired, trimmed_reads.reverse_paired, params.kraken2_db)
    binning_results = metawrapBinning(assembly.contigs, trimmed_reads.forward_paired, trimmed_reads.reverse_paired)
    dastool_output_bins = DAStool(binning_results.concoct_bins_tsv, binning_results.maxbin2_bins_tsv, binning_results.metabat2_bins_tsv, assembly.contigs)
    checkm_output = CheckM(dastool_output_bins)
    metaphlan_profile_reads = Metaphlan_reads(trimmed_reads.forward_paired, trimmed_reads.reverse_paired)
    merged_bins = PrependBinNames(dastool_output_bins)
    vsearch_output = Vsearch(merged_bins, params.utax_reference_db)
    amrfinder_output = AMRfinder(merged_bins)
    plasmidfinder_output = PlasmidFinder(merged_bins)
    virsorter2_output = VirSorter2(merged_bins, params.virsorter2_db)
    gtdbtk_output = GTDBTK(dastool_output_bins, params.gtdbtk_db)
}









