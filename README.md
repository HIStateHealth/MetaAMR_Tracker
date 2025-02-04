# MetaAMR_Tracker

MetaAMR-Tracker is a comprehensive metagenomics pipeline designed for tracking Antimicrobial Resistance (AMR) genes and fungi genes in environmental wastewater samples. This multifaceted pipeline supports processing of Illumina short reads, Nanopore long reads, and hybrid reads, catering to a diverse range of sequencing technologies.

Users must specify the type of reads (Illumina, Nanopore, or Hybrid) they wish to process in the UI frontend interface. This input is crucial for the pipeline to select the appropriate path and tools for analysis.

From quality control to contig assembly and resistance gene identification, MetaAMR-Tracker integrates state-of-the-art tools like FastQC, MEGAHIT, MetaPhlAn, and ResFinder. The pipeline delivers a holistic view of the microbial resistome, providing valuable insights into the prevalence and distribution of AMR genes in environmental samples.



![image](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/0f595650-5f7e-403d-b6b3-a4c25f85e44d)



![image](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/f834b920-49aa-4d3e-9dec-4e2fcb19e25c)



![image](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/df06b4d7-76af-40c0-814e-65d983f45890)







![pipeline](https://github.com/HIStateHealth/MetaAMR_Tracker/assets/138935158/ca9602d0-0d8f-44e6-ad3d-5764e41adeb2)




![image](https://github.com/user-attachments/assets/aa8467ff-9c88-4e81-aa3d-94a352e83898)





# Metagenomic Analysis Pipeline with Nextflow

## Introduction
This pipeline serves as a comprehensive solution for metagenomic data analysis. It is designed to run in a Nextflow environment, leveraging Docker containers for seamless software management. The pipeline incorporates a variety of bioinformatics tools for tasks ranging from quality control to specialized analyses like antibiotic resistance and viral identification.

## Prerequisites
Before you begin, ensure that you have installed the following:
- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/products/docker-desktop)
- Python 3.0 or higher
- git
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

# Metagenomic Analysis Programs for Antimicrobial Gene Surveillance ILLUMINA SIDE****

## FastQC
FastQC provides a quick and intuitive overview of data quality, assessing sequencing reads for various quality metrics, including base quality, sequence duplication levels, and overrepresented sequences. This step is crucial for identifying potential issues early in the pipeline, ensuring that subsequent analyses are based on high-quality data.

## Trimmomatic
Trimmomatic performs read trimming and filtering, removing adapters and low-quality bases. This step is essential for eliminating sequencing artifacts that could interfere with downstream analyses, such as assembly and gene prediction, thus improving the accuracy of antimicrobial resistance gene detection.

## MegaHit
MegaHit assembles high-throughput sequencing reads into longer contiguous sequences (contigs). In metagenomics, it enables the reconstruction of microbial genomes from mixed community samples. This step is foundational for identifying the microbial composition and potential antimicrobial resistance genes within a sample.

## MetaQUAST
MetaQUAST evaluates the quality of metagenomic assemblies, providing metrics like N50, L50, and the number of misassemblies, offering insights into the assembly's completeness and accuracy. High-quality assemblies are critical for reliable identification and annotation of genes related to antimicrobial resistance.

## Metaphlan4
Metaphlan4 specializes in profiling the composition of microbial communities from metagenomic sequencing data, using unique clade-specific marker genes for precise taxonomic assignment. This information is pivotal for understanding microbial diversity and detecting potential reservoirs of antimicrobial resistance.

## Kraken2
Kraken2 rapidly assigns taxonomic labels to short DNA sequences, offering insights into the microbial composition of a sample. This step is vital for antimicrobial gene surveillance as it helps identify which microbes carry resistance genes, contributing to our understanding of resistance spread and evolution.

## Metawrap (using Maxbin2, Metabat2, CONCOCT)
Metawrap integrates tools like Maxbin2, Metabat2, and CONCOCT for metagenomic binning, grouping contigs into bins representing potential genomes. This process is crucial for resolving community structure and identifying microbes that may harbor antimicrobial resistance genes.

## DAS TOOL
DAS TOOL integrates binning results from multiple sources to optimize the quality and completeness of metagenomic bins. By selecting the best bins across different algorithms, it enhances the accuracy of identifying microbes and their antimicrobial resistance genes within complex samples.

## CheckM
CheckM assesses metagenomic bins for completeness and contamination, ensuring that genomic reconstructions used in antimicrobial resistance gene analysis are reliable.

## GTDB-tk
GTDB-tk assigns taxonomy to microbial genomes based on the Genome Taxonomy Database (GTDB), facilitating accurate tracking of antimicrobial resistance genes across different taxa.

## Prepend Python Program

The Prepend Program is designed to enhance the traceability and manageability of contigs across different bins in metagenomic analyses. It does so by modifying the identifiers (IDs) of the contigs to include the name of the bin they belong to. This modification aids in downstream analyses where the origin of each contig is crucial for understanding the distribution of features, such as antibiotic resistance genes or specific metabolic capabilities, across the microbial communities present in a sample.

## How It Works
The program takes as input:
- A directory containing FASTA files, where each file represents a bin (a collection of contigs belonging to a putative genome) produced by binning tools like CONCOCT, MaxBin2, or MetaBAT2.
- The desired output directory for the modified FASTA files.

For each FASTA file (representing a bin), the program reads through the file and prepends the bin name to each contig's ID. For example, if a contig originally has the ID `>k141_393177` and it's contained in a file named `CONCOCT_bin.3.fa`, the modified ID will be `>CONCOCT_bin.3_k141_393177`. This process is repeated for all contigs across all bins.

## Output
The output of the Prepend Program is a set of modified FASTA files, one for each input bin, where the IDs of the contigs are enriched with their corresponding bin names. These modified files are stored in the specified output directory. Optionally, the program can merge all the modified bins into a single FASTA file, facilitating analyses that require the combined dataset rather than individual bins. This merged file simplifies the input for various genomic analysis tools that may not support multi-file input or that benefit from a holistic view of the metagenome.

This approach ensures that each contig can be easily traced back to its original bin, enhancing the interpretability of results from downstream analytical processes.


# Kaiju

Kaiju is a bioinformatics program designed for the taxonomic classification of metagenomic sequences. It uses a protein-level alignment approach, making it particularly effective for identifying and classifying microbial species in environmental samples, including those with incomplete or fragmented genomes.

## How It Works
- **Protein-Level Alignment**: Kaiju translates DNA sequences into protein sequences and compares them against a reference database of known protein sequences. This approach improves accuracy in identifying microbial species, especially for sequences that are not well-represented in nucleotide databases.
- **Reference Database**: Kaiju uses a reference database that includes protein sequences from various microbial species. This database can be customized to include specific groups of interest, such as pathogens or antibiotic resistance genes.
- **Classification**: Kaiju assigns taxonomic labels to sequences based on the best matches found in the reference database. It can identify bacteria, archaea, fungi, viruses, and other microorganisms.



## Applications in AMR Detection
- **Microbial Identification**: Kaiju helps identify the presence of specific microbial species in metagenomic samples. This is crucial for understanding the microbial community structure and detecting potential pathogens.
- **Resistance Gene Classification**: By identifying the species carrying antibiotic resistance genes, Kaiju provides insights into the spread and origin of these genes within microbial communities.
- **Environmental Monitoring**: Kaiju is used in monitoring environmental samples such as wastewater, soil, and water to detect and track the presence of antimicrobial resistance genes and their host organisms.

## Example Usage
```bash
# Run Kaiju for taxonomic classification
kaiju -t nodes.dmp -f kaiju_db.fmi -i input.fastq -o output.kaiju

# Generate a summary report
kaijuReport -t nodes.dmp -n names.dmp -i output.kaiju -o report.out


## AMR Finder Plus and Enhanced Python Script Integration

### AMR Finder Plus
This essential tool scans our samples, matching sequences against a rigorously maintained database to spot antimicrobial resistance genes. It's a whiz at detecting not just the genes themselves but also their variants and the mechanisms behind their resistance superpowers. Think of AMR Finder Plus as our detective, uncovering the hidden secrets of microbial resistance within our communities.

### Python Script Enhancements
We've boosted our AMR process by coupling it with a suite of customized Python scripts. These scripts have been fine-tuned to work seamlessly with AMR Finder Plus, automating and refining the steps from gene extraction to visualization.

- **Sequence Extraction**: Our `extract_amr_contigs.py` script meticulously carves out the AMR gene sequences from larger genetic canvases, setting the stage for precise analysis.
- **KMA Alignment**: With the power of the KMA tool, we align raw sequence data to our focused library of AMR genes, ensuring that every potential match is explored and noted.
- **Abundance Calculation**: The `kma_to_amr_abundance_for_amrfindergenes.py` script takes the baton and calculates how much of each resistance gene is present, giving us quantitative insight into their prevalence.
- **Data Integration**: We then bring everything together, merging our findings into the AMR Finder Plus reports using `update_amr_finder_report_with_abundance.py`, creating a detailed map of resistance across our samples.
- **Visualization with Heatmaps**: Finally, `heatmap_from_amr_finder_report_abundance.py` transforms our data into vivid heatmaps, making complex patterns of resistance easily digestible at a glance.

In short, we've not only implemented AMR Finder Plus but also customized additional steps that extract, align, quantify, and visualize AMR genes. This fortified process should give the user with a more comprehensive understanding of antimicrobial resistance.


## Virsorter2
Virsorter2 detects viral sequences within metagenomic data, including prophages within microbial genomes. Understanding these elements is crucial for grasping horizontal gene transfer mechanisms, which can contribute to the spread of antimicrobial resistance.

## Plasmidfinder
Plasmidfinder identifies plasmid sequences in bacterial genomes, playing a pivotal role in understanding and monitoring the dissemination of resistance within and across microbial communities, as plasmids are key vectors for the horizontal transfer of antimicrobial resistance genes.

NANOPORE SIDE****




# Nextflow Pipeline Guide

## Overview
This guide provides instructions on setting up and running a Nextflow pipeline for metagenomic analysis. The pipeline uses Docker containers for each process, ensuring consistency and reproducibility.
All Docker images included in the Nextflow code are available and can be searched on "Docker Hub" located here: **https://hub.docker.com/**  

## Pipeline Setup

### 1. Nextflow Configuration
Create a `nextflow.config` file in your pipeline's root directory  OR download the config file here in the github repository above with the following content to enable Docker and set the DSL version:



```docker {
enabled = true
}

//executor {
  // queueSize = 1
//}

process {
    memory = '56 GB'
}

params {


### DOWNLOAD ALL NECESSSARY DATABASE FILES AND KEEP THEM IN THE DATABASES FOLDER STORED LOCALLY ON YOUR SYSTME SOMEWHERE...

**https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa**
adapters = "${baseDir}/databases/TruSeq3-PE.fa"

**https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20240605.tar.gz**
kraken2_db = "${baseDir}/databases/k2_plus_small_db_8gb"


virsorter2_db = "${baseDir}/databases/virsorter2_DOCKER/db"

**https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz**
gtdbtk_db = "${baseDir}/databases/gtdbtk_r214_data/release214"


}


## Running the Pipeline

Ensure Docker is installed and configured with a Nextflow config file.

To run the pipeline, use the following command structure, replacing `<placeholders>` with appropriate paths or values:


## COMMAND LINE RUN EXAMPLE...


##################  *******************  Command line for ILLUMINA  **************************##########################


nextflow run tracker_test_pipeline_2.nf --first_fastq /home/user01/milestones_6/metagenomics_custom_pipeline/fastq_files/DRR274968_1.fastq.gz --second_fastq /home/user01/milestones_6/metagenomics_custom_pipeline/fastq_files/DRR274968_2.fastq.gz --output_dir /home/user01/milestones_6/metagenomics_custom_pipeline/nextflow_workflow/Meta_Pipeline/MetaAMR_Tracker/Nextflow_output_directory -resume



##################  *******************  Command line for NANOPORE  **************************##########################



nextflow run NANOPORE_Nexflow_Pipeline.nf --nanopore_fastq /home/user01/milestones_6/CAMISIM_PROGRAM/genomes_fasta_files/NANO_SIM_PROGRAM/NanoSim_1k_OUTPUT_TEST_sample0_aligned_reads.fastq.gz --output_dir /home/user01/milestones_6/metagenomics_custom_pipeline/nextflow_workflow/Meta_Pipeline/MetaAMR_Tracker/NANOPORE_NEXTFLOW_OUTPUT -resume



##################  *******************  Command line for HYBRID  **************************##########################

nextflow run HYBRID_Meta_Tracker_Pipeline.nf --first_fastq /home/user01/milestones_6/metagenomics_custom_pipeline/HYBRID_Pipeline/Illumina/SRR23926929_1.fastq.gz --second_fastq /home/user01/milestones_6/metagenomics_custom_pipeline/HYBRID_Pipeline/Illumina/SRR23926929_2.fastq.gz --output_dir /home/user01/milestones_6/metagenomics_custom_pipeline/nextflow_workflow/Meta_Pipeline/MetaAMR_Tracker/HYBRID_Nextflow_output -resume




```



## Important Note on Nextflow DSL Version

This pipeline is written using Nextflow DSL 2 syntax. To run this pipeline, your environment needs to support Nextflow DSL 2.

Please ensure that your Nextflow configuration includes the line `nextflow.enable.dsl=2`. This line is necessary to enable DSL 2 features in the script. If you are using an older version of Nextflow, you might need to update it to a version that supports DSL 2.

To update or confirm the DSL version in your Nextflow environment, you can create or modify the `nextflow.config` file in your project directory with the following content:

```nextflow
nextflow.enable.dsl=2

```Also you need to make sure to create a 


## Example Command
Here's an example command to help you get started. @@@####****REPLACE THE PATHS WITH THOSE
CORRESPONDING TO YOUR OWN SETUP IF YOU WANT TO MODIFY THE PIPELINE. Or to run this pipeline, just clone this a run it as is with
your different fastq file for it to work, as per instructions in the command line in the previous step....*****####

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



                   #####*****                   Documentation and Performance Metrics for Established Workflows      ######********

## Introduction

This documentation outlines the standard operating procedures (SOPs), periodic updates, upgrades, and performance metrics for the established workflows of our project. Our goal is to ensure transparency, maintain high-quality standards, and continuously improve the efficiency and reliability of our workflows.

## Table of Contents

- [Standard Operating Procedures (SOPs)](#standard-operating-procedures-sops)
- [Periodic Updates and Upgrades](#periodic-updates-and-upgrades)
- [Quantitative Performance Metrics](#quantitative-performance-metrics)
- [Evaluation Reports](#evaluation-reports)

## Standard Operating Procedures (SOPs)

The SOPs for each workflow are documented to provide clear instructions on their execution, maintenance, and troubleshooting. These documents are crucial for onboarding new team members and ensuring consistency in workflow execution.

- **Workflow Documentation**: Located in `docs/workflows/`, each workflow has a dedicated documentation file outlining its purpose, inputs, outputs, and step-by-step instructions.

## Periodic Updates and Upgrades

To ensure our workflows remain up-to-date with the latest software developments and security standards, we have established a schedule for periodic updates and upgrades.

- **Update Schedule**: Updates are scheduled on a quarterly basis, with emergency patches applied as needed.
- **Upgrade Process**: The upgrade process is documented in `docs/updates/UPGRADE_GUIDE.md`, detailing the steps for safely updating packages, tools, and functions.

## Quantitative Performance Metrics

Performance metrics are established to quantitatively evaluate the efficiency, accuracy, and reliability of our workflows. For comprehensive evaluation, we utilize simulated datasets generated by CAMISIM and NanoSim to mimic real-world fastq files. This approach allows us to assess how well the Meta AMR Tracker Pipeline performs under controlled yet realistic conditions.

### Simulation Tools Used:

- **CAMISIM**: Generates simulated metagenomic sequencing data, providing a diverse set of microbial communities. This tool helps us understand the pipeline's capability to handle complex metagenomic samples.
- **NanoSim**: Simulates Nanopore sequencing reads, offering insights into the pipeline's performance with long-read technologies.

### Metrics Employed:

- **Accuracy**: Measures the correctness of the workflow output against the known data from simulated datasets.
- **Efficiency**: Evaluated by the runtime and resource usage, providing insights into the computational demands of the pipeline.
- **Scalability**: Assessed based on the workflow's performance as the data volume and complexity scale.
- **Reliability**: Frequency of failures or errors encountered during workflow execution, ensuring the pipeline's stability.
- **Time Scale and Temporal Analysis**: Evaluates the pipeline's performance over time, considering factors such as data growth and algorithmic efficiency.
- **Space Complexity Scale**: Measures the amount of storage space required by all programs in the pipeline, offering insights into its scalability and efficiency in resource utilization.

### Evaluation Reports

Detailed evaluation reports are generated for each major release, incorporating analysis based on the above metrics. These reports delve into the pipeline's performance across different simulated datasets, highlighting areas for potential improvement.

- **Report Location**: Reports are stored in `docs/evaluation_reports/`, with each report named according to the workflow and version evaluated.
- **Contents**: Each report includes an executive summary, methodology, detailed metric evaluations, and recommendations for future enhancements.

The integration of simulation tools and advanced performance metrics ensures a thorough assessment of the Meta AMR Tracker Pipeline, guiding ongoing optimizations and improvements.



                     



       # **MetaAMR Tracker Pipeline Evaluation Report**




                                                    ## **Overview**
This report outlines the testing and evaluation of the MetaAMR Tracker pipeline, which is designed to handle comprehensive analysis of antimicrobial resistance (AMR) genes from environmental samples. This document serves to demonstrate the pipeline's effectiveness, based on rigorous testing metrics.

                                                ## **Testing Environment**
The pipeline was executed using Nextflow, which provides robust data handling capabilities to manage complex workflows. The pipeline was tested under various computational loads to simulate real-world scenarios and to ensure its performance across different systems.

## **Metrics Employed**
- **Accuracy**: Verification against known datasets to ensure the pipeline's outputs are accurate.
- **Efficiency**: Monitored by runtime and computational resource utilization.
- **Scalability**: Tested by incrementally increasing the data volume to evaluate performance.
- **Reliability**: Assessed through the frequency of failures or errors during executions.
- **Time Scale and Temporal Analysis**: Performance monitored over time to ensure stability and efficiency.
- **Space Complexity Scale**: Evaluated the storage requirements to optimize space utilization.

## **Testing Commands**
The pipeline was executed with the following Nextflow command to enable detailed resource monitoring:

```bash
nextflow run ILLUMINA_Meta_Tracker_Pipeline_tested_final.nf --first_fastq /home/chris/RESEARCH/metagenome/GUI/DRR274968_1.fastq.gz --second_fastq /home/chris/RESEARCH/metagenome/GUI/DRR274968_2.fastq.gz --output_dir /home/chris/RESEARCH/metagenome/GUI/illumina_output_metrics -with-report -with-trace -with-timeline -resume


| Task ID | Process            | Status     | CPU Usage | Peak RAM  | Duration  |
|---------|--------------------|------------|-----------|-----------|-----------|
| 1       | fastQC             | CACHED     | 102.9%    | 562.2 MB  | 5m 10s    |
| 2       | Trimmomatic        | CACHED     | 100.7%    | 1.3 GB    | 21m 23s   |
| 3       | Kraken2            | CACHED     | 107.5%    | 8.2 GB    | 9m 52s    |
| 4       | Metaphlan_reads    | CACHED     | 952.5%    | 18.5 GB   | 28m 37s   |
| 5       | MegaHit            | CACHED     | 4538.5%   | 12.3 GB   | 1h 28m 10s|
| 6       | metawrapBinning    | COMPLETED  | 854.9%    | 9.9 GB    | 1h 3m 46s |
| 7       | QUAST              | COMPLETED  | 127.1%    | 650.6 MB  | 46.9s     |
| 8       | DAStool            | COMPLETED  | 1049.7%   | 2.4 GB    | 4m 6s     |
| 9       | CheckM             | COMPLETED  | 116.5%    | 70.2 GB   | 1h 9m 48s |
| 10      | PrependBinNames    | COMPLETED  | 99.5%     | 1.2 GB    | 4.6s      |
| 11      | GTDBTK             | COMPLETED  | 142.2%    | 729.1 GB  | 47m 57s   |
| 12      | PlasmidFinder      | COMPLETED  | 100.0%    | 4.4 GB    | 1h 36m 59s|
| 13      | VirSorter2         | COMPLETED  | 1468.4%   | 2.2 GB    | 19h 33m 24s|
| 14      | AMRfinder          | COMPLETED  | 238.8%    | 1.2 GB    | 7m 47s    |







