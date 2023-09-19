#!/usr/bin/env nextflow

params.first_fastq = "path to file/first_fastq"
params.second_fastq = "path to file/second_fastq"
params.output_dir = "path to output/output_dir"

process fastQC {
        output:
        file "*_fastqc.zip" into fastqc_results

        script:
        """
        fastqc -o ${params.output_dir}/fastqc_output ${params.first_fastq} ${params.second_fastq}

        """
     }


process Trimmomatic {

       input:
        file fastqc_results from fastqc_results.collect()

        output:
        file "*_paired.fq" into trimmed_reads}

        script:
        """
            trimmomatic PE ${params.first_fastq} ${params.second_fastq} ${params.output_dir}/trimmed_first_paired.fq ${params.output_dir}/trimmed_first_unpaired.fq \
             ${params.output_dir}/trimmed_second_paired.fq ${params.output_dir}/trimmed_second_unpaired.fq \
             ILLUMINACLIP:/path/to/TruSeq-PE.fa:2:30:10:2:true LEADING:3 TRAILING:3 MINLEN:36

                """
