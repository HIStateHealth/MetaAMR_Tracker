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
process Metaphlan4 {

        input:
        file trimmed_reads from trimmed_reads.colect()

        output:
        file  "metaphlan_output.txt" into metaphlan_results

        script:
        """
        metaphlan ${params.first_fastq} ${params.second_fastq} --input_type fastq > ${params.output_dir}/metaphlan_output.txt
        """

}

process Kraken2 {
        input:
        file trimmed_reads from trimmed_reads.collect()

        output:
        file "kraken2output.txt into kraken2_results

        script:

        """
        kraken2 --db /path to kraken2db.fa --paired ${params.first_fastq} ${params.second_fastq} --output ${params.output_dir}/kraken2_output.txt
        """
        }

process Megahit {
        input:
        file trimmed_reads from trimmed_reads.collect()

        output:
        file 'megahit_output' into megahit_results

        script:
        """
        megahit -1 ${params.first_fastq} -2 ${params.second_fastq} -o ${params.output_dir}/megahit_output
        """
}
process QUAST {
        input:
        file megahit_output from megahit_results

        output:
         file 'quast_output' into quast_results

        scripts:
        """
        quast ${params.output_dir}/megahit_output -o ${params.output_dir}/quast_output
        """

}

process AMRFinder {

        input:
        file megahit_output from megahit_results

        output:
        file 'amr_output.txt' into amr_results

        script:
        """
        plasmidfinder -i ${params.output_dir}/megahit_output -o ${params.output_dir}/plasmidfinder_output
}

process VirSorter2 {
          input:
           file megahit_output from megahit_results

           output:
           file 'virsorter_output' into virsorter2_results

            script:
            """
            virsorter run -w ${params.output_dir}/virsorter_output -i ${params.output_dir}/megahit_output --min-length 1500 -j 4 all

            """
}

workflow  {
    fastQC()
    Trimmomatic()
    Metaphlan4()
    Kraken2()
    Megahit()
    QUAST()
    AMRFinder()
    PlasmidFinder()
    VirSorter2()
}
