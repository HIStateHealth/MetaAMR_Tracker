#!/bin/bash

# Check if the correct number of arguments are passed
## USers have to put the location of their fastq files....
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <path_to_first_fastq> <path_to_second_fastq> <input_directory> <output_directory>"
    exit 1
fi

# Assign command line arguments to variables
first_fastq="$1"
second_fastq="$2"
input_dir="$3"
output_dir="$4"

##Running fastQC....
echo "FastQC program is running now..."
fastqc -o ${output_dir}/fastqc_output $first_fastq $second_fastq

## Running Trimmomatic program to trim reads and adapters off...
echo "Trimmomatic Program is Running...."
trimmomatic PE $first_fastq $second_fastq   ${output_dir}/trimmed_$(basename $first_fastq .fastq)_paired.fq ${output_dir}/trimmed_$(basename $first_fastq .fastq)_unpaired.fq   ${output_dir}/trimmed_$(basename $second_fastq .fastq)_paired.fq ${output_dir}/trimmed_$(basename $second_fastq .fastq)_unpaired.fq   ILLUMINACLIP:/path/to/TruSeq-PE.fa:2:30:10:2:true   LEADING:3   TRAILING:3   MINLEN:36

echo "Trimmomatic has finished processing samples $(basename $first_fastq) and $(basename $second_fastq)."

### Next is taxonomy classification....
## We want to run Metaphlan and Kraken2 in parallel...
echo "Running Metaphla4 and Kraken2..."
metaphlan $first_fastq $second_fastq   --input_type fastq > ${output_dir}/metaphlan_output.txt &

kraken2 --db /path/to/kraken2db   --paired $first_fastq $second_fastq   --output ${output_dir}/kraken2_output.txt &

## wait for them to finish...
wait

## Now we run Megahit for assembly...
echo "...running Megahit now for assembly...."
megahit -1 ${output_dir}/trimmed_$(basename $first_fastq .fastq)_paired.fq -2 ${output_dir}/trimmed_$(basename $second_fastq .fastq)_paired.fq   -o ${output_dir}/megahit_output

## Run Abricate here....
echo "...running Abricate now..."
abricate --db card --threads 4 ${output_dir}/megahit_output/final.contigs.fa > ${output_dir}/abricate_output.txt

PlasmidVerify -i ${output_dir}/megahit_output/final.contigs.fa -o ${output_dir}/plasmidverify_output

echo  "...running VirSorter2..."
virsorter run -w ${output_dir}/virsorter_output -i ${output_dir}/megahit_output/final.contigs.fa --min-length 1500 -j 4 all

echo "Pipeline finished...."







