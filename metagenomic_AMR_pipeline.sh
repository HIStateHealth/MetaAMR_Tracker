#!/bin/bash


## This is the input and also the Output directories to where
## everything will go...

input_dir = ""
output_dir = ""

##Running fastQC....
echo "FastQC program is running now..."

fastqc -o ${output_dir}/fastqc_output ${input_dir}/DRR274968_1.fastq/DRR274968_2.fastq

## Running Trimmomatic program to trim reads and adapters off...

echo "Trimmomatic Program is Running...."

trimmomatic PE ${input_dir}/DRR274968_1.fastq ${input_dir}/DRR274968_2.fastq \
  ${output_dir}/trimmed_DRR274968_1_paired.fq ${output_dir}/trimmed_DRR274968_1_unpaired.fq \
  ${output_dir}/trimmed_DRR274968_2_paired.fq ${output_dir}/trimmed_DRR274968_2_unpaired.fq \
  ILLUMINACLIP:/path/to/TruSeq-PE.fa:2:30:10:2:true \
  LEADING:3 \
  TRAILING:3 \
  MINLEN:36


echo "Trimmomatic has finished processing samples DRR274968_1.fastq and DRR274968_2.fastq."


### Next is taxonomy classification....
## We want to run Metaphlan and Kraken2 in parallel...
echo "Running Metaphla4 and Kraken2..."

metaphlan ${input_dir}/DRR274968_1.fastq ${input_dir}/DRR274968_2.fastq \
  --input_type fastq > ${output_dir}/metaphlan_output.txt &

kraken2 --db /path/to/kraken2db \
  --paired ${input_dir}/DRR274968_1.fastq ${input_dir}/DRR274968_2.fastq \
  --output ${output_dir}/kraken2_output.txt &

## wait for them to finish...

wait

## Now we run Megahit for assembly...

echo "...running Megahit now for assembly...."

megahit -1 ${output_dir}/trimmed_DRR274968_1_paired.fq -2 ${output_dir}/trimmed_DRR274968_2_paired.fq \
  -o ${output_dir}/megahit_output

## Run Abricate here....

echo "...running Abricate now..."

abricate --db card --threads 4 ${output_dir}/megahit_output/final.contigs.fa > ${output_dir}/abricate_output.txt

PlasmidVerify -i ${output_dir}/megahit_output/final.contigs.fa -o ${output_dir}/plasmidverify_output


echo  "...running VirSorter2..."

virsorter run -w ${output_dir}/virsorter_output -i ${output_dir}/megahit_output/final.contigs.fa --min-length 1500 -j 4 all

echo "Pipeline finished...."







