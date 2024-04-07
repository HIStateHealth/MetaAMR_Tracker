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
    publishDir "${params.output_dir}/trimmomatic_output", mode: 'copy'

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
    """
}

process MegaHit {

        container 'nanozoo/megahit:latest'
        publishDir "${params.output_dir}/megahit_output", mode: 'copy'
	
	input:
	path forward_paired
	path reverse_paired



	output:
	path("megahit_output_${task.process.hashCode()}/final.contigs.fa"), emit: contigs



	script:
	"""
	mkdir -p megahit_output
	megahit -1 ${forward_paired} -2 ${reverse_paired} -o megahit_output_${task.process.hashCode()}
	"""
	}
	
process QUAST {

    container 'staphb/quast:latest'
    publishDir "${params.output_dir}/quast_output", mode: 'copy'

    input:
    path contigs

    output:
    path("./quast_output"), emit: quast_report

    script:
    """
    mkdir -p quast_output
    metaquast.py ${contigs} -o ./quast_output
    """
}	

process Kraken2 {

    container 'staphb/kraken2:latest'
    publishDir "${params.output_dir}/kraken2_output", mode: 'copy'


    input:
    path forward_paired
    path reverse_paired
    path kraken2_db
    
    output:
    path "kraken2_output", emit: kraken2_report

    script:
    """
    mkdir -p kraken2_output
    kraken2 --db ${kraken2_db} --paired ${forward_paired} ${reverse_paired} --report kraken2_output/kraken2_report.txt
    """
}

process Metaphlan_reads {
    
    container 'stang/metaphlan4:v1'
    publishDir "${params.output_dir}/metaphlan_output", mode: 'copy'

    input:
    path forward_paired
    path reverse_paired

    output:
    path "metaphlan_output/metaphlan_profile_reads.txt", emit: metaphlan_profile_reads
    path "metaphlan_output/metaphlan_bowtie2.txt", emit: bowtie2_output

    script:
    """
    mkdir -p metaphlan_output
    metaphlan ${forward_paired},${reverse_paired} --input_type fastq --bowtie2out metaphlan_output/metaphlan_bowtie2.txt --nproc 10 > metaphlan_output/metaphlan_profile_reads.txt
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
    ls dastool_output
    """
}

process CheckM {
    
    container 'nanozoo/checkm:1.1.3--c79a047'
    publishDir "${params.output_dir}/checkm_output", mode: 'copy'

    input:
    path dastool_output_bins

    output:
    path "checkm_output/*", emit: checkm_output_files

    script:
    """
    mkdir -p checkm_output
    checkm lineage_wf -x fa ${dastool_output_bins} checkm_output
    """
}


process GTDBTK {

	container 'nanozoo/gtdbtk:2.3.2--641ec99'
	publishDir "${params.output_dir}/gtdbtk_output", mode: 'copy'

	input:
	path dastool_output_bins
	path gtdbtk_db

	output:
	path "gtdbtk_output/*", emit: gtdbtk_output_files

	script:
	"""
	export GTDBTK_DATA_PATH=${gtdbtk_db}
	gtdbtk classify_wf --genome_dir ${dastool_output_bins} --out_dir gtdbtk_output --skip_ani_screen --cpus 20 --extension fa
	"""
	}

process PrependBinNames {
    
    container 'bladerunner2945/prepend_bin_names:latest'
    publishDir "${params.output_dir}/prepend_bin_names_output", mode: 'copy'

    input:
    path dastool_output_bins

    output:
    path "merged_bins.fa", emit: merged_bins_fa

    script:
    """
    prepend_bin_names.py ${dastool_output_bins} merged_bins.fa
    """
}

process Vsearch  {

	container 'manutamminen/vsearch:v1.28'
        publishDir "${params.output_dir}/vsearch_output", mode: 'copy'
	
	input:
	path merged_bins_fa
	path utax_reference_db
	
	output:
	path "output_fungi_taxa_bins.txt", emit: vsearch_fungi_output
	
	script:
	"""
	vsearch --usearch_global ${merged_bins_fa} --db ${utax_reference_db} --id 0.97 --blast6out output_fungi_taxa_bins.txt --top_hits_only
        """}


process AMRfinder {
	container 'pnatarajan84/amr_kma_fpkm_heatmap:latest'
	publishDir "${params.output_dir}/amrfinder_output", mode: 'copy'

	input:
	path merged_bins_fa
	path forward_paired
	path reverse_paired

	output:
	path "amrfinder_output/*", emit: amr_finder_output_files

	script:
	"""
	mkdir -p amrfinder_output
	amrfinder -n ${merged_bins_fa} -o amrfinder_output/amr_finder_report.txt
	kma index -i /data/AMR_CDS.fasta -o amr_kma_index

	# Execute kma and capture its exit status
	kma -ipe ${forward_paired} ${reverse_paired} -o amr_kma_OUTPUT -t_db amr_kma_index -ef || {
	    exit_status=\$?
	    if [[ \$exit_status -ne 95 ]]; then
		echo "kma failed with exit status \$exit_status"
		exit \$exit_status
	    fi
	    echo "Ignoring kma exit status 95, continuing..."
	}

	kma_to_amr_abundance.py amr_kma_OUTPUT.mapstat amr_kma_OUTPUT.res amrfinder_output/amr_abundance.csv
	amr_heatmap.py amrfinder_output/amr_abundance.csv amrfinder_output/amr_abundance_heatmap_top_genes.jpg
	"""
	}
process Deep_Arg {

	container 'bladerunner2945/deep_arg:latest'
	publishDir "${params.output_dir}/deep_arg_output", mode: 'copy'

	input:
	path merged_bins_fa

	output:
	path "deep_arg_output", emit: deep_arg_report

	script:
	"""
	mkdir -p deep_arg_output
	deeparg predict --model LS --input ${merged_bins_fa} --output deep_arg_output/deep_arg_results --d /data/deeparg_data
	"""
	}	

	

process PlasmidFinder {



	container 'staphb/plasmidfinder:latest'
	publishDir "${params.output_dir}/plasmidfinder_output", mode: 'copy'

	input:
	path merged_bins_fa


	output:
	path "plasmidfinder_output/*", emit: plasmidfinder_output_files


	script:
	"""
	mkdir -p plasmidfinder_output
	plasmidfinder.py -i ${merged_bins_fa} -o plasmidfinder_output -x -q
	"""
	}

process VirSorter2 {

	container 'staphb/virsorter2:latest'
	publishDir "${params.output_dir}/virsorter2_output", mode: 'copy'

	input:
	path merged_bins_fa
	path virsorter2_db

	output:
	path "virsorter2_output", emit: virsorter2_output_files


	script:
	"""
	mkdir -p virsorter2_output
	virsorter run -i ${merged_bins_fa} -w virsorter2_output --provirus-off --max-orf-per-seq 10 all -j 12 --db-dir ${virsorter2_db} --min-score 0.9
	"""
	}


workflow  {

    fastQC(params.first_fastq, params.second_fastq)
    trimmed_reads = Trimmomatic(params.first_fastq, params.second_fastq, params.adapters)
    assembly = MegaHit(trimmed_reads.forward_paired, trimmed_reads.reverse_paired)
    QUAST(assembly.contigs)
    Kraken2(trimmed_reads.forward_paired, trimmed_reads.reverse_paired, params.kraken2_db)
    binning_results = metawrapBinning(assembly.contigs, trimmed_reads.forward_paired, trimmed_reads.reverse_paired)
    dastool_output_bins = DAStool(binning_results.concoct_bins_tsv, binning_results.maxbin2_bins_tsv, binning_results.metabat2_bins_tsv, assembly.contigs)
    checkm_output = CheckM(dastool_output_bins)
    metaphlan_profile_reads = Metaphlan_reads(trimmed_reads.forward_paired, trimmed_reads.reverse_paired)
    merged_bins = PrependBinNames(dastool_output_bins)
    vsearch_output = Vsearch(merged_bins, params.utax_reference_db)
    amrfinder_output = AMRfinder(merged_bins,trimmed_reads.forward_paired, trimmed_reads.reverse_paired)
    deeparg_results = Deep_Arg(merged_bins)
    plasmidfinder_output = PlasmidFinder(merged_bins)
    virsorter2_output = VirSorter2(merged_bins, params.virsorter2_db)
    gtdbtk_output = GTDBTK(dastool_output_bins, params.gtdbtk_db)
    
     }

