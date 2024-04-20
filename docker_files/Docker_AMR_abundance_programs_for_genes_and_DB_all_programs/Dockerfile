# Start with Ubuntu as the base image
FROM ubuntu:latest

# Install dependencies for AMRFinderPlus, KMA, and Python scripts
# Added HMMER for hmmpress
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    git \
    perl \
    ncbi-blast+ \
    hmmer \
    python3 \
    python3-pip \
    unzip \
    libcurl4-gnutls-dev

# Install AMRFinderPlus from GitHub
RUN git clone https://github.com/ncbi/amr.git AMRFinderPlus && \
    cd AMRFinderPlus && \
    make && \
    make install

# Download the AMRFinderPlus database
RUN amrfinder_update

# Install Python packages needed for the scripts
RUN pip3 install pandas seaborn matplotlib numpy biopython

# Assume you've downloaded and prepared KMA source code beforehand
# Adjust this COPY command according to your actual KMA source code location
COPY genomicepidemiology-kma-c10a387d8f43 /kma_source
WORKDIR /kma_source
RUN make && \
    cp kma /usr/local/bin/ && \
    cp kma_index /usr/local/bin/

# Add the Python scripts and AMR_CDS.fasta to the container
COPY amr_heatmap_from_db_genes.py /scripts/
COPY heatmap_from_amr_finder_report_abundance.py /scripts/
COPY update_amr_finder_report_with_abundance.py /scripts/
COPY extract_amr_contigs.py /scripts/
COPY kma_to_amr_abundance_from_db_genes.py /scripts/
COPY kma_to_amr_abundance_for_amrfindergenes.py /scripts/
COPY Modified_AMR_CDS_3.12.fasta /data/

# Make scripts in /scripts executable from anywhere
RUN chmod +x /scripts/*.py
ENV PATH="/scripts:${PATH}"

# Adjust the working directory as needed
WORKDIR /

# Set the default command or entrypoint
CMD ["/bin/bash"]
