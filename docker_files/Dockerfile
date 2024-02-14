# Use bioinforpi/semibin2:v0.0.1 as the base image
FROM bioinforpi/semibin2:v0.0.1

# Update package lists and install minimap2 and samtools
RUN apt-get update && \
    apt-get install -y minimap2 samtools

# (Optional) Clean up to reduce image size
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*
