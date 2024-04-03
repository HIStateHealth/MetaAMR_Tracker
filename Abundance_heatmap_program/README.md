# AMR Abundance Heatmap Generator

## Overview
The AMR Abundance Heatmap Generator is a Python tool designed to visualize the top 20 genes by antimicrobial resistance (AMR) abundance using a heatmap. The tool processes input data from a CSV file, sorts genes by their AMR abundance, and generates a heatmap with a custom green-to-red color scale to represent the data visually. This tool is useful for researchers and scientists who want to quickly identify and visualize key genes associated with AMR in microbial genomes.

## Installation

### Prerequisites
- Python 3.x
- pandas
- numpy
- seaborn
- matplotlib

Ensure Python is installed on your system. You can download Python from [python.org](https://www.python.org/downloads/). After installing Python, install the required libraries using pip:

### Example
This command processes `gene_data.csv`, generates a heatmap for the top 20 genes by AMR abundance, and saves the heatmap as `gene_heatmap.jpeg`.
