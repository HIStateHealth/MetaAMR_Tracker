#!/usr/bin/env python3
import pandas as pd
import re
import argparse

def partial_match(query, target_list):
    """Check for partial match of query in a list of targets using regular expressions."""
    pattern = re.compile(re.escape(query), re.IGNORECASE)
    return any(pattern.search(target) for target in target_list)

def find_classification_for_bin(bin_query, gtdbtk_data):
    """Find classification based on bin query in GTDBtk data."""
    return gtdbtk_data.get(bin_query, "Unknown")

def update_amr_gene_location(row, primary_viral_contigs, primary_plasmid_contigs):
    """Determine AMR gene location based on viral and plasmid data."""
    contig_part = '_'.join(row['Contig id'].split('_')[:4])
    if partial_match(contig_part, primary_viral_contigs) or partial_match(contig_part, primary_plasmid_contigs):
        return "Bac chromosome"
    elif row['Location_identity'] == "Unknown":
        return "Unknown"
    else:
        return "Bac chromosome"

def main(args):
    """Main function to load data, process, and save to an Excel file."""
    # Load data
    amr_finder_df = pd.read_csv(args.amr_finder, delimiter='\t')
    gtdbtk_df = pd.read_csv(args.gtdbtk, delimiter='\t')
    plasmid_df = pd.read_csv(args.plasmid, delimiter='\t')
    viral_df = pd.read_csv(args.viral, delimiter='\t')

    # Clean and extract primary contig parts for plasmid and viral data
    plasmid_df['primary_contig'] = plasmid_df['Contig'].apply(lambda x: x.split(' ')[0])
    primary_plasmid_contigs = plasmid_df['primary_contig'].unique()
    viral_df['primary_viral_contig'] = viral_df['seqname'].apply(lambda x: x.split('||')[0])
    primary_viral_contigs = viral_df['primary_viral_contig'].unique()

    # Extract bin numbers correctly by splitting at '_k'
    amr_finder_df['Bin Number'] = amr_finder_df['Contig id'].apply(lambda x: x.split('_k')[0])

    # Create a map from the GTDBtk data using 'user_genome' as keys for classifications
    gtdbtk_map = gtdbtk_df.set_index('user_genome')['classification'].to_dict()

    # Apply mapping functions
    amr_finder_df['Location_identity'] = amr_finder_df['Bin Number'].apply(lambda x: find_classification_for_bin(x, gtdbtk_map))
    amr_finder_df['AMR_gene_location (bac_chromosome, plasmid, phage)'] = amr_finder_df.apply(
        lambda row: update_amr_gene_location(row, primary_viral_contigs, primary_plasmid_contigs), axis=1
    )

    # Rearrange columns to put new columns at the beginning
    cols = ['AMR_gene_location (bac_chromosome, plasmid, phage)', 'Location_identity'] + [col for col in amr_finder_df.columns if col not in ['AMR_gene_location (bac_chromosome, plasmid, phage)', 'Location_identity']]
    amr_finder_df = amr_finder_df[cols]

    # Save to Excel
    amr_finder_df.to_excel(args.output, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process AMR Finder data with associated microbial classification and resistance gene location.")
    parser.add_argument('--amr_finder', type=str, required=True, help='Path to the AMR Finder report with abundance file')
    parser.add_argument('--gtdbtk', type=str, required=True, help='Path to the GTDBtk report file')
    parser.add_argument('--plasmid', type=str, required=True, help='Path to the plasmid report file')
    parser.add_argument('--viral', type=str, required=True, help='Path to the viral report file')
    parser.add_argument('--output', type=str, required=True, help='Output path for the Excel file')

    args = parser.parse_args()
    main(args)
