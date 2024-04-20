#!/usr/bin/env python3

"""
Update AMR Finder Report with AMR Abundance

This script appends AMR abundance data to the AMR Finder Report based on the accession numbers.
It takes the path to an AMR Finder Report (a TSV file), the path to an AMR Abundance data file (CSV),
and the path for the output file as command-line arguments.

Usage:
Run the script from the command line with the following syntax:
python update_amr_finder_report_with_abundance.py --amr_report_path <path_to_amr_report.txt> --amr_abundance_path <path_to_amr_abundance.csv> --output_path <output_file_path.txt>
"""

import pandas as pd
import argparse

def main(args):
    # Load the AMR Finder Report
    amr_report = pd.read_csv(args.amr_report_path, sep='\t', header=0)

    # Load the AMR Abundance file
    amr_abundance = pd.read_csv(args.amr_abundance_path)

    # Merge the data based on 'Accession of closest sequence' and '#Template' columns
    merged_data = pd.merge(amr_report, amr_abundance, how='left', left_on='Accession of closest sequence', right_on='#Template')

    # Optionally drop the duplicate accession column from amr_abundance
    merged_data.drop(columns=['#Template'], inplace=True)

    # Save the merged data back to a new file
    merged_data.to_csv(args.output_path, sep='\t', index=False)

    print(f"Updated AMR Finder Report has been saved to {args.output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Append AMR abundance data to AMR Finder Report based on accession numbers.")
    parser.add_argument("--amr_report_path", required=True, help="Path to the AMR Finder Report file (TSV format).")
    parser.add_argument("--amr_abundance_path", required=True, help="Path to the AMR Abundance file (CSV format).")
    parser.add_argument("--output_path", required=True, help="Path for the output file to save the updated AMR Finder Report (TSV format).")
    args = parser.parse_args()
    main(args)
