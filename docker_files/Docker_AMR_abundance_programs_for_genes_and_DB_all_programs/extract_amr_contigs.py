#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

# Function to reverse complement a sequence if needed
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def main(args):
    # Load the table data, assuming tab-separated format
    df = pd.read_csv(args.table_path, sep='\t')

    # Load the contig FASTA file
    contigs = SeqIO.to_dict(SeqIO.parse(args.contigs_path, "fasta"))

    # Process each row and extract sequences
    extracted_sequences = []
    for index, row in df.iterrows():
        contig_id = row['Contig id']
        start = row['Start'] - 1  # Convert to 0-based index
        stop = row['Stop']
        strand = row['Strand']
        accession = row['Accession of closest sequence']
        
        if contig_id in contigs:
            sequence = contigs[contig_id].seq[start:stop]
            if strand == '-':
                sequence = reverse_complement(sequence)
            
            # Create a new SeqRecord
            new_record = SeqRecord(
                Seq(sequence),
                id=accession,
                description=""
            )
            extracted_sequences.append(new_record)

    # Write the extracted sequences to a new FASTA file
    SeqIO.write(extracted_sequences, args.output_path, "fasta")

    print(f"Extracted sequences have been saved to {args.output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract specific sequences from a FASTA file based on a TSV file.")
    parser.add_argument("--table_path", required=True, help="Path to the tab-separated values (TSV) table file.")
    parser.add_argument("--contigs_path", required=True, help="Path to the contig FASTA file.")
    parser.add_argument("--output_path", required=True, help="Path to output the extracted sequences in FASTA format.")
    args = parser.parse_args()
    main(args)

