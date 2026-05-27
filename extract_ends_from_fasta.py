#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def extract_ends(input_fasta, output_fasta, flank_size=1000):
    output_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        desc = record.description.strip()
        try:
            name = desc.split()[1]  # ID<TAB>name
        except IndexError:
            name = desc

        sequence = record.seq
        seq_len = len(sequence)

        if flank_size >= seq_len:
            print(f"Skipping {name}: flank size {flank_size} >= sequence length {seq_len}")
            continue

        # First flank_size bp or less
        start_seq = sequence[:flank_size]
        start_record = SeqRecord(start_seq, id=f"{name}_start", description="")
        output_records.append(start_record)

        # Last flank_size bp or less
        end_seq = sequence[-flank_size:]
        end_record = SeqRecord(end_seq, id=f"{name}_end", description="")
        output_records.append(end_record)

    SeqIO.write(output_records, output_fasta, "fasta")
    print(f"Wrote: {output_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Extract start and end flanks from each FASTA record.")
    parser.add_argument("-i", "--input_fasta", help="Multi-FASTA input file (e.g., example_regions.fasta)")
    parser.add_argument("-o", "--output_fasta", help="Output FASTA file for flanking regions")
    parser.add_argument("-f", "--flank", type=int, default=1000, help="Flank size in bp (default: 1000)")
    args = parser.parse_args()

    extract_ends(args.input_fasta, args.output_fasta, args.flank)

if __name__ == "__main__":
    main()

