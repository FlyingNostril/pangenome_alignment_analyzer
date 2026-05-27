#!/usr/bin/env python3

import argparse
import os
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def open_fasta(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def load_genome_index(genome_path):
    with open_fasta(genome_path) as handle:
        genome = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return genome

def parse_bed_line(line):
    fields = line.strip().split()
    if len(fields) < 3:
        raise ValueError(f"Invalid BED line: {line.strip()}")
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    name = fields[3] if len(fields) > 3 else f"{chrom}_{start+1}_{end}"
    return chrom, start, end, name

def extract_region_record(genome, chrom, start, end, name):
    if chrom not in genome:
        raise ValueError(f"Chromosome {chrom} not found in genome.")
    seq = genome[chrom][start:end]
    record_id = f"{chrom}:{start+1}-{end}"  ### 1-based BED-style
    record_desc = name
    return SeqRecord(seq.seq, id=record_id, description=record_desc)

def main():
    parser = argparse.ArgumentParser(description="Extract regions from a genome FASTA using a BED file into a single multi-FASTA.")
    parser.add_argument("-b", "--bed", help="BED file with regions to extract (with optional name column)")
    parser.add_argument("-g", "--genome", help="Reference genome FASTA (can be .gz)")
    parser.add_argument("-o", "--out", default="extracted_regions.fasta", help="Output FASTA file name (default: extracted_regions.fasta)")
    args = parser.parse_args()

    print("Loading genome...")
    genome = load_genome_index(args.genome)

    records = []
    print("Processing BED file...")
    with open(args.bed) as bed:
        for line in bed:
            if line.startswith("#") or line.strip() == "":
                continue
            try:
                chrom, start, end, name = parse_bed_line(line)
                record = extract_region_record(genome, chrom, start, end, name)
                records.append(record)
            except Exception as e:
                print(f"Skipping line: {line.strip()} - {e}")

    if not records:
        print("No valid regions extracted.")
        return

    with open(args.out, "w") as outfile:
        SeqIO.write(records, outfile, "fasta")
    print(f"Multi-FASTA written to {args.out}")

if __name__ == "__main__":
    main()

