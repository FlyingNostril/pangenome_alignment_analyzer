#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ThreadPoolExecutor, as_completed

def open_fasta(path):
    return open(path, "r")

def load_genome_index(genome_path):
    with open_fasta(genome_path) as handle:
        genome = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return genome

def parse_bed_line(line):
    fields = line.strip().split()
    if len(fields) < 4:
        raise ValueError(f"Invalid BED line: {line.strip()}")
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    name = fields[3]
    return chrom, start, end, name

def process_bed_file(bed_file, bed_dir, genome_dir):
    region_sequences = defaultdict(list)
    genome_name = bed_file.replace("_best_hits.bed", "")

    genome_fasta_path = os.path.join(genome_dir, f"{genome_name}.fasta")
    if not os.path.exists(genome_fasta_path):
        genome_fasta_path = os.path.join(genome_dir, f"{genome_name}.fa")
    if not os.path.exists(genome_fasta_path):
        return region_sequences

    genome = load_genome_index(genome_fasta_path)

    with open(os.path.join(bed_dir, bed_file)) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, region_name = parse_bed_line(line)
            if chrom not in genome:
                continue
            seq = genome[chrom].seq[start:end]
            region_sequences[region_name].append(
                SeqRecord(seq, id=f"{genome_name}_{region_name}", description = "")
                )
    print(f"processing {bed_file} is complete")
    return region_sequences

def merge_region_maps(region_maps):
    merged = defaultdict(list)
    for m in region_maps:
        for region_name, records in m.items():
            merged[region_name].extend(records)
    return merged

def write_region_fasta(region_sequences, output_dir):
    # Write one FASTA per region
    for region_name, records in region_sequences.items():
        out_path = os.path.join(output_dir, f"{region_name}.fasta")
        SeqIO.write(records, out_path, "fasta")
        print(f"Wrote: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate one multi-FASTA file per region from BED and genome files.")
    parser.add_argument("-b", "--bed_dir", help="Directory containing *_best_hits.bed files")
    parser.add_argument("-g", "--genome_dir", help="Directory containing genome FASTA files (named genome1.fasta, etc.)")
    parser.add_argument("-o", "--output_dir", help="Directory to write per-region FASTA files")
    parser.add_argument("-t", "--max_threads", type=int, default=1, help="Number of threads allotted to the pipline. For this part of the pipeline, concurrent jobs = max_threads")
    args = parser.parse_args()
   
    os.makedirs(args.output_dir, exist_ok=True)

    bed_files = [f for f in os.listdir(args.bed_dir) if f.endswith("_best_hits.bed")]

    with ThreadPoolExecutor(max_workers=args.max_threads) as ex:
        futures = [ex.submit(process_bed_file, bed_file, args.bed_dir, args.genome_dir)
                for bed_file in bed_files]
        region_maps = [f.result() for f in as_completed(futures)]
        
    all_region_sequences = merge_region_maps(region_maps)
    write_region_fasta(all_region_sequences, args.output_dir)



if __name__ == "__main__":
    main()

