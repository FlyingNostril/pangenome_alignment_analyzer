#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_gff3_mrnas(gff3_file):
    mrnas = defaultdict(list)
    with open(gff3_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            chrom, _, feature_type, start, end, _, strand, _, attributes = parts
            if feature_type != "mRNA":
                continue
            attr_dict = dict(kv.split("=", 1) for kv in attributes.split(";") if "=" in kv)
            mrna_id = attr_dict.get("ID")
            if not mrna_id:
                continue
            mrnas[chrom].append({
                "id": mrna_id,
                "start": int(start),
                "end": int(end),
                "strand": strand
            })
    return mrnas

def parse_bed_line(line):
    fields = line.strip().split()
    if len(fields) < 4:
        raise ValueError(f"Invalid BED line: {line.strip()}")
    return fields[0], int(fields[1]), int(fields[2]), fields[3]

def extract_mrnas_to_regions(bed_file, genome_seq, mrnas, label="reference"):
    region_seqs = defaultdict(list)

    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            try:
                chrom, bed_start, bed_end, region_name = parse_bed_line(line)
                if chrom not in mrnas:
                    continue
                for mrna in mrnas[chrom]:
                    if mrna["end"] < bed_start or mrna["start"] > bed_end:
                        continue  # no overlap
                    if chrom not in genome_seq:
                        continue
                    seq = genome_seq[chrom].seq[mrna["start"]-1:mrna["end"]]
                    if mrna["strand"] == "-":
                        seq = seq.reverse_complement()
                    seq_record = SeqRecord(seq, id=f"{label}_{mrna['id']}_{region_name}", description="")
                    region_seqs[region_name].append(seq_record)
            except Exception as e:
                print(f"Skipping line in BED: {line.strip()} - {e}")
    return region_seqs

def create_gene_id_list(region_name, region_records, output_dir):
    gene_ids=set()
    for record in region_records:
        parts = record.id.split("_")
        if len(parts) >= 2:
            gene_id = parts[1]
            gene_ids.add(gene_id)
    output_file = os.path.join(output_dir, f"{region_name}_gene_IDs.txt")
    with open(output_file, "w") as out:
        for gid in sorted(gene_ids):
            out.write(gid + "\n")
    print(f"Wrote gene ID list to {output_file}")
    

def main():
    parser = argparse.ArgumentParser(description="Extract mRNA transcripts overlapping BED regions from a reference genome.")
    parser.add_argument("-b", "--bed", help="BED file with regions (4 columns: chrom, start, end, name)")
    parser.add_argument("-g", "--genome", help="Reference genome FASTA file")
    parser.add_argument("-G", "--gff3", help="Reference GFF3 annotation file")
    parser.add_argument("-o", "--out", default="reference_transcripts", help="Output directory")
    args = parser.parse_args()

    print("Loading genome...")
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

    print("Parsing GFF3 for mRNA entries...")
    mrna_dict = parse_gff3_mrnas(args.gff3)

    print("Extracting region transcripts...")
    region_records = extract_mrnas_to_regions(args.bed, genome, mrna_dict, label="reference")
    
    print("Writing reference transcripts")
    os.makedirs(args.out, exist_ok=True)
    for region, records in region_records.items():
        out_path = os.path.join(args.out, f"{region}_reference_transcripts.fasta")
        SeqIO.write(records, out_path, "fasta")
        print(f"Wrote: {out_path}")

    print("creating gene ID list")
    for region, records in region_records.items():
        create_gene_id_list(region, records, args.out)

if __name__ == "__main__":
    main()

