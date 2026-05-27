#!/usr/bin/env python3

import os
import csv
import argparse
import re
from collections import defaultdict, OrderedDict

def sniff_delimiter(filepath):
    with open(filepath, newline="") as csvfile:
        sample = csvfile.read(2048)
        sniffer = csv.Sniffer()
        return sniffer.sniff(sample).delimiter

def csv_rows_to_intervals(input_csv, flank):
    """
    Parse CSV/TSV -> chrom -> list of (start, end, name), preserving first-seen order.
    Deduplicates identical region_name rows.
    """
    intervals = defaultdict(list)
    seen_regions = set()
    seen_order = OrderedDict()  # name -> first index appearance
    order_counter = 0

    delimiter = sniff_delimiter(input_csv)
    with open(input_csv, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=delimiter)
        for row in reader:
            try:
                region_name = row[reader.fieldnames[0]].strip()
                if region_name in seen_regions:
                    # skip duplicates by name (exact same label)
                    continue
                seen_regions.add(region_name)
                if region_name not in seen_order:
                    seen_order[region_name] = order_counter
                    order_counter += 1

                chrom_number = int(row[reader.fieldnames[1]])
                position = int(row[reader.fieldnames[2]])
                chrom = f"CHR_{chrom_number:02d}"

                start = max(0, position - flank)  # 0-based inclusive
                end = position + flank            # 0-based exclusive
                intervals[chrom].append((start, end, region_name))
            except Exception as e:
                print(f"Skipping row {row} due to error: {e}")

    return intervals, seen_order

def merge_intervals(intervals_by_chrom, name_order):
    """
    Merge overlapping/touching intervals per chromosome.
    For each merged block, join all original names (in first-seen input order) with underscores.
    Returns: chrom -> list of (start, end, merged_name)
    """
    merged = {}
    for chrom, items in intervals_by_chrom.items():
        if not items:
            merged[chrom] = []
            continue

        # Sort by start, then end
        items_sorted = sorted(items, key=lambda x: (x[0], x[1]))

        cur_start, cur_end = items_sorted[0][0], items_sorted[0][1]
        cur_names = {items_sorted[0][2]}

        merged_list = []
        for s, e, nm in items_sorted[1:]:
            if s <= cur_end:  # treat touching as overlap
                # extend
                if e > cur_end:
                    cur_end = e
                cur_names.add(nm)
            else:
                # flush current
                merged_name = "_".join(sorted(cur_names, key=lambda n: name_order.get(n, 10**12)))
                merged_list.append((cur_start, cur_end, merged_name))
                # start new
                cur_start, cur_end, cur_names = s, e, {nm}

        # flush last
        merged_name = "_".join(sorted(cur_names, key=lambda n: name_order.get(n, 10**12)))
        merged_list.append((cur_start, cur_end, merged_name))

        merged[chrom] = merged_list

    return merged

def load_gene_ids(path):
    """
    take a list of gene IDs and make sure they fit the criteria for our gff3 file
    remove duplicate IDs
    """
    valid_ids = []
    duplicate_ids = []
    bad_format_ids = []
    seen = set()
    pattern = re.compile(r"-R[A-Za-z]$")
    with open(path) as g:
        for line in g:
            if line.startswith("#") or not line.strip():
                continue
            gene_id = line.strip()
            if not pattern.search(gene_id):
                bad_format_ids.append(gene_id)
            elif gene_id in seen:
                duplicate_ids.append(gene_id)
            else:
                seen.add(gene_id)
                valid_ids.append(gene_id)

    return valid_ids, duplicate_ids, bad_format_ids

def parse_gff3_mrnas_by_id(gff3_file, gene_ids, bed_path, flank):
    """
    Search the input gff3 file for mRNA entries matching the input mRNA ids
    return the information needed to generate a bedfile, one line per mRNA
    """
    mrnas = defaultdict(list)
    genelist = set(gene_ids)
    with open(gff3_file) as gf:
        for line in gf:
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
            
            if mrna_id not in genelist:
                continue
            
            start_init = int(start) - 1
            end_init = int(end)

            if strand == '-':
                adj_start = max(0, start_init - 1000)
                adj_end = end_init + flank
            else:
                adj_start = max(0, start_init - flank)
                adj_end = end_init + 1000

            mrnas[chrom].append({
                "start": adj_start,
                "end": adj_end,
                "name": mrna_id
                })
    with open(bed_path, "w") as bedfile:
        for chrom in sorted(mrnas):
            rows = mrnas[chrom]
            for row in sorted(rows, key=lambda x: (x["start"], x["end"])):
                start = row["start"]
                end = row["end"]
                out_name = row["name"]
                bedfile.write(f"{chrom}\t{start}\t{end}\t{out_name}\n")
    


def write_bed(intervals_by_chrom, bed_path, max_name_len=None):
    """
    intervals_by_chrom: chrom -> list of (start, end, name)
    If max_name_len is set and exceeded, truncate and warn once per write.
    """
    warned = False
    with open(bed_path, "w") as bedfile:
        for chrom in sorted(intervals_by_chrom):
            rows = intervals_by_chrom[chrom]
            for start, end, name in sorted(rows, key=lambda x: (x[0], x[1])):
                out_name = name
                if max_name_len and len(out_name) > max_name_len:
                    if not warned:
                        print(f"Warning: one or more merged names exceed --max-name-len ({max_name_len}); truncating.")
                        warned = True
                    out_name = out_name[:max_name_len]
                bedfile.write(f"{chrom}\t{start}\t{end}\t{out_name}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert a CSV/TSV of markers to a BED file; optionally merge overlapping intervals and annotate merged names in the BED name field.\nAlternatively, use a list of genes to generate a bed file for downstream use")
    # the bedfile and gene ID options are mutually exclusive
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-c", "--csvfile", help="Input CSV/TSV file (first 3 cols: name, chromosome_number, position)")
    group.add_argument("-g", "--gene_id_list", help="Input gene ID list, one ID per line. Requires -G (--gff3)")
    
    parser.add_argument("-b", "--bedfile", help="BED output filename (default: <input basename>.bed)")
    parser.add_argument("-G", "--gff3", help="gff3 file used with -g (--gene_id_list) to generate the genelist bedfile")
    parser.add_argument("-f", "--flank", type=int, default=50000, help="Flanking size for each side of the position (default: 50000)")
    parser.add_argument("-F", "--gene_id_flank", type=int, default=2000, help="Upstream flanking BP for gene id based bedfile. Default 2000.")
    parser.add_argument("--merge", action="store_true", help="Merge overlapping/touching intervals and join names with underscores")
    parser.add_argument("--max-name-len", type=int, default=None, help="Optional max length for BED name field; long names get truncated")
    args = parser.parse_args()

    if args.gene_id_list and not args.gff3:
      parser.error("-G (--gff3) is required when using -g (--gene_id_list)")

    if args.gene_id_list:
        if not args.bedfile:
            output_bed = os.path.splitext(args.gene_id_list)[0] + ".bed"
        else:
            output_bed = os.path.join(args.bedfile)
        print("parsing the gene ID list and writing a bedfile for it.")
        gene_ids, double_ids, bad_ids = load_gene_ids(args.gene_id_list)
        print("gene ID list parsed:")
        print("these gene ids failed to parse:")
        print(f"Duplicated: {double_ids}")
        print(f"bad format: {bad_ids}")
        parse_gff3_mrnas_by_id(args.gff3, gene_ids, output_bed, args.gene_id_flank)
        print(f"wrote gene-id based bedfile to {output_bed}")
    else:
        if not args.bedfile:
            output_bed = os.path.splitext(args.csvfile)[0] + ".bed"
        else:
            output_bed = os.path.join(args.bedfile)
        intervals_by_chrom, name_order = csv_rows_to_intervals(args.csvfile, args.flank)
        if not args.merge:
            # Unmerged: just write (start, end, original name)
            write_bed(intervals_by_chrom, output_bed, max_name_len=args.max_name_len)
            print(f"Wrote BED (unmerged): {output_bed}")
        else:
            merged_by_chrom = merge_intervals(intervals_by_chrom, name_order)
            write_bed(merged_by_chrom, output_bed, max_name_len=args.max_name_len)
            print(f"Wrote BED (merged, names joined by '_'): {output_bed}")

if __name__ == "__main__":
    main()

