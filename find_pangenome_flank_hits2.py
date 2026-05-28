#!/usr/bin/env python3

import os
import subprocess
import argparse
from Bio import SeqIO
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_blast(query_fasta, subject_fasta, db_prefix, out_file, db_loc, threads, max_targets):
    blast_db_loc = os.path.join(db_loc, db_prefix)
    required = [f"{blast_db_loc}.nhr", f"{blast_db_loc}.nin", f"{blast_db_loc}.nsq"]
    if not all(os.path.exists(p) and os.path.getsize(p) > 0 for p in required):
        subprocess.run(["makeblastdb", "-in", subject_fasta, "-dbtype", "nucl", "-out", blast_db_loc], check=True)

    blast_cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", blast_db_loc,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-max_target_seqs", str(max_targets),
        "-word_size", "25",
        "-evalue", "1e-10",
        "-num_threads", str(threads),
        "-out", out_file
    ]
    subprocess.run(blast_cmd, check=True)

def parse_best_hits(blast_output_path):
    best_starts = {}
    best_ends = {}

    with open(blast_output_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 12:
                continue
            # set up the BLAST hit variables:
            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            mismatch = int(fields[4])
            gapopen = int(fields[5])
            qstart = int(fields[6])
            qend = int(fields[7])
            sstart = int(fields[8])
            send = int(fields[9])
            evalue = float(fields[10])
            score = float(fields[11])
            
            key = qseqid.replace("_start", "").replace("_end", "")
            
            # collect all of the start and end for each chromosome


            if qseqid.endswith("_start"):
                # prefer a higher bitscore, break ties with a lower sstart
                if key not in best_starts or score > float(best_starts[key][11]) or (score == float(best_starts[key][11] and sstart < int(best_starts[key][8]))):
                    best_starts[key] = fields
            elif qseqid.endswith("_end"):
                if key not in best_ends or score > float(best_ends[key][11] or (score == float(best_ends[key][11]) and send > int(best_ends[key][9]))):
                    best_ends[key] = fields

    return best_starts, best_ends

def find_pangenome_regions(genome_file, outdir, flank_fasta, threads, blast_db, max_targets):
    genome_base = os.path.splitext(os.path.basename(genome_file))[0]
    db_prefix = genome_base + "_blastdb"
    blast_out = os.path.join(outdir, genome_base + "_flanks_blast.tsv")

    run_blast(flank_fasta, genome_file, db_prefix, blast_out, blast_db, threads, max_targets)

    best_starts, best_ends = parse_best_hits(blast_out)
    print("Start keys:", best_starts.keys())
    print("End Keys:", best_ends.keys())
    bed_file = os.path.join(outdir, genome_base + "_best_hits.bed")
    with open(bed_file, "w") as out:
        for name in sorted(set(best_starts) & set(best_ends)):
            start_hit = best_starts[name]
            end_hit = best_ends[name]
            chr_start = start_hit[1]
            chr_end   = end_hit[1]
            if chr_start != chr_end:
                print(f"Warning: {name} start and end on different chromosomes. Skipping.")
                continue

            start_pos = min(int(start_hit[8]), int(start_hit[9])) - 1  # 0-based
            end_pos   = max(int(end_hit[8]), int(end_hit[9]))          # end-exclusive
            out.write(f"{chr_start}\t{start_pos}\t{end_pos}\t{name}\n")

    print(f"Wrote best hits to {bed_file}")

# functions meant to use with the gene ID mode
def load_gene_ids(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip() and not line.startswith("#")}

def find_pangenome_genes(genome_file, outdir, flank_fasta, threads, blast_db, gene_id_file, max_targets):
    genome_base = os.path.splitext(os.path.basename(genome_file))[0]
    db_prefix = genome_base + "_blastdb"
    blast_out = os.path.join(outdir, genome_base + "flanks_blast.tsv")

    requested_ids = load_gene_ids(gene_id_file)

    run_blast(flank_fasta, genome_file, db_prefix, blast_out, blast_db, threads, max_targets)

    best_starts, best_ends = parse_best_hits(blast_out)
    bed_file = os.path.join(outdir, genome_base + "_best_hits.bed")

    with open(bed_file, "w") as out:
        for name in sorted(requested_ids):
            if name not in best_starts or name not in best_ends: 
                print(f"Warning, missing start or end hit for gene {name} in {genome_base}")
                continue
            start_hit = best_starts[name]
            end_hit = best_ends[name]

            chr_start = start_hit[1]
            chr_end = end_hit[1]
            if chr_start != chr_end:
                print(f"Warning: blast hit for gene {name} starts and ends on different chromosomes in {genome_base}, skipping.")
                continue
            start_pos = min(int(start_hit[8]), int(start_hit[9])) - 1
            end_pos = max(int(end_hit[8]), int(end_hit[9]))

            out.write(f"{chr_start}\t{start_pos}\t{end_pos}\t{name}\n")

    print(f"Wrote gene-specific best hits to {bed_file}")


def main():
    parser = argparse.ArgumentParser(description="BLAST region flanks against multiple genome FASTA files and report best hit for each.")
    parser.add_argument("-i", "--flank_fasta", help="Input FASTA with _start and _end sequences")
    parser.add_argument("-l", "--genome_list", help="Text file listing genome FASTA files (one per line)")
    parser.add_argument("-g", "--genome_dir", help="directory where the pangenome fasta files are located")
    parser.add_argument("-o", "--outdir", default="blast_results", help="Output directory (default: blast_results/)")
    parser.add_argument("-b", "--blast_db", required=True, help="FULL PATH to blast database locations")
    parser.add_argument("-t", "--threads", default = 1, type = int, help="Number of threads allowed for BLAST to find pangenome flank hits")
    parser.add_argument("-T", "--max_threads", default =1, type = int, help="Maximum number of threads alloted to the program")
    parser.add_argument("-G", "--gene_id_file", help="The gene id file used with gene ID mode. Using this restricts output to requested genes only")
    parser.add_argument("-m", "--max_targets", type=int, default=5, help="The number of matches blast will search for. Default 5. \nIf this step runs slowly then reduce this number")
    args = parser.parse_args()

    max_concurrent = max(1, args.max_threads // args.threads)
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.blast_db, exist_ok=True)

    with open(args.genome_list) as f:
        genome_files = [os.path.join(args.genome_dir, line.strip()) for line in f if line.strip()]

    worker_fn = find_pangenome_genes if args.gene_id_file else find_pangenome_regions

    with ThreadPoolExecutor(max_workers=max_concurrent) as ex:
        if args.gene_id_file:
            futures = [ex.submit(worker_fn, genome_file, args.outdir, args.flank_fasta, args.threads, args.blast_db, args.gene_id_file, args.max_targets)
                for genome_file in genome_files]
        else:
            futures = [ex.submit(worker_fn, genome_file, args.outdir, args.flank_fasta, args.threads, args.blast_db, args.max_targets)
                    for genome_file in genome_files]

        for f in as_completed(futures):
            f.result()
        

if __name__ == "__main__":
    main()

