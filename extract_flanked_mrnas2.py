#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess

def ensure_blast_db(fasta_path, db_prefix, force_ignore):
    db_dir = os.path.dirname(db_prefix)
    if db_dir:
        os.makedirs(db_dir, exist_ok=True)

    required = [f"{db_prefix}.nhr", f"{db_prefix}.nin", f"{db_prefix}.nsq"]
    
    db_exists = all(os.path.exists(p) and os.path.getsize(p) > 0 for p in required)

    if db_exists and not force_ignore:
        fasta_mtime = os.path.getmtime(fasta_path)
        db_mtime = min(os.path.getmtime(p) for p in required)
        if db_mtime >= fasta_mtime:
            print(f"using existing BLAST database: {db_prefix}")
            return db_prefix

    # lightweight FASTA sanity check (no full load)
    try:
        has_seq = any(len(rec.seq) > 0 for rec in SeqIO.parse(fasta_path, "fasta"))
        if not has_seq:
            print(f"Warning: No sequences found in {fasta_path}. Skipping BLAST DB creation.")
            return None
    except Exception as e:
        print(f"Warning: Failed to parse {fasta_path}: {e}. Skipping BLAST DB creation.")
        return None

    print(f"creating the BLAST database for {fasta_path}")
    subprocess.run(["makeblastdb", "-in", fasta_path, "-dbtype", "nucl", "-out", db_prefix], check=True)
    return db_prefix

def run_blast(query_fasta, db_prefix, out_file, threads):
    blast_cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", db_prefix,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-evalue", "1e-10",
        "-num_threads", str(threads),
        "-out", out_file
    ]
    subprocess.run(blast_cmd, check=True)

# gotta look at this one too
def parse_blast_results(blast_file, genome_prefixes, query_lengths, min_cov, max_extra_len=100):
    best_hits = defaultdict(dict)

    with open(blast_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 12:
                continue
            qseqid, sseqid = fields[0], fields[1]
            aligned_len = int(fields[3])
            sstart, send = int(fields[8]), int(fields[9])
            bitscore = float(fields[11])
            query_len = query_lengths.get(qseqid)
            if not query_len:
                continue
            coverage = aligned_len / query_len
            if coverage < min_cov or aligned_len > (query_len + max_extra_len):
                continue
            # determine genome by matching prefix
            genome_hit=None
            for prefix in genome_prefixes:
                if sseqid.startswith(prefix + "_"):
                    genome_hit = prefix
                    break
            if not genome_hit:
                continue
            prev_hit = best_hits[qseqid].get(genome_hit)
            if not prev_hit or bitscore > prev_hit["bitscore"] or (
                bitscore == prev_hit["bitscore"] and aligned_len > prev_hit["aligned_len"]):
                best_hits[qseqid][genome_hit] = {
                    "sseqid": sseqid,
                    "sstart": sstart,
                    "send": send,
                    "bitscore": bitscore,
                    "aligned_len": aligned_len
                }

    return best_hits

def load_sequences_by_id(pangenome_path):
    return SeqIO.to_dict(SeqIO.parse(pangenome_path, "fasta"))

def group_reference_sequence(reference_fasta):
    ref_records = {}
    query_lengths = {}
    strands = {}
    for record in SeqIO.parse(reference_fasta, "fasta"):
        qid = record.id
        fields = qid.split("_")
        strand = fields[-1] if fields[-1] in {"+", "-"} else None
        ref_records[qid] = record
        query_lengths[qid] = len(record.seq)
        strands[qid] = strand
    return ref_records, query_lengths, strands

def parse_gff3_mrnas(gff3_file):
    mrnas = defaultdict(list)
    with open(gff3_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                continue
            chrom, _, _, start, end, _, strand, _, attributes = parts
            attr_dict = dict(kv.split("=", 1) for kv in attributes.split(";") if "=" in kv)
            mrna_id = attr_dict.get("ID")
            if mrna_id:
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
        raise ValueError(f"Invalid Bed Line: {line.strip()}")
    return fields[0], int(fields[1]), int(fields[2]), fields[3]

def extract_flanked_sequences(bed_path, mrnas, genome, flank_up, flank_down):
    flanked_records = defaultdict(list)
    with open(bed_path) as nostril:  # nostril is a funny word
        for line in nostril:
            if line.startswith('#') or not line.strip():
                continue
            try: 
                chrom, bed_start, bed_end, region_name = parse_bed_line(line)
                if chrom not in mrnas or chrom not in genome:
                    continue
                seq_len = len(genome[chrom])
                for mrna in mrnas[chrom]:
                    if mrna["end"] < bed_start or mrna["start"] > bed_end:
                        continue
                    if mrna["strand"] == "+":
                        start = max(0, mrna["start"] - flank_up - 1)
                        end = min(seq_len, mrna["end"] + flank_down)
                    else:
                        start = max(0, mrna["start"] - flank_down - 1)
                        end = min(seq_len, mrna["end"] + flank_up)

                    seq = genome[chrom].seq[start:end]
                    if mrna["strand"] == "-":
                        seq = seq.reverse_complement()

                    record_id = f"{mrna['id']}_{chrom}_{start+1}_{end}_{mrna['strand']}"
                    flanked_records[region_name].append(SeqRecord(seq, id=record_id, description=""))

            except Exception as mouse:
                print(f"Error Error, +++insert more cheese+++ Skipping line in {bed_path}: {line.strip()} - {mouse}")
    return flanked_records

def assemble_fasta_by_gene(region_outdir, ref_records, transcripts, genome_names, pangenome_seqs, reference_strands, best_hits):
        with open(os.path.join(region_outdir, "no_matches.txt"), "w") as no_match_log:
            for qid, r_record in ref_records.items():
                out_fasta = os.path.join(region_outdir, f"{qid}_best_hits.fasta")

                # initialize records for this gene:
                records = [SeqRecord(r_record.seq, id=qid, description="")]

                # add in the reference transcript:
                transcript_id = qid.split("_")[0]
                matched_transcript = next((transcripts[tid] for tid in transcripts if tid.startswith(transcript_id)), None)
                if matched_transcript:
                    records.append(matched_transcript)

                # add in the matching transcripts from the pangenome
                missing_genomes = []
                hits = best_hits.get(qid, {})
                for genome in genome_names:
                    hit = hits.get(genome)
                    if not hit:
                        missing_genomes.append(genome)
                        continue
                    seq_id = hit["sseqid"]
                    if seq_id not in pangenome_seqs:
                        missing_genomes.append(genome)
                        continue
                    full_seq = pangenome_seqs[seq_id].seq
                    hit_start = min(hit["sstart"], hit["send"]) -1 # -1 to account for 0-based coords
                    hit_end = max(hit["sstart"], hit["send"])
                    hit_seq = full_seq[hit_start:hit_end]
                    if reference_strands[qid] == "-":
                        hit_seq = hit_seq.reverse_complement()   
                    records.append(SeqRecord(hit_seq, id=f"{seq_id}_{hit_start+1}_{hit_end}", description="BLAST HIT"))

                if len(records) > 1:
                    SeqIO.write(records, out_fasta, "fasta")
                    print(f"Wrote: {out_fasta}")

                if missing_genomes:
                    no_match_log.write(f"{qid},{','.join(missing_genomes)}\n")


def process_region(outdir, region, pangenome_dir, blast_db, threads, min_cov, transcripts, genome_names, ignore):
    region_name = region.replace("_flanked_mRNAs.fasta", "")
    ref_path = os.path.join(outdir, region)
    region_outdir = os.path.join(outdir, region_name)
    os.makedirs(region_outdir, exist_ok=True)
    # find matching pangenome region
    pangenome_path= None
    for ext in [".fasta", ".fa"]:
        candidate = os.path.join(pangenome_dir, f"{region_name}{ext}")
        if os.path.exists(candidate):
            pangenome_path = candidate
            break
    if pangenome_path is None:
        print(f"Error: pangenome region for {region_name} not found.")
        return

    # use ensure_blast_db:
    db_prefix = ensure_blast_db(pangenome_path, os.path.join(blast_db, region_name), ignore)
    if db_prefix is None:
        print(f"Skipping region {region_name} due to empty or invalid FASTA.")
        return

    # run the blast:
    ref_records, q_lengths, ref_strands = group_reference_sequence(ref_path)
    temp_blast_output = os.path.join(region_outdir, "blast_hits.tsv")
    run_blast(ref_path, db_prefix, temp_blast_output, threads)

    # Parse the output:
    best_hits = parse_blast_results(temp_blast_output, genome_names, q_lengths, min_cov)
    pangenome_seqs = load_sequences_by_id(pangenome_path)

    # use the new assemble_fasta_by_gene function
    assemble_fasta_by_gene(region_outdir, ref_records, transcripts, genome_names, pangenome_seqs, ref_strands, best_hits)
    


def main():
    parser = argparse.ArgumentParser(description="Extract mRNA sequences with flanking sequence.")
    parser.add_argument("-g", "--genome", required=True, help="Reference genome FASTA.  Required")
    parser.add_argument("-l", "--pangenome_list", required=True, help="list of genomes inside of the pangenome.  Required")
    parser.add_argument("-G", "--gff3", required=True, help="GFF3 annotation file of the reference genome.  Required")
    parser.add_argument("-o", "--outdir", default="flanked_transcripts_out", help="Output FASTA file directory of the flanked_fastas from the pangenome")
    parser.add_argument("-u", "--flank_up", type=int, default=2000, help="number of bases to select upstream of the transcription start site, IE the promoter, default is 2000bp")
    parser.add_argument("-d", "--flank_down", type=int, default=1000, help="number of bases to select downstream of the transcription end site, default is 1000bp")
    parser.add_argument("-B", "--bedfile_dir", required=True, help="Directory with bed files of the reference genome regions to extract transcripts with flanking sequnce from.  Required")
    parser.add_argument("-p", "--pangenome_dir", required=True, help="directory with the pangenome region fasta files.")
    parser.add_argument("-D", "--blast_db", default='blast_db', help="blast db directory that holds the pangenome blast databases")
    parser.add_argument("-m", "--min_cov", type=float, default=0.8, help="float, 0-1, Minimum coverage required for a blast hit to be considered")
    parser.add_argument("-t", "--transcripts", required=True, help="The fasta file of reference transcripts.  Make sure that the IDs match! Required")
    parser.add_argument("--threads", type=int, default = 1, help="number of threads to pass to BLAST, default 1")
    parser.add_argument("--max_threads", type=int, default = 1, help="maximum number of threads to pass to BLAST, default 1.\nBLAST will divide max_threads by threads to determine the number of concurrent blast searches to run.")
    parser.add_argument("--force_ignore", action="store_true", help = "Rebuild blast databases even if they already exist.")
    args = parser.parse_args()
    
    # total CPU budget = max_threads, BLAST per-job CPUs = threads
    max_concurrent = max(1, args.max_threads // args.threads)
    
    # make the transcript file a dictionary
    transcript_dict = SeqIO.to_dict(SeqIO.parse(args.transcripts, "fasta"))

    os.makedirs(args.outdir, exist_ok=True)
    # generate the query files
    print("Loading the genome....")
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    print("parsing the gff3 file for mRNAs")
    mrnas = parse_gff3_mrnas(args.gff3)
    print("extracting mRNAs with flanking sequence from bedfiles")
    for bedfile in sorted(os.listdir(args.bedfile_dir)):
        if not bedfile.endswith(".bed"):
            continue
        bed_path = os.path.join(args.bedfile_dir, bedfile)
        print(f"Processing {bedfile} for mRNAs. . . .")
        region_seqs = extract_flanked_sequences(bed_path, mrnas, genome, args.flank_up, args.flank_down)
        
        for region_name, records in region_seqs.items():
            output_path = os.path.join(args.outdir, f"{region_name}_flanked_mRNAs.fasta")
            SeqIO.write(records, output_path, "fasta")
            print(f"wrote: {output_path} yay!")

    # load the genome list
    print("loading the genome prefix list")
    with open(args.pangenome_list) as pangenome_list:
        genome_names = [os.path.splitext(os.path.basename(line.strip()))[0]
            for line in pangenome_list if line.strip()]

    # calculating max targets
    print("calculating max_targets for BLAST to report")
    max_targets = len(genome_names) * 2
    
    # send the contents of args.outdir to a list so that each entry can have a separate job scheduled
    region_files = sorted(f for f in os.listdir(args.outdir) if f.endswith("_flanked_mRNAs.fasta"))

    print(f"running analysis on {max_concurrent} regions at a time")
    with ThreadPoolExecutor(max_workers=max_concurrent) as ex:
        futures = [ex.submit(process_region, args.outdir, region_file, args.pangenome_dir, args.blast_db, 
            args.threads, args.min_cov, transcript_dict, genome_names, args.force_ignore)
                for region_file in region_files]
        for f in as_completed(futures):
            f.result()


if __name__ == "__main__":
    main()

