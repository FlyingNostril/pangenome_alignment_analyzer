#!/usr/bin/env python3

import os
import argparse
from extract_CDS import coding_regions
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from shutil import move

def classify_alignments(master_dir, gene_id_dir):
    four_or_more_protein_changes = 0
    three_or_fewer_protein_changes = 0
    silent_mutations = 0
    noncoding_mutations = 0
    no_mutations = 0
    for region_folder in os.listdir(master_dir):
        region_path = os.path.join(master_dir, region_folder)
        if not os.path.isdir(region_path):
            continue
        region_name = region_folder.replace("_alignments", "")
        gene_id_file = os.path.join(gene_id_dir, f"{region_name}_gene_IDs.txt")
        if not os.path.exists(gene_id_file):
            print(f"Gene ID list missing for {region_name}, skipping.")
            continue
        with open(gene_id_file) as f:
            gene_ids = set(line.strip() for line in f if line.strip())

        dna_changed_dir = os.path.join(region_path, "changed_DNA")
        protein_changed_dir = os.path.join(region_path, "changed_protein")
        protein_minor_dir = os.path.join(protein_changed_dir, "3_or_fewer")
        protein_major_dir = os.path.join(protein_changed_dir, "4_or_more")
        unchanged_dir = os.path.join(region_path, "unchanged")
        os.makedirs(dna_changed_dir, exist_ok=True)
        os.makedirs(protein_changed_dir, exist_ok=True)
        os.makedirs(unchanged_dir, exist_ok=True)
        os.makedirs(protein_minor_dir, exist_ok=True)
        os.makedirs(protein_major_dir, exist_ok=True)

        for fname in os.listdir(region_path):
            if not (fname.endswith(".fasta") or fname.endswith(".fa")):
                continue
            aln_path = os.path.join(region_path, fname)
            if os.path.getsize(aln_path) == 0:
                print(f"Alignment file is empty in {fname}, skipping.")
                continue
            try:
                alignment = AlignIO.read(aln_path, "fasta")
            except ValueError as err:
                print(f"Could not read alignment in {fname}, skipping. {err}")
                continue
            ref_record_cds = None
            ref_record_full = None
            for record in alignment:
                if record.id in gene_ids:
                    ref_record_cds = record
                elif any(gid in record.id for gid in gene_ids) and ref_record_full is None:
                    ref_record_full = record
            if ref_record_cds is None:
                print(f"No reference match for CDS in {fname}, skipping.")
                continue
            if ref_record_full is None:
                print(f"No reference record match for CDS with flanking sequence in {fname}, skipping.")

            ref_seq = str(ref_record_cds.seq.upper())
            ref_seq_full = str(ref_record_full.seq.upper())
            
            CDS_log_path = os.path.join(region_path, fname.replace(".fasta", "_CDS_errors.log"))

            coding_positions, noncoding_positions, CDS_errors = coding_regions(ref_seq, ref_seq_full)
            with open(CDS_log_path, "a") as log:
                if CDS_errors != []:
                    log.write(f"positions are relative to the alignment FASTA\n")
                for error_pos,errors in CDS_errors:
                    log.write(f"{record.id} has an error at: {error_pos}   {errors}\n")
            
            actual_CDS = ''.join(ref_seq[pos] for pos in coding_positions if ref_seq[pos] !='-')
            
            if len(actual_CDS) % 3 != 0:
                print(f"Reference coding sequence in {fname} is not divisable by 3, skipping! Hah!")
                print(f"here is the length of the coding sequence that failed this test: {len(actual_CDS)}")
                print(f"here is the ref_seq:\n{ref_seq}")
                print(f"here is the ref_coding:\n{actual_CDS}")
                continue

            ref_protein = str(Seq(actual_CDS).translate(to_stop=True))

            protein_changed_count = 0
            noncoding_CDS_differences = 0
            untranslated_differences = 0
            protein_records = [SeqRecord(Seq(ref_protein), id=ref_record_cds.id, description = "reference")]

            for record in alignment:
                diff_log_path = os.path.join(region_path, fname.replace(".fasta", "_differences.log"))
                if record.id == ref_record_cds.id:
                    continue
                elif ref_record_full and record.id == ref_record_full.id:
                    continue
                with open(diff_log_path, "a") as log:
                    other_seq = str(record.seq)
                    other_positions = [i for i, base in enumerate(other_seq)]
                    # testing the record vs the two references.  
                    other_coding = ''.join(other_seq[pos] for pos in coding_positions if other_seq[pos] !='-')
                    other_noncoding = ''.join(other_seq[pos] for pos in noncoding_positions)
                    other_protein = str(Seq(other_coding).translate(to_stop=True))
                    protein_records.append(SeqRecord(Seq(other_protein), id=record.id, description="translated"))

                    if other_coding != actual_CDS:
                        for i, (a1, a2) in enumerate(zip(actual_CDS, other_coding)):
                            if a1 != a2:
                                log.write(f"Nucleotide difference in: gene={ref_record_cds.id}    \
                                        pangenome_id={record.id}    pos={i} ref={a1}    alt={a2}\n")
                        noncoding_CDS_differences += 1

                    if other_protein != ref_protein:
                        for i, (a1, a2) in enumerate(zip(ref_protein, other_protein)):
                            if a1 != a2:
                                log.write(f"Protein difference in: gene={ref_record_cds.id}    \
                                        pangenome_id={record.id}    codon={i}   codon_pos={i*3}   ref={a1}    alt={a2}\n")
                        protein_changed_count += 1

                    if other_noncoding != noncoding_positions:
                        untranslated_differences += 1
                        for i, (a1, a2) in enumerate(zip(noncoding_positions, other_noncoding)):
                            if a1 != a2:
                                log.write(f"Nucleotide difference in non-coding region: gene={ref_record_cds.id}    \
                                        pangenome_id={record.id}    pos={i} ref={a1}    alt={a2}\n")
            
            if protein_changed_count >= 4:
                dest = protein_major_dir
                four_or_more_protein_changes += 1
            elif protein_changed_count <= 3 and protein_changed_count >= 1:
                dest = protein_minor_dir
                three_or_fewer_protein_changes += 1
            elif noncoding_CDS_differences >= 1:
                dest = dna_changed_dir
                silent_mutations += 1
            elif untranslated_differences >=1:
                dest = dna_changed_dir
                noncoding_mutations += 1
            else:
                dest = unchanged_dir
                no_mutations += 1
            
            prot_out_path = os.path.join(dest, fname.replace(".fasta", "_proteins.fasta"))
            with open(prot_out_path, "w") as prot_out:
                SeqIO.write(protein_records, prot_out, "fasta")

            move(aln_path, os.path.join(dest, fname))
            print(f"Moved {fname} to {dest}")
            print(f"Wrote translated proteins to {prot_out_path}")
    return four_or_more_protein_changes, three_or_fewer_protein_changes, silent_mutations, noncoding_mutations, no_mutations
def main():
    parser = argparse.ArgumentParser(description="Classify alignments by synonymous and non-synonymous coding sequence differences.")
    parser.add_argument("-m", "--master_dir", help="Directory with region alignment subdirectories")
    parser.add_argument("-g", "--gene_id_dir", help="Directory with region_name_gene_IDs.txt files")
    args = parser.parse_args()

    four_or_more, three_or_less, silent, noncoding, no_mutations = classify_alignments(args.master_dir, args.gene_id_dir)

    print(f"{four_or_more} genes have 4 or more pangenome members with different protein sequences\n \
            {three_or_less} genes have 3 or less pangenome members with different protein sequences\n \
            {silent} genes have silent CDS changes\n \
            {noncoding} genes have changes only in non-coding UTR or Introns\n \
            {no_mutations} genes have no differences at all.\n\n\n\nyay!")


if __name__ == "__main__":
    main()
