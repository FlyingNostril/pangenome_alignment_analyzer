import argparse
import os
import subprocess
import time
import sys
import shutil

# big wrapper script to run the marker to gene alignment pipeline
def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Pangenome pipeline driver script")

    parser.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    parser.add_argument("-R", "--pangenome_region_flanking", type=int, default=2500, help="flanking bp to help identify the regions in the pangenome that correspond to the reference regions/markers.")
    parser.add_argument("-p", "--pangenome_list", required=True, help="Text file listing pangenome FASTA file names (one per line)")
    parser.add_argument("-g", "--gff3", required=True, help="Reference GFF3 annotation file")
    parser.add_argument("-o", "--outdir", default="pipeline_out", help="Top-level output directory")
    parser.add_argument("-D", "--blast_db_loc", required=True, help="FULL PATH to the location of blast data bases for pangenome members")
    parser.add_argument("-u", "--flank_up", type=int, default=2000, help="Flank size for 5' region edges (default 2000bp)\nIs also used in gene-id input mode")
    parser.add_argument("-d", "--flank_down", type=int, default=1000, help="Flank size for 3' region edges (default 1000bp)")
    parser.add_argument("-G", "--pangenome_dir", required=True, help="directory containing the genome files for the pangenome")
    parser.add_argument("-m", "--min_cov", type=float, default=0.8, help="minimum coverage required for blast hit matches")
    parser.add_argument("-t", "--transcripts", required=True, help="transcripts of the genes.  The CDS file of the reference genome")
    parser.add_argument("-T", "--threads", type=int, default=1, help="threads for BLAST jobs")
    parser.add_argument("-M", "--max_threads", type=int, default = 1, help = "Maximum number of threads alloted to the program, if -M > -T then multiple BLAST jobs will be run concurrently.")
    parser.add_argument("-F", "--region_flank", type=int, default=50000, help="Flanking for the region/marker position. Default is 50000 bp, only used with --csvfile")
    parser.add_argument("-I", "--force_ignore", action="store_true", help="set this flag to force the mrna extractor to use stale blast databases for pangenome members instead of computing new ones.")
    parser.add_argument("--max_targets", type=int, default=5, help="max targets allowed for the pangenome region blast run. default 5. reduce this amount if this part of the script runs slowly.")

    # Muscle5 options
    parser.add_argument("--muscle_jobs", type=int, default=1, help="maximum jobs that muscle is allowed to run concurrently. Set with care.\nMuscle can draw large amounts of RAM without warning.")
    parser.add_argument("--muscle_threads", type=int, default=1, help="Maximum number of threads allowed per muscle job. Set with care.\nAdditional threads can lead to drastic increases in RAM draw.")
    parser.add_argument("--super5", action="store_true", help="Use when aligning large amounts of sequence with muscle. reduces memory use.")

    # mutually exclusive pathways: supplied csv file, bed file, or gene id file
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-c", "--csvfile", help="The csv file with the regions/markers.  Must be at least three columns, column 1 is region/marker name, column two is the chromosome in numeric format, column three is the position in a single base pair. The file MUST also have a header line.")
    group.add_argument("-b", "--bedfile", help="bed file to use directly, skips bedfile generation")
    group.add_argument("-i", "--gene_id_file", help="File with gene IDs that matches IDs of the reference genome. One gene ID per line.\nSkips CSV mode and bedfile mode")

    args = parser.parse_args()

    if args.max_threads < args.muscle_jobs * args.muscle_threads:
        parser.error("--max_threads cannot be less than --muscle_jobs * --muscle_threads")

    os.makedirs(args.outdir, exist_ok=True)

    logfile = open(os.path.join(args.outdir, "pipeline.log"), "w")
    sys.stdout = sys.stderr = logfile

    print("Pipeline driver initialized.")
    print(f"Reference genome: {args.reference}")
    print(f"Pangenome list: {args.pangenome_list}")
    print(f"GFF3: {args.gff3}")
    print(f"Output directory: {args.outdir}")
    
    # set base name
    if args.bedfile:
        base = os.path.splitext(os.path.basename(args.bedfile))[0]
    elif args.csvfile: 
        base = os.path.splitext(os.path.basename(args.csvfile))[0]
    else:
        base = os.path.splitext(os.path.basename(args.gene_id_file))[0]

    bed_dir = os.path.join(args.outdir) # removed , base to make the bed_dir the main dir. 
    os.makedirs(bed_dir, exist_ok=True)
    bedfile_path = args.bedfile or os.path.join(bed_dir, f"{base}.bed")

    #copy the bedfile to the bedfile directory
    if args.bedfile: 
        shutil.copy(args.bedfile, bed_dir)

    # create the bedfiles from the csv file of regions/markers (optional)
    if args.bedfile:
        print(f"using the provided bed file: {args.bedfile}")
    elif args.csvfile:
        print(f"creating the bed files from {args.csvfile}")
        subprocess.run([
            "python3", "generate_bedfile3.py", "-c", args.csvfile, "-b", bedfile_path, "-f", str(args.region_flank), "--merge"], check=True)
    else:
        print(f"creating a bed file from the provided gene ID list: {args.gene_id_file}")
        subprocess.run([
            "python3", "generate_bedfile3.py", "-g", args.gene_id_file, "-b", bedfile_path, "-F", str(args.flank_up), "-G", args.gff3], check=True)

    
    # initializing the name of the fasta output file for extract_bed_regions.py
    output_fasta = os.path.join(args.outdir, f"{base}.fasta")

    # create the fasta file of reference sequence
    print(f"creating the fasta file of reference sequence: {output_fasta}")
    subprocess.run([
        "python3", "extract_bed_regions.py", "-b", bedfile_path, "-g", args.reference, "-o", output_fasta],
        check=True)

    # extract the flanking regions for each region in the fasta file
    print("Extracting region flanking sequences from the reference")

    # initialize the variables and give them names/locations
    flank_fasta = os.path.join(args.outdir, f"{base}_flanking_sequence.fasta")

    subprocess.run([
        "python3", "extract_ends_from_fasta.py", "-i", output_fasta, "-o", flank_fasta, 
        "-f", str(args.pangenome_region_flanking)], check=True)

    # set up the blast hits location
    blast_location = os.path.join(args.outdir, "blast_hits_for_pangenome_regions")
    
    #initialize the log file for the pangenome region construction
    flank_hits_log = os.path.join(args.outdir, "find_pangenome_flank_hits.log")
    
    find_pangenome_flank_cmd = ["python3", "find_pangenome_flank_hits2.py", "-i", flank_fasta, "-l", 
            args.pangenome_list, "-g", args.pangenome_dir, "-o", blast_location, "-b", args.blast_db_loc, 
            "-t", str(args.threads), "-T", str(args.max_threads), "-m", str(args.max_targets)]

    # now to run the blast to get the pangenome regions
    print(f"running blast to get the pangenome regions")
    with open(flank_hits_log, "w") as lf:
        subprocess.run(find_pangenome_flank_cmd, check = True, stdout=lf, stderr=subprocess.STDOUT, text=True) 
    print(f"wrote the find_pangenome_flank_hits log to {flank_hits_log}")

    # initializing the pangenome region fasta output directory
    pangenome_fasta_outdir = os.path.join(args.outdir, "pangenome_regions")

    # creating the pangenome region fastas
    print("creating the pangenome region fastas")
    subprocess.run([
        "python3", "generate_pangenome_regions2.py", "-b", blast_location, 
        "-g", args.pangenome_dir, "-o", pangenome_fasta_outdir, "-t", str(args.max_threads)], check=True)
    
    # initialize the blast db and flanked mrnas outputs for extract_flanked_mrnas.py
    blast_db_path = os.path.join(args.outdir, "blast_db")
    flanked_mrnas = os.path.join(args.outdir, "flanked_mrnas")
    
    # run the extract_flanked_mrnas.py script
    print("extracting flanked mrnas from pangenome")
    extract_cmd = ["python3", "extract_flanked_mrnas2.py", "-g", args.reference, "-l", args.pangenome_list,
        "-G", args.gff3, "-o", flanked_mrnas, "-u", str(args.flank_up), "-d", str(args.flank_down), "-B",  bed_dir,
        "-p", pangenome_fasta_outdir, "-D", blast_db_path, "-m", str(args.min_cov), "-t", args.transcripts,
        "--threads", str(args.threads), "--max_threads", str(args.max_threads)]

    if args.force_ignore:
        extract_cmd.append("--force_ignore")

    subprocess.run(extract_cmd, check=True)

    # define the output directory for the alignment program
    aligned_path = os.path.join(args.outdir, "aligned_mrnas")

    # run the alignment script
    print("generate the alignments for each region/marker")
    align_cmd = ["python3", "generate_alignments2.py", "-d", flanked_mrnas, "-o",
        aligned_path, "-t", str(args.max_threads), "--max_jobs", str(args.muscle_jobs), 
        "--max_threads", str(args.muscle_threads)]

    if args.gene_id_file: 
        align_cmd.extend(["-g", args.gene_id_file])
    if args.super5: 
        align_cmd.append("--super5")

    subprocess.run(align_cmd, check=True) 
    
    # define the gene id text files output directory
    ref_transcripts = os.path.join(args.outdir, f"{base}_reference_transcripts")

    # create the gene id txt files
    print("create the gene id lists for each marker/region")
    subprocess.run([
        "python3", "extract_reference_transcripts.py", "-b", bedfile_path, "-g", args.reference, 
        "-G", args.gff3, "-o", ref_transcripts], check=True)

    # run the alignment classifer
    print("classifying the alignments")
    classifier_result = subprocess.run([
        "python3", "evaluate_alignments.py", "-m", aligned_path, "-g", ref_transcripts], check=True, capture_output=True, text=True)
    
    print(classifier_result.stdout.strip())
    if classifier_result.stderr:
        print(classifier_result.stderr.strip())

    print("Alignments created and classified.")
    
    # this step may not be needed. 
    if args.gene_id_file: 
        print("condensing gene ID based results.")
        subprocess.run([
            "python3", "condense_gene_id_results2.py", "-i", aligned_path], check=True)

    end_time = time.time()
    elapsed_sec = end_time - start_time
    print(f"Total time elapsed: {elapsed_sec:.2f} seconds")


    print("Pipeline complete")
    logfile.close()

    # next steps to come here

if __name__ == "__main__":
    main()
