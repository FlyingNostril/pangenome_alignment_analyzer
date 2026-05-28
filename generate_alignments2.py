#!/usr/bin/env python3

import os
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_muscle(input_fasta, output_fasta, threads, log_dir, super5=False):

    base_name = os.path.splitext(os.path.basename(output_fasta))[0]
    stdout_log = os.path.join(log_dir, base_name + ".stdout.log")
    stderr_log = os.path.join(log_dir, base_name + ".stderr.log")
    if not super5:
        cmd = ["muscle", "-align", input_fasta, "-output", output_fasta, "-threads", str(threads)]
    else:
        cmd = ["muscle", "-super5", input_fasta, "-output", output_fasta, "-threads", str(threads)]
    
    with open(stdout_log, "w") as out_handle, open(stderr_log, "w") as err_handle:
        try:
            subprocess.run(cmd, check=True, stdout=out_handle, stderr=err_handle)
            return {"status": "Success", "input": input_fasta, "output": output_fasta,
                    "returncode": 0, "stdout_log": stdout_log, "stderr_log": stderr_log,
                    "cmd": " ".join(cmd)
                    }
        except subprocess.CalledProcessError as e:
            return {"status": "Failed", "input": input_fasta, "output": output_fasta, 
                    "returncode":e.returncode, "stdout_log": stdout_log, 
                    "stderr_log": stderr_log, "cmd": " ".join(cmd)
                    }

def gather_fasta_jobs(master_dir, output_root, extension, gene_id_file=None):
    jobs = []

    if not gene_id_file:
        for subdir in sorted(os.listdir(master_dir)):
            subdir_path = os.path.join(master_dir, subdir)
            if not os.path.isdir(subdir_path):
                continue

            output_subdir = os.path.join(output_root, f"{subdir}_alignments")
            os.makedirs(output_subdir, exist_ok=True)

            for file in sorted(os.listdir(subdir_path)):
                if file.endswith(".fasta") or file.endswith(".fa"):
                    input_fasta = os.path.join(subdir_path, file)
                    base_name = os.path.splitext(file)[0]
                    output_fasta = os.path.join(output_subdir, base_name + extension)
                    jobs.append((input_fasta, output_fasta))
    else:
        input_set = load_gene_ids(gene_id_file)
        
        # putting all of the outputs dont work good.
        #outdir_name = os.path.splitext(os.path.basename(gene_id_file))[0] 
        #output_subdir = os.path.join(output_root, outdir_name)
        #os.makedirs(output_subdir, exist_ok=True)
        
        for subdir in sorted(os.listdir(master_dir)):
            subdir_path = os.path.join(master_dir, subdir)
            if not os.path.isdir(subdir_path):
                continue
            if subdir not in input_set:
                continue

            output_subdir = os.path.join(output_root, f"{subdir}_alignments")
            os.makedirs(output_subdir, exist_ok=True)

            for file in sorted(os.listdir(subdir_path)):
                if not (file.endswith(".fasta") or file.endswith(".fa")):
                    continue

                name = file.split("_CHR_")[0]
                
                if name != subdir:
                    continue
                
                input_fasta = os.path.join(subdir_path, file)
                base_name = os.path.splitext(file)[0]
                output_fasta = os.path.join(output_subdir, base_name + extension)
                jobs.append((input_fasta, output_fasta))
    return jobs

def load_gene_ids(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip() and not line.startswith("#")}
    

def main():
    parser = argparse.ArgumentParser(description="Align all FASTA files in region-specific subdirectories using MUSCLE.")
    parser.add_argument("-d", "--directory", required=True, help="Master directory of region-grouped FASTA records")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory to store alignment subdirectories")
    parser.add_argument("-g", "--gene_id_file", help="Use this option when providing a gene id file to prevent subdirectory bloat and excessive duplication of alignments")
    parser.add_argument("-t", "--total_threads", type=int, default=1, help="max_threads passed to the pipeline or total threads available for this job.\nMust be at least equal to max_jobs * max_threads")
    parser.add_argument("--max_jobs", type=int, default=1, help="maximum number of concurrent muscle jobs allowed. Set with care, Muscle 5 draws huge amounts of RAM.")
    parser.add_argument("--max_threads", type=int, default=1, help="Maximum number of threads allotted to each muscle job. Set with car.\nMuscle 5 RAM draw for each extra thread is non-linear and large.")
    parser.add_argument("--super5", action="store_true", help="Muscle5's option to reduce memory use on large alignments")

    args = parser.parse_args()
    
    if args.total_threads < args.max_jobs * args.max_threads:
        parser.error("--total_threads cannot be less than --max_jobs * --max_threads")

    max_threads = args.max_threads
    max_jobs = args.max_jobs

    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a valid directory.")
        exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    log_dir = os.path.join(args.outdir, "muscle_logs")
    os.makedirs(log_dir, exist_ok=True)
    failure_log = os.path.join(args.outdir, "alignment_failures.log")

    file_extension = ".fasta"

    jobs = gather_fasta_jobs(args.directory, args.outdir, file_extension, args.gene_id_file)
    if not jobs:
        print("No FASTA files found in the provided directory structure.")
        exit(1)

    print(f"Found {len(jobs)} FASTA files to align using {args.max_threads} thread(s)...")
    with open(failure_log, "w") as fail_log:
        fail_log.write("input\toutput\treturncode\tstderr_log\tstdout_log\tcmd\n")

        with ProcessPoolExecutor(max_workers=max_jobs) as executor:
            future_to_job = {executor.submit(
                run_muscle, inp, out, max_threads, log_dir, args.super5
                ): (inp, out) for inp, out in jobs}

            for future in as_completed(future_to_job):
                inp, out = future_to_job[future]
                result = future.result()
                print(result)

                if result["status"] == "Failed":
                    fail_log.write(
                        f"{result['input']}\t"
                        f"{result['output']}\t"
                        f"{result['returncode']}\t"
                        f"{result['stderr_log']}\t"
                        f"{result['stdout_log']}\t"
                        f"{result['cmd']}\n"
                    )
                    fail_log.flush()

if __name__ == "__main__":
    main()

