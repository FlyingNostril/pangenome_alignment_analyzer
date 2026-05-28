#!/usr/bin/env python3

import argparse
import os
import shutil


def classify_source_path(path_parts, filename):
    if "changed_DNA" in path_parts:
        return "changed_dna"
    if "unchanged" in path_parts:
        return "unchanged"
    if "changed_protein" in path_parts and "3_or_fewer" in path_parts:
        return "protein_3_or_less"
    if "changed_protein" in path_parts and "4_or_more" in path_parts:
        return "protein_4_or_more"
    if filename.endswith(".log"):
        return "logs"
    return None


def transfer_file(src, dst, move_files):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if move_files:
        shutil.move(src, dst)
    else:
        shutil.copy2(src, dst)


def condense(aligned_dir, condensed_dir, move_files=False):
    categories = [
        "changed_dna",
        "protein_3_or_less",
        "protein_4_or_more",
        "unchanged",
        "logs",
        "duplicates",
    ]
    for cat in categories:
        os.makedirs(os.path.join(condensed_dir, cat), exist_ok=True)

    copied = 0
    duplicates = 0

    for entry in sorted(os.listdir(aligned_dir)):
        region_dir = os.path.join(aligned_dir, entry)
        if not os.path.isdir(region_dir) or not entry.endswith("_alignments"):
            continue

        region_name = entry[:-11]  # remove "_alignments"
        for root, _, files in os.walk(region_dir):
            path_parts = set(root.split(os.sep))
            for fname in files:
                source = os.path.join(root, fname)
                category = classify_source_path(path_parts, fname)
                if category is None:
                    continue

                target = os.path.join(condensed_dir, category, fname)
                if os.path.exists(target):
                    duplicates += 1
                    dup_name = f"{region_name}__{fname}"
                    dup_target = os.path.join(condensed_dir, "duplicates", dup_name)
                    transfer_file(source, dup_target, move_files)
                else:
                    transfer_file(source, target, move_files)
                copied += 1

    return copied, duplicates


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Condense gene-ID mode classified alignment outputs into a flat summary "
            "layout with duplicate quarantine."
        )
    )
    parser.add_argument(
        "-i",
        "--input_dir",
        required=True,
        help="Path to aligned_mrnas directory containing *_alignments outputs.",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        default=None,
        help="Condensed output directory (default: <input_dir>/condensed).",
    )
    parser.add_argument(
        "--move",
        action="store_true",
        help="Move files instead of copying them.",
    )
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = (
        os.path.abspath(args.output_dir)
        if args.output_dir
        else os.path.join(input_dir, "condensed")
    )
    os.makedirs(output_dir, exist_ok=True)

    copied, duplicates = condense(input_dir, output_dir, move_files=args.move)
    action = "Moved" if args.move else "Copied"
    print(f"{action} files: {copied}")
    print(f"Duplicate-name files quarantined: {duplicates}")
    print(f"Condensed output: {output_dir}")


if __name__ == "__main__":
    main()
